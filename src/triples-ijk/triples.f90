
! =============================================================================
!   Copyright (C) 2010.  Ryan M. Olson
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
! =============================================================================

program cc_triples_restart
use common_cc
implicit none

  integer :: input
  integer :: d_t2, d_vm, d_voe, d_vei, d_vej, d_vek, d_t3, d_v3
  integer :: lt2, lvm, lvoe, lvei, lvej, lvek, lt3, lv3

  double precision, allocatable :: eh(:), ep(:), t1(:), v1(:)
  double precision :: mpi_wtime
  double precision start_wall, stop_wall

  character*1 addr(1)


  call ddi_init
  call ddi_memory(0,1000,0)
  call ddi_nproc(ddi_np,ddi_me)
  call ddi_nnode(ddi_nn,ddi_my)
  call ddi_smp_nproc(smp_np,smp_me)

  if(ddi_me.eq.0) then
    input = 80
    open(unit=input, file='triples.restart', status='old', action='read',form='unformatted',access='sequential')
    read(input) no
    read(input) nu
  end if
  call ddi_bcast(123,'f',no,1,0)
  call ddi_bcast(124,'f',nu,1,0)

  call ddi_sync(1)

  call common_cc_init(no,nu)
 
  allocate( eh(no) )
  allocate( ep(nu) )
  allocate( v1(nou) )
  allocate( t1(nou) )

  call ddi_smp_create(no2u2,d_t2)
  call ddi_smp_create(no3u ,d_vm)
  call ddi_smp_create(no2u2,d_voe)
  call ddi_smp_create(nu3  ,d_t3)
  call ddi_smp_create(nu3  ,d_v3)
  call ddi_smp_create(nu3  ,d_vei)
  call ddi_smp_create(nu3  ,d_vej)
  call ddi_smp_create(nu3  ,d_vek)

  call ddi_smp_offset(d_t2, addr, lt2)
  call ddi_smp_offset(d_vm, addr, lvm)
  call ddi_smp_offset(d_voe,addr, lvoe)
  call ddi_smp_offset(d_t3 ,addr, lt3) 
  call ddi_smp_offset(d_v3 ,addr, lv3)
  call ddi_smp_offset(d_vei,addr, lvei)
  call ddi_smp_offset(d_vej,addr, lvej)
  call ddi_smp_offset(d_vek,addr, lvek)

  lt2  = lt2+1
  lvm  = lvm+1
  lvoe = lvoe+1
  lt3  = lt3+1
  lv3  = lv3+1
  lvei = lvei+1
  lvej = lvej+1
  lvek = lvek+1

  call ddi_create(nutr,nou,d_vvvo)
! call ddi_create(nou,nu2,d_vovv)

  if(ddi_me.eq.0) then
     call cc_triples_readinp(input,eh,ep,t1,addr(lt2),addr(lvm),addr(lvoe),addr(lv3))
  endif

  call ddi_sync(1234)

  call ddi_bcast(1,'F',eh,no,0)
  call ddi_bcast(2,'F',ep,nu,0)
  call ddi_bcast(3,'F',t1,nou,0)
  
  if(smp_me.eq.0) then
    call ddi_masters_bcast(4,'F',addr(lt2),no2u2,0)
    call ddi_masters_bcast(5,'F',addr(lvm),no3u,0)
    call ddi_masters_bcast(6,'F',addr(lvoe),no2u2,0)
  end if

  call ddi_sync(2)

  if(ddi_me.eq.0) then
    write(6,*) 'restarting triples'
  end if

! if(ddi_my.eq.0 .and. smp_me.eq.1) then
!   call cc_test(eh,ep,t1,addr(lt2),addr(lvm),addr(lvoe))
! end if

  call ddi_sync(3)

! if(ddi_me.eq.0) call cc_convert_ijk_to_abc(addr(lvm),addr(lt3),addr(lv3))
  start_wall = mpi_wtime()

  call ddi_sync(4)

  call cc_triples(eh,ep,v1,t1,addr(lt2),addr(lv3),addr(lt3), &
               addr(lvm),addr(lvoe),addr(lvei),addr(lvej),addr(lvek))

  stop_wall = mpi_wtime()

  call ddi_sync(4)

  if(ddi_me.eq.0) write(6,9000) (stop_wall-start_wall)
9000 format('walltime=',F15.5)

  call ddi_destroy(d_vvvo)
  call ddi_smp_destroy(d_vek)
  call ddi_smp_destroy(d_vej)
  call ddi_smp_destroy(d_vei)
  call ddi_smp_destroy(d_v3)
  call ddi_smp_destroy(d_t3)
  call ddi_smp_destroy(d_voe)
  call ddi_smp_destroy(d_vm)
  call ddi_smp_destroy(d_t2)

  call ddi_finalize
end program cc_triples_restart

subroutine cc_test(eh,ep,t1,t2,vm,voe)
use common_cc
implicit none
double precision :: eh(no), ep(nu), t1(nou), t2(no2u2), vm(no3u), voe(no2u2)
      write(6,*) 'check eh(1)=',eh(1)
      write(6,*) 'check eh(no-10)=',eh(no-10)
      write(6,*) 'check ep(1)=',ep(1)
      write(6,*) 'check ep(nu-10)=',ep(nu-10)
      write(6,*) 'check o1(1)=',t1(1)
      write(6,*) 'check o1(nou/2)=',t1(nou/2)
      write(6,*) 'check t2(1)=',t2(1)
      write(6,*) 'check t2(no2u2/2)=',t2(no2u2/2)
      write(6,*) 'check t2(no2u2-1)=',t2(no2u2-1)
      write(6,*) 'check vm(1)=',vm(1)
      write(6,*) 'check vm(no3u/2)=',vm(no3u/2)
      write(6,*) 'check voe(1)=',voe(1)
      write(6,*) 'check voe(no2u2/2)=',voe(no2u2/2)
      flush(6)
return
end subroutine cc_test


subroutine cc_triples_readinp(input,eh,ep,t1,t2,vm,voe,tmp)
use common_cc
implicit none

integer :: i, j, input
double precision :: eh(no), ep(nu), t1(nou), t2(no2u2), vm(no3u), voe(no2u2), tmp(nu3)

read(input,err=911) (eh(i),i=1,no)
read(input,err=911) (ep(i),i=1,nu)
read(input,err=911) (t1(i),i=1,nou)
read(input,err=911) (t2(i),i=1,no2u2)
read(input,err=911) (vm(i),i=1,no3u)
read(input,err=911) (voe(i),i=1,no2u2)

! generate check values
! call cc_test(eh,ep,t1,t2,vm,voe)

do j = 1,nou
  read(input,err=911) (tmp(i),i=1,nutr)
  if(j.eq.1 .or. j.eq.nou) then
    write(6,*) 'check ve_ia(1)=',tmp(1)
  end if
  call ddi_put(d_vvvo,1,nutr,j,j,tmp)
end do
return

911 write(6,*) "error reading restart file"
stop
end subroutine cc_triples_readinp


subroutine cc_triples(eh,ep,v1,t1,t2,v3,t3,vm,voe,ve_i,ve_j,ve_k)
use common_cc
implicit none

double precision, allocatable :: tmp(:)
double precision :: eh(no),ep(nu),v1(nou),t1(*),t2(*),v3(*),t3(*)
double precision :: vm(*),voe(*),ve_i(*),ve_j(*),ve_k(*)
integer :: ntuples, nr, sr, iwrk, i, j, k, mytask, divisor, partial, icntr
integer :: iold, jold, kold, comm_core

double precision ddot
double precision :: mpi_wtime
double precision :: ijk_start, ijk_stop

allocate( tmp(nu2) )

ntuples = (no*(no-1)*(no-2))/6
call div_even(ntuples,ddi_nn,ddi_my,nr,sr)
sr = sr-1

divisor = 10
if(nr.lt.10) divisor = nr
partial = nr/divisor
if(mod(nr,divisor).ne.0) partial = partial+1

etd = zero

iold = -1
jold = -1
kold = -1

v1(1:nou) = 0.0D+00
if(smp_me.eq.0) then
   v3(1:nu3) = 0.0D+00
endif
call smp_sync()
call ddi_dlbreset()
call ddi_sync(1234)

if(ddi_me.eq.0) ijk_start = mpi_wtime()

do iwrk = sr, sr+nr-1
  mytask = iwrk
  call ddcc_t_task(mytask,no,i,j,k)

  comm_core = 0

  if(i.ne.iold) then
    if(smp_me.eq.comm_core) call ddcc_t_getve(nu,i,tmp,ve_i)
    comm_core = comm_core+1
    if(comm_core.eq.smp_np) comm_core=0
  end if

  if(j.ne.jold) then
    if(smp_me.eq.comm_core) call ddcc_t_getve(nu,j,tmp,ve_j)
    comm_core = comm_core+1
    if(comm_core.eq.smp_np) comm_core=0
  end if

  if(k.ne.kold) then
    if(smp_me.eq.comm_core) call ddcc_t_getve(nu,k,tmp,ve_k)
    comm_core = comm_core+1
    if(comm_core.eq.smp_np) comm_core=0
  end if

  if(i.ne.iold) then
    call trant3_1(nu,ve_i)
    iold = i
  end if

  if(j.ne.jold) then
    call trant3_1(nu,ve_j)
    jold = j
  end if

  if(k.ne.kold) then
    call trant3_1(nu,ve_k)
    kold = k
  end if

  call ddcc_t_ijk_big(no,nu,i,j,k,v1,t2,vm,v3,t3,voe,t1,tmp,eh,ep,ve_i,ve_j,ve_k)
  call smp_sync()

  if(ddi_me.eq.0) then
!    write(6,*) 'mytask',mytask
     if(mod(mytask,partial).eq.0) then
        if(mytask/partial .gt. 0) then
          write(6,*) '% complete = ', mytask/partial*divisor
          if(mod(nr,divisor).ne.0) then
            if(mytask/partial .eq. mod(nr,divisor)) partial = partial-1
          end if
        end if
     end if
  end if
end do

call ddi_sync(1234)
if(ddi_me.eq.0) then
   ijk_stop = mpi_wtime()
   if(ddi_me.eq.0) write(6,9001) (ijk_stop-ijk_start)
 9001 format('walltime=',F15.5)
end if

! counters and load-balancing for iij and ijj tuples
icntr = 0

if(smp_me.eq.0) call ddi_dlbnext(mytask)
call ddi_smp_bcast(1234,'I',mytask,1,0)

! iij and ijj tuples
do i=1,no
  do j=1,i-1
  ! ijj tuple
    if(mytask.eq.icntr) then
       if(i.ne.iold) then
         if(smp_me.eq.0) call ddcc_t_getve(nu,i,tmp,ve_i)
         call tranmd_23(ve_i,nu,nu,nu,1)
         iold = i
       end if
       if(j.ne.jold) then
         if(smp_me.eq.0) call ddcc_t_getve(nu,j,tmp,ve_j)
         call tranmd_23(ve_j,nu,nu,nu,1)
         jold = j
       end if
       call ddcc_t_ijj_big(no,nu,i,j,v1,t2,vm,v3,t3,voe,t1,eh,ep,tmp,ve_i,ve_j)
       if(smp_me.eq.0) call ddi_dlbnext(mytask)
       call ddi_smp_bcast(1235,'I',mytask,1,0)
    end if
    icntr = icntr + 1
  ! iij tuple
    if(mytask.eq.icntr) then
       if(i.ne.iold) then
         if(smp_me.eq.0) call ddcc_t_getve(nu,i,tmp,ve_i)
         call tranmd_23(ve_i,nu,nu,nu,1)
         iold = i
       end if
       if(j.ne.jold) then
         if(smp_me.eq.0) call ddcc_t_getve(nu,j,tmp,ve_j)
         call tranmd_23(ve_j,nu,nu,nu,1)
         jold = j
       end if
       call ddcc_t_iij_big(no,nu,i,j,v1,t2,vm,v3,t3,voe,t1,eh,ep,tmp,ve_i,ve_j)
       if(smp_me.eq.0) call ddi_dlbnext(mytask)
       call ddi_smp_bcast(1236,'I',mytask,1,0)
    end if
    icntr = icntr + 1
  end do
end do

call ddi_gsumf(123,eh,no)
call ddi_gsumf(124,ep,nu)
call ddi_gsumf(125,v1,nou)
call dd_t3squa_gsum

call ddi_sync(1)

! there is a problem calculating ets
! etd is correct; ets is wrong - compared against triples-nersc
if(ddi_me.eq.0) then
  call trpose(t1,tmp,no,nu,1)
  ets = 2.0D+00*ddot(nou,t1,1,v1,1)
end if

if(ddi_me.eq.0) then
  write(6,9000) ets,etd,ets+etd
end if

call ddi_sync(2)
! gsum the common block from call dd_t3squa_gsum

! if(ddi_me.eq.0) then
!   call trt1(nu,no,t3,v1)
!   call t1sq(no,nu,t3,v1)
!   call addden1(no,nu,v1,eh,ep)
! end if

return
9000 format(1x,'ets/etd/ets+etd=',3F15.10)
end subroutine cc_triples




subroutine cc_convert_ijk_to_abc(vm,vm_tmp,v3)
use common_cc
implicit none

integer :: i,j,k,a,b,c
integer :: ilo,ihi,jlo,jhi
double precision :: vm(no,nu,no,no), vm_tmp(no,no,no,nu), v3(nu,nu,nu)
double precision, allocatable :: tmp(:,:)

! ----------------------------
! transform vm to look like ve
! ----------------------------

! vm(i,j,2,1) == vm(1,j,2,i)
! vm(i,a,j,k) == vm(k,a,j,i) 
! vm(h1,p?,h3,h2) == vm(h2,p?,h3,h1)
! store as (h1,h2,h3,p?) to be equivalent in storage to v3
! but note, ve requires a 23 sym operation prior to use
! so that we can done on vm just once after the order is arranged.

do k = 1,no       ! h3
  do j = 1,no     ! h2
    do a = 1,nu   ! p?
      do i = 1,no ! h1
        vm_tmp(i,j,k,a) = vm(i,a,k,j)
      end do
    end do
  end do
end do

call dcopy(no3u,vm_tmp,1,vm,1)
call tranmd_23(vm,no,no,no,nu)

! vm is ready to be used like ve
! vm_a = vm(1,1,1,a) and can be used like ve_i
! no fetching or manipulating is required


! ----------------------------
! transform ve to look like vm
! prepare d_vovv from d_vvvo
! ----------------------------

! ve is stored as ve(p4,p5,p6,h?) == ve(p5,p4,p6,h?)
! needs to be transformed to ve(p4,h?,p6,p5)

! d_vvvo is stored as ve(p4,p5,p6,h?) [p4 p5 | p6 h?]
! d_vovv is stored as ve(p4,h?,p6,p5) [p4 h? | p6 p5] = [p4 h? | p5 p6] = <p4 p5 | h? p6>

! d_vovv is not entirely needed; however the operation
! to get a ve(p4,i,a,b) would require NO gets of size
! NU from patch (nu*(a-1)+1:nu*a,b+no*(i-1)+b:b+no*i

allocate( tmp(nu,nu) )

do i = 1,no   ! h?
  
  call ddcc_t_getve(nu,1,tmp,v3) ! v3(p4,p5,p6)

  ilo = (i-1)*nu + 1
  ihi = ilo + nu

  do b = 1, nu  ! p5

    jlo = (b-1)*nu + 1
    jhi = jlo + nu

    do c = 1,nu   ! p6
      do a = 1,nu ! p4
         tmp(a,c) = v3(a,b,c)
      end do
    end do

    call ddi_put(d_vovv,ilo,ihi,jlo,jhi,tmp)

  end do ! end p5
end do ! end h?


deallocate( tmp )

! this is the test loops to used to determine how vm was being used
! test vm do i = 1, no
! test vm   do j = 1, nu
! test vm       write(6,9000) 1,vm(i,j,2,1), vm(2,j,i,1)
! test vm       write(6,9000) 2,vm(i,j,2,1), vm(1,j,i,2)
! test vm       write(6,9000) 3,vm(i,j,2,1), vm(1,j,2,i)
! test vm       write(6,9000) 4,vm(i,j,2,1), vm(2,j,1,i)
! test vm   end do
! test vm end do

return
9000 format(I5,2F20.15)
end subroutine cc_convert_ijk_to_abc

! we can run tests on just the [t], i.e. we can simply calculate etd without doing
! the (t) part or ets part.

! we have to get the ijj and iij routines into shape, but that should not be hard,
! because instead of 6 sets of dgemm transposes, we only have 3 sets.

