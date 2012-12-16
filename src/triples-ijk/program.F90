
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
  integer :: nu3gpu

  double precision, allocatable :: eh(:), ep(:), t1(:), v1(:)
  double precision :: mpi_wtime
  double precision start_wall, stop_wall

  character*1 addr(1)


  call ddi_init
  call ddi_memory(0,1000,0)
  call ddi_nproc(ddi_np,ddi_me)
  call ddi_nnode(ddi_nn,ddi_my)
  call ddi_smp_nproc(smp_np,smp_me)
  
! read in the restart file
  if(ddi_me.eq.0) then
    input = 80
    open(unit=input, file='triples.restart', status='old', action='read',form='unformatted',access='sequential')
    read(input) no
    read(input) nu
    write(6,*) 'no=',no
    write(6,*) 'nu=',nu
  end if
  call ddi_bcast(123,'f',no,1,0)
  call ddi_bcast(124,'f',nu,1,0)

  call ddi_sync(1)

  call common_cc_init(no,nu)

! allocation replicated storage
  allocate( eh(no) )
  allocate( ep(nu) )
  allocate( v1(nou) )
  allocate( t1(nou) )

! the following shared quanties are shared between
! both the CPU code and GPU code
  nu3gpu = nu3+nu3*gpu_nd
  if(ddi_me.eq.0) then
     write(6,*) 'nu3   =',nu3
     write(6,*) 'nu3gpu=',nu3gpu
  endif
  call ddi_smp_create(no2u2 ,d_t2)
  call ddi_smp_create(no3u  ,d_vm)
  call ddi_smp_create(no2u2 ,d_voe)
  call ddi_smp_create(nu3   ,d_t3)  ! only needed for iij/ijj tuples
  call ddi_smp_create(nu3   ,d_v3)
  call ddi_smp_create(nu3gpu,d_vei)
  call ddi_smp_create(nu3gpu,d_vej)
  call ddi_smp_create(nu3   ,d_vek)

  call ddi_smp_offset(d_t2, addr, lt2)
  call ddi_smp_offset(d_vm, addr, lvm)
  call ddi_smp_offset(d_voe,addr, lvoe)
  call ddi_smp_offset(d_t3 ,addr, lt3) 
  call ddi_smp_offset(d_v3 ,addr, lv3)
  call ddi_smp_offset(d_vei,addr, lvei)
  call ddi_smp_offset(d_vej,addr, lvej)
  call ddi_smp_offset(d_vek,addr, lvek)

! set offsets
  lt2  = lt2+1
  lvm  = lvm+1
  lvoe = lvoe+1
  lt3  = lt3+1
  lv3  = lv3+1
  lvei = lvei+1
  lvej = lvej+1
  lvek = lvek+1

! each rank that drive a gpu will have it's own offset to the 
! non-cpu arrays; otherwise, gve{ijk} will point to the cpu storage
if(smp_me.lt.gpu_nd) then
   if(ddi_my.eq.0) write(6,*) 'me=',ddi_me,' lvei=',lvei
   lvei = lvei + nu3*(smp_me+1)
   lvej = lvej + nu3*(smp_me+1)
   if(ddi_my.eq.0) write(6,*) 'me=',ddi_me,' gvei=',lvei
endif

  call ddi_create(nutr,nou,d_vvvo)

! if we want to conserve shared memory on the node we can put
! t2 and voe into distributed memory - this will require more 
! remote memory accesses ==> give and take scenario

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

  call cc_triples_acc(eh,ep,v1,t1,addr(lt2),addr(lv3),addr(lt3), &
               addr(lvm),addr(lvoe),addr(lvei),addr(lvej),addr(lvek),d_vvvo)
      

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

!----------------------------------------------------------

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

!----------------------------------------------------------

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


