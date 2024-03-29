subroutine cc_triples_acc(eh,ep,v1,t1,t2,v3,t3,vm,voe,vvvo_hnd,gms_ets,gms_etd)
use common_cc
implicit none

double precision, allocatable :: ve_i(:), ve_j(:), ve_k(:)
double precision, allocatable :: tmp_i(:), tmp_j(:), tmp_k(:)
double precision, allocatable :: tmp(:), vei(:), vej(:)
double precision :: eh(no),ep(nu),v1(nou),t1(*),t2(*),v3(*),t3(*)
double precision :: vm(*),voe(*)
double precision :: gms_ets, gms_etd
integer :: vvvo_hnd

integer iii

integer :: nr, sr, iwrk, i, j, k, mytask, divisor, partial, icntr
integer :: iold, jold, kold, comm_core
integer :: ilo,ihi


integer :: n_ijk_tuples, n_iij_tuples, n_ijj_tuples
integer :: ierr

double precision ddot
double precision :: mpi_wtime
double precision :: ijk_start, ijk_stop

integer ioff, joff, koff, iloop, jloop, kloop, ij

d_vvvo = vvvo_hnd

etd = zero
v1(1:nou) = 0.0D+00
call ddi_dlbreset()
call ddi_sync(1234)

if(smp_me.ne.0) goto 999

#ifdef USE_OPEN_ACC
  allocate(ve_i(nu3),ve_j(nu3),ve_k(nu3))
  allocate(tmp(nutr*nu),tmp_i(nutr*nu),tmp_j(nutr*nu),tmp_k(nutr*nu))
#else
  allocate(tmp(nutr*nu),tmp_i(nutr*nu),tmp_j(nutr*nu),tmp_k(nutr*nu))
! allocate( tmp(nu2) )
#endif


! load-balancing strategy
! =======================


n_ijk_tuples = (no*(no-1)*(no-2))/6
n_ijj_tuples = (no*(no-1))/2
n_iij_tuples = n_ijj_tuples

if(ddi_me.eq.0) then
   write(6,*) 'ijk tuples = ',n_ijk_tuples
   write(6,*) 'iij tuples = ',n_iij_tuples
   write(6,*) 'ijj tuples = ',n_ijj_tuples
endif

! write(6,/30C=10I/) 'IJK Tuples',n_ijk_tuples
! write(6,/30C=10I/) 'IIJ Tuples',n_iij_tuples

! flops_per_ijk_tuple = some_value
! flops_per_ijj_tuple = some_value
! flops_per_iij_tuple = some_value

! write(6,/30C=10I/) 'IJK Est. Flops',flops_per_ijk_tuple
! write(6,/30C=10I/) 'IIJ Est. Flops',flops_per_iij_tuple

! work ratios between ijk and (iij+ijj)
! note: infate the ratio by 1.25 due to extra communication overhead 
! with respect to flops for iij/ijj
! count_ratio = 1.0 * (n_iij_tuples + n_ijj_tuples) / n_ijk_tuples
! flop_ratio = 1.0 * (flops_per_iij_tuple + flops_per_ijj_tuple) / flops_per_ijk_tuple
! ratio = 1.0 + (count_ratio * flop_ratio *. 1.5)

! ijk_gpu_to_cpu_efficiency is an efficency ratio of 
! 1 node of CPU cores vs 1 GPU
! ijk_efficiency is scaled by the ratio of CPU nodes to GPU devices
! ijk_gpu_to_cpu_efficiency is a user runtime parameter
! ijk_gpu_efficiency = (ijk_gpu_to_cpu_efficiency * gpu_nd) / ( 1.0*ddi_nn ) 


! tuples per gpu ratio (tpgr)
! tpgr = 1.0 * n_ijk_tuples / gnu_nd

! any remaining tumples represent a load-IMBALANCE
! they can be assigned to CPU nodes or GPU nodes depending
! on the GPU:CPU efficency ratio
! tuples_per_gpu = div(n_ijk_tuples,gpu_nd)
! tuples_remaining = mod(n_ijk_tuples,gpu_nd)

! if the tuples per gpu ratio exceeds the ijk_gpu_to_cpu_efficiency 
! ratio, then we need to assign at least some portion of the
! ijk tuples to some subset of the CPU nodes

! if the gpu iterations per cpu iteration (gipci) is equal to 1
! gipci = tpgr / (ijk_gpu_efficiency+1)

call div_even(n_ijk_tuples,ddi_nn,ddi_my,nr,sr)
sr = sr-1


divisor = 10
if(nr.lt.10) divisor = nr
partial = nr/divisor
if(mod(nr,divisor).ne.0) partial = partial+1


iold = -1
jold = -1
kold = -1

!if(smp_np.gt.1) call smp_sync()

#ifndef USE_OPEN_ACC
if(smp_me.eq.0) then
   v3(1:nu3) = 0.0D+00
endif
#endif


!$acc data copyout(v1(1:nou))  &
!$acc& copyin(eh,ep,t2(1:nu2*no*no),vm(1:no*nu*no*no),voe(1:no2u2),t3(1:nu3),etd)  &
!$acc& create(ve_i(1:nu3),ve_j(1:nu3),ve_k(1:nu3),v3(1:nu3),tmp_i,tmp_j,tmp_k,x3)
!$acc wait

do iwrk = 1,10
!$acc update device( ve_i(iwrk:iwrk) ) async(iwrk)
end do
!$acc wait

call dgemm_async_setup(10, 1)

#ifdef USE_OPEN_ACC
!$acc kernels
  v1(1:nou) = 0.0D+00   ! initialize on GPU only
!$acc end kernels
#endif // end USE_OPEN_ACC

if(ddi_me.eq.0) ijk_start = mpi_wtime()


! ---------- ijk-tuples loop ---------------

do iwrk = sr, sr+nr-1
  mytask = iwrk
  call ddcc_t_task(mytask,no,i,j,k)

  comm_core = 0

   ! OpenACC CODE ----------------------------

#ifdef USE_OPEN_ACC
     if(j.ne.jold) then
!      if(smp_me.eq.comm_core) then
          ilo = nu*(j-1) + 1
          ihi = ilo + nu - 1
          call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_j)
!$acc update device(tmp_j(1:nutr*nu)) async(2)
          call ddcc_t_getve_acc(2,nu,j,tmp_j,ve_j)
          call trant3_1_async(2,nu,ve_j)
!      end if
!      comm_core = comm_core+1
!      if(comm_core.eq.smp_np) comm_core=0
     end if

     if(k.ne.kold) then
!      if(smp_me.eq.comm_core) then
          ilo = nu*(k-1) + 1
          ihi = ilo + nu - 1
          call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_k)
!$acc update device(tmp_k(1:nutr*nu)) async(3)
          call ddcc_t_getve_acc(3,nu,k,tmp_k,ve_k)
          call trant3_1_async(3,nu,ve_k)
!      end if
!      comm_core = comm_core+1
!      if(comm_core.eq.smp_np) comm_core=0
     end if

     if(i.ne.iold) then
!      if(smp_me.eq.comm_core) then
          ilo = nu*(i-1) + 1
          ihi = ilo + nu - 1
          call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_i)
!$acc update device(tmp_i(1:nutr*nu)) async(1)
          call ddcc_t_getve_acc(1,nu,i,tmp_i,ve_i)
          call trant3_1_async(1,nu,ve_i)
!      end if
!      comm_core = comm_core+1
!      if(comm_core.eq.smp_np) comm_core=0
     end if

!    if(smp_np.gt.1) call smp_sync()
#endif

     call ddcc_t_ijk_acc(no,nu,i,j,k,v1,t2,vm,v3,voe,eh,ep,ve_i,ve_j,ve_k)


  iold = i
  jold = j
  kold = k

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

! ---------- end of ijk-tuples loop ---------------

998 continue
!call ddi_sync(1234)
if(ddi_me.eq.0) then 
   ijk_stop = mpi_wtime()
   write(6,*) 'ijk time on rank 0 = ',(ijk_stop-ijk_start)
endif
!if(smp_me.ne.0) goto 999

! counters and load-balancing for iij and ijj tuples
icntr = 0

if(smp_me.eq.0) call ddi_dlbnext(mytask)
!call ddi_smp_bcast(1237,'I',mytask,1,0)

!$acc wait

! ----------- iij and ijj tuples -------------
do i=1,no
  do j=1,i-1
  ! ijj tuple
    if(mytask.eq.icntr) then
       ! if(smp_me.eq.0) write(6,*) ddi_me,' tasked with ',mytask,' iij/ijj tuple'
       if(i.ne.iold) then
         if(smp_me.eq.0) then
           ilo = nu*(i-1) + 1
           ihi = ilo + nu - 1
           call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_i)
!$acc update device(tmp_i(1:nutr*nu)) async(1)
           call ddcc_t_getve_acc(1,nu,i,tmp_i,ve_i)
         end if
         call tranmd_23_acc(1,ve_i,nu,nu,nu,1)
         iold = i
       end if
       if(j.ne.jold) then
         if(smp_me.eq.0) then
           ilo = nu*(j-1) + 1
           ihi = ilo + nu - 1
           call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_j)
!$acc update device(tmp_j(1:nutr*nu)) async(2)
           call ddcc_t_getve_acc(2,nu,j,tmp_j,ve_j)
         end if
         call tranmd_23_acc(2,ve_j,nu,nu,nu,1)
         jold = j
       end if
       call ddcc_t_ijj_acc(no,nu,i,j,v1,t2,vm,v3,t3,voe,t1,eh,ep,tmp,ve_i,ve_j)
       if(smp_me.eq.0) call ddi_dlbnext(mytask)
!      call ddi_smp_bcast(1235,'I',mytask,1,0)
    end if
    icntr = icntr + 1
  ! iij tuple
    if(mytask.eq.icntr) then
       ! if(smp_me.eq.0) write(6,*) ddi_me,' tasked with ',mytask,' iij/ijj tuple'
       if(i.ne.iold) then
         if(smp_me.eq.0) then
           ilo = nu*(i-1) + 1
           ihi = ilo + nu - 1
           call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_i)
!$acc update device(tmp_i(1:nutr*nu)) async(1)
           call ddcc_t_getve_acc(1,nu,i,tmp_i,ve_i)
         end if
         call tranmd_23_acc(1,ve_i,nu,nu,nu,1)
         iold = i
       end if
       if(j.ne.jold) then
         if(smp_me.eq.0) then
           ilo = nu*(j-1) + 1
           ihi = ilo + nu - 1
           call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_j)
!$acc update device(tmp_j(1:nutr*nu)) async(2)
           call ddcc_t_getve_acc(2,nu,j,tmp_j,ve_j)
         end if
         call tranmd_23_acc(2,ve_j,nu,nu,nu,1)
         jold = j
       end if
       call ddcc_t_iij_acc(no,nu,i,j,v1,t2,vm,v3,t3,voe,t1,eh,ep,tmp,ve_i,ve_j)
       if(smp_me.eq.0) call ddi_dlbnext(mytask)
!      call ddi_smp_bcast(1234,'I',mytask,1,0)
    end if
    icntr = icntr + 1
  end do
end do
! ----------- end of iij and ijj tuples -------------

!$acc end data       

999 continue
call ddi_sync(1235)
call ddi_gsumf(125,v1,nou)
call ddi_gsumf(126,etd,1)

if(ddi_me.eq.0) then
  call trpose(t1,tmp,no,nu,1)
  ets = 2.0D+00*ddot(nou,t1,1,v1,1)
end if

if(ddi_me.eq.0) then
  write(6,9000) ets,etd,ets+etd
end if

gms_ets = ets
gms_etd = etd

call ddi_sync(2)
! gsum the common block from call dd_t3squa_gsum

! if(ddi_me.eq.0) then
!   call trt1(nu,no,t3,v1)
!   call t1sq(no,nu,t3,v1)
!   call addden1(no,nu,v1,eh,ep)
! end if

return
9000 format(1x,'ets/etd/ets+etd=',3F15.10)
end subroutine cc_triples_acc
