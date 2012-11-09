subroutine cc_triples_acc(eh,ep,v1,t1,t2,v3,t3,vm,voe,vvvo_hnd)
use common_cc
implicit none

double precision, allocatable :: ve_i(:), ve_j(:), ve_k(:)
double precision, allocatable :: tmp_i(:), tmp_j(:), tmp_k(:)
double precision, allocatable :: tmp(:), vei(:), vej(:)
double precision :: eh(no),ep(nu),v1(nou),t1(*),t2(*),v3(*),t3(*)
double precision :: vm(*),voe(*)
integer :: vvvo_hnd

integer iii

integer :: nr, sr, iwrk, i, j, k, mytask, divisor, partial, icntr
integer :: iold, jold, kold, comm_core
integer :: ilo,ihi

#ifdef USE_CUDA
integer :: gpu_driver
#endif

integer :: n_ijk_tuples, n_iij_tuples, n_ijj_tuples
integer :: ierr

double precision ddot
double precision :: mpi_wtime
double precision :: ijk_start, ijk_stop

integer ioff, joff, koff, iloop, jloop, kloop, ij

d_vvvo = vvvo_hnd

#ifdef USE_OPEN_ACC
  allocate(ve_i(nu3),ve_j(nu3),ve_k(nu3))
  allocate(tmp(nutr*nu),tmp_i(nutr*nu),tmp_j(nutr*nu),tmp_k(nutr*nu))
#else
  allocate(tmp(nutr*nu),tmp_i(nutr*nu),tmp_j(nutr*nu),tmp_k(nutr*nu))
! allocate( tmp(nu2) )
#endif

! gpu
#ifdef USE_CUDA
gpu_driver = 0
if(smp_me.lt.gpu_nd) gpu_driver=1
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

etd = zero

iold = -1
jold = -1
kold = -1

if(smp_np.gt.1) call smp_sync()
call ddi_dlbreset()

#ifndef USE_OPEN_ACC
v1(1:nou) = 0.0D+00
if(smp_me.eq.0) then
   v3(1:nu3) = 0.0D+00
endif
#endif

! switch scopes
#ifdef USE_CUDA
call ddi_sync(1234)
working_smp_comm = hybrid_smp_comm
working_compute_comm = hybrid_compute_comm
call mpi_comm_rank(working_smp_comm, smp_me, ierr)
call mpi_comm_size(working_smp_comm, smp_np, ierr)
call mpi_comm_rank(working_compute_comm, ddi_me, ierr)
call mpi_comm_size(working_compute_comm, ddi_np, ierr)
call ddi_sync(1234)

if(gpu_driver.eq.1) then
   call ijk_gpu_init(no,nu,eh,ep,v1,t2,vm,voe)
endif
#endif // ifdef USE_CUDA


call ddi_sync(1234)

!$acc data copyout(v1(1:nou))  &
!$acc& copyin(eh,ep,t2(1:nu2*no*no),vm(1:no*nu*no*no),voe(1:no2u2),t3(1:nu3))  &
!$acc& create(ve_i(1:nu3),ve_j(1:nu3),ve_k(1:nu3),v3(1:nu3),tmp_i,tmp_j,tmp_k)

#ifdef USE_OPEN_ACC
!$acc kernels
  v1(1:nou) = 0.0D+00   ! initialize on GPU only
!$acc end kernels
#endif // end USE_OPEN_ACC

if(ddi_me.eq.0) ijk_start = mpi_wtime()

#ifdef USE_CUDA
if(gpu_driver.eq.1) then
allocate( vei(nu3), vej(nu3) )
#endif


! ---------- ijk-tuples loop ---------------

do iwrk = sr, sr+nr-1
  mytask = iwrk
  call ddcc_t_task(mytask,no,i,j,k)

  comm_core = 0

# ifdef USE_CUDA
  if(gpu_driver.eq.0) then
# endif

   ! OpenACC CODE ----------------------------

#ifdef USE_OPEN_ACC
     if(i.ne.iold) then
       if(smp_me.eq.comm_core) then
          ilo = nu*(i-1) + 1
          ihi = ilo + nu - 1
          call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_i)
!$acc update device(tmp_i(1:nutr*nu)) async(1)
          call ddcc_t_getve_acc(1,nu,i,tmp_i,ve_i)
          call trant3_1_async(1,nu,ve_i)
       end if
       comm_core = comm_core+1
       if(comm_core.eq.smp_np) comm_core=0
     end if

     if(j.ne.jold) then
       if(smp_me.eq.comm_core) then
          ilo = nu*(j-1) + 1
          ihi = ilo + nu - 1
          call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_j)
!$acc update device(tmp_j(1:nutr*nu)) async(2)
          call ddcc_t_getve_acc(2,nu,j,tmp_j,ve_j)
          call trant3_1_async(2,nu,ve_j)
       end if
       comm_core = comm_core+1
       if(comm_core.eq.smp_np) comm_core=0
     end if

     if(k.ne.kold) then
       if(smp_me.eq.comm_core) then
          ilo = nu*(k-1) + 1
          ihi = ilo + nu - 1
          call ddi_get(d_vvvo,1,nutr,ilo,ihi,tmp_k)
!$acc update device(tmp_k(1:nutr*nu)) async(3)
          call ddcc_t_getve_acc(3,nu,k,tmp_k,ve_k)
          call trant3_1_async(3,nu,ve_k)
       end if
       comm_core = comm_core+1
       if(comm_core.eq.smp_np) comm_core=0
     end if

     if(smp_np.gt.1) call smp_sync()
#endif

# ifdef USE_CUDA
  else ! gpu_driver.eq.1

   ! CUDA CODE ------------------------------
     if(i.ne.iold) then
!      if(smp_me.eq.comm_core) call ddcc_t_getve(nu,i,tmp,vei)
       ilo = nu*(i-1) + 1
       ihi = ilo + nu       ! should one be subtracted?
       if(smp_me.eq.comm_core) call ddi_get(d_vvvo,1,nutr,ilo,ihi,vei)
       comm_core = comm_core+1
       if(comm_core.eq.smp_np) comm_core=0
     end if
!   
     if(j.ne.jold) then
       ilo = nu*(j-1) + 1
       ihi = ilo + nu       ! should one be subtracted?
       if(smp_me.eq.comm_core) call ddi_get(d_vvvo,1,nutr,ilo,ihi,vej)
       if(comm_core.eq.smp_np) comm_core=0
     end if
   
     if(k.ne.kold) then
       ilo = nu*(k-1) + 1
       ihi = ilo + nu       ! should one be subtracted?
       if(smp_me.eq.comm_core) call ddi_get(d_vvvo,1,nutr,ilo,ihi,ve_k)
       comm_core = comm_core+1
       if(comm_core.eq.smp_np) comm_core=0
     end if
     ! NOTE: ddcc_t_ijk_gpu must expand and perfrom a trant3_1 operation
     ! on vei, vej, vek prior to use!

  end if ! gpu_driver.eq.1
# endif

# ifdef USE_CUDA
  if(gpu_driver.eq.0) then
# endif

     call ddcc_t_ijk_acc(no,nu,i,j,k,v1,t2,vm,v3,voe,eh,ep,ve_i,ve_j,ve_k)

# ifdef USE_CUDA
  else
     call ddcc_t_ijk_gpu(no,nu,i,j,k,v1,t2,vm,voe,eh,ep,vei,vej,ve_k) !  call ijk_gpu_driver
  end if
# endif
  if(smp_np.gt.1) call smp_sync()

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


#ifdef USE_CUDA
call ijk_gpu_finalize(no,nu,etd)   ! deallocate device-arrays, copy back "v1, etd"
deallocate(vei,vej)
end if ! gpu_driver == 1
#endif

! call ddi_sync(1234)

! remote sync at this point in hybrid code
! this was used to measure the ijk tuple time
!
! call ddi_sync(1234)
! if(ddi_me.eq.0) then
!    ijk_stop = mpi_wtime()
!    print *,"IJK-tuples time=",(ijk_stop-ijk_start),"  etd=",etd
! end if


#ifdef USE_CUDA
!if(gpu_driver .eq. 0) then    ! <--- vjg commented out so iij/ijj tuples would run
#endif

! counters and load-balancing for iij and ijj tuples
icntr = 0

if(smp_me.eq.0) call ddi_dlbnext(mytask)
call ddi_smp_bcast(1237,'I',mytask,1,0)

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
       call ddi_smp_bcast(1235,'I',mytask,1,0)
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
       call ddi_smp_bcast(1234,'I',mytask,1,0)
    end if
    icntr = icntr + 1
  end do
end do
! ----------- end of iij and ijj tuples -------------

!$acc end data       

! if(smp_me.eq.0) write(6,*) 'cpu node ',ddi_my,' finshed iij/ijj'
#ifdef USE_CUDA
!end if  ! gpu_driver == 0   ! vjg commented out 
#endif

! switch scopes
#ifdef USE_CUDA
call ddi_sync(1234)
working_smp_comm = global_smp_comm
working_compute_comm = global_compute_comm
call mpi_comm_rank(working_smp_comm, smp_me, ierr)
call mpi_comm_size(working_smp_comm, smp_np, ierr)
call mpi_comm_rank(working_compute_comm, ddi_me, ierr)
call mpi_comm_size(working_compute_comm, ddi_np, ierr)
call ddi_sync(1234)
#endif


call ddi_gsumf(123,eh,no)
call ddi_gsumf(124,ep,nu)
call ddi_gsumf(125,v1,nou)
call ddi_gsumf(126,etd,1)

call ddi_sync(1)

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
end subroutine cc_triples_acc
