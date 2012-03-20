program gpu_driver
implicit none

integer :: world_comm, gpu_comm
integer :: ddi_np, ddi_me
integer :: ddi_nn, ddi_my
integer :: smp_np, smp_me, gpu_count
integer :: mpi_me

call ddi_init
call ddi_memory(0,1000,0)

call ddi_get_workingcomm( world_comm )

call ddi_nproc(ddi_np,ddi_me)
call ddi_nnode(ddi_nn,ddi_my)
call ddi_smp_nproc(smp_np,smp_me)
write(6,1000) ddi_np,ddi_me,ddi_nn,ddi_my,smp_np,smp_me

mpi_me = ddi_me

call ddi_gpu_createcomm( world_comm, gpu_comm )
call ddi_scope( gpu_comm )

call ddi_nproc(ddi_np,ddi_me)
call ddi_nnode(ddi_nn,ddi_my)
call ddi_smp_nproc(smp_np,smp_me)
call ddi_gpu_device_count(gpu_count)
write(6,1001) mpi_me,ddi_np,ddi_me,ddi_nn,ddi_my,smp_np,smp_me,gpu_count



if(smp_me.eq.0) then
   write(6,*) 'testing smp sync',mpi_me
endif
call ddi_smp_sync()
if(smp_me.eq.0) then
   write(6,*) 'smp sync completed',mpi_me
endif

call ddi_scope( world_comm )
call ddi_finalize()

1000 format(6I6)
1001 format(8I6)
end program gpu_driver
