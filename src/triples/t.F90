module cc_triples
implicit none



contains

subroutine triples_init()
end subroutine triples_init()


subroutine triples_finalize()
end subroutine triples_finalize()


subroutine triples_calc()
implicit none
integer :: ntuples
integer :: cuda_driver, ngpus

ngpus = 0
cuda_driver = 0
ntuples = no*(no-1)*(no-2)/6
call div_even(ntuples, ddi_nn, ddi_my, nr, sr)

#if HAVE_CUDA
    if(smp_me .lt. cuda_device_count()) cuda_driver = 1
    call split(acc_driver, ddi_me, sg)
    call nprocs(nn, my, sg)
   
    ! sync with ranks not compiled with cuda
    ngpus = cuda_driver
    call ddi_gsumi(ngpus)

    ! divide up ijk work between gpus
    call div_even(ntuples, nn, my, nr, sr)

    if(cuda_driver .eq. 1) then
       call triples_cuda_init( ... )
       call triples_cuda_driver(no, nu, sr, nr, d_vvvo)
       call triples_cuda_finalize( ... )
       return
    end if
#endif

! sync with ranks compiled with cuda
call ddi_gsumi(ngpus)
if(ngpus.gt.0) then
   sr = 0 ! adjust ijk distribution
   nr = 0 ! for now, we assign all ijk work to the gpu
endif

call triples_driver(no, nu, sr, nr)
end subroutine triples_calc()


end module cc_triples
