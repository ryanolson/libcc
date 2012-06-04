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
#else
    ! sync with ranks compiled with cuda
    call ddi_gsumi(ngpus)
#endif

if(ngpus.gt.0) then
   sr = 0 ! adjust ijk distribution
   nr = 0 ! for now, we assign all ijk work to the gpu
endif

call triples_cpu_driver(no, nu, sr, nr)
end subroutine triples_calc()


subroutine ijk_task(ijk,no,i,j,k)
implicit none
integer, intent(in)  :: ijk, no
integer, intent(out) :: i, j, k

integer :: icntr
integer :: II,JJ,KK,I1,J1,ICNTR
      ICNTR=0
      DO II=1,NO
         I1=II-1
         DO JJ=1,I1
            J1=JJ-1
            DO KK=1,J1
               IF(ICNTR.EQ.MYTASK) THEN
                  I=II
                  J=JJ
                  K=KK
                  GO TO 10
               END IF
               ICNTR=ICNTR+1
            END DO
         END DO
      END DO
   10 CONTINUE
      RETURN
end subroutine ijk_task


end module cc_triples
