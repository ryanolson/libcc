      SUBROUTINE DDCC_T_IJK_GPU(NO,NU,I,J,K,T1,T2,VM,VOE,
     *                          EH,EP,ve_i,ve_j,ve_k)
      use common_cc, only: smp_np, smp_me, nu2, etd
      IMPLICIT NONE
C
      INTEGER NO,NU,I,J,K
      DOUBLE PRECISION T1(NU,*),VM(NO,NU,NO,NO),
     &                 VOE(NU,NU,NO,NO),
     &                 EH(NO),EP(NU),T2(NU*NU,NO,NO),
     &                 ve_i(*),ve_j(*),ve_k(*) 
C
      call ddcc_t_ijk_big_cuda_wrapper(nu,no,i,j,k,
     &               t2(1,1,i),t2(1,1,j),t2(1,1,k),
     &               vm(1,1,i,j),vm(1,1,j,i),vm(1,1,i,k),
     &               vm(1,1,k,i),vm(1,1,j,k),vm(1,1,k,j),
     &               ve_i,ve_j,ve_k,
     &               voe(1,1,i,j),voe(1,1,j,i),voe(1,1,i,k),
     &               voe(1,1,k,i),voe(1,1,j,k),voe(1,1,k,j),
     &               t1, eh, ep, etd )
      RETURN
      END

