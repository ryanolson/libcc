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
      call ijk_gpu_driver(nu,no,i,j,k,
     &               ve_i,ve_j,ve_k,
     &               t1, eh, ep, etd )
      RETURN
      END

