      SUBROUTINE DDCC_T_IJK_BIG(NO,NU,I,J,K,T1,T2,VM,V3,VOE,
     *                          EH,EP,ve_i,ve_j,ve_k)
      use common_cc, only: smp_np, smp_me, nu2, nu3
      IMPLICIT NONE
C
      INTEGER NO,NU,I,J,K
C
      DOUBLE PRECISION T1(*),VM(NO,NU,NO,NO),
     &                 V3(NU3),VOE(NU,NU,NO,NO),
     &                 EH(NO),EP(NU),T2(NU*NU,NO,NO),
     &                 ve_i(*),ve_j(*),ve_k(*) 
C
      call ijk_tuple(i,j,k,T2(1,1,i),T2(1,1,j),T2(1,1,k),
     &               VM(1,1,i,j),VM(1,1,j,i),VM(1,1,i,k),
     &               VM(1,1,k,i),VM(1,1,j,k),VM(1,1,k,j),
     &               ve_i,ve_j,ve_k,V3)
      if(smp_np.gt.1) call smp_sync()
      call t1wt3_ijk(i,j,k,no,nu,V3,
     &              VOE(1,1,i,j),VOE(1,1,j,i),VOE(1,1,i,k),VOE(1,1,k,i),
     &              VOE(1,1,j,k),VOE(1,1,k,j),T1,EH,EP)
      RETURN
      END

