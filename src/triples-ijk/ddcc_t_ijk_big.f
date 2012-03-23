      SUBROUTINE DDCC_T_IJK_BIG(NO,NU,I,J,K,T1,T2,VM,V3,VOE,
     *                          EH,EP,ve_i,ve_j,ve_k)
      use common_cc, only: smp_np, smp_me, nu2
      IMPLICIT NONE
C
      INTEGER NO,NU,I,J,K
C
      DOUBLE PRECISION T1(*),VM(NO,NU,NO,NO),
     &                 V3(*),VOE(NU,NU,NO,NO),
     &                 EH(NO),EP(NU),T2(NU*NU,NO,NO),
     &                 ve_i(*),ve_j(*),ve_k(*) 
C
      call ijk_tuple(i,j,k,t2(1,1,i),t2(1,1,j),t2(1,1,k),
     &               vm(1,1,i,j),vm(1,1,j,i),vm(1,1,i,k),
     &               vm(1,1,k,i),vm(1,1,j,k),vm(1,1,k,j),
     &               ve_i,ve_j,ve_k,v3)
      call smp_sync()
      call t1wt3_ijk(i,j,k,no,nu,v3,
     &              voe(1,1,i,j),voe(1,1,j,i),voe(1,1,i,k),voe(1,1,k,i),
     &              voe(1,1,j,k),voe(1,1,k,j),t1,eh,ep)
      RETURN
      END

