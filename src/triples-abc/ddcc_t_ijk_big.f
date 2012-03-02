      SUBROUTINE DDCC_T_IJK_BIG(NO,NU,I,J,K,T1,T2,VM,V3,T3,VOE,O1,TMP,
     *                          EH,EP,ve_i,ve_j,ve_k)
      IMPLICIT NONE
C
      INTEGER NO,NU,I,J,K
C
      DOUBLE PRECISION T1(*),VM(NO,NU,NO,NO),
     *                 V3(*),T3(*),TMP(NU*NU),
     *                 VOE(1),O1(*),EH(NO),EP(NU),T2(NU*NU,NO,NO),
     &                 ve_i(*),ve_j(*),ve_k(*) 
C
      DOUBLE PRECISION DEH
      INTEGER ITMP,JTMP,KTMP
C
      DOUBLE PRECISION ZERO,ONE,OM
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,OM=-1.0D+00)
C
      call ijk_tuple_formv3(i,j,k,t2(1,1,i),t2(1,1,j),t2(1,1,k),
     &                      vm(1,1,i,j),vm(1,1,j,i),vm(1,1,i,k),
     &                      vm(1,1,k,i),vm(1,1,j,k),vm(1,1,k,j),
     &                      ve_i,ve_j,ve_k,v3)
C
      CALL T3SQUA_SMP(I,J,K,NO,NU,O1,T2,V3,EH,EP)
      DEH=EH(I)+EH(J)+EH(K)
      CALL ADT3DEN_SMP(NU,DEH,V3,EP)
C
      ITMP=I
      JTMP=J
      KTMP=K
      CALL DRT1WT3IJK_SMP(ITMP,JTMP,KTMP,NO,NU,T1,VOE,V3,T3)
      RETURN
      END

