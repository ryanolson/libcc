C =============================================================================
C   Copyright (C) 2010.  Ryan M. Olson
C
C   This program is free software: you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation, either version 3 of the License, or
C   (at your option) any later version.

C   This program is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.

C   You should have received a copy of the GNU General Public License
C   along with this program.  If not, see <http://www.gnu.org/licenses/>.
C =============================================================================
C
      SUBROUTINE DDCC_T_TASK(MYTASK,NO,I,J,K)
      IMPLICIT NONE
C
      INTEGER MYTASK,NO,I,J,K
      INTEGER II,JJ,KK,I1,J1,ICNTR
C
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
      END

      SUBROUTINE DIV_EVEN(N,NP,ME,NR,SR)
      IMPLICIT NONE
      INTEGER N,NP,ME,NR,SR
      INTEGER NE
C
      NR = N / NP
      NE = MOD(N,NP)
C
      IF(ME.LT.NE) THEN
         NR = NR + 1
         SR = NR*ME + 1
      ELSE
         SR = (NR+1)*NE + NR*(ME-NE) + 1
      END IF
C
      RETURN
      END

      SUBROUTINE DDCC_T_GETVE(NUINP,INDX,TEMP,VE)
      use common_cc
      INTEGER INDX,I,IOFF
C
C     TEMP SHOULD BE NU^2 AND VE SHOULD BE NU^3
      DOUBLE PRECISION TEMP(NU*NU),VE(NU,NU,NU)
      INTEGER LOOP
C
      IOFF = 1
      LOOP = INDX - 1
C
      DO I = LOOP*NU, (LOOP*NU+NU-1)
        CALL DDI_GET(D_VVVO,1,NUTR,I+1,I+1,TEMP)
        CALL CPYTSQ(TEMP,VE(1,1,IOFF),NU,1)
        IOFF = IOFF + 1
      END DO
C
      RETURN
      END

      SUBROUTINE CPYTSQ(A,B,NA,INCA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),B(NA,NA)
C
C     ---- COPY TRIANGULAR A TO SQUARE B (NA BY NA) ----
C     THE INCREMENT BETWEEN ELEMENTS OF A WILL USUALLY BE 1.
C
      IJ=1
      DO 200 I=1,NA
         DO 100 J=1,I
            B(I,J) = A(IJ)
            B(J,I) = A(IJ)
            IJ = IJ + INCA
  100    CONTINUE
  200 CONTINUE
      RETURN
      END


      SUBROUTINE ZEROT3(T3,NU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A
      DIMENSION T3(NU,NU,NU)
      DATA ZERO/0.0D+00/
C
      DO  A=1,NU
      T3(A,A,A)=ZERO
      ENDDO           
      RETURN 
      END 


      SUBROUTINE ADT3DEN_SMP(NU,DEH,T3,EP)
      use common_cc, only: smp_np, smp_me
      IMPLICIT NONE
      INTEGER A,B,C,NU
      DOUBLE PRECISION T3(NU,NU,NU),EP(NU),DEH,DEN
C
      CALL smp_sync()
      DO 11 A=1,NU
        IF(MOD(A,SMP_NP).NE.SMP_ME) GOTO 11
      DO 10 B=1,NU
      DO 10 C=1,NU
      DEN=DEH-EP(A)-EP(B)-EP(C)
      T3(C,B,A)=T3(C,B,A)/DEN
 10   CONTINUE
 11   CONTINUE
      CALL smp_sync()
      RETURN
      END

      SUBROUTINE DRT1WT3IJK_SMP(I,J,K,NO,NU,T1,VOE,TI,T3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T1(1),VOE(1),T3(1),TI(1)
C
      CALL T1WT3IJK_SMP(I,J,K,NO,NU,T1,VOE,TI,T3)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_SMP(I,K,J,NO,NU,T1,VOE,TI,T3)
      call trant3_4(nu,ti)
      CALL T1WT3IJK_SMP(J,I,K,NO,NU,T1,VOE,TI,T3)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_SMP(J,K,I,NO,NU,T1,VOE,TI,T3)
      call trant3_5(nu,ti)
      CALL T1WT3IJK_SMP(K,I,J,NO,NU,T1,VOE,TI,T3)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_SMP(K,J,I,NO,NU,T1,VOE,TI,T3)
      RETURN
      END

      SUBROUTINE T1WT3IJK_SMP(I,J,K,NO,NU,T1,VOE,TI,T3)
      use common_cc, only: smp_np, smp_me
      IMPLICIT NONE
      INTEGER NR,SR,T3OFF,VOEOFF,I,J,K,NO,NU,NU2
      DOUBLE PRECISION T1(NU,NO),VOE(*),T3(1),TI(1),ONE
      DATA ONE/1.0D+00/
C
      NU2 = NU*NU
      CALL SMT3FOUR_SMP(NU,T3,TI)
      CALL DIV_EVEN(NU2,SMP_NP,SMP_ME,NR,SR)
      T3OFF = (SR - 1)*NU + 1
      VOEOFF = (K - 1)*NU2*NO + (J - 1)*NU2  + SR
      CALL DGEMM('N','N',NU,1,NR,ONE,T3(T3OFF),NU,VOE(VOEOFF),NU2,ONE,
     *           T1(1,I),NU)
      CALL smp_sync()
      RETURN
      END

      SUBROUTINE SMT3FOUR_SMP(NU,T3,V3)
      use common_cc, only: smp_np, smp_me
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION T3(NU,NU,NU),V3(NU,NU,NU)
      DATA TWO/2.0D+00/
      DATA ZERO/0.0D+00/
      DATA OM/-1.0D+00/
C
      CALL smp_sync()
C
      DO 2 C=1,NU
        IF(MOD(C,SMP_NP).NE.SMP_ME) GOTO 2

      do b=1,nu
         t3(b,b,c) = zero
      end do

      DO B=2,NU
      DO A=1,B-1
         T3(A,B,C)=(V3(A,B,C)-V3(B,A,C))*TWO-V3(A,C,B)+V3(B,C,A)
         T3(B,A,C)=T3(A,B,C)*OM
      END DO
      END DO

 2    CONTINUE
      CALL smp_sync()
      RETURN
      END

C*MODULE CCDDI   *DECK DRT1WT3IJ_SMP
      SUBROUTINE DRT1WT3IJ_SMP(I,J,NO,NU,T1,VOE,TI,T3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T1(1),VOE(1),T3(1),TI(1)
C
      CALL T1WT3IJK_SMP(I,J,J,NO,NU,T1,VOE,TI,T3)
      CALL TRANT3_SMP(TI,NU,2)
      CALL T1WT3IJK_SMP(J,I,J,NO,NU,T1,VOE,TI,T3)
c     CALL TRANT3_SMP(TI,NU,1)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_SMP(J,J,I,NO,NU,T1,VOE,TI,T3)
      RETURN
      END

      SUBROUTINE DRT1WT3JK_SMP(J,K,NO,NU,T1,VOE,TI,T3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T1(1),VOE(1),T3(1),TI(1)
C
      CALL T1WT3IJK_SMP(J,J,K,NO,NU,T1,VOE,TI,T3)
C     CALL TRANT3_SMP(TI,NU,1)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_SMP(J,K,J,NO,NU,T1,VOE,TI,T3)
      CALL TRANT3_SMP(TI,NU,2)
      CALL T1WT3IJK_SMP(K,J,J,NO,NU,T1,VOE,TI,T3)
      RETURN
      END


      SUBROUTINE SYMT311(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION V(NU,NU,NU)
C     
      IF (ID.EQ.23) GO TO 23
      IF (ID.EQ.12) GO TO 12
      IF (ID.EQ.13) GO TO 13
 23   CONTINUE
      DO 100 B=1,NU
      DO 100 C=1,B
      DO 100 A=1,NU
      X=V(A,B,C)+V(A,C,B)
      V(A,B,C)=X
      V(A,C,B)=X
  100 CONTINUE
      GO TO 1000
 12   CONTINUE
      DO 101 C=1,NU
      DO 101 A=1,NU
      DO 101 B=1,A
      X=V(A,B,C)+V(B,A,C)  
      V(A,B,C)=X
      V(B,A,C)=X
 101  CONTINUE
      GO TO 1000
 13   CONTINUE
      DO 102 A=1,NU
      DO 102 B=1,NU
      DO 102 C=1,A
      X=V(A,B,C)+V(C,B,A)
      V(A,B,C)=X
      V(C,B,A)=X
 102  CONTINUE
      GO TO 1000
 1000 CONTINUE
      END

      SUBROUTINE DD_T3SQUA_GSUM
      use common_cc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     COMMON /CCRENO/ OSS,ODS,ODD,OTS,OTD,OTT,ODS_S,ODS_D,ODS_T,
c    *                OQS,OQDS,OQDD,OQTS,ESD,ETD,ETS,ETTM,ESD_TM
C
c     CALL DDI_GSUMF(1234,OTS,1)
c     CALL DDI_GSUMF(1234,OTD,1)
      CALL DDI_GSUMF(1234,ETD,1)
      RETURN
      end


C*MODULE MTHLIB  *DECK TRPOSE
      SUBROUTINE TRPOSE(A,B,N,M,KIND)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION A(N,M), B(M,N)
C
C* 14 JAN 1983 - STE * 8 MAR 1980
C*
C*    AUTHOR: S. T. ELBERT (AMES LABORATORY-USDOE)
C*
C*    PURPOSE -
C*       STORE TRANSPOSE OF N BY M MATRIX A IN MATRIX B OR A
C*             **   ****
C*
C*    ON ENTRY -
C*       A     - W.P. REAL (N,M)
C*               MATRIX TO BE TRANSPOSED
C*       N      - INTEGER
C*                ROWS OF INPUT MATRIX, COLUMNS OF OUTPUT MATRIX
C*       M      - INTEGER
C*                COLUMNS OF INPUT MATRIX, ROWS OF OUTPUT MATRIX
C*       KIND   - INTEGER
C*                IF NOT ZERO, TRANSPOSED MATRIX IS COPYIED BACK INTO A
C*
C*    ON EXIT -
C*       B      - W.P. REAL (M,N)
C*                TRANSPOSED COPY OF INPUT MATRIX
C*       A (OPTIONAL) - W.P. REAL (M,N)
C*                TRANSPOSED COPY OF INPUT MATRIX
C
      IF(N.LE.0 .OR. M.LE.0) RETURN
      DO 120 J=1,M
         DO 110 I=1,N
            B(J,I) = A(I,J)
  110    CONTINUE
  120 CONTINUE
      IF(KIND.NE.0) CALL DCOPY(M*N,B,1,A,1)
      RETURN
      END

