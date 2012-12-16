
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

      SUBROUTINE CPYTSQ(A,B,NA,INCA)
      use common_cc
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

      END

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


      SUBROUTINE TRANMD_SMP(A,N1,N2,N3,N4,IJ)
      use common_cc, only: smp_np, smp_me
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ::  A(*)
     
      if(smp_np.gt.1) CALL smp_sync()
 
        N12=N1*N2
        N123=N12*N3
        IF(IJ.EQ.12) GO TO 12
        IF(IJ.EQ.13) GO TO 13
        IF(IJ.EQ.14) GO TO 14
        IF(IJ.EQ.23) GO TO 23
        IF(IJ.EQ.24) GO TO 24
        IF(IJ.EQ.34) GO TO 34
        IF(IJ.EQ.231) GO TO 231
        IF(IJ.EQ.312) GO TO 312
        IF(IJ.EQ.341) GO TO 341
        IF(IJ.EQ.413) GO TO 413
        IF(IJ.EQ.1234) GO TO 1234
        GOTO 100
   12 CONTINUE
        DO 11 L=1,N4
           N123L=N123*(L-1)
        DO 9 K=1,N3
           IF(MOD(K,SMP_NP).NE.SMP_ME) GOTO 9
           N12K=(K-1)*N12+N123L
        DO 10 I=1,N1
           N1I=(I-1)*N1
        DO 10 J=1,I
           IJKL=N12K+(J-1)*N1+I
           JIKL=N12K+N1I+J
        X=A(IJKL)
        A(IJKL)=A(JIKL)
        A(JIKL)=X
   10 CONTINUE
    9 CONTINUE
   11 CONTINUE
        GO TO 100
   13 CONTINUE
        DO 22 L=1,N4
           N123L=N123*(L-1)
        DO 21 I=1,N1
           IF(MOD(I,SMP_NP).NE.SMP_ME) GOTO 21
           N12I=N123L+N12*(I-1)
        DO 20 J=1,N2
           N1J=N1*(J-1)
        DO 20 K=1,I
           IJKL=N123L+(K-1)*N12+N1J+I
           KJIL=N12I+N1J+K
        X=A(IJKL)
        A(IJKL)=A(KJIL)
        A(KJIL)=X
   20 CONTINUE
   21 CONTINUE
   22 CONTINUE
        GO TO 100
   14 CONTINUE
        DO 26 I=1,N1
           IF(MOD(I,SMP_NP).NE.SMP_ME) GOTO 26
           N123I=N123*(I-1)
        DO 25 J=1,N2
           N1J=N1*(J-1)
        DO 25 K=1,N3
           N12K=N12*(K-1)+N1J
        DO 25 L=1,I
           IJKL=(L-1)*N123+N12K+I
           LJKI=N123I     +N12K+L
        X=A(IJKL)
        A(IJKL)=A(LJKI)
        A(LJKI)= X
   25 CONTINUE
   26 CONTINUE
        GO TO 100
   23 CONTINUE
        DO 31 J=1,N2
           IF(MOD(J,SMP_NP).NE.SMP_ME) GOTO 31
           J1=J-1
           N1J=N1*J1
           N12J=N12*J1
        DO 30 K=1,J
           K1=K-1
           N1K=N1*K1
           N12K=N12*K1
        DO 30 L=1,N4
           N123L=N123*(L-1)
        DO 30 I=1,N1
           N123LI=N123L+I
           IJKL=N123LI+N12K+N1J
           IKJL=N123LI+N12J+N1K
        X=A(IJKL)
        A(IJKL)=A(IKJL)
        A(IKJL)=X
   30 CONTINUE
   31 CONTINUE
        GO TO 100
   24 CONTINUE
        DO 41 J=1,N2
           IF(MOD(J,SMP_NP).NE.SMP_ME) GOTO 41
           J1=J-1
           N123J=N123*J1
           N1J  =N1*J1
        DO 40 L=1,J
           L1=L-1
           N123L=N123*L1
           N1L=N1*L1
        DO 40 K=1,N3
           N12K=N12*(K-1)
        DO 40 I=1,N1
           N12KI=N12K+I
           IJKL=N123L+N12KI+N1J
           ILKJ=N123J+N12KI+N1L
        X=A(IJKL)
        A(IJKL)=A(ILKJ)
        A(ILKJ)=X
   40 CONTINUE
   41 CONTINUE
        GO TO 100
   34 CONTINUE
        DO 51 K=1,N3
           IF(MOD(K,SMP_NP).NE.SMP_ME) GOTO 51
           K1=K-1
           N12K=N12*K1
           N123K=N123*K1
        DO 50 L=1,K
           L1=L-1
           N12L=N12*L1
           N123L=N123*L1
        DO 50 J=1,N2
           J1=J-1
           N1J=N1*(J-1)
        DO 50 I=1,N1
           N1JI=N1J+I
           IJKL=N123L+N12K+N1JI
           IJLK=N123K+N12L+N1JI
        X=A(IJKL)
        A(IJKL)=A(IJLK)
        A(IJLK)=X
   50 CONTINUE
   51 CONTINUE
        GO TO 100
 231  CONTINUE
        DO 61 L=1,N4
           IF(MOD(L,SMP_NP).NE.SMP_ME) GOTO 61
           N123L=N123*(L-1)
        DO 60 J=1,N1
           J1=J-1
           N12J=N12*J1
           N1J=N1*J1
        DO 60 K=1,J
           N12JK=N12J+K
           K1=K-1
           N12K=N12*K1
           N1K=N1*K1
        DO 60 I=1,K
           I1=I-1
           IJKL=N123L+N12K+N1J+I
           JKIL=N123L+I1*N12+N1K+J
           KIJL=N123L+N12JK +I1*N1
        X=A(IJKL)
        A(IJKL)=A(JKIL)
        A(JKIL)=A(KIJL)
        A(KIJL)=X
        IF(J.EQ.K.OR.K.EQ.I) GOTO 60
           JIKL=N123L+N12K  +N1*I1+J
           IKJL=N123L+N12J  +N1K  +I
           KJIL=N123L+N12*I1+N1J  +K
        X=A(JIKL)
        A(JIKL)=A(IKJL)
        A(IKJL)=A(KJIL)
        A(KJIL)=X
 60   CONTINUE
 61   CONTINUE
        GOTO 100
 312  CONTINUE
        DO 71 L=1,N4
           IF(MOD(L,SMP_NP).NE.SMP_ME) GOTO 71
           N123L=N123*(L-1)
        DO 70 I=1,N1
           I1=I-1
           N1I=N1*I1
           N12I=N12*I1
        DO 70 J=1,I
           J1=J-1
           N1J =N1*J1
           N12J=N12*J1
        DO 70 K=1,J
           K1=K-1
           N1K=N1*K1
           N12K=N12*K1
           IJKL=N123L+N12K+N1J+I
           JKIL=N123L+N12I+N1K+J
           KIJL=N123L+N12J+N1I+K
        X=A(JKIL)
        A(JKIL)=A(IJKL)
        A(IJKL)=A(KIJL)
        A(KIJL)=X
        IF (I.EQ.J.OR.J.EQ.K) GOTO 70
           IKJL=N123L+N12J+N1K+I
           JIKL=N123L+N12K+N1I+J
           KJIL=N123L+N12I+N1J+K
        X=A(IKJL)
        A(IKJL)=A(JIKL)
        A(JIKL)=A(KJIL)
        A(KJIL)=X
 70   CONTINUE
 71   CONTINUE
        GOTO 100
 341  CONTINUE
        DO 81 L=1,N2
           IF(MOD(L,SMP_NP).NE.SMP_ME) GOTO 81
        DO 80 J=1,N1
        DO 80 K=1,J
        DO 80 I=1,K
           ILJK=(K-1)*N123+(J-1)*N12+(L-1)*N1+I
           JLKI=(I-1)*N123+(K-1)*N12+(L-1)*N1+J
           KLIJ=(J-1)*N123+(I-1)*N12+(L-1)*N1+K
        X=A(ILJK)
        A(ILJK)=A(JLKI)
        A(JLKI)=A(KLIJ)
        A(KLIJ)=X
        IF(J.EQ.K.OR.K.EQ.I) GOTO 80
           ILKJ=(J-1)*N123+(K-1)*N12+(L-1)*N1+I
           JLIK=(K-1)*N123+(I-1)*N12+(L-1)*N1+J
           KLJI=(I-1)*N123+(J-1)*N12+(L-1)*N1+K
        X=A(JLIK)
        A(JLIK)=A(ILKJ)
        A(ILKJ)=A(KLJI)
        A(KLJI)=X
 80   CONTINUE
 81   CONTINUE
        GOTO 100
 413  CONTINUE
        DO 91 L=1,N2
           IF(MOD(L,SMP_NP).NE.SMP_ME) GOTO 91
        DO 90 I=1,N1
        DO 90 J=1,I
        DO 90 K=1,J
           JLKI=(I-1)*N123+(K-1)*N12+(L-1)*N1+J
           ILJK=(K-1)*N123+(J-1)*N12+(L-1)*N1+I
           KLIJ=(J-1)*N123+(I-1)*N12+(L-1)*N1+K
        X=A(JLKI)
        A(JLKI)=A(ILJK)
        A(ILJK)=A(KLIJ)
        A(KLIJ)=X
        IF (I.EQ.J.OR.J.EQ.K) GOTO 90
           ILKJ=(J-1)*N123+(K-1)*N12+(L-1)*N1+I
           JLIK=(K-1)*N123+(I-1)*N12+(L-1)*N1+J
           KLJI=(I-1)*N123+(J-1)*N12+(L-1)*N1+K
        X=A(ILKJ)
        A(ILKJ)=A(JLIK)
        A(JLIK)=A(KLJI)
        A(KLJI)=X
 90   CONTINUE
 91   CONTINUE
        GO TO 100
 1234 CONTINUE
!--      WRITE(6,76) A
!-- 76   FORMAT(4F15.10)
        DO 96 I=1,N1
           IF(MOD(I,SMP_NP).NE.SMP_ME) GOTO 96
        DO 95 J=1,N2
        DO 95 K=1,J
        DO 95 L=1,I
           IJKL=(L-1)*N123+(K-1)*N12+(J-1)*N1+I
           LKJI=(I-1)*N123+(J-1)*N12+(K-1)*N1+L
        X=A(IJKL)
        A(IJKL)=A(LKJI)
        A(LKJI)=X
        IF (I.EQ.L.OR.K.EQ.J) GOTO 95
           LJKI=(I-1)*N123+(K-1)*N12+(J-1)*N1+L
           IKJL=(L-1)*N123+(J-1)*N12+(K-1)*N1+I
        X=A(LJKI)
        A(LJKI)=A(IKJL)
        A(IKJL)=X
 95   CONTINUE
 96   CONTINUE
 100  CONTINUE
        if(smp_np.gt.1) CALL smp_sync()
        RETURN
        END SUBROUTINE TRANMD_SMP

