      SUBROUTINE TRANT3_SMP(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU)
      INTEGER SMP_ME,SMP_NP
C
      CALL DDI_SMP_NPROC(SMP_NP,SMP_ME)
      CALL DDI_SMP_SYNC()
C
      GO TO (1,2,3,4,5),ID
    1 CONTINUE
      DO 101 B=1,NU
        IF(MOD(B,SMP_NP).NE.SMP_ME) GOTO 101
      DO 100 C=1,B
      DO 100 A=1,NU
      X=V(A,B,C)
      V(A,B,C)=V(A,C,B)
      V(A,C,B)=X
  100 CONTINUE
  101 CONTINUE
      GO TO 1000
    2 CONTINUE
      DO 201 C=1,NU
        IF(MOD(C,SMP_NP).NE.SMP_ME) GOTO 201
      DO 200 A=1,NU
      DO 200 B=1,A
      X=V(A,B,C)
      V(A,B,C)=V(B,A,C)
      V(B,A,C)=X
  200 CONTINUE
  201 CONTINUE
      GO TO 1000
    3 CONTINUE
      DO 301 B=1,NU
        IF(MOD(B,SMP_NP).NE.SMP_ME) GOTO 301
      DO 300 A=1,NU
      DO 300 C=1,A
      X=V(A,B,C)
      V(A,B,C)=V(C,B,A)
      V(C,B,A)=X
  300 CONTINUE
  301 CONTINUE
      GO TO 1000
    4 CONTINUE
      DO 401 B=1,NU
        IF(MOD(B,SMP_NP).NE.SMP_ME) GOTO 401
      DO 400 C=1,B
      DO 400 A=1,C
      X=V(A,B,C)
      V(A,B,C)=V(B,C,A)
      V(B,C,A)=V(C,A,B)
      V(C,A,B)=X
      IF(B.EQ.C.OR.C.EQ.A) GO TO 400
      X=V(B,A,C)
      V(B,A,C)=V(A,C,B)
      V(A,C,B)=V(C,B,A)
      V(C,B,A)=X
  400 CONTINUE
  401 CONTINUE
      GO TO 1000
    5 CONTINUE
      DO 501 A=1,NU
        IF(MOD(A,SMP_NP).NE.SMP_ME) GOTO 501
      DO 500 C=1,A
      DO 500 D=1,C
      X=V(C,D,A)
      V(C,D,A)=V(A,C,D)
      V(A,C,D)=V(D,A,C)
      V(D,A,C)=X
      IF (A.EQ.C.OR.C.EQ.D) GO TO 500
      X=V(A,D,C)
      V(A,D,C)=V(C,A,D)
      V(C,A,D)=V(D,C,A)
      V(D,C,A)=X
  500 CONTINUE
  501 CONTINUE
      GO TO 1000
 1000 CONTINUE
      CALL DDI_SMP_SYNC()
      RETURN
      END
