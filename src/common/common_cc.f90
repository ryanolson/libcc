module common_cc
implicit none

integer :: no, nu
integer :: nou, no2, no3, no4, nu2, nu3, nu4, no2u, no3u, nou2, nou3, no2u2
integer :: notr, nutr ! number of elements in triagular matrices of order no and nu

! MO integrals stored in d_vvvo are of [p4 p5 | p6 h*] with permutational symmetry
! only in the bra [ p4 p5 | = [ p5 p4 | ... therefore the # of rows is nutr
integer :: d_vvvo  ! created as (nutr,nou)

! MO integrals stored in d_vovv are of [p4 h* | p6 p5].  Like d_vvvo the only symmetry
! is exchange between p4 and p5.  this matrix must be stored in square block form
! note: it is possible to extract integrals in this order from d_vvvo with multiple
! small gets.  however, unless distributed memory storage is at a premium, it is more
! efficient to reorder and store the quantity in expanded block form.  it seems likely
! that an alternative code path could be taken if d_vovv.eq.-1 ... 
integer :: d_vovv ! created as (nou,nu2) .. only needed for abc tuples


integer :: ddi_np, ddi_me
integer :: ddi_nn, ddi_my
integer :: smp_np, smp_me
integer :: gpu_nd    ! on the node

integer :: global_comm, global_smp_comm, global_compute_comm
integer :: hybrid_comm, hybrid_smp_comm, hybrid_compute_comm
integer :: working_comm, working_smp_comm, working_compute_comm

integer :: flops

double precision :: ets, etd

double precision :: zero, one, two, four, eight, om, half
parameter(zero=0.0D+00,one=1.0D+00,two=2.0D+00,four=4.0D+00,eight=8.0D+00,om=-1.0D+00,half=0.5D+00)

! --------------------------------------------------------------------------------------------------
! common_cc subroutine definitions
! --------------------------------------------------------------------------------------------------
  contains

  subroutine common_cc_init(nocc,nvir)
  implicit none
  integer :: nocc, nvir
  no = nocc
  nu = nvir
  nou = no*nu
  no2 = no*no
  no3 = no*no2
  no4 = no*no3
  nu2 = nu*nu
  nu3 = nu*nu2
  nu4 = nu*nu3
  no2u = no2*nu
  no3u = no3*nu
  nou2 = no*nu2
  nou3 = no*nu3
  no2u2 = no2*nu2
  notr = (no2+no)/2
  nutr = (nu2+nu)/2
  ets = zero
  etd = zero
  flops = 0
  call ddi_nproc(ddi_np,ddi_me)
  call ddi_nnode(ddi_nn,ddi_my)
  call ddi_smp_nproc(smp_np,smp_me)
  return
  end subroutine common_cc_init


  subroutine sync(comm)
  implicit none
  integer :: comm
  integer :: ok
  integer :: ierr
  ok = 0
  if(ok.eq.1 .or. comm.eq.global_smp_comm)     ok=1
  if(ok.eq.1 .or. comm.eq.global_compute_comm) ok=1
  if(ok.eq.1 .or. comm.eq.hybrid_smp_comm)     ok=1
  if(ok.eq.1 .or. comm.eq.hybrid_compute_comm) ok=1
  if(ok.ne.1) then
     write(6,*) 'unknown comm'
     stop
  end if
  call mpi_barrier(comm,ierr)
  end subroutine sync
end module


  subroutine trant3_1(n,v)
  use common_cc, only: smp_np, smp_me
  implicit none
  integer :: n, a, b, c
  double precision :: v(n,n,n), x
  integer :: icntr,nr,sr,ltr

  icntr = 0
  ltr = (n*n-n)/2
  call div_even(ltr,smp_np,smp_me,nr,sr)
  if(smp_np.gt.1) call smp_sync()

      DO B=1,N
      DO C=1,B-1
         icntr = icntr+1
         if(icntr.lt.sr) cycle
         if(icntr.ge.sr+nr) cycle
      DO A=1,N
         X=V(A,B,C)
         V(A,B,C)=V(A,C,B)
         V(A,C,B)=X
      end do
      end do
      end do

  if(smp_np.gt.1) call smp_sync()
  return
  end subroutine trant3_1

  subroutine trant3_4(n,v)
  use common_cc, only: smp_np, smp_me
  implicit none
  integer :: n, a, b, c 
  double precision :: v(n,n,n), x

  if(smp_np.gt.1) call smp_sync()
      DO 401 B=1,N
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

  if(smp_np.gt.1) call smp_sync()
  return
  end subroutine trant3_4

  subroutine trant3_5(n,v)
  use common_cc, only: smp_np, smp_me
  implicit none
  integer :: n, a, c, d
  double precision :: v(n,n,n), x
  if(smp_np.gt.1) call smp_sync()

      DO 501 A=1,N
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

  if(smp_np.gt.1) call smp_sync()
  return
  end subroutine trant3_5


  subroutine tranmd_23(a,n1,n2,n3,n4)
  use common_cc, only: smp_np, smp_me
  implicit double precision(a-h,o-z)
  integer :: n1,n2,n3,n4,n12,n123,icntr
  double precision :: a(*)

  n12=n1*n2
  n123=n12*n3
  icntr=0

  if(smp_np.gt.1) call smp_sync()

   23 CONTINUE
      DO J=1,N2
         J1=J-1
         N1J=N1*J1
         N12J=N12*J1
      DO K=1,J
         K1=K-1
         N1K=N1*K1
         N12K=N12*K1
         icntr=icntr+1
         if(icntr.eq.smp_np) icntr=0
         if(icntr.ne.smp_me) cycle
      DO L=1,N4
         N123L=N123*(L-1)
      DO I=1,N1
         N123LI=N123L+I
         IJKL=N123LI+N12K+N1J
         IKJL=N123LI+N12J+N1K
         X=A(IJKL)
         A(IJKL)=A(IKJL)
         A(IKJL)=X
      end do
      end do
      end do
      end do

  if(smp_np.gt.1) call smp_sync()

  return
  end subroutine tranmd_23
  
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


      SUBROUTINE TRANT3_SMP(V,NU,ID)
      use common_cc, only: smp_np, smp_me
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU)

      if(smp_np.gt.1) CALL smp_sync()

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
      if(smp_np.gt.1) CALL smp_sync()
      RETURN
      END


subroutine smp_sync()
use common_cc
implicit none
call sync(working_smp_comm)
return
end subroutine smp_sync
