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
  implicit none
  integer :: n, a, b, c
  double precision :: v(n,n,n), x
  integer :: icntr,nr,sr,ltr

!$acc parallel loop private(x)
      DO B=1,N
      DO C=1,B-1
      DO A=1,N
         X=V(A,B,C)
         V(A,B,C)=V(A,C,B)
         V(A,C,B)=X
      end do
      end do
      end do
!$acc end parallel loop

  return
  end subroutine trant3_1

  subroutine trant3_3(n,v)
  implicit none
  integer :: n, a, b, c
  double precision :: v(n,n,n), x
  integer :: icntr,nr,sr,ltr

!$acc parallel loop private(x)
      DO B=1,N
      DO A=1,N
      DO C=1,A
         X=V(A,B,C)
         V(A,B,C)=V(C,B,A)
         V(C,B,A)=X
      end do
      end do
      end do
!$acc end parallel loop

  return
  end subroutine trant3_3

  subroutine trant3_4(n,v)
  implicit none
  integer :: n, a, b, c 
  double precision :: v(n,n,n), x


!$acc parallel loop private(x)
      DO 401 B=1,N
!$acc loop
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
!$acc end parallel loop

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

  subroutine tranmd_23_acc(acc_sync,a,n1,n2,n3,n4)
  implicit double precision(a-h,o-z)
  integer :: n1,n2,n3,n4,n12,n123,icntr,acc_sync
  double precision :: a(n1*n2*n3)

  n12=n1*n2
  n123=n12*n3
  icntr=0

!$acc parallel loop private(x) async(acc_sync)
      DO J=1,N2
         J1=J-1
         N1J=N1*J1
         N12J=N12*J1
      DO K=1,J
         K1=K-1
         N1K=N1*K1
         N12K=N12*K1
      DO L=1,N4
         N123L=N123*(L-1)
!$acc loop
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
!$acc end parallel loop

  end subroutine tranmd_23_acc
  

      SUBROUTINE TRANT3_ACC(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C,D
      DIMENSION V(NU,NU,NU)

      GO TO (1,2,3,4,5),ID
    1 CONTINUE

!$acc parallel loop private(x) present(v)
      DO 101 B=1,NU
      DO 100 C=1,B
      DO 100 A=1,NU
      X=V(A,B,C)
      V(A,B,C)=V(A,C,B)
      V(A,C,B)=X
  100 CONTINUE
  101 CONTINUE
!$acc end parallel loop

      GO TO 1000
    2 CONTINUE

!$acc parallel loop private(x) present(v)
      DO 201 C=1,NU
!$acc loop
      DO 200 A=1,NU
      DO 200 B=1,A
      X=V(A,B,C)
      V(A,B,C)=V(B,A,C)
      V(B,A,C)=X
  200 CONTINUE
  201 CONTINUE
!$acc end parallel loop

      GO TO 1000
    3 CONTINUE

!$acc parallel loop private(x) present(v)
      DO 301 B=1,NU
!$acc loop
      DO 300 A=1,NU
      DO 300 C=1,A
      X=V(A,B,C)
      V(A,B,C)=V(C,B,A)
      V(C,B,A)=X
  300 CONTINUE
  301 CONTINUE
!$acc end parallel loop

      GO TO 1000
    4 CONTINUE

!$acc parallel loop private(x) present(v)
      DO 401 B=1,NU
!$acc loop
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
!$acc end parallel loop

      GO TO 1000
    5 CONTINUE

!$acc parallel loop private(x) present(v)
      DO 501 A=1,NU
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
!$acc end parallel loop

      GO TO 1000
 1000 CONTINUE
      RETURN
      END


subroutine smp_sync()
implicit none
call ddi_smp_sync()
return
end subroutine smp_sync

  subroutine common_cc_init(nocc,nvir)
  use common_cc
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

