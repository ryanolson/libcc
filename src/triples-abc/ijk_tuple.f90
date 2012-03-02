subroutine ijk_tuple_formv3(i,j,k,t2_i,t2_j,t2_k,vm_ij,vm_ji,vm_ik,vm_ki,vm_jk,vm_kj,ve_i,ve_j,ve_k,v3)
use common_cc
implicit none

integer :: i, j, k
double precision :: t2_i(nu2,no), t2_j(nu2,no), t2_k(nu2,no)
double precision :: vm_ij(nou), vm_ji(nou), vm_ik(nou), vm_ki(nou), vm_jk(nou), vm_kj(nou)
double precision :: ve_i(nu3), ve_j(nu3), ve_k(nu3)
double precision :: v3(nu3)

double precision :: om,one,zero
parameter(zero=0.0D+00, om=-1.0D+00, one=1.0D+00)

integer nr, sr, veoff

call div_even(nu2,smp_np,smp_me,nr,sr)
veoff = (sr-1)*nu + 1

! #0 & #11
call dgemm('n','n',nr,nu,no,om,t2_i(sr,1),nu2,vm_kj,no,zero,v3(sr),nu2)       ! #0: Type A
call dgemm('t','t',nr,nu,nu,one,ve_i(veoff),nu,t2_k(1,j),nu,one,v3(sr),nu2)   ! #11: Type BT
call ddi_smp_sync()
! #8 & #3
call dgemm('n','n',nu,nr,nu,one,t2_i(1,k),nu,ve_j(veoff),nu,one,v3(veoff),nu) ! #8: Type B
call dgemm('t','t',nu,nr,no,om,vm_ik,no,t2_j(sr,1),nu2,one,v3(veoff),nu)      ! #3: Type AT
! transform v3
call trant3_1(v3)
! #1 & #10
call dgemm('n','n',nr,nu,no,om,t2_i(sr,1),nu2,vm_jk,no,one,v3(sr),nu2)        ! #1: Type A
call dgemm('t','t',nr,nu,nu,one,ve_i(veoff),nu,t2_j(1,k),nu,one,v3(sr),nu2)   ! #10: Type BT
call ddi_smp_sync()
! #6 & #5
call dgemm('n','n',nu,nr,nu,one,t2_i(1,j),nu,ve_k(veoff),nu,one,v3(veoff),nu) ! #6: Type B
call dgemm('t','t',nu,nr,no,om,vm_ij,no,t2_k(sr,1),nu2,one,v3(veoff),nu)      ! #5: Type AT
! transform v3
call trant3_4(v3)
! #2 & #9
call dgemm('n','n',nr,nu,no,om,t2_j(sr,1),nu2,vm_ki,no,one,v3(sr),nu2)        ! #2: Type A
call dgemm('t','t',nr,nu,nu,one,ve_j(veoff),nu,t2_k(1,i),nu,one,v3(sr),nu2)   ! #9: Type BT
call ddi_smp_sync()
! transform v3
call trant3_1(v3) 
! #4 & #7
call dgemm('t','t',nu,nr,no,om,vm_ji,no,t2_k(sr,1),nu2,one,v3(veoff),nu)      ! #4: Type AT
call dgemm('n','n',nu,nr,nu,one,t2_j(1,i),nu,ve_k(veoff),nu,one,v3(veoff),nu) ! #7: Type B
! transform v3
call trant3_4(v3)

if(smp_me.eq.0) call zerot3(v3,nu)
call ddi_smp_sync()

return
end subroutine ijk_tuple_formv3

subroutine trant3_1(v)
use common_cc
implicit none
double precision :: v(nu,nu,nu)
integer :: a, b, c
double precision :: x

call ddi_smp_sync()

      DO 101 B=1,NU
        IF(MOD(B,SMP_NP).NE.SMP_ME) GOTO 101
      DO 100 C=1,B
      DO 100 A=1,NU
      X=V(A,B,C)
      V(A,B,C)=V(A,C,B)
      V(A,C,B)=X
  100 CONTINUE
  101 CONTINUE

call ddi_smp_sync()

return
end subroutine trant3_1

subroutine trant3_4(v)
use common_cc
implicit none
double precision :: v(nu,nu,nu)
integer :: a, b, c
double precision :: x

call ddi_smp_sync()

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

call ddi_smp_sync()

return
end subroutine trant3_4

subroutine trant3_5(v)
use common_cc
implicit none
double precision :: v(nu,nu,nu)
integer :: a, c, d
double precision :: x

call ddi_smp_sync()

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

call ddi_smp_sync()

return
end subroutine trant3_5

