subroutine abc_tuple_formv3(a,b,c,t2_a,t2_b,t2_c,ve_ab,ve_ba,ve_ac,ve_ca,ve_bc,ve_cb,vm_a,vm_b,vm_c,v3)
use common_cc
implicit none

integer :: a, b, c
double precision :: t2_a(no2,nu), t2_b(no2,nu), t2_c(no2,nu)
double precision :: ve_ab(no2), ve_ba(no2), ve_ac(no2), ve_ca(no2), ve_bc(no2), ve_cb(no2)
double precision :: vm_a(no3), vm_b(no3), vm_c(no3)
double precision :: v3(no3)

double precision :: om,one,zero
parameter(zero=0.0D+00, om=-1.0D+00, one=1.0D+00)

integer nr, sr, veoff

! old call div_even(nu2,smp_np,smp_me,nr,sr)
  call div_even(no2,smp_np,smp_me,nr,sr)
  veoff = (sr-1)*no + 1

! #0
! call dgemm('n','n',nr,nu,no,om,t2_i(sr,1),nu2,vm_kj,no,zero,v3(sr),nu2)
  call dgemm('n','n',nr,no,nu,om,t2_a(sr,1),no2,ve_cb,nu,zero,v3(sr),no2)
  call ddi_smp_sync()
! #8
! call dgemm('n','n',nu,nr,nu,one,t2_i(1,k),nu,ve_j(veoff),nu,one,v3(veoff),nu)
  call dgemm('n','n',no,nr,no,one,t2_a(1,c),no,vm_b(veoff),no,one,v3(veoff),no)
! transform v3
  call trant3_1(no,v3)
! #1
! call dgemm('n','n',nr,nu,no,om,t2_i(sr,1),nu2,vm_jk,no,one,v3(sr),nu2)
  call dgemm('n','n',nr,no,nu,om,t2_a(sr,1),no2,ve_bc,nu,one,v3(sr),no2)
  call ddi_smp_sync()
! #6
! call dgemm('n','n',nu,nr,nu,one,t2_i(1,j),nu,ve_k(veoff),nu,one,v3(veoff),nu)
  call dgemm('n','n',no,nr,no,one,t2_a(1,b),no,vm_c(veoff),no,one,v3(veoff),no)
! transform v3
  call trant3_4(no,v3)
! #2
! call dgemm('n','n',nr,nu,no,om,t2_j(sr,1),nu2,vm_ki,no,one,v3(sr),nu2)
  call dgemm('n','n',nr,no,nu,om,t2_b(sr,1),no2,ve_ca,nu,one,v3(sr),no2)
  call ddi_smp_sync()
! #10
! call dgemm('n','n',nu,nr,nu,one,t2_j(1,k),nu,ve_i(veoff),nu,one,v3(veoff),nu)
  call dgemm('n','n',no,nr,no,one,t2_b(1,c),no,vm_a(veoff),no,one,v3(veoff),no)
! transform v3
  call trant3_1(no,v3) 
! #3
! call dgemm('n','n',nr,nu,no,om,t2_j(sr,1),nu2,vm_ik,no,one,v3(sr),nu2)
  call dgemm('n','n',nr,no,nu,om,t2_b(sr,1),no2,ve_ac,nu,one,v3(sr),no2)
  call ddi_smp_sync()
! #7
! call dgemm('n','n',nu,nr,nu,one,t2_j(1,i),nu,ve_k(veoff),nu,one,v3(veoff),nu)
  call dgemm('n','n',no,nr,no,one,t2_b(1,a),no,vm_c(veoff),no,one,v3(veoff),no)
! transform v3
  call trant3_5(no,v3)
! #4
! call dgemm('n','n',nr,nu,no,om,t2_k(sr,1),nu2,vm_ji,no,one,v3(sr),nu2)
  call dgemm('n','n',nr,no,nu,om,t2_c(sr,1),no2,ve_ba,nu,one,v3(sr),no2)
  call ddi_smp_sync()
! #11
! call dgemm('n','n',nu,nr,nu,one,t2_k(1,j),nu,ve_i(veoff),nu,one,v3(veoff),nu)
  call dgemm('n','n',no,nr,no,one,t2_c(1,b),no,vm_a(veoff),no,one,v3(veoff),no)
! transform v3
  call trant3_1(no,v3)
! #5
! call dgemm('n','n',nr,nu,no,om,t2_k(sr,1),nu2,vm_ij,no,one,v3(sr),nu2)
  call dgemm('n','n',nr,no,nu,om,t2_c(sr,1),no2,ve_ab,nu,one,v3(sr),no2)
  call ddi_smp_sync()
! #9 
! call dgemm('n','n',nu,nr,nu,one,t2_k(1,i),nu,ve_j(veoff),nu,one,v3(veoff),nu)
  call dgemm('n','n',no,nr,no,one,t2_c(1,a),no,vm_b(veoff),no,one,v3(veoff),no)

! return to original ordering
  call trant3_4(no,v3)
  call trant3_1(no,v3)

  if(smp_me.eq.0) call zerot3(v3,no)
  call ddi_smp_sync()

return
end subroutine abc_tuple_formv3

subroutine trant3_1(n,v)
use common_cc
implicit none
integer :: n
double precision :: v(n,n,n)
integer :: a, b, c
double precision :: x

call ddi_smp_sync()

      DO 101 B=1,N
        IF(MOD(B,SMP_NP).NE.SMP_ME) GOTO 101
      DO 100 C=1,B
      DO 100 A=1,N
      X=V(A,B,C)
      V(A,B,C)=V(A,C,B)
      V(A,C,B)=X
  100 CONTINUE
  101 CONTINUE

call ddi_smp_sync()

return
end subroutine trant3_1

subroutine trant3_4(n,v)
use common_cc
implicit none
integer :: n
double precision :: v(n,n,n)
integer :: a, b, c
double precision :: x

call ddi_smp_sync()

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

call ddi_smp_sync()

return
end subroutine trant3_4

subroutine trant3_5(n,v)
use common_cc
implicit none
integer :: n
double precision :: v(n,n,n)
integer :: a, c, d
double precision :: x

call ddi_smp_sync()

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

call ddi_smp_sync()

return
end subroutine trant3_5

