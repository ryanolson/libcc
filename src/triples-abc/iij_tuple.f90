subroutine iij_tuple_formv3(i,j,t2_i,t2_j,vm_ij,vm_ji,vm_ii,ve_i,ve_j,v3)
use common_cc
implicit none

integer :: i, j
double precision :: t2_i(nu2,no), t2_j(nu2,no)
double precision :: vm_ij(nou), vm_ji(nou), vm_ii(nou)
double precision :: ve_i(nu3), ve_j(nu3)
double precision :: v3(nu3)

integer :: nr,sr,t3off
double precision :: zero,one,om
parameter(zero=0.0D+00,one=1.0D+00,om=-1.0D+00)

call div_even(nu2,smp_np,smp_me,nr,sr)


call dgemm('n','n',nr,nu,no,om,t2_i(sr,1),nu2,vm_ji,no,zero,v3(sr),nu2)
call trant3_1(v3)
!     T2OFF = (I-1)*NU2*NO + SR
!     CALL DCOPY(NOU,VM(1,1,J,I),1,TEMP,1)
!     CALL DGEMM('N','N',NR,NU,NO,OM,T2(T2OFF),NU2,TEMP,
!    *NO,ZERO,V3(SR),NU2)
!     CALL TRANT3_SMP(V3,NU,1)


call dgemm('n','n',nr,nu,no,om,t2_i(sr,1),nu2,vm_ij,no,one,v3(sr),nu2)
call trant3_smp(v3,nu,2)
!     CALL DCOPY(NOU,VM(1,1,I,J),1,TEMP,1)
!     CALL DGEMM('N','N',NR,NU,NO,OM,T2(T2OFF),NU2,TEMP,
!    *NO,ONE,V3(SR),NU2)
!     CALL TRANT3_SMP(V3,NU,2)


call dgemm('n','n',nr,nu,no,om,t2_j(sr,1),nu2,vm_ii,no,one,v3(sr),nu2)
call trant3_smp(v3,nu,2)
!     T2OFF = (J-1)*NU2*NO + SR
!     CALL DCOPY(NOU,VM(1,1,I,I),1,TEMP,1)
!     CALL DGEMM('N','N',NR,NU,NO,OM,T2(T2OFF),NU2,TEMP,
!    *NO,ONE,V3(SR),NU2)
!     CALL TRANT3_SMP(V3,NU,2)

! ve_j is already present
!     IF(SMP_ME.EQ.0) THEN
!     CALL DDCC_T_GETVE(NU,J,TEMP,T3)
!     END IF
!     CALL TRANMD_SMP(T3,NU,NU,NU,1,23)
!     CALL DDI_SMP_SYNC()
!
      T3OFF = (SR-1)*NU + 1


call dgemm('n','n',nu,nr,nu,one,t2_i(1,i),nu,ve_j(t3off),nu,one,v3(t3off),nu)
call trant3_1(v3)
!     T2OFF = (I-1)*NU2*NO + (I-1)*NU2 + 1
!     CALL DCOPY(NU2,T2(T2OFF),1,TEMP,1)
!     CALL DGEMM('N','N',NU,NR,NU,ONE,TEMP,NU,T3(T3OFF),
!    *NU,ONE,V3(T3OFF),NU)
!     CALL TRANT3_SMP(V3,NU,1)

! ve_i is already present
!     IF(SMP_ME.EQ.0) THEN
!     CALL DDCC_T_GETVE(NU,I,TEMP,T3)
!     END IF
!     CALL TRANMD_SMP(T3,NU,NU,NU,1,23)
!     CALL DDI_SMP_SYNC()


call dgemm('n','n',nu,nr,nu,one,t2_i(1,j),nu,ve_i(t3off),nu,one,v3(t3off),nu)
call trant3_smp(v3,nu,3)
!     T2OFF = (I-1)*NU2*NO + (J-1)*NU2 + 1
!     CALL DCOPY(NU2,T2(T2OFF),1,TEMP,1)
!     CALL DGEMM('N','N',NU,NR,NU,ONE,TEMP,NU,T3(T3OFF),
!    *NU,ONE,V3(T3OFF),NU)
!     CALL TRANT3_SMP(V3,NU,3)


call dgemm('n','n',nu,nr,nu,one,t2_j(1,i),nu,ve_i(t3off),nu,one,v3(t3off),nu)
call trant3_smp(v3,nu,3)
!     T2OFF = (J-1)*NU2*NO + (I-1)*NU2 + 1
!     CALL DCOPY(NU2,T2(T2OFF),1,TEMP,1)
!     CALL DGEMM('N','N',NU,NR,NU,ONE,TEMP,NU,T3(T3OFF),
!    *NU,ONE,V3(T3OFF),NU)
!     CALL TRANT3_SMP(V3,NU,3)


return
end subroutine iij_tuple_formv3

