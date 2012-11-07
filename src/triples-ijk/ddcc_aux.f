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
      SUBROUTINE DDCC_T_GETVE_acc(acc_sync,NUINP,INDX,TEMP,VE)
      use common_cc
      INTEGER INDX,I,IOFF,acc_sync
C
C     TEMP SHOULD BE NU^2 AND VE SHOULD BE NU^3

      DOUBLE PRECISION TEMP(NUTR,NU),VE(NU,NU,NU)
      INTEGER ij_1,ij_2

!$acc parallel loop async(acc_sync)
      DO IOFF = 1, NU
        ij_1 = 1
        ij_2 = 1
        do i = 1,NU
          do j = 1,i
             VE(i,j,IOFF) = TEMP(ij_1,IOFF)
             ij_1 = ij_1 + 1
          end do
          do j = 1,i
             VE(j,i,IOFF) = TEMP(ij_2,IOFF)
             ij_2 = ij_2 + 1
          end do
        end do

      END DO
!$acc end parallel loop

      END


      SUBROUTINE ZEROT3_ACC(T3,NU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER A
      DIMENSION T3(NU,NU,NU)
      DATA ZERO/0.0D+00/

!$acc parallel loop 
      DO  A=1,NU
      T3(A,A,A)=ZERO
      ENDDO           
!$acc end parallel loop

      END 


      SUBROUTINE ADT3DEN_ACC(NU,DEH,T3,EP)
      use common_cc, only: smp_np, smp_me
      IMPLICIT NONE
      INTEGER A,B,C,NU
      DOUBLE PRECISION T3(NU,NU,NU),EP(NU),DEH,DEN
C
      if(smp_np.gt.1) CALL smp_sync()

!$acc parallel loop 
      DO 11 A=1,NU
#ifndef USE_OPEN_ACC
        IF(MOD(A,SMP_NP).NE.SMP_ME) GOTO 11
#endif
      DO 10 B=1,NU
      DO 10 C=1,NU
      DEN=DEH-EP(A)-EP(B)-EP(C)
      T3(C,B,A)=T3(C,B,A)/DEN
 10   CONTINUE
 11   CONTINUE
!$acc end parallel loop

      if(smp_np.gt.1) CALL smp_sync()

      END

      SUBROUTINE DRT1WT3IJK_ACC(I,J,K,NO,NU,T1,VOE,TI,T3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T1(1),VOE(1),T3(1),TI(1)
C
      CALL T1WT3IJK_ACC(I,J,K,NO,NU,T1,VOE,TI,T3)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_ACC(I,K,J,NO,NU,T1,VOE,TI,T3)
      call trant3_4(nu,ti)
      CALL T1WT3IJK_ACC(J,I,K,NO,NU,T1,VOE,TI,T3)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_ACC(J,K,I,NO,NU,T1,VOE,TI,T3)
      call trant3_5(nu,ti)
      CALL T1WT3IJK_ACC(K,I,J,NO,NU,T1,VOE,TI,T3)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_ACC(K,J,I,NO,NU,T1,VOE,TI,T3)
      RETURN
      END

      SUBROUTINE T1WT3IJK_ACC(I,J,K,NO,NU,T1,VOE,TI,T3)
      use common_cc, only: smp_np, smp_me, no2u2, nu3
      IMPLICIT NONE
      INTEGER NR,SR,T3OFF,VOEOFF,I,J,K,NO,NU,NU2
      DOUBLE PRECISION T1(NU,NO),VOE(no2u2),T3(nu3),TI(nu3),ONE
      DATA ONE/1.0D+00/
C
      NU2 = NU*NU
      CALL SMT3FOUR_ACC(NU,T3,TI)
      CALL DIV_EVEN(NU2,SMP_NP,SMP_ME,NR,SR)
      T3OFF = (SR - 1)*NU + 1 
      VOEOFF = (K - 1)*NU2*NO + (J - 1)*NU2  + SR

!$acc data present(t3,voe,t1)
!$acc host_data use_device(t3,voe,t1)
      CALL DGEMM('N','N',NU,1,NR,ONE,T3(T3OFF),NU,VOE(VOEOFF),NU2,ONE,
     *           T1(1,I),NU)
!$acc end host_data
!$acc end data
      RETURN
      END

      SUBROUTINE SMT3FOUR_ACC(NU,T3,V3)
      use common_cc, only: smp_np, smp_me
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      DIMENSION T3(NU,NU,NU),V3(NU,NU,NU)
      DATA TWO/2.0D+00/
      DATA ZERO/0.0D+00/
      DATA OM/-1.0D+00/
C
!$acc parallel loop
      DO 2 C=1,NU
      do b=1,nu
         T3(b,b,c) = zero
      end do

      DO B=2,NU
      DO A=1,B-1
         T3(A,B,C)=(V3(A,B,C)-V3(B,A,C))*TWO-V3(A,C,B)+V3(B,C,A)
         T3(B,A,C)=T3(A,B,C)*OM
      END DO
      END DO

 2    CONTINUE
!$acc end parallel loop
      RETURN
      END

      SUBROUTINE DRT1WT3IJ_ACC(I,J,NO,NU,T1,VOE,TI,T3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T1(1),VOE(1),T3(1),TI(1)
C
      CALL T1WT3IJK_ACC(I,J,J,NO,NU,T1,VOE,TI,T3)
      CALL TRANT3_ACC(TI,NU,2)
      CALL T1WT3IJK_ACC(J,I,J,NO,NU,T1,VOE,TI,T3)
c     CALL TRANT3_ACC(TI,NU,1)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_ACC(J,J,I,NO,NU,T1,VOE,TI,T3)
      RETURN
      END

      SUBROUTINE DRT1WT3JK_ACC(J,K,NO,NU,T1,VOE,TI,T3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T1(1),VOE(1),T3(1),TI(1)
C
      CALL T1WT3IJK_ACC(J,J,K,NO,NU,T1,VOE,TI,T3)
C     CALL TRANT3_ACC(TI,NU,1)
      call trant3_1(nu,ti)
      CALL T1WT3IJK_ACC(J,K,J,NO,NU,T1,VOE,TI,T3)
      CALL TRANT3_ACC(TI,NU,2)
      CALL T1WT3IJK_ACC(K,J,J,NO,NU,T1,VOE,TI,T3)
      RETURN
      END

      SUBROUTINE SYMT311_ACC(V,NU,ID)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,B,C
      double precision V(NU,NU,NU),x

      IF (ID.EQ.23) go to 23
      IF (ID.EQ.12) go to 12
      IF (ID.EQ.13) go to 13

 23   CONTINUE

!$acc parallel loop private(x)
      DO B=1,NU
      DO C=1,B
      DO A=1,NU
        X=V(A,B,C)+V(A,C,B)
        V(A,B,C)=X
        V(A,C,B)=X
      end do
      end do
      end do
!$acc end parallel loop
      go to 1000

 12   CONTINUE
!$acc parallel loop private(x)
      DO C=1,NU
!$acc loop
      DO A=1,NU
      DO B=1,A
        X=V(A,B,C)+V(B,A,C)
        V(A,B,C)=X
        V(B,A,C)=X
      end do
      end do
      end do
!$acc end parallel loop
      go to 1000

 13   CONTINUE
!$acc parallel loop private(x)
      DO A=1,NU
!$acc loop
      DO B=1,NU
      DO C=1,A
        X=V(A,B,C)+V(C,B,A)
        V(A,B,C)=X
        V(C,B,A)=X
      end do
      end do
      end do
!$acc end parallel loop

 1000 CONTINUE

      END

!     previous version:
!
!     SUBROUTINE SYMT311(V,NU,ID)
!     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     INTEGER A,B,C
!     DIMENSION V(NU,NU,NU)
!     
!     IF (ID.EQ.23) GO TO 23
!     IF (ID.EQ.12) GO TO 12
!     IF (ID.EQ.13) GO TO 13
!23   CONTINUE
!     DO 100 B=1,NU
!     DO 100 C=1,B
!     DO 100 A=1,NU
!     X=V(A,B,C)+V(A,C,B)
!     V(A,B,C)=X
!     V(A,C,B)=X
! 100 CONTINUE
!     GO TO 1000
!12   CONTINUE
!     DO 101 C=1,NU
!     DO 101 A=1,NU
!     DO 101 B=1,A
!     X=V(A,B,C)+V(B,A,C)  
!     V(A,B,C)=X
!     V(B,A,C)=X
!101  CONTINUE
!     GO TO 1000
!13   CONTINUE
!     DO 102 A=1,NU
!     DO 102 B=1,NU
!     DO 102 C=1,A
!     X=V(A,B,C)+V(C,B,A)
!     V(A,B,C)=X
!     V(C,B,A)=X
!102  CONTINUE
!     GO TO 1000
!1000 CONTINUE
!     END
