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
      SUBROUTINE DDCC_T_IJJ_ACC(NO,NU,I,J,T1,T2,VM,
     *                      V3,T3,VOE,O1,EH,EP,TEMP,ve_i,ve_j)
      use common_cc, only: smp_np, smp_me, ddi_me, nu3
      IMPLICIT NONE
C
      INTEGER I,J,NO,NU
C
      DOUBLE PRECISION VM(NO,NU,NO,NO),VOE(1),TEMP(1),
     *                 T1(1),T3(*),V3(*),O1(1),EH(NO),
     *                 EP(NU),T2(NU*NU,NO,NO),
     &                 VE_I(*),VE_J(*)
C
      INTEGER NU2,NOU,T2OFF,T3OFF
      INTEGER ITMP,JTMP,NR,SR
      DOUBLE PRECISION DEH
C
      DOUBLE PRECISION ZERO,ONE,OM
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,OM=-1.0D+00)
C
      NU2 = NU*NU
      NOU = NO*NU
C
      call ijj_tuple_formv3(i,j,t2(1,1,i),t2(1,1,j),
     &                      vm(1,1,i,j),vm(1,1,j,i),vm(1,1,j,j),
     &                      ve_i,ve_j,V3)


      CALL SYMT311_ACC(V3,NU,23)      ! modifies V3
      CALL ZEROT3_ACC(V3,NU)          ! modifies V3
      CALL T3SQUA_ACC(I,J,J,NO,NU,O1,T2,V3,EH,EP)  ! no array change
      DEH=EH(I)+EH(J)+EH(J)
      CALL ADT3DEN_ACC(NU,DEH,V3,EP)  ! modifies V3
      ITMP=I
      JTMP=J
      CALL DRT1WT3IJ_ACC(ITMP,JTMP,NO,NU,T1,VOE,V3,T3) ! modifies V3, T1
      RETURN
      END
C
C*MODULE DDICC   *DECK DDCC_T_IIJ
      SUBROUTINE DDCC_T_IIJ_ACC(NO,NU,I,J,T1,T2,VM,
     *                      V3,T3,VOE,O1,EH,EP,TEMP,ve_i,ve_j)
      use common_cc, only: smp_np, smp_me, ddi_me, nu3, no2u2
      IMPLICIT NONE
C
C
      INTEGER I,J,NO,NU
C
      DOUBLE PRECISION VM(NO,NU,NO,NO),VOE(1),TEMP(1),
     *                 T1(1),T3(*),V3(*),O1(1),EH(NO),
     *                 EP(NU),T2(NU*NU,NO,NO),
     &                 ve_i(*),ve_j(*)
C
      INTEGER NU2,NOU,SR,NR,T2OFF,T3OFF
      INTEGER ITMP,JTMP
      DOUBLE PRECISION DEH
C
      DOUBLE PRECISION ZERO,ONE,OM
      PARAMETER(ZERO=0.0D+00,ONE=1.0D+00,OM=-1.0D+00)
C
      NU2 = NU*NU
      NOU = NO*NU

      call iij_tuple_formv3(i,j,t2(1,1,i),t2(1,1,j),
     &                      vm(1,1,i,j),vm(1,1,j,i),vm(1,1,i,i),
     &                      ve_i,ve_j,v3)


      CALL SYMT311_ACC(V3,NU,12)  ! modifies V3
      CALL ZEROT3_ACC(V3,NU)      ! modifies V3
      CALL T3SQUA_ACC(I,I,J,NO,NU,O1,T2,V3,EH,EP)  ! no array mods
      DEH=EH(I)+EH(I)+EH(J)
      CALL ADT3DEN_ACC(NU,DEH,V3,EP)  ! modifies V3
      ITMP=I
      JTMP=J
      CALL DRT1WT3JK_ACC(ITMP,JTMP,NO,NU,T1,VOE,V3,T3) ! modifies V3, T1, T3
      RETURN
      END
