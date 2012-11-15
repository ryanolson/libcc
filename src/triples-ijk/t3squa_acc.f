C =============================================================================
C   Copyright (C) 2010.  Ryan M. Olson
C
C   This program is free software: you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation, either version 3 of the License, or
C   (at your option) any later version.
C
C   This program is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program.  If not, see <http://www.gnu.org/licenses/>.
C =============================================================================

      SUBROUTINE T3SQUA_ACC(I,J,K,NORM,NURM,T1,T2,T3,EH,EP)
      use common_cc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IEJ,JEK
      INTEGER A,B,C
      DIMENSION T1(NO,NU),T2(NU,NU,NO,NO),T3(NU,NU,NU),EH(NO),EP(NU)
C -RESTART-      common/lame  / my_first_time

      DIJK=EH(I)+EH(J)+EH(K)
      X3=ZERO

!$acc parallel loop private(dc,dbc,DABC,DENOM,d1,d2,d3,f) 
!$acc&         reduction(+:x3)
         do c = 1,nu
            dc = ep(c)
         do b = 1,nu
            dbc = ep(b) + dc
         do a = 1,nu
            if(a.eq.b .and. b.eq.c) cycle
            DABC = EP(A) + DBC
            DENOM=DIJK-DABC
            DENOM=1/DENOM
            D1=  T3(A,B,C)
            D2=  T3(A,C,B)+T3(C,B,A)+T3(B,A,C)
            D3=  T3(B,C,A)+T3(C,A,B)
            F=D1*EIGHT-FOUR*D2+D3*TWO
            X3=X3+F*D1*DENOM
         end do
         end do
         end do
!$acc end parallel loop

      CF=ONE
      IEJ=I.EQ.J
      JEK=J.EQ.K
      IF(IEJ.OR.JEK) CF=HALF
      ETD=ETD+CF*X3

      END
