
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

      SUBROUTINE T3SQUA_SMP(I,J,K,NORM,NURM,T1,T2,T3,EH,EP)
      use common_cc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IEJ,JEK
      INTEGER A,B,C
      DIMENSION T1(NO,NU),T2(NU,NU,NO,NO),T3(NU,NU,NU),EH(NO),EP(NU)
      integer blocksize,aa,bb,cc
      parameter(blocksize=16)
C -RESTART-      common/lame  / my_first_time

      CALL smp_sync()

      nblocks = nu/blocksize
      DIJK=EH(I)+EH(J)+EH(K)
      X1=ZERO
      X2=ZERO
      X3=ZERO
      icntr=0

      do cc = 1,nu,16
         icntr = icntr+1
         if(icntr.eq.smp_np) icntr=0
         if(icntr.ne.smp_me) cycle
      do bb = 1,nu,16
      do aa = 1,nu,16
         do c = cc,min(cc+15,nu)
            dc = ep(c)
         do b = bb,min(bb+15,nu)
            dbc = ep(b) + dc
         do a = aa,min(aa+15,nu)
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
      end do
      end do
      end do

      CF=ONE
      IEJ=I.EQ.J
      JEK=J.EQ.K
      IF(IEJ.OR.JEK) CF=HALF
      ETD=ETD+CF*X3
      CALL smp_sync()

C -RESTART-      if(ddi_me.eq.0 .and. my_first_time.eq.0) then
C -RESTART-         ifile=80
C -RESTART-         open(unit=ifile, file='t3squa.restart', status='new', 
C -RESTART-     &   action='write', form='unformatted', access='sequential')
C -RESTART-
C -RESTART-         write(ifile) nu
C -RESTART-         write(ifile) dijk
C -RESTART-         write(ifile) ep
C -RESTART-         write(ifile) t3
C -RESTART-         write(ifile) x3
C -RESTART-
C -RESTART-         close(ifile)
C -RESTART-         my_first_time = 1
C -RESTART-      end if

      RETURN
      END
