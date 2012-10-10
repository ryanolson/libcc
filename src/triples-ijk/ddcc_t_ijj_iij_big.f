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
      SUBROUTINE DDCC_T_IJJ_BIG(NO,NU,I,J,T1,T2,VM,
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


C
C-IJJ-      CALL DIV_EVEN(NU2,SMP_NP,SMP_ME,NR,SR)
C-IJJ-C
C-IJJ-      T2OFF = (I-1)*NU2*NO + SR
C-IJJ-      CALL DDI_SMP_SYNC()
C-IJJ-      CALL DCOPY(NOU,VM(1,1,J,J),1,TEMP,1)
C-IJJ-      CALL DGEMM('N','N',NR,NU,NO,OM,T2(T2OFF),NU2,TEMP,
C-IJJ-     *NO,ZERO,V3(SR),NU2)
C-IJJ-C
C-IJJ-      CALL TRANT3_SMP(V3,NU,2)
C-IJJ-C
C-IJJ-      T2OFF = (J-1)*NU2*NO + SR
C-IJJ-      CALL DCOPY(NOU,VM(1,1,J,I),1,TEMP,1)
C-IJJ-      CALL DGEMM('N','N',NR,NU,NO,OM,T2(T2OFF),NU2,TEMP,
C-IJJ-     *NO,ONE,V3(SR),NU2)
C-IJJ-C
C-IJJ-      CALL TRANT3_SMP(V3,NU,1)
C-IJJ-C
C-IJJ-      CALL DCOPY(NOU,VM(1,1,I,J),1,TEMP,1)
C-IJJ-      CALL DGEMM('N','N',NR,NU,NO,OM,T2(T2OFF),NU2,TEMP,
C-IJJ-     *NO,ONE,V3(SR),NU2)
C-IJJ-C
C-IJJ-      CALL TRANT3_SMP(V3,NU,4)
C-IJJ-      IF(SMP_ME.EQ.0) THEN
C-IJJ-      CALL DDCC_T_GETVE(NU,J,TEMP,T3)
C-IJJ-      END IF
C-IJJ-      CALL TRANMD_SMP(T3,NU,NU,NU,1,23)
C-IJJ-      CALL DDI_SMP_SYNC()
C-IJJ-C
C-IJJ-      T3OFF = (SR-1)*NU + 1
C-IJJ-C
C-IJJ-      T2OFF = (I-1)*NU2*NO + (J-1)*NU2 + 1
C-IJJ-      CALL DCOPY(NU2,T2(T2OFF),1,TEMP,1)
C-IJJ-      CALL DGEMM('N','N',NU,NR,NU,ONE,TEMP,NU,T3(T3OFF),
C-IJJ-     *NU,ONE,V3(T3OFF),NU)
C-IJJ-C
C-IJJ-      CALL TRANT3_SMP(V3,NU,3)
C-IJJ-C
C-IJJ-      T2OFF = (J-1)*NU2*NO + (I-1)*NU2 + 1
C-IJJ-      CALL DCOPY(NU2,T2(T2OFF),1,TEMP,1)
C-IJJ-      CALL DGEMM('N','N',NU,NR,NU,ONE,TEMP,NU,T3(T3OFF),
C-IJJ-     *NU,ONE,V3(T3OFF),NU)
C-IJJ-C
C-IJJ-      CALL TRANT3_SMP(V3,NU,1)
C-IJJ-      IF(SMP_ME.EQ.0) THEN
C-IJJ-      CALL DDCC_T_GETVE(NU,I,TEMP,T3)
C-IJJ-      END IF
C-IJJ-      CALL TRANMD_SMP(T3,NU,NU,NU,1,23)
C-IJJ-      CALL DDI_SMP_SYNC()
C-IJJ-C
C-IJJ-      T2OFF = (J-1)*NU2*NO + (J-1)*NU2 + 1
C-IJJ-      CALL DCOPY(NU2,T2(T2OFF),1,TEMP,1)
C-IJJ-      CALL DGEMM('N','N',NU,NR,NU,ONE,TEMP,NU,T3(T3OFF),
C-IJJ-     *NU,ONE,V3(T3OFF),NU)
C-IJJ-C
C-IJJ-      CALL TRANT3_SMP(V3,NU,2)
C-IJJ-
      IF(SMP_ME.EQ.0) THEN
      CALL SYMT311(V3,NU,23)      ! modifies V3
      CALL ZEROT3(V3,NU)          ! modifies V3
      END IF
      CALL T3SQUA_SMP(I,J,J,NO,NU,O1,T2,V3,EH,EP)  ! no array change
      DEH=EH(I)+EH(J)+EH(J)
      CALL ADT3DEN_SMP(NU,DEH,V3,EP)  ! modifies V3
      ITMP=I
      JTMP=J
C
      CALL DRT1WT3IJ_SMP(ITMP,JTMP,NO,NU,T1,VOE,V3,T3) ! modifies V3, T1
      if(smp_np.gt.1) CALL smp_sync()
      RETURN
      END
C
C*MODULE DDICC   *DECK DDCC_T_IIJ
      SUBROUTINE DDCC_T_IIJ_BIG(NO,NU,I,J,T1,T2,VM,
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

C-IIJ-      CALL DIV_EVEN(NU2,SMP_NP,SMP_ME,NR,SR)
C-IIJ-C
C-IIJ-      CALL DDI_SMP_SYNC()
C-IIJ-C
C-IIJ-      T2OFF = (I-1)*NU2*NO + SR
C-IIJ-      CALL DCOPY(NOU,VM(1,1,J,I),1,TEMP,1)
C-IIJ-      CALL DGEMM('N','N',NR,NU,NO,OM,T2(T2OFF),NU2,TEMP,
C-IIJ-     *NO,ZERO,V3(SR),NU2)
C-IIJ-C
C-IIJ-      CALL TRANT3_SMP(V3,NU,1)
C-IIJ-C
C-IIJ-      CALL DCOPY(NOU,VM(1,1,I,J),1,TEMP,1)
C-IIJ-      CALL DGEMM('N','N',NR,NU,NO,OM,T2(T2OFF),NU2,TEMP,
C-IIJ-     *NO,ONE,V3(SR),NU2)
C-IIJ-C
C-IIJ-      CALL TRANT3_SMP(V3,NU,2)
C-IIJ-C
C-IIJ-      T2OFF = (J-1)*NU2*NO + SR
C-IIJ-      CALL DCOPY(NOU,VM(1,1,I,I),1,TEMP,1)
C-IIJ-      CALL DGEMM('N','N',NR,NU,NO,OM,T2(T2OFF),NU2,TEMP,
C-IIJ-     *NO,ONE,V3(SR),NU2)
C-IIJ-C
C-IIJ-      CALL TRANT3_SMP(V3,NU,2)
C-IIJ-      IF(SMP_ME.EQ.0) THEN
C-IIJ-      CALL DDCC_T_GETVE(NU,J,TEMP,T3)
C-IIJ-      END IF
C-IIJ-      CALL TRANMD_SMP(T3,NU,NU,NU,1,23)
C-IIJ-      CALL DDI_SMP_SYNC()
C-IIJ-C
C-IIJ-      T3OFF = (SR-1)*NU + 1
C-IIJ-C
C-IIJ-      T2OFF = (I-1)*NU2*NO + (I-1)*NU2 + 1
C-IIJ-      CALL DCOPY(NU2,T2(T2OFF),1,TEMP,1)
C-IIJ-      CALL DGEMM('N','N',NU,NR,NU,ONE,TEMP,NU,T3(T3OFF),
C-IIJ-     *NU,ONE,V3(T3OFF),NU)
C-IIJ-C
C-IIJ-      CALL TRANT3_SMP(V3,NU,1)
C-IIJ-      IF(SMP_ME.EQ.0) THEN
C-IIJ-      CALL DDCC_T_GETVE(NU,I,TEMP,T3)
C-IIJ-      END IF
C-IIJ-      CALL TRANMD_SMP(T3,NU,NU,NU,1,23)
C-IIJ-      CALL DDI_SMP_SYNC()
C-IIJ-C
C-IIJ-      T2OFF = (I-1)*NU2*NO + (J-1)*NU2 + 1
C-IIJ-      CALL DCOPY(NU2,T2(T2OFF),1,TEMP,1)
C-IIJ-      CALL DGEMM('N','N',NU,NR,NU,ONE,TEMP,NU,T3(T3OFF),
C-IIJ-     *NU,ONE,V3(T3OFF),NU)
C-IIJ-C
C-IIJ-      CALL TRANT3_SMP(V3,NU,3)
C-IIJ-C
C-IIJ-      T2OFF = (J-1)*NU2*NO + (I-1)*NU2 + 1
C-IIJ-      CALL DCOPY(NU2,T2(T2OFF),1,TEMP,1)
C-IIJ-      CALL DGEMM('N','N',NU,NR,NU,ONE,TEMP,NU,T3(T3OFF),
C-IIJ-     *NU,ONE,V3(T3OFF),NU)
C-IIJ-C
C-IIJ-      CALL TRANT3_SMP(V3,NU,3)

      IF(SMP_ME.EQ.0) THEN
      CALL SYMT311(V3,NU,12)  ! modifies V3
      CALL ZEROT3(V3,NU)      ! modifies V3
C     IF(IDISC.EQ.0.AND.MET.GT.4) THEN
C       CALL WRT3(KK,NU,V3)
C     END IF
      END IF

      CALL T3SQUA_SMP(I,I,J,NO,NU,O1,T2,V3,EH,EP)  ! no array mods

      DEH=EH(I)+EH(I)+EH(J)
      CALL ADT3DEN_SMP(NU,DEH,V3,EP)  ! modifies V3

      ITMP=I
      JTMP=J

      CALL DRT1WT3JK_SMP(ITMP,JTMP,NO,NU,T1,VOE,V3,T3) ! modifies V3, T1, T3
      if(smp_np.gt.1) CALL smp_sync()
      RETURN
      END
