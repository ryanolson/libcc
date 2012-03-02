      SUBROUTINE T3SQUA_SMP(I,J,K,NORM,NURM,T1,T2,T3,EH,EP)
      use common_cc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL IEJ,JEK
      INTEGER A,B,C
      DIMENSION T1(NO,NU),T2(NU,NU,NO,NO),T3(NU,NU,NU),EH(NO),EP(NU)
      DATA TWO/2.0D+00/,FOUR/4.0D+00/,EIGHT/8.0D+00/,ZERO/0.0D+00/,
     *     HALF/0.5D+00/,ONE/1.0D+00/
      integer blocksize
      parameter(blocksize=16)
C -RESTART-      common/lame  / my_first_time

      CALL DDI_SMP_SYNC()

      nblocks = nu/blocksize
      DIJK=EH(I)+EH(J)+EH(K)
      X1=ZERO
      X2=ZERO
      X3=ZERO

! cache blocked loops for blocks not divisible by 16
      if(mod(nu,blocksize).gt.0) then
         nblocks=nblocks+1

      do 151 icblock = 1,nblocks
        if(mod(icblock,smp_np).ne.smp_me) goto 151
        icstart = blocksize*(icblock-1) + 1
      do 150 ibblock = 1,nblocks
        ibstart = blocksize*(ibblock-1) + 1
      do 150 iablock = 1,nblocks
        iastart = blocksize*(iablock-1) + 1
      DO 150 C=icstart, icstart+blocksize-1
      DC = EP(C)
      DO 150 B=ibstart, ibstart+blocksize-1
      DBC = EP(B) + DC
      DO 150 A=iastart, iastart+blocksize-1
      IF (A.EQ.B.AND.B.EQ.C) GOTO 150
      IF (A.GT.NU .OR. B.GT.NU .OR. C.GT.NU) GOTO 150
c
c performance counters for the inner loop
c 8 adds, 5 mults, 1 divide = 14 flops
c passes: nu3-nu
c
      DABC = EP(A) + DBC
      DENOM=DIJK-DABC
      DENOM=1/DENOM
      D1=  T3(A,B,C)
      D2=  T3(A,C,B)+T3(C,B,A)+T3(B,A,C)
      D3=  T3(B,C,A)+T3(C,A,B)
      F=D1*EIGHT-FOUR*D2+D3*TWO
      X3=X3+F*D1*DENOM
 150  CONTINUE
 151  CONTINUE

      else
 
! cache blocked loops for nu evenly divisible by 16
      do 161 icblock = 1,nblocks
        if(mod(icblock,smp_np).ne.smp_me) goto 161
        icstart = blocksize*(icblock-1) + 1
      do 160 ibblock = 1,nblocks
        ibstart = blocksize*(ibblock-1) + 1
      do 160 iablock = 1,nblocks
        iastart = blocksize*(iablock-1) + 1
!dir$ noblocking
      DO 160 C=icstart, icstart+blocksize-1
      DC = EP(C)
      DO 160 B=ibstart, ibstart+blocksize-1
      DBC = EP(B) + DC
      DO 160 A=iastart, iastart+blocksize-1
      IF (A.EQ.B.AND.B.EQ.C) GOTO 160
c
c performance counters for the inner loop
c 8 adds, 5 mults, 1 divide = 14 flops
c passes: nu3-nu
c
      DABC = EP(A) + DBC
      DENOM=DIJK-DABC
      DENOM=1/DENOM
      D1=  T3(A,B,C)
      D2=  T3(A,C,B)+T3(C,B,A)+T3(B,A,C)
      D3=  T3(B,C,A)+T3(C,A,B)
      F=D1*EIGHT-FOUR*D2+D3*TWO
      X3=X3+F*D1*DENOM
 160  CONTINUE
 161  CONTINUE

      end if

      CF=ONE
      IEJ=I.EQ.J
      JEK=J.EQ.K
      IF(IEJ.OR.JEK) CF=HALF
      ETD=ETD+CF*X3
      CALL DDI_SMP_SYNC()

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
