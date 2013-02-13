! =============================================================================
!   Copyright (C) 2010.  Ryan M. Olson
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
! =============================================================================

subroutine t1wt3_ijk(i,j,k,no,nu,v3,voe_ij,voe_ji,voe_ik,voe_ki,voe_jk,voe_kj,t1,eh,ep)
use common_cc, only: smp_np, smp_me, nu3, nu2, no2, one, om, zero, two, ddi_me, etd, half, eight, four
implicit none

integer,intent(in) :: i,j,k,no,nu
double precision :: v3(nu,nu,nu), t1(nu,no)
double precision :: voe_ij(nu,nu), voe_ji(nu,nu), voe_ik(nu,nu), voe_ki(nu,nu), voe_jk(nu,nu), voe_kj(nu,nu)
double precision :: eh(no), ep(nu)

double precision :: t3_ab1, t3_ab2, t3_ab3, t3_ab4, t3_ab5, t3_ab6, denom, dijk, dabc, x3, f, d1, d2, d3

integer :: a,b,c,aa,bb,cc,icntr

if(smp_np.gt.1) call smp_sync()

icntr = 0
x3 = zero
denom = one
dijk = eh(i) + eh(j) + eh(k)

do cc = 1,nu,16
   icntr = icntr+1
   if(icntr.eq.smp_np) icntr=0
   if(icntr.ne.smp_me) cycle

   do bb = 1,nu,16
   do aa = 1,nu,16

   do c = cc, min(cc+15,nu)
   do b = bb, min(bb+15,nu)
   do a = aa, min(aa+15,nu)

      if(a.gt.b) cycle
      if(a.eq.b .and. b.eq.c) cycle

      dabc    = ep(a) + ep(b) + ep(c)
      denom   = 1/(dijk-dabc)

! etd
      d1      = v3(a,b,c)
      d2      = v3(a,c,b) + v3(c,b,a) + v3(b,a,c)
      d3      = v3(b,c,a) + v3(c,a,b) 
      f       = d1*eight  - d2*four   + d3*two
      x3      = x3 + f*d1*denom

      if(a.eq.b) cycle

      d1      = v3(b,a,c)
      d2      = v3(b,c,a) + v3(c,a,b) + v3(a,b,c)
      d3      = v3(a,c,b) + v3(c,b,a)
      f       = d1*eight  - d2*four   + d3*two
      x3      = x3 + f*d1*denom

! ets
      t3_ab1  = (V3(A,B,C)-V3(B,A,C))*TWO-V3(A,C,B)+V3(B,C,A) ! ijk; abc => abc
      t3_ab2  = (v3(a,c,b)-v3(b,c,a))*two-v3(a,b,c)+v3(b,a,c) ! ikj; abc -> acb _1? ==> acb
      t3_ab3  = (v3(b,a,c)-v3(a,b,c))*two-v3(c,a,b)+v3(c,b,a) ! jik; abc -> cab _4? ==> bac
      t3_ab4  = (v3(b,c,a)-v3(a,c,b))*two-v3(c,b,a)+v3(c,a,b) ! jki; abc -> acb _1? ==> bca
      t3_ab5  = (v3(c,a,b)-v3(c,b,a))*two-v3(b,a,c)+v3(a,b,c) ! kij; abc -> bca _5? ==> cab
      t3_ab6  = (v3(c,b,a)-v3(c,a,b))*two-v3(b,c,a)+v3(a,c,b) ! kji; abc -> acb _1? ==> cba


      t1(a,i) = t1(a,i) + ( t3_ab1*voe_jk(b,c) + t3_ab2*voe_kj(b,c) )*denom
      t1(b,i) = t1(b,i) + ( t3_ab1*voe_jk(a,c) + t3_ab2*voe_kj(a,c) )*denom*om

      t1(a,j) = t1(a,j) + ( t3_ab3*voe_ik(b,c) + t3_ab5*voe_ki(b,c) )*denom
      t1(b,j) = t1(b,j) + ( t3_ab3*voe_ik(a,c) + t3_ab5*voe_ki(a,c) )*denom*om

      t1(a,k) = t1(a,k) + ( t3_ab4*voe_ij(b,c) + t3_ab6*voe_ji(b,c) )*denom
      t1(b,k) = t1(b,k) + ( t3_ab4*voe_ij(a,c) + t3_ab6*voe_ji(a,c) )*denom*om

   end do
   end do
   end do
   
   end do
   end do

end do

if(i.eq.j .or. j.eq.k) then
  etd = etd + x3*half
else
  etd = etd + x3
end if


if(smp_np.gt.1) call smp_sync()

return
9000 format(3I5,1F20.15)
end subroutine t1wt3_ijk

subroutine t1wt3_ijk_temp(i,j,k,no,nu,v3,voe_ij,voe_ji,voe_ik,voe_ki,voe_jk,voe_kj,t1,eh,ep)
use common_cc, only: smp_np, smp_me, nu3, nu2, no2, one, om, zero, two, ddi_me, etd, half, eight, four
implicit none

integer,intent(in) :: i,j,k,no,nu
double precision :: v3(nu,nu,nu), t1(nu,no)
double precision :: voe_ij(nu,nu), voe_ji(nu,nu), voe_ik(nu,nu), voe_ki(nu,nu), voe_jk(nu,nu), voe_kj(nu,nu)
double precision :: eh(no), ep(nu)

double precision :: t3_ab1, t3_ab2, t3_ab3, t3_ab4, t3_ab5, t3_ab6, denom, dijk, dabc, x3, f, d1, d2, d3, x3temp, x3accum

integer :: a,b,c,aa,bb,cc,icntr

if(smp_np.gt.1) call smp_sync()

icntr = 0
x3 = zero
denom = one
dijk = eh(i) + eh(j) + eh(k)

!   do a = 1, nu
!   do b = 1, nu
!   x3accum = 0.0
!   do c = 1, nu
!   x3temp = 0.d0
!

!      if(a.gt.b) cycle
!      if(a.eq.b .and. b.eq.c) cycle
!
!      dabc    = ep(a) + ep(b) + ep(c)
!      denom   = 1/(dijk-dabc)
!
!      if(a.eq.b) cycle

! ets
!      t3_ab1  = (V3(A,B,C)-V3(B,A,C))*TWO-V3(A,C,B)+V3(B,C,A) ! ijk; abc => abc
!      t3_ab2  = (v3(a,c,b)-v3(b,c,a))*two-v3(a,b,c)+v3(b,a,c) ! ikj; abc -> acb _1? ==> acb
!      t3_ab3  = (v3(b,a,c)-v3(a,b,c))*two-v3(c,a,b)+v3(c,b,a) ! jik; abc -> cab _4? ==> bac
!      t3_ab4  = (v3(b,c,a)-v3(a,c,b))*two-v3(c,b,a)+v3(c,a,b) ! jki; abc -> acb _1? ==> bca
!      t3_ab5  = (v3(c,a,b)-v3(c,b,a))*two-v3(b,a,c)+v3(a,b,c) ! kij; abc -> bca _5? ==> cab
!      t3_ab6  = (v3(c,b,a)-v3(c,a,b))*two-v3(b,c,a)+v3(a,c,b) ! kji; abc -> acb _1? ==> cba


!      t1(a,i) = t1(a,i) + ( t3_ab1*voe_jk(b,c) + t3_ab2*voe_kj(b,c) )*denom
!      t1(b,i) = t1(b,i) + ( t3_ab1*voe_jk(a,c) + t3_ab2*voe_kj(a,c) )*denom*om

!      t1(a,j) = t1(a,j) + ( t3_ab3*voe_ik(b,c) + t3_ab5*voe_ki(b,c) )*denom
!      t1(b,j) = t1(b,j) + ( t3_ab3*voe_ik(a,c) + t3_ab5*voe_ki(a,c) )*denom*om

!      t1(a,k) = t1(a,k) + ( t3_ab4*voe_ij(b,c) + t3_ab6*voe_ji(b,c) )*denom
!      t1(b,k) = t1(b,k) + ( t3_ab4*voe_ij(a,c) + t3_ab6*voe_ji(a,c) )*denom*om
!   end do
!   end do
!   end do

   do b = 1, nu
   do a = 1, nu
   x3accum = 0.0
   do c = 1, nu
   x3temp = 0.d0


      if(a.gt.b) cycle
      if(a.eq.b .and. b.eq.c) cycle

      dabc    = ep(a) + ep(b) + ep(c)
      denom   = 1/(dijk-dabc)

      if(a.eq.b) cycle

! ets
      t3_ab1  = (V3(A,B,C)-V3(B,A,C))*TWO-V3(A,C,B)+V3(B,C,A) ! ijk; abc => abc
      t3_ab2  = (v3(a,c,b)-v3(b,c,a))*two-v3(a,b,c)+v3(b,a,c) ! ikj; abc -> acb _1? ==> acb
      t3_ab3  = (v3(b,a,c)-v3(a,b,c))*two-v3(c,a,b)+v3(c,b,a) ! jik; abc -> cab _4? ==> bac
      t3_ab4  = (v3(b,c,a)-v3(a,c,b))*two-v3(c,b,a)+v3(c,a,b) ! jki; abc -> acb _1? ==> bca
      t3_ab5  = (v3(c,a,b)-v3(c,b,a))*two-v3(b,a,c)+v3(a,b,c) ! kij; abc -> bca _5? ==> cab
      t3_ab6  = (v3(c,b,a)-v3(c,a,b))*two-v3(b,c,a)+v3(a,c,b) ! kji; abc -> acb _1? ==> cba


!      t1(a,i) = t1(a,i) + ( t3_ab1*voe_jk(b,c) + t3_ab2*voe_kj(b,c) )*denom
!      t1(b,i) = t1(b,i) + ( t3_ab1*voe_jk(a,c) + t3_ab2*voe_kj(a,c) )*denom*om

!      t1(a,j) = t1(a,j) + ( t3_ab3*voe_ik(b,c) + t3_ab5*voe_ki(b,c) )*denom
!      t1(b,j) = t1(b,j) + ( t3_ab3*voe_ik(a,c) + t3_ab5*voe_ki(a,c) )*denom*om

!      t1(a,k) = t1(a,k) + ( t3_ab4*voe_ij(b,c) + t3_ab6*voe_ji(b,c) )*denom
!      t1(b,k) = t1(b,k) + ( t3_ab4*voe_ij(a,c) + t3_ab6*voe_ji(a,c) )*denom*om
   end do
   end do
   end do
   
   do b = 1, nu
   do a = 1, nu
   x3accum = 0.0
   do c = 1, nu
   x3temp = 0.d0


      if(a.gt.b) cycle
      if(a.eq.b .and. b.eq.c) cycle

      dabc    = ep(a) + ep(b) + ep(c)
      denom   = 1/(dijk-dabc)

! etd
      d1      = v3(a,b,c)
      d2      = v3(a,c,b) + v3(c,b,a) + v3(b,a,c)
      d3      = v3(b,c,a) + v3(c,a,b) 
      f       = d1*eight  - d2*four   + d3*two
      x3      = x3 + f*d1*denom
      x3temp  = x3temp + f*d1*denom
      x3accum = x3accum + x3temp

      if(a.eq.b) cycle

      d1      = v3(b,a,c)
      d2      = v3(b,c,a) + v3(c,a,b) + v3(a,b,c)
      d3      = v3(a,c,b) + v3(c,b,a)
      f       = d1*eight  - d2*four   + d3*two
      x3      = x3 + f*d1*denom
      x3temp  = x3temp + f*d1*denom
      x3accum = x3accum + x3temp
  enddo
  enddo
  enddo

!do b = 1, no
!do a = 1, nu
! write(6,*)'a,b, t1',a,b,t1(a,b)
!enddo
!enddo

!write(6,*)'Fortran etd x3 ', etd, x3
! if(i.eq.j .or. j.eq.k) then
!   etd = etd + x3*half
! else
!   etd = etd + x3
! end if

if(smp_np.gt.1) call smp_sync()

return
9000 format(3I5,1F20.15)
end subroutine t1wt3_ijk_temp
