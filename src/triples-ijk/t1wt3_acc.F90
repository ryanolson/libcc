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

#define ETD_X3_ON_DEVICE 1

subroutine t1wt3_ijk(i,j,k,no,nu,v3,voe_ij,voe_ji,voe_ik,voe_ki,voe_jk,voe_kj,t1,eh,ep)
use common_cc, only: nu3, nu2, no2, one, om, zero, two, ddi_me, etd, half, eight, four, x3
implicit none

integer,intent(in) :: i,j,k,no,nu
double precision :: v3(nu,nu,nu), t1(nu,no)
double precision :: voe_ij(nu,nu), voe_ji(nu,nu), voe_ik(nu,nu), voe_ki(nu,nu), voe_jk(nu,nu), voe_kj(nu,nu)
double precision :: eh(no), ep(nu)

double precision :: t3_ab1, t3_ab2, t3_ab3, t3_ab4, t3_ab5, t3_ab6, denom, dijk, dabc, f, d1, d2, d3
double precision :: t1ai, t1bi, t1aj, t1bj, t1ak, t1bk

integer :: a,b,c,aa,bb,cc,icntr

x3 = zero
denom = one
dijk = eh(i) + eh(j) + eh(k)

#if ETD_X3_ON_DEVICE
!$acc data present( etd, x3 )
!$acc kernels present( x3 ) async(5)
  x3 = zero
!$acc end kernels
#endif

!$acc parallel loop private(t1ai,t1aj,t1ak) reduction(+:x3) private(dabc,denom,d1,d2,d3,f,t3_ab1,t3_ab2,t3_ab3,t3_ab4,t3_ab5,t3_ab6) async(5)
   do a = 1, nu
    t1ai = t1(a,i)
    t1aj = t1(a,j)
    t1ak = t1(a,k)
   do b = 1, nu
   do c = 1, nu

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
      t3_ab1  = (v3(a,b,c)-v3(b,a,c))*two-v3(a,c,b)+v3(b,c,a) ! ijk; abc => abc
      t3_ab2  = (v3(a,c,b)-v3(b,c,a))*two-v3(a,b,c)+v3(b,a,c) ! ikj; abc -> acb _1? ==> acb
      t3_ab3  = (v3(b,a,c)-v3(a,b,c))*two-v3(c,a,b)+v3(c,b,a) ! jik; abc -> cab _4? ==> bac
      t3_ab4  = (v3(b,c,a)-v3(a,c,b))*two-v3(c,b,a)+v3(c,a,b) ! jki; abc -> acb _1? ==> bca
      t3_ab5  = (v3(c,a,b)-v3(c,b,a))*two-v3(b,a,c)+v3(a,b,c) ! kij; abc -> bca _5? ==> cab
      t3_ab6  = (v3(c,b,a)-v3(c,a,b))*two-v3(b,c,a)+v3(a,c,b) ! kji; abc -> acb _1? ==> cba

      t1ai = t1ai + ( t3_ab1*voe_jk(b,c) + t3_ab2*voe_kj(b,c) )*denom
      t1aj = t1aj + ( t3_ab3*voe_ik(b,c) + t3_ab5*voe_ki(b,c) )*denom
      t1ak = t1ak + ( t3_ab4*voe_ij(b,c) + t3_ab6*voe_ji(b,c) )*denom

   end do
   end do
      t1(a,i) = t1ai
      t1(a,j) = t1aj
      t1(a,k) = t1ak
   end do
!$acc end parallel loop

!$acc parallel loop private(t1bi,t1bj,t1bk) private(dabc,denom,t3_ab1,t3_ab2,t3_ab3,t3_ab4,t3_ab5,t3_ab6) async(5)
   do b = 1, nu
    t1bi = t1(b,i)
    t1bj = t1(b,j)
    t1bk = t1(b,k)
   do c = 1, nu

   do a = 1, nu

      if(a.gt.b) cycle
      if(a.eq.b .and. b.eq.c) cycle

      dabc    = ep(a) + ep(b) + ep(c)
      denom   = 1/(dijk-dabc)

      if(a.eq.b) cycle

! ets
      t3_ab1  = (v3(a,b,c)-v3(b,a,c))*two-v3(a,c,b)+v3(b,c,a) ! ijk; abc => abc
      t3_ab2  = (v3(a,c,b)-v3(b,c,a))*two-v3(a,b,c)+v3(b,a,c) ! ikj; abc -> acb _1? ==> acb
      t3_ab3  = (v3(b,a,c)-v3(a,b,c))*two-v3(c,a,b)+v3(c,b,a) ! jik; abc -> cab _4? ==> bac
      t3_ab4  = (v3(b,c,a)-v3(a,c,b))*two-v3(c,b,a)+v3(c,a,b) ! jki; abc -> acb _1? ==> bca
      t3_ab5  = (v3(c,a,b)-v3(c,b,a))*two-v3(b,a,c)+v3(a,b,c) ! kij; abc -> bca _5? ==> cab
      t3_ab6  = (v3(c,b,a)-v3(c,a,b))*two-v3(b,c,a)+v3(a,c,b) ! kji; abc -> acb _1? ==> cba

      t1bi = t1bi + ( t3_ab1*voe_jk(a,c) + t3_ab2*voe_kj(a,c) )*denom*om
      t1bj = t1bj + ( t3_ab3*voe_ik(a,c) + t3_ab5*voe_ki(a,c) )*denom*om
      t1bk = t1bk + ( t3_ab4*voe_ij(a,c) + t3_ab6*voe_ji(a,c) )*denom*om

   end do
   end do
      t1(b,i) = t1bi
      t1(b,j) = t1bj
      t1(b,k) = t1bk
   end do
!$acc end parallel loop

#if ETD_X3_ON_DEVICE
!$acc kernels present( etd, x3) async(5)
  etd = etd + x3
!$acc end kernels
!$acc update host( etd ) async(5)
!$acc end data
#else
!$acc wait(5)
etd = etd + x3
#endif

return
9000 format(3I5,1F20.15)
end subroutine t1wt3_ijk

