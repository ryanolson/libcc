
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
!
subroutine iij_tuple_formv3(i,j,t2_i,t2_j,vm_ij,vm_ji,vm_ii,ve_i,ve_j,v3)
use common_cc
implicit none

integer :: i, j
double precision :: t2_i(nu2,no), t2_j(nu2,no)
double precision :: vm_ij(nou), vm_ji(nou), vm_ii(nou)
double precision :: ve_i(nu3), ve_j(nu3)
double precision :: v3(nu3)

integer :: nr,sr,t3off

call div_even(nu2,smp_np,smp_me,nr,sr)
      T3OFF = (SR-1)*NU + 1

if(smp_np.gt.1)  call smp_sync()

! #1, #3 - 12563478
call dgemm('n','n',nr,nu,no,om,t2_i(sr,1),nu2,vm_ij,no,zero,v3(sr),nu2)              ! #1: Type A
if(smp_np.gt.1) call smp_sync()
call dgemm('n','n',nu,nr,nu,one,t2_i(1,i),nu,ve_j(t3off),nu,one,v3(t3off),nu)        ! #3: Type B

! transform
call trant3_acc(v3,nu,2)

! #2 - 15263748 
call dgemm('n','n',nr,nu,no,om,t2_j(sr,1),nu2,vm_ii,no,one,v3(sr),nu2)               ! #2: Type A

! transform
call trant3_1(nu,v3)

! #5 - 15372648
call dgemm('n','n',nu,nr,nu,one,t2_j(1,i),nu,ve_i(t3off),nu,one,v3(t3off),nu)        ! #5: Type B

! transform
call trant3_acc(v3,nu,3)

! #0, #4 - 12345678
call dgemm('n','n',nr,nu,no,om,t2_i(sr,1),nu2,vm_ji,no,one,v3(sr),nu2)               ! #0: Type A
if(smp_np.gt.1) call smp_sync()
call dgemm('n','n',nu,nr,nu,one,t2_i(1,j),nu,ve_i(t3off),nu,one,v3(t3off),nu)        ! #4: Type B

if(smp_np.gt.1) call smp_sync()
end subroutine iij_tuple_formv3
