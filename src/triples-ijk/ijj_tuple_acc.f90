
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

subroutine ijj_tuple_formv3(i,j,t2_i,t2_j,vm_ij,vm_ji,vm_jj,ve_i,ve_j,v3)
use common_cc
implicit none

integer :: i, j
double precision :: t2_i(nu2,no), t2_j(nu2,no)
double precision :: vm_ij(nou), vm_ji(nou), vm_jj(nou)
double precision :: ve_i(nu3), ve_j(nu3)
double precision :: v3(nu3)

integer :: nr,sr,t3off

call div_even(nu2,smp_np,smp_me,nr,sr)
      T3OFF = (SR-1)*NU + 1

if(smp_np.gt.1) call smp_sync()

!$acc data present(t2_i,t2_j,vm_ij,vm_ji,vm_jj,ve_i,ve_j,v3)

!$acc host_data use_device(t2_i,t2_j,vm_ij,vm_ji,vm_jj,ve_i,ve_j)

! starting arrangement of v3
! #2, #3 - 12345678 == T(13572468)

!$acc host_data use_device(v3)
call dgemm('t','t',nu,nr,no,om,vm_ij,no,t2_j(sr,1),nu2,zero,v3(t3off),nu)           ! #2: Type AT

!$acc wait  ! for ve_j

call dgemm('n','n',nu,nr,nu,one,t2_i(1,j),nu,ve_j(t3off),nu,one,v3(t3off),nu)       ! #3: Type B
!$acc end host_data 

! transform 12345678 -> 15372648
call trant3_acc(v3,nu,3)

! #4, #1 - 15372648 == T(13245758)
!$acc host_data use_device(v3)
call dgemm('n','n',nu,nr,nu,one,t2_j(1,i),nu,ve_j(t3off),nu,one,v3(t3off),nu)       ! #4: Type B
call dgemm('t','t',nu,nr,no,om,vm_ji,no,t2_j(sr,1),nu2,one,v3(t3off),nu)            ! #1: Type AT
!$acc end host_data 

! transform 15372648 -> 15263748
call trant3_1(nu,v3)

! #5, #0 - 15263748 == T(12345678)
!$acc host_data use_device(v3)
call dgemm('n','n',nu,nr,nu,one,t2_j(1,j),nu,ve_i(t3off),nu,one,v3(t3off),nu)       ! #5: Type B
call dgemm('t','t',nu,nr,no,om,vm_jj,no,t2_i(sr,1),nu2,one,v3(t3off),nu)            ! #0: Type AT
!$acc end host_data 

! transform 15263748 -> 12563478
call trant3_acc(v3,nu,2)

! ending arrangement of v3
! #6 - 12563478

!$acc end host_data 
!$acc end data

if(smp_np.gt.1) call smp_sync()

end subroutine ijj_tuple_formv3
