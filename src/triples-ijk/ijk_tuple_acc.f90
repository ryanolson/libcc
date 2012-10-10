
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
! FLOPS Analysis
!  2 Types of DGEMMs
!  - Type A: nr,nu,no contracting t2 and vm ==> 2*nr*nou flops per dgemm
!  - Type B: nu,nr,nu contracting t2 and ve ==> 2*nr*nu2 flops per dgemm
!
!  There are 6 DGEMMs of each type:
!  - FLOPS per RANK = 6 * 2 * nr * (nou + nu2)
!                   = 12 * nr  * (nou + nu2)
!
!  - FLOPS per NODE = 12 * nu2 * (nou + nu2)
!
! MEMORY Analysis
!  All arrays that are passed in are "shared" arrays. That is these arrays
!  are in shared memory and can be accessed by ANY PE/rank on the node.
!
subroutine ijk_tuple(i,j,k,t2_i,t2_j,t2_k,vm_ij,vm_ji,vm_ik,vm_ki,vm_jk,vm_kj, &
                     ve_i,ve_j,ve_k,v3)
use common_cc
implicit none

integer,intent(in) :: i, j, k
double precision,intent(in) :: t2_i(nu2,no), t2_j(nu2,no), t2_k(nu2,no)
double precision,intent(in) :: vm_ij(nou), vm_ji(nou), vm_ik(nou), vm_ki(nou), vm_jk(nou), vm_kj(nou)
double precision,intent(in) :: ve_i(nu3), ve_j(nu3), ve_k(nu3)
double precision,intent(out) :: v3(nu3)

integer nr

nr = nu2

! note: the syncs between dgemms are required because each dgemm is updating a 
! different disjoint portion of v3.

!$acc data present(v3,t2_i,t2_j,t2_k,ve_i,ve_j,ve_k,vm_ij,vm_ji,vm_ik,vm_ki,vm_jk,vm_kj)


!$acc host_data use_device(t2_i,t2_j,t2_k,ve_i,ve_j,ve_k,vm_ij,vm_ji,vm_ik,vm_ki,vm_jk,vm_kj)
! #2 & #9
!$acc host_data use_device(v3)
call dgemm('n','n',nr,nu,no,om,t2_j,nu2,vm_ki,no,zero,v3,nu2)       ! #2: Type A

!$acc wait  ! for ve_j

call dgemm('t','t',nr,nu,nu,one,ve_j,nu,t2_k(1,i),nu,one,v3,nu2)    ! #9: Type BT
!$acc end host_data

! transform v3
call trant3_1(nu,v3) 

! #4 & #7
!$acc host_data use_device(v3)
call dgemm('t','t',nu,nr,no,om,vm_ji,no,t2_k,nu2,one,v3,nu)         ! #4: Type AT
call dgemm('n','n',nu,nr,nu,one,t2_j(1,i),nu,ve_k,nu,one,v3,nu)     ! #7: Type B
!$acc end host_data

! transform v3
call trant3_4(nu,v3)

! #0 & #11
!$acc host_data use_device(v3)
call dgemm('n','n',nr,nu,no,om,t2_i,nu2,vm_kj,no,one,v3,nu2)        ! #0: Type A
call dgemm('t','t',nr,nu,nu,one,ve_i,nu,t2_k(1,j),nu,one,v3,nu2)    ! #11: Type BT
! #8 & #3
call dgemm('n','n',nu,nr,nu,one,t2_i(1,k),nu,ve_j,nu,one,v3,nu)     ! #8: Type B
call dgemm('t','t',nu,nr,no,om,vm_ik,no,t2_j,nu2,one,v3,nu)         ! #3: Type AT
!$acc end host_data

! transform v3
call trant3_1(nu,v3)

! #1 & #10
!$acc host_data use_device(v3)
call dgemm('n','n',nr,nu,no,om,t2_i,nu2,vm_jk,no,one,v3,nu2)        ! #1: Type A
call dgemm('t','t',nr,nu,nu,one,ve_i,nu,t2_j(1,k),nu,one,v3,nu2)    ! #10: Type BT
! #6 & #5
call dgemm('n','n',nu,nr,nu,one,t2_i(1,j),nu,ve_k,nu,one,v3,nu)     ! #6: Type B
call dgemm('t','t',nu,nr,no,om,vm_ij,no,t2_k,nu2,one,v3,nu)         ! #5: Type AT
!$acc end host_data

! transform v3
call trant3_1(nu,v3)

!$acc end host_data
!$acc end data

end subroutine ijk_tuple

!-------------------------------------

  subroutine trant3_1_async(acc_sync,n,v)
  implicit none
  integer :: n, a, b, c, acc_sync
  double precision :: v(n,n,n), x

!$acc parallel loop private(x) async(acc_sync)
  DO B=1,N
  DO C=1,B-1
  DO A=1,N
     X=V(A,B,C)
     V(A,B,C)=V(A,C,B)
     V(A,C,B)=X
  end do
  end do
  end do
!$acc end parallel loop

  end subroutine trant3_1_async
