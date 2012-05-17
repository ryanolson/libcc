#define TEST_NO  5
#define TEST_NU 11

program unit_tests
implicit none

call test_4d_from_2d()
call test_2d_from_4d()
call test_increment_4d()
call test_swap()

end program unit_tests


#define TEST_II 3
#define TEST_JJ 7
#define TEST_KK 7
#define TEST_LL 3

subroutine test_4d_from_2d
use cc_transformations
implicit none

integer :: ii, jj
integer :: i, j, k, l

do jj = 1,TEST_JJ
do ii = 1,TEST_II
   call translate_4d_from_2d(ii,jj,TEST_II,TEST_JJ,i,j,k,l,TEST_II,TEST_JJ,1,1)

   if(ii.eq.i .and. jj.eq.j) cycle

   write(6,911) ii,jj,i,j
911 format('failure : ',4I4)
   stop

end do
end do

write(6,*) 'pass: test_4d_from_2d'

end subroutine test_4d_from_2d

subroutine test_2d_from_4d
use cc_transformations
implicit none

integer :: ii, jj
integer :: i, j, k, l
integer :: it, jt, kt, lt
integer :: nr, nc

nr = TEST_II*TEST_JJ
nc = TEST_KK*TEST_LL

do l = 1, TEST_LL
do k = 1, TEST_KK
do j = 1, TEST_JJ
do i = 1, TEST_II

   call translate_2d_from_4d(i,j,k,l,TEST_II,TEST_JJ,TEST_KK,TEST_LL,ii,jj,nr,nc)
   call translate_4d_from_2d(ii,jj,nr,nc,it,jt,kt,lt,TEST_II,TEST_JJ,TEST_KK,TEST_LL)

   if(i.eq.it .and. j.eq.jt .and. k.eq.kt .and. l.eq.lt) cycle

   write(6,911) i,j,k,l,ii,jj
911 format('failure : ',6I4)

end do
end do
end do
end do

write(6,*) 'pass: test_2d_from_4d'

end subroutine test_2d_from_4d


subroutine test_increment_4d()
use cc_transformations
implicit none

integer :: i, j, k, l
integer :: it, jt, kt, lt

it = 1
jt = 1
kt = 1
lt = 1

do l = 1,TEST_LL-1
do k = 1,TEST_KK
do j = 1,TEST_JJ
do i = 1,TEST_II
   if(i.ne.it .or. j.ne.jt .or. k.ne.kt .or. l.ne.lt) then
      write(6,9000) i,j,k,l,it,jt,kt,lt
 9000 format(' fail : ',8I4)
      return
   endif
   call increment_4d(it,jt,kt,lt,TEST_II,TEST_JJ,TEST_KK,TEST_LL,1)
end do
end do
end do
end do

write(6,*) 'pass: increment_4d'
end subroutine test_increment_4d

subroutine test_swap
use cc_transformations
implicit none
integer :: i, j

i=32
j=46
call swap(i,j)
if(i.eq.46 .and. j.eq.32) then
   write(6,*) 'pass: swap'
   return
endif

write(6,*) 'fail: swap'
end subroutine test_swap


subroutine test_indexed_swap12
use cc_transformations
implicit none

integer :: no, nu
integer :: tmp(TEST_NU,TEST_NU,TEST_NU)
integer :: a_src(TEST_NO,TEST_NU,TEST_NU,TEST_NO)
integer :: a_dst(TEST_NU,TEST_NO,TEST_NU,TEST_NO)
integer :: cntr, length
integer :: i, a, b, j
integer :: i1, i2, i3, i4
integer :: n1, n2, n3, n4

no = TEST_NO
nu = TEST_NU
cntr = 0
length = no*no*nu*nu

do j = 1, no
do b = 1, nu
do a = 1, nu
do i = 1, no
   cntr = cntr+1
   a_src(i,a,b,j) = cntr
end do
end do
end do
end do

! copy src ==> dst
call dcopy(length, a_src, 1, a_dst, 1)

! reorder dst
call insi12(no,nu,nu,tmp,a_dst)

do j = 1, no
do b = 1, nu
do i = 1, no
do a = 1, nu
   cntr = a_dst(a,i,b,j)

   i1 = a
   n1 = no

   i2 = i
   n2 = no

   i3 = b
   n3 = nu

   i4 = j
   n4 = no

   call transform_4d_by_swap12(i1, i2, i3, i4, n1, n2, n3, n4)

   if(a_src(i1,i2,i3,i4).ne.cntr) then
      write(6,9000) a, i, b, j, i1, i2, i3, i4
 9000 format(' fail : ',8I4)
      stop
   endif

end do
end do
end do
end do

end subroutine test_indexed_swap12
