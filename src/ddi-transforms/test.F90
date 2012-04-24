program unit_tests
implicit none

call test_4d_from_2d()
call test_2d_from_4d()

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
