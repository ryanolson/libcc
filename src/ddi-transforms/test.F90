program unit_tests
implicit none

call test()

end program unit_tests


#define TEST_JJ 7
#define TEST_II 3

subroutine test
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

write(6,*) 'pass'

end subroutine test
