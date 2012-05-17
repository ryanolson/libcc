#define DEBUG_ASSERT(assertion)

module cc_transformations
implicit none
contains

subroutine translate_2d_to_4d(ii,jj,nr,nc,i,j,k,l,n1,n2,n3,n4)
implicit none
integer :: ii, jj, nr, nc, n1, n2, n3, n4
integer :: i, j, k, l
call translate_4d_from_2d(ii,jj,nr,nc,i,j,k,l,n1,n2,n3,n4)
end subroutine translate_2d_to_4d

subroutine translate_4d_from_2d(ii,jj,nr,nc,i,j,k,l,n1,n2,n3,n4)
implicit none

integer, intent(in)  :: ii, jj, nr, nc, n1, n2, n3, n4
integer, intent(out) :: i, j, k, l

integer :: iijj, n12, n123

DEBUG_ASSERT(nr*nc .eq. n1*n2*n3*n4)
DEBUG_ASSERT(ii.gt.0 .and. ii.le.nr)
DEBUG_ASSERT(jj.gt.0 .and. jj.le.nc)

n12  = n1*n2
n123 = n12*n3

iijj = nr*(jj-1) + ii

l = (iijj-1) / n123

iijj = iijj - n123*l
k = (iijj-1) / n12

iijj = iijj - n12*k
j = (iijj-1) / n1

iijj = iijj - n1*j
i = iijj

l = l + 1
k = k + 1
j = j + 1
i = i

DEBUG_ASSERT(i.gt.0 .and. i.le.n1)
DEBUG_ASSERT(j.gt.0 .and. j.le.n2)
DEBUG_ASSERT(k.gt.0 .and. k.le.n3)
DEBUG_ASSERT(l.gt.0 .and. l.le.n4)
end subroutine translate_4d_from_2d


subroutine translate_4d_to_2d(i,j,k,l,n1,n2,n3,n4,ii,jj,nr,nc)
implicit none
integer :: i, j, k, l, nr, nc, n1, n2, n3, n4
integer :: ii, jj
call translate_2d_from_4d(i,j,k,l,n1,n2,n3,n4,ii,jj,nr,nc)
subroutine translate_4d_to_2d

subroutine translate_2d_from_4d(i,j,k,l,n1,n2,n3,n4,ii,jj,nr,nc)
implicit none

integer, intent(in)  :: i, j, k, l, nr, nc, n1, n2, n3, n4
integer, intent(out) :: ii, jj

integer :: ijkl, n12, n123

DEBUG_ASSERT(nr*nc .eq. n1*n2*n3*n4)
DEBUG_ASSERT(i.gt.0 .and. i.le.n1)
DEBUG_ASSERT(j.gt.0 .and. j.le.n2)
DEBUG_ASSERT(k.gt.0 .and. k.le.n3)
DEBUG_ASSERT(l.gt.0 .and. l.le.n4)

n12  = n1*n2
n123 = n12*n3

ijkl = n123*(l-1) + n12*(k-1) + n1*(j-1) + i

jj = (ijkl-1) / nr

ijkl = ijkl - nr*jj
ii = ijkl

jj = jj+1

DEBUG_ASSERT(ii.gt.0 .and. ii.le.nr)
DEBUG_ASSERT(jj.gt.0 .and. jj.le.nc)
end subroutine translate_2d_from_4d



subroutine increment_4d(i,j,k,l,n1,n2,n3,n4,incr)
implicit none

integer :: i, j, k, l
integer, intent(in) :: n1, n2, n3, n4, incr

DEBUG_ASSERT(incr .gt. 0)
DEBUG_ASSERT(i.gt.0 .and. i.le.n1)
DEBUG_ASSERT(j.gt.0 .and. j.le.n2)
DEBUG_ASSERT(k.gt.0 .and. k.le.n3)
DEBUG_ASSERT(l.gt.0 .and. l.le.n4)

i = i + incr

do while(i .gt. n1)
  j = j + 1
  i = i - n1
end do

do while(j .gt. n2)
  k = k + 1
  j = j - n2
end do

do while(k .gt. n3)
  l = l + 1
  k = k - n3
end do

DEBUG_ASSERT(i.gt.0 .and. i.le.n1)
DEBUG_ASSERT(j.gt.0 .and. j.le.n2)
DEBUG_ASSERT(k.gt.0 .and. k.le.n3)
DEBUG_ASSERT(l.gt.0 .and. l.le.n4)
end subroutine increment_4d


#if 0

subroutine tranmd_12_using_puts(d_in, n1, n2, n3, n4, d_out)
implicit none

call ddi_dimensions(d_in, nrows_in, ncols_in)
call ddi_dimensions(d_out, nrows_out, ncols_out)

call ddi_distrib(d_in, ddi_me, ilo, ihi, jlo, jhi)

call ddi_sync_on_array(d_in)

do jj = jlo,jhi
do ii = ilo,ihi

   call translate_4d_from_2d(ii,jj,nrows_in,ncols_in, &
                             i,j,k,l,n1,n2,n3,n4)

   call translate_2d_from_4d(j,i,k,l,n2,n1,n3,n4, &
                             kk,ll,nrows_out,ncols_out)

   call ddi_nbput(d_out, kk, kk, ll, ll, local(ii,jj))

end do
end do

call ddi_get_communicator_for_array(d_out, comm)

call ddi_sync(comm)

end subroutine tranmd_12_using_puts


subroutine tranmd_12_using_gets(d_src, n1, n2, n3, n4, d_dst)
implicit none

integer :: nrows_src, ncols_src
integer :: nrows_dst, ncols_dst

call ddi_dimensions(d_src, nrows_src, ncols_src)
call ddi_dimensions(d_dst, nrows_dst, ncols_dst)

call ddi_distrib(d_out, ddi_me, ilo, ihi, jlo, jhi)

call ddi_sync_on_array(d_src)

do jj = jlo, jhi
do ii = ilo, ihi

 ! for each element of the local patch of d_dst
 ! determine what the [i,j,k,l] coordinates are
   call translate_4d_from_2d(ii, jj, nrows_dst, ncols_dst, &
                             i, j, k, l, n1, n2, n3, n4)
   
   call transform_4d_by_trans12(i,j,k,l)
 ! prior to the transformation [i,j,k,l] came from [j,i,k,l]
   call translate_2d_from_4d(j, i, k, l, n1, n2, n3, n4, &
                             kk, ll, nrows_src, ncols_src)

   call ddi_get(d_src, kk, kk, ll, ll, local(ii,jj)

end subroutine tranmd_12_using_gets

#endif


subroutine swap(i,j)
implicit none
integer :: i,j,tmp
tmp = i
i = j
j = tmp
end subroutine swap

subroutine transform_4d_by_tran12(i,j,k,l)
implicit none
integer :: i,j,k,l,tmp
call swap(i,j)
end subroutine transform_4d_by_tran12

subroutine transform_4d_by_tran13(i,j,k,l)
implicit none
integer :: i,j,k,l,tmp
call swap(i,k)
end subroutine transform_4d_by_tran13

subroutine transform_4d_by_tran14(i,j,k,l)
implicit none
integer :: i,j,k,l,tmp
call swap(i,l)
end subroutine transform_4d_by_tran14

subroutine transform_4d_by_tran23(i,j,k,l)
implicit none
integer :: i,j,k,l,tmp
call swap(j,k)
end subroutine transform_4d_by_tran23

subroutine transform_4d_by_tran24(i,j,k,l)
implicit none
integer :: i,j,k,l,tmp
call swap(j,l)
end subroutine transform_4d_by_tran24

subroutine transform_4d_by_tran34(i,j,k,l)
implicit none
integer :: i,j,k,l,tmp
call swap(k,l)
end subroutine transform_4d_by_tran34


subroutine transform_4d_by_swap12(i,j,k,l,n1,n2,n3,n4)
implicit none
integer :: i,j,k,l,n1,n2,n3,n4,tmp
call swap(i,j)
call swap(n1,n2)
end subroutine transform_4d_by_swap12

subroutine transform_4d_by_swap13(i,j,k,l,n1,n2,n3,n4)
implicit none
integer :: i,j,k,l,n1,n2,n3,n4,tmp
call swap(i,k)
call swap(n1,n3)
end subroutine transform_4d_by_swap13

subroutine transform_4d_by_swap23(i,j,k,l,n1,n2,n3,n4)
implicit none
integer :: i,j,k,l,n1,n2,n3,n4,tmp
call swap(j,k)
call swap(n2,n3)
end subroutine transform_4d_by_swap23

end module cc_transformations
