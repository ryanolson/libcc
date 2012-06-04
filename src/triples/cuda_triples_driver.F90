! =====================================================================
! Notes:
!
! Prior to entering the triples_cuda_driver, the application has
! formed disjoint communicators for the gpu drivers and the remaining
! cpu ranks
! =====================================================================

subroutine triples_cuda_driver(no, nu, sr, nr, vei, vej, vek, d_vvvo)
implicit none

integer :: no, nu, sr, nr, d_vvvo
double precision :: vei, vej, vek
integer :: ijk, i, j, k, ilo, ihi
integer :: iold, jold, kold, nutr

iold = -1
jold = -1
kold = -1
nutr = (nu*nu+nu)/2

do ijk = sr, sr+nr-1
   
   call ddcc_t_task(ijk,no,i,j,k)

   if(i.ne.iold) then
     ilo = nu*(i-1) + 1
     ihi = ilo + nu
     call ddi_get(d_vvvo,1,nutr,ilo,ihi,vei)
   end if
  
   if(j.ne.jold) then
     ilo = nu*(j-1) + 1
     ihi = ilo + nu
     call ddi_get(d_vvvo,1,nutr,ilo,ihi,vej)
   end if

   if(k.ne.kold) then
     ilo = nu*(k-1) + 1
     ihi = ilo + nu
     call ddi_get(d_vvvo,1,nutr,ilo,ihi,vek)
   end if

   call ijk_gpu_driver(nu,no,i,j,k,vei,vej,vek)

   iold = i
   jold = j
   kold = k

end do

end subroutine triples_cuda_driver
