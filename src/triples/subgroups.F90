module subgroups
implicit none

#include "mpif.h"

type subgroup
  integer :: ddi_comm, compute_comm, smp_comm
  type (subgroup), pointer :: next => null ()
end type

type (subgroup), private :: root

contains

subroutine split(color, key, sg)
implicit none
integer, intent(in) :: color, key
type (subgroup), intent(out), pointer :: sg

integer :: working_comm, compute_comm, smp_comm

call ddi_get_workingcomm(working_comm)
call ddi_get_computecomm(working_comm, compute_comm)
call ddi_get_smpcomm(working_comm, smp_comm)

allocate( sg )

call mpi_comm_split(compute_comm, color, key, sg%compute_comm)
call mpi_comm_split(smp_comm, color, key, sg%smp_comm)

end subroutine split

subroutine sync(sg)
implicit none
type (subgroup) :: sg
integer :: ierr
call mpi_barrier(sg % compute_comm, ierr)
end subroutine sync

subroutine smp_sync(sg)
implicit none
type (subgroup) :: sg
integer :: ierr
call mpi_barrier(sg % smp_comm, ierr)
end subroutine smp_sync

subroutine bcast(buffer, cnt, datatype, root, sg)
implicit none
type (subgroup) :: sg
double precision buffer(*)
integer cnt, datatype, root, comm, ierr
call mpi_bcast(buffer, cnt, datatype, root, sg % compute_comm, ierr)
end subroutine bcast

subroutine smp_bcast(buffer, cnt, datatype, root, sg)
implicit none
type (subgroup) :: sg
double precision buffer(*)
integer cnt, datatype, root, comm, ierr
call mpi_bcast(buffer, cnt, datatype, root, sg % smp_comm, ierr)
end subroutine smp_bcast

end module subgroups
