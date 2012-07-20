subroutine triples_cuda_init( no, nu, eh, ep, v1, t2, vm, voe )
implicit none
integer :: no, nu
double precision :: eh(*), ep(*), v1(*), t2(*), vm(*), voe(*)
return
end subroutine

subroutine triples_cuda_driver(no, nu, sr, nr, vei, vej, vek, d_vvvo)
implicit none
integer :: no, nu, sr, nr, d_vvvo
double precision :: vei(*), vej(*), vek(*)
return
end subroutine

subroutine triples_cuda_finalize( no, nu, etd )
implicit none
integer :: no, nu
double precision :: etd
return
end subroutine
