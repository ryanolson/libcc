#ifndef HAS_VE_EXPANSION
#define HAS_VE_EXPANSION 0
#endif

#define INDX(a,b,c,ld) ( ( (c) * (ld) * (ld) ) \
                       + ( (b) * (ld) ) \
                       + ( (a) ) )

#define SHARED_REDUCTION_SIZE 128

#ifndef HAVE_VE_EXPANSION_KERNEL
#define HAVE_VE_EXPANSION_KERNEL 1
#endif

extern "C" {

__global__ void exp_trsq_kernel( long int n, double *v, double *a );

__global__ void expand_tr_kernel( long int n, double *v );

__global__ void expand_trsq_kernel( long int n, double *v );

__global__ void trant3_1_kernel( long int n, double *v );

__global__ void trant3_4_kernel( long int n, double *v );

__global__ void etd_cuda_kernel( const int i, const int j, const int k,
        const int no, const int nu,
        const double *v3, const double *voe_ij, const double *voe_ji,
        const double *voe_ik, const double *voe_ki,
        const double *voe_jk, const double *voe_kj,
        double *t1, const double *eh, const double *ep, double *etd_reduce );

__global__ void t1a_cuda_kernel( const int i, const int j, const int k,
        const int no, const int nu,
        const double *v3, const double *voe_ij, const double *voe_ji,
        const double *voe_ik, const double *voe_ki,
        const double *voe_jk, const double *voe_kj,
        double *t1, const double *eh, const double *ep, double *etd_reduce );

__global__ void reduce_etd_kernel( const long int, const double *, double * );


} /* end extern C */
