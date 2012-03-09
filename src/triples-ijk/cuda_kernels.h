#define INDX(a,b,c,ld) ( ( (c) * (ld) * (ld) ) \
                       + ( (b) * (ld) ) \
                       + ( (a) ) )

extern "C" {


__global__ void trant3_1_kernel( long int n, double *v );
__global__ void trant3_4_kernel( long int n, double *v );
__global__ void t1wt3_cuda_kernel( const int i, const int j, const int k,
        const int no, const int nu,
        const double *v3, const double *voe_ij, const double *voe_ji,
        const double *voe_ik, const double *voe_ki,
        const double *voe_jk, const double *voe_kj,
        double *t1, const double *eh, const double *ep, double *etd_reduce );

__global__ void reduce_etd_kernel( const long int, const double *, double * );


} /* end extern C */
