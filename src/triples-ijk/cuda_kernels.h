#define INDX(a,b,c,ld) ( ( (c) * (ld) * (ld) ) \
                       + ( (b) * (ld) ) \
                       + ( (a) ) )

extern "C" {


__global__ void trant3_1_kernel( long int n, double *v );
__global__ void trant3_4_kernel( long int n, double *v );


} /* end extern C */
