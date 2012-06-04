#include <stdio.h>
#include "cuda_runtime.h"
#include "cuda_kernels.h"

#define INDX(a,b,c,ld) ( ( (c) * (ld) * (ld) ) \
                       + ( (b) * (ld) ) \
                       + ( (a) ) )

extern "C" {

__device__ void warpReduce( volatile double *sdata )
{
#if 1
  if( blockDim.x >= 64 && blockDim.x % 32 == 0 )
  {
    for( unsigned int s = blockDim.x / 2; s > 32; s >>= 1)
    {
      if( threadIdx.x < s ) 
         sdata[threadIdx.x] += sdata[threadIdx.x + s];
      __syncthreads();
    } /* end for */
    if( threadIdx.x < 32 ) 
    {
      sdata[threadIdx.x] += sdata[threadIdx.x + 32];
      sdata[threadIdx.x] += sdata[threadIdx.x + 16];
      sdata[threadIdx.x] += sdata[threadIdx.x +  8];
      sdata[threadIdx.x] += sdata[threadIdx.x +  4];
      sdata[threadIdx.x] += sdata[threadIdx.x +  2];
      sdata[threadIdx.x] += sdata[threadIdx.x +  1];
    } /* end if */
  }
  else
#endif
  {
    if( threadIdx.x == 0 )
    {
      for( int idx = 1; idx < blockDim.x; idx++ )
      {
        sdata[0] += sdata[idx];
      } /* end for */
    } /* end if */
  } /* end */
} /* end warp Reduce */

__global__ void etd_cuda_kernel( const int i, const int j, const int k, 
        const int no, const int nu,
        const double *v3, const double *voe_ij, const double *voe_ji, 
        const double *voe_ik, const double *voe_ki, 
        const double *voe_jk, const double *voe_kj,
        double *t1, const double *eh, const double *ep, double *etd_reduce )
{

__shared__ double etd_shared[SHARED_REDUCTION_SIZE];

  int a = blockIdx.x;
  int b = blockIdx.y;
  double dijk = eh[i-1] + eh[j-1] + eh[k-1];
  double x3 = 0.0;
  const double two = 2.0, four = 4.0, eight = 8.0, om = -1.0;
  double d1,d2,d3,f, t1ai = 0.0, t1bi = 0.0;
  double t1aj = 0.0, t1bj = 0.0, t1ak = 0.0, t1bk = 0.0;

  for( int idx = 0; idx < nu; idx += blockDim.x )
  {
    int c = idx + threadIdx.x;

/*
 * don't do the loop if my id is outside the bounds of nu
 */

    if( c < nu )
    {

      if( a > b ) goto loop_end;
      if( a == b && b == c ) goto loop_end;
      double dabc = ep[a] + ep[b] + ep[c];
      double denom = 1.0 / ( dijk - dabc );

      d1 = v3[INDX(a, b, c, nu)];
      d2 = v3[INDX(a, c, b, nu)] + v3[INDX(c, b, a, nu)] 
                + v3[INDX(b, a, c, nu)];
      d3 = v3[INDX(b, c, a, nu)] + v3[INDX(c, a, b, nu)];
      f  = d1*eight - d2*four + d3*two;
      x3        += f*d1*denom;

      if( a == b ) goto loop_end;

      d1 = v3[INDX(b, a, c, nu)];
      d2 = v3[INDX(b, c, a, nu)] + v3[INDX(c, a, b, nu)] 
         + v3[INDX(a, b, c, nu)];
      d3 = v3[INDX(a, c, b, nu)] + v3[INDX(c, b, a, nu)];
      f  = d1*eight - d2*four + d3*two;
      x3 += f*d1*denom;

    } /* end if */
    
loop_end:
  } /* end for */


  etd_shared[threadIdx.x] = 0.0;
  etd_shared[threadIdx.x] = x3;

  __syncthreads();

  int offset = INDX(a, b, 0, gridDim.x );

  double temp = 0.0;

#if 0
  if( threadIdx.x == 0 )
  {
    for( int idx = 0; idx < blockDim.x; idx++ )
    {
      temp += etd_shared[idx];
    } /* end for */
    etd_reduce[offset] = temp;
  } /* end if */
#endif
#if 0
  for( unsigned int s = blockDim.x / 2; s > 32; s >>= 1)
  {
    if( threadIdx.x < s ) 
       etd_shared[threadIdx.x] += etd_shared[threadIdx.x + s];
    __syncthreads();
  } /* end for */

  if( threadIdx.x < 32 ) warpReduce( etd_shared, threadIdx.x );

  if( threadIdx.x == 0 ) etd_reduce[offset] = etd_shared[0];
#endif
#if 1
  warpReduce( etd_shared );
#endif
  __syncthreads();
  if( threadIdx.x == 0 ) etd_reduce[offset] = etd_shared[0];



//  etd_shared[threadIdx.x] = t1ai;

} /* end kernel */


__global__ void t1a_cuda_kernel( const int i, const int j, const int k, 
        const int no, const int nu,
        const double *v3, const double *voe_ij, const double *voe_ji, 
        const double *voe_ik, const double *voe_ki, 
        const double *voe_jk, const double *voe_kj,
        double *t1, const double *eh, const double *ep, double *etd_reduce )
{

__shared__ double etd_shared[SHARED_REDUCTION_SIZE];

  int a = blockIdx.x;
  int b;
  double dijk = eh[i-1] + eh[j-1] + eh[k-1];
  const double two = 2.0, four = 4.0, eight = 8.0, om = -1.0;
  double t1ai = 0.0, t1aj = 0.0, t1ak = 0.0;
  int c_major, c_minor;
  __shared__ double BUFFER_cab[16][128];		// BUFFER_cab[c_minor][threadIdx.x] <=> v3[INDX(c,a,b,nu)]
  __shared__ double BUFFER_cba[16][128];		// BUFFER_cba[c_minor][threadIdx.x] <=> v3[INDX(c,b,a,nu)]

  for( int idx = 0; idx < nu; idx += blockDim.x )
  {
        b = idx + threadIdx.x;
	  int temp = threadIdx.x - (threadIdx.x%16);
      for( c_major = 0; c_major < nu; c_major += 16 ) {
         for( c_minor = 0; c_minor < 16; c_minor++ ) {
	      BUFFER_cab[threadIdx.x%16][temp + c_minor] = v3[INDX(c_major + (threadIdx.x %16),a,idx + c_minor + temp,nu)];
	      BUFFER_cba[threadIdx.x%16][temp + c_minor] = v3[INDX(c_major + (threadIdx.x %16),idx + c_minor + temp,a,nu)];
         }
		 __syncthreads();
    if( b < nu && a < b )
    {
      for( c_minor = 0; c_minor < 16; c_minor++ )
      {
	    int c = c_major + c_minor;
		if ( c >= nu ) break;						// Break to c_major loop, if syncthreads allow.
        double dabc = ep[a] + ep[b] + ep[c];
        double denom = 1.0 / ( dijk - dabc );

                double abcbac = v3[INDX(a,b,c,nu)] - v3[INDX(b,a,c,nu)];
                double acbbca = v3[INDX(a,c,b,nu)] - v3[INDX(b,c,a,nu)];
//              double cabcba = v3[INDX(c,a,b,nu)] - v3[INDX(c,b,a,nu)];
              double cabcba = BUFFER_cab[c_minor][threadIdx.x] - BUFFER_cba[c_minor][threadIdx.x];

        double t3_ab1 = abcbac * two - acbbca;
        double t3_ab2 = acbbca * two - abcbac;
        double t3_ab3 = - ( abcbac * two + cabcba );
                double t3_ab4 = ( - acbbca * two ) + cabcba;
                double t3_ab5 = cabcba * two + abcbac;
                double t3_ab6 = ( - cabcba * two ) + acbbca;
/*
        double t3_ab1 = ( v3[INDX(a,b,c,nu)] - v3[INDX(b,a,c,nu)] ) * two
                      -   v3[INDX(a,c,b,nu)] + v3[INDX(b,c,a,nu)];

        double t3_ab2 = ( v3[INDX(a,c,b,nu)] - v3[INDX(b,c,a,nu)] ) * two
                      -   v3[INDX(a,b,c,nu)] + v3[INDX(b,a,c,nu)];

        double t3_ab3 = ( v3[INDX(b,a,c,nu)] - v3[INDX(a,b,c,nu)] ) * two
                      -   v3[INDX(c,a,b,nu)] + v3[INDX(c,b,a,nu)];

        double t3_ab5 = ( v3[INDX(c,a,b,nu)] - v3[INDX(c,b,a,nu)] ) * two
                      -   v3[INDX(b,a,c,nu)] + v3[INDX(a,b,c,nu)];

        double t3_ab4 = ( v3[INDX(b,c,a,nu)] - v3[INDX(a,c,b,nu)] ) * two
                      -   v3[INDX(c,b,a,nu)] + v3[INDX(c,a,b,nu)];

        double t3_ab6 = ( v3[INDX(c,b,a,nu)] - v3[INDX(c,a,b,nu)] ) * two
                      -   v3[INDX(b,c,a,nu)] + v3[INDX(a,c,b,nu)];
*/
//      if( a == 0 && b == 4 ) printf("c %d t3_ab1 %22.17e t3_ab2 %22.17e\n",
//         c,t3_ab1,t3_ab2);

        t1ai += ( t3_ab1*voe_jk[INDX(b,c,0,nu)] 
             +    t3_ab2*voe_kj[INDX(b,c,0,nu)] ) * denom;

        t1aj += ( t3_ab3*voe_ik[INDX(b,c,0,nu)] 
             +    t3_ab5*voe_ki[INDX(b,c,0,nu)] ) * denom;

        t1ak += ( t3_ab4*voe_ij[INDX(b,c,0,nu)] 
             +    t3_ab6*voe_ji[INDX(b,c,0,nu)] ) * denom;

      } /* end c_minor loop */
    } /* end if */
      } /* end c_major loop */
  } /* end idx loop */


  etd_shared[threadIdx.x] = t1ai;

  __syncthreads();

  double tempi = 0.0, tempj = 0.0, tempk = 0.0;

  int offi = INDX(a,i-1,0,nu);
#if 0
  if( threadIdx.x == 0 ) 
  {
    for( int idx = 0; idx < blockDim.x; idx++ )
    {
      tempi += etd_shared[idx];
    } /* end for */
    t1[offi] += tempi;
  } /* end if */
#endif
  warpReduce( etd_shared );
  __syncthreads();
  if( threadIdx.x == 0 ) t1[offi] += etd_shared[0];

  __syncthreads();

  etd_shared[threadIdx.x] = t1aj;

  __syncthreads();

  int offj = INDX(a,j-1,0,nu);
#if 0
  if( threadIdx.x == 0 ) 
  {
    for( int idx = 0; idx < blockDim.x; idx++ )
    {
      tempj += etd_shared[idx];
    } /* end for */
    t1[offj] += tempj;
  } /* end if */
#endif
  warpReduce( etd_shared );
  __syncthreads();
  if( threadIdx.x == 0 ) t1[offj] += etd_shared[0];

  __syncthreads();

  etd_shared[threadIdx.x] = t1ak;

  __syncthreads();

  int offk = INDX(a,k-1,0,nu);
#if 0
  if( threadIdx.x == 0 ) 
  {
    for( int idx = 0; idx < blockDim.x; idx++ )
    {
      tempk += etd_shared[idx];
    } /* end for */
    t1[offk] += tempk;
  } /* end if */
#endif
  warpReduce( etd_shared );
  __syncthreads();
  if( threadIdx.x == 0 ) t1[offk] += etd_shared[0];

//    Next calculation:

  double t1bi = 0.0, t1bj = 0.0, t1bk = 0.0;
  dijk = eh[i-1] + eh[j-1] + eh[k-1];

  b = blockIdx.x;
  for( int idx = 0; idx < nu; idx += blockDim.x )
  {
        a = idx + threadIdx.x;
/*
 * don't do the loop if my id is outside the bounds of nu
 */
	  int temp = threadIdx.x - (threadIdx.x%16);
      for( c_major = 0; c_major < nu; c_major += 16 ) {
         for( c_minor = 0; c_minor < 16; c_minor++ ) {
	      BUFFER_cab[threadIdx.x%16][temp + c_minor] = v3[INDX(c_major + (threadIdx.x %16),idx + c_minor + temp,b,nu)];
	      BUFFER_cba[threadIdx.x%16][temp + c_minor] = v3[INDX(c_major + (threadIdx.x %16),b,idx + c_minor + temp,nu)];
/*
int c = c_major + c_minor;
BUFFER_cab[c_minor][threadIdx.x] = v3[INDX(c,a,b,nu)];
BUFFER_cba[c_minor][threadIdx.x] = v3[INDX(c,b,a,nu)];
*/
         }
		 __syncthreads();
    if( a < nu && a < b )
    {
      for( c_minor = 0; c_minor < 16; c_minor++ )
      {
	    int c = c_major + c_minor;
		if ( c >= nu ) break;						// Break to c_major loop, if syncthreads allow.
        double dabc = ep[a] + ep[b] + ep[c];
        double denom = 1.0 / ( dijk - dabc );

                double abcbac = v3[INDX(a,b,c,nu)] - v3[INDX(b,a,c,nu)];
                double acbbca = v3[INDX(a,c,b,nu)] - v3[INDX(b,c,a,nu)];
//                double cabcba = v3[INDX(c,a,b,nu)] - v3[INDX(c,b,a,nu)];
              double cabcba = BUFFER_cab[c_minor][threadIdx.x] - BUFFER_cba[c_minor][threadIdx.x];
/*
if ( BUFFER_cab[c_minor][threadIdx.x] != v3[INDX(c,a,b,nu)] )
	printf("ThreadIdx.x=%d c_minor=%d c_major=%d c=%d a=%d b=%d nu=%d BUFFER_cab[c_minor][threadIdx.x]=%g v3[INDX(c,a,b,nu)]=%g\n",
	       threadIdx.x, c_minor, c_major, c, a, b, nu, BUFFER_cab[c_minor][threadIdx.x], v3[INDX(c,a,b,nu)] );
*/
        double t3_ab1 = abcbac * two - acbbca;
        double t3_ab2 = acbbca * two - abcbac;
        double t3_ab3 = - ( abcbac * two + cabcba );
                double t3_ab4 = ( - acbbca * two ) + cabcba;
                double t3_ab5 = cabcba * two + abcbac;
                double t3_ab6 = ( - cabcba * two ) + acbbca;
/*
        double t3_ab1 = ( v3[INDX(a,b,c,nu)] - v3[INDX(b,a,c,nu)] ) * two
                      -   v3[INDX(a,c,b,nu)] + v3[INDX(b,c,a,nu)];

        double t3_ab2 = ( v3[INDX(a,c,b,nu)] - v3[INDX(b,c,a,nu)] ) * two
                      -   v3[INDX(a,b,c,nu)] + v3[INDX(b,a,c,nu)];

        double t3_ab3 = ( v3[INDX(b,a,c,nu)] - v3[INDX(a,b,c,nu)] ) * two
                      -   v3[INDX(c,a,b,nu)] + v3[INDX(c,b,a,nu)];

        double t3_ab5 = ( v3[INDX(c,a,b,nu)] - v3[INDX(c,b,a,nu)] ) * two
                      -   v3[INDX(b,a,c,nu)] + v3[INDX(a,b,c,nu)];

        double t3_ab4 = ( v3[INDX(b,c,a,nu)] - v3[INDX(a,c,b,nu)] ) * two
                      -   v3[INDX(c,b,a,nu)] + v3[INDX(c,a,b,nu)];

        double t3_ab6 = ( v3[INDX(c,b,a,nu)] - v3[INDX(c,a,b,nu)] ) * two
                      -   v3[INDX(b,c,a,nu)] + v3[INDX(a,c,b,nu)];
*/
//      if( a == 0 && b == 4 ) printf("c %d t3_ab1 %22.17e t3_ab2 %22.17e\n",
//         c,t3_ab1,t3_ab2);

        t1bi += ( t3_ab1*voe_jk[INDX(a,c,0,nu)] 
             +    t3_ab2*voe_kj[INDX(a,c,0,nu)] ) * denom * om;

        t1bj += ( t3_ab3*voe_ik[INDX(a,c,0,nu)] 
             +    t3_ab5*voe_ki[INDX(a,c,0,nu)] ) * denom * om;

        t1bk += ( t3_ab4*voe_ij[INDX(a,c,0,nu)] 
             +    t3_ab6*voe_ji[INDX(a,c,0,nu)] ) * denom * om;

      } /* end c_minor loop */
    } /* end if */
      } /* end c_major loop */
  } /* end idx loop */


  __syncthreads();


  etd_shared[threadIdx.x] = t1bi;

  __syncthreads();

  tempi = 0.0, tempj = 0.0, tempk = 0.0;

  offi = INDX(b,i-1,0,nu);
#if 0
  if( threadIdx.x == 0 ) 
  {
    for( int idx = 0; idx < blockDim.x; idx++ )
    {
      tempi += etd_shared[idx];
    } /* end for */
    t1[offi] += tempi;
  } /* end if */
#endif
  warpReduce( etd_shared );
  __syncthreads();
  if( threadIdx.x == 0 ) t1[offi] += etd_shared[0];

  __syncthreads();

  etd_shared[threadIdx.x] = t1bj;

  __syncthreads();

  offj = INDX(b,j-1,0,nu);
#if 0
  if( threadIdx.x == 0 ) 
  {
    for( int idx = 0; idx < blockDim.x; idx++ )
    {
      tempj += etd_shared[idx];
    } /* end for */
    t1[offj] += tempj;
  } /* end if */
#endif
  warpReduce( etd_shared );
  __syncthreads();
  if( threadIdx.x == 0 ) t1[offj] += etd_shared[0];

  __syncthreads();

  etd_shared[threadIdx.x] = t1bk;

  __syncthreads();

  offk = INDX(b,k-1,0,nu);
#if 0
  if( threadIdx.x == 0 ) 
  {
    for( int idx = 0; idx < blockDim.x; idx++ )
    {
      tempk += etd_shared[idx];
    } /* end for */
    t1[offk] += tempk;
  } /* end if */
#endif
  warpReduce( etd_shared );
  __syncthreads();
  if( threadIdx.x == 0 ) t1[offk] += etd_shared[0];
} /* end t1a_cuda_kernel_b */



__global__ void reduce_etd_kernel( const long int size, const double *a,
     double *result )
{
  long int i;
  for( i = 0; i < size; i++ ) 
  {
    result[0] += a[i];
  } /* end for */
  return;
} /* end kernel */

} /* end extern C */

