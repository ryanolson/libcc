#include <stdio.h>
#include "cuda_runtime.h"

#define INDX(a,b,c,ld) ( ( (c) * (ld) * (ld) ) \
                       + ( (b) * (ld) ) \
                       + ( (a) ) )

extern "C" {

__global__ void t1wt3_cuda_kernel( const int i, const int j, const int k, 
        const int no, const int nu,
        const double *v3, const double *voe_ij, const double *voe_ji, 
	const double *voe_ik, const double *voe_ki, 
	const double *voe_jk, const double *voe_kj,
        double *t1, const double *eh, const double *ep, double *etd_reduce )
{
#if 0
  myindx = ( blockIdx.z * ( gridDim.x * gridDim.y ) 
           + blockIdx.y * ( gridDim.x ) 
	   + blockIdx.x ) 
           * blockDim.x
	   + threadIdx.x;

  if( myindx >= ( nu * nu * nu )  ) return;
  
  a = ( myindx % ( nu * nu ) ) % nu;
  b = ( myindx % ( nu * nu ) ) / nu;
  c =   myindx / ( nu * nu );
#endif

__shared__ double etd_shared[1024];

  int a = blockIdx.x;
  int b = blockIdx.y;
  double dijk = eh[i-1] + eh[j-1] + eh[k-1];
  double x3 = 0.0;
  const double two = 2.0, four = 4.0, eight = 8.0;
  double d1,d2,d3,f;

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
//    if( a == 32 && b == 33 ) printf("tid %d c %d x3: %e\n",threadIdx.x,c,x3);
  } /* end for */


  etd_shared[threadIdx.x] = 0.0;
  etd_shared[threadIdx.x] = x3;

  __syncthreads();

//  if( a == 0 && b == 0 ) printf("tid %d x3: %e\n",threadIdx.x,x3);
  int offset = INDX(a, b, 0, gridDim.x );

  double temp = 0.0;

  if( threadIdx.x == 0 )
  {
//    for( int idx = 0; idx < min(nu,blockDim.x); idx++ )
    for( int idx = 0; idx < blockDim.x; idx++ )
    {
      temp += etd_shared[idx];
//      if(a ==0 && b == 0 ) printf("idx %d val %e\n",idx,etd_shared[idx]);
    } /* end for */
 //   if( a == 0 && b == 0 ) printf("offset %d x3accum: %e\n",offset,temp);
    etd_reduce[offset] = temp;
  } /* end if */

} /* end kernel */



__global__ void reduce_etd_kernel( const long int size, const double *a,
     double *result )
{
  long int i;
  for( i = 0; i < size; i++ ) 
  {
    result[0] += a[i];
  //  printf("i %d a %e\n",i,a[i] );
  } /* end for */
  return;
} /* end kernel */

} /* end extern C */
