#include <stdio.h>
#include "cuda_runtime.h"

#define INDX(a,b,c,ld) ( ( (c) * (ld) * (ld) ) \
                       + ( (b) * (ld) ) \
                       + ( (a) ) )

extern "C" {

__global__ void t1wt3_cuda_kernel( long int n, double *v )
{
  int ltr, n32, a, b, c, myindx;
  double temp;

  n32 = (int) n;
  ltr = ( n*n - n ) / 2;

  myindx = ( blockIdx.z * ( gridDim.x * gridDim.y ) 
           + blockIdx.y * ( gridDim.x ) 
	   + blockIdx.x ) 
           * blockDim.x
	   + threadIdx.x;

  if( myindx >= ( n32 * n32 * n32 )  ) return;
  
  a = ( myindx % ( n32 * n32 ) ) % n32;
  b = ( myindx % ( n32 * n32 ) ) / n32;
  c =   myindx / ( n32 * n32 );

//  printf("tidx %d bidx %d bidy %d myindex %d a %d b %d c %d\n",
 //     threadIdx.x, blockIdx.x, blockIdx.y, myindx, a, b, c );

  if ( c >= b ) return;

  temp = v[INDX(a, b, c, n32)];
  v[INDX(a, b, c, n32)] = v[INDX(a, c, b, n32)];
  v[INDX(a, c, b, n32)] = temp;

  return;

} /* end trant3_1_kernel */

} /* end extern C */
