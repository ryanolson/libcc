#include "cuda_runtime.h"
#include "stdio.h"

#define INDX(a,b,c,ld) ( ( (c) * (ld) * (ld) ) \
                       + ( (b) * (ld) ) \
                       + ( (a) ) )

extern "C" {

__global__ void exp_trsq_kernel( long int n, double *src, double *dest )
{

  int n32, a, b, c, myindx;

  n32 = (int) n;
  int nu = n32;

  int nu2 = nu * nu;
  int nutr = (nu2 + nu) / 2;

  myindx = ( blockIdx.z * ( gridDim.x * gridDim.y ) 
           + blockIdx.y * ( gridDim.x ) 
	   + blockIdx.x ) 
           * blockDim.x
	   + threadIdx.x;

  if( myindx >= ( n32 * n32 * n32 )  ) return;
  
  a = ( myindx % ( n32 * n32 ) ) % n32;
  b = ( myindx % ( n32 * n32 ) ) / n32;
  c =   myindx / ( n32 * n32 );

  int row = min( a, b );
  int col = max( a, b );

  int offset = c * nutr;

  offset += ( col * (col + 1) / 2);
  offset += row;

    dest[INDX(a, b, c, n32)] = src[offset];

} /* end exp_trsq_kernel */

__global__ void expand_tr_kernel( long int n, double *v )
{
  int nu = (int) n;
  int nu2 = nu * nu;
  int ij = ( ( ( nu2 + nu ) / 2 ) * nu ) - 1;

  for( int kloop = nu - 1; kloop >= 0; kloop-- )
  {
    int koff = kloop * nu2;
    for( int iloop = nu - 1; iloop >= 0; iloop-- )
    {
      for( int jloop = iloop; jloop >= 0; jloop-- )
      {
        int joff = iloop*nu + jloop;
//        v[koff + joff] = v[ij];
        v[INDX(jloop, iloop, kloop, nu)] = v[ij];
        ij--;
      } /* end for */
    } /* end for */
  } /* end for */

  return; 
} /* end expand_tr */

__global__ void expand_trsq_kernel( long int n, double *v )
{
  int n32, a, b, c, myindx;

  n32 = (int) n;

  myindx = ( blockIdx.z * ( gridDim.x * gridDim.y ) 
           + blockIdx.y * ( gridDim.x ) 
	   + blockIdx.x ) 
           * blockDim.x
	   + threadIdx.x;

  if( myindx >= ( n32 * n32 * n32 )  ) return;
  
  a = ( myindx % ( n32 * n32 ) ) % n32;
  b = ( myindx % ( n32 * n32 ) ) / n32;
  c =   myindx / ( n32 * n32 );

  if( b > a ) return;

  v[INDX(a, b, c, n32)] = v[INDX(b, a, c, n32)];

  return;

} /* end expand_trsq */

__global__ void trant3_1_kernel( long int n, double *v )
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

__global__ void trant3_4_kernel( long int n, double *v )
{
  int n32, a, b, c, myindx;
  double temp;

  n32 = (int) n;

  myindx = ( blockIdx.z * ( gridDim.x * gridDim.y )
           + blockIdx.y * ( gridDim.x )
           + blockIdx.x )
           * blockDim.x
           + threadIdx.x;

#if 1
  if( myindx >= ( n32 * n32 * n32 )  ) return;

  a = ( myindx % ( n32 * n32 ) ) % n32;
  b = ( myindx % ( n32 * n32 ) ) / n32;
  c =   myindx / ( n32 * n32 );

//  printf("tidx %d bidx %d bidy %d myindex %d a %d b %d c %d\n",
 //     threadIdx.x, blockIdx.x, blockIdx.y, myindx, a, b, c );

  if( c > b ) return;
  if( a > c ) return;

  temp = v[INDX(a, b, c, n32)];
  v[INDX(a, b, c, n32)] = v[INDX(b, c, a, n32)];
  v[INDX(b, c, a, n32)] = v[INDX(c, a, b, n32)];
  v[INDX(c, a, b, n32)] = temp;

  if( ( b == c ) || ( c == a ) ) return;

  temp = v[INDX(b, a, c, n32)];
  v[INDX(b, a, c, n32)] = v[INDX(a, c, b, n32)];
  v[INDX(a, c, b, n32)] = v[INDX(c, b, a, n32)];
  v[INDX(c, b, a, n32)] = temp;

  return;
#endif
#if 0
  if( myindx > 0 ) return;
  printf("tidx %d bidx %d bidy %d myindex %d a %d b %d c %d\n",
      threadIdx.x, blockIdx.x, blockIdx.y, myindx, a, b, c );
  for( b = 0; b < n32; b++ )
  {
    for( c = 0; c <= b; c++ )
    {
      for( a = 0; a <= c; a++ )
      {
        temp = v[INDX(a, b, c, n32)];
        v[INDX(a, b, c, n32)] = v[INDX(b, c, a, n32)];
        v[INDX(b, c, a, n32)] = v[INDX(c, a, b, n32)];
        v[INDX(c, a, b, n32)] = temp;
        if( b == c || c == a ) {}
	else
        {
          temp = v[INDX(b, a, c, n32)];
          v[INDX(b, a, c, n32)] = v[INDX(a, c, b, n32)];
          v[INDX(a, c, b, n32)] = v[INDX(c, b, a, n32)];
          v[INDX(c, b, a, n32)] = temp;
        }
      }
    }
  }
  return;
#endif
} /* end trant3_4_kernel */

} /* end extern C */
