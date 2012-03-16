//#define INDX(a,b,c,ld) ( ( (c) * (ld) * (ld) ) \
                       + ( (b) * (ld) ) \
                       + ( (a) ) )

//extern "C" {

__device__ void warpReduce( volatile double *sdata )
{
#if 0
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
    if( threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0 )
    {
      for( int idx = 1; idx < blockDim.x * blockDim.y * blockDim.z; idx++ )
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

#ifdef NOTEX
      double abc = v3[ INDX(a, b, c, nu) ];
      double acb = v3[ INDX(a, c, b, nu) ];
      double bac = v3[ INDX(b, a, c, nu) ];
      double bca = v3[ INDX(b, c, a, nu) ];
      double cba = v3[ INDX(c, b, a, nu) ];
      double cab = v3[ INDX(c, a, b, nu) ];

      d1 = abc;
      d2 = acb + cba + bac;
      d3 = bca + cab;
      f  = d1*eight - d2*four + d3*two;
      x3        += f*d1*denom;
#if 0
      d1 = v3[INDX(a, b, c, nu)];
      d2 = v3[INDX(a, c, b, nu)] + v3[INDX(c, b, a, nu)] 
	        + v3[INDX(b, a, c, nu)];
      d3 = v3[INDX(b, c, a, nu)] + v3[INDX(c, a, b, nu)];
      f  = d1*eight - d2*four + d3*two;
      x3        += f*d1*denom;
#endif
#else
      double abc = fetch_x_v3( INDX(a, b, c, nu) );
      double acb = fetch_x_v3( INDX(a, c, b, nu) );
      double bac = fetch_x_v3( INDX(b, a, c, nu) );
      double bca = fetch_x_v3( INDX(b, c, a, nu) );
      double cba = fetch_x_v3( INDX(c, b, a, nu) );
      double cab = fetch_x_v3( INDX(c, a, b, nu) );

      d1 = abc;
      d2 = acb + cba + bac;
      d3 = bca + cab;
      f  = d1*eight - d2*four + d3*two;
      x3        += f*d1*denom;

//      d1 = fetch_x_v3( INDX(a, b, c, nu) );
//      d2 = fetch_x_v3( INDX(a, c, b, nu) ) + fetch_x_v3( INDX(c, b, a, nu) )
//	        + fetch_x_v3( INDX(b, a, c, nu) );
//      d3 = fetch_x_v3( INDX(b, c, a, nu) ) + fetch_x_v3( INDX(c, a, b, nu) );
//      f  = d1*eight - d2*four + d3*two;
//      x3        += f*d1*denom;
#endif

      if( a == b ) goto loop_end;

#ifdef NOTEX
      d1 = bac;
      d2 = bca + cab + abc;
      d3 = acb + cba;
      f  = d1*eight - d2*four + d3*two;
      x3 += f*d1*denom;
#if 0
      d1 = v3[INDX(b, a, c, nu)];
      d2 = v3[INDX(b, c, a, nu)] + v3[INDX(c, a, b, nu)] 
	 + v3[INDX(a, b, c, nu)];
      d3 = v3[INDX(a, c, b, nu)] + v3[INDX(c, b, a, nu)];
      f  = d1*eight - d2*four + d3*two;
      x3 += f*d1*denom;
#endif
#else
      d1 = bac;
      d2 = bca + cab + abc;
      d3 = acb + cba;
      f  = d1*eight - d2*four + d3*two;
      x3 += f*d1*denom;

//      d1 = fetch_x_v3( INDX(b, a, c, nu) );
//      d2 = fetch_x_v3( INDX(b, c, a, nu) ) + fetch_x_v3( INDX(c, a, b, nu) )
//	 + fetch_x_v3( INDX(a, b, c, nu) );
//      d3 = fetch_x_v3( INDX(a, c, b, nu) ) + fetch_x_v3( INDX(c, b, a, nu) );
//      f  = d1*eight - d2*four + d3*two;
//      x3 += f*d1*denom;
#endif

    } /* end if */
    
loop_end:
  } /* end for */


  etd_shared[threadIdx.x] = 0.0;
  etd_shared[threadIdx.x] = x3;

  __syncthreads();

  int offset = INDX(a, b, 0, gridDim.x );

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
  const double two = 2.0, om = -1.0;
  double t1ai = 0.0, t1aj = 0.0, t1ak = 0.0;
  double t1bi = 0.0, t1bj = 0.0, t1bk = 0.0;
  int ti1d = threadIdx.z * ( blockDim.x * blockDim.y )
           + threadIdx.y * ( blockDim.x ) 
	   + threadIdx.x;

  for( b = 0; b < nu; b++ )
  { 
    for( int idx = 0; idx < nu; idx += ( blockDim.x * blockDim.y * blockDim.z) )
    {
      int c = idx + ti1d;

/*
 * don't do the loop if my id is outside the bounds of nu
 */

      if( c < nu )
      {
        if( a > b ) goto loop_end;
        if( a == b && b == c ) goto loop_end;
        double dabc = ep[a] + ep[b] + ep[c];
        double denom = 1.0 / ( dijk - dabc );
        if( a == b ) goto loop_end;

#ifdef NOTEX
        double abc = v3[ INDX(a, b, c, nu) ];
        double acb = v3[ INDX(a, c, b, nu) ];
        double bac = v3[ INDX(b, a, c, nu) ];
        double bca = v3[ INDX(b, c, a, nu) ];
        double cba = v3[ INDX(c, b, a, nu) ];
        double cab = v3[ INDX(c, a, b, nu) ];

        double t3_ab1 = ( abc - bac ) * two
   	              -   acb + bca;

        double t3_ab2 = ( acb - bca ) * two
	              -   abc + bac;

        double t3_ab3 = ( bac - abc ) * two
   	              -   cab + cba;

        double t3_ab5 = ( cab - cba )  * two
	              -   bac + abc;

        double t3_ab4 = ( bca - acb ) * two
	              -   cba + cab;

        double t3_ab6 = ( cba - cab ) * two
	              -   bca + acb;
#if 0
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
#endif
#else
        double abc = fetch_x_v3( INDX(a, b, c, nu) );
        double acb = fetch_x_v3( INDX(a, c, b, nu) );
        double bac = fetch_x_v3( INDX(b, a, c, nu) );
        double bca = fetch_x_v3( INDX(b, c, a, nu) );
        double cba = fetch_x_v3( INDX(c, b, a, nu) );
        double cab = fetch_x_v3( INDX(c, a, b, nu) );

        double t3_ab1 = ( abc - bac ) * two
   	              -   acb + bca;

        double t3_ab2 = ( acb - bca ) * two
	              -   abc + bac;

        double t3_ab3 = ( bac - abc ) * two
   	              -   cab + cba;

        double t3_ab5 = ( cab - cba )  * two
	              -   bac + abc;

        double t3_ab4 = ( bca - acb ) * two
	              -   cba + cab;

        double t3_ab6 = ( cba - cab ) * two
	              -   bca + acb;
#endif



//      if( a == 0 && b == 4 ) printf("c %d t3_ab1 %22.17e t3_ab2 %22.17e\n",
//	   c,t3_ab1,t3_ab2);

        t1ai += ( t3_ab1*voe_jk[INDX(b,c,0,nu)] 
	     +    t3_ab2*voe_kj[INDX(b,c,0,nu)] ) * denom;

        t1aj += ( t3_ab3*voe_ik[INDX(b,c,0,nu)] 
	     +    t3_ab5*voe_ki[INDX(b,c,0,nu)] ) * denom;

        t1ak += ( t3_ab4*voe_ij[INDX(b,c,0,nu)] 
	     +    t3_ab6*voe_ji[INDX(b,c,0,nu)] ) * denom;

      } /* end if */

loop_end:

    } /* end idx loop */

  } /* end b loop */


  etd_shared[ti1d] = t1ai;

  __syncthreads();

  int offi = INDX(a,i-1,0,nu);
#if 0
  if( ti1d == 0 ) 
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
  if( ti1d == 0 ) t1[offi] += etd_shared[0];

  __syncthreads();

  etd_shared[ti1d] = t1aj;

  __syncthreads();

  int offj = INDX(a,j-1,0,nu);
#if 0
  if( ti1d == 0 ) 
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
  if( ti1d == 0 ) t1[offj] += etd_shared[0];

  __syncthreads();

  etd_shared[ti1d] = t1ak;

  __syncthreads();

  int offk = INDX(a,k-1,0,nu);
#if 0
  if( ti1d == 0 ) 
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
  if( ti1d == 0 ) t1[offk] += etd_shared[0];

#if 1

  b = blockIdx.x;

  for( a = 0; a < nu; a++ )
  { 
    for( int idx = 0; idx < nu; idx += ( blockDim.x * blockDim.y * blockDim.z ) )
    {
      int c = idx + ti1d;

/*
 * don't do the loop if my id is outside the bounds of nu
 */

      if( c < nu )
      {
        if( a > b ) goto loop_end1;
        if( a == b && b == c ) goto loop_end1;
        double dabc = ep[a] + ep[b] + ep[c];
        double denom = 1.0 / ( dijk - dabc );
        if( a == b ) goto loop_end1;

#ifdef NOTEX
        double abc = v3[ INDX(a, b, c, nu) ];
        double acb = v3[ INDX(a, c, b, nu) ];
        double bac = v3[ INDX(b, a, c, nu) ];
        double bca = v3[ INDX(b, c, a, nu) ];
        double cba = v3[ INDX(c, b, a, nu) ];
        double cab = v3[ INDX(c, a, b, nu) ];

        double t3_ab1 = ( abc - bac ) * two
   	              -   acb + bca;

        double t3_ab2 = ( acb - bca ) * two
	              -   abc + bac;

        double t3_ab3 = ( bac - abc ) * two
   	              -   cab + cba;

        double t3_ab5 = ( cab - cba )  * two
	              -   bac + abc;

        double t3_ab4 = ( bca - acb ) * two
	              -   cba + cab;

        double t3_ab6 = ( cba - cab ) * two
	              -   bca + acb;
#if 0
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
#endif
#else
        double abc = fetch_x_v3( INDX(a, b, c, nu) );
        double acb = fetch_x_v3( INDX(a, c, b, nu) );
        double bac = fetch_x_v3( INDX(b, a, c, nu) );
        double bca = fetch_x_v3( INDX(b, c, a, nu) );
        double cba = fetch_x_v3( INDX(c, b, a, nu) );
        double cab = fetch_x_v3( INDX(c, a, b, nu) );

        double t3_ab1 = ( abc - bac ) * two
   	              -   acb + bca;

        double t3_ab2 = ( acb - bca ) * two
	              -   abc + bac;

        double t3_ab3 = ( bac - abc ) * two
   	              -   cab + cba;

        double t3_ab5 = ( cab - cba )  * two
	              -   bac + abc;

        double t3_ab4 = ( bca - acb ) * two
	              -   cba + cab;

        double t3_ab6 = ( cba - cab ) * two
	              -   bca + acb;
#endif

        t1bi += ( t3_ab1*voe_jk[INDX(a,c,0,nu)] 
	     +    t3_ab2*voe_kj[INDX(a,c,0,nu)] ) * denom * om;

        t1bj += ( t3_ab3*voe_ik[INDX(a,c,0,nu)] 
	     +    t3_ab5*voe_ki[INDX(a,c,0,nu)] ) * denom * om;

        t1bk += ( t3_ab4*voe_ij[INDX(a,c,0,nu)] 
	     +    t3_ab6*voe_ji[INDX(a,c,0,nu)] ) * denom * om;

      } /* end if */

loop_end1:

    } /* end idx loop */

  } /* end a loop */


  __syncthreads();


  etd_shared[ti1d] = t1bi;

  __syncthreads();

  offi = INDX(b,i-1,0,nu);
#if 0
  if( ti1d == 0 ) 
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
  if( ti1d == 0 ) t1[offi] += etd_shared[0];

  __syncthreads();

  etd_shared[ti1d] = t1bj;

  __syncthreads();

  offj = INDX(b,j-1,0,nu);
#if 0
  if( ti1d == 0 ) 
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
  if( ti1d == 0 ) t1[offj] += etd_shared[0];

  __syncthreads();

  etd_shared[ti1d] = t1bk;

  __syncthreads();

  offk = INDX(b,k-1,0,nu);
#if 0
  if( ti1d == 0 ) 
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
  if( ti1d == 0 ) t1[offk] += etd_shared[0];
#endif
} /* end t1a_cuda_kernel */



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

//} /* end extern C */
