#include <stdio.h>
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "cuda_kernels.h"

//#define NO_CUDA_DEBUG
#ifndef NO_CUDA_DEBUG
#define CUDA_ERROR_CHECK()                                                                     \
{                                                                                              \
    cudaError_t err = cudaGetLastError();                                                      \
     if ( err != cudaSuccess && err != cudaErrorSetOnActiveProcess ) { \
     printf(" +----------------------------------------\n"); \
     printf(" | ** CUDA ERROR! ** \n"); \
     printf(" | Error: \n"); \
     printf(" | Msg: %s\n", cudaGetErrorString(err) ); \
     printf(" | File: %s\n", __FILE__ );  \
     printf(" | Line: %d\n", __LINE__ );  \
     printf(" +----------------------------------------\n"); \
                 exit(-1);  } \
}
#else
#define CUDA_ERROR_CHECK() {}
#endif

#define VM_INDEX(i,j) ( (no*nu*no)*(j-1) + (no*nu)*(i-1) )
#define T2_INDEX(i)   ( (nu*nu*no)*(i-1) )

extern "C" {

typedef long Integer;

static long iold = -1;
static long jold = -1;
static long kold = -1;

static double *d_eh = NULL;
static double *d_ep = NULL;
static double *d_vm = NULL;
static double *d_v3 = NULL;
static double *d_ve_i = NULL;
static double *d_ve_j = NULL;
static double *d_ve_k = NULL;
static double *d_t2 = NULL;

void gpu_arrays_init_(
        Integer *f_no,
        Integer *f_nu,
        double *f_eh,
        double *f_ep,
        double *f_vm)
{

        long no = (long) *f_no;
        long nu = (long) *f_nu;
        long no3 = no*no*no;
        long nou2 = no*nu*nu;
        long nu3 = nu*nu*nu;
        size_t numbytes, gpu_bytes=0;

        cudaError_t cudaStat;

     // d_eh
        numbytes = sizeof(double) * no; gpu_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_eh, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_eh, f_eh, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

     // d_ep
        numbytes = sizeof(double) * nu; gpu_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_ep, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_ep, f_ep, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

     // d_vm
        numbytes = sizeof(double) * no3 * nu; gpu_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_vm, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_vm, f_vm, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

     // d_v3
        numbytes = sizeof(double) * nu3; gpu_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_v3, numbytes );
        CUDA_ERROR_CHECK();

     // d_ve_i
        numbytes = sizeof(double) * nu3; gpu_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_ve_i, numbytes );
        CUDA_ERROR_CHECK();

     // d_ve_j
        numbytes = sizeof(double) * nu3; gpu_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_ve_j, numbytes );
        CUDA_ERROR_CHECK();

     // d_ve_k
        numbytes = sizeof(double) * nu3; gpu_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_ve_k, numbytes );
        CUDA_ERROR_CHECK();

     // d_t2
        numbytes = sizeof(double) * nou2 * 3; gpu_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_t2, numbytes );
        CUDA_ERROR_CHECK();

}

void gpu_arrays_finalize_()
{
        cudaFree( d_eh );
        CUDA_ERROR_CHECK();
        cudaFree( d_ep );
        CUDA_ERROR_CHECK();
        cudaFree( d_vm );
        CUDA_ERROR_CHECK();
        cudaFree( d_v3 );
        CUDA_ERROR_CHECK();
        cudaFree( d_ve_i );
        CUDA_ERROR_CHECK();
        cudaFree( d_ve_j );
        CUDA_ERROR_CHECK();
        cudaFree( d_ve_k );
        CUDA_ERROR_CHECK();
        cudaFree( d_t2 );
        CUDA_ERROR_CHECK();
}

void ddcc_t_ijk_big_cuda_wrapper_(
    long int *p_nu, 
    long int *p_no,
    long int *p_i,
    long int *p_j,
    long int *p_k,
    double *t2_i,
    double *t2_j,
    double *t2_k,
    double *vm_ij,
    double *vm_ji,
    double *vm_ik,
    double *vm_ki,
    double *vm_jk,
    double *vm_kj,
    double *ve_i,
    double *ve_j,
    double *ve_k,
    double *voe_ij,
    double *voe_ji,
    double *voe_ik,
    double *voe_ki,
    double *voe_jk,
    double *voe_kj,
    double *t1,
    double *eh,
    double *ep,
    double *etd)
{

  cublasStatus_t stat;
  cublasHandle_t handle;
  cudaError_t cudaStat; 

  const double om = -1.0, zero = 0.0, one = 1.0;

  long int i = (*p_i); //fortran pointer offset
  long int j = (*p_j); //fortran pointer offset
  long int k = (*p_k); //fortran pointer offset
  long int no = *p_no;
  long int nu = *p_nu;
  long int nu2 = nu * nu;
  long int nu3 = nu2 * nu;
  long int nou2 = no * nu2;

  double *d_t2_i, *d_t2_j, *d_t2_k;
  double *d_vm_ij, *d_vm_ji, *d_vm_ik, *d_vm_ki, *d_vm_kj, *d_vm_jk;

  size_t numbytes;

/**
 * Determine VM offsets
 */
  d_vm_ij = d_vm + VM_INDEX(i,j);
  d_vm_ji = d_vm + VM_INDEX(j,i);
  d_vm_ik = d_vm + VM_INDEX(i,k);
  d_vm_ki = d_vm + VM_INDEX(k,i);
  d_vm_jk = d_vm + VM_INDEX(j,k);
  d_vm_kj = d_vm + VM_INDEX(k,j);

/**
 * Determine T2 offsets
 */
  d_t2_i = d_t2;
  d_t2_j = d_t2_i + nou2; 
  d_t2_k = d_t2_j + nou2;

  if(j != jold) {
     numbytes = sizeof(double) * no * nu2;
     cudaMemcpy( d_t2_j, t2_j, numbytes, cudaMemcpyHostToDevice );
     CUDA_ERROR_CHECK();
  }

// this copy is unnecessary because d_v3 is the product with a beta=0.0
// cudaMemcpy( d_v3, v3, numbytes, cudaMemcpyHostToDevice );
// CUDA_ERROR_CHECK();

  stat = cublasCreate( &handle );

  stat = cublasDgemm( handle,
           CUBLAS_OP_N, CUBLAS_OP_N,
	   nu2, nu, no, &om,
	   d_t2_j, nu2,
	   d_vm_ki, no, &zero,
	   d_v3, nu2 );

  if(j != jold) {
     numbytes = sizeof(double) * nu3;
     cudaMemcpy( d_ve_j, ve_j, numbytes, cudaMemcpyHostToDevice );
     CUDA_ERROR_CHECK();
  }

  if(k != kold) {
     numbytes = sizeof(double) * no * nu2;
     cudaMemcpy( d_t2_k, t2_k, numbytes, cudaMemcpyHostToDevice );
     CUDA_ERROR_CHECK();
  }

  stat = cublasDgemm( handle,
           CUBLAS_OP_T, CUBLAS_OP_T,
	   nu2, nu, nu, &one,
	   d_ve_j, nu,
	   &d_t2_k[nu2*(i-1)], nu, &one,
	   d_v3, nu2 );

  int blockx = 512;

  long int numblocks = ( nu3 / blockx ) + 1;

  dim3 block(blockx,1,1);
  long int gridx = 1;
  long int gridy = 1;

  if( numblocks <= 65535 )
  {
    gridx = numblocks;
  } else
  if( numblocks > 65535 && numblocks < (long int) 65535 * (long int )65535 )
  {
    gridx =  (long int) ceil( sqrt( (double) numblocks ) );
    gridy = gridx;
  } else
  {
    printf("too large grid requested...exiting\n");
    exit( 911 );
  } /* end if */

  dim3 grid( gridx, gridy, 1 );

  trant3_1_kernel<<< grid, block >>>( nu, d_v3 );
  CUDA_ERROR_CHECK();

  stat = cublasDgemm( handle,
           CUBLAS_OP_T, CUBLAS_OP_T,
	   nu, nu2, no, &om,
	   d_vm_ji, no,
	   d_t2_k, nu2, &one,
	   d_v3, nu );
  
  if(k != kold) {
     numbytes = sizeof(double) * nu3;
     cudaMemcpy( d_ve_k, ve_k, numbytes, cudaMemcpyHostToDevice );
     CUDA_ERROR_CHECK();
  }
  
  stat = cublasDgemm( handle,
           CUBLAS_OP_N, CUBLAS_OP_N,
	   nu, nu2, nu, &one,
	   &d_t2_j[nu2*(i-1)], nu,
	   d_ve_k, nu, &one,
	   d_v3, nu );


  trant3_4_kernel<<< grid, block >>>( nu, d_v3 );
  CUDA_ERROR_CHECK();

  if(i != iold) {
     numbytes = sizeof(double) * no * nu2;
     cudaMemcpy( d_t2_i, t2_i, numbytes, cudaMemcpyHostToDevice );
     CUDA_ERROR_CHECK();
   }

  stat = cublasDgemm( handle,
           CUBLAS_OP_N, CUBLAS_OP_N,
	   nu2, nu, no, &om,
	   d_t2_i, nu2,
	   d_vm_kj, no, &one,
	   d_v3, nu2 );

  if(i != iold) {
     numbytes = sizeof(double) * nu3;
     cudaMemcpy( d_ve_i, ve_i, numbytes, cudaMemcpyHostToDevice );
     CUDA_ERROR_CHECK();
  }


  stat = cublasDgemm( handle,
           CUBLAS_OP_T, CUBLAS_OP_T,
	   nu2, nu, nu, &one,
	   d_ve_i, nu,
	   &d_t2_k[nu2*(j-1)], nu, &one,
	   d_v3, nu2 );


  stat = cublasDgemm( handle,
           CUBLAS_OP_N, CUBLAS_OP_N,
	   nu, nu2, nu, &one,
	   &d_t2_i[nu2*(k-1)], nu,
	   d_ve_j, nu, &one,
	   d_v3, nu );

  stat = cublasDgemm( handle,
           CUBLAS_OP_T, CUBLAS_OP_T,
	   nu, nu2, no, &om,
	   d_vm_ik, no,
	   d_t2_j, nu2, &one,
	   d_v3, nu );


  trant3_1_kernel<<< grid, block >>>( nu, d_v3 );
  CUDA_ERROR_CHECK();


  stat = cublasDgemm( handle,
           CUBLAS_OP_N, CUBLAS_OP_N,
	   nu2, nu, no, &om,
	   d_t2_i, nu2,
	   d_vm_jk, no, &one,
	   d_v3, nu2 );


  stat = cublasDgemm( handle,
           CUBLAS_OP_T, CUBLAS_OP_T,
	   nu2, nu, nu, &one,
	   d_ve_i, nu,
	   &d_t2_j[nu2*(k-1)], nu, &one,
	   d_v3, nu2 );


  stat = cublasDgemm( handle,
           CUBLAS_OP_N, CUBLAS_OP_N,
	   nu, nu2, nu, &one,
	   &d_t2_i[nu2*(j-1)], nu,
	   d_ve_k, nu, &one,
	   d_v3, nu );


  stat = cublasDgemm( handle,
           CUBLAS_OP_T, CUBLAS_OP_T,
	   nu, nu2, no, &om,
	   d_vm_ij, no,
	   d_t2_k, nu2, &one,
	   d_v3, nu );


  trant3_1_kernel<<< grid, block >>>( nu, d_v3 );
  CUDA_ERROR_CHECK();

/* 
 * final copy back of v3
 */

//  numbytes = sizeof(double) * nu3;
//  cudaMemcpy( v3, d_v3, numbytes, cudaMemcpyDeviceToHost );
//  CUDA_ERROR_CHECK();

//  cudaFree( d_v3 );
//  CUDA_ERROR_CHECK();

  double x3;

  double *d_t1;
  double *d_voe_ij, *d_voe_ji, *d_voe_ik, *d_voe_ki, *d_voe_kj, *d_voe_jk;
  double *d_x3, *d_etd_reduce;


//  numbytes = sizeof(double) * nu3;
//  cudaStat = cudaMalloc( (void **) &d_v3, numbytes );
//  CUDA_ERROR_CHECK();
//  cudaMemcpy( d_v3, v3, numbytes, cudaMemcpyHostToDevice );
//  CUDA_ERROR_CHECK();

  numbytes = sizeof(double) * nu2;
  cudaStat = cudaMalloc( (void **) &d_voe_ij, numbytes );
  CUDA_ERROR_CHECK();
  cudaMemcpy( d_voe_ij, voe_ij, numbytes, cudaMemcpyHostToDevice );
  CUDA_ERROR_CHECK();

  numbytes = sizeof(double) * nu2;
  cudaStat = cudaMalloc( (void **) &d_voe_ji, numbytes );
  CUDA_ERROR_CHECK();
  cudaMemcpy( d_voe_ji, voe_ji, numbytes, cudaMemcpyHostToDevice );
  CUDA_ERROR_CHECK();
  
  numbytes = sizeof(double) * nu2;
  cudaStat = cudaMalloc( (void **) &d_voe_ik, numbytes );
  CUDA_ERROR_CHECK();
  cudaMemcpy( d_voe_ik, voe_ik, numbytes, cudaMemcpyHostToDevice );
  CUDA_ERROR_CHECK();

  numbytes = sizeof(double) * nu2;
  cudaStat = cudaMalloc( (void **) &d_voe_ki, numbytes );
  CUDA_ERROR_CHECK();
  cudaMemcpy( d_voe_ki, voe_ki, numbytes, cudaMemcpyHostToDevice );
  CUDA_ERROR_CHECK();

  numbytes = sizeof(double) * nu2;
  cudaStat = cudaMalloc( (void **) &d_voe_jk, numbytes );
  CUDA_ERROR_CHECK();
  cudaMemcpy( d_voe_jk, voe_jk, numbytes, cudaMemcpyHostToDevice );
  CUDA_ERROR_CHECK();

  numbytes = sizeof(double) * nu2;
  cudaStat = cudaMalloc( (void **) &d_voe_kj, numbytes );
  CUDA_ERROR_CHECK();
  cudaMemcpy( d_voe_kj, voe_kj, numbytes, cudaMemcpyHostToDevice );
  CUDA_ERROR_CHECK();

  numbytes = sizeof(double) * nu * no;
  cudaStat = cudaMalloc( (void **) &d_t1, numbytes );
  CUDA_ERROR_CHECK();
  cudaMemcpy( d_t1, t1, numbytes, cudaMemcpyHostToDevice );
  CUDA_ERROR_CHECK();

  numbytes = sizeof(double);
  cudaStat = cudaMalloc( (void **) &d_x3, numbytes );
  CUDA_ERROR_CHECK();
  cudaMemset( d_x3, 0, numbytes );
  CUDA_ERROR_CHECK();

  int device = 0;
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties( &deviceProp, device );

//  const int blockx = deviceProp.warpSize * 6;
  block.x = SHARED_REDUCTION_SIZE;
//    block.x = 128;

//  printf("warpSize is %d\n",blockx);
  block.y = 1;


  gridx = 1;
  gridy = 1;

  if( nu <= 65535 )
  {
    gridx = nu;
    gridy = nu;
  } else
  {
    printf("too large grid requested...exiting\n");
    exit( 911 );
  } /* end if */

  grid.x = gridx;
  grid.y = gridy;
  grid.z = 1;

  numbytes = sizeof(double) * (gridx * gridy) ;
  cudaStat = cudaMalloc( (void **) &d_etd_reduce, numbytes );
  CUDA_ERROR_CHECK();

/*
 * set the temporary array to zero it will be used for the reduction
 */

  cudaMemset( d_etd_reduce, 0, numbytes );
  CUDA_ERROR_CHECK();

//  printf("nu %d\n", nu);

//  printf("block x y z %d %d %d\n",block.x,block.y,block.z);
//  printf("grid x y z %d %d %d\n",grid.x,grid.y,grid.z);

  etd_cuda_kernel<<< grid, block >>>( i, j, k, no, nu, d_v3,
       d_voe_ij, d_voe_ji, d_voe_ik, d_voe_ki, d_voe_jk, d_voe_kj, 
       d_t1, d_eh, d_ep, d_etd_reduce );
  CUDA_ERROR_CHECK();

  reduce_etd_kernel<<<1,1>>>( gridx * gridy, d_etd_reduce, d_x3 );
  CUDA_ERROR_CHECK();

  grid.x = nu;
  grid.y = 1;

//  printf("block x y z %d %d %d\n",block.x,block.y,block.z);
//  printf("grid x y z %d %d %d\n",grid.x,grid.y,grid.z);

  t1a_cuda_kernel<<< grid, block >>>( i, j, k, no, nu, d_v3,
       d_voe_ij, d_voe_ji, d_voe_ik, d_voe_ki, d_voe_jk, d_voe_kj, 
       d_t1, d_eh, d_ep, d_etd_reduce );
  CUDA_ERROR_CHECK();

/* 
 * final copy back of v3 and t1
 */

#if 1
  numbytes = sizeof(double) * nu * no;
  cudaMemcpy( t1, d_t1, numbytes, cudaMemcpyDeviceToHost );
  CUDA_ERROR_CHECK();
#endif
  numbytes = sizeof(double);
  cudaMemcpy( &x3, d_x3, numbytes, cudaMemcpyDeviceToHost );
  CUDA_ERROR_CHECK();

/*
 * no need to copy v3 back to host
 */
// numbytes = sizeof(double) * nu3;
// cudaMemcpy( v3, d_v3, numbytes, cudaMemcpyDeviceToHost );
// CUDA_ERROR_CHECK();

//  printf("C etd %e x3 %e\n",*etd,x3);

  if( i == j || j == k ) 
  {
    *etd = (*etd) + x3 * 0.5;
  } /* end if */
  else
  {
    *etd = (*etd) + x3;
  } /* end else */

  cudaFree( d_voe_ij );
  CUDA_ERROR_CHECK();
  cudaFree( d_voe_ji );
  CUDA_ERROR_CHECK();
  cudaFree( d_voe_ik );
  CUDA_ERROR_CHECK();
  cudaFree( d_voe_ki );
  CUDA_ERROR_CHECK();
  cudaFree( d_voe_jk );
  CUDA_ERROR_CHECK();
  cudaFree( d_voe_kj );
  CUDA_ERROR_CHECK();
  cudaFree( d_t1 );
  CUDA_ERROR_CHECK();
  cudaFree( d_x3 );
  CUDA_ERROR_CHECK();
  cudaFree( d_etd_reduce );
  CUDA_ERROR_CHECK();

  iold = i;
  jold = j;
  kold = k;

  return;
   
} /* end void */
} /* end extern C */
