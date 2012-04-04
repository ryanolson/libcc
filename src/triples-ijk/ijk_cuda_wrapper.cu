#define VM_INDEX(i,j) ( (no*nu*no)*(j-1) + (no*nu)*(i-1) )


static int ijk_gpu_initialized = 0;

/**
 *
 */
static long int no = -1;
static long int nu = -1;
static long int d_vvvo = -1;

/**
 * pointer to shared fortran arrays
 */
static double *t1 = NULL;
static double *t2 = NULL;
static double *vm = NULL;
static double *voe = NULL;
static double *ep = NULL;
static double *eh = NULL;

/**
 * cpu driver  storage for the [vv|vo] integrals
 * on the driver side, these can be as small as nutr*nu provided
 * there is a kernal and temp storage to expand the arrays to nu3
 */
static double *ve_i = NULL;
static double *ve_j = NULL;
static double *ve_k = NULL;

/**
 * iold, jold, kold indicate the ijk tuple from the previous
 * iteration.  if any are teh same, the ve_ data can be reused
 */
static long int iold = -1;
static long int jold = -1;
static long int kold = -1;

/**
 * GPU array pointers
 */
static double *d_eh; // IN
static double *d_ep; // IN
static double *d_vm; // IN
static double *d_t1; // IN/OUT
static double *d_x3; // OUT

/**
 * static local functions
 */
static void ijk_task(long mytask, int *i, int *j, int *k);
static void get_ve(int i, void *buff);


/**
 * Main IJK tuples driver
 */
void
ijk_gpu_driver(Integer *f_no, Integer *f_nu, double *f_t1, double *f_t2,
               double *f_vm, double *f_voe, double *f_eh, double *f_ep,
               Integer *f_vvvo, Integer *comm)
{
        size_t numbytes, gpu_bytes=0;
        
        long no2, no3;
        long nu2, nu3;
        long nou;
        
        MPI_Comm comm = (MPI_Comm) *f_comm;

        no = (int)*f_no;
        nu = (int)*f_nu;
        d_vvvo = (int)*f_vvvo;

        nou = no*nu;
        no2 = no*no;
        nu2 = nu*nu;
        no3 = no*no2;
        nu3 = nu*nu2;

     // Device Resident Arrays
     // ======================

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

     // Determine if t2 can be resident, otherwise register it
     // Determine if voe can be resident, otherwise register it

     // CPU Memory Resident Arrays
     // ==========================
      # if HAS_VE_EXPANSION
        numbytes = sizeof(double) * nutr * nu; cpu_bytes += (3*numbytes);
      # else
        numbytes = sizeof(double) * nu3; cpu_bytes += (3*numbytes);
      # endif
        ve_i = malloc(numbyte);
        ve_j = malloc(numbyte);
        ve_k = malloc(numbyte);

     // Set up basic load-balancing parameters
        n_ijk_tuples = (no*(no-1)*(no-2)) / 6;

     // Evenly divide the ijk-tuples over the GPU nodes
     // This is where the load-balancing "smarts" need to be improved
        div_even(n_ijk_tuples, np, me, &nr, &sr);

        for(ijk=sr; ijk<nr; ijk++)
        {
           ijk_task(ijk, &i, &j, &k);
           get_ve(i, ve_i);
           get_ve(j, ve_j);
           get_ve(k, ve_k);
           formv3( ... );
           t1wt3( ... );
        }

        return;
}


static void
ijk_task(long mytask, int *i, int *j, int *k)
{
        long cntr = 0;
        int ii,jj,kk;
        for(ii=1; ii<=no; ii++)
        for(jj=1; jj<=(ii-1); jj++)
        for(kk=1; kk<=(jj-1); kk++) 
        {
            if(cntr == mytask)
            {
               *i = ii;
               *j = jj;
               *k = kk;
               return;
            }
            ++cntr;
        }
        assert(0);
}

static void
get_ve(int i, void *buff)
{
        if(buff == ve_i && i == iold) return;
        if(buff == ve_j && j == jold) return;
        if(buff == ve_k && k == kold) return;

	DDI_Patch patch;
	patch.ilo = 0;
	patch.ihi = nutr-1;
	patch.jlo = i-1;
	patch.jhi = i-1;
        DDI_GetP(d_vvvo, &patch, buff);
      # if HAS_VE_EXPANSION == 0
        assert(0);
      # endif
        return;
}

static void
formv3(
        long int i,
        long int j,
        long int k
)
{
        cublasStatus_t stat;
        cublasHandle_t handle;
        cudaError_t cudaStat; 

        const double om = -1.0, zero = 0.0, one = 1.0;

        long int nou = no * nu;
        long int nu2 = nu * nu;
        long int nu3 = nu2 * nu;
        long int nutr = (nu2 + nu) / 2;

        double *d_t2_i, *d_t2_j, *d_t2_k;
        double *d_ve_i, *d_ve_j, *d_ve_k;
        double *d_vm_ij, *d_vm_ji, *d_vm_ik, *d_vm_ki, *d_vm_kj, *d_vm_jk;
        double *d_v3;

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




        numbytes = sizeof(double) * nu2 * no;
        cudaStat = cudaMalloc( (void **) &d_t2_j, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_t2_j, t2_j, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

//      numbytes = sizeof(double) * no * nu;
//      cudaStat = cudaMalloc( (void **) &d_vm_ki, numbytes );
//      CUDA_ERROR_CHECK();
//      cudaMemcpy( d_vm_ki, vm_ki, numbytes, cudaMemcpyHostToDevice );
//      CUDA_ERROR_CHECK();

        numbytes = sizeof(double) * nu3;
        cudaStat = cudaMalloc( (void **) &d_v3, numbytes );
        CUDA_ERROR_CHECK();

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

/**
* The first time d_ve_j is access we may have to expand it
*/
# if HAS_VE_EXPANSION
        numbytes = sizeof(double) * nutr * nu;
        cudaMemcpy( d_ve_j, ve_j, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

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

// final location for expansion of ve_i
/*
* expand from packed triangular to expanded upper triangular
*/
        expand_tr_kernel<<< 1, 1 >>>( nu, d_ve_j );
        CUDA_ERROR_CHECK();

/* 
* expand from upper triangular to full matrix
*/
        expand_trsq_kernel<<< grid, block >>>( nu, d_ve_j );
        CUDA_ERROR_CHECK();
# else // if HAS_VE_EXPANSION
        numbytes = sizeof(double) * nu3;
        cudaStat = cudaMalloc( (void **) &d_ve_j, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_ve_j, ve_j, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();
# endif // if HAS_VE_EXPANSION

/**
* Regardless of whether or not we have to expand d_ve_j,
* the responsibility of transforming d_ve_j by trant3_1 has been
* moved to the GPU
*/
        trant3_1_kernel<<< grid, block >>>( nu, d_ve_j );
        CUDA_ERROR_CHECK();

        numbytes = sizeof(double) * nu2 * no;
        cudaStat = cudaMalloc( (void **) &d_t2_k, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_t2_k, t2_k, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

        stat = cublasDgemm( handle,
                 CUBLAS_OP_T, CUBLAS_OP_T,
	   nu2, nu, nu, &one,
	   d_ve_j, nu,
	   &d_t2_k[nu2*(i-1)], nu, &one,
	   d_v3, nu2 );



//  printf("nu3 %d\n", nu3);

//  printf("block x y z %d %d %d\n",block.x,block.y,block.z);
//  printf("grid x y z %d %d %d\n",grid.x,grid.y,grid.z);
        trant3_1_kernel<<< grid, block >>>( nu, d_v3 );
        CUDA_ERROR_CHECK();


//      numbytes = sizeof(double) * nu * no;
//      cudaStat = cudaMalloc( (void **) &d_vm_ji, numbytes );
//      CUDA_ERROR_CHECK();
//      cudaMemcpy( d_vm_ji, vm_ji, numbytes, cudaMemcpyHostToDevice );
//      CUDA_ERROR_CHECK();

        stat = cublasDgemm( handle,
                 CUBLAS_OP_T, CUBLAS_OP_T,
	   nu, nu2, no, &om,
	   d_vm_ji, no,
	   d_t2_k, nu2, &one,
	   d_v3, nu );
        

        numbytes = sizeof(double) * nu3;
        cudaStat = cudaMalloc( (void **) &d_ve_k, numbytes );
        CUDA_ERROR_CHECK();
        numbytes = sizeof(double) * nutr * nu;
        cudaMemcpy( d_ve_k, ve_k, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

        expand_tr_kernel<<< 1, 1 >>>( nu, d_ve_k );
        CUDA_ERROR_CHECK();

        expand_trsq_kernel<<< grid, block >>>( nu, d_ve_k );
        CUDA_ERROR_CHECK();

        trant3_1_kernel<<< grid, block >>>( nu, d_ve_k );
        CUDA_ERROR_CHECK();

        
        stat = cublasDgemm( handle,
                 CUBLAS_OP_N, CUBLAS_OP_N,
	   nu, nu2, nu, &one,
	   &d_t2_j[nu2*(i-1)], nu,
	   d_ve_k, nu, &one,
	   d_v3, nu );

#if 1

        trant3_4_kernel<<< grid, block >>>( nu, d_v3 );
        CUDA_ERROR_CHECK();

        numbytes = sizeof(double) * nu2 * no;
        cudaStat = cudaMalloc( (void **) &d_t2_i, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_t2_i, t2_i, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

//      numbytes = sizeof(double) * no * nu;
//      cudaStat = cudaMalloc( (void **) &d_vm_kj, numbytes );
//      CUDA_ERROR_CHECK();
//      cudaMemcpy( d_vm_kj, vm_kj, numbytes, cudaMemcpyHostToDevice );
//      CUDA_ERROR_CHECK();

        stat = cublasDgemm( handle,
                 CUBLAS_OP_N, CUBLAS_OP_N,
	   nu2, nu, no, &om,
	   d_t2_i, nu2,
	   d_vm_kj, no, &one,
	   d_v3, nu2 );

        numbytes = sizeof(double) * nu3;
        cudaStat = cudaMalloc( (void **) &d_ve_i, numbytes );
        CUDA_ERROR_CHECK();
        numbytes = sizeof(double) * nutr * nu;
        cudaMemcpy( d_ve_i, ve_i, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

// final location for expansion of ve_i
        expand_tr_kernel<<< 1, 1 >>>( nu, d_ve_i );
        CUDA_ERROR_CHECK();

        expand_trsq_kernel<<< grid, block >>>( nu, d_ve_i );
        CUDA_ERROR_CHECK();

        trant3_1_kernel<<< grid, block >>>( nu, d_ve_i );
        CUDA_ERROR_CHECK();

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

//      numbytes = sizeof(double) * no * nu;
//      cudaStat = cudaMalloc( (void **) &d_vm_ik, numbytes );
//      CUDA_ERROR_CHECK();
//      cudaMemcpy( d_vm_ik, vm_ik, numbytes, cudaMemcpyHostToDevice );
//      CUDA_ERROR_CHECK();

        stat = cublasDgemm( handle,
                 CUBLAS_OP_T, CUBLAS_OP_T,
	   nu, nu2, no, &om,
	   d_vm_ik, no,
	   d_t2_j, nu2, &one,
	   d_v3, nu );


        trant3_1_kernel<<< grid, block >>>( nu, d_v3 );
        CUDA_ERROR_CHECK();

        
//      numbytes = sizeof(double) * no * nu;
//      cudaStat = cudaMalloc( (void **) &d_vm_jk, numbytes );
//      CUDA_ERROR_CHECK();
//      cudaMemcpy( d_vm_jk, vm_jk, numbytes, cudaMemcpyHostToDevice );
//      CUDA_ERROR_CHECK();

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

//      numbytes = sizeof(double) * no * nu;
//      cudaStat = cudaMalloc( (void **) &d_vm_ij, numbytes );
//      CUDA_ERROR_CHECK();
//      cudaMemcpy( d_vm_ij, vm_ij, numbytes, cudaMemcpyHostToDevice );
//      CUDA_ERROR_CHECK();


        stat = cublasDgemm( handle,
                 CUBLAS_OP_T, CUBLAS_OP_T,
	   nu, nu2, no, &om,
	   d_vm_ij, no,
	   d_t2_k, nu2, &one,
	   d_v3, nu );


        trant3_1_kernel<<< grid, block >>>( nu, d_v3 );
        CUDA_ERROR_CHECK();
#endif
/* 
* final copy back of v3
*/

//  numbytes = sizeof(double) * nu3;
//  cudaMemcpy( v3, d_v3, numbytes, cudaMemcpyDeviceToHost );
//  CUDA_ERROR_CHECK();

       cudaFree( d_vm_ij );
       CUDA_ERROR_CHECK();
       cudaFree( d_vm_ji );
       CUDA_ERROR_CHECK();
       cudaFree( d_vm_ik );
       CUDA_ERROR_CHECK();
       cudaFree( d_vm_ki );
       CUDA_ERROR_CHECK();
       cudaFree( d_vm_jk );
       CUDA_ERROR_CHECK();
       cudaFree( d_vm_kj );
       CUDA_ERROR_CHECK();
       cudaFree( d_t2_i );
       CUDA_ERROR_CHECK();
       cudaFree( d_t2_j );
       CUDA_ERROR_CHECK();
       cudaFree( d_t2_k );
       CUDA_ERROR_CHECK();
       cudaFree( d_ve_i );
       CUDA_ERROR_CHECK();
       cudaFree( d_ve_j );
       CUDA_ERROR_CHECK();
       cudaFree( d_ve_k );
       CUDA_ERROR_CHECK();
//  cudaFree( d_v3 );
//  CUDA_ERROR_CHECK();

