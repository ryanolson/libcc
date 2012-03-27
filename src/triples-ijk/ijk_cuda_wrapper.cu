
static int ijk_gpu_initialized = 0;

/**
 *
 */
static int no = -1;
static int nu = -1;

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
static int iold = -1;
static int jold = -1;
static int kold = -1;

/**
 * GPU array pointers
 */
static double *d_eh; // IN
static double *d_ep; // IN
static double *d_vm; // IN
static double *d_t1; // IN/OUT
static double *d_x3; // OUT

void
ijk_gpu_driver(Integer *f_no, Integer *f_nu, double *f_t1, double *f_t2,
               double *f_vm, double *f_voe, double *f_eh, double *f_ep)
{
        size_t numbytes, total_bytes=0;
        
        long no2, no3;
        long nu2, nu3;
        long nou;

        no = (int)*f_no;
        nu = (int)*f_nu;

        nou = no*nu;
        no2 = no*no;
        nu2 = nu*nu;
        no3 = no*no2;
        nu3 = nu*nu3;

     // Device Resident Arrays
     // ======================

     // d_eh
        numbytes = sizeof(double) * no; total_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_eh, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_eh, f_eh, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

     // d_ep
        numbytes = sizeof(double) * nu; total_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_ep, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_ep, f_ep, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

     // d_vm
        numbytes = sizeof(double) * no3 * nu; total_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_vm, numbytes );
        CUDA_ERROR_CHECK();
        cudaMemcpy( d_vm, f_vm, numbytes, cudaMemcpyHostToDevice );
        CUDA_ERROR_CHECK();

     // d_v3
        numbytes = sizeof(double) * nu3; total_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_v3, numbytes );
        CUDA_ERROR_CHECK();

     // d_ve_i
        numbytes = sizeof(double) * nu3; total_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_ve_i, numbytes );
        CUDA_ERROR_CHECK();

     // d_ve_j
        numbytes = sizeof(double) * nu3; total_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_ve_j, numbytes );
        CUDA_ERROR_CHECK();

     // d_ve_k
        numbytes = sizeof(double) * nu3; total_bytes += numbytes;
        cudaStat = cudaMalloc( (void **) &d_ve_k, numbytes );
        CUDA_ERROR_CHECK();

     // Determine if t2 can be resident, otherwise register it

     // Determine if voe can be resident, otherwise register it
     
     // Set up basic load-balancing

     // 
