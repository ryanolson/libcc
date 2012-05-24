#include <stdio.h>
#include <stdlib.h>
#include "cuda_runtime.h"
#include "cublas_v2.h"

extern "C" {

#include "ddi.h"

extern double MPI_Wtime();

#define CUDA_RC_CHECK(a) \
        do { a; \
        cudaError_t err = cudaGetLastError();                                                      \
        if ( err != cudaSuccess && err != cudaErrorSetOnActiveProcess ) { \
           printf(" +----------------------------------------\n"); \
           printf(" | ** CUDA ERROR! ** \n"); \
           printf(" | Error: \n"); \
           printf(" | Msg: %s\n", cudaGetErrorString(err) ); \
           printf(" | File: %s\n", __FILE__ );  \
           printf(" | Line: %d\n", __LINE__ );  \
           printf(" +----------------------------------------\n"); \
                 exit(-1);  \
           } \
        } while(0)

#define CUBLAS_RC_CHECK(a) do { a; } while(0)

#define SD_SHIFT(ptr,shift) \
        if(ptr == ptr ## _tail) ptr = ptr ## _head; \
        else                    ptr += shift;

#define CUDA_EVENT_COUNT        2
#define CUDA_D2H_EVENT_COUNT    3

static  cublasStatus_t stat;
static  cublasHandle_t cublas_hnd;
static  cudaError_t cudaStat;

typedef struct {
        int * ddi_handle;
} sd_t;

void streaming_dgemm(int m, int n, int k, double alpha, const sd_t *a, const sd_t *b,
                     double beta, const sd_t *c, int tm, int tn, int tk);

int main(int argc, char *argv[])
{
        int me, np, my, nn;
        int ddi_a, ddi_b, ddi_c;

     // read in m, n, k
        int m = atoi(argv[1]);
        int n = atoi(argv[2]);
        int k = atoi(argv[3]);

     // read in tile dimensions
        int tm = atoi(argv[4]); 
        int tn = atoi(argv[5]);
        int tk = atoi(argv[6]);

     // determine distriubted matrix requirements
        size_t dm = ( (long)m*k + (long)k*n + (long)m*n ) * 1.5;
        dm /= 1000000;

     // initialized ddi / mpi
        DDI_Init(argc,argv);
        DDI_NProc(&np, &me);
        DDI_NNode(&nn, &my);
        DDI_Memory(dm);

     // Create Full A, B and C
        DDI_Create(m, k, &ddi_a);
        DDI_Create(k, n, &ddi_b);
        DDI_Create(m, n, &ddi_c);

        sd_t a, b, c;
        a.ddi_handle = &ddi_a;
        b.ddi_handle = &ddi_b;
        c.ddi_handle = &ddi_c;

     // load cublas
        stat = cublasCreate( &cublas_hnd );

     // call streaming dgemm
        streaming_dgemm( m, n, k, 1.0, &a, &b, 1.0, &c, tm, tn, tk );

     // Clean up memory - distributed
        DDI_Destroy( ddi_c );
        DDI_Destroy( ddi_b );
        DDI_Destroy( ddi_a );

     // Destory cublas handle
        stat = cublasDestroy( cublas_hnd );

     // Finalize
        DDI_Finalize();
        return 0;
}

void streaming_dgemm(int m, int n, int k, double alpha, const sd_t *a, const sd_t *b, 
                     double beta, const sd_t *c, int tm, int tn, int tk)
{
     // id
        int me, np, my, nn;

     // timing
        double start_time, stop_time;

     // ddi handles
        int ddi_a = *(a->ddi_handle);
        int ddi_b = *(b->ddi_handle);
        int ddi_c = *(c->ddi_handle);

     // id
        DDI_NProc(&np, &me);
        DDI_NNode(&nn, &my);

     // determine the number of patches
        int m_patch_count = (m + tm - 1) / tm;
        int n_patch_count = (n + tn - 1) / tn;
        int k_patch_count = (k + tk - 1) / tk;
        int tiled_dgemm_count = m_patch_count * n_patch_count * k_patch_count;

     // Host Memory - define
        size_t ta_count = (long)tm * (long)tk;
        size_t tb_count = (long)tk * (long)tn;
        size_t tc_count = (long)tm * (long)tn;

        size_t ta_size = ta_count * sizeof(double);
        size_t tb_size = tb_count * sizeof(double);
        size_t tc_size = tc_count * sizeof(double);

        size_t h_a_size = ta_size * 3;
        size_t h_b_size = tb_size * 3;
        size_t h_c_size = tc_size * 3;

        double h_a_size_in_mb = (double) h_a_size / (1024*1024);
        double h_b_size_in_mb = (double) h_b_size / (1024*1024);
        double h_c_size_in_mb = (double) h_c_size / (1024*1024);
        double h_size_in_mb = h_a_size_in_mb + h_b_size_in_mb + h_c_size_in_mb;

     // Host Memory = alloc
        double * h_a = (double *) malloc( h_a_size );
        double * h_b = (double *) malloc( h_b_size );
        double * h_c = (double *) malloc( h_c_size );

     // Host Memory - register
        CUDA_RC_CHECK( cudaHostRegister( h_a, h_a_size, 0 ) );
        CUDA_RC_CHECK( cudaHostRegister( h_b, h_b_size, 0 ) );
        CUDA_RC_CHECK( cudaHostRegister( h_c, h_c_size, 0 ) );

     // Device Memory - define
        size_t d_a_size = ta_count * sizeof(double) * 2;
        size_t d_b_size = tb_count * sizeof(double) * 2;
        size_t d_c_size = tc_count * sizeof(double) * 2;

        double d_a_size_in_mb = (double) d_a_size / (1024*1024);
        double d_b_size_in_mb = (double) d_b_size / (1024*1024);
        double d_c_size_in_mb = (double) d_c_size / (1024*1024);
        double d_size_in_mb = d_a_size_in_mb + d_b_size_in_mb + d_c_size_in_mb;

     // Device Memory - alloc
        double * d_a = NULL;
        double * d_b = NULL;
        double * d_c = NULL;
        CUDA_RC_CHECK( cudaMalloc( (void **) &d_a, d_a_size ) );
        CUDA_RC_CHECK( cudaMalloc( (void **) &d_b, d_b_size ) );
        CUDA_RC_CHECK( cudaMalloc( (void **) &d_c, d_c_size ) );
        CUDA_RC_CHECK( cudaMemset( d_c, 0, d_c_size ) );

     // Create CUDA Streams
        cudaStream_t * stream = (cudaStream_t *) malloc( sizeof(cudaStream_t) * CUDA_STREAM_COUNT );

     // Initialize Streams & Record Initial events
        for(int i=0; i<CUDA_STREAM_COUNT; i++)  CUDA_RC_CHECK( cudaStreamCreate( &stream[i] ) );

     // Create CUDA Events
        cudaEvent_t * event     = (cudaEvent_t *) malloc( sizeof(cudaEvent_t) * CUDA_EVENT_COUNT );
        cudaEvent_t * d2h_event = (cudaEvent_t *) malloc( sizeof(cudaEvent_t) * CUDA_D2H_EVENT_COUNT );
        
     // Initialzie Cuda Events
        for(int i=0; i<CUDA_EVENT_COUNT; i++)     CUDA_RC_CHECK( cudaEventCreate( &event[i] ) );
        for(int i=0; i<CUDA_D2H_EVENT_COUNT; i++) CUDA_RC_CHECK( cudaEventCreate( &d2h_event[i] ) );
     

     // Prepare buffering pointers
        double * h_a_head = h_a;
        double * h_a_tail = h_a + ta_count*2;
        double * h_b_head = h_b;
        double * h_b_tail = h_b + tb_count*2;
        double * h_c_head = h_c;
        double * h_c_tail = h_c + tc_count*2;

        double * d_a_head = d_a;
        double * d_a_tail = d_a + ta_count;
        double * d_b_head = d_b;
        double * d_b_tail = d_b + tb_count;
        double * d_c_head = d_c;
        double * d_c_tail = d_c + tc_count;

        cudaStream_t * stream_head = stream;
        cudaStream_t * stream_tail = stream + 1;

        if(me == 0) 
        {
           printf("\n");
           printf("---------------- Streaming DGEMM using DDI ---------------- \n"); 
           printf("global dimensions: %8d %8d %8d\n",m,n,k);
           printf("tile dimensions  : %8d %8d %8d\n",tm,tn,tk);
           printf("\n");
           printf("dividing %d work packets over %d nodes\n", (m_patch_count*n_patch_count),nn);
           printf("each work packet consists of %d tiled dgemms\n",k_patch_count);
           printf("\n");
           printf("host requirements   = %.2lf + %.2lf + %.2lf = %.2lf\n", h_a_size_in_mb,
                   h_b_size_in_mb, h_c_size_in_mb, h_size_in_mb);
           printf("device requirements = %.2lf + %.2lf + %.2lf = %.2lf\n", d_a_size_in_mb,
                   d_b_size_in_mb, d_c_size_in_mb, d_size_in_mb);
           printf("\n");
           fflush(stdout);
        }

        DDI_Sync(1234);

     // DDOT Algorithm - Start Up
     // Dynamically load-balance all incoming patches
     // Stream patches of A & B from distributed memory 
     // C remains device resident over the summing index

        size_t patch_count = m_patch_count * n_patch_count;
        DDI_Patch a_patch, b_patch;
        DDI_Patch * c_patch = (DDI_Patch *) malloc( sizeof(DDI_Patch) * H_C_PATCH_COUNT );
        DDI_Patch * c_patch_head = c_patch;
        DDI_Patch * c_patch_tail = c_patch + H_C_PATCH_COUNT - 1;
        DDI_Patch * put_patch = c_patch;
        DDI_Patch * put_patch_head = c_patch_head;
        DDI_Patch * put_patch_tail = c_patch_tail;
        
        int * remaining_iterations = (int *) malloc( sizeof(int) * H_C_PATCH_COUNT );
        int * remaining_iterations_head = remaining_iterations;
        int * remaining_iterations_tail = remaining_iterations + H_C_PATCH_COUNT - 1; 

        double * hc_put   = hc;

        long ip = 0;            // number of patches streamed from the network
        size_t dlb_counter;     // dynamic load balancer

        double _one  = 1.0;
        double _zero = 0.0;
        const double *one = &_one;
        const double *zero = &_zero;

        double *palpha = &alpha;
        double *pbeta  = &beta;

        DDI_DLBNext(&dlb_counter);

        start_time = MPI_Wtime();
        while(dlb_counter < patch_count)
        {
           int patch_coord_i = dlb_counter / m_patch_count;
           int patch_coord_j = dlb_counter % m_patch_count;
        // Constant over these patch coordinates (i,j)
           a_patch.ilo = tm * patch_coord_i;
           a_patch.ihi = a_patch.ilo + tm - 1;
           b_patch.jlo = tn * patch_coord_j;
           b_patch.jhi = b_patch.jlo + tn - 1;
           for(int ik=0; ik<k_patch_count; ik++,ip++)
           {
               if( ip > k_patch_count ) {
                   if( *remaining_iterations == 1) CUDA_RC_CHECK( cudaEventSynchronize( d2h_event ) ); // it has to be done this time
                   cudaQuery = cudaEventQuery( d2h_event );  // does this return true every - i think yes (it is a value < the current event counter val)
                   if( *remaining_iterations == 1) assert( cudaQuery == cudaSuccess );
                   if(cudaQuery != cudaSuccess) break;
                   ncols_remaining = put_patch->jhi - put_patch.jlo + 1; 
                   if(*remaining_iterations > k_patch_count) *remaining_iterations = k_patch_count; // this allows use to finish early!
                   cols_per_iteration = ncols_remaining / *remaining_iterations; // if remaining_iterations is > k_patch_count, then we are ahead of 
                   memcpy(patch, put_patch, sizeof(DDI_Patch));                  // schedule and we are try to finish ahead of schedule
                   patch.jhi = patch.jlo + cols_per_iteration;
                   put_patch.jlo = patch.jhi + 1;
                   DDI_Put(ddi_c, &patch, hc_put);
                   hc_put += patch->size; // ensure ddi sets size
                   (*remaining_iterations)--;
                   if(*remaining_iterations == 0) {
                      SD_SHIFT( put_patch, 1 );
                      SD_SHIFT( d2h_event, 1 );
                      SD_SHIFT( h_c, tc_count );
                      SD_SHIFT( remaining_iterations, 1 );
                      hc_put = h_c;
               }   }
               a_patch.jlo = tk * ik;
               a_patch.jhi = a_patch.jlo + tk - 1;
               b_patch.ilo = a_patch.jlo;
               b_patch.ihi = a_patch.jhi;
               DDI_GetP(ddi_a, &a_patch, h_a);
               DDI_GetP(ddi_b, &b_patch, h_b);
               if(ip > 1) CUDA_RC_CHECK( cudaEventSynchronize( *event ) );
               CUDA_RC_CHECK( cudaMemcpyAsync( d_a, h_a, ta_size, cudaMemcpyHostToDevice, *stream ) );
               CUDA_RC_CHECK( cudaMemcpyAsync( d_b, h_b, tb_size, cudaMemcpyHostToDevice, *stream ) );
               CUBLAS_RC_CHECK( cublasSetStream( cublas_hnd, *stream ) );
               CUBLAS_RC_CHECK( cublasDgemm( cublas_hnd, CUBLAS_OP_N, CUBLAS_OP_N,
                                             tn, tm, tk, palpha, d_a, tn, d_b, tk, pbeta, d_c, tn ) );
               CUDA_RC_CHECK( cudaEventRecord( *event, stream ) );
               if(ik == patch_count-1) {
                  CUDA_RC_CHECK( cudaMemcpyAsync( h_c, d_c, tc_size, cudaMemcpyDeviceToHost, *stream ) );
                  CUDA_RC_CHECK( cudaEventRecord( d2h_event, stream ) );
                  *remaining_iterations = 2 * k_patch_count;
                  SD_SHIFT( d_c, tc_count );
               }
               SD_SHIFT( h_a, ta_count ); 
               SD_SHIFT( h_b, tb_count );
               SD_SHIFT( d_a, ta_count );
               SD_SHIFT( d_b, tb_count );
               SD_SHIFT( stream , 1 );
               SD_SHIFT( event , 1 );
           } // end ik loop over k_patch_count
           c_patch->ilo = a_patch.ilo;
           c_patch->ihi = a_patch.ihi;
           c_patch->jlo = b_patch.jlo;
           c_patch->jhi = b_patch.jhi;
           SD_SHIFT( c_patch, 1 );
           DDI_DLBNext(&dlb_counter);
        } // end while loop

        // the only thing we need to do outside the loops is process putting hc back into ddi_c!

        // finished loops ==> drain buffers - step 1 of 2
        // moves last get on to the device & executes the dgemm for the data already there
        if(ip > 1)
        {
            CUBLAS_RC_CHECK( cublasSetStream( cublas_hnd, *stream ) );
            CUBLAS_RC_CHECK( cublasDgemm( cublas_hnd, CUBLAS_OP_N, CUBLAS_OP_N,
                                          tn, tm, tk, palpha, d_a, tn, d_b, tk, pbeta, d_c, tn ) );
            SD_SHIFT( stream, 1 );
            SD_SHIFT( d_a, ta_count );
            SD_SHIFT( d_b, tb_count );
        }
  
        cudaStreamSynchronize( *stream );
        // add copy c off logic in here
        CUDA_RC_CHECK( cudaMemcpyAsync( d_a, h_a, ta_size, cudaMemcpyHostToDevice, *stream ) );
        CUDA_RC_CHECK( cudaMemcpyAsync( d_b, h_b, tb_size, cudaMemcpyHostToDevice, *stream ) );
        SD_SHIFT( h_a, ta_count );
        SD_SHIFT( h_b, tb_count );

        ++ip;

        // finished loops ==> drain buffers - step 2 of 2
        CUBLAS_RC_CHECK( cublasSetStream( cublas_hnd, *stream ) );
        CUBLAS_RC_CHECK( cublasDgemm( cublas_hnd, CUBLAS_OP_N, CUBLAS_OP_N,
                                      tn, tm, tk, palpha, d_a, tn, d_b, tk, pbeta, d_c, tn ) );

        // move last patch of c off the device
        CUDA_RC_CHECK( cudaMemcpyAsync( h_c, d_c, tc_size, cudaMemcpyDeviceToHost, *stream ) );
        CUDA_RC_CHECK( cudaStreamSynchronize( *stream ) );

        DDI_Sync(1234);
        stop_time = MPI_Wtime();

        if(me == 0) {
           printf("walltime = %.6lf\n", (stop_time - start_time));
           printf("\n");
           printf("---------------- Streaming DGEMM using DDI ---------------- \n"); 
           fflush(stdout);
        }

     // Destroy Streams
        CUDA_RC_CHECK( cudaStreamSynchronize( c_stream ) );
        CUDA_RC_CHECK( cudaStreamDestroy( c_stream ) );
        // for(int i=0; i<2; i++) {
        //    CUDA_RC_CHECK( cudaStreamSynchronize( stream[i] ) );
        //    CUDA_RC_CHECK( cudaStreamDestroy( stream[i] ) );
       //  }

     // Clean up memory - host
        free( h_a_head );
        free( h_b_head );
        free( h_c_head );

     // Clean up memory - device
        CUDA_RC_CHECK( cudaFree( d_a_head ) );
        CUDA_RC_CHECK( cudaFree( d_b_head ) );
        CUDA_RC_CHECK( cudaFree( d_c_head ) );
}
}
