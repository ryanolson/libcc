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

int main(int argc, char *argv[])
{
        int me, np, my, nn;
        int ddi_a, ddi_b, ddi_c;
        double start_time, stop_time;
        cublasStatus_t stat;
        cublasHandle_t cublas_hnd;
        cudaError_t cudaStat;

     // load cublas
        stat = cublasCreate( &cublas_hnd );

     // read in m, n, k
        printf("argc=%d\n",argc);
        fflush(stdout);
        int m = atoi(argv[1]);
        int n = atoi(argv[2]);
        int k = atoi(argv[3]);

     // read in tile dimensions
        int tm = atoi(argv[4]);
        int tn = atoi(argv[5]);
        int tk = atoi(argv[6]);

     // determine the number of patches
        int m_patch_count = (m + tm - 1) / tm;
        int n_patch_count = (n + tn - 1) / tn;
        int k_patch_count = (k + tk - 1) / tk;
        int tiled_dgemm_count = m_patch_count * n_patch_count * k_patch_count;

     // determine distriubted matrix requirements
        printf("m,n,k = %d, %d, %d\n",m,n,k);
        fflush(stdout);
        size_t dm = ( (long)m*k + (long)k*n + (long)m*n ) * 3;
        printf("dm = %ld\n",dm);
        fflush(stdout);
        dm /= 1000000;
        printf("dm = %ld\n",dm);
        fflush(stdout);

     // initialized ddi / mpi
        DDI_Init(argc,argv);
        DDI_NProc(&np, &me);
        DDI_NNode(&nn, &my);
        DDI_Memory(dm);

     // Create Full A, B and C
        DDI_Create(m, k, &ddi_a);
        DDI_Create(k, n, &ddi_b);
        DDI_Create(m, n, &ddi_c);

     // Host Memory - define
        size_t ta_count = tm * tk;
        size_t tb_count = tk * tn;
        size_t tc_count = tm * tn;

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


     // Create CUDA Streams
        cudaStream_t * stream = (cudaStream_t *) malloc( sizeof(cudaStream_t) * 2 );
        for(int i=0; i<2; i++) CUDA_RC_CHECK( cudaStreamCreate( &stream[i] ) );

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

        long ip = 0;            // number of patches streamed from the network
        size_t dlb_counter;       // dynamic load balancer

        if(me == 0) 
        {
           printf("\n");
           printf("---------------- Streaming DGEMM using DDI ---------------- \n"); 
           printf("global dimensions: %8d %8d %8d\n",m,n,k);
           printf("tile dimensions  : %8d %8d %8d\n",tm,tn,tk);
           printf("\n");
           printf("performing %d tiled dgemms over %d nodes\n", tiled_dgemm_count, nn);
           printf("host requirements   = %.2lf + %.2lf + %.2lf = %.2lf\n", h_a_size_in_mb,
                   h_b_size_in_mb, h_c_size_in_mb, h_size_in_mb);
           printf("device requirements = %.2lf + %.2lf + %.2lf = %.2lf\n", d_a_size_in_mb,
                   d_b_size_in_mb, d_c_size_in_mb, d_size_in_mb);
           fflush(stdout);
        }

        DDI_Sync(1234);

     // DDOT Algorithm - Start Up
     // Dynamically load-balance patches of C
     // Stream A & B from Distributed Memory 
     // C remains device resident

     // Start Up
     // Get load-balance counter
     // Convert counter to determine patch of C
        DDI_Patch a_patch, b_patch;
        size_t patch_count = m_patch_count * n_patch_count;
       
        double _alpha = 1.0;
        double _beta  = 1.0;
        double *alpha = &_alpha;
        double *beta  = &_beta;

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
               a_patch.jlo = tk * ik;
               a_patch.jhi = a_patch.jlo + tk - 1;
               b_patch.ilo = a_patch.jlo;
               b_patch.ihi = a_patch.jhi;
   
               if(ip > 1)
               {
                   CUBLAS_RC_CHECK( cublasSetStream( cublas_hnd, *stream ) );
                   CUBLAS_RC_CHECK( cublasDgemm( cublas_hnd, CUBLAS_OP_N, CUBLAS_OP_N,
                                                 tn, tm, tk, alpha, d_a, tn, d_b, tk, beta, d_c, tn ) );
                   if(stream == stream_tail) stream = stream_head;
                   else                      stream++;
                   if(d_a == d_a_tail) d_a = d_a_head;
                   else                d_a += ta_count;
                   if(d_b == d_b_tail) d_b = d_b_head;
                   else                d_b += tb_count;
               }
   
               if(ip > 0)
               {
                   cudaStreamSynchronize( *stream );
                   /*
                   if(++dgemms_completed == k_patch_count) {
                      cudaMemcpyAsync( h_c, d_c, c_size, cudaMemcpyDeviceToHost, stream );
                      if(h_c == h_c_tail) h_c = h_c_head;
                      else                h_c += tc_count;
                      if(d_c == d_c_tail) d_c = d_c_head;
                      else                d_c += tc_count;
                      // zero out new d_c
                      dgemms_completed = 0;
                   }
                   */
                   CUDA_RC_CHECK( cudaMemcpyAsync( d_a, h_a, ta_size, cudaMemcpyHostToDevice, *stream ) );
                   CUDA_RC_CHECK( cudaMemcpyAsync( d_b, h_b, tb_size, cudaMemcpyHostToDevice, *stream ) );
                   if(h_a == h_a_tail) h_a = h_a_head;
                   else                h_a += ta_count;
                   if(h_b == h_b_tail) h_b = h_b_head;
                   else                h_b += tb_count;
               }

               DDI_GetP(ddi_a, &a_patch, h_a);
               DDI_GetP(ddi_b, &b_patch, h_b);

               // printf("end of loop %d\n", ip);
   
           } // end loop on ik

           DDI_DLBNext(&dlb_counter);

        } // end while loop


        // finished loops ==> drain buffers - step 1 of 2
        // moves last get on to the device & executes the dgemm for the data already there
        if(ip > 1)
        {
            CUBLAS_RC_CHECK( cublasSetStream( cublas_hnd, *stream ) );
            CUBLAS_RC_CHECK( cublasDgemm( cublas_hnd, CUBLAS_OP_N, CUBLAS_OP_N,
                                          tn, tm, tk, alpha, d_a, tn, d_b, tk, beta, d_c, tn ) );
            if(stream == stream_tail) stream = stream_head;
            else                      stream++;
            if(d_a == d_a_tail) d_a = d_a_head;
            else                d_a += ta_count;
            if(d_b == d_b_tail) d_b = d_b_head;
            else                d_b += tb_count;
        }
  
        if(ip > 0)
        {
            cudaStreamSynchronize( *stream );
            /*
            if(++dgemms_completed == k_patch_count) {
               cudaMemcpyAsync( h_c, d_c, c_size, cudaMemcpyDeviceToHost, stream );
               if(h_c == h_c_tail) h_c = h_c_head;
               else                h_c += tc_count;
               if(d_c == d_c_tail) d_c = d_c_head;
               else                d_c += tc_count;
               // zero out new d_c
               dgemms_completed = 0;
            }
            */
            CUDA_RC_CHECK( cudaMemcpyAsync( d_a, h_a, ta_size, cudaMemcpyHostToDevice, *stream ) );
            CUDA_RC_CHECK( cudaMemcpyAsync( d_b, h_b, tb_size, cudaMemcpyHostToDevice, *stream ) );
            if(h_a == h_a_tail) h_a = h_a_head;
            else                h_a += ta_count;
            if(h_b == h_b_tail) h_b = h_b_head;
            else                h_b += tb_count;
        }

        // finished loops ==> drain buffers - step 2 of 2
        ++ip;

        if(ip > 1)
        {
            CUBLAS_RC_CHECK( cublasSetStream( cublas_hnd, *stream ) );
            CUBLAS_RC_CHECK( cublasDgemm( cublas_hnd, CUBLAS_OP_N, CUBLAS_OP_N,
                                          tn, tm, tk, alpha, d_a, tn, d_b, tk, beta, d_c, tn ) );
        }

        CUDA_RC_CHECK( cudaMemcpyAsync( h_c, d_c, tc_size, cudaMemcpyDeviceToHost, *stream ) );
        cudaStreamSynchronize( *stream );

        DDI_Sync(1234);
        stop_time = MPI_Wtime();

        if(me == 0) {
           printf("walltime = %.6lf\n", (stop_time - start_time));
           fflush(stdout);
        }

     // Clean up memory
        DDI_Destroy(ddi_c);
        DDI_Destroy(ddi_b);
        DDI_Destroy(ddi_a);
        DDI_Finalize();
        return 0;
}
}
