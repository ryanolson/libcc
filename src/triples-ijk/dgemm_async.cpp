#include <openacc.h>
#include <cstdio>
#include <cstdlib>

#define __no_return__ 
#define __forceinline__ 
#define __align__(n) 
#define __thread__ 
#define __import__ 
#define __export__ 
#define __cdecl 
#define __annotate__(a) \
	        __attribute__((a))
#define __location__(a) \
	        __annotate__(a)
#define CUDARTAPI 

#include <cublas_v2.h>

/* ===== Function Prototypes ===== */
extern "C" void dgemm_async_setup_(long *numStreams, long *firstStreamId);
extern "C" void dgemm_async_(long *openaccStreamId, char *transa, char *transb, 
		long *m, long *n, long *k,
		double *alpha, double *a, long *lda, double *b, long *ldb,
		double *beta, double *c, long *ldc);
extern "C" void dgemm_async_shutdown_(void);


/* ===== Global Variables ===== */
cublasHandle_t *__da_cuhandles;
int __da_startStreamId, __da_numStreams;


/* ===== Function Definitions ===== */
extern "C" void dgemm_async_setup_(long *numStreams, long *firstStreamId) {

	int num = (*numStreams);
	__da_numStreams = num;
	__da_startStreamId = (int)(*firstStreamId);

	__da_cuhandles = (cublasHandle_t *) malloc(sizeof(cublasHandle_t)*num);
	for(int i = 0; i < num; ++i) {
		cudaStream_t custream;
		int streamId = (int)(*firstStreamId) + i;

		cublasCreate_v2(&__da_cuhandles[i]);
		cray_acc_get_async_info(streamId, &custream);
		cublasSetStream(__da_cuhandles[i], custream);
	}

        // printf("dgemm_async: initialized with %d streams with the first stream correspending to %d async_id\n",__da_numStreams, __da_startStreamId);
        // fflush(stdout);
}

extern "C" void dgemm_async_(long *openaccStreamId, char *transa, char *transb, 
		long *m, long *n, long *k,
		double *alpha, double *a, long *lda, double *b, long *ldb,
		double *beta, double *c, long *ldc) {

	cublasOperation_t transA = CUBLAS_OP_T, transB = CUBLAS_OP_T;

	if ((*transa) == 'n' || (*transa) == 'N') {
		transA = CUBLAS_OP_N;
	}
	if ((*transb) == 'n' || (*transb) == 'N') {
		transB = CUBLAS_OP_N;
	}

        // printf("dgemm_async %8d %8d %8d\n", (int)*m, (int)*n, (int)*k); fflush(stdout);

	cublasDgemm_v2(__da_cuhandles[(*openaccStreamId) - __da_startStreamId], transA, transB, 
			*m, *n, *k, 
			alpha, a, *lda, b, *ldb,
			beta, c, *ldc);
}

extern "C" void dgemm_async_shutdown_() {

	for(int i = 0; i <  __da_numStreams; ++i) {
		cublasDestroy_v2(__da_cuhandles[i]);
	}

	free(__da_cuhandles);
}
