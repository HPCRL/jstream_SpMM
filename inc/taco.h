#ifndef TACO_C_HEADERS
#define TACO_C_HEADERS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "../inc/reader.h"
#include "../inc/cacheflush.h"

/*
#include ""
#ifdef PAPI 
	#include <papi.h>
	#define PAPI_ERROR_CHECK(X) if((X)!=PAPI_OK) std::cerr<<X<<" Error \n";
#endif
*/
#define TACO_MIN(_a,_b) ((_a) < (_b) ? (_a) : (_b))
#define TACO_MAX(_a,_b) ((_a) > (_b) ? (_a) : (_b))
#define TACO_DEREF(_a) (((___context___*)(*__ctx__))->_a)
#ifndef TACO_TENSOR_T_DEFINED
#define TACO_TENSOR_T_DEFINED
#define rem(X) if(X!= NULL) {free(X); X=NULL;}
#define rem_index(X) if(X!= NULL) {free(X[1][0]); free(X[1][1]); X=NULL;}
#define REM_TENSOR(X) { rem(X->vals); rem_index(X->indices);}

typedef enum { taco_mode_dense, taco_mode_sparse } taco_mode_t;
typedef struct {
  int32_t      order;         // tensor order (number of modes)
  int32_t*     dimensions;    // tensor dimensions
  int32_t      csize;         // component size
  int32_t*     mode_ordering; // mode storage ordering
  taco_mode_t* mode_types;    // mode storage types
  uint32_t***   indices;       // tensor index data (per mode)
  TYPE*     vals;          // tensor values
  int32_t      vals_size;     // values array size
} taco_tensor_t;
#endif
/*
int cmp(const void *a, const void *b) {
  return *((const int*)a) - *((const int*)b);
}
*/
int taco_spmm(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B); 
int taco_sddmm(taco_tensor_t *D, taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C);
int taco_spmm_papi(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B);
int taco_sddmm_papi(taco_tensor_t *D, taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C);

taco_tensor_t* create_dense(int N, int M, int mode);
taco_tensor_t* create_sparse(taco_tensor_t* A, int mode );
taco_tensor_t* read_sparse();



#endif
