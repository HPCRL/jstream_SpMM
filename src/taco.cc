#ifndef TACO_CC
#define TACO_CC

#include  "../inc/taco.h"
#include <algorithm>


struct timeval starttime, midtime, endtime,timediff;
double elapsed;


/*
double findMedian(double a[], int n) 
{ 
    // First we sort the array 
    //cout<<"here"<<endl;
    std::sort(a, a + n); 
  
    //cout<<"here"<<endl;
    // check for even case 
    if (n % 2 != 0) 
        return (double)a[n / 2]; 
  
    return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0; 
} 
*/

int taco_spmm(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) 
{
	int C1_dimension = (int)(C->dimensions[0]);
	int C2_dimension = (int)(C->dimensions[1]);
	double* restrict C_vals2 = (double*)(C->vals);
	int A1_dimension = (int)(A->dimensions[0]);
	int* restrict A2_pos = (int*)(A->indices[1][0]);
	int* restrict A2_crd = (int*)(A->indices[1][1]);
	double* restrict A_vals = (double*)(A->vals);
	int B1_dimension = (int)(B->dimensions[0]);
	int B2_dimension = (int)(B->dimensions[1]);
	double* restrict B_vals2 = (double*)(B->vals);

	#pragma omp parallel for schedule(static)
	for (int32_t pC = 0; pC < (C1_dimension * C2_dimension); pC++) {
		C_vals2[pC] = 0.0;
	}

	uint64_t layers = 16;
	uint64_t mem_capacity = 4000000000; // 4 billion doubles = 32 GB

	if ( ((uint64_t) layers) * (((uint64_t) B1_dimension)* ((uint64_t) B2_dimension) + ((uint64_t) C1_dimension)* ((uint64_t) C2_dimension) ) > mem_capacity)
		layers = mem_capacity / (((uint64_t) B1_dimension)* ((uint64_t) B2_dimension) + ((uint64_t) C1_dimension)* ((uint64_t) C2_dimension) );
	if (layers < 1)
		layers = 1;

//	cout << "layers "<<layers<<endl;

	uint64_t sizeB = ((uint64_t) layers )*((uint64_t) B1_dimension)* ((uint64_t) B2_dimension);
	uint64_t sizeC = ((uint64_t) layers )*((uint64_t) C1_dimension)* ((uint64_t) C2_dimension);
	
	double* B_vals_all = (TYPE* ) malloc( ((uint64_t) sizeof(TYPE))*sizeB );
	double* C_vals_all = (TYPE* ) malloc( ((uint64_t) sizeof(TYPE))*sizeC );

	uint64_t singleB = ((uint64_t) B1_dimension)* ((uint64_t) B2_dimension);
	uint64_t singleC = ((uint64_t) C1_dimension)* ((uint64_t) C2_dimension);


	for (uint64_t i = 0 ; i< sizeB ; i++)
	{
		B_vals_all[i] = B_vals2[i%singleB];
	}


	for (uint64_t i = 0 ; i< sizeC; i ++)
	{
		C_vals_all[i] = C_vals2[i%singleC];
	}	


	double *elapsed_arr = new double[REP];

	for(int repeat = 0; repeat < REP ; repeat ++ )
	{
		cacheFlush(X,Y);
		//gettimeofday(&starttime,NULL); 
		uint64_t jumpB = sizeB / layers * (repeat % layers);
		uint64_t jumpC = sizeC / layers * (repeat % layers);
		double* B_vals = B_vals_all + jumpB;
		double* C_vals = C_vals_all + jumpC;
		auto start = std::chrono::high_resolution_clock::now();
		#pragma omp parallel for schedule(runtime)
		for (int32_t i = 0; i < A1_dimension; i++) {
			for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++) {
				int32_t j = A2_crd[pA2];
				for (int32_t k = 0; k < B2_dimension; k++) {
					int32_t pC2 = i * C2_dimension + k;
					int32_t pB2 = j * B2_dimension + k;
					C_vals[pC2] = C_vals[pC2] + A_vals[pA2] * B_vals[pB2];
				}
			}
		}
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		elapsed = diff.count();
		elapsed_arr[repeat] = elapsed;
		//gettimeofday(&endtime,NULL); 
		//elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
		//cout<<"TACO SpMM K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000 ;
		//cout<< endl;
	}

	elapsed = findMedian(elapsed_arr,REP);
	cout<<"TACO SpMM median K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000 ;
		cout<< endl;

	elapsed = elapsed_arr[0];
	cout<<"TACO SpMM max K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000 ;
		cout<< endl;

	elapsed = elapsed_arr[REP-1];
	cout<<"TACO SpMM min K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000 ;
		cout<< endl;

	//gettimeofday(&endtime,NULL); 
	//elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	//cout<<"TACO SpMM K= "<<C2_dimension<<", time= " << elapsed/REP << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000 * REP;
	//cout<< endl;
	return 0;
}


int taco_sddmm(taco_tensor_t *D, taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C)
{
	int D1_dimension = (int)(D->dimensions[0]);
	double* restrict D_vals = (double*)(D->vals);
	int A1_dimension = (int)(A->dimensions[0]);
	int* restrict A2_pos = (int*)(A->indices[1][0]);
	int* restrict A2_crd = (int*)(A->indices[1][1]);
	double* restrict A_vals = (double*)(A->vals);
	int B1_dimension = (int)(B->dimensions[0]);
	int B2_dimension = (int)(B->dimensions[1]);
	double* restrict B_vals2 = (double*)(B->vals);
	int C1_dimension = (int)(C->dimensions[0]);
	int C2_dimension = (int)(C->dimensions[1]);
	double* restrict C_vals2 = (double*)(C->vals);

	int32_t pD2 = 0;

	//gettimeofday(&starttime,NULL);
	double* elapsed_arr = new double[REP];

	uint64_t layers = 16;
	uint64_t mem_capacity = 4000000000; // 4 billion doubles = 32 GB

	if ( ((uint64_t) layers) * (((uint64_t) B1_dimension)* ((uint64_t) B2_dimension) + ((uint64_t) C1_dimension)* ((uint64_t) C2_dimension) ) > mem_capacity)
		layers = mem_capacity / (((uint64_t) B1_dimension)* ((uint64_t) B2_dimension) + ((uint64_t) C1_dimension)* ((uint64_t) C2_dimension) );
	if (layers < 1)
		layers = 1;

//	cout << "layers "<<layers<<endl;

	uint64_t sizeB = ((uint64_t) layers )*((uint64_t) B1_dimension)* ((uint64_t) B2_dimension);
	uint64_t sizeC = ((uint64_t) layers )*((uint64_t) C1_dimension)* ((uint64_t) C2_dimension);
	
	double* B_vals_all = (TYPE* ) malloc( ((uint64_t) sizeof(TYPE))*sizeB );
	double* C_vals_all = (TYPE* ) malloc( ((uint64_t) sizeof(TYPE))*sizeC );

	uint64_t singleB = ((uint64_t) B1_dimension)* ((uint64_t) B2_dimension);
	uint64_t singleC = ((uint64_t) C1_dimension)* ((uint64_t) C2_dimension);

	
	for (uint64_t i = 0 ; i< sizeB ; i++)
	{
		B_vals_all[i] = B_vals2[i%singleB];
	}


	for (uint64_t i = 0 ; i< sizeC; i ++)
	{
		C_vals_all[i] = C_vals2[i%singleC];
	}	





	for(int repeat = 0; repeat < REP ; repeat ++ )
	{
		cacheFlush(X,Y);
		//gettimeofday(&starttime,NULL); 
		uint64_t jumpB = sizeB / layers * (repeat % layers);
		uint64_t jumpC = sizeC / layers * (repeat % layers);
		double* B_vals = B_vals_all + jumpB;
		double* C_vals = C_vals_all + jumpC;
		//cout<<"here"<<endl;
		pD2 = 0;
		auto start = std::chrono::high_resolution_clock::now();
		#pragma omp parallel for 
		for (int32_t i = 0; i < C1_dimension; i++) 
		{
			for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++) 
			{
				int32_t j = A2_crd[pA2];
				D_vals[pD2] = 0.0;
				for (int32_t k = 0; k < C2_dimension; k++) 
				{
					int32_t pB2 = j * B2_dimension + k;
					int32_t pC2 = i * C2_dimension + k;
					D_vals[pD2] = D_vals[pD2] + A_vals[pA2] * B_vals[pB2] * C_vals[pC2];
				}
				pD2++;
			}
			
		}
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = end-start;
		elapsed = diff.count();

		elapsed_arr[repeat] = elapsed;



		//gettimeofday(&endtime,NULL); 
		//elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
		//cout<<"TACO SDDMM K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000 ;
		//cout<< endl;
	}
	//gettimeofday(&endtime,NULL); 

	TYPE sum = 0;
	for (int32_t i = 0; i < C1_dimension; i++) 
	{
		for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++) 
		{
			sum += D_vals[pA2];
		}
	}

	
	cout<<sum<<endl;


	elapsed = findMedian(elapsed_arr,REP);
	cout<<"TACO SDDMM median K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000 ;
		cout<< endl;

	elapsed = elapsed_arr[0];
	cout<<"TACO SDDMM max K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000 ;
		cout<< endl;

	elapsed = elapsed_arr[REP-1];
	cout<<"TACO SDDMM min K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000 ;
		cout<< endl;
	
	//elapsed += ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	//cout<<"TACO SDDMM K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000;
	//cout<< endl;

	return 0;
}

int taco_spmm_papi(taco_tensor_t *C, taco_tensor_t *A, taco_tensor_t *B) 
{
	int C1_dimension = (int)(C->dimensions[0]);
	int C2_dimension = (int)(C->dimensions[1]);
	double* restrict C_vals = (double*)(C->vals);
	int A1_dimension = (int)(A->dimensions[0]);
	int* restrict A2_pos = (int*)(A->indices[1][0]);
	int* restrict A2_crd = (int*)(A->indices[1][1]);
	double* restrict A_vals = (double*)(A->vals);
	int B1_dimension = (int)(B->dimensions[0]);
	int B2_dimension = (int)(B->dimensions[1]);
	double* restrict B_vals = (double*)(B->vals);

	#pragma omp parallel for schedule(static)
	for (int32_t pC = 0; pC < (C1_dimension * C2_dimension); pC++) {
		C_vals[pC] = 0.0;
	}

	gettimeofday(&starttime,NULL); 
	{
		#include "../inc/papi_init_reg.cc"
		for (int32_t i = 0; i < A1_dimension; i++) {
			for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++) {
				int32_t j = A2_crd[pA2];
				for (int32_t k = 0; k < B2_dimension; k++) {
					int32_t pC2 = i * C2_dimension + k;
					int32_t pB2 = j * B2_dimension + k;
					C_vals[pC2] = C_vals[pC2] + A_vals[pA2] * B_vals[pB2];
				}
			}
		}
		#include "../inc/papi_end.cc"
	}
	gettimeofday(&endtime,NULL); 
	{
		#include "../inc/papi_init_lx.cc"
		for (int32_t i = 0; i < A1_dimension; i++) {
			for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++) {
				int32_t j = A2_crd[pA2];
				for (int32_t k = 0; k < B2_dimension; k++) {
					int32_t pC2 = i * C2_dimension + k;
					int32_t pB2 = j * B2_dimension + k;
					C_vals[pC2] = C_vals[pC2] + A_vals[pA2] * B_vals[pB2];
				}
			}
		}
		#include "../inc/papi_end.cc"
	}
	{
		#include "../inc/papi_init_stall.cc"
		for (int32_t i = 0; i < A1_dimension; i++) {
			for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++) {
				int32_t j = A2_crd[pA2];
				for (int32_t k = 0; k < B2_dimension; k++) {
					int32_t pC2 = i * C2_dimension + k;
					int32_t pB2 = j * B2_dimension + k;
					C_vals[pC2] = C_vals[pC2] + A_vals[pA2] * B_vals[pB2];
				}
			}
		}
		#include "../inc/papi_end.cc"
	}
	elapsed += ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	cout<<"TACO SpMM K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000;
	cout<< endl;
	return 0;
}


int taco_sddmm_papi(taco_tensor_t *D, taco_tensor_t *A, taco_tensor_t *B, taco_tensor_t *C)
{
	int D1_dimension = (int)(D->dimensions[0]);
	double* restrict D_vals = (double*)(D->vals);
	int A1_dimension = (int)(A->dimensions[0]);
	int* restrict A2_pos = (int*)(A->indices[1][0]);
	int* restrict A2_crd = (int*)(A->indices[1][1]);
	double* restrict A_vals = (double*)(A->vals);
	int B1_dimension = (int)(B->dimensions[0]);
	int B2_dimension = (int)(B->dimensions[1]);
	double* restrict B_vals = (double*)(B->vals);
	int C1_dimension = (int)(C->dimensions[0]);
	int C2_dimension = (int)(C->dimensions[1]);
	double* restrict C_vals = (double*)(C->vals);

	int32_t pD2 = 0;

	gettimeofday(&starttime,NULL);
	{
		#include "../inc/papi_init_reg.cc"
		for (int32_t i = 0; i < C1_dimension; i++) {
			for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++) {
				int32_t j = A2_crd[pA2];
				D_vals[pD2] = 0.0;
				for (int32_t k = 0; k < C2_dimension; k++) {
					int32_t pB2 = j * B2_dimension + k;
					int32_t pC2 = i * C2_dimension + k;
					D_vals[pD2] = D_vals[pD2] + A_vals[pA2] * B_vals[pB2] * C_vals[pC2];
				}
				pD2++;
			}
		}
		#include "../inc/papi_end.cc"
	}
	gettimeofday(&endtime,NULL); 
	pD2 = 0;
	{
		#include "../inc/papi_init_lx.cc"
		for (int32_t i = 0; i < C1_dimension; i++) {
			for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++) {
				int32_t j = A2_crd[pA2];
				D_vals[pD2] = 0.0;
				for (int32_t k = 0; k < C2_dimension; k++) {
					int32_t pB2 = j * B2_dimension + k;
					int32_t pC2 = i * C2_dimension + k;
					D_vals[pD2] = D_vals[pD2] + A_vals[pA2] * B_vals[pB2] * C_vals[pC2];
				}
				pD2++;
			}
		}
		#include "../inc/papi_end.cc"
	}
	pD2 = 0;
	{
		#include "../inc/papi_init_stall.cc"
		for (int32_t i = 0; i < C1_dimension; i++) {
			for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++) {
				int32_t j = A2_crd[pA2];
				D_vals[pD2] = 0.0;
				for (int32_t k = 0; k < C2_dimension; k++) {
					int32_t pB2 = j * B2_dimension + k;
					int32_t pC2 = i * C2_dimension + k;
					D_vals[pD2] = D_vals[pD2] + A_vals[pA2] * B_vals[pB2] * C_vals[pC2];
				}
				pD2++;
			}
		}
		#include "../inc/papi_end.cc"
	}
	elapsed += ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	cout<<"TACO SDDMM K= "<<C2_dimension<<", time= " << elapsed << ", GFLOPS=  " << (double)2*(double)gold_ne*(double)C2_dimension/elapsed/1000000000;
	cout<< endl;

	return 0;
}


taco_tensor_t* create_dense(int M, int N, int mode) 
{
	taco_tensor_t* C = (taco_tensor_t*) malloc(sizeof(taco_tensor_t));
	int32_t* C_dims = new int32_t[2];
	C -> dimensions = C_dims; 
	C -> dimensions[0] = M;
	C -> dimensions[1] = N;
	C -> indices = NULL;
	int C1_dimension = (int) N;
	int C2_dimension = (int) M;
	double* restrict C_vals = (double*)(C->vals);

	C_vals = (double*)malloc(sizeof(double) * (C1_dimension * C2_dimension));

	

	for (int i =0 ; i < C1_dimension * C2_dimension ; i ++)
	{
		if (mode == 0)
		{
			C_vals[i] = 0.0;
		}
		else if (mode == 1)
		{
			C_vals[i] = ((TYPE) (i%100))/100;
		}
	}

	C->vals = C_vals;

	return C;
}

taco_tensor_t* create_sparse(taco_tensor_t* A, int mode )
{
	taco_tensor_t* C = (taco_tensor_t*) malloc(sizeof(taco_tensor_t));
	C -> dimensions = new int32_t[1];
	int C1_dimension = A->dimensions[0];
	C -> dimensions[0] = C1_dimension;
	C -> vals_size = A->vals_size;
	C -> indices = new uint32_t**[2];
	for (int i =0; i<2; i++)
		C -> indices[i] = new uint32_t*[2];
	int A1_dimension = (int)(A->dimensions[0]);
	int* restrict A2_pos = (int*)(A->indices[1][0]);
	int* restrict A2_crd = (int*)(A->indices[1][1]);

	TYPE* C_vals = (TYPE*) malloc(sizeof(TYPE)*(C->vals_size));
	uint32_t* C_col_ind = (uint32_t*) malloc(sizeof(uint32_t)*(C->vals_size));
	uint32_t* C_row_ptr = (uint32_t*) malloc(sizeof(uint32_t)*((C1_dimension+1))); 

	C_row_ptr[0] = 0;

	for (int32_t i = 0; i < C1_dimension; i++)
	{
		uint32_t nnz_row = 0; 
		for (int32_t pA2 = A2_pos[i]; pA2 < A2_pos[(i + 1)]; pA2++)
		{
			nnz_row ++ ;
			int pos = C_row_ptr[i] + nnz_row;
			C_col_ind[pos] = A2_crd[pA2];
			if (mode == 0)
			{
				C_vals[pos] = 0;
			}
			else if(mode == 1)
			{
				C_vals[pos] = ((TYPE) (i%100))/100;
			}
			else if(mode == 2)
			{
				C_vals[pos] = A->vals[pos];
			}
		}

		C_row_ptr[i+1] = C_row_ptr[i] + nnz_row;
	}

	C->indices[1][0] = (uint32_t*) (C_row_ptr);
	C->indices[1][1] = (uint32_t*) (C_col_ind);
	C->vals = (TYPE *) C_vals;


	return C;
}

taco_tensor_t* read_sparse()
{
	taco_tensor_t* C = (taco_tensor_t*) malloc(sizeof(taco_tensor_t));
	int32_t* C_dim = new int32_t[2];
	C -> dimensions = C_dim; //(int32_t*) malloc(sizeof(int32_t)); //
	int C1_dimension = nr;
	C -> dimensions[0] = C1_dimension;
	C -> vals_size = ne;
	C -> indices = new uint32_t**[2];
	for (int i =0; i<2; i++)
		C -> indices[i] = new uint32_t*[2];
	
	
	TYPE* C_vals = (TYPE*) malloc(sizeof(TYPE)*(C->vals_size));
	uint32_t* C_col_ind = (uint32_t*) malloc(sizeof(uint32_t)*(C->vals_size));
	uint32_t* C_row_ptr = (uint32_t*) malloc(sizeof(uint32_t)*((C1_dimension+1))); 


	int old_row = 0;
	C_row_ptr[0] = 0;
	for (int32_t i = 0; i < ne; i++)
	{
		C_vals[i] = gold_temp_v[i].val;
		C_col_ind[i] = gold_temp_v[i].col;
		while (old_row != gold_temp_v[i].row)
		{
			old_row ++;
			C_row_ptr[old_row] = i;
		}
	}
	while (old_row < C1_dimension)
	{
		old_row++;
		C_row_ptr[old_row] = C_row_ptr[old_row-1];
	}

	C->indices[1][0] = (uint32_t*) (C_row_ptr);
	C->indices[1][1] = (uint32_t*) (C_col_ind);
	C->vals = (TYPE *) C_vals;

	return C;
}

#endif
