#ifndef J_STREAM_CC
#define J_STREAM_CC
#include "../inc/j_stream.h"
#include <algorithm>

//double **X =(double**) malloc(sizeof(double*)), **Y=(double**) malloc(sizeof(double*));

//int initFlush = 0;

TYPE *B=NULL, *C=NULL, *D=NULL;
TYPE *vout_gold = NULL;

/*

void cacheFlush(double **X, double **Y) 
{
	int size_flush = 8000000;
	if (initFlush == 0)
	{
		//printf("Init X\n");
		X[0] = (double*)malloc(size_flush*sizeof(double));
		Y[0] = (double*)malloc(size_flush*sizeof(double));

		

		initFlush = 1;
	}
	double * XX = X[0];
	double * YY = Y[0];

	for(int i=0; i<size_flush; i++) 
	{
		XX[i]=YY[i]+rand() % 5;
		YY[i] += rand() % 7;
	}
	//printf("Finished FLushing\n");
}
*/
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

void spmmkernel(int *segment, int *column_number, int *column_index, int *row_number, TYPE *A_value, TYPE* B, TYPE* C ,int tile_size, int K)
{
	int Nk = K;
	#include "../inc/j_stream_spmm.cc"
}

void sddmmkernel(int *segment, int *column_number, int *column_index, int *row_number, TYPE *A_value, TYPE* B, TYPE* C , TYPE* D, int tile_size, int K)
{
	#include "../inc/j_stream_sddmm.cc"
}


void spmm(int Nk, int Tk, int CacheSize, int kij, int order)
{
	//gettimeofday(&starttime,NULL);
	#include "../inc/init_test.cc"
	vout_gold = (TYPE *)malloc(sizeof(TYPE)*nr*(sc+PADC));
	
	#include "../inc/malloc_check.cc"
	//gettimeofday(&endtime,NULL); 
	//elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	               
	//cout<<"Initialization time = " << elapsed <<endl;

	elapsed = 0;

	// Warm-up
	#include "../inc/j_stream_spmm.cc"

	double *elapsed_arr = new double[REP];

	double min = -1, max = -1;
	if(kij)
	{

		#define KIJ
		for(int rep = 0; rep<REP ; rep++)
		{
			cacheFlush(X,Y);
			//gettimeofday(&starttime,NULL); 
			//#include "../inc/j_stream_spmm.cc"
			uint64_t bjump = sizeB / layers * (rep % layers);
			uint64_t cjump = sizeC / layers * (rep % layers);
			auto start = std::chrono::high_resolution_clock::now();
			spmmkernel(segment, column_number, column_index, row_number, A_value, B + bjump, C + cjump, tile_size, K);
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			elapsed = diff.count();
			elapsed_arr[rep] = elapsed;
			//gettimeofday(&endtime,NULL); 
			//elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;

			//cout<<"J Stream SpMM KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
			//cout<< endl;

		}
		#undef KIJ
	}
	else
	{

		for(int rep = 0; rep<REP ; rep++)
		{
			cacheFlush(X,Y);
			uint64_t bjump = sizeB / layers * (rep % layers);
			uint64_t cjump = sizeC / layers * (rep % layers);
			//gettimeofday(&starttime,NULL); 
			//#include "../inc/j_stream_spmm.cc"
			auto start = std::chrono::high_resolution_clock::now();
			spmmkernel(segment, column_number, column_index, row_number, A_value, B + bjump, C + cjump, tile_size, K);
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			elapsed = diff.count();
			elapsed_arr[rep] = elapsed;
			//gettimeofday(&endtime,NULL); 
			//elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
			//cout<<"J Stream SpMM KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
			//cout<< endl;
		}
	}


	elapsed /= REP;
		               
	

	gettimeofday(&starttime,NULL); 
	#include "../inc/correctness.cc"
	gettimeofday(&endtime,NULL); 
	//elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	
	elapsed = findMedian(elapsed_arr,REP);
	cout<<"J Stream SpMM median KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
	cout<<endl;
	elapsed = elapsed_arr[0];
	cout<<"J Stream SpMM max KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
	cout<<endl;
	elapsed = elapsed_arr[REP-1];
	cout<<"J Stream SpMM min KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
	cout<<endl;
	//cout<<"Correctness check time Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000 << endl;

	#include "../inc/end_test.cc"
	
	//vout_gold = NULL;

}

void spmm_papi(int Nk, int Tk, int CacheSize, int kij, int order)
{
	#define SINGLE
	#include "../inc/init_test.cc"
	vout_gold = (TYPE *)malloc(sizeof(TYPE)*nr*(sc+PADC));
	#include "../inc/malloc_check.cc"

	if(kij)
	{
		#define KIJ
		{
		
			#include "../inc/papi_init_reg.cc"
			#include "../inc/j_stream_spmm.cc"	
			#include "../inc/papi_end.cc"
		
		}
		{
		
			#include "../inc/papi_init_lx.cc"
			#include "../inc/j_stream_spmm.cc"	
			#include "../inc/papi_end.cc"
		
		}
		{

			#include "../inc/papi_init_stall.cc"
			#include "../inc/j_stream_spmm.cc"	
			#include "../inc/papi_end.cc"
		
		}		
		#undef KIJ	
	}
	else
	{
		{
		
			#include "../inc/papi_init_reg.cc"
			#include "../inc/j_stream_spmm.cc"	
			#include "../inc/papi_end.cc"
		
		}
		{
		
			#include "../inc/papi_init_lx.cc"
			#include "../inc/j_stream_spmm.cc"	
			#include "../inc/papi_end.cc"
		
		}
		{

			#include "../inc/papi_init_stall.cc"
			#include "../inc/j_stream_spmm.cc"	
			#include "../inc/papi_end.cc"
		
		}		
	}


	//gettimeofday(&starttime,NULL);
	
	auto start = std::chrono::high_resolution_clock::now();
	#include "../inc/j_stream_spmm.cc"	
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end-start;
	elapsed = diff.count();
	
	//gettimeofday(&endtime,NULL); 
	//elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	cout<<"J Stream SpMM KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
	cout<< endl;

	#include "../inc/correctness.cc"

	#include "../inc/end_test.cc"
	#undef SINGLE
}


void sddmm(int Nk, int Tk, int CacheSize, int kij, int order)
{
	#include "../inc/init_test.cc"
	vout_gold = (TYPE *)malloc(sizeof(TYPE)*ne);

	D = (TYPE *)_mm_malloc(sizeof(TYPE) * ne ,64);
	for(i = 0;i<ne;i++)
	{
		D[i] = 0.0;
	}
	for(i = 0; i < nr * (Nk+PADC) ; i++)
	{
		C[i] = (TYPE) (rand()%10)/10;
	}
	#include "../inc/malloc_check.cc"
	if (D == NULL)
	{
		printf("Malloc for D failed exiting...\n");
		exit(1);
	}
	
	double *elapsed_arr = new double[REP];

	elapsed = 0;

	// Warm-up
	#include "../inc/j_stream_sddmm.cc"

	if(kij)
	{
		#define KIJ
		for(int rep = 0; rep<REP ; rep++)
		{
			cacheFlush(X,Y);
			//gettimeofday(&starttime,NULL); 
			//#include "../inc/j_stream_sddmm.cc"
			uint64_t bjump = sizeB / layers * (rep % layers);
			uint64_t cjump = sizeC / layers * (rep % layers);
			auto start = std::chrono::high_resolution_clock::now();
			sddmmkernel(segment, column_number, column_index, row_number, A_value, B + bjump, C + cjump, D, tile_size, K);
			//gettimeofday(&endtime,NULL); 
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			elapsed = diff.count();
			elapsed_arr[rep] = elapsed;
			//elapsed += ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
		}
		#undef KIJ
	}
	else
	{

		for(int rep = 0; rep<REP ; rep++)
		{
			cacheFlush(X,Y);
			//gettimeofday(&starttime,NULL); 
			//#include "../inc/j_stream_sddmm.cc"
			uint64_t bjump = sizeB / layers * (rep % layers);
			uint64_t cjump = sizeC / layers * (rep % layers);
			auto start = std::chrono::high_resolution_clock::now();
			sddmmkernel(segment, column_number, column_index, row_number, A_value, B + bjump, C + cjump, D, tile_size, K);
			//gettimeofday(&endtime,NULL); 
			auto end = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> diff = end-start;
			elapsed = diff.count();
			elapsed_arr[rep] = elapsed;
			//elapsed += ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
		}
	}
	elapsed /= REP;


	elapsed = findMedian(elapsed_arr,REP);
	cout<<"J Stream SDDMM median KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
	cout<<endl;
	elapsed = elapsed_arr[0];
	cout<<"J Stream SDDMM max KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
	cout<<endl;
	elapsed = elapsed_arr[REP-1];
	cout<<"J Stream SDDMM min KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
	cout<<endl;

	#include "../inc/correctness_sddmm.cc"
	               
	//cout<<"J Stream SDDMM KIJ "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
	//cout<< endl;

	#include "../inc/end_test.cc"
}


void sddmm_papi(int Nk, int Tk, int CacheSize, int kij, int order)
{
	#define SINGLE
	#include "../inc/init_test.cc"
	vout_gold = (TYPE *)malloc(sizeof(TYPE)*ne);

	D = (TYPE *)_mm_malloc(sizeof(TYPE) * ne ,64);
	for(i = 0;i<ne;i++)
	{
		D[i] = 0.0;
	}
	for(i = 0; i < nr * (Nk+PADC) ; i++)
	{
		C[i] = (TYPE) (rand()%10)/10;
	}
	#include "../inc/malloc_check.cc"
	if (D == NULL)
	{
		printf("Malloc for D failed exiting...\n");
		exit(1);
	}
	
	if(kij)
	{
		#define KIJ
		{
			#include "../inc/papi_init_reg.cc"
			#include "../inc/j_stream_sddmm.cc"	
			#include "../inc/papi_end.cc"	
		}
		{
			#include "../inc/papi_init_lx.cc"
			#include "../inc/j_stream_sddmm.cc"	
			#include "../inc/papi_end.cc"	
		}
		{
			#include "../inc/papi_init_stall.cc"
			#include "../inc/j_stream_sddmm.cc"	
			#include "../inc/papi_end.cc"	
		}
		#undef KIJ
	}
	else
	{
		
		{
			#include "../inc/papi_init_reg.cc"
			#include "../inc/j_stream_sddmm.cc"	
			#include "../inc/papi_end.cc"	
		}
		{
			#include "../inc/papi_init_lx.cc"
			#include "../inc/j_stream_sddmm.cc"	
			#include "../inc/papi_end.cc"	
		}
		{
			#include "../inc/papi_init_stall.cc"
			#include "../inc/j_stream_sddmm.cc"	
			#include "../inc/papi_end.cc"	
		}	
		
	}


	//gettimeofday(&starttime,NULL);
	
	auto start = std::chrono::high_resolution_clock::now();
	#include "../inc/j_stream_sddmm.cc"	
	//gettimeofday(&starttime,NULL); 
	//#include "../inc/j_stream_sddmm.cc"
	//gettimeofday(&endtime,NULL); 
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end-start;
	elapsed = diff.count();
	//elapsed += ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;

	
	//gettimeofday(&endtime,NULL); 
	//elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	cout<<"J Stream SDDMM "<<kij<<" order "<<order<<" K= "<< K <<" Ti="<<SM_WIDTH<< ", Tk= "<< tile_size << ", Elapsed= " << elapsed << " sec, GFLOPS=  " << (double)2*(double)gold_ne*(double)K/elapsed/1000000000;
	cout<< endl;

	#include "../inc/correctness_sddmm.cc"


	#include "../inc/end_test.cc"
	#undef SINGLE
}



#endif