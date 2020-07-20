#ifndef CACHEFLUSH_CC
#define CACHEFLUSH_CC

#include "../inc/cacheflush.h"
#include <algorithm>

double **X =(double**) malloc(sizeof(double*)), **Y=(double**) malloc(sizeof(double*));

int initFlush = 0;
double *Amkl, *Bmkl, *Cmkl;
int m, n, k;
double alpha, beta;

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


void cacheFlush(double **X, double **Y) 
{
	int size_flush = 6*1000000;

	int mode = 0;




	if (initFlush == 0)
	{
		//printf("Init X\n");
		X[0] = (double*)malloc(size_flush*sizeof(double));
		Y[0] = (double*)malloc(size_flush*sizeof(double));


		m = 2000, k = 200, n = 1000;
		Amkl = (double *)mkl_malloc( m*k*sizeof( double ), 64 );
		Bmkl = (double *)mkl_malloc( k*n*sizeof( double ), 64 );
		Cmkl = (double *)mkl_malloc( m*n*sizeof( double ), 64 );
		

		initFlush = 1;
	}
	double * XX = X[0];
	double * YY = Y[0];
	double total_sum = 0;


	if(mode > 0)
	{
		for(int i=0; i<size_flush; i++) 
		{
			XX[i] = 0.0;
			YY[i] = 0.0;
		}

		#pragma omp parallel for 
		for(int i=0; i<size_flush; i++) 
		{
			XX[i] = (double) (i % 3213123);
			YY[i] += (double) (i % 323);
		}

		for(int i=0; i<size_flush; i++) 
		{
			total_sum = XX[i] + YY[i];
		}
		

		printf("flush %lf\n",total_sum);
		//printf("Finished FLushing\n");
	}




	// WARM Up
	    


    /*
    printf ("\n This example computes real matrix C=alpha*A*B+beta*C using \n"
            " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
            " alpha and beta are double precision scalars\n\n");
	*/

	if(mode > 1)
	{
		m = 2000, k = 200, n = 1000;

		int i,j;

		//printf (" Initializing data for matrix multiplication C=A*B for matrix \n"
		//        " A(%ix%i) and matrix B(%ix%i)\n\n", m, k, k, n);
		alpha = 1.0; beta = 1.0;

		//printf (" Allocating memory for matrices aligned on 64-byte boundary for better \n"
		//        " performance \n\n");

		if (Amkl == NULL || Bmkl == NULL || Cmkl == NULL) {
		  printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
		  mkl_free(Amkl);
		  mkl_free(Bmkl);
		  mkl_free(Cmkl);
		  return ;
		}

		//printf (" Intializing matrix data \n\n");
		for (i = 0; i < (m*k); i++) {
		    Amkl[i] = (double)(i+1);
		}

		for (i = 0; i < (k*n); i++) {
		    Bmkl[i] = (double)(-i-1);
		}

		for (i = 0; i < (m*n); i++) {
		    Cmkl[i] = 0.0;
		}

		//printf (" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");
		for(i = 0; i < 10 ; i++)
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
		            m, n, k, alpha, Amkl, k, Bmkl, n, beta, Cmkl, n);
	}
	//printf ("\n Computations completed.\n\n");

	/*
	printf (" Top left corner of matrix A: \n");
	for (i=0; i<min(m,6); i++) {
	  for (j=0; j<min(k,6); j++) {
	    printf ("%12.0f", A[j+i*k]);
	  }
	  printf ("\n");
	}

	printf ("\n Top left corner of matrix B: \n");
	for (i=0; i<min(k,6); i++) {
	  for (j=0; j<min(n,6); j++) {
	    printf ("%12.0f", B[j+i*n]);
	  }
	  printf ("\n");
	}

	printf ("\n Top left corner of matrix C: \n");
	for (i=0; i<min(m,6); i++) {
	  for (j=0; j<min(n,6); j++) {
	    printf ("%12.5G", C[j+i*n]);
	  }
	  printf ("\n");
	}
	*/
	//  printf ("\n Deallocating memory \n\n");

	/*
	mkl_free(A);
	mkl_free(B);
	mkl_free(C);
	*/
	//printf (" Example completed. \n\n");


    //return 0;


}
#endif