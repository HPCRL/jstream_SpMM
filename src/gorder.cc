#ifndef GORDER_CC
#define GORDER_CC

#include "../inc/gorder.h"


int gorder(int* rowPtr, int* colInd,  int window, int n, int m, int* result)
{
	int* count = new int[m];
	int* used = new int[n];


	#pragma omp parallel
	for(int i=0; i<m; i++)
	{
		count[i] = 0;
	}


	#pragma omp parallel
	for(int i=0; i<n; i++)
	{
		used[i] = 0;
	}


	srand (time(NULL));

	result[0] = 0; 
	// result[0] = rand()%n; // try this too
	int seed = 0;

/*

	for (int i=0; i<n; i++)
		for(int j=rowPtr[i] ; j<rowPtr[i+1] ; j++)
			if(colInd[j]<0 || colInd[j]>=m)
				printf("Error in CSR\n");
*/
	// we got a seed. Now pick the best matches for a given window
	while ( seed < n - 1)
	{
//		if (result[seed] > n || result[seed]<0) printf("Error here. seed  %d result[seed] %d n is %d\n",seed,result[seed],n);
		used[result[seed]] = 1;
		#pragma omp parallel for
		for (int i = rowPtr[result[seed]] ; i < rowPtr[result[seed]+1] ; i++)
		{
//			if (colInd[i] > m) printf("Error here\n");
			count[colInd[i]]++;
		}
		if (seed >= window)
		{
			#pragma omp parallel for
			for (int i = rowPtr[result[seed - window]] ; i < rowPtr[result[seed - window]+1] ; i++)
				count[colInd[i]]--;
		}
		int max = 0;
		int max_id = -1;
		for(int i=n-1; i>=0 ; i--)
		{
			int match = 0;
			if(used[i] == 0)
			{
				for (int j = rowPtr[result[i]] ; j < rowPtr[result[i]+1] ; j++)
					match += count[colInd[j]];
				if (match >= max)
				{
					max = match;
					max_id = i;
				}
			}
		}
		seed++;
		result[seed] = max_id;
	}

	delete [] count;
	delete [] used;

	return 0;
}

#endif