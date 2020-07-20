#ifndef PREPROCESSOR_CC
#define PREPROCESSOR_CC
#define rem_mm(X) if(X!= NULL) {_mm_free(X); X=NULL;}
#define rem(X) if(X!= NULL) {free(X); X=NULL;}

#include "../inc/preprocessor.h"


int SM_K = 8;
int SM_WIDTH = 64;
bool HIGH_SD = false;

int* pb = NULL;
int* didx = NULL;
int* d = NULL;
TYPE* dv = NULL;
int* dc = NULL;


void gen_structure()
{
	rem_mm(pb);
	rem_mm(didx);
	rem_mm(d);
	rem_mm(dv);
	rem_mm(dc);


	int i,j;
	
	nseg = CEIL(nr, SM_WIDTH);
	d = (int *)_mm_malloc(sizeof(int)*(ne+128), 64); // can be char
	dv = (TYPE *)_mm_malloc(sizeof(TYPE)*(ne+128), 64); // can be char
	int * new_pb = (int *)_mm_malloc(sizeof(int)*(nseg+1), 64);
	pb = (int *)_mm_malloc(sizeof(int)*(nseg+1), 64);
	int * new_pb3 = (int *)_mm_malloc(sizeof(int)*(nseg+1), 64);
	
	if(d == NULL || dv == NULL || new_pb == NULL || pb == NULL || new_pb3 == NULL)
	{
		
		printf("Memory allocation problem in preprocessor\n");
		exit(1);
	}

	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	int m = nr;
	int k = nc;

	double mean = 0;
	double sd = 0;
	//if(SM_WIDTH == 256)
	//cout<<"Average ele in rows : "<<ne/nc<<endl;

	
	mean = ne/nc;
	//if(SM_WIDTH == 256)
	//cout<<"\n Mean : "<<mean<<endl;
	int temp_colx=temp_v[0].col;
	double sumOfSquare = 0;
	int rows_in_column = 0;
	int count=0;
	/*
	for (int i = 0; i < ne; i++) 
	{
			
			if(temp_colx != temp_v[i].col){
				double square_dist = (count ) - mean ;
				if (square_dist < 0)
					square_dist = square_dist * (-1);

				sumOfSquare += pow(square_dist,2);
				count=1;
				temp_colx=temp_v[i].col;
			}
			else
			{
				count++;
			}
			
	}

	sd = sumOfSquare/nc;
	sd = sqrt(sd);

	*/
	//if(SM_WIDTH == 256)
	//cout<<" \nstandard deviation = "<< sd<<endl;
	/*
	if(sd > 200)
	{
		SM_WIDTH = 16;
		HIGH_SD = true;
	}
	*/
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	//change
	// cout<<"nseg "<<nseg<<endl;
	// pb = (int *)_mm_malloc(sizeof(int)*(nseg+1), 64);
	// memset(pb, -1, sizeof(int)*(nseg+1));
	//    dc = (int *)malloc(sizeof(int)*ne); // loose UB
	//    didx = (int *)malloc(sizeof(int)*ne); // loose UB

	//#pragma vector aligned
	struct timeval starttime, midtime, endtime,timediff;
	double elapsed = 0;
	gettimeofday(&starttime,NULL);

	#pragma omp parallel for
	for(i=0;i<ne;i++) {
		temp_v[i].grp = temp_v[i].row / SM_WIDTH;//change      
	}
	
	
	

	//#pragma vector aligned
	#pragma omp parallel for 
	for(int i=0;i<nseg+1;i++ )
		new_pb[i]=0;
	
	int temp_panel = -2;
	int temp_col = -2;
	for(int i=0;i < ne ; i++)
	{
		if(temp_v[i].grp != temp_panel || temp_v[i].col != temp_col)
		{
			new_pb[temp_v[i].grp+1]++;
			temp_panel=temp_v[i].grp;
			temp_col = temp_v[i].col;
		}
	}

	for(int i=0;i<nseg;i++ )
	{
		new_pb[i+1]+=new_pb[i];
		
	}
	
	pb[0]=0;
	
	//#pragma vector aligned
	#pragma omp parallel for num_threads(CORE)
	for(int i=0;i<nseg+1;i++ )
	{
		pb[i]=new_pb[i];
		new_pb3[i]=new_pb[i];
		
	}
	dcnt = new_pb[nseg];
	// cout<<"\n dcnt x = "<<dcnt<<endl;

	gettimeofday(&endtime,NULL);
	elapsed += ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;


	dc = (int *)_mm_malloc(sizeof(int)*(dcnt+128), 64); // loose UB
	didx = (int *)_mm_malloc(sizeof(int)*(dcnt+128), 64); // loose UB
	int * new_didx2 = (int *)_mm_malloc(sizeof(int)*(dcnt+128), 64);

	gettimeofday(&starttime,NULL);
	if (dc == NULL || didx == NULL || new_didx2 == NULL)
	{
		
		printf("Memory allocation problem in preprocessor\n");
		exit(1);
	}


	//#pragma vector aligned
	#pragma omp parallel for
	for(int i=0;i<dcnt+128;i++)
	{
		dc[i]=0;
		didx[i]=0;
	}
	
	temp_panel = -2; 
	temp_col = -2;
	didx[0]=-1;
	//#pragma vector aligned
	for(int i=0;i < ne ; i++)
	{
		if(temp_v[i].grp != temp_panel || temp_v[i].col != temp_col)
		{
			dc[new_pb[temp_v[i].grp]] = temp_v[i].col;
			didx[new_pb[temp_v[i].grp]]++;
			new_pb[temp_v[i].grp]++;
			temp_panel=temp_v[i].grp;
			temp_col = temp_v[i].col;
			
		}
		else
		{
			didx[new_pb[temp_v[i].grp]]++;
		}
	}

	for(int i=0;i<dcnt+127;i++ )
		didx[i+1]+=didx[i];
	didx[dcnt]++; 

	//#pragma vector aligned
	#pragma omp parallel for 
	for(int i=0;i<dcnt+128;i++)
	{
		new_didx2[i]=didx[i];

	}

	
	
	temp_panel = pb[temp_v[0].grp];
	temp_col = temp_v[0].col;
	
	//#pragma vector aligned
	#pragma omp parallel for num_threads(CORE)
	for(int i=0;i<ne;i++)
	{
		d[i]=0;
		dv[i]=0.0;
	}

	int range = nseg/CORE;
	range++;

	#pragma vector aligned
	#pragma omp parallel num_threads(CORE)
	{
		int threadID = omp_get_thread_num();
		for(int i=0; i< ne ; i++)
		{
			if( ( temp_v[i].grp >= range*threadID ) &&  ( temp_v[i].grp < range*(threadID+1) ) )
			{
				int didx_index = new_pb3[temp_v[i].grp] ;
				int row = new_didx2[didx_index];
				d[row] = temp_v[i].row;
				dv[row] = temp_v[i].val;
				new_didx2[didx_index]+=1;
				gold_temp_v[i].grp = row;
				if(new_didx2[didx_index] == new_didx2[didx_index+1])
					new_pb3[temp_v[i].grp]+=1;
			}
		}
	}
	gettimeofday(&endtime,NULL);
	elapsed += ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
//	cout<<"BCSC Preprocess overhead w/o memory operations "<<" Elapsed: "<<elapsed<<endl;

	_mm_free(new_pb);
	_mm_free(new_pb3);
	_mm_free(new_didx2);

}

#endif
