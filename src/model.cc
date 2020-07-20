#ifndef MODEL_CC
#define MODEL_CC
#include "../inc/model.h"
//#define LOGC(X) ceil(log2(X))
/*
#define LOGC(X) logc(X)
#define REVERSE_LOGC(X) pow(2,X)
#define AVERAGE_BIN(X) X*3/4
*/

#define FACTOR 1
#define LOGC(X) X/FACTOR
#define REVERSE_LOGC(X) X*FACTOR
#define AVERAGE_BIN(X) (2*X)/2



uint64_t ** intervalsp = new uint64_t*[1];
uint64_t ** notsmallerp = new uint64_t*[1];
uint64_t ** notsmallerweightedp = new uint64_t*[1];



unsigned int logc(unsigned int n)
{
	if (n == 0)
	return 0;

	unsigned int msb = 0;
	while (n != 0) 
	{
		n = n >> 1;
		msb++;
	}

	return (msb);
}


int init_model(bool print)
{
	uint64_t* intervals = new uint64_t[nr + 1];
	uint64_t* notsmaller = new uint64_t[nr + 1];
	uint64_t* notsmallerweighted = new uint64_t[nr+1];
	for(uint64_t i = 0; i< nc+1 ; i++)
	{
		intervals [i] = 0;
		notsmaller [i] = 0;
		notsmallerweighted[i] = 0;
	}
	struct timeval starttime, midtime, endtime,timediff;
	double elapsed=0;

	gettimeofday(&starttime,NULL); 
	uint64_t i=0;
	for (uint64_t j=0; j<nc ; j++)
	{
		uint64_t lasti = -1;
		while(gold_temp_v[i].col == j)
		{
			uint64_t d = gold_temp_v[i].row - lasti - 1;
			intervals[d] += 1;
			lasti = gold_temp_v[i].row;
			i++;
		}
		uint64_t d = nr - lasti - 1;
		intervals[d] += 1;
		
	}

	gettimeofday(&midtime,NULL); 
	elapsed = ((midtime.tv_sec-starttime.tv_sec)*1000000 + midtime.tv_usec-starttime.tv_usec)/1000000.0;
	if(print)
		cout<<"Preprocessing midpoint time: " << elapsed <<endl;

	if(print)
		cout << ne << " "<<i<<endl;

	notsmaller[nr] = intervals[nr];
	notsmallerweighted[nr] = nr*intervals[nr];

	for(uint64_t Ti = nc-1; Ti > 0 ; Ti--)
	{
		notsmaller [Ti] = notsmaller [Ti + 1] + intervals[Ti];
		notsmallerweighted[Ti] += notsmallerweighted [Ti + 1] + Ti*intervals[Ti];
	}

	/*
	for(uint64_t Ti = 1; Ti <= nc ; Ti++)
	{
		notsmaller [Ti] = notsmaller [Ti -1] + intervals[Ti];
		notsmallerweighted[Ti] += notsmallerweighted [Ti - 1] + Ti*intervals[Ti];
	}
	*/


	gettimeofday(&endtime,NULL); 
	elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	cout<<"Preprocessing for the model time: " << elapsed <<endl;
	*intervalsp = intervals;
	*notsmallerp = notsmaller;
	*notsmallerweightedp = notsmallerweighted;

	return 0;
}

int init_model_parallel(bool print)
{
	uint64_t* intervals = new uint64_t[nr + 1];
	uint64_t* notsmaller = new uint64_t[nr + 1];
	uint64_t* notsmallerweighted = new uint64_t[nr+1];
	double* d_notsmaller = new double[nr + 1];
	double* d_notsmallerweighted = new double[nr+1];


	int lbin_size = LOGC(nc)+1;
	int core = omp_get_max_threads();
	omp_set_num_threads(core);
	uint64_t* log_intervals = new uint64_t[core*lbin_size];
	//cout<<"log func returned "<<lbin_size<<endl;

	for(uint64_t i = 0; i< nc+1 ; i++)
	{
		intervals [i] = 0;
		notsmaller [i] = 0;
		notsmallerweighted[i] = 0;
	}

	for(uint64_t i = 0; i< core*lbin_size ; i++)
	{
		log_intervals[i] = 0;
	}


	struct timeval starttime, midtime, endtime,timediff;
	double elapsed;
/*
	int* s_pos = new int[core];
	int* e_pos = new int[core];
	for(int i =0; i<core; i++)
	{
		start_pos = 
	}
*/


	gettimeofday(&starttime,NULL); 
	
	#pragma omp parallel
	{

		int thread_id = omp_get_thread_num();
		int start_pos= (thread_id * ne)/core;
		int end_pos =  ((thread_id+1) * ne)/core; 

		uint64_t lasti = -1;
		if(start_pos > 0 && gold_temp_v[start_pos].col == gold_temp_v[start_pos-1].col)
			lasti = gold_temp_v[start_pos-1].row;
		int j = gold_temp_v[start_pos].col;
		//int i = start_pos;
		for(int i=start_pos ; i<end_pos ; i++)
		{
			//uint64_t lasti = -1;
			while(gold_temp_v[i].col == j && i<end_pos)
			{
				uint64_t dist = gold_temp_v[i].row - lasti - 1;
				uint64_t d =  LOGC(dist);
				log_intervals[lbin_size*thread_id + d] += 1;
				lasti = gold_temp_v[i].row;
				i++;
			}
			if(i<end_pos)
			{
				uint64_t dist = nr - lasti - 1;
				uint64_t d = LOGC(dist);
				log_intervals[lbin_size*thread_id + d] += 1;
				j++;
				lasti = -1;
			}	
		}

	}

	// REDUCE the sums can be parallelized later on

	for (int th = 1; th < core ; th++)
	{
		#pragma novector
		for (int i = 0; i<lbin_size ; i++)
		{
			log_intervals[i] += log_intervals[th*lbin_size + i];
		}
	}

	log_intervals[lbin_size] = 0; 
	for (int i = 0; i < lbin_size ; i++)
	{
		log_intervals[lbin_size] += log_intervals[i];
		//cout<< "log_intervals["<<i<<"]: "<<log_intervals[i]<<endl;
	}

	if(print)
		cout<<ne<<" "<<log_intervals[lbin_size]<<endl;


	gettimeofday(&midtime,NULL); 
	elapsed = ((midtime.tv_sec-starttime.tv_sec)*1000000 + midtime.tv_usec-starttime.tv_usec)/1000000.0;
	if(print)
		cout<<"Preprocessing midpoint time: " << elapsed <<" "<<lbin_size<<endl;

	notsmaller[lbin_size] = 0;// log_intervals[lbin_size];
	notsmallerweighted[lbin_size] = 0;// log_intervals[lbin_size];

	for(int lTi=lbin_size-1; lTi>=0 ; lTi--)
	{
		uint64_t Ti = (uint64_t)  AVERAGE_BIN(REVERSE_LOGC(lTi));
		notsmaller [lTi] = notsmaller [lTi + 1] + log_intervals[lTi];
		notsmallerweighted[lTi] = notsmallerweighted [lTi + 1] + Ti*log_intervals[lTi];
	}

	for(int Ti = nc; Ti >= 0 ; Ti --)
	{
		int lTi = LOGC(Ti);
		notsmaller[Ti] = notsmaller[lTi];
		notsmallerweighted[Ti] = notsmallerweighted[lTi];
	}

	gettimeofday(&endtime,NULL); 



	elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	cout<<"Preprocessing for the model time: " << elapsed <<endl;
	*intervalsp = intervals;
	*notsmallerp = notsmaller;
	*notsmallerweightedp = notsmallerweighted;

	return 0;
}


uint64_t active_cols(int Ti, bool print)
{
	if (intervalsp == NULL || notsmallerp == NULL || notsmallerweightedp == NULL)
	{
		init_model();
	}
	uint64_t* intervals = *intervalsp;
	uint64_t* notsmaller = *notsmallerp;
	uint64_t* notsmallerweighted = *notsmallerweightedp;
	if(Ti > nc)
		Ti = nc;


	uint64_t div = Ti;
	uint64_t res = notsmallerweighted[Ti] - ((Ti - 1)*notsmaller[Ti]);
	uint64_t all_tiles = (nr-Ti+1);
	all_tiles *= nc;

	bool overflow = (notsmallerweighted[Ti] < ((Ti - 1)*notsmaller[Ti])) || (all_tiles < res);

	if (overflow)
		res = 0;
	else
		res = all_tiles - res;
	
	if (print)
		cout << "raw res "<< res  <<" overflow "<< (overflow ? "true" : "false") 
			<< " "<< nr<< " nsw " << notsmallerweighted[Ti] << " ns " << notsmaller[Ti] << " intervals[Ti] " << intervals[Ti] << " " << (nr-Ti+1)*nc <<" " << all_tiles<< endl;

	//res /= Ti;
	double ratio = ((double) (res)/all_tiles);  //((nr-1)/Ti+1)*nc - res;
	res = ((nr-1)/Ti+1);
	res *= nc;
	res *= ratio;

	if(print)
		cout <<"Ti: "<<Ti << ", model_ac_cols: "<< res <<endl;
	return res;
}

int pick_tile(int &Ti, int &Tk , int Nk , int cache)
{
	uint64_t lowest = -1;
	//for(int k = 32; k <= 128 ; k*=2)
	{
		for (int i = 2; i< 2*nc /*CACHE/k */ && cache/i > 0; i*=2)
		{
			uint64_t k = MIN(Nk,cache/i);
			k = k - (k %32);
			if (k <= 0)
				continue;
			uint64_t test =  active_cols(i,false);

			uint64_t model = ((uint64_t )ne )* ((uint64_t) Nk) / k *5 / 2; // sparse part
			model += test * ((uint64_t) Nk) ; // input part
			model += 2* ((uint64_t) nr)*((uint64_t) Nk); // output part
//			cout<<i<<"\t"<<test<<'\t'<<model/8<<endl;
			if (lowest == -1 || model < lowest)
			{
				lowest = model;
				Ti = i;
				Tk = k;
			}

		}
	}

	return 0;
}

int pick_tile_sddmm(int &Ti, int &Tk , int Nk , int cache)
{
	uint64_t lowest = -1;
	//for(int k = 32; k <= 128 ; k*=2)
	{
		for (int i = 2; i< 2*nc /*CACHE/k */ && cache/i > 0; i*=2)
		{
			uint64_t k = MIN(Nk,cache/i);
			k = k - (k %32);
			if (k <= 0)
				continue;
			uint64_t test =  active_cols(i,false);

			uint64_t model = ((uint64_t )ne )* ((uint64_t) Nk) / k *7 / 2; // sparse part
			model += test * ((uint64_t) Nk) ; // input part
			model += 2* ((uint64_t) nr)*((uint64_t) Nk); // output part
//			cout<<i<<"\t"<<test<<'\t'<<model/8<<endl;
			if (lowest == -1 || model < lowest)
			{
				lowest = model;
				Ti = i;
				Tk = k;
			}

		}
	}

	return 0;
}


#endif
