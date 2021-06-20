#include "../inc/model.h"
#include "../inc/j_stream.h"


int main(int argc, char **argv)
{
	ready(argc, argv);
	struct timeval starttime, midtime, endtime,timediff;
	double total_time=0.0, avg_time=0.0;
	int Ti,Tk;
	//qsort(gold_temp_v, ne, sizeof(struct v_struct), compare1);

	/*
	init_model_parallel();

	
	pick_tile(Ti,Tk);
	init_model();
	pick_tile(Ti,Tk);
	*/

	cout<<"================================STATS================================="<<endl;

	int kij = 0;	
	int order = 0;
	uint64_t cache = CACHE;
	int Nk = 1000;
	SM_K = Nk;
	Tk = Nk;

	int argcc = 1;
	while(argcc < argc-1)
	{
		argcc ++;
		printf("New option %s\n",argv[argcc]);
		string opt = argv[argcc];
		if(opt.compare("kij")==0)
		{
			argcc ++;
			kij = atoi( argv[argcc]);
		}
		if(opt.compare("order")==0)
		{
			argcc ++;
			order = atoi( argv[argcc]);
		}
		if(opt.compare("cache")==0)
		{
			argcc ++;
			cache = atoi( argv[argcc]);
		}
		if(opt.compare("Nk")==0)
		{
			argcc ++;
			SM_K = Nk = atoi( argv[argcc]);
		}
		if(opt.compare("Tk")==0)
		{
			argcc ++;
			Tk = atoi( argv[argcc]);
			if (Tk > Nk)
				return 1;
		}
	}

	init_model();
	pick_tile(Ti,Tk,Nk,cache);

//	cout<<nr/28<<endl;
	int ti_est = nr/28;

//	cout<<"ti est is "<<ti_est<<endl;
	uint64_t accol = active_cols(ti_est,false);
//	cout<<"ti est is "<<ti_est<<" "<<accol<<endl;



	cache = 35*1024*1024/8/2;
	uint64_t sp_data = 2*ne;
	uint64_t d_data = 0;
	uint64_t nnz = ne;
	uint64_t best_tk = 0;
	
	uint64_t model_min=0;
	for(uint64_t Tkk = 32; Tkk < Nk + Nk/10 ; Tkk *= 2)
	{
		if(Tkk > Nk)
		{
			Tkk = Nk;
		}	
		Ti = ne/nr;
		uint64_t accol_tile = 0;
		do
		{
			accol = active_cols(Ti,false);
			accol_tile  = accol/((nr - 1)/Ti + 1);
			Ti *= 2;
//			cout << Ti << endl;
		}
		while(accol_tile*Tkk < cache && Ti < nr);
//		cout<<"Tk "<<Tkk<<"  ";	
		if( 2*nnz + Tkk*nc < cache)
			sp_data = 2*nnz;
		else
			sp_data = 2*nnz*Nk/Tkk;

		if( Tkk*nc < cache)
			d_data = nc * Nk;
		else
			d_data = accol * Nk/8;
	
		cout<<sp_data + d_data << " "<<sp_data<<" "<<d_data<<endl;	
		if(model_min == 0 || model_min >= sp_data + d_data)
		{
			model_min = sp_data + d_data;
			best_tk = Tkk;
		}
	}

	cout<<"Ti "<<nc<<" Tk "<<best_tk<<endl;
	
	//cout<<"Ti, Tk, cache, "<<Ti << ", "<<Tk<<", "<<cache<<endl;
	/*
	if (Ti*CORE > nr)
	{
		Ti = (nr-1)/CORE +1;
		Tk = MIN(Nk,cache/Ti);
		Tk = Tk - (Tk % 32);
	}
	*/
/*
	SM_WIDTH = Ti;
	gettimeofday(&starttime,NULL);
	gen_structure();
	gettimeofday(&endtime,NULL);
	double elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
//	cout<<"BCSC Preprocess overhead "<<" Elapsed: "<<elapsed<<endl;
	spmm(Nk,Tk,cache,kij,order);
	cout<<"================================STATS================================="<<endl;
*/
}
