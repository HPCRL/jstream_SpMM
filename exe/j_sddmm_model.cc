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
	int cache = CACHE;
	int Nk = 128;
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
	pick_tile_sddmm(Ti,Tk,Nk,cache);

	//cout<<"Ti, Tk, cache, "<<Ti << ", "<<Tk<<", "<<cache<<endl;

	if (Ti*CORE > nr)
	{
		Ti = (nr-1)/CORE +1;
		Tk = MIN(Nk,cache/Ti);
		Tk = Tk - (Tk % 32);
	}

	SM_WIDTH = Ti;
	gettimeofday(&starttime,NULL);
	gen_structure();
	gettimeofday(&endtime,NULL);
	double elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
//	cout<<"BCSC Preprocess overhead "<<" Elapsed: "<<elapsed<<endl;
	sddmm(Nk,Tk,cache,kij,order);
	cout<<"================================STATS================================="<<endl;
}
