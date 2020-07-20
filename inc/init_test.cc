int i,j;
double tot_ms;
// cout<<endl<<"DENSE SPMM"<<endl;
sc = Nk;

struct timeval starttime, midtime, endtime,timediff;
double elapsed;

int tile_size = Tk;
//        if(sc >=1024)
//            tile_size = 512;   
//cacheFlush(X, Y);
uint64_t layers = 16;
uint64_t mem_capacity = 4000000000; // 4 billion doubles = 32 GB

if ( ((uint64_t) layers) * (((uint64_t) nc)* ((uint64_t) (Nk+PADB)) + ((uint64_t) nr)* ((uint64_t) (Nk+PADC)) ) > mem_capacity)
	layers = mem_capacity / (((uint64_t) nc)* ((uint64_t) (Nk+PADB)) + ((uint64_t) nr)* ((uint64_t) (Nk+PADC)) );
if (layers < 1)
	layers = 1;

//cout << "layers "<<layers<<endl;

uint64_t sizeB = ((uint64_t) layers )*((uint64_t) nc*(Nk+PADB));
uint64_t sizeC = ((uint64_t) layers )*((uint64_t) nr*(Nk+PADC));

if (B == NULL)
{
	//printf("B init\n");
	srand(unsigned(time(0)));
	B = (TYPE *)_mm_malloc(((uint64_t) sizeof(TYPE))* sizeB,64);
	//RandGen G;
	#pragma omp parallel for
	for(i=0;i<sizeB;i++) 
	{
		//:qB[i] = (TYPE) (rand()%10)/10;
		//B[i] = G.RandReal();
		//B[i] = rand()/((TYPE) RAND_MAX );
		B[i] = ((TYPE) (i%100))/100;
		//B[i] = 1.3;
	}  
	
}

int* ordervec = (int*) malloc(sizeof(int)*nseg);	


for (int i=0; i< nseg ; i++)
	ordervec[i] = i;


D = NULL;


if(order)
{
	//printf("ORDER\n");
	gorder(pb,dc,CORE,nseg,nc,ordervec);
}


if (C == NULL)
{
	C = (TYPE *)_mm_malloc(((uint64_t) sizeof(TYPE))* sizeC,64);
//	memset(C, 0, sizeof(TYPE)*(Nk+PADC)*nr);

	#pragma omp parallel for
	for(i=0;i<sizeC;i++) 
	{
		C[i] = 0.0f;
	}

}




//vout_gold = (TYPE *)malloc(sizeof(TYPE)*ne);



double total_time=0.0, avg_time=0.0;

int K=Nk;

int *segment=pb;
int *column_number=dc;
int *column_index=didx;
int *row_number=d;
TYPE *A_value=dv;
//===========================Dense MM============================================================================================
//===========================Dense MM============================================================================================
//===========================Dense MM============================================================================================clea

__assume_aligned(row_number, 64);
__assume_aligned(A_value, 64);
__assume_aligned(C, 64);
__assume_aligned(B, 64);
__assume_aligned(column_index, 64);        
