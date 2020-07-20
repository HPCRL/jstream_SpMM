#include "../inc/reader.h"
#include "../inc/cacheflush.h"




int main(int argc, char **argv)
{
	ready(argc, argv);
	qsort(gold_temp_v, ne, sizeof(struct v_struct), compare1);

	//struct timeval starttime, midtime, endtime,timediff;
	//double total_time=0.0, avg_time=0.0;
	
	int k = nc;
	int m = nr;
	TYPE         alpha = 1.0, beta = 0.0;
	char          transa, uplo, nonunit;
	char          matdescra[6];
	MKL_INT       i, j, is;

	transa = 'N';
	matdescra[0] = 'G';
	matdescra[1] = 'L';
	matdescra[2] = 'N';
	matdescra[3] = 'C';
	TYPE *Aval = (TYPE *)malloc(sizeof(TYPE)*ne);
	int *rowPtr = (int*)malloc(sizeof(int)*(nr+1));
	int *colInd = (int*)malloc(sizeof(int)*(ne));


	rowPtr[0] = 0;
	int row_id = 0;
	for(int i = 0 ; i< ne ; i++)
	{
		Aval[i] = gold_temp_v[i].val;

		colInd[i] = gold_temp_v[i].col;
		while(gold_temp_v[i].row > row_id)
		{
			row_id ++;
			rowPtr[row_id] = i;
			//cout<<row_id << " " << i<<endl;
		}
		
	}

	while(gold_temp_v[ne].row > row_id)
	{
		row_id ++;
		rowPtr[row_id] = ne;
		//cout<<row_id << " " << i<<endl;
	}
	
		

	struct timeval starttime, midtime, endtime,timediff;
	long double total_time=0.0, avg_time=0.0;

	int n = 0;
	
	if(argc==3)
		n = atoi(argv[2]);
	else 
	{
		cout<<"Usage is ./mkl.exe file_name K"<<endl;
		exit(1);
	}
	

	

	

	//for (int base = 32 ; base <= 50 ; base += 18)
	{
	//	for (int n= base ; n<= base * 32 ; n*= 2)
		{

			//n = 1;
			

			uint64_t singleB = ((uint64_t) k)* ((uint64_t) n);
			uint64_t singleC = ((uint64_t) m)* ((uint64_t) n);
			

			uint64_t layers = 16;
			uint64_t mem_capacity = 4000000000; // 4 billion doubles = 32 GB

			if ( ((uint64_t) layers) * (singleB + singleC) > mem_capacity)
				layers = mem_capacity / (singleB + singleC);
			if (layers < 1)
				layers = 1;

			cout << "layers "<<layers<<endl;
			uint64_t sizeB = ((uint64_t) layers )* singleB;
			uint64_t sizeC = ((uint64_t) layers )* singleC;
			

			TYPE *B_all = (TYPE *)malloc( ((uint64_t) sizeof(TYPE)) * sizeB);
			TYPE *C_all = (TYPE *)malloc( ((uint64_t) sizeof(TYPE)) * sizeC);
			for(int i=0;i<sizeB; i++) 
			{
				B_all[i] = ((TYPE)   ( (i % singleB) %100))/100;
			}
			for(int i=0;i< sizeC; i++) 
			{
				C_all[i] = 0.0;
			}
			TYPE* B = B_all;
			TYPE* C = C_all;

			mkl_dcsrmm (&transa, &m, &n, &k, &alpha, matdescra, Aval, colInd, rowPtr,  &(rowPtr[1]), &(B[0]), &n,  &beta, &(C[0]), &n);

			int rep = REP;

			cout<<"warmup finished"<<endl;		
			double* elapsed_arr = new double[REP];


			for(int i=0; i<rep ; i++)
			{
				cacheFlush(X,Y);
				uint64_t jumpB = sizeB / layers * (i % layers);
				uint64_t jumpC = sizeC / layers * (i % layers);
				B = B_all + jumpB;
				C = C_all + jumpC;

				//gettimeofday(&starttime,NULL);
				auto start = std::chrono::high_resolution_clock::now();
				mkl_dcsrmm (&transa, &m, &n, &k, &alpha, matdescra, Aval, colInd, rowPtr,  &(rowPtr[1]), &(B[0]), &n,  &beta, &(C[0]), &n);
				//gettimeofday(&endtime,NULL);
				auto end = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> diff = end-start;
				double elapsed = diff.count();
				elapsed_arr[i] = elapsed;
				TYPE sum = 0;
				for(int i=0;i<m*n;i++) 
				sum += C[i];
				//cout<<sum<<endl;
				//double elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
				//cout<<"MKL "<<argv[1] << ", K, "<<n <<", time, "<<elapsed<<", GFLOPS:, "<<((double)2*ne*n)/elapsed/1e9<<endl;	
			}

			
			double elapsed = findMedian(elapsed_arr,REP);
			cout<<"MKL median"<<argv[1] << ", K, "<<n <<", time, "<<elapsed<<", GFLOPS:, "<<((double)2*ne*n)/elapsed/1e9<<endl;	

			elapsed = elapsed_arr[0];
			cout<<"MKL max"<<argv[1] << ", K, "<<n <<", time, "<<elapsed<<", GFLOPS:, "<<((double)2*ne*n)/elapsed/1e9<<endl;	


			elapsed = elapsed_arr[REP-1];
			cout<<"MKL min"<<argv[1] << ", K, "<<n <<", time, "<<elapsed<<", GFLOPS:, "<<((double)2*ne*n)/elapsed/1e9<<endl;	
			

						

			for(int i=0;i<m*n;i++) 
			{
				C[i] = 0.0;
			}
			mkl_dcsrmm (&transa, &m, &n, &k, &alpha, matdescra, Aval, colInd, rowPtr,  &(rowPtr[1]), &(B[0]), &n,  &beta, &(C[0]), &n);

			#define ERROR_RATE 10
			int err_cnt[ERROR_RATE] ;
			TYPE threshold = 0.1;

			for(int i=0; i<ERROR_RATE ; i++)
				err_cnt[i] = 0;

			for(int i = 0 ; i < m ; i++)
			{
				for(int kk = 0 ; kk < n ; kk++)
				{
					TYPE val = 0;
					for (int jp = rowPtr[i] ; jp < rowPtr[i+1] ; jp++)
					{
						int j = colInd[jp];
						val += Aval[jp]*B[j*n + kk];
						//val += Aval[jp]*B[kk*n + j];
					}
					TYPE diff = C[i*n + kk] - val;
					//TYPE diff = C[kk*n + i] - val;
					//cout<<val<<" ";
					diff /= val;
					if(diff < 0)
						diff = -diff;
					
					for(int err_no = 0; err_no < ERROR_RATE ; err_no ++)
					{
						TYPE thresh = threshold*(1+err_no);
						if ( diff > thresh )
						{
							//if (err_cnt[err_no] < 10)	cout << diff << " "<< val <<endl;
							err_cnt[err_no] ++;
						}
					}
				}
				//cout<<endl;
			}

			//if (err_cnt[0]>0)
			{
				cout<< "Error % for (";
				for(int err_no = 0; err_no < ERROR_RATE ; err_no ++)
					cout<< threshold*(1+err_no) <<  ((err_no == (ERROR_RATE -1)) ? "":",");
				cout<<") : ";
				for(int err_no = 0; err_no < ERROR_RATE ; err_no ++)
					cout<< (int) (( 100*((double)err_cnt[err_no])) / ((double)n) /((double)k) )<<" ";
			}
			cout << endl;
			
		}
	}
	
						
	
}
