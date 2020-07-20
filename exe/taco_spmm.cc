#include "../inc/taco.h"
//#include "../inc/j_stream.h"


int main(int argc, char **argv)
{
	ready(argc, argv);
	struct timeval starttime, midtime, endtime,timediff;
	double total_time=0.0, avg_time=0.0;
	qsort(gold_temp_v, ne, sizeof(struct v_struct), compare1);
	taco_tensor_t *A,*B,*C,*D;

	int Nk = 128;
	A = read_sparse();

	//int Nk = 0;
	if(argc==3)
		Nk = atoi(argv[2]);
	else 
	{
		cout<<"Usage is ./taco_spmm.exe file_name K"<<endl;
		exit(1);
	}
	//for (int base = 32 ; base <= 50 ; base += 18)
	{
		//for (int Nk = base ; Nk <= base * 32 ; Nk*= 2)
		{
			B = create_dense(nr,Nk,1); // Randomly initialized matrix
			C = create_dense(nr,Nk,0); // Zero initialized matrix
			//D = create_sparse(A,0);
			taco_spmm(C,A,B);
			REM_TENSOR(B);
			REM_TENSOR(C);
		}
	}

	

	
	//taco_sddmm(D,A,B,C);

}
