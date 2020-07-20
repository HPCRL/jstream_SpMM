#include "../inc/taco.h"


int main(int argc, char **argv)
{
	ready(argc, argv);
	struct timeval starttime, midtime, endtime,timediff;
	double total_time=0.0, avg_time=0.0;
	qsort(gold_temp_v, ne, sizeof(struct v_struct), compare1);
	taco_tensor_t *A,*B,*C,*D;

	int Nk = 128;
	A = read_sparse();
	D = create_sparse(A,0);
	/*
	B = create_dense(nr,Nk,1); // Randomly initialized matrix
	C = create_dense(nr,Nk,1); // Randomly initialized matrix
	D = create_sparse(A,0);

	//taco_spmm(C,A,B);
	taco_sddmm(D,A,B,C);


	Nk = 1024;
	B = create_dense(nr,Nk,1); // Randomly initialized matrix
	C = create_dense(nr,Nk,1); // Randomly initialized matrix
	D = create_sparse(A,0);
	taco_sddmm(D,A,B,C);
	*/


	if(argc==3)
		Nk = atoi(argv[2]);
	else 
	{
		cout<<"Usage is ./taco_sddmm.exe file_name K"<<endl;
		exit(1);
	}
	//for (int base = 32 ; base <= 50 ; base += 18)
	{
		//for (int Nk = base ; Nk <= base * 32 ; Nk*= 2)
		{
			B = create_dense(nr,Nk,1); // Randomly initialized matrix
			C = create_dense(nr,Nk,1); // Zero initialized matrix
			//D = create_sparse(A,0);
			taco_sddmm(D,A,B,C);
			REM_TENSOR(B);
			REM_TENSOR(C);
		}
	}

}
