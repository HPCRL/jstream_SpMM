#pragma omp parallel for
for(i=0;i<ne;i++) 
{
	D[i] = 0.0f;
}

#include "../inc/j_stream_sddmm.cc"

for(i = 0; i < ne; i++)
{
	vout_gold[i] = 0.0f;
}

for(i = 0; i < ne; i++)
{
	for(j = 0; j < Nk; j++)
	{
		vout_gold[i] += C[gold_temp_v[i].row * (sc+PADC) + j] *B[(sc+PADB) * gold_temp_v[i].col + j] * gold_temp_v[i].val;
		// cout<<" C " << vout_gold[gold_temp_v[i].row*sc+j] <<"  "<<B[sc*gold_temp_v[i].col+j] << "   "<<gold_temp_v[i].val<<endl;
	}
}

// for(int y=0;y<nr*sc;y++)
// {
//     if(vout_gold[y]!=0)
//         cout<<" non zeros in vout_dold " << vout_gold[y]<<endl;
// }
int num_diff = 0;

for(i = 0; i < ne; i++)
{
	//if(i%nr != i%nr0) continue;
	TYPE p1 = vout_gold[i];
	TYPE p2 = D[gold_temp_v[i].grp];

	if(p1 < 0)
	{
		p1 *= -1;
	}

	if(p2 < 0)
	{
		p2 *= -1;
	}

	TYPE diff;
	diff = p1 - p2;

	if(diff < 0)
	{
		diff *= -1;
	}

	if(diff / MAX(p1, p2) > 0.01)
	{
						if(num_diff < 20*1) fprintf(stdout, "row %d col %d %lf %lf\n",i/sc, i%sc, D[i], vout_gold[i]);
		num_diff++;
	}

	//      cout<<p1<<" ";
	//      if((i+1)%sc == 0)
	//          cout<<endl;
}

//    fprintf(stdout, "num_diff : %d\n", num_diff);
if(num_diff != 0)
{
	fprintf(stdout, " DIFF %lf, %d\n", (double)num_diff / (ne) * 100, num_diff);
}