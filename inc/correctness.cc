{
	#pragma omp parallel for
	for(i=0;i<nr*(Nk+PADC);i++) 
	{
		C[i] = 0.0f;
	}

	#include "../inc/j_stream_spmm.cc"

	for(i=0;i<nr*sc;i++) {
		vout_gold[i] = 0.0f;
	}
	for(i=0;i<gold_ne;i++) {
		for(j=0;j<sc;j++) {
			vout_gold[gold_temp_v[i].row*(K+PADC)+j] += B[(K+PADB)*gold_temp_v[i].col+j] * gold_temp_v[i].val;
			// cout<<" C " << vout_gold[gold_temp_v[i].row*sc+j] <<"  "<<B[sc*gold_temp_v[i].col+j] << "   "<<gold_temp_v[i].val<<endl;
		}
	}
	// for(int y=0;y<nr*sc;y++)
	// {
	//     if(vout_gold[y]!=0)
	//         cout<<" non zeros in vout_dold " << vout_gold[y]<<endl;
	// }
	int num_diff=0;
	for(i=0;i<nr*sc;i++) {
		//if(i%nr != i%nr0) continue;
		TYPE p1 = vout_gold[i]; TYPE p2 = C[i];
		if(p1 < 0) p1 *= -1;
		if(p2 < 0) p2 *= -1;
		TYPE diff;
		diff = p1 - p2;
		if(diff < 0) diff *= -1;
		if(diff / MAX(p1,p2) > 0.01) {
			if(num_diff < 20*1) fprintf(stdout, "row %d col %d %lf %lf\n",i/sc, i%sc, C[i], vout_gold[i]);
			num_diff++;
		}
	}
	//    fprintf(stdout, "num_diff : %d\n", num_diff);
	if(num_diff != 0)
	fprintf(stdout, " DIFF %lf, %d\n", (double)num_diff/(nr*sc)*100, num_diff);
	//    fprintf(stdout, "ne : %d\n", gold_ne);
	//#endif
	//fprintf(stdout, "\n");
}