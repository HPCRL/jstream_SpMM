#ifdef KIJ
	for (int k = 0; k < K; k+=tile_size) //slices
	{
		#ifndef SINGLE
			#pragma ivdep//k
			#pragma vector aligned//k
			//#pragma temporal (C)
			#pragma omp parallel for num_threads(CORE) //schedule(dynamic, 1) 
		#endif
		for(int row_panel=0; row_panel<nseg; row_panel++)
		{
			if(segment[row_panel]==segment[row_panel+1]) continue;
#else
	#ifndef SINGLE
		#pragma ivdep//k
		#pragma vector aligned//k
		//#pragma temporal (C)
		#pragma omp parallel for num_threads(CORE) //schedule(dynamic, 1) 
	#endif
	for(int row_panel=0; row_panel<nseg; row_panel++)
	{
		if(segment[row_panel]==segment[row_panel+1]) continue;
		for (int k = 0; k < K; k+=tile_size) //slices
		{
#endif


			for(int j= segment[row_panel]; j< segment[row_panel+1] ; j++) //seperate dCSC matrices 
			{
				int colNumber = column_number[j];
				int numOfRows = column_index[j+1] - column_index[j];
				int i = 0;
				for( ; i < numOfRows ; i ++) 
				{ 
					for (int kk = k; kk < k + tile_size; kk++)
					{
						D[ i + 0 + column_index[j] ] += B[ column_number[j] * (K+PADB) + kk + 0 ] * C[ row_number[i + column_index[j]] * (K+PADC) + kk + 0];
					}
					//if (k + Tk == Nk)       D[ i + 0 + column_index[j] ] *= A_value[ i + 0 + column_index[j] ];
			    }

			}
		}
	}

#pragma omp parallel for num_threads(CORE)
for (int i=0; i<ne ; i++)
	D[i] *= A_value[i];
