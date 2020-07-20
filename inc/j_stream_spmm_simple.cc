#ifdef KIJ
	//printf("KIJ\n");
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


				for(; i < ((int)numOfRows/8) * 8 ; i+=8) 
				{
					#pragma ivdep
					#pragma unroll(8)
					#pragma vector nontemporal (A_value)
					#pragma prefetch C:_MM_HINT_T1
					//#pragma prefetch B:_MM_HINT_T1
					//#pragma temporal (C)
					for (int kk = k; kk < MIN(Nk,k+tile_size); kk++)// i, rows of 
					{

						C[ row_number[i + 0 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 0 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 0 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 1 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 1 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 1 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 2 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 2 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 2 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 3 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 3 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 3 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 4 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 4 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 4 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 5 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 5 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 5 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 6 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 6 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 6 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 7 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 7 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 7 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						//C[ rowInd * K + kk + kkk ] = C[ rowInd * K + kk + kkk ] + (tx * B[ col * K + kk + kkk ] );
					}
				}
				if( numOfRows - i >= 4)
				{

					#pragma ivdep
					#pragma vector nontemporal (A_value)
					#pragma prefetch C:_MM_HINT_T1
					#pragma unroll(8)
					//#pragma temporal (C)
					for (int kk = k; kk < MIN(Nk,k+tile_size); kk++)// i, rows of 
					{

						C[ row_number[i + 0 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 0 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 0 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 1 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 1 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 1 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 2 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 2 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 2 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						C[ row_number[i + 3 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 3 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 3 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
						

					}
					i+=4;
				}

				for( ; i < numOfRows ; i ++) 
				{ 
					#pragma ivdep//k
					#pragma vector nontemporal (A_value)//k
					#pragma prefetch C:_MM_HINT_T1
					#pragma prefetch B:_MM_HINT_T2
					#pragma unroll(8)
					//#pragma temporal (C)
					for (int kk = k; kk < MIN(Nk,k+tile_size); kk++)
					{
						C[ row_number[i + 0 + column_index[j]] * (K+PADC) + kk  ] = C[ row_number[i + 0 + column_index[j]] * (K+PADC) + kk  ] + (A_value[i + 0 + column_index[j] ] * B[ colNumber * (K+PADB) + kk  ] );
					}
				}
			}
		}
	}