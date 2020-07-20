#define vsetzero _mm256_setzero_pd 
#define vload _mm256_loadu_pd
#define vntload(addr) \
    _mm256_stream_load_si256((__m256i*)addr)
#define vstore _mm256_storeu_pd
#define vntstore _mm256_stream_pd
#define rtype __m256d
#define vfma _mm256_fmadd_pd
#define vbroadcast(dest, addr) dest = _mm256_broadcast_sd(addr)
#define vset(val)  _mm256_set1_pd(val)
#define vmul(a,b) _mm256_mask_mul_pd(a,b)
#define vadd(a,b) _mm256_add_pd(a,b)
#define vhreduce(a)\
    a=_mm256_hadd_pd(a,a); \
    a=_mm256_hadd_pd(a,a);
#define vlower(a) _mm256_cvtsd_f64(a)

#ifdef KIJ
	for (int k = 0; k < K; k+=tile_size) //slices
	{
		#ifndef SINGLE
			#pragma ivdep//k
			#pragma vector aligned//k
			//#pragma temporal (C)
			#pragma omp parallel for  //schedule(dynamic, 1) 
		#endif
		for(int row_panel=0; row_panel<nseg; row_panel++)
		{
			if(segment[row_panel]==segment[row_panel+1]) continue;
#else
	#ifndef SINGLE
		#pragma ivdep//k
		#pragma vector aligned//k
		//#pragma temporal (C)
		#pragma omp parallel for  //schedule(dynamic, 1) 
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

				int bindex = column_number[j] * (K+PADB);
				for( ; i < numOfRows ; i ++) 
				{ 
					rtype sum1 = vsetzero();
					rtype sum2 = vsetzero();
					rtype sum3 = vsetzero();
					rtype sum4 = vsetzero();
					
					int cindex = row_number[i + column_index[j]] * (K+PADC);	
					int dindex = i + 0 + column_index[j];
					
					for (int kk = k; kk < MIN(k + tile_size ,K); kk+=16)
					{
						rtype bval1 = vload(&B[ bindex + kk + 0 ]);
						rtype bval2 = vload(&B[ bindex + kk + 4 ]);
						rtype bval3 = vload(&B[ bindex + kk + 8 ]);
						rtype bval4 = vload(&B[ bindex + kk + 12 ]);

						rtype cval1 = vload(&C[ cindex + kk + 0]);
						rtype cval2 = vload(&C[ cindex + kk + 4]);
						rtype cval3 = vload(&C[ cindex + kk + 8]);
						rtype cval4 = vload(&C[ cindex + kk + 12]);

						sum1 = vfma(bval1, cval1, sum1);
						sum2 = vfma(bval2, cval2, sum2);
						sum3 = vfma(bval3, cval3, sum3);
						sum4 = vfma(bval4, cval4, sum4);
					}

					sum1 = vadd(sum1, sum2);
					sum3 = vadd(sum3, sum4);
					sum1 = vadd(sum1, sum3);

					double tmp = hsum_double_avx(sum1);

					D[dindex] +=  tmp;


					//for (int kk = k; kk < k + tile_size; kk++)
					//{
					//D[ i + 0 + column_index[j] ] += B[ column_number[j] * (K+PADB) + kk + 0 ] * C[ row_number[i + column_index[j]] * (K+PADC) + kk + 0];
					//}
				}

			}
		}
	}

#pragma omp parallel for 
for (int i=0; i<ne ; i++)
	D[i] *= A_value[i];
