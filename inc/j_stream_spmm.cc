
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
#define vhreduce(a)\
    a=_mm256_hadd_pd(a,a); \
    a=_mm256_hadd_pd(a,a);
#define vlower(a) _mm256_cvtsd_f64(a)


//printf("test\n");
#ifdef KIJ
//printf("KIJ\n");
for (int k = 0; k < K; k += tile_size)   //slices
{
#ifndef SINGLE
#pragma ivdep//k
#pragma vector aligned//k
    //#pragma temporal (C)
    #pragma omp parallel for schedule(static, 1)
#endif

    for(int row_panel = 0; row_panel < nseg; row_panel++)
    {
        if(segment[row_panel] == segment[row_panel + 1])
        {
            continue;
        }

#else
#ifndef SINGLE
#pragma ivdep//k
#pragma vector aligned//k
//#pragma temporal (C)
#pragma omp parallel for schedule(static, 1)
#endif

for(int row_panel = 0; row_panel < nseg; row_panel++)
{
    if(segment[row_panel] == segment[row_panel + 1])
    {
        continue;
    }

    for (int k = 0; k < K; k += tile_size)   //slices
    {
#endif

        for(int j = segment[row_panel]; j < segment[row_panel + 1] ; j++)   //seperate dCSC matrices
        {

            int colNumber = column_number[j];
            int numOfRows = column_index[j + 1] - column_index[j];
            int bindex = colNumber * (K + PADB) ;

            /// 32 "B" elements will be kept in 8 vec registers
            for (int kk = k; kk < MIN(Nk, k + tile_size); kk += 32)
            {
                rtype Breg1 = vload(&B[ bindex + kk + 0 ]);
                rtype Breg2 = vload(&B[ bindex + kk + 4 ]);
                rtype Breg3 = vload(&B[ bindex + kk + 8 ]);
                rtype Breg4 = vload(&B[ bindex + kk + 12]);
                rtype Breg5 = vload(&B[ bindex + kk + 16]);
                rtype Breg6 = vload(&B[ bindex + kk + 20]);
                rtype Breg7 = vload(&B[ bindex + kk + 24]);
                rtype Breg8 = vload(&B[ bindex + kk + 28]);

		//_mm_prefetch(&B[ bindex + kk + 32 + 0 ], _MM_HINT_T0);
		//_mm_prefetch(&B[ bindex + kk + 32 + 4 ], _MM_HINT_T0);
		//_mm_prefetch(&B[ bindex + kk + 32 + 8 ], _MM_HINT_T0);
		//_mm_prefetch(&B[ bindex + kk + 32 + 12], _MM_HINT_T0);
		//_mm_prefetch(&B[ bindex + kk + 32 + 16], _MM_HINT_T0);
		//_mm_prefetch(&B[ bindex + kk + 32 + 20], _MM_HINT_T0);
		//_mm_prefetch(&B[ bindex + kk + 32 + 24], _MM_HINT_T0);
		//_mm_prefetch(&B[ bindex + kk + 32 + 28], _MM_HINT_T0);

                for( int i = 0 ; i < numOfRows ; i ++)
                {
                    int cindex = row_number[i + 0 + column_index[j]] * (K + PADC);
                    TYPE aval = A_value[i + 0 + column_index[j] ] ;
                    rtype Areg = vset(aval);

                    /// The following code is fully unrolled
                    /// we are only going to use 4 vec regs for C to prevent spills
                    rtype Creg1 = vload(&C[ cindex + kk + 0 ]);
                    Creg1 = vfma(Areg, Breg1, Creg1);
                    vstore(&C[ cindex + kk + 0 ], Creg1);

                    rtype Creg2 = vload(&C[ cindex + kk + 4 ]);
                    Creg2 = vfma(Areg, Breg2, Creg2);
                    vstore(&C[ cindex + kk + 4  ], Creg2);

                    rtype Creg3 = vload(&C[ cindex + kk + 8 ]);
                    Creg3 = vfma(Areg, Breg3, Creg3);
                    vstore(&C[ cindex + kk + 8 ], Creg3);

                    rtype Creg4 = vload(&C[ cindex + kk + 12 ]);
                    Creg4 = vfma(Areg, Breg4, Creg4);
                    vstore(&C[ cindex + kk + 12  ], Creg4);

                    /// C reuse beg

                    Creg1 = vload(&C[ cindex + kk + 16 ]);
                    Creg1 = vfma(Areg, Breg5, Creg1);
                    vstore(&C[ cindex + kk + 16 ], Creg1);

                    Creg2 = vload(&C[ cindex + kk + 20 ]);
                    Creg2 = vfma(Areg, Breg6, Creg2);
                    vstore(&C[ cindex + kk + 20  ], Creg2);

                    Creg3 = vload(&C[ cindex + kk + 24 ]);
                    Creg3 = vfma(Areg, Breg7, Creg3);
                    vstore(&C[ cindex + kk + 24 ], Creg3);

                    Creg4 = vload(&C[ cindex + kk + 28 ]);
                    Creg4 = vfma(Areg, Breg8, Creg4);
                    vstore(&C[ cindex + kk + 28  ], Creg4);

                    /// C reuse beg


                }
            }
        }
    }
}

