#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include "mkl.h"
#include <time.h>
#include <ctime>
#include <cmath>
#include <climits>  
#include <omp.h>
#include <unistd.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include "rdtsc.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <numeric>
#include <random>
#ifdef PAPI
	#include <papi.h>
#endif


#define MFACTOR (32)
#define LOG_MFACTOR (5)
#define MFACTOR16 (16)
#define LOG_MFACTOR16 (4)
#define MFACTOR8 (8)
#define LOG_MFACTOR8 (3)
#define CORE 28
#define TYPE double
#define CACHE 32*1024//35*1024*128 // RI2 specs
#define BC (64)
#define BR (64)
#define PADB 8
#define PADC 8
#define BF 4
#define REP 100
#define BSIZE (BF*32)
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define CEIL(a,b) (((a)+(b)-1)/(b))

using namespace std;


#ifdef PAPI
	#define PAPI_ERROR_CHECK(X) if((X)!=PAPI_OK) std::cerr<<X<<" Error \n";
#endif


struct v_struct {
    int row, col;
    int grp;
    short idx;
    TYPE val;
};

extern int pcnt[100000];
extern int nr0, nr, nc, ne, np, nn, upper_size;
extern int gold_ne;
extern struct v_struct *temp_v;
extern struct v_struct *gold_temp_v;
extern struct v_struct *new_temp_v;
extern int *pb;
//int *dc, *dc;
extern int dcnt;
extern double **X, **Y;
extern int nseg;
extern int *dc, *didx;
extern int *d;
extern TYPE *dv;

extern int sc;


extern int SM_K;
extern int SM_WIDTH;
extern bool HIGH_SD;

#endif
