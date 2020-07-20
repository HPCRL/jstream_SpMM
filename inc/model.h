#ifndef MODEL_H
#define MODEL_H

#include "../inc/reader.h"
#include <cmath>
#include <math.h> 

int active_col();

extern uint64_t** intervalsp;
extern uint64_t** notsmallerp;
extern uint64_t** notsmallerweightedp;

typedef unsigned __int128 uint128_t;
int init_model(bool print=false);
int init_model_parallel(bool print=false);
int pick_tile(int &Ti, int &Tk , int Nk = 1024, int cache = CACHE);
int pick_tile_sddmm(int &Ti, int &Tk , int Nk = 1024, int cache = CACHE);
uint64_t active_cols(int Ti, bool print);

#endif