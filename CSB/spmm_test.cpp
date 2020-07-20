#define NOMINMAX
#include <iostream> 
#include <algorithm>
#include <numeric>
#include <functional>
#include <fstream>
#include <ctime>
#include <cmath>
#include <string>
#include <array>

#include "timer.gettimeofday.c"
#include "cilk_util.h"
#include "aligned.h"

#ifdef PAPI
#include <papi.h>
#define PAPI_ERROR_CHECK(X) if((X)!=PAPI_OK) std::cerr<<X<<" Error \n";
#endif


#define INDEXTYPE uint32_t
#ifdef SINGLEPRECISION
	#define VALUETYPE float
#else
	#define VALUETYPE double
#endif

#ifndef RHSDIM
	#define RHSDIM 16
#endif
#define ALIGN 32

#include "utility.h"
#include "triple.h"
#include "csc.h"
#include "bicsb.h"
#include "bmcsb.h"
#include "spvec.h"
#include "Semirings.h"

#define SWAPPY std::swap(x,y_bicsb);
#define SWAPPYT std::swap(xt,yt_bicsb);


using namespace std;


template <typename NT, typename ALLOC, int DIM>
void fillzero (vector< array<NT,DIM>, ALLOC > & vecofarr)
{
	for(auto itr = vecofarr.begin(); itr != vecofarr.end(); ++itr)
	{
		itr->fill(static_cast<NT> (0));
	}
}



void cacheFlush(double *X, double *Y) {

}



template <typename NT, typename ALLOC, int DIM>
void fillrandom (vector< array<NT,DIM>, ALLOC > & vecofarr)
{
	for(auto itr = vecofarr.begin(); itr != vecofarr.end(); ++itr)
	{
	#if (__GNUC__ == 4 && (__GNUC_MINOR__ < 7) )
		RandGen G;
		for(auto refarr = itr->begin(); refarr != itr->end(); ++refarr)
		{
			*refarr = G.RandReal();
		}
	#else
		std::uniform_real_distribution<NT> distribution(0.0f, 1.0f); //Values between 0 and 1
		std::mt19937 engine; // Mersenne twister MT19937
		auto generator = std::bind(distribution, engine);
		std::generate_n(itr->begin(), DIM, generator); 
	#endif	
	}
}

template <typename NT, typename ALLOC, int DIM>
void VerifyMM (vector< array<NT,DIM>, ALLOC > & control, vector< array<NT,DIM>, ALLOC > & test)
{
    	vector< array<NT,DIM> > error;
	pair<size_t, size_t> maxerrloc;
	NT prevmax = 0.0;
	
	for(auto itr1 = control.begin(), itr2 = test.begin(); itr1 != control.end(); ++itr1, ++itr2)
	{
		array<NT,DIM> entry;
		transform(itr1->begin(), itr1->end(), itr2->begin(), entry.begin(), absdiff<NT>());
		auto maxelement = max_element(entry.begin(), entry.end());
		size_t maxcol = maxelement-entry.begin();
		if(*maxelement > prevmax)
		{
			maxerrloc = make_pair(itr1- control.begin(), maxcol);
			prevmax = *maxelement;
		}
		error.emplace_back(entry);
	}
	
    	cout << "Max error is: " << prevmax << " on y[" << maxerrloc.first <<"][" << maxerrloc.second << "]=";
	cout << test[maxerrloc.first][maxerrloc.second] << endl;
    
	NT machEps = machineEpsilon<NT>();
    	cout << "Absolute machine epsilon is: " << machEps <<" and y[" << maxerrloc.first <<"][" << maxerrloc.second;
	cout << "]*EPSILON becomes " << machEps * test[maxerrloc.first][maxerrloc.second] << endl;
    
    	NT sqrtm = sqrt(static_cast<NT>(control.size()));
    	cout << "sqrt(n) * relative error is: " << abs(machEps * test[maxerrloc.first][maxerrloc.second]) * sqrtm << endl;
    	if ( (abs(machEps * test[maxerrloc.first][maxerrloc.second]) * sqrtm) < abs(prevmax))
    	{
        	cout << "*** ATTENTION ***: error is more than sqrt(n) times the relative machine epsilon" << endl;
    	}
}



int main(int argc, char* argv[])
{
#ifndef CILK_STUB
	int gl_nworkers = __cilkrts_get_nworkers();
#else
	int gl_nworkers = 0;
#endif
	bool syminput = false;
	bool binary = false;
	bool iscsc = false;
	INDEXTYPE m = 0, n = 0, nnz = 0, forcelogbeta = 0;
	string inputname;
	if(argc < 2)
	{
		cout << "Normal usage: ./a.out inputmatrix.mtx sym/nosym binary/text triples/csc" << endl;
		cout << "Assuming matrix.txt is the input, matrix is unsymmetric, and stored in text(ascii) file" << endl;
		inputname = "matrix.txt";
	}
	else if(argc < 3)
	{
		cout << "Normal usage: ./a.out inputmatrix.mtx sym/nosym binary/text triples/csc" << endl;
		cout << "Assuming that the matrix is unsymmetric, and stored in text(ascii) file" << endl;
		inputname =  argv[1];
	}
	else if(argc < 4)
	{
		cout << "Normal usage: ./a.out inputmatrix.mtx sym/nosym binary/text triples/csc" << endl;
		cout << "Assuming matrix is stored in text(ascii) file" << endl;
		inputname =  argv[1];
		string issym(argv[2]);
		if(issym == "sym")
			syminput = true;
		else if(issym == "nosym")
			syminput = false;
		else
			cout << "unrecognized option, assuming nosym" << endl;
	}
	else
	{
		inputname =  argv[1];
		string issym(argv[2]);
		if(issym == "sym")
			syminput = true;
		else if(issym == "nosym")
			syminput = false;
		else
			cout << "unrecognized option, assuming unsymmetric" << endl;

		string isbinary(argv[3]);
		if(isbinary == "text")
			binary = false;
		else if(isbinary == "binary")
			binary = true;
		else
			cout << "unrecognized option, assuming text file" << endl;
	
		if(argc > 4)
		{
			string type(argv[4]);
			if(type == "csc")
			{
				iscsc = true;
				cout << "Processing CSC binary" << endl;
			}
		}
			
		if(argc == 6)
			forcelogbeta = atoi(argv[5]);
	}

	Csc<VALUETYPE, INDEXTYPE> * csc;
	if(binary)
	{
		FILE * f = fopen(inputname.c_str(), "r");
		if(!f)
		{
			cerr << "Problem reading binary input file\n";
			return 1;
		}
		if(iscsc)
		{
			fread(&n, sizeof(INDEXTYPE), 1, f);
			fread(&m, sizeof(INDEXTYPE), 1, f);
			fread(&nnz, sizeof(INDEXTYPE), 1, f);
		}
		else
		{
			fread(&m, sizeof(INDEXTYPE), 1, f);
			fread(&n, sizeof(INDEXTYPE), 1, f);
			fread(&nnz, sizeof(INDEXTYPE), 1, f);
		}
		if (m <= 0 || n <= 0 || nnz <= 0)
		{
			cerr << "Problem with matrix size in binary input file\n";	
			return 1;		
		}
		long tstart = cilk_get_time();	// start timer
		cout << "Reading matrix with dimensions: "<< m << "-by-" << n <<" having "<< nnz << " nonzeros" << endl;
		INDEXTYPE * rowindices = new INDEXTYPE[nnz];
		VALUETYPE * vals = new VALUETYPE[nnz];
		INDEXTYPE * colindices;
		INDEXTYPE * colpointers;
		if(iscsc)
		{
			colpointers = new INDEXTYPE[n+1];
			size_t cols = fread(colpointers, sizeof(INDEXTYPE), n+1, f);
			if(cols != n+1)
			{
				cerr << "Problem with FREAD, aborting... " << endl;
                        	return -1;
			}
		}
		else
		{
			colindices = new INDEXTYPE[nnz];
			size_t cols = fread(colindices, sizeof(INDEXTYPE), nnz, f);
			if(cols != nnz)
			{
				cerr << "Problem with FREAD, aborting... " << endl;
                        	return -1;
			}
		}
		size_t rows = fread(rowindices, sizeof(INDEXTYPE), nnz, f);
		size_t nums = fread(vals, sizeof(VALUETYPE), nnz, f);

		if(rows != nnz || nums != nnz)
		{
			cerr << "Problem with FREAD, aborting... " << endl;
			return -1;
		}
		long tend = cilk_get_time();	// end timer	
		cout<< "Reading matrix in binary took " << ((VALUETYPE) (tend-tstart)) /1000 << " seconds" <<endl;
		fclose(f);
		if(iscsc)
		{
			csc = new Csc<VALUETYPE, INDEXTYPE>();
			csc->SetPointers(colpointers, rowindices, vals , nnz, m, n, true);	// do the fortran thing
			// csc itself will manage the data in this case (shallow copy)
		}
		else
		{
			csc = new Csc<VALUETYPE, INDEXTYPE>(rowindices, colindices, vals , nnz, m, n);
			delete [] colindices;
			delete [] rowindices;
			delete [] vals;
		}
	}
	else
	{
		cout << "reading input matrix in text(ascii)... " << endl;
		ifstream infile(inputname.c_str());
		char line[256];
		char c = infile.get();
		while(c == '%')
		{
			infile.getline(line,256);
			c = infile.get();
		}
		infile.unget();
		infile >> m >> n >> nnz;	// #{rows}-#{cols}-#{nonzeros}

		long tstart = cilk_get_time();	// start timer	
		Triple<VALUETYPE, INDEXTYPE> * triples = new Triple<VALUETYPE, INDEXTYPE>[nnz];
	
		if (infile.is_open())
		{
			INDEXTYPE cnz = 0;	// current number of nonzeros
			while (! infile.eof() && cnz < nnz)
			{
				infile >> triples[cnz].row >> triples[cnz].col >> triples[cnz].val;	// row-col-value
				triples[cnz].row--;
				triples[cnz].col--;
				++cnz;
			}
			assert(cnz == nnz);	
		}
		long tend = cilk_get_time();	// end timer	
		cout<< "Reading matrix in ascii took " << ((double) (tend-tstart)) /1000 << " seconds" <<endl;
	
		cout << "converting to csc ... " << endl;
		csc= new Csc<VALUETYPE,INDEXTYPE>(triples, nnz, m, n);
		delete [] triples;
	}

	cout << "# workers: "<< gl_nworkers << endl;
	BiCsb<VALUETYPE, INDEXTYPE> bicsb(*csc, gl_nworkers, forcelogbeta);
		
	double mflops = (2.0 * static_cast<double>(nnz) * RHSDIM) / 1000000.0;
	cout << "generating " << RHSDIM << " multi vectors... " << endl;
	typedef array<VALUETYPE, RHSDIM> PACKED;


	uint64_t layers = 16;

	/*
	//vector< PACKED, aligned_allocator<PACKED, ALIGN> >  z[layers]; //(n);
	//vector< PACKED, aligned_allocator<PACKED, ALIGN> > tempz(n);
	//fillrandom<VALUETYPE, aligned_allocator<PACKED, ALIGN>, RHSDIM> (x);
	for(int i =0 ; i< layers; i++)
	{
		z[i] = tempz;

	}
	//PACKED ptest = 
	//z[0][0][0] = 10;
	//z[1][0][0] = 20;

	//cout << z[0][0][0] << " " << z[1][0][0] << endl;
	*/
	//	vector< PACKED, aligned_allocator<PACKED, ALIGN> > z[4];
	vector< PACKED, aligned_allocator<PACKED, ALIGN> > x(n);
	vector< PACKED, aligned_allocator<PACKED, ALIGN> > y_bicsb(m);
	vector< PACKED, aligned_allocator<PACKED, ALIGN> > y_csc(m);
	vector< PACKED, aligned_allocator<PACKED, ALIGN> > xs[layers];
	vector< PACKED, aligned_allocator<PACKED, ALIGN> > y_bicsbs[layers];

	fillzero<VALUETYPE, aligned_allocator<PACKED, ALIGN>, RHSDIM>(y_csc);
	fillzero<VALUETYPE, aligned_allocator<PACKED, ALIGN>, RHSDIM>(y_bicsb);
	fillrandom<VALUETYPE, aligned_allocator<PACKED, ALIGN>, RHSDIM>(x);
	
	for (int i=0; i<layers ; i++)
	{
		xs[i] = x;
		y_bicsbs[i] = y_bicsb;
	}



	typedef PTSRArray<VALUETYPE,VALUETYPE, RHSDIM> PTARR;		
	cout << "starting SpMV ... " << endl;
	cout << "Row imbalance is: " << RowImbalance(bicsb) << endl;
	cout << "Col imbalance is: " << ColImbalance(bicsb) << endl;
	timer_init();

	double *X = NULL ; // (double*)malloc(20 * 1000000 * sizeof(double));
	double *Y =  NULL ; //(double*)malloc(20 * 1000000 * sizeof(double));

	
	bicsb_gespmv<PTARR>(bicsb, &(x[0]), &(y_bicsb[0]));
	double t0,t1,sum = 0;// = timer_seconds_since_init();

	for(int i=0; i < REPEAT; ++i)
	{
		cacheFlush(X, Y);
		//std::swap()
		x = xs[i%layers];
		y_bicsb = y_bicsbs[i%layers];

		t0 = timer_seconds_since_init();

		bicsb_gespmv<PTARR>(bicsb, &(x[0]), &(y_bicsb[0]));

		t1 = timer_seconds_since_init();
		double single_time = t1-t0;
		
		//SWAPPY
		cout<< "CSB" << " time: " << single_time  << " Gflop/sec: " << mflops  / single_time / 1000 <<endl;

		sum += t1 - t0;
	}

	#ifdef PAPI
		{
			#include "papi_init_reg.cc"
			bicsb_gespmv<PTARR>(bicsb, &(x[0]), &(y_bicsb[0]));		
			#include "papi_end.cc"
		}
		{
			#include "papi_init_lx.cc"
			bicsb_gespmv<PTARR>(bicsb, &(x[0]), &(y_bicsb[0]));		
			#include "papi_end.cc"
		}
		{
			#include "papi_init_stall.cc"
			bicsb_gespmv<PTARR>(bicsb, &(x[0]), &(y_bicsb[0]));		
			#include "papi_end.cc"
		}
	#endif
//	double t1 = timer_seconds_since_init();

	double time = (sum)/REPEAT;
	cout<< "Flushed BiCSB" << " time: " << time << " seconds" <<endl;
	cout<< "Flushed BiCSB" << " mflop/sec: " << mflops  / time <<endl;


	cout << "starting SpMV_T" << endl;
	vector< PACKED, aligned_allocator<PACKED, ALIGN> > xt(m);
	vector< PACKED, aligned_allocator<PACKED, ALIGN> > yt_bicsb(n);
	vector< PACKED, aligned_allocator<PACKED, ALIGN> > yt_csc(n);
	
	fillzero<VALUETYPE, aligned_allocator<PACKED, ALIGN>, RHSDIM>(yt_csc);
	fillzero<VALUETYPE, aligned_allocator<PACKED, ALIGN>, RHSDIM>(yt_bicsb);
	fillrandom<VALUETYPE, aligned_allocator<PACKED, ALIGN>, RHSDIM>(xt);
	
	bicsb_gespmvt<PTARR>(bicsb, &(xt[0]), &(yt_bicsb[0]));		// warm-up computation	
	sum = 0;
	//t0 = timer_seconds_since_init();
	for(int i=0; i < REPEAT; ++i)
	{
		cacheFlush(X, Y);
		t0 = timer_seconds_since_init();
		
		bicsb_gespmvt<PTARR>(bicsb, &(xt[0]), &(yt_bicsb[0]));

		t1 = timer_seconds_since_init();
		//SWAPPYT
		sum += t1 - t0;
	//	bicsb_gespmvt<PTARR>(bicsb, &(xt[0]), &(yt_bicsb[0]));
	}
	#ifdef PAPI
		{
			#include "papi_init_reg.cc"
			bicsb_gespmv<PTARR>(bicsb, &(x[0]), &(y_bicsb[0]));		
			#include "papi_end.cc"
		}
		{
			#include "papi_init_lx.cc"
			bicsb_gespmv<PTARR>(bicsb, &(x[0]), &(y_bicsb[0]));		
			#include "papi_end.cc"
		}
		{
			#include "papi_init_stall.cc"
			bicsb_gespmv<PTARR>(bicsb, &(x[0]), &(y_bicsb[0]));		
			#include "papi_end.cc"
		}
	#endif
	//t1 = timer_seconds_since_init();
	
	double totaltime = time + (sum)/REPEAT;
	time = (t1-t0)/REPEAT;
	cout<< "Flushed BiCSB Trans" << " time: " << time << " seconds" <<endl;
	cout<< "Flushed BiCSB Trans" << " mflop/sec: " << mflops  / time <<endl;
	

	cout<< "Flushed BiCSB Total" << " time: " << totaltime << " seconds" <<endl;
	cout<< "Flushed BiCSB Total" << " mflop/sec: " << 2*mflops  / totaltime <<endl;
	
	// Verify with CSC (serial)
	csc_gaxpy_mm<RHSDIM>(*csc, &(x[0]), &(y_csc[0]));
	t0 = timer_seconds_since_init();
	for(int i=0; i < REPEAT; ++i)
	{
	    	csc_gaxpy_mm<RHSDIM>(*csc, &(x[0]), &(y_csc[0]));
	    	//SWAPPY
	}
	t1 = timer_seconds_since_init();
	double csctime = (t1-t0)/REPEAT;
	cout<< "CSC" << " time: " << csctime << " seconds" <<endl;
	cout<< "CSC" << " mflop/sec: " << mflops / csctime <<endl;

	VerifyMM<VALUETYPE, aligned_allocator<PACKED, ALIGN>, RHSDIM>(y_csc, y_bicsb);
	
	
	csc_gaxpy_mm_trans<RHSDIM> ( *csc, &(xt[0]), &(yt_csc[0]));
	t0 = timer_seconds_since_init();
	for(int i=0; i < REPEAT; ++i)
	{
		csc_gaxpy_mm_trans<RHSDIM> ( *csc, &(xt[0]), &(yt_csc[0]));
		//SWAPPYT
	}
	t1 = timer_seconds_since_init();
	time = (t1-t0)/REPEAT;
	cout <<"Transposed CSC time: " << time << " seconds" << endl;
	cout <<"Transposed CSC mflop/sec: " << mflops/ time << endl;
	
	VerifyMM<VALUETYPE, aligned_allocator<PACKED, ALIGN>, RHSDIM>(yt_csc, yt_bicsb);	
	
	delete csc;
	
}

