#define rem_mm(X) if(X!= NULL) {_mm_free(X); X=NULL;}
#define rem(X) if(X!= NULL) {free(X); X=NULL;}

rem_mm(B)
rem_mm(C)
rem_mm(D)
rem(vout_gold)
rem(ordervec)

