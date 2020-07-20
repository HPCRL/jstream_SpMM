//#include "../inc/common.h"
#ifndef READER_CC
#define READER_CC
#include "../inc/reader.h"


int pcnt[100000];
int nr0, nr, nc, ne, np, nn, upper_size;
int gold_ne;
struct v_struct *temp_v;
struct v_struct *gold_temp_v;
struct v_struct *new_temp_v;
//int *pb;
//int *dc, *dc;
int dcnt;

int nseg;
//int *dc, *didx;
//int *d;
//TYPE *dv;

int sc;

double time_in_mill_now() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double time_in_mill =
	(tv.tv_sec) * 1000.0 + (tv.tv_usec) / 1000.0;
  return time_in_mill;
}

int compare1(const void *a, const void *b)
{
	if (((struct v_struct *)a)->grp - ((struct v_struct *)b)->grp > 0) return 1;
	if (((struct v_struct *)a)->grp - ((struct v_struct *)b)->grp < 0) return -1;
	if (((struct v_struct *)a)->row - ((struct v_struct *)b)->row > 0) return 1;
	if (((struct v_struct *)a)->row - ((struct v_struct *)b)->row < 0) return -1;
	return ((struct v_struct *)a)->col - ((struct v_struct *)b)->col;
}
//temp_v[i].grp = temp_v[i].row / SM_WIDTH;
int compare2(const void *a, const void *b)
{
	if (((struct v_struct *)a)->grp - ((struct v_struct *)b)->grp > 0) return 1;
	if (((struct v_struct *)a)->grp - ((struct v_struct *)b)->grp < 0) return -1;
	if (((struct v_struct *)a)->col - ((struct v_struct *)b)->col > 0) return 1;
	if (((struct v_struct *)a)->col - ((struct v_struct *)b)->col < 0) return -1;
	return ((struct v_struct *)a)->row - ((struct v_struct *)b)->row;
}


void ready(int argc, char **argv)
{
	FILE *fp;
	int *loc;
	char buf[300];
	int nflag, sflag;
	int dummy, pre_count=0, tmp_ne;
	int i;
	
	srand(time(NULL));
	
	//sc = atoi(argv[2]);
	fp = fopen(argv[1], "r");
	fgets(buf, 300, fp);
	if(strstr(buf, "symmetric") != NULL || strstr(buf, "Hermitian") != NULL) sflag = 1; // symmetric
	else sflag = 0;
	if(strstr(buf, "pattern") != NULL) nflag = 0; // non-value
	else if(strstr(buf, "complex") != NULL) nflag = -1;
	else nflag = 1;
	
	#ifdef SYM
		sflag = 1;
	#endif
	
	while(1) {
		pre_count++;
		fgets(buf, 300, fp);
		if(strstr(buf, "%") == NULL) break;
	}
	fclose(fp);
	
	fp = fopen(argv[1], "r");
	for(i=0;i<pre_count;i++)
		fgets(buf, 300, fp);
	
	fscanf(fp, "%d %d %d", &nr, &nc, &ne);
	nr0 = nr;
	ne *= (sflag+1);
	nr = CEIL(nr, BF)*BF;
	nc = CEIL(nc, BF)*BF;
	
//	cout<<endl;
	cout<<"Matrix "<<argv[1]<<" nr="<<nr<<" nc="<<nc<<" ne="<<ne<<endl;
	
	temp_v = (struct v_struct *)malloc(sizeof(struct v_struct)*(ne+1));
	gold_temp_v = (struct v_struct *)malloc(sizeof(struct v_struct)*(ne+1));
	
	for(i=0;i<ne;i++) {
		fscanf(fp, "%d %d", &temp_v[i].row, &temp_v[i].col);
		temp_v[i].grp = 0;
		temp_v[i].row--; temp_v[i].col--;
		
		if(temp_v[i].row < 0 || temp_v[i].row >= nr || temp_v[i].col < 0 || temp_v[i].col >= nc) {
			fprintf(stdout, "A vertex id is out of range %d %d\n", temp_v[i].row, temp_v[i].col);
			exit(0);
		}
		if(nflag == 0) temp_v[i].val = (TYPE)(rand()%1048576)/1048576;
		else if(nflag == 1) {
			TYPE ftemp;
			fscanf(fp, " %lf ", &ftemp);
			temp_v[i].val = ftemp;
		} else { // complex
			TYPE ftemp1, ftemp2;
			fscanf(fp, " %lf %lf ", &ftemp1, &ftemp2);
			temp_v[i].val = ftemp1;
		}
		if(sflag == 1) {
			i++;
			temp_v[i].row = temp_v[i-1].col;
			temp_v[i].col = temp_v[i-1].row;
			temp_v[i].val = temp_v[i-1].val;
			temp_v[i].grp = 0;
		}
	}
	qsort(temp_v, ne, sizeof(struct v_struct), compare2);
	
	loc = (int *)malloc(sizeof(int)*(ne+1));
	
	memset(loc, 0, sizeof(int)*(ne+1));
	loc[0]=1;
	for(i=1;i<ne;i++) {
		if(temp_v[i].row == temp_v[i-1].row && temp_v[i].col == temp_v[i-1].col)
			loc[i] = 0;
		else loc[i] = 1;
	}
	for(i=1;i<=ne;i++)
		loc[i] += loc[i-1];
	for(i=ne; i>=1; i--)
		loc[i] = loc[i-1];
	loc[0] = 0;
	
	for(i=0;i<ne;i++) {
		temp_v[loc[i]].row = temp_v[i].row;
		temp_v[loc[i]].col = temp_v[i].col;
		temp_v[loc[i]].val = temp_v[i].val;
		temp_v[loc[i]].grp = temp_v[i].grp;
	}
	ne = loc[ne];
	temp_v[ne].row = nr;
	gold_ne = ne;
	for(i=0;i<=ne;i++) {
		gold_temp_v[i].row = temp_v[i].row;
		gold_temp_v[i].col = temp_v[i].col;
		gold_temp_v[i].val = temp_v[i].val;
		gold_temp_v[i].grp = temp_v[i].grp;
	}
	free(loc);

	// for(int i =0; i< 100 ; i++)
	// {
	//     cout<<" gold_temp_v[i].val "<<gold_temp_v[i].val<<" temp_v[i].val "<<temp_v[i].val<<endl;
	// }
}

#endif
