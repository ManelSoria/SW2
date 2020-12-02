// sppde: structured parallel pde solver
// Manel Soria 2017 - UPC - ESEIAAT - TUAREG


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/times.h>
#include <stdio.h>
#include <stdarg.h>
#include "mpi.h"
#include "sppde.h"
#include "sppde_tests.h" // prototypes of test functions

#define F2d(i,j) 1+abs(i)+abs(j)*1000 // non zero value

#define XX(i,j) ac(x,(i),(j),M)
#define YY(i,j) ac(y,(i),(j),M)


// NEW MACROS to replace the previous


// number of nodes (excluding halos)
#define numn2d(M) (M->gl1[1]-M->gl0[1]+1 ) * ( M->gl1[2]-M->gl0[2]+1 )

// if where y is not zero it is equal to x, then returns 1, else 0
static int partial_test(double *x,double *y,map *M) {
	int i,j;
	int r=1;
	int rg;
	int info=0;
	forallh(i,j,M) {
		if (YY(i,j)!=0)
			if (XX(i,j)!=YY(i,j)) r=0;
	}
	if (info) pprintf("partial_test r=%d \n",r);
    MPI_Allreduce(&r, &rg, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

	return(rg);
}

// total number of non-null positions (including halos)
static int count_non_null(double *x,map *M) {
	int i,j;
	int s=0;
	int sg;
	int info=0;
	forallh(i,j,M) {
		if (XX(i,j)!=0)
			s+=1;
	}
	if (info) pprintf("count_non_null s=%d \n",s);
    MPI_Allreduce(&s, &sg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	return(sg);

}

int test_sweeps_2d(map *M,int TESTtheTEST) {

	if (M->nd==3) CRASH("sorry has to be 2d");

	int info=0;
  	double *x,*y;
  	int i,j,pt,pta=1;

	if (info>1) printmap("test_sweeps_2d",M);
	x=dmem(M);
	y=dmem(M);

  	setzero_scaf(x,M); // x will be filled in with a single sweep
	setzero_scaf(y,M); // y will be filled in with the partial sweeps

	forall(i,j,M) XX(i,j)=F2d(i,j);

	// We accumulate Y values, as we should pass exactly once at each point,
	// at the end it has to be equal to X

    // sweep the inner (non BC) areas
	foralli(i,j,M)   YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M);  pta&=pt; if (!pta) pprintf("test_sweeps_2d i !\n");

	// sweep the edges (where there are BC)
	fornorth(i,j,M) YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M); pta&=pt; if (!pta) pprintf("test_sweeps_2d n !\n");
 	forsouth(i,j,M) YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M); pta&=pt; if (!pta) pprintf("test_sweeps_2d s !\n");
	foreast(i,j,M)  YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M); pta&=pt; if (!pta) pprintf("test_sweeps_2d e !\n");
	forwest(i,j,M)  YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M); pta&=pt; if (!pta) pprintf("test_sweeps_2d w !\n");

	// sweep the vertices (where there are BC)
	forne(i,j,M)    YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M); pta&=pt; if (!pta) pprintf("test_sweeps_2d ne !\n");
	fornw(i,j,M)    YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M); pta&=pt; if (!pta) pprintf("test_sweeps_2d nw !\n");
	forse(i,j,M)    YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M); pta&=pt; if (!pta) pprintf("test_sweeps_2d se !\n");
	forsw(i,j,M)    YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M); pta&=pt; if (!pta) pprintf("test_sweeps_2d sw !\n");

	if (TESTtheTEST) {
	foralli(i,j,M)   YY(i,j)=YY(i,j)+F2d(i,j);   pt=partial_test(x,y,M);  pta&=pt;
	}

	if (info) pprintf("test_sweeps_2d pta=%d \n",pta);

	// with pt we check that no point has been sweept more than once

	if (info>1) print_scaf(y,"y",0,M,"%+10.3e ");

	// with nn we check that all points have been sweept
	// if domain is periodic, limiting nodes are included in foralli

	int nn=count_non_null(y,M);
	if (info) pprintf("test_sweeps_2d non_null in y = %d should be: %d \n",nn,numn2d(M));

  	free(x);
  	free(y);

  	if (pta && nn==numn2d(M))
	  	return(1);
	else
	  	return(0);
}
