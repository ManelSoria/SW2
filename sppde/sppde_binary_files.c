/*
This file is part of SW2 code

2018-2020
Manel Soria, Arnau Prat, Arnau Sabates, Marc Andres-Carcasona, Arnau Miro, Enrique Garcia-Melendo
UPC - ESEIAAT - TUAREG

(c) Manel Soria, Enrique Garcia-Melendo 2018-2020

LICENSED UNDER: Attribution 4.0 International

*/
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

//#define TESTtheTEST // see below

int pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,map *M); // todo: this has to be in sppde.h
int un_pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,map *M);

#define F2d(i,j) (i)+(j)*1000 // better keep the test macros at local file level

// proc 0 gathers x and then saves as binary; file must be open
// all processors MUST call this function

#define MAGIC_SPPDE 15

void fwrite_scaf(double *x, FILE *f,map *M) {
	int info=0;
  int magic=MAGIC_SPPDE;

	if (M->nd==3) CRASH("3d not implemented yet sorry.. ");

	map MJ;
	double *xg;
	xg= alloc_gather(x,M,&MJ,0);

	if (info>1) printmap("global",&MJ);
	if (quisoc()==0 && info>1)
		print_scaf(xg,"global",1,&MJ,"%+10.3e ");
	int nw;
	if (quisoc()==0) {
    nw=fwrite( (void *)&magic,sizeof(int)   ,1   ,f); // save magic
    if (nw!=1)
      CRASH("oops problems writing binary file magic number");
		nw=fwrite( (void *)xg    ,sizeof(double),MJ.S,f);
		if (nw!=MJ.S)
			CRASH("oops problems writing binary file MJ.S=%d nw=%d",MJ.S,nw);
	}

	free(xg);
}

// proc 0 gathers x and then prints field x to file fname, that is created
// all processors MUST call this function

void text_write_scaf(double *x, char *fname,map *M) {
  int info=0;

  if (M->nd==3) CRASH("3d not implemented yet sorry.. ");

  map MJ;
  double *xg;
  xg= alloc_gather(x,M,&MJ,0);

  if (info>1) printmap("global",&MJ);
  if (quisoc()==0 && info>1)
    print_scaf(xg,"global",1,&MJ,"%+10.3e ");
  int nw;
  if (quisoc()==0) {
    FILE *fo;
    int i,j;
    fo=fopen(fname,"w");
    if (fo==NULL) CRASH("Can't open <%s>",fname);
    fprintf(fo,"%d %d %d %d\n",gsx(&MJ),gex(&MJ),gsy(&MJ),gey(&MJ)); // write dimensions
    forall(i,j,&MJ)
      fprintf(fo,"%.20e\n",ac(x,i,j,&MJ)); // write all the field
    fclose(fo);
  }

  free(xg);
}

// proc 0 loads x (file must be open) and scatters to all proc
// halo is updated
// all processors MUST call this function

void text_read_scaf(double *x, char *fname, map *M) {
  int info = 1;

  map MJ;
  joinM(M,&MJ); // create global map
  double *xg;

  if (quisoc()==0) {

    FILE *fo;
    int i,j;
    fo=fopen(fname,"r");
    if (fo==NULL) CRASH("Can't open <%s>",fname);
    int vgsx, vgex, vgsy, vgey;
    int matched;
    matched = fscanf(fo, "%d %d %d %d\n", &vgsx, &vgex, &vgsy, &vgey);
    if(matched!=4) CRASH("text_read_scaf incorrect or corrupt file");

    if(vgex<vgsx) CRASH("text_read_scaf global x axis starting index %d is greater than end index %d", vgsx, vgex);
    if(vgey<vgsy) CRASH("text_read_scaf global y axis starting index %d is greater than end index %d", vgsy, vgey);

    if (vgsx!=gsx(&MJ)) CRASH("vgsx=%d should be %d",vgsx,gsx(&MJ));
    if (vgex!=gex(&MJ)) CRASH("vgex=%d should be %d",vgex,gex(&MJ));
    if (vgsy!=gsy(&MJ)) CRASH("vgsy=%d should be %d",vgsy,gsy(&MJ));
    if (vgey!=gey(&MJ)) CRASH("vgey=%d should be %d",vgey,gey(&MJ));

    xg=dmem(&MJ); // proc 0 allocs memory for the global map

    int value;
    forall(i,j,&MJ) {
      value = fscanf(fo,"%lf\n",&(ac(xg,i,j,&MJ))); // write all the field
      if(value!=1) CRASH("text_read_scaf incorrect or corrupt file");
    }
    fclose(fo);

  }

  gs_field(1,0,x,M,xg,&MJ); // scatter xg (MJ) to x (M); don't transfer external halos

  if (quisoc()==0)
    free(xg);

  halo_update(x,M);

}


// proc 0 loads x (file must be open) and scatters to all proc
// halo is updated
// all processors MUST call this function

void fread_scaf(double *x, FILE *f,map *M) {
	int info=0;
  int magic=0;
	if (M->nd==3) CRASH("3d not implemented yet sorry.. ");
	map MJ;
	joinM(M,&MJ); // create global map
	double *xg;


	int nw;

	if (quisoc()==0) {

    nw=fread( (void *)&magic, sizeof(int), 1 , f);
    if (nw!=1)
      CRASH("oops problems reading binary file magic number");
    if (magic!=MAGIC_SPPDE) CRASH("Wrong binary field magic, is %d and should be %d",magic,MAGIC_SPPDE);

		xg=dmem(&MJ); // proc 0 allocs memory for the global map
		nw=fread( (void *)xg,sizeof(double),MJ.S,f);
		if (nw!=MJ.S)
			CRASH("oops problems reading binary file MJ.S=%d nw=%d",MJ.S,nw);
    if (info) pprintf("fread_scaf read nw=%d doubles ok \n",nw);
		if (info>1) print_scaf(xg,"global_read",1,&MJ,"%+10.3e ");

	}

	gs_field(1,0,x,M,xg,&MJ); // scatter xg (MJ) to x (M); don't transfer external halos

	if (quisoc()==0)
		free(xg);

	halo_update(x,M);
}


/// ******************************************************************************

// generates a field x, saves it and reads as y
// tets is ok if x==y
int test_binaryfiles_2d(map *M,int TESTtheTEST) { // 2d
  int info=0;
  if (M->nd==3) CRASH("3d not implemented yet sorry.. ");

  double *x,*y;
  x=dmem(M);
  y=dmem(M);

  setzero_scaf(x,M);
  setzero_scaf(y,M);

  int i,j;

  forall(i,j,M)
  	ac(x,i,j,M)=F2d(i,j);

  halo_update(x,M);
  if (info>1) print_scaf(x,"x",1,M,"%+10.3e ");

  FILE *fo;
  if (quisoc()==0) {
	fo=fopen("bintest.bin","w"); // must be called only by 0 !!
  	if (fo==NULL) CRASH("can't open bintest.bin to write");
  }

  fwrite_scaf(x,fo,M); // must be called by all !!

  if (quisoc()==0) // must be called only by 0 !!
  	 fclose(fo);

  if (quisoc()==0) {
	fo=fopen("bintest.bin","r");
  	if (fo==NULL) CRASH("can't open bintest.bin to read");
  }

  fread_scaf(y,fo,M); // must be called by all !!

  if (quisoc()==0) // must be called only by 0 !!
  	 fclose(fo);

  if (info>1) {
    print_scaf(x,"x",1,M,"%+10.3e ");
    print_scaf(y,"y",1,M,"%+10.3e ");
  }

  // IMPORTANT: IT SHOULD BE A DELIBERATE ERROR TO TEST THE TEST
  if (TESTtheTEST )
	  ifowned(3,3,M) ac(x,3,3,M)=33;

  forallh(i,j,M) {
  	ac(x,i,j,M)=ac(x,i,j,M)-ac(y,i,j,M);
  }
  double maxabs=norm_scaf(x,0,M);
  if (info) pprintf("maxabs=%e \n",maxabs);
  free(x);
  free(y);

  if (maxabs==0)
	  return(1); // okk :)
  else
  	  return(0); // not ok :(
}
