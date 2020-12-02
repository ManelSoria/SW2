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


int kkfu(map *M,int TESTtheTEST) {

	if (M->nd==3) CRASH("sorry has to be 2d");

	int info=0;
  	double *x,*y;
  	int i,j,pt,pta=1;

	if (info>1) printmap("test_sweeps_2d",M);
	x=dmem(M);
	y=dmem(M);

  	setzero_scaf(x,M); // x will be filled in with a single sweep
	setzero_scaf(y,M); // y will be filled in with the partial sweeps

	forall(i,j,M)
		XX(i,j)=F2d(i,j);

	print_scaf(x,"x",0,M,"%+10.3e ");

	pprintf("lsy=%d ley=%d \n",lsy(M),ley(M));

	for (j=lsy(M);j<=ley(M);j++)
		for (i=lsx(M);i<=lex(M);i++) {
			YY(i,j)=F2d(i,j);
		}

	print_scaf(y,"y",0,M,"%+10.3e ");

	pprintf("normdif=%f \n",norm_dif_scaf(x,y,0,M));

  	free(x);
  	free(y);


  	int px,py; // processor coordinates
  	for (py=1;py<=NPY(M);py++)
  		for (px=1;px<=NPX(M);px++)
  			pprintf("px=%d py=%d \n",px,py);

  	if (1)
	  	return(1);
	else
	  	return(0);
}


double rZon(double rE,double epsilon2,double lat) {
	double rZon;
    rZon = ((rE)/sqrt(1.0+(tan(lat)*tan(lat))/(epsilon2)));
    return (rZon);
}

double rMed(double rE,double epsilon2,double lat) {
	double rMed;
    rMed = (rE)/(epsilon2)*pow( rZon(rE,epsilon2,lat) / (rE*cos(lat)),3 );
    return(rMed);
}



double lon0=0; // angular coordinate of whatever origin... explicar be !!
double lat0=0; // angular coordinate of whatever

double dlon=1*2*3.1415927/360;
double dlat=1*2*3.1415927/360;

#define LONE(i,j) lon0+(i-1)*dlon
#define LATE(i,j) lat0+(j-1)*dlat

int test_sweeps2_2d(map *M,int TESTtheTEST) {

	if (M->nd==3) CRASH("sorry has to be 2d");

	int info=0;
  	double *x,*y;
  	int i,j;

	if (info>1) printmap("test_coords",M);
	x=dmem(M);
	y=dmem(M);

  	setzero_scaf(x,M); // x will be filled in with a single sweep
	setzero_scaf(y,M); // y will be filled in with the partial sweeps


	double epsilon2=1.10;
	double rE=71492000; // eq radius, m
	double lat;
	double lon;
	double rZ;
	double rM;
	double posx,posy;
	double deltax,deltay;

	double startx,starty;

	posy=0;
	for (j=asy(M);j<=aey(M);j++) { // ho correm absolutament tot
		posx=0;
		for (i=asx(M);i<=aex(M);i++) {
			lat=LATE(i,j);
			lon=LONE(i,j);
			// ara som a lon,lat
			// calculem els radis aqui
			rZ=rZon(rE,epsilon2,lat); // tendeix a 0 a mida que ens acostem al pol
			rM=rMed(rE,epsilon2,lat);
			deltax=dlon*rZ; // distancies de la nostra cel.la
			deltay=dlat*rM;
			posx=posx+deltax;
			pprintf("i=%d j=%d lat=%e lon=%e posx=%e posy=%e \n",i,j,lat,lon,posx,posy);
			// save the position of nodes at 0,0 per despres restar
			if (i==0) startx=posx;
			if (j==0) starty=posy;
			ifownedh(i,j,M) { // si es meva la guardo
 				XX(i,j)=posx;
				YY(i,j)=posy;
			}

		}
		posy=posy+deltay;
	}

	pprintf("startx=%e starty=%e \n",startx,starty);

	forallh(i,j,M) {
		XX(i,j)=XX(i,j)-startx;
		YY(i,j)=YY(i,j)-starty;
	}

	print_scaf(x,"x",1,M,"%+10.3e ");
	//print_scaf(y,"y",1,M,"%+10.3e ");


  	free(x);
  	free(y);


  	int px,py; // processor coordinates
  	for (py=1;py<=NPY(M);py++)
  		for (px=1;px<=NPX(M);px++)
  			pprintf("px=%d py=%d \n",px,py);

  	if (1)
	  	return(1);
	else
	  	return(0);
}
