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
#include "sppde_internals.h"

// given two fields x and y
//   x is distributed x according to map M
//   y is global, according to MJ (consistent with M), only available to proc 0
// m==0: gather  proc 0 receives all the field x and stores in the preallocated global field y
// m==1: scatter proc 0 sends global field y to the procs, that store it in their prealloc areas of x


// eh controls the behaviour of the halos
// if eh==0, halos are not gathered nor scattered
// if eh==1, external halos are gathered to 0, and **all** the halos are scattered

void gs_field(int m,int eh,double *x,map *M,double *y,map *MJ) {
	int info = 0;

	if (info) pprintf("gs_field entering mode m=%d \n",m);
	if (M->nd==3) CRASH("3d not implemented yet sorry.. ");
	if (m!=0 && m!=1) CRASH("gs_field m=%d and should be 0 or 1",m);
	if (quisoc()==0) {
		if (M->sh != MJ->sh) CRASH("gs_field halo sizes should be equal M->sh=%d MJ->sh=%d",M->sh,MJ->sh);
	}

	if (info>1) {
		printmap("gs_field M ",M);
		if (quisoc()==0) printmap("gs_field MJ ",MJ);
	}

	// transfer from 0 to 0
	if (quisoc()==0) {
		double *bb; // buffer
		int H=0;
		if (eh) H=M->sh; // we assume that both halo sizes are equal
		bb=(double *)mallocc(sizeof(double)*(M->l1[1]-M->l0[1]+1+2*H)*(M->l1[2]-M->l0[2]+1+2*H) );
		if (info) pprintf("gs_field 0/1\n");
		if (m==0) {
			if (info) pprintf("gs_field 0/2\n");
			pack(   bb,0,x, M->l0[1]-H ,M->l1[1]+H ,M->l0[2]-H ,M->l1[2]+H ,M);  // x to bb
			if (info) pprintf("gs_field 0/21\n");
			un_pack(bb,0,y, M->l0[1]-H ,M->l1[1]+H ,M->l0[2]-H ,M->l1[2]+H ,MJ); // bb to y
		} else {
			if (info) pprintf("gs_field 0/3\n");
			pack(   bb,0,y, M->l0[1]-H,M->l1[1]+H,M->l0[2]-H,M->l1[2]+H,MJ); // y to bb
			if (info) pprintf("gs_field 0/31\n");
			un_pack(bb,0,x, M->l0[1]-H,M->l1[1]+H,M->l0[2]-H,M->l1[2]+H ,M); // bb to x
		}
		free(bb);

		if (info) pprintf("gs_field /00\n");

	}

	if (info) pprintf("gs_field/0 done\n");
	if (quisoc()==0) {
		for (int pp=1;pp<=M->gpc[1]*M->gpc[2]-1;pp++) {
			if (info) pprintf("gs_field receive pp=%d \n",pp);
			int coord[4];
			MPI_Status st;
			int inix,fix,iniy,fiy;
			inix=M->al0[pp][1]; fix=M->al1[pp][1]; iniy=M->al0[pp][2]; fiy=M->al1[pp][2]; // area to be transfered
			if (eh) { // if halos have to be included
				if (m==0) {  // for gathering, include only external halos
					if (M->pco[pp][1]==1)         inix-=M->sh;
					if (M->pco[pp][1]==M->gpc[1]) fix +=M->sh;
					if (M->pco[pp][2]==1)         iniy-=M->sh;
					if (M->pco[pp][2]==M->gpc[2]) fiy +=M->sh;
				}
				else { // for scattering, include all halos
					inix-=M->sh;
					fix +=M->sh;
					iniy-=M->sh;
					fiy +=M->sh;
				}
			}
			if (info) pprintf("gs_field/part1 pp=%d area= x:%d %d y:%d %d \n",pp,inix,fix,iniy,fiy);
			int sz=(fix-inix+1)*(fiy-iniy+1);
			if (info) pprintf("gs_field inix=%d fix=%d iniy=%d fiy=%d sz=%d \n",inix,fix,iniy,fiy,sz);
			double *bb;
			bb=mallocc(sizeof(double)*sz);
			if (m==0) { // reveive from proc pp and unpack to global y
				checkr(MPI_Recv(bb,sz,MPI_DOUBLE,pp,2,MPI_COMM_WORLD,&st),"gs_2");
				un_pack(bb,0,y,inix,fix,iniy,fiy,MJ);
			} else { // pack from global y and send to proc pp
				pack   (bb,0,y,inix,fix,iniy,fiy,MJ);
				checkr(MPI_Ssend(bb,sz,MPI_DOUBLE,pp,2,MPI_COMM_WORLD),"gs_3");
			}
			free(bb);
			MPI_Barrier(MPI_COMM_WORLD);
		}
	} else {
		for (int pp=1;pp<=M->gpc[1]*M->gpc[2]-1;pp++) {
			if (quisoc()==pp) {
				MPI_Status st;
				if (info) pprintf("gather_field send pp=%d \n",pp);
				int inix,fix,iniy,fiy;
				inix=M->al0[pp][1]; fix=M->al1[pp][1]; iniy=M->al0[pp][2]; fiy=M->al1[pp][2]; // area to be transfered
				if (eh) { // if halos have to be included
					if (m==0) {  // for gathering, include only external halos
						if (M->pco[pp][1]==1)         inix-=M->sh;
						if (M->pco[pp][1]==M->gpc[1]) fix +=M->sh;
						if (M->pco[pp][2]==1)         iniy-=M->sh;
						if (M->pco[pp][2]==M->gpc[2]) fiy +=M->sh;
					}
					else { // for scattering, include all halos
						inix-=M->sh;
						fix +=M->sh;
						iniy-=M->sh;
						fiy +=M->sh;
					}
				}
				if (info) pprintf("gs_field/part2 pp=%d area= x:%d %d y:%d %d \n",pp,inix,fix,iniy,fiy);
				int sz=(fix-inix+1)*(fiy-iniy+1);
				if (info) pprintf("gs_field inix=%d fix=%d iniy=%d fiy=%d sz=%d \n",inix,fix,iniy,fiy,sz);
				double *bb;
				bb=mallocc(sizeof(double)*sz);
				if (m==0) { // pack from distributed x and send to 0
					pack(bb   ,0,x,inix,fix,iniy,fiy,M);
					checkr(MPI_Ssend(bb,sz,MPI_DOUBLE,0,2,MPI_COMM_WORLD),"gs_5");
				} else { // receive from 0 and unpack in global y
					checkr(MPI_Recv(bb,sz,MPI_DOUBLE,0,2,MPI_COMM_WORLD,&st),"gs_6");
					un_pack(bb,0,x,inix,fix,iniy,fiy,M);
				}
				free(bb);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	if (info) pprintf("gs_field end\n");

}




#define F2d(i,j) 1+abs(i)+abs(j)*100 // better keep the test macros at local file level

// generates a distributed field x, gathers it to gx, and then scatters it to y
// tets is ok if x==y
// also tests alloc_gather
// tests it with eh=1 (halos transferred) or eh==0

int test_gs_field_2d_eh(int eh,map *M,int TESTtheTEST) { // 2d
  int info=1;
  if (M->nd==3) CRASH("3d not implemented yet sorry.. ");
  map MJ;

  if (info) pprintf("test_gs_field_2d eh=%d \n",eh);

  if (info>1) printmap("test_gs_field_2d ",M);

  double *x,*gx,*y;
  x=dmem(M);
  y=dmem(M);

  setzero_scaf(x,M);
  setzero_scaf(y,M);

  int i,j;

  forallh(i,j,M)
  	ac(x,i,j,M)=F2d(i,j);

  gx=alloc_gather(x,M,&MJ,eh); // this tests gs_field(0) (gather from all procs to 0, external halos gathered also)

  // remember: now only proc 0 can acess gx !!

  if (info>1) print_scaf(x,"x",1,M,"%+10.3e ");

  if (info>1) {
  	if (quisoc()==0)
	  	print_scaf(gx,"gx",1,&MJ,"%+10.3e ");
  }


  gs_field(1,eh,y,M,gx,&MJ); // from gx (local to 0) to y (distributed)

  if (TESTtheTEST) {
	  ifowned(1,1,M) ac(y,1,1,M)=ac(y,1,1,M)+1;
  }

  if (info>1) print_scaf(y,"y",1,M,"%+10.3e ");

  int ndif=0,ndifa;
  forallh(i,j,M) {
  	if (i>=M->l0[1] && i<=M->l1[1] && j>=M->l0[2] && j<=M->l1[2] ) { // in the inner areas..
		if (ac(x,i,j,M)!=ac(y,i,j,M)) ndif++; // y should be equal to x
  	} else { // in the halo areas..
  		if (eh==1) {
  			if (ac(x,i,j,M)!=ac(y,i,j,M)) ndif++; // y should be equal to x if eh==1
  		} else {
  			if (ac(y,i,j,M)!=0) ndif++; // else, y should be 0
  		}

  	}
  }
  MPI_Allreduce(&ndif, &ndifa, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


  if (info) pprintf("ndif=%d ndifa=%d \n",ndif,ndifa);

  free(x);
  free(gx);
  free(y);

  if (ndifa==0)
	  return(1); // okk :)
  else
  	  return(0); // not ok :(
}


int test_gs_field_2d(map *M,int TESTtheTEST) {
	int r1,r2,r;
	r1=test_gs_field_2d_eh(0,M,TESTtheTEST);
	r2=test_gs_field_2d_eh(1,M,TESTtheTEST);
	if (r1 && r2) r=1; else r=0;
	return(r);
}
