// sppde: structured parallel pde solver
// Manel Soria 2017 - UPC - ESEIAAT - TUAREG

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "sppde.h"

#include "sppde_internals.h"



// buffer b should be already allocated
// returns q, the position where the next field has to be packed
int pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,map *M) {
	int i,j,q;
	int info=0;
	if (M->nd==3) CRASH("has to be 2d");
	if (info) pprintf("pack entering q0=%d  i1=%d i2=%d j1=%d j2=%d \n",q0,i1,i2,j1,j2);
	q=q0;
	for (j=j1;j<=j2;j++)
		for (i=i1;i<=i2;i++) {
			b[q]=X(i,j,M);
			if (info>1) pprintf("pack i=%d j=%d %e \n",i,j,X(i,j,M));
			q++;
		}
	if (info) pprintf("pack leaving q=%d\n",q);
	return(q);
}

int un_pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,map *M) {
	int i,j,q;
	int info=0;
	if (M->nd==3) CRASH("has to be 2d");
	if (info) pprintf("un_pack entering q0=%d  i1=%d i2=%d j1=%d j2=%d \n",q0,i1,i2,j1,j2);
	q=q0;
	for (j=j1;j<=j2;j++)
		for (i=i1;i<=i2;i++)
			X(i,j,M)=b[q++];
	return(q);
}

// ====================0 3d

// buffer b should be already allocated
// returns q, the position where the next field has to be packed
int Pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,int k1,int k2,map *M) {
	int i,j,k,q;
	int info=0;
	if (M->nd==2) CRASH("has to be 3d");
	if (info) pprintf("Pack entering (3d) q0=%d  i1=%d i2=%d j1=%d j2=%d k1=%d k2=%d \n",q0,i1,i2,j1,j2,k1,k2);
	q=q0;
	for (k=k1;k<=k2;k++)
		for (j=j1;j<=j2;j++)
			for (i=i1;i<=i2;i++) {
				b[q]=X3(i,j,k,M);
				if (info>1) pprintf("pack i=%d j=%d k=%d %e \n",i,j,X3(i,j,k,M));
				q++;
			}
	if (info) pprintf("pack leaving q=%d\n",q);
	return(q);
}

int Un_pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,int k1,int k2,map *M) {
	int i,j,k,q;
	int info=0;
	if (M->nd==2) CRASH("has to be 3d");
	if (info) pprintf("Un_pack entering (3d) q0=%d  i1=%d i2=%d j1=%d j2=%d k1=%d k2=%d \n",q0,i1,i2,j1,j2,k1,k2);
	q=q0;
	for (k=k1;k<=k2;k++)
		for (j=j1;j<=j2;j++)
			for (i=i1;i<=i2;i++)
				X3(i,j,k,M)=b[q++];
	return(q);
}

int g_pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,int k1,int k2,map *M) {
	int r;
	if (M->nd==2)
		r=pack(b,q0,x,i1,i2,j1,j2,M);
	else
		r=Pack(b,q0,x,i1,i2,j1,j2,k1,k2,M);
	return(r);
}

int g_un_pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,int k1,int k2,map *M) {
	int r;
	if (M->nd==2)
		r=un_pack(b,q0,x,i1,i2,j1,j2,M);
	else
		r=Un_pack(b,q0,x,i1,i2,j1,j2,k1,k2,M);
	return(r);
}

// =======================
void easy_s(int nb, double *b, int ndata) {
		int info=1;
		if (info) pprintf("easy_s nb=%d ndata=%d \n",nb,ndata);
		if (info>1) {
			for (int i=0;i<=ndata-1;i++)
				pprintf("easy_s %d %f \n",i,b[i]);
		}
 		if (nb==-1) return;
		checkr(MPI_Ssend(b,ndata,MPI_DOUBLE,nb,1,MPI_COMM_WORLD),"easy_s");
}

void easy_r(int nb,double *b,int ndata) {
	MPI_Status st;

	int info=1;
	if (info) pprintf("easy_r nb=%d ndata=%d \n",nb,ndata);
	if (nb==-1) { for (int i=0;i<=ndata-1;i++) b[i]=0.0; return;}
	checkr(MPI_Recv(b,ndata,MPI_DOUBLE,nb,1,MPI_COMM_WORLD,&st),"easy_r");
	if (info>1) {
		for (int i=0;i<=ndata-1;i++)
			pprintf("easy_r %d %f \n",i,b[i]);
	}
}


void easy_sr(int nb,double *bs,double *br,int ndata){

/*

Replaces the old easy_s and easy_r at the same time so processors can exchange information in a full-duplex network. This means that the processors can be senders and receivers at the same time.

- double nb (input) The neighbour to exchange information.
- double *bs (input/output) The sending [previously allocatted] buffer.
- double *br (input/output) The receiving [previously allocatted] buffer.
- int ndata (input) The number of elements to be sent or received. The number of elements that are exchanged are always the same for the sending and receiving buffers.

This function makes use of MPI_Sendrecv.

*/

    MPI_Status st;

    int info=0;
    int i;
    if (info) pprintf("easy_sr nb=%d ndata=%d \n",nb,ndata);
    if (nb==-1) {
        for (int i=0;i<=ndata-1;i++) br[i]=0.0;
        return;
    }
    if (info>1) {
        for (i=0;i<=ndata-1;i++)
            pprintf("easy_sr send_op %d %f \n",i,bs[i]);
    }
    checkr(MPI_Sendrecv(bs,ndata,MPI_DOUBLE,nb,1,br,ndata,MPI_DOUBLE,nb,1,MPI_COMM_WORLD,&st),"easy_sr");
    if (info>1) {
        for (i=0;i<=ndata-1;i++)
            pprintf("easy_sr recv_op %d %f \n",i,br[i]);
    }
}

void halo_update_y(double *x,map *M) {

/*

Updates the processor's north and south halos when the number of processors in the column is >1.

- double *x (input/output) The matrix/field to be updated.
- map *M (input) The scaf's initrinsic properties.

*/

    double *bbs, *bbr;
    int info=0;
    int nd=M->nd; // QRAT -> int nd=M->nd;??

    int inix=M->l0[1]; // only inner halo (NOT diagonal)
    int fix=M->l1[1];
    int iniz=M->l0[3];
    int fiz=M->l1[3];
    int SH=M->sh*(fix-inix+1);
    if (nd==3) SH=SH*(fiz-iniz+1);

    if (info) pprintf("halo_updatey M->pc[2]=%d inix=%d fix=%d SH=%d\n",M->pc[2],inix,fix,iniz,fiz,SH);

    bbs=(double *)mallocc(sizeof(double)*SH);
    bbr=(double *)mallocc(sizeof(double)*SH);

    if (M->pc[2] % 2 ==1) {
        if (info) pprintf("halo_updatey 1: odd sends top receives top");
        g_pack(bbs,0,x,inix,fix,M->l1[2]-M->sh+1,M->l1[2],iniz,fiz,M);
        easy_sr(M->nb_n,bbs,bbr,SH);
        g_un_pack(bbr,0,x,inix,fix,M->l1[2]+1,M->l1[2]+M->sh,iniz,fiz,M);

        if (info) pprintf("halo_updatey 2: odd sends bottom receives bottom");
        g_pack(bbs,0,x,inix,fix,M->l0[2],M->l0[2]+M->sh-1,iniz,fiz,M);
        easy_sr(M->nb_s,bbs,bbr,SH);
        g_un_pack(bbr,0,x,inix,fix,M->l0[2]-M->sh,M->l0[2]-1,iniz,fiz,M);

    } else {
        if (info) pprintf("halo_updatey 1: odd sends bottom receives bottom");
        g_pack(bbs,0,x,inix,fix,M->l0[2],M->l0[2]+M->sh-1,iniz,fiz,M);
        easy_sr(M->nb_s,bbs,bbr,SH);
        g_un_pack(bbr,0,x,inix,fix,M->l0[2]-M->sh,M->l0[2]-1,iniz,fiz,M);

        if (info) pprintf("halo_updatey 2: odd sends top receives top");
        g_pack(bbs,0,x,inix,fix,M->l1[2]-M->sh+1,M->l1[2],iniz,fiz,M);
        easy_sr(M->nb_n,bbs,bbr,SH);
        g_un_pack(bbr,0,x,inix,fix,M->l1[2]+1,M->l1[2]+M->sh,iniz,fiz,M);

    }

    free(bbs);
    free(bbr);
}

void halo_update_x(double *x,map *M) {

/*

Updates the processor's east and west halos when the number of processors in the row is >1.

- double *x (input/output) The matrix/field to be updated.
- map *M (input) The scaf's initrinsic properties.

*/

    double *bbs, *bbr;
    int info=0;
    int nd=M->nd; // QRAT -> int nd=M->nd;??

    int iniy=M->l0[2]-M->sh; // also diagonal halos !!!
    int fiy=M->l1[2]+M->sh;
    int iniz=M->l0[3];
    int fiz=M->l1[3];
    int SH=M->sh*(fiy-iniy+1);
    if (nd==3) SH=SH*(fiz-iniz+1);

    if (info) pprintf("halo_update M->pc[1]=%d iniy=%d fiy=%d iniz=%d fiz=%d SH=%d\n",M->pc[1],iniy,fiy,iniz,fiz,SH);

    bbs=(double *)mallocc(sizeof(double)*SH);
    bbr=(double *)mallocc(sizeof(double)*SH);

    if (M->pc[1] % 2 ==1 ) {
        if (info) pprintf("halo_updatey 1: odd sends right receives right");
        g_pack(bbs,0,x,M->l1[1]-M->sh+1,M->l1[1],iniy,fiy,iniz,fiz,M);
        easy_sr(M->nb_e,bbs,bbr,SH);
        g_un_pack(bbr,0,x,M->l1[1]+1,M->l1[1]+M->sh,iniy,fiy,iniz,fiz,M);

        if (info) pprintf("halo_updatey 2: odd sends left receives left");
        g_pack(bbs,0,x,M->l0[1],M->l0[1]+M->sh-1,iniy,fiy,iniz,fiz,M);
        easy_sr(M->nb_w,bbs,bbr,SH);
        g_un_pack(bbr,0,x,M->l0[1]-M->sh,M->l0[1]-1,iniy,fiy,iniz,fiz,M);

    } else {
        if (info) pprintf("halo_updatey 1: odd sends left receives left");
        g_pack(bbs,0,x,M->l0[1],M->l0[1]+M->sh-1,iniy,fiy,iniz,fiz,M);
        easy_sr(M->nb_w,bbs,bbr,SH);
        g_un_pack(bbr,0,x,M->l0[1]-M->sh,M->l0[1]-1,iniy,fiy,iniz,fiz,M);

        if (info) pprintf("halo_updatey 2: odd sends right receives right");
        g_pack(bbs,0,x,M->l1[1]-M->sh+1,M->l1[1],iniy,fiy,iniz,fiz,M);
        easy_sr(M->nb_e,bbs,bbr,SH);
        g_un_pack(bbr,0,x,M->l1[1]+1,M->l1[1]+M->sh,iniy,fiy,iniz,fiz,M);

    }
    free(bbs);
    free(bbr);
}

void halo_update_single_map_y(double *x,map *M) {

/*

Updates the processor's north and south halos if in the column there is only one processor.

- double * (input/output) The matrix/field to be updated.
- map *M (input) The scaf's initrinsic properties.

*/

    double *bb;
    int info=0;
    int nd=M->nd; // QRAT -> this has been changed in the other halo functions! before it was int nd;

    int inix=M->l0[1]; // only inner halo (NOT diagonal)
    int fix=M->l1[1];
    int iniz=M->l0[3];
    int fiz=M->l1[3];
    int SH=M->sh*(fix-inix+1);
    if (nd==3) SH=SH*(fiz-iniz+1);

    bb=(double *)mallocc(sizeof(double)*SH);

    // Halo update in y axis

    g_pack(bb,0,x,inix,fix,M->l1[2]-M->sh+1,M->l1[2],iniz,fiz,M);
    g_un_pack(bb,0,x,inix,fix,M->l0[2]-M->sh,M->l0[2]-1,iniz,fiz,M);

    g_pack(bb,0,x,inix,fix,M->l0[2],M->l0[2]+M->sh-1,iniz,fiz,M);
    g_un_pack(bb,0,x,inix,fix,M->l1[2]+1,M->l1[2]+M->sh,iniz,fiz,M);

    free(bb);
} // QRAT

void halo_update_single_map_x(double *x,map *M) {

/*

Updates the processor's east and west halos if in the row there is only one processor.

- double *x (input/output) The matrix/field to be updated.
- map *M (input) The scaf's initrinsic properties.

*/

    double *bb;
    int info=0;
    int nd=M->nd; // QRAT -> this has been changed in the other halo functions! before it was int nd;

    int iniy=M->l0[2]-M->sh; // also diagonal halos !!!
    int fiy=M->l1[2]+M->sh;
    int iniz=M->l0[3];
    int fiz=M->l1[3];
    int SH=M->sh*(fiy-iniy+1);
    if (nd==3) SH=SH*(fiz-iniz+1);

    bb=(double *)mallocc(sizeof(double)*SH);

    // Halo update in x axis

    g_pack(bb,0,x,M->l1[1]-M->sh+1,M->l1[1],iniy,fiy,iniz,fiz,M);
    g_un_pack(bb,0,x,M->l0[1]-M->sh,M->l0[1]-1,iniy,fiy,iniz,fiz,M);

    g_pack(bb,0,x,M->l0[1],M->l0[1]+M->sh-1,iniy,fiy,iniz,fiz,M);
    g_un_pack(bb,0,x,M->l1[1]+1,M->l1[1]+M->sh,iniy,fiy,iniz,fiz,M);

    free(bb);
} // QRAT

void halo_update(double *x,map *M) { // 2d or 3d QRAT

/*

Updates the processor's halos.

- double *x (input/output) The matrix/field to be updated.
- map *M (input) The scaf's initrinsic properties.

This function distinguishes between np=1 and np>1. It could, however, also discriminate cases where npx=1 or npy=1 because currently, periodic boundary conditions do not work for when in the given axis there is only 1 processor.

*/
    if(M->gpc[1]*M->gpc[2]*M->gpc[3]==1) { // QRAT
        halo_update_single_map_y(x,M); // QRAT
        halo_update_single_map_x(x,M); // QRAT
    } else { // QRAT
        halo_update_y(x,M); // has to be the first !!
        halo_update_x(x,M);
    }
}



/*

The function below may be integrated with ac to improve the usage of ShallowWorlds. If it is integrated with ac, the parallel code could be theoretically [almost] identical to sequential code. However, loads of network communications would be carried out and therefore it currently is not a good idea.

*/

double get_val(double *scaf, int i, int j, map *M) {

/*

Returns the value of (i,j) of a scaf stored in memory and defined by M, even if the point is not owned by the processor.

- double *scaf (input) The scalar field stored in memory of which the value is desired.
- int i (input) The row index of the matrix scaf.
- int j (input) The column index of the matrix scaf.
- map *M (input) The scaf's initrinsic properties.

This function makes two network communications. The first one, MPI_Allreduce is needed to know which processor owns the point and therefore, to know the processor that will broadcast the value to the other processors through MPI_Bcast. If this function is integrated into halo_update, halo updating and process comunication will not be seen by the user.

*/

    int info=0;
    int nprocs=quants();
    int *has_value;
    double val;
    int root;
    has_value=(int *)callocc(nprocs, sizeof(int));
    ifowned(i,j,M) {
        *(has_value+quisoc())=1;
        val=ac(scaf,i,j,M);
    }
    MPI_Allreduce(MPI_IN_PLACE, has_value, nprocs, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    int proc;
    for(proc=0;proc<nprocs;proc++)
        if(*(has_value+proc)==1) {
            root=proc;
            break;
        }
    if(info) pprintf("proc: %d\n",proc);
    if(info) pprintf("val: %e\n",val);

    MPI_Bcast(&val, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

    if(info) pprintf("%d\n", root);
    if(info) pprintf("%e\n", val);
    free(has_value);
    return val;

}


// ============================== TEST


 // aux function to test halo update with periodic BC
static int pertest(int i,int ini,int fi,int per) { // function local to this file (TODO)
  int r;
  r=i;
  if (per==1) {
    // make it periodic if needed
   if (i<ini) r=fi-(ini-i)+1;
   if (i>fi)  r=ini+i-fi-1;
  }
//  pprintf("pertest i=%d ini=%d fi=%d per=%d r=%d\n",i,ini,fi,per,r);
  return(r);
}

// aux
static int funfunij(int i,int j,map *M) {
  int r;
  r=pertest(i,M->gl0[1],M->gl1[1],M->per[1])+
    pertest(j,M->gl0[2],M->gl1[2],M->per[2])*10;
  return(r);
}


int test_halo_update_2d(map *M,int TESTtheTEST) {
  int info=0;
	int i,j;
  //CRASH("Epps atencio que cal verificar test_halo_update !!!");

  if (M->nd==3) CRASH("sorry has to be 2d");

	double *x,*y,*z; // three different fields
	x=dmem(M); // allocate the memory for x
	y=dmem(M);
	z=dmem(M);
  // we set all fields to zero (including the halos)
	setzero_scaf(x,M);
  setzero_scaf(y,M);
  setzero_scaf(z,M);

  // X(2,32,M)=33;
// only inside our domain, we define in X a field
  forall(i,j,M) {
    //pprintf("forall i=%d j=%d \n",i,j);
    X(i,j,M)=funfunij(i,j,M); // Arbitray (periodic in x / y) function
  }

  if (TESTtheTEST)
    ifowned(2,2,M)
      X(2,2,M)=456;


  forallh(i,j,M) { // Y is equal, bus also defined in the halos
    //pprintf("forallh i=%d j=%d \n",i,j);
    Y(i,j,M)=funfunij(i,j,M);
  }
  if (info>1) print_scaf(x,"x:forall",1,M,"%+10.3e ");

  if (info>1) print_scaf(y,"y:forallh",1,M,"%+10.3e ");

  // we update the halo of X
  double t0=MPI_Wtime();
  halo_update(x,M);
  double t=MPI_Wtime()-t0;

   // now they should be equal

  //ifowned(4,7,M) X(4,7,M)=333; // introduce an error


  if (info>1) print_scaf(x,"test linop_updated",1,M,"%+10.3e ");

  setzero_scaf(z,M);

  forallh(i,j,M) { // mark error positions in the inner halos with 1
    if (  ( (j>=M->gl0[2] && j<=M->gl1[2]) || M->per[2]==1) &&  // if we are in inner nodes in j .. or periodic
          ( (i>=M->gl0[1] && i<=M->gl1[1]) || M->per[1]==1) && // idem in i .. or periodic
          X(i,j,M) != Y(i,j,M) )
      Z(i,j,M)=1;
  }
  if (info>1) print_scaf(z,"test error",1,M,"%.0f ");

  double ss=0,gss;

  forallh(i,j,M) { // mark error positions in the inner halos with 1
    ss=ss+Z(i,j,M);
  }

   if (info) pprintf("test_halo_update ss=%f \n",ss);

  MPI_Allreduce(&ss, &gss, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // add all errors

  if (info) pprintf("test_halo_update gss=%f time=%f \n",gss,t);

  free(x);
  free(y);
  free(z);

  if (gss==0)
    return(1); // okk!
  else
    return(0);

}
