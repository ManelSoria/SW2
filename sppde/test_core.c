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


#define F2d(i,j) (i)+(j)*1000

#define F3d(i,j,k) (i)+(j)*10+(k)*100

// 2d or 3d
// sets testing field f
void set_field(double *x,int f,map *M) {
  int i,j,k;
  switch(M->nd) {
    case 2: // ================== 2d
      switch(f) {
        case 1: forall(i,j,M) X(i,j,M)=F2d(i,j); break;
        case 2: forall(i,j,M) X(i,j,M)=i+0.01*j; break;
        default: CRASH("uhhh f=%d",f);
      }
    break;

    case 3: // ================== 3d
      switch(f) {
        case 1: Forall(i,j,k,M)  X3(i,j,k,M)=F3d(i,j,k); break;
        default: CRASH("uhhh f=%d",f);
      }
    break;
    default: CRASH("uhhh nd=%d",M->nd);
  }
}


void setvalue(int i,int j,int k,double val,double *x,map *M) {
  if (gOwned(i,j,k,M))
      gac(x,i,j,k,M)=val;
}

void test_operators(map *M) { // 2d or 3d
  double *x,*y,*z;
  x=dmem(M);
  y=dmem(M);
  z=dmem(M);

  setzero_scaf(x,M);
  setzero_scaf(y,M);
  setzero_scaf(z,M);

  setvalue(M->gl0[1],M->gl0[2],M->gl0[3],1,x,M);
  setvalue(M->gl0[1],M->gl1[2],M->gl0[3],2,x,M);
  setvalue(M->gl1[1],M->gl0[2],M->gl0[3],3,x,M);
  setvalue(M->gl1[1],M->gl1[2],M->gl0[3],4,x,M);

  double n2=norm_scaf(x,2,M); // L2
  double ni=norm_scaf(x,0,M); // L inf (maxabs)

  pprintf("test_operators n2=%e \n",n2);
  pprintf("test_operators ni=%e \n",ni);

  free(x);
  free(y);
  free(z);
}


int check_partial_s(char *msg,double *x,double *y,int p,map *M) {
    int i,j,k,r;
    double ss=0,gss;
    Forall(i,j,k,M) {
      if (Y3(i,j,k,M)!=0) {
          if (X3(i,j,k,M)==Y3(i,j,k,M)) // if non zero, it should be equal to the X value
            ss++;
          else
            CRASH("check_partial %s ijk=%d %d %d X=%f Y=%f",msg,i,j,k,X3(i,j,k,M),Y3(i,j,k,M));
      }
    }
    MPI_Allreduce(&ss, &gss, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // add all non-zeros
    r=(int)gss;
    pprintf("check_partial_s %s p=%d r=%d r-p=%d\n",msg,p,r,r-p);
    return(r);
}

void test_partial_sweeps_3d(map *M) {
  double *x,*y;
  int i,j,k,p;

  if (M->nd!=3) CRASH("only 3d maps");
  x=dmem(M); // allocate memory for x
  y=dmem(M);

  setzero_scaf(y,M); // y will be filled in with the partial sweeps
  setzero_scaf(x,M);

  Forallh(i,j,k,M) // loops in all our area, including halos
    X3(i,j,k,M)=1+F3d(i,j,k);

  Foralli(i,j,k,M) Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("Foralli",x,y,0,M);

  ForE(i,j,k,M)    Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForE",x,y,p,M);
  ForW(i,j,k,M)    Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForW",x,y,p,M);
  ForN(i,j,k,M)    Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForN",x,y,p,M);
  ForS(i,j,k,M)    Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForS",x,y,p,M);
  ForT(i,j,k,M)    Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForT",x,y,p,M);
  ForB(i,j,k,M)    Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForB",x,y,p,M);

  ForEN(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForEN",x,y,p,M);
  ForES(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForES",x,y,p,M);
  ForWN(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForWN",x,y,p,M);
  ForWS(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForWS",x,y,p,M);
  ForNT(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForNT",x,y,p,M);
  ForNB(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForNB",x,y,p,M);
  ForST(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForST",x,y,p,M);
  ForSB(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForSB",x,y,p,M);
  ForET(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForET",x,y,p,M);
  ForEB(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForEB",x,y,p,M);
  ForWT(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForWT",x,y,p,M);
  ForWB(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForWB",x,y,p,M);

  ForWSB(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForWSB",x,y,p,M);
  ForWST(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForWST",x,y,p,M);
  ForWNB(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForWNB",x,y,p,M);
  ForWNT(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForWNT",x,y,p,M);
  ForESB(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForESB",x,y,p,M);
  ForEST(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForEST",x,y,p,M);
  ForENB(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("ForENB",x,y,p,M);
  ForENT(i,j,k,M)   Y3(i,j,k,M)=Y3(i,j,k,M)+1+F3d(i,j,k); p=check_partial_s("For",x,y,p,M);

  Forall(i,j,k,M) {
    if (Y3(i,j,k,M)==0) CRASH("uhh? point/s not sweept ijk=%d %d %d ",i,j,k);
  }

  int t=( M->gl1[1]-M->gl0[1] +1 ) * ( M->gl1[2]-M->gl0[2] +1 ) * ( M->gl1[3]-M->gl0[3] +1) ;

  pprintf("test_partial_sweeps_3d t=%d \n",t);
  if (p != t)
    CRASH("uhhh ? p=%d should be %d ",p,t);

  free(y);
  free(x);

  pprintf("test_partial_sweeps_3d OK \n");


}


void test_distributed_arrays_3d(map *M) {
    double *x;
    int i,j,k;

    if (M->nd!=3) CRASH("only 3d maps");

    pprintf("test_distributed_arrays_3d\n");
    printmap("test_distributed_arrays_3d",M);

    test_partial_sweeps_3d(M);

    x=dmem(M); // allocate memory for x

if (0) {
    setzero_scaf(x,M);
    Forallh(i,j,k,M) // loops in all our area, including halos
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"forallh_3d",1,M,"%+10.3e ");

    setzero_scaf(x,M);
    Forall(i,j,k,M) // loops in all our area, NOT including halos
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"forall_3d",1,M,"%+10.3e ");

    setzero_scaf(x,M);
    Foralli(i,j,k,M) // loops in our area, EXCLUDING BOUNDARY CONDITIONS and HALOS
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"foralli_3d",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    Foreast(i,j,k,M)  {// loops east boundary (if owned)
      // pprintf("Foreast %d %d %d\n",i,j,k);
      X3(i,j,k,M)=F3d(i,j,k);
    }
    print_scaf(x,"foreast_3d",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    Forwest(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"forwest_3d",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    Fornorth(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"fornorth_3d",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    Forsouth(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"forsouth_3d",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    Fortop(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"fortop_3d",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    Forbottom(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"forbottom_3d",0,M,"%+10.3e ");
}

if (0) {
    setzero_scaf(x,M);
    ForE(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"ForE",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    ForWS(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"ForWS",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    ForWSB(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"ForWSB",0,M,"%+10.3e ");

}

if (0) {

    setzero_scaf(x,M);
    ForEN(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"ForEN",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    ForES(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"ForES",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    ForWS(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"ForWS",0,M,"%+10.3e ");

    setzero_scaf(x,M);
    ForWN(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"ForWN",0,M,"%+10.3e ");


    setzero_scaf(x,M);
    ForNT(i,j,k,M)
      X3(i,j,k,M)=F3d(i,j,k);
    print_scaf(x,"ForNT",0,M,"%+10.3e ");

}

// all the rest are untested !!!

    free(x);
}

void test_distributed_arrays_2d(map *M) {
    double *x;
    int i,j;

    if (M->nd!=2) CRASH("only 2d maps");


    x=dmem(M); // allocate memory for x

    pprintf("test_distributed_arrays_2d\n");
    printmap("test_distributed_arrays_2d",M);

    setzero_scaf(x,M);
    forallh(i,j,M) // loops in all our area, including halos
      X(i,j,M)=F2d(i,j);
    print_scaf(x,"forallh_2d",1,M,"%+10.3e ");

    setzero_scaf(x,M);
    forall(i,j,M) // loops in all our area, NOT including halos
      X(i,j,M)=F2d(i,j);
    print_scaf(x,"forall_2d",1,M,"%+10.3e ");

    setzero_scaf(x,M);
    foralli(i,j,M) // loops in our area, EXCLUDING BOUNDARY CONDITIONS and HALOS
      X(i,j,M)=F2d(i,j);
    print_scaf(x,"foralli_2d",1,M,"%+10.3e ");

    setzero_scaf(x,M);
    foreast(i,j,M) // loops east boundary (if owned)
      X(i,j,M)=F2d(i,j);
    print_scaf(x,"foreast_2d",1,M,"%+10.3e ");

    setzero_scaf(x,M);
    forwest(i,j,M)
      X(i,j,M)=F2d(i,j);
    print_scaf(x,"forwest_2d",1,M,"%+10.3e ");

    setzero_scaf(x,M);
    fornorth(i,j,M)
      X(i,j,M)=F2d(i,j);
    print_scaf(x,"fornorth_2d",1,M,"%+10.3e ");

    setzero_scaf(x,M);
    forsouth(i,j,M)
      X(i,j,M)=F2d(i,j);
    print_scaf(x,"forsouth_2d",1,M,"%+10.3e ");

    free(x);
}

void test_join_map(map *M) {
  map MJ;
  map *MJ_=&MJ;
  int i,j,nerror=0;
  int info=0;
  double *x,*y;
  x=dmem(M); // field to be gathered

  forall(i,j,M) X(i,j,M)=F2d(i,j); // put some data

  //ifowned(7,8,M) X(7,8,M)=11;
  //ifowned(2,8,M) X(2,8,M)=11;

  y= alloc_gather(x,M,&MJ,0); // only proc 0 allocs the map !!

  printmap("global",&MJ);


  if (quisoc()==0) {
    forall(i,j,MJ_) {
      if (Y(i,j,MJ_) != F2d(i,j) ) {
        nerror++;
        pprintf("error_join %d %d \n",i,j);
      }
      else {
        if (info) pprintf("ok_join %d %d \n",i,j);
      }
    }
    print_scaf(y,"gathered field",0,&MJ,"%+13.5e ");
  }

  free(x);
  if (quisoc()==0) free(y);
  if (nerror>0) CRASH("oops falla alloc_gather nerror=%d",nerror);
}


 // aux function to test halo update with periodic BC
int pertest(int i,int ini,int fi,int per) {
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



int funfunij(int i,int j,map *M) {
  int r;
  r=pertest(i,M->gl0[1],M->gl1[1],M->per[1])+
    pertest(j,M->gl0[2],M->gl1[2],M->per[2])*10;
  return(r);
}




void test_distributed_arrays() { // should run with NN procs
    map M;

    // test2d
    int gl0[4]={0,1,1,0};
    int gl1[4]={0,10,10,0};
    int np[4]={0,2,2,1};
    int per[4]={0,1,1,0};

/* test 2d */
    if (np[1]*np[2]*np[3]!=quants()) CRASH("Please run with %d procs",np[1]*np[2]*np[3]);

    createM(2,gl0,gl1,np,2,per,&M);
    printmap("test_distributed_arrays 2d",&M);

    //test_distributed_arrays_2d(&M);  // check basic functions and loops
    //test_halo_update(&M,0);  // check halo update
    test_join_map(&M); // check join
    //test_get_region(&M);
    //test_operators(&M);
    //pprintf("test 2d finished \n");

if (0) { //  test3d
 // note that map should not be deallocated, all memory is static

    gl0[1]=1; gl0[2]=1;  gl0[3]=1;
    gl1[1]=10; gl1[2]=10;  gl1[3]=10;
    np[1]=3;  np[2]=3;   np[3]=1;
    per[1]=0; per[2]=0;  per[3]=0;
    createM(3,gl0,gl1,np,2,per,&M);

    test_distributed_arrays_3d(&M);  // check basic functions and loops
}
//    test_halo_update_3d(&M);
//    test_operators(&M);
    pprintf("test 3d finished \n");
}


void test_halo_update_3d(map *M) {
  int i,j,k;
  double *x;

  CRASH("Epps atencio que cal verificar test_halo_update_3d !!!");
  x=dmem(M); // allocate the memory for x

  if (M->nd==2) CRASH("sorry has to be 3d");
  setzero_scaf(x,M);
  Forall(i,j,k,M) {
    X3(i,j,k,M)=F3d(pertest(i,M->gl0[1],M->gl1[1],M->per[1]),j,k); // Arbitray (periodic in x) function
  }
  halo_update(x,M);
  print_scaf(x,"test_halo_update_3d",1,M,"%+10.3e ");

  Ifowned(2,2,1,M) X3(2,2,1,M)=33;

  Forallh(i,j,k,M) {
      double val=F3d(pertest(i,M->gl0[1],M->gl1[1],M->per[1]),j,k) ;
      if (  k>=M->gl0[3] && k<=M->gl1[3] &&  // we are in inner nodes in k
            j>=M->gl0[2] && j<=M->gl1[2] &&  // if we are in inner nodes in j
            ( ( i>=M->gl0[1] && i<=M->gl1[1] ) || M->per[1]==1) && // idem in i .. or periodic
            X3(i,j,k,M) != val ) // if the result doesn't agree
        pprintf("error test_halo_update_3d ijk=%d %d %d is=%e should be=%e \n",i,j,k,X3(i,j,k,M),val );
  }

  free(x);
}


/// ************



void benchmark_halo_update_2d(int NX,int NY,int NPX,int NPY) { // should run with NPX*NPY procs
    map M;
    int info=1;

    // test2d
    int gl0[4]={0,1,1,0};
    int gl1[4]={0,NX,NY,0};
    int np[4]={0,NPX,NPY,1};
    int per[4]={0,1,1,0};

/* test 2d */
    if (np[1]*np[2]*np[3]!=quants()) CRASH("Please run with %d procs",np[1]*np[2]*np[3]);

    createM(2,gl0,gl1,np,2,per,&M);
    if (info) printmap("benchmark_halo_update 2d",&M);

    //double r=test_halo_update(&M,1);  // check halo update


    pprintf("benchmark_halo_update\n");
}
