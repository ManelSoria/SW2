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
#include "mpi.h"

#include "sppde.h"

#include "sppde_tests.h" // prototypes of test functions

#include "sppde_parser.h"

#include "sw.h"

// test in normal mode then in TEST the TEST mode
#define CTEST(fun) pprintf("TESTING %s ... ",#fun); \
   oneok=fun(M,0);  allok&=oneok; pprintf("%d  ",oneok); \
   pprintf("TtT.."); oneok=fun(M,1); pprintf("%d\n",oneok); allok&=!oneok;

int calltests(map *M) {

    int oneok,allok=1;

    // test_sweeps2_2d(M,0);
    test_gs_field_2d(M,0);
    // call them from the most simple to the most complex.
    // eg. binaryfiles uses halo_update, gather ... so it has to be about the last
    if (0) {
    CTEST(test_sweeps_2d);
    CTEST(test_gs_field_2d);
    CTEST(test_halo_update_2d);
    CTEST(test_binaryfiles_2d);
    }

    return(allok);
}

int test_all(int SX,int SY,int perx,int pery,int npx,int npy,int npz) {
    int ok;
    map M;
    // Cal posar mapes diferents en els dos eixos i amb nombres primers pq surtin combinacions rares
    // Tambe un mapa no distribuit... donara problemes pero els anirem trobant

    int gl0[4]={0,1,1,0}; // global limits (start)
    int gl1[4]={0,SX,SY,0}; // global limits (end)
    int np[4]={0,npx,npy,npz}; // processor grid
    //int per[4]={0,1,0,0}; // periodicity (only in x)
    int per[4]={0,perx,pery,0}; //  periodicity


    if (np[1]*np[2]*np[3]!=quants()) CRASH("Please run with %d procs",np[1]*np[2]*np[3]);

    createM(2,gl0,gl1,np,2,per,&M);

    ok=calltests(&M);

    return(ok);

}

void pbdhalo(double *t,map *M) {
    int i,j;
    int info=1;
    if (info) pprintf("obN=%d obS=%d obE=%d obW=%d \n",obN(M),obS(M),obE(M),obW(M));
    if (obN(M)) {
        for (i=lsx(M)-HS(M);i<=lex(M)+HS(M);i++)
            for (j=ley(M)+1;j<=ley(M)+HS(M);j++) {
                if (info) pprintf("pdbhaloN %d %d \n",i,j);
                ac(t,i,j,M)=ac(t,i,j-1,M);
            }
    }
    if (obS(M)) {
        for (i=lsx(M)-HS(M);i<=lex(M)+HS(M);i++)
            for (j=lsy(M)-1;j>=lsy(M)-HS(M);j--) {
                if (info) pprintf("pdbhaloS %d %d \n",i,j);
                ac(t,i,j,M)=ac(t,i,j+1,M);
            }
    }

}

void test_easy() { // MANEL
    int SX=10,SY=10;
    int npx=3,npy=3,npz=1;
    int perx=1,pery=0;
    map M;
    // Cal posar mapes diferents en els dos eixos i amb nombres primers pq surtin combinacions rares
    // Tambe un mapa no distribuit... donara problemes pero els anirem trobant

    int gl0[4]={0,1,1,0}; // global limits (start)
    int gl1[4]={0,SX,SY,0}; // global limits (end)
    int np[4]={0,npx,npy,npz}; // processor grid
    //int per[4]={0,1,0,0}; // periodicity (only in x)
    int per[4]={0,perx,pery,0}; //  periodicity


    if (np[1]*np[2]*np[3]!=quants()) CRASH("Please run with %d procs",np[1]*np[2]*np[3]);

    createM(2,gl0,gl1,np,2,per,&M);

    double *t=dmem(&M);

    set_field(t,2,&M);
    halo_update(t,&M);
    printmap("test",&M);
    print_scaf(t,"before",1,&M,"%6.2f ");

    pbdhalo(t,&M);

    print_scaf(t,"after",1,&M,"%6.2f ");

    free(t);

    //test_gs_field_2d_eh(1,&M,0);
 //   int ok=test_binaryfiles_2d(&M,0);
    int ok=1;
    if (ok==0) CRASH("uhhh ? test failed!!");
}

int main(int argc, char **argv) { // test
    checkr(MPI_Init(&argc,&argv),"init"); // MPI should be started before calling p_startinput


    init_pprintf(0,1,"/tmp/");

    pprintf("demopetita STARTED VERSION: %s\n",VERSION); // MANEL
    pprintf("COMPILED: %s %s\n",__DATE__,__TIME__); // MANEL

    pprintf("Hola\n");



//    FILE *fin;

//    fin=p_startinput("testparser.txt",1); // MANEL
//    p_setdebug(1); // turn on debug, extra level
 //   double radius=p_getdouble(fin,"planet_radius");
 //   fclose(fin);

    test_easy();
    end_pprintf(); // End print jobs

    MPI_Finalize();

    return(0);
}

// main de tests..
// main(..)
        // create_map(&M) creem un map d'una certa mida, especificada desde comandes
        //
        // pprintf("test_something = %d \n",test_something(M)); 1: ok; 0: error
        // pprintf("test_other     = %d \n",test_other(M)); 1: ok; 0: error
        // ...
//        test_cumsum();

int main2(int argc, char **argv) {
    int ok;

    checkr(MPI_Init(&argc,&argv),"init"); // Initialize MPI execution environment
    init_pprintf(1,1,"/tmp/");

    int SX,SY,perx,pery,npx,npy,npz;

    if (argc!=8) {
        CRASH("test SX SY perx pery npx npy npz");
    }

    sscanf(argv[1],"%d",&SX);
    sscanf(argv[2],"%d",&SY);
    sscanf(argv[3],"%d",&perx);
    sscanf(argv[4],"%d",&pery);
    sscanf(argv[5],"%d",&npx);
    sscanf(argv[6],"%d",&npy);
    sscanf(argv[7],"%d",&npz);

    pprintf("SX=%d SY=%d perx=%d pery=%d npx=%d npy=%d npz=%d \n",SX,SY,perx,pery,npx,npy,npz);

    if (npx*npy*npz!=quants()) CRASH("number of division and number of processors must be equal");
    if (npz!=1) CRASH("sorry npz=%d and should be 1",npz);

    ok=test_all(SX,SY,perx,pery,npx,npy,npz);
    if (ok==1)
        pprintf("SPPDE: ALL SEEMS OK !! \n");
    else
        pprintf("Ooops there is at least one error !! \n");
    end_pprintf(); // End print jobs
    MPI_Finalize(); // Terminate MPI environment

    return !ok;
}



#define XX(i,j) ac(x,i,j,&M)

int main3(int argc, char **argv) {
    checkr(MPI_Init(&argc,&argv),"init"); // Initialize MPI execution environment
    init_pprintf(0,1,"/tmp/");
    pprintf("hola quisoc=%d quants=%d \n",quisoc(),quants());
    int info=0;

    map M;
    //         procs  xlim  ylim  hs  per
    createM_2d(2,2,   1,13, 1,8,  2,  0,0,&M);

    printmap("prova",&M);

    pprintf("HS=%d\n",HS(&M));

    double *x,*gx;

    x=dmem(&M);

    int i,j,i0,j0;

    i0=2; j0=3;

    forallh(i,j,&M) XX(i,j)=0; // init to 0

    ifowned(i0,j0,&M) XX(i0,j0)=1;

    for (int ite=1;ite<=100;ite++) {
        halo_update(x,&M);

        foralli(i,j,&M) XX(i,j)=(1./4.)*(XX(i+1,j)+XX(i-1,j)+XX(i,j+1)+XX(i,j-1));

        fornorth(i,j,&M) {
            if (info) pprintf("north i=%d j=%d \n",i,j);
            XX(i,j)=(1./3.)*(XX(i+1,j)+XX(i-1,j)+XX(i,j-1));
        }

        forsouth(i,j,&M) {
            if (info) pprintf("south i=%d j=%d \n",i,j);
            XX(i,j)=(1./3.)*(XX(i+1,j)+XX(i-1,j)+XX(i,j+1));
        }

        foreast(i,j,&M) {
            if (info) pprintf("east i=%d j=%d \n",i,j);
            XX(i,j)=(1./3.)*(XX(i-1,j)+XX(i,j-1)+XX(i,j+1));
        }

        forwest(i,j,&M) {
            if (info) pprintf("west i=%d j=%d \n",i,j);
            XX(i,j)=(1./3.)*(XX(i+1,j)+XX(i,j-1)+XX(i,j+1));
        }

        forne(i,j,&M) {
            if (info) pprintf("ne i=%d j=%d \n",i,j);
            XX(i,j)=(1./2.)*(XX(i-1,j)+XX(i,j-1));
        }

        fornw(i,j,&M) {
            if (info) pprintf("nw i=%d j=%d \n",i,j);
            XX(i,j)=(1./2.)*(XX(i+1,j)+XX(i,j-1));
        }

        forsw(i,j,&M) {
            if (info) pprintf("sw i=%d j=%d \n",i,j);
            XX(i,j)=(1./2.)*(XX(i+1,j)+XX(i,j+1));
        }

        forse(i,j,&M) {
            if (info) pprintf("se i=%d j=%d \n",i,j);
            XX(i,j)=(1./2.)*(XX(i-1,j)+XX(i,j+1));
        }

        pprintf("ite=%d max abs= %e \n",ite,norm_scaf(x,0,&M));
    }

    print_scaf(x,"x",1,&M,"%+10.3e ");

    map MJ;

    gx=alloc_gather(x,&M,&MJ,0);

    if (quisoc()==0)
        print_scaf(gx,"gx",0,&MJ,"%+10.3e ");

    free(x); free(gx);


    end_pprintf(); // End print jobs
    MPI_Finalize(); // Terminate MPI environment
}
