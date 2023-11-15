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
#include <string.h>
#include "mpi.h"

#include "sppde.h"
#include "sppde_parser.h"
#include "sw.h"

int initialize(int npx,int npy,int force1proc,char *bname, int copyinput,sw *SW) {

/*

Initializes the Shallow World by reading data from a file.

- sw *SW (input/output) The Shallow World.

*/

    int info = 0;

    int i,j;

    FILE *fin;

    char tmpstr[MAXOUTLEN];
    snprintf(tmpstr,MAXOUTLEN,"%s.sw",bname);

    if (info) pprintf("initialize calling p_startinput with <%s> npx=%d npy=%d\n",tmpstr,npx,npy);

    fin = p_startinput(tmpstr,copyinput);

    p_setdebug(1);

    int sch = p_getint(fin,"sch");
    SW->sch = sch;

    SW->polar = p_getint(fin,"polar");
        if (SW->polar!=0 && SW->polar!=1 && SW->polar!=-1) CRASH("polar=%d and should be 0, 1 (north pole) or -1 (south pole)",SW->polar);

        if (SW->polar==0) { // "normal" geometries
            // Domain dimensions
            SW->lon0 = deg2rad(p_getdouble(fin,"lon0"));
            SW->lon1 = deg2rad(p_getdouble(fin,"lon1"));
            SW->lat0 = deg2rad(p_getdouble(fin,"lat0"));
            SW->lat1 = deg2rad(p_getdouble(fin,"lat1"));


            // Domain cells
            SW->nx = p_getint(fin,"nx");
            SW->ny = p_getint(fin,"ny");
            SW->nz = 1; // p_getint(fin,"nz");
        } else { // polar geometries
            double degSim;
            degSim = deg2rad(p_getdouble(fin,"angleSim")); // ie, 10 would be from 80 degrees to the pole
            SW->lon0 = -degSim;
            SW->lon1 = degSim;
            SW->lat0 = -degSim;
            SW->lat1 = degSim;


            // Domain cells
            SW->nx = p_getint(fin,"n");
            SW->ny = SW->nx;
            SW->nz = 1; // p_getint(fin,"nz");

        }

    int npz=1;

    pprintf("initialize npx= %d npy= %d kN= %.1f kN/P= %.1f\n",npx,npy,1e-3*(SW->nx)*(SW->ny),1e-3*(SW->nx)*(SW->ny) / (npx*npy) );

    if (force1proc==1) {npx=1; npy=1; npz=1;}

    if (npx*npy*npz!=quants()) CRASH("oops quants=%d doesnt agree with npx npy npz = %d %d %d. Run with %d procs",quants(),npx,npy,npz,npx*npy*npz);

    int y_per = 1;
    if(SW->polar == 0) y_per = 0;


    int gl0[4] = {0,1,1,0};                // Start index
    int gl1[4] = {0,SW->nx,SW->ny,SW->nz}; // End index
    int np[4]  = {0,npx,npy,npz};          // Number of procs per axis
    int per[4] = {0,IS_PERIODIC_X,y_per,0};

    SW->M = &(SW->m);
    createM(2,gl0,gl1,np,2,per,SW->M);     // Initialize struct map

    pprintf("Reading/1\n");
    // Time params
    SW->t0 = p_getdouble(fin,"t0");
    SW->t1 = p_getdouble(fin,"t1");
    SW->Dt = p_getdouble(fin,"Dt");
    SW->Courant_max = p_getdouble(fin,"CFL");

    SW->N = (int)(((SW->t1)-(SW->t0))/(SW->Dt)); // Check and /or change
    double tau = p_getdouble(fin,"Tau");
    if (tau>1e10)
            SW->sigma_dis = 0;
        else
            SW->sigma_dis = 1.0/tau;

    //Hyperviscosity
    double auxnu[3];
    p_getndoubles(fin, 3, auxnu,"Hypernu");
    SW->nu2 =  auxnu[0];
    SW->nu4 =  -auxnu[1];//-1e3; // should be entered as positive !!
    SW->nu6 =  auxnu[2];//1e2;
    if (SW->nu2<0 || SW->nu4>0 || SW->nu6<0) CRASH("Hyperviscosity coefs should be entered as >0 and they are %e %e %e",SW->nu2,-SW->nu4,SW->nu6);

    SW->IteInfo=p_getint(fin,"IteInfo");
    SW->saveEvery=p_getint(fin,"SaveEvery");
    SW->loadFrom=p_getint(fin,"LoadFrom");

    // reference files
    SW->reftype=p_getint(fin,"reftype");
    char *buf;
    switch (SW->reftype) { // REMEMBER TO FREE MEMORY IN finish.c !!
        case 0:
            SW->reftable=(double *)NULL;
            SW->refname[0]='\0';
        break;
        case 1: // compare with text files
            buf=p_getstring(fin,"refname");
            strcpy(SW->refname,buf);
            if (info) pprintf("refname=<%s>\n",SW->refname);
            SW->reftable=p_alloc_gettable(fin,&(SW->reftable_nr),&(SW->reftable_nc),0);
            if (info) pprintf("reftable nr=%d nc=%d\n",SW->reftable_nr,SW->reftable_nc);
        break;
        default: CRASH("wrong reftype=%d",SW->reftype);

    }
    // Physical variables
    SW->Omega = p_getdouble(fin,"Omega");
    SW->rE = p_getdouble(fin,"rE");
    SW->rP = p_getdouble(fin,"rP");
    SW->epsilon2 = POW2((SW->rE)/(SW->rP));
    SW->gE = p_getdouble(fin,"gE");
    SW->depth = p_getdouble(fin,"depth");

    if (SW->polar != 0){
    SW->beta_sponge=deg2rad(p_getdouble(fin,"BetaSp"));
    SW->sigma_sponge=p_getdouble(fin,"SigSp");
    SW->tau_sponge=-log(0.01)/(sqrt(2.0)*SW->lat1-SW->beta_sponge); // (99% of the sigma at the end of the sponge)
    } else {
    SW->beta_sponge=1e5;
    SW->sigma_sponge=0;
    SW->tau_sponge=0;
    }

    // CV dimesions
    SW->Dxe = dmem(SW->M);
    SW->Dxc = dmem(SW->M);
    SW->Dxv = dmem(SW->M);
    SW->Dyn = dmem(SW->M);
    SW->Dyc = dmem(SW->M);
    SW->Dyv = dmem(SW->M);
    SW->xe = dmem(SW->M);
    SW->xc = dmem(SW->M);
    SW->xv = dmem(SW->M);
    SW->yn = dmem(SW->M);
    SW->yc = dmem(SW->M);
    SW->yv = dmem(SW->M);

    // Gravity and Coriolis parameter matrices
    SW->ge = dmem(SW->M);
    SW->gn = dmem(SW->M);
    SW->gc = dmem(SW->M);
    SW->gv = dmem(SW->M);
    SW->fe = dmem(SW->M);
    SW->fn = dmem(SW->M);
    SW->fc = dmem(SW->M);
    SW->fv = dmem(SW->M);


    if(SW->polar == 0){
      init_channel_region(SW);
    }else if(SW->polar == 1 || SW->polar == -1){
      init_polar_region(SW);
    }else{
      CRASH("Error with SW->polar, it's not 0, 1 or -1");
    }

    init_coords_spheroid(SW);

    SW->tracer = dmem(SW->M);
    SW->perturbs = dmem(SW->M);
    SW->u = dmem(SW->M);
    SW->v = dmem(SW->M);
    SW->windx = dmem(SW->M);
    SW->windy = dmem(SW->M);
    SW->eta = dmem(SW->M);
    SW->hB = dmem(SW->M);

    setzero_scaf(SW->tracer, SW->M);
    setzero_scaf(SW->perturbs, SW->M);
    setzero_scaf(SW->u, SW->M);
    setzero_scaf(SW->v, SW->M);
    setzero_scaf(SW->eta, SW->M);
    setzero_scaf(SW->hB, SW->M);

    SW->geoeq=p_getint(fin,"geoeq");
    if (SW->geoeq!=0 && SW->geoeq!=1) CRASH("geoeq=%d and should be 0 or 1",SW->geoeq);

    pprintf("Reading/2 (Vortices)\n");
    read_vortices(fin,SW);

    pprintf("Reading/3 (Gaussians)\n");
    read_Gaussians(fin,SW);

    SW->tracer_dis = p_getdouble(fin,"TracerDis");

    pprintf("Reading/4 (wind). Echo disabled for this part.\n");
    p_stopdebug();

    read_wind(fin,SW);

    // setzero_scaf(SW->windx, SW->M);
    if(SW->polar == 0)halo_update(SW->windx,SW->M);

    //setzero_scaf(SW->windy, SW->M);
    if(SW->polar == 0) halo_update(SW->windy,SW->M);

    p_endinput(fin);

    // Other params
    return 0;
}

void init_polar_region(sw *SW){
  /*
    This function fills the gravity and Coriolis matrices either for north
    and south poles
  */
  int i,j;
  int info = 0;
  double ns = ((double)(SW->polar));

  forallh(i,j,SW->M){
    GE(i,j) = g_eff(lat2pclat(POLAR_LATE(i,j)));
    GN(i,j) = g_eff(lat2pclat(POLAR_LATN(i,j)));
    GC(i,j) = g_eff(lat2pclat(POLAR_LATC(i,j)));
    GV(i,j) = g_eff(lat2pclat(POLAR_LATV(i,j)));

    FE(i,j) = f(POLAR_LATE(i,j));
    FN(i,j) = f(POLAR_LATN(i,j));
    FC(i,j) = f(POLAR_LATC(i,j));
    FV(i,j) = f(POLAR_LATV(i,j));
  }

  if(info){
    forall(i,j,SW->M){
      if (j == round(SW->ny/2.0)){
        printf("CORIOLIS %8.4f %8.6e \n", rad2deg(POLAR_LATC(i,j)),FC(i,j));
        printf("GRAVITY %8.4f %8.6e \n",  rad2deg(POLAR_LATC(i,j)),GC(i,j));
      }
    }
  }

  double new_radius = (SW->rE)*(SW->rE)/(SW->rP);

  SW->rE = new_radius;
  SW->rP = new_radius;
}

void init_channel_region(sw *SW){
  /*
    This function computes and fills the scalar fields of gravity and Coriolis parameter
  */

  int i,j;

  forallh(i,j,SW->M){
    GE(i,j) = g_eff(lat2pclat(LATE(i,j)));
    GN(i,j) = g_eff(lat2pclat(LATN(i,j)));
    GC(i,j) = g_eff(lat2pclat(LATC(i,j)));
    GV(i,j) = g_eff(lat2pclat(LATV(i,j)));

    FE(i,j) = f(LATE(i,j));
    FN(i,j) = f(LATN(i,j));
    FC(i,j) = f(LATC(i,j));
    FV(i,j) = f(LATV(i,j));
  }
}
