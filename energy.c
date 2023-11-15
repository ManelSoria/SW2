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
#include "mpi.h"

#include "sppde.h"
#include "sw.h"

void print_energy(double t,double ke,double ape) {
        pprintf("ENERGY %e %e %e %e\n",t,ke,ape,ke+ape);
}

void calc_energy_balance(double *ke,double *ape,sw *SW){

    int i, j;

    double avegh;
    double surf;
    double myavegh=0;
    double mysurf=0;

    double myke=0,myape=0;
    mysurf = 0;
    double uc, vc;
    forall(i,j,SW->M) {
        uc      = 0.5*(UWINDX(i,j)+UWINDX(i-1,j)); // wind+perturbation
        vc      = 0.5*(VWINDY(i,j)+VWINDY(i,j-1));
        mysurf  = DXC(i,j)*DYC(i,j); // correcte ??? aproximar x un trapezi ?
        myke   += 0.5*H(i,j)*(uc*uc+vc*vc)*mysurf-0.5*(POW2(WINDX(i,j))+POW2(WINDY(i,j)))*mysurf;
        myape  +=  0.5*GC(i,j)*POW2(H(i,j))*mysurf - 0.5*GC(i,j)*POW2(SW->depth)*mysurf;
    }

    double myene[2],ene[2];
    myene[0]=myke;
    myene[1]=myape;
    checkr(MPI_Allreduce(myene, ene, 2, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD),"Energy");
    *ke=ene[0];
    *ape=ene[1];

//    checkr(MPI_Allreduce(&myke, ke, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD),"KE");
//    checkr(MPI_Allreduce(&myape, ape, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD),"APE");

  }
