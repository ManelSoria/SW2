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

void init_coords_spheroid(sw *SW) {

/*

Un-Bizantine simplified version
Initializes the coordinate information of the Shallow World according to an ellipsoidal coordinate system with uniformly distributed latitudes and longitudes.

- sw *SW (input/output) The Shallow World.

*/

    int info=0;


	// angular distances
    double dlon = (SW->lon1-SW->lon0)/((double)SW->nx);
    double dlat = (SW->lat1-SW->lat0)/((double)SW->ny);
    double dlone = dlon, dlonc = dlon, dlonv = dlon; // A macro e.g. DLONE(i,j) is commented out in sw.h but because the angular grid is uniform it is better to fix a value instead of calling functions every time
    double dlatn = dlat, dlatc = dlat, dlatv = dlat;

    double xe, xc, xv;
    double dxe, dxc, dxv;
    double yn, yc, yv;
    double dyn, dyc, dyv;

    double deltax,deltay;
    int i,j;

    yn = 0; yc = 0; yv = 0;
    for (j=asy(SW->M);j<=aey(SW->M);j++) { // loop all the domain, including external halos
        xe = 0; xc = 0; xv = 0;

        // get distances (m), only depend on the latitude (j)
        dxe = dlone*rZ(latc(j));
        dxc = dlonc*rZ(latc(j));
        dxv = dlonv*rZ(latn(j));
        dyn = dlatn*rM(latc(j));
        dyc = dlatc*rM(latn(j-1));
        dyv = dlatv*rM(latc(j));

        // accumulate vertical distances
        yn += dyn;
        yc += dyc;
        yv += dyv;

        for (i=asx(SW->M);i<=aex(SW->M);i++) { // loop all the domain, including external halos
        	// accumulate horizontal distances
            xe += dxe;
            xc += dxc;
            xv += dxv;

            if(i==0 && j==0 ){ // save coordinates of the starting point
                deltax= xe;
                deltay= yn;
            }


			// when we reach our area (including halos), we save the data
      // we don't do a halo update because in the periodic direction coordinates are messed up
            ifownedh(i,j,SW->M) {
                DXE(i,j) = dxe;
                DXC(i,j) = dxc;
                DXV(i,j) = dxv;
                DYN(i,j) = dyn;
                DYC(i,j) = dyc;
                DYV(i,j) = dyv;

                XE(i,j) = xe;
                XC(i,j) = xc;
                XV(i,j) = xv;
                YN(i,j) = yn;
                YC(i,j) = yc;
                YV(i,j) = yv;
            }
        }
    }

    // substract the values so that zero position is where it should
    forallh(i,j,SW->M) {
        XE(i,j) -= deltax;
        XC(i,j) -= deltax;
        XV(i,j) -= deltax;
        YN(i,j) -= deltay;
        YC(i,j) -= deltay;
        YV(i,j) -= deltay;
    }

    if(info>1) {
        print_scaf(SW->Dxe, "Dxe", 1, SW->M,"%8.2e ");
        print_scaf(SW->Dxc, "Dxc", 1, SW->M,"%8.2e ");
        print_scaf(SW->Dxv, "Dxv", 1, SW->M,"%8.2e ");
        print_scaf(SW->Dyn, "Dyn", 1, SW->M,"%8.2e ");
        print_scaf(SW->Dyc, "Dyc", 1, SW->M,"%8.2e ");
        print_scaf(SW->Dyv, "Dyv", 1, SW->M,"%8.2e ");
    }
    if(info>1) {
        print_scaf(SW->xe, "xe", 1, SW->M,"%10.6e ");
        print_scaf(SW->xc, "xc", 1, SW->M,"%8.2e ");
        print_scaf(SW->xv, "xv", 1, SW->M,"%8.2e ");
        print_scaf(SW->yn, "yn", 1, SW->M,"%8.2e ");
        print_scaf(SW->yc, "yc", 1, SW->M,"%8.2e ");
        print_scaf(SW->yv, "yv", 1, SW->M,"%8.2e ");
    }

    //CRASH("poc a poc");
}


//test_init_coords_spheroid(sw *SW)
