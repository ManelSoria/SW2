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
#include <string.h>
#include <math.h>
#include <libgen.h>
#include "mpi.h"

#include "sppde.h"
#include "sw.h"
#include "ensight.h"

// POSTPROCESSING
#define UC(i,j) ((U(i-1,j) + U(i,j)) / 2.) // The horizontal velocity at a centered point
#define VC(i,j) ((ac(SW->v,i,j-1,SW->M) + ac(SW->v,i,j,SW->M)) / 2.) // The vertical velocity at a centered point
#define UWINDXC(i,j) ((UWINDX(i-1,j) + UWINDX(i,j)) / 2.)
#define DVDXV(i,j) ((VWINDY(i+1,j)-VWINDY(i,j))/DXV(i,j))
#define DUWINDXDYV(i,j) ((UWINDX(i,j+1)-UWINDX(i,j))/DYC(i,j+1))
#define VORT(i,j) (DVDXV(i,j) - DUWINDXDYV(i,j)) // The vorticity at a vertex point
#define PVORT(i,j) ((VORT(i,j)+FV(i,j))/(0.25*(H(i,j)+H(i+1,j)+H(i,j+1)+H(i+1,j+1)))) // The potential vorticity at a vertex point


void save_ens(int outCoordsMode,char *bname, int Nstart, int Nend, char *binaries_folder, char *output_folder, sw *SW);


int main(int argc, char **argv) {
    int info=0;

    char binaries_folder[MAXOUTLEN]; // Where binary files are stored
    char output_folder[MAXOUTLEN]; // Where the EnSight Gold and stdout files will be stored
    char bname[MAXOUTLEN]; // Case name

    checkr(MPI_Init(&argc,&argv),"init"); // Initialize MPI execution environment

    if (quants()!=1) CRASH("Please run postproc with a single processor");

    MPI_Barrier(MPI_COMM_WORLD); // allow all procs to finish Init

    if (argc!=7) CRASH("Postproc VERSION=%s COMPILED %s %s\nRun as:\npostproc (1:lon-lat/2:2d/3:3d/4:text/5:only energy) input_file Nstart Nend binaries_folder output_folder\nUse Nend=-1 to process all files\n",VERSION,__DATE__,__TIME__);

    if (strlen(argv[6])>MAXOUTLEN) CRASH("Output path=<%s> is too long",argv[6]);
    if (strlen(argv[5])>MAXOUTLEN) CRASH("Output path=<%s> is too long",argv[5]);
    strncpy(output_folder,argv[6],MAXOUTLEN);
    strncpy(binaries_folder,argv[5],MAXOUTLEN);
    int Nend = atoi(argv[4]); // Start case number
    int Nstart = atoi(argv[3]); // End case number
    strncpy(bname,argv[2],MAXOUTLEN);
    int outCoordsMode = atoi(argv[1]); // Output coordinate mode 1:lon/lat/2:2d/3:3d /4:save text files

    if (outCoordsMode<0 || outCoordsMode>5) CRASH("wrong outCoordsMode=%d",outCoordsMode);

    if (Nstart<0) CRASH("oops Nstart=%d should be >=0 ",Nstart);
    if (Nend<-1 ) CRASH("oops Nend=%d shoud be >=-1",Nend);
    // mode:
    // 0: all processors to file
    // 1: processor 0 to stdout, the rest to the file
    // 2: processor 0 to stdout, the rest ignore output
    // set flushall to 0 for production runs !!!
    int mode=2;
    #ifdef CHECK
    int flushall=1;
    mode=0;
    #else
    int flushall=0;
    #endif

    init_pprintf(mode,flushall,output_folder);


    if (info) {
        pprintf("Nstart=%d Nend=%d \n",Nstart,Nend);
    }

    cr_reset(); // Initialize  profiling tools

    pprintf("postproc STARTED VERSION: %s\n",VERSION);
    pprintf("%d processes started\n",quants());
    cr_start("TOTAL",0);

    sw myworld;
    sw *SW=&myworld;

    initialize(1,1,1,bname,0,SW); // npx=1,npy=1, force 1..; DONT COPY INPUT FILE TO DESTINATION
    save_ens(outCoordsMode,bname,Nstart,Nend,binaries_folder,output_folder,SW);
    destroy(SW);

    cr_end("TOTAL",0);
    pprintf("PROGRAM ENDED\n");

    cr_info(); // Print profiling data
    end_pprintf(); // End print jobs
    MPI_Finalize(); // Terminate MPI environment

    return 0;
}

void save_ens(int outCoordsMode,char *bname, int Nstart, int Nend, char *binaries_folder, char *output_folder, sw *SW) {

    int info=0;

    char *var_names[N_OUT_VARIABLES] = { // Depends on N_OUT_VARIABLES
        "ETA",
        "U",
        "V",
        "WINDX",
        "TRACER",
        "PVORT",
        "PERTS"
    };
    double *pvort;
    pvort = dmem(SW->M);
    double *var_fields[N_OUT_VARIABLES] = { // Depends on N_OUT_VARIABLES
        SW->eta,
        SW->u,
        SW->v,
        SW->windx,
        SW->tracer,
        pvort,
        SW->perturbs
    };
    double var_gridpos[N_OUT_VARIABLES][2] = { // Depends on N_OUT_VARIABLES
        /*
        {cx, cy} if not centered it is assumed as staggered (cx=1 => sx=0, cy=1 => sy=0)
        {0, 0} for true value instead of linearly interpolated to nodal/vertex position
        */
        {1, 1}, // SW->eta
        {0, 1}, // SW->u
        {1, 0}, // SW->v
        {0, 1}, // SW->zwx
        {1, 1}, // SW->tracer
        {0, 0}, // pvort
        {1, 1}, // perturbs
    };


    double t;
    double *timevals;
    // Memory for the time values
    timevals = (double *)mallocc(sizeof(double)*(1));

    char outname[MAXOUTLEN];
    int i,j,k;
    int nstep; // dummy variable
    int load_error;
    int q=0; // position of our frame in the timevals vector
    k=Nstart;
    while (1) {
        if (info) pprintf("k=%d q=%d \n",k,q);
        load_error = load_binary_sw(k,&nstep,&t,bname,binaries_folder,SW); // load_binary_sw updates SW->loadFrom
        if(load_error==1){
            pprintf("file %d not found, stopping post-process\n",k);
            break;
        }

        do_channel(SW);

        // energy::

        forallh(i,j,SW->M) {
            // if((PVORT(i,j))<(-2.4e-07)) {
                // pprintf(
                // "PVORT=%e, FV=%e, UC=%e, VC=%e, UWINDXC=%e, DVDXV=%e, DUWINDXDYV=%e, VORT=%e\n",
                // PVORT(i,j), FV(i,j), UC(i,j), VC(i,j), UWINDXC(i,j), DVDXV(i,j), DUWINDXDYV(i,j), VORT(i,j)
                // );
            // }
            ac(pvort,i,j,SW->M) = PVORT(i,j); // PVORT(i,j) does not access a map SW->M
        }

        halo_update(pvort,SW->M);
        forsouth(i,j,SW->M){
          ac(pvort,i,j-1,SW->M) = PVORT(i,j);
        }
        ifownedh(0,0,SW->M){
          ac(pvort,0,0,SW->M) = PVORT(1,1);
        }

        double ke,ape;
        calc_energy_balance(&ke,&ape,SW);
        print_energy(t,ke,ape);

        if (outCoordsMode!=5) {
            if (outCoordsMode==4) { // special case 4, just write text files
                snprintf(outname,MAXOUTLEN,"%s/%s.%04d.ETA",output_folder,basename(bname),q);
                text_write_scaf(SW->eta,outname,SW->M);
            }  else {
                snprintf(outname,MAXOUTLEN,"%s/%s",output_folder,basename(bname));
                save_ens_step(outname, var_names, var_fields, var_gridpos, k,SW);
            }

            pprintf("postproc generated nstep=%d t=%f as frame %d\n",nstep,t,q);
        }

        timevals=(double *)realloc(timevals,sizeof(double)*(q+1));
        if (timevals==NULL) CRASH("Realloc failed.");
        timevals[q]=t;
        q++; // now q contains the number of frames
        if (k==Nend) break; // stop if reached last frame
        k++;
    }


    free(pvort);

    /*
    The following cannot be put at the top and before the loop.
    */
    if (outCoordsMode!=4 && outCoordsMode!=5) { // for normal ensight files..
        save_ens_geo_sw(outname, outCoordsMode, SW);
        save_ens_case(outname, var_names, Nstart, q, timevals);
    }

    free(timevals);
}
