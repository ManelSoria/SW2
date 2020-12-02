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

// run as:

int main(int argc, char **argv) {
    int info=0;

    char binaries_folder[MAXOUTLEN]; // Where binary files are stored
    char output_folder[MAXOUTLEN]; // Where the EnSight Gold and stdout files will be stored
    char bname[MAXOUTLEN]; // Case name

    checkr(MPI_Init(&argc,&argv),"init"); // Initialize MPI execution environment

    if (quants()!=1) CRASH("Please run casegen with a single processor");

    MPI_Barrier(MPI_COMM_WORLD); // allow all procs to finish Init

    if (argc!=6) CRASH("Run as:\ncasegen input_file Nstart Nend binaries_folder output_folder\n");

    if (strlen(argv[1])>MAXOUTLEN) CRASH("<%s> is too long",argv[1]);
    if (strlen(argv[2])>MAXOUTLEN) CRASH("<%s> is too long",argv[2]);
    if (strlen(argv[3])>MAXOUTLEN) CRASH("<%s> is too long",argv[3]);

    strncpy(output_folder,argv[5],MAXOUTLEN);
    strncpy(binaries_folder,argv[4],MAXOUTLEN);
    int Nend = atoi(argv[3]); // Start case number
    int Nstart = atoi(argv[2]); // End case number
    strncpy(bname,argv[1],MAXOUTLEN);

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

    char *var_names[N_OUT_VARIABLES] = { // Depends on N_OUT_VARIABLES
        "ETA",
        "U",
        "V",
        "WINDX",
        "TRACER",
        "PVORT",
        "PERTURBS"
    };

    double t;
    double *timevals=(double *)malloc(sizeof(double)*1);

    char outname[MAXOUTLEN];
    int i,j,k;
    int nstep; // dummy variable
    int load_error;
    int q=0; // position of our frame in the timevals vector
    k=Nstart;
    while (1) {
        if (info) pprintf("k=%d q=%d \n",k,q);
        load_error = load_binary_header_sw(k,&nstep,&t,bname,binaries_folder);
        if(load_error==1){
            pprintf("file %d not found, stopping post-process\n",k);
            break;
        }
        pprintf("casegen generated nstep=%d t=%f as frame %d\n",nstep,t,q);
        timevals=(double *)realloc(timevals,sizeof(double)*(q+1));
        if (timevals==NULL) CRASH("Realloc failed.");
        timevals[q]=t;
        q++;
        if (k==Nend) break; // stop if reached last frame
        k++;
    }

    snprintf(outname,MAXOUTLEN,"%s/%s",output_folder,basename(bname));
    save_ens_case(outname, var_names, Nstart, q, timevals);

    end_pprintf(); // End print jobs
    MPI_Finalize(); // Terminate MPI environment

    return 0;
}
