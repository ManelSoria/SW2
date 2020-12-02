// SW2
// 2018-2020
// Manel Soria, Arnau Prat, Arnau Sabates, Marc Andres, Enrique Garcia-Melendo
// UPC - ESEIAAT - TUAREG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#include "sppde.h"

#include "sw.h"

// run as: SW input_file output_folder
// the data generated will be saved in output_folder and named after input_file

int main(int argc, char **argv) {
    char output_folder[MAXOUTLEN];
    char bname[MAXOUTLEN];

    checkr(MPI_Init(&argc,&argv),"init"); // Initialize MPI execution environment

    MPI_Barrier(MPI_COMM_WORLD); // allow all procs to finish Init

//                               0  1   2   3          4
    if (argc!=5) CRASH("SW VERSION=%s COMPILED %s %s\nRun as\nSW npx npy input_file output_folder",VERSION,__DATE__,__TIME__); //MARC

    if (strlen(argv[4])>MAXOUTLEN) CRASH("Output path=<%s> is too long",argv[2]);
    strncpy(output_folder,argv[4],MAXOUTLEN);
    strncpy(bname,argv[3],MAXOUTLEN);


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


    cr_reset(); // Initialize  profiling tools

    pprintf("SW STARTED VERSION: %s\n",VERSION);
    pprintf("COMPILED: %s %s\n",__DATE__,__TIME__);
    pprintf("%d processes started\n",quants());
    cr_start("TOTAL",0);

    sw myworld;
    sw *SW=&myworld;

    int npx,npy;
    sscanf(argv[1],"%d",&npx);
    sscanf(argv[2],"%d",&npy);
    pprintf("npx=%d npy=%d \n",npx,npy);

    initialize(npx,npy,
        0, // don't force 1 cpu
        bname,
        1, // copy input file to destination folder
        SW);
    loop(bname,output_folder,SW);
    destroy(SW);

    cr_end("TOTAL",0);
    pprintf("PROGRAM ENDED\n");

    cr_info(); // Print profiling data
    end_pprintf(); // End print jobs
    MPI_Finalize(); // Terminate MPI environment

    return 0;
}
