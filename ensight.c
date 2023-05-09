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
#include <libgen.h>
#include <math.h>
#include "mpi.h"

#include "sppde.h"
#include "sppde_parser.h"
#include "sw.h"
#include "ensight.h"

void save_ens_case(char *bname, char **var_names, int nstart, int nsteps, double *timevals) {

/*

Create EnSight Gold Case file with all variables in var_names.

- char *bname (input) The name of the file.
- char **var_names (input) The name of the variables in the output file.
- double nsteps (input) Number of steps.
- double *timevals (input) The time values.

*/

    map MJ; // distributed map
    map *MJ_=&MJ;
    FILE *fout;
    char nn[MAXOUTLEN]; // complete file name

    int i,j,k;

    if(quisoc()==0) {
        sprintf(nn,"%s.ensi.case",bname);
        fout=fopen(nn,"w");
        if (fout==NULL) CRASH("file cannot be opened <%s>",nn);
    } else
        fout= (FILE *)NULL; // only proc 0 writes

    if(quisoc()==0) {
        fprintf(fout, "FORMAT\n");
        fprintf(fout, "type:\tensight gold\n");
        fprintf(fout, "GEOMETRY\n");
        fprintf(fout, "model:%10d\t%s.ensi.geo\n",1,basename(bname));
        fprintf(fout, "VARIABLES\n");
        k = 0;
        for(k; k < N_OUT_VARIABLES; k++) {
            fprintf(fout, "scalar per node:%10d\t%s\t%s.ensi.%s-****\n",1,var_names[k],basename(bname),var_names[k]);
        }
        fprintf(fout, "TIME\n");
        fprintf(fout, "time set:\t%d\n",1);
        fprintf(fout, "number of steps:\t%d\n",nsteps);
        fprintf(fout, "filename start number:\t%d\n",nstart);
        fprintf(fout, "filename increment:\t%d\n",1);
        fprintf(fout, "time values:\n");
        k = 0;
        for(k ; k < nsteps; k++) {
            fprintf(fout, "%012.5e\n", timevals[k]);
        }
    }

    if(quisoc()==0) {
        fclose(fout);
    }

    checkr(MPI_Barrier(MPI_COMM_WORLD),"Barrier");

}

void save_ens_geo_sw(char *bname, int mode, sw *SW) {

/*

Create ShallowWorlds-specific EnSight Gold Geometry file. Calls generic save_ens_geo(...).

- char *bname (input) The name of the file.
- int mode (input) See below.
- sw *SW (input) The Shallow World.

*/

    int info = 1;

    int i,j;

    double *coords[4];

    switch(mode) {
        case 1:
            coords[1] = dmem(SW->M); forallh(i,j,SW->M) ac(coords[1],i,j,SW->M) =rad2deg(LONV(i,j));
            coords[2] = dmem(SW->M); forallh(i,j,SW->M) ac(coords[2],i,j,SW->M) =rad2deg(LATV(i,j));
            coords[3] = dmem(SW->M); setzero_scaf(coords[3], SW->M);
         break;

        case 2: // 2d (Z=0)
            coords[1] = (SW -> xv);
            coords[2] = (SW -> yv);
            coords[3] = dmem(SW->M); setzero_scaf(coords[3], SW->M); // Set z coordinate to 0
        break;

        case 3: // 3d
            coords[1] = dmem(SW->M);
            coords[2] = dmem(SW->M);
            coords[3] = dmem(SW->M);
            forallh(i,j,SW->M) ac(coords[1],i,j,SW->M) = XV3D(i,j);
            forallh(i,j,SW->M) ac(coords[2],i,j,SW->M) = YV3D(i,j);
            forallh(i,j,SW->M) ac(coords[3],i,j,SW->M) = ZV3D(i,j);
        break;
        default: CRASH("unknown mode=%d",mode);
    }


    if(info>1) {
        print_scaf(coords[1], "coords[1]", 1, SW->M,"%8.2e ");
        print_scaf(coords[2], "coords[2]", 1, SW->M,"%8.2e ");
        print_scaf(coords[3], "coords[3]", 1, SW->M,"%8.2e ");
    }

    int nelems[4] = {0, SW->nx, SW->ny, SW->nz}; // First element (current 0) is not used

    save_ens_geo(bname, coords, nelems, SW->M);

    switch(mode) {
       case 1:
            free(coords[1]);
            free(coords[2]);
            free(coords[3]);
        break;

        case 2:
            free(coords[3]);
        break;

        case 3:
            free(coords[1]);
            free(coords[2]);
            free(coords[3]);
        break;
        default: CRASH("unknown mode=%d",mode);
    }


}

void save_ens_geo(char *bname, double **coords, int *nelems, map *M) { //x,y are nodal coordintates (xe, ye) but z is not used in shallow worlds until 3d is implemented

/*

Create EnSight Gold Geometry file (generic function).

- char *bname (input) The name of the file.
- double **coords (input) coords[1] are the x coordinates, coords[2] are the y coordinates,...
- double *nelems(input) nelems[1] is the number of elements in the x axis, nelems[2] is the number of elements in the y axis
- map *M (input) The sppde's map.

*/

    // Should easily be ported to parallell by avoiding allgather and changing MJ_ by M
    // nelemsx is SW->nx...

    map MJ; // distributed map
    map *MJ_=&MJ;
    FILE *fout;
    char nn[MAXOUTLEN]; // complete file name

    double *gscaf; // Scalar field

    double *x = coords[1];
    double *y = coords[2];
    double *z = coords[3];

    int nelemsx = nelems[1];
    int nelemsy = nelems[2];
    int nelemsz = 1; // nelems[3];

    int i,j,k;
    int info=1;

    if(quisoc()==0) {
        sprintf(nn,"%s.ensi.geo",bname);
        fout=fopen(nn,"w");
        if (fout==NULL) CRASH("file cannot be opened <%s>",bname);
    } else
        fout= (FILE *)NULL; // only proc 0 writes

    if(quisoc()==0) {
        // printmap("MJ_",MJ_);
        fprintf(fout, "Problem name: %s\n", basename(bname));
        fprintf(fout, "Geometry file\n");
        fprintf(fout, "node id given\n");
        fprintf(fout, "element id given\n");
        fprintf(fout, "part\n");
        fprintf(fout, "%10d\n", 1);
        fprintf(fout, "Volume Mesh\n");
        fprintf(fout, "coordinates\n");
        fprintf(fout, "%10d\n", (nelemsx+1)*(nelemsy+1));
    }
    joinM(M,MJ_); // GEt global map MJ_ from local M
    if(quisoc()==0) {
        k = 1;
        forallX(i,j,MJ_) {
            fprintf(fout, "%10d\n", k);
            k++;
        }
    }
    // Node x coordinates
    gscaf=alloc_gather(x,M,MJ_,1);
    if(quisoc()==0) {
        forallX(i,j,MJ_) {
            fprintf(fout, "%12.5e\n", ac(gscaf,i,j,MJ_));
        }
        if(info>1) print_scaf(gscaf, "X", 1, MJ_,"%+8.2e ");
        free(gscaf);
    }
    // Node y coordinates
    gscaf=alloc_gather(y,M,MJ_,1);
    if(quisoc()==0) {
        forallX(i,j,MJ_) {
            fprintf(fout, "%12.5e\n", ac(gscaf,i,j,MJ_));
        }
        if(info>1) print_scaf(gscaf, "Y", 1, MJ_,"%+8.2e ");
        free(gscaf);
    }
    // Node z coordinates
    gscaf=alloc_gather(z,M,MJ_,1);
    if(quisoc()==0) {
        forallX(i,j,MJ_) {
            fprintf(fout, "%12.5e\n", ac(gscaf,i,j,MJ_));
        }
        if(info>1) print_scaf(gscaf, "Z", 1, MJ_,"%+8.2e ");
        free(gscaf);
    }


    if(quisoc()==0) {
        fprintf(fout, "quad4\n");
        fprintf(fout, "%10d\n", nelemsx*nelemsy);
        // Element Id
        k = 1;
        forall(i,j,MJ_){
            fprintf(fout, "%10d\n",k);
            k++;
        }
        // Element nodes
        k = 1;
        forall(i,j,MJ_){
            fprintf(fout, "%10d%10d%10d%10d\n", i+(j-1)*(nelemsx+1), i+1+(j-1)*(nelemsx+1), (i+1)+(nelemsx+1)+(j-1)*(nelemsx+1), i+(nelemsx+1)+(j-1)*(nelemsx+1)); // nelsx + 1 is the number of nodes in the x direction. For a nelsx-by-nelsy = 3-by-3, the node numbers for the first element (el 1) are 1 2 6 5 (order is important as filetype assumes that for quad4 numbering is anticlockwise)
            k++;
        }
    }

    if(quisoc()==0)
        fclose(fout);


    checkr(MPI_Barrier(MPI_COMM_WORLD),"Barrier");

}

void save_ens_step(char *bname, char **var_names, double **var_fields, double (*var_gridpos)[2], int frame, sw *SW) {

/*

Generates Ensight Gold Step file for all variables by calling save_ens_var(...).

- char *bname (input) The name of the file.
- char **var_names (input) The name of the variables in the output file.
- double **var_fields (inputs) Reference to the addresses of the fields.
- double (*var_gridpos)[2] Contains info of the position on the grid of each variable. {{1,1}, {1,0}, {0,0}, ..., {cx, cy}}
- int frame (input) The frame number.
- sw *SW (input) The Shallow World.

*/

    int info = 0;
    int k;

    for(k=0; k <= N_OUT_VARIABLES-1; k++) {
        save_ens_var(bname, var_names[k], var_fields[k], frame, var_gridpos[k][0], var_gridpos[k][1], SW->M);
    }


}

void save_ens_var(char *bname, char *scafname, double *scaf, int step, int cx, int cy, map *M) {

/*

Create EnSight Gold Step file for a single variable.

- char *bname (input) The name of the file.
- char *scafname (input) Variable name for EnSight variable file.
- double *scaf (ipnut) The scalar field to be printed in file.
- double step (input) current step.
- int cx (input) Is centered in x? 1 = yes, 0 = no (assumes staggered)
- int cy (input) Is centered in y? 1 = yes, 0 = no (assumes staggered)
- map *M (input) The sppde's map.

*/

    map MJ; // distributed map
    map *MJ_=&MJ;
    FILE *fout;
    char nn[MAXOUTLEN]; // complete file name POSAR UNA MIDA SEGURA

    double *gscaf;
    // print_scaf(scaf, "scaf", 1, M,"%8.2e ");
    gscaf=alloc_gather(scaf,M,MJ_,1);

    int i,j,k;
    int iscrash = 0;

    if(quisoc()==0) {
        sprintf(nn,"%s.ensi.%s-%04d",bname,scafname,step);
        fout=fopen(nn,"w");
        if (fout==NULL) CRASH("file cannot be opened <%s>",bname);
    } else
        fout= (FILE *)NULL; // only proc 0 writes

    if(quisoc()==0) {
        // print_scaf(gscaf, "gscaf", 1, MJ_,"%8.2e ");
        fprintf(fout,"EnSight Gold --- Scalar per-node variables file\n");
        fprintf(fout,"part\n");
        fprintf(fout, "%10d\n", 1);
        fprintf(fout,"coordinates\n");
        if(cx==1&&cy==1) {
            forallX(i,j,MJ_) {
                fprintf(fout, "%12.5e\n", (ac(gscaf,i,j,MJ_) + ac(gscaf,i+1,j,MJ_) + ac(gscaf,i,j+1,MJ_) + ac(gscaf,i+1,j+1,MJ_))/4); // move to vertex point
            }
        }
        if(cx==1&&cy==0) {
            forallX(i,j,MJ_) {
                fprintf(fout, "%12.5e\n", (ac(gscaf,i,j,MJ_) + ac(gscaf,i+1,j,MJ_))/2);  // move to vertex point
            }
        }
        if(cx==0&&cy==1) {
            forallX(i,j,MJ_) {
                fprintf(fout, "%12.5e\n", (ac(gscaf,i,j,MJ_) + ac(gscaf,i,j+1,MJ_))/2);  // move to vertex point
            }
        }
        if(cx==0&&cy==0) {
            forallX(i,j,MJ_) {
                fprintf(fout, "%12.5e\n", ac(gscaf,i,j,MJ_));  // already at vertex point
            }
        }
    }

    checkr(MPI_Allreduce(MPI_IN_PLACE, &iscrash, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD),"Allreduce"); // Update iscrash to all procs.

    if(iscrash) CRASH("Case undefined"); // If iscrash == 1, all procs crash

    if(quisoc()==0) {
        free(gscaf); // IMPORTANT: ONLY 0 DEALLOCATES !!
        fclose(fout);
    }

    checkr(MPI_Barrier(MPI_COMM_WORLD),"Barrier");

}
