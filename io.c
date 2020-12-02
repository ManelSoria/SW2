// SW2
// 2018-2020
// Manel Soria, Arnau Prat, Arnau Sabates, Marc Andres, Enrique Garcia-Melendo
// UPC - ESEIAAT - TUAREG

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <math.h>
#include "mpi.h"

#include "sppde.h"
#include "sppde_parser.h"
#include "sw.h"

#define MAGICSW1 27 // eta, u, v, tracer
#define MAGICSW2 28 // eta, u, v, tracer, perturbs


void save_binary_sw(int frame, int nstep,double timedouble,char *bname, char *output_folder,sw *SW) { // bname: name without output folder
    FILE *fout;
    int info=0;

    if (quisoc()==0) { // open and write header
        char outname[MAXOUTLEN];
        snprintf(outname,MAXOUTLEN,"%s/%s.%04d.bin",output_folder,basename(bname),frame);
        fout=fopen(outname,"w");
        if (fout==NULL) CRASH("file cannot be opened <%s>",outname);
        int magic=MAGICSW2;
        fwrite( (void *)&magic,sizeof(int),1,fout );
        fwrite( (void *)&nstep,sizeof(int),1,fout );
        fwrite( (void *)&timedouble,sizeof(double),1,fout );
        fwrite( (void *)&(SW->nx),sizeof(int),1,fout );
        fwrite( (void *)&(SW->ny),sizeof(int),1,fout );
        pprintf("save_binary_sw saved frame=%d nstep=%d timedouble=%e outname=<%s> \n",frame,nstep,timedouble,outname);
    }

    // MAGICSW2:
    fwrite_scaf(SW->eta    , fout,SW->M);
    fwrite_scaf(SW->u      , fout,SW->M);
    fwrite_scaf(SW->v      , fout,SW->M);
    fwrite_scaf(SW->tracer , fout,SW->M);
    fwrite_scaf(SW->perturbs, fout,SW->M);

    if(info>1){
        pprintf("save_binary_sw saved binary file with: \n");
        stats_scaf(SW->eta,SW->M,"eta"); // Print field information
        stats_scaf(SW->u,SW->M,"u"); // Print field information
        stats_scaf(SW->v,SW->M,"v"); // Print field information
        stats_scaf(SW->tracer,SW->M,"tracer"); // Print field information
        stats_scaf(SW->perturbs,SW->M,"perturbs"); // Print field information
    }

    if (quisoc()==0) {
        fclose(fout);
    }
}

// loads only the header of a binary sw file
// frame: frame to be loaded
// bfname, binaries_folder: name of the file and folder
// output:
// nstep: writes the step loaded
// timedouble: idem time
// returns the magic number

int load_binary_header_core(FILE *fin, int *nstep, double *timedouble, int *nxr, int *nyr) {
    int nr;
    int magic;
    nr=fread( (void *)&magic,sizeof(int),1,fin );
    if (nr!=1) CRASH("File to be loaded is corrupted");
    if (magic<MAGICSW1 || magic>MAGICSW2) CRASH("Wrong magic = %d should be %d or %d",magic,MAGICSW1,MAGICSW2);

    nr=fread( (void *)nstep,sizeof(int),1,fin);
    if (nr!=1) CRASH("File to be loaded is corrupted");

    nr=fread( (void *)timedouble,sizeof(double),1,fin );
    if (nr!=1) CRASH("File to be loaded is corrupted");

    nr=fread( (void *)nxr,sizeof(int),1,fin );
    if (nr!=1) CRASH("File to be loaded is corrupted");

    nr=fread( (void *)nyr,sizeof(int),1,fin );
    if (nr!=1) CRASH("File to be loaded is corrupted");

    return(magic);
}

// returns 0 if no errors
int load_binary_header_sw(int frame,int *nstep,double *timedouble, char *bname, char *binaries_folder) {
    FILE *fin;
    int info=0;
    if (quisoc()==0) {
         char inname[MAXOUTLEN];
        snprintf(inname,MAXOUTLEN,"%s/%s.%04d.bin",binaries_folder,basename(bname),frame);
        fin=fopen(inname,"r");
        if (fin==NULL)
            return(1);
        int nxr, nyr; // variables not currently used (dummy variables)
        int magic=load_binary_header_core(fin,nstep,timedouble, &nxr, &nyr);
        if (info) pprintf("load_binary_header_sw found magic=%d frame=%d nstep=%d timedouble=%e inname=<%s> \n",magic,frame,*nstep,*timedouble,inname);
        fclose(fin);
   }
    return(0);
}

// returns 0 if file was loaded, 1 if the file does not exist or is corrupted
int load_binary_sw(int frame, int *nstep,double *timedouble,char *bname, char *binaries_folder,sw *SW) { // bname: name without binaries folder
    FILE *fin;
    int info=0;
    int nr;
    int magic,magicG;

    if (info) pprintf("load_binary_sw stating\n");
    if (quisoc()==0) { // proc 0 open and write header
        char inname[MAXOUTLEN];
        snprintf(inname,MAXOUTLEN,"%s/%s.%04d.bin",binaries_folder,basename(bname),frame);
        fin=fopen(inname,"r");
        if (fin==NULL)
            return(1);
        int nxr, nyr;
        magic=load_binary_header_core(fin,nstep,timedouble, &nxr, &nyr);
        if (nxr!=SW->nx) CRASH("oops nxr=%d should be %d",nxr,SW->nx);
        if (nyr!=SW->ny) CRASH("oops nyr=%d should be %d",nyr,SW->ny);
        if (info) pprintf("load_binary_sw loaded frame=%d nstep=%d timedouble=%e inname=<%s> \n",frame,*nstep,*timedouble,inname);
    } else magic=0; // the rest have a 0 here

    checkr(MPI_Allreduce(&magic,&magicG,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD),"Allreduce");

    if (info) pprintf("load_binary_sw got header magicG=%d \n",magicG);

    switch(magicG) { // binary file legacy compatible
        case MAGICSW1:
            if (info) pprintf("reading binary magic=%d \n",MAGICSW1);
            fread_scaf(SW->eta    , fin,SW->M);
            fread_scaf(SW->u      , fin,SW->M);
            fread_scaf(SW->v      , fin,SW->M);
            fread_scaf(SW->perturbs , fin,SW->M); // tracer field
            setzero_scaf(SW->tracer,SW->M); // no tracer
        break;
        case MAGICSW2:
            if (info) pprintf("reading binary magic=%d \n",MAGICSW2);
            fread_scaf(SW->eta    , fin,SW->M);
            fread_scaf(SW->u      , fin,SW->M);
            fread_scaf(SW->v      , fin,SW->M);
            fread_scaf(SW->tracer , fin,SW->M);
            fread_scaf(SW->perturbs , fin,SW->M);
        break;
    }

    if(info>1){
        pprintf("load_binary_sw loaded binary file with: \n");
        stats_scaf(SW->tracer,SW->M,"tracer"); // Print field information
        stats_scaf(SW->eta,SW->M,"eta"); // Print field information
        stats_scaf(SW->u,SW->M,"u"); // Print field information
        stats_scaf(SW->v,SW->M,"v"); // Print field information
    }
    checkr(MPI_Bcast(nstep,1,MPI_INT,0,MPI_COMM_WORLD),"Bcast"); // Add Bcast to sppde_core
    checkr(MPI_Bcast(timedouble,1,MPI_DOUBLE,0,MPI_COMM_WORLD),"Bcast"); // Add Bcast to sppde_core
    if (quisoc()==0) {
        fclose(fin);
    }
    return(0);
}

void read_wind(FILE *fin, sw *SW) {

/*

Reads a zonal wind file and updates the array of the Shallow World dedicated to it.

- char *fname (input) The name of the file.
- sw *SW (input/output) The Shallow World.

Requires sppde_tinyexpr and sppde_parser.

*/

    int info=0;
    double *intable, *outtable;

    int nr,nc;
    intable=p_alloc_gettable(fin, &nr, &nc,1);
    if(info>1)
        for (int i=0;i<=nr*nc-1;i++)
              pprintf("i=%d intable=%e\n",i,intable[i]);

    if(info) pprintf("Table allocated.\n");

    int i,j;

    if(SW->polar == 0){
      if(info) pprintf("Limits of table = %e and %e\n",intable[0],intable[nr-1]);
      if(rad2deg(latc(-1))<intable[0] || rad2deg(latc(SW->ny-1+2))>intable[nr-1]) CRASH("Latitude value (%e, %e) out of bounds (%e, %e).",rad2deg(latc(-1)),intable[0],rad2deg(latc(SW->ny-1+2)),intable[SW->ny-1+2]);
      forallh(i,j,SW->M){
          WINDX(i,j) = linterp(rad2deg(latc(j)),intable, intable+nr, nr); // SW->wind of size the number of points
          // There is the `linterp_vec(...)` function which may be more elegant than linterp, as it can deal with vectors and interpolations altogether
      }
    }
    if (SW->polar == 1 || SW->polar == -1){
      double angle;
      double mid_LON, mid_LAT;
      double dist;
      mid_LON  = MEAN(SW->lon0,SW->lon1);
      mid_LAT  = MEAN(SW->lat0,SW->lat1);
      double k = log(0.01)/(sqrt(2.0)*SW->lat1 - SW->beta_sponge);
      forallh(i,j,SW->M){
        if(rad2deg(POLAR_LATC(i,j))<intable[0] || rad2deg(POLAR_LATC(i,j))>intable[nr-1]){
          WINDX(i,j) = 0;
          WINDY(i,j) = 0;
        }else{
          angle = acos((LATC(i,j)-mid_LAT)/sqrt(POW2(LATC(i,j)-mid_LAT)+POW2(LONC(i,j)-mid_LON)));
          if ((LONC(i,j)-mid_LON)>=0){
            angle =  angle;
          }else{
            angle = -angle;
          }
          if (SW->polar == 1){
            WINDX(i,j) = -cos(angle)*linterp(rad2deg(POLAR_LATC(i,j)),intable, intable+nr, nr);
            WINDY(i,j) =  sin(angle)*linterp(rad2deg(POLAR_LATC(i,j)),intable, intable+nr, nr);
          }else{
            WINDX(i,j) =  cos(angle)*linterp(rad2deg(POLAR_LATC(i,j)),intable, intable+nr, nr);
            WINDY(i,j) = -sin(angle)*linterp(rad2deg(POLAR_LATC(i,j)),intable, intable+nr, nr);
          }
          dist = fabs(get_great_angle(MEAN(SW->lon0,SW->lon1), MEAN(SW->lat0,SW->lat1), LONC(i,j), LATC(i,j)))-SW->beta_sponge;
          if (dist >= 0){
            WINDX(i,j) = WINDX(i,j)*exp(k*dist);
            WINDY(i,j) = WINDY(i,j)*exp(k*dist);
          }
        }
      }
    }


    free(intable);
    pprintf("Zonal winds added.\n");
}


void read_Gaussians(FILE *fin, sw *SW) {

/*

Reads a file containing the parameters of the Gaussian perturbations to be carried out on the surface of the Shallow World.

- char *fname (input) The name of the file.
- sw *SW (input/output) The Shallow World.

Requires sppde_tinyexpr and sppde_parser.

*/


    int info=0;
    double *table;


    int nr,nc;
    SW->Gaussians = p_alloc_gettable(fin, &nr, &nc,0);
    if(nc != N_COLS_GAUSS) CRASH("Each SW perturbation should be defined with %d columns. The file containing the perturbations is not %d columns long.", N_COLS_GAUSS, N_COLS_GAUSS); // ARNAU

    SW->num_Gaussians = nr; // Update number of perturbations
    pprintf("Gaussians table read. There is/are %d gaussian perturbations in the domain.\n",SW->num_Gaussians);
}


void read_vortices(FILE *fin, sw *SW) {

/*

Reads a file containing the parameters of the vortices.

- char *fname (input) The name of the file.
- sw *SW (input/output) The Shallow World.

Requires sppde_tinyexpr and sppde_parser.

*/


    int info=0;
    double *table;


    int nr,nc;
    SW->vortices = p_alloc_gettable(fin, &nr, &nc,0);
    if(nc != N_COLS_VORTICES) CRASH("Each SW vortex should be defined with %d columns. The file containing the perturbations is not %d columns long.", N_COLS_VORTICES, N_COLS_VORTICES); // ARNAUMARC

    SW->num_vortices = nr; // Update number of perturbations
    pprintf("Vortices table read. There is/are %d vortices in the domain.\n",SW->num_vortices);
} // ARNAUMARC
