// SW2
// 2018-2020
// Manel Soria, Arnau Prat, Arnau Sabates, Marc Andres, Enrique Garcia-Melendo
// UPC - ESEIAAT - TUAREG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libgen.h>
#include "mpi.h"

#include "sppde.h"
#include "sw.h"

static int firstcall=1;
static FILE *fo;
static char fon[MAXOUTLEN];

double compare_ref(sw *SW,int n,double t, char *bname,char *output_folder) {
	int info=0;
	double r=-1; // no comparison made
    if (firstcall==1) {
        sprintf(fon,"%s/%s.ref.txt",output_folder,bname);
        fo=fopen(fon,"w");
        if (fo==NULL) CRASH("Can't open %s",fon);
        fclose(fo);
        firstcall=0;
    }
    switch (SW->reftype) {
        case 0:
        break;
        case 1: // compare with text files
 			for (int i=0;i<=SW->reftable_nr-1;i++) {
 				int frame=(int)(SW->reftable[i*SW->reftable_nc]);
 				double reftime=SW->reftable[i*SW->reftable_nc+1];
 				if (info) pprintf("compare_ref i=%d frame=%d reftime=%f\n",i,frame,reftime);
 				double epsilon = fabs(t-reftime);
 				if (epsilon < SW->Dt/2 || (epsilon <= SW-> Dt/2 && t<reftime))  { // should we compare ?
 					char bname[MAXOUTLEN],fname[MAXOUTLEN];
 					sprintf(bname,SW->refname,frame);
 					sprintf(fname,"%s/%s",output_folder,bname);
 					if (info) pprintf("compare_ref COMPARE t=%f reftime=%f DT=%f frame=%d fname=<%s>\n",t,reftime,SW->Dt,frame,fname);
 					double *fread;
 					fread=dmem(SW->M); // alloc memory for the reference field
 					text_read_scaf(fread,fname,SW->M);
 					check_stats_scaf(fread,SW->M,"fread",n,-1e10,1e10,1);
 					r=norm_dif_scaf(SW->eta,fread,2,SW->M);
                    if (info) pprintf("compare_ref t %f reftime %f frame %d fname %s delta %e\n",t,reftime,frame,fname,r);
 					if (info) pprintf("compare_ref COMPARE r=%e\n",r);
 					free(fread);
                    if (quisoc()==0) {
                        fo=fopen(fon,"a");
                        if (fo==NULL) CRASH("Can't open %s",fon);
     					fprintf(fo,"frame %d time %e delta %e\n",frame,t,r);
                        fflush(fo);
                        fclose(fo);
                    }
 					break; // exit loop (we already found it)
 				}
 			}
        break;
        default: CRASH("wrong reftype=%d",SW->reftype);
    }

    return(r);
}
