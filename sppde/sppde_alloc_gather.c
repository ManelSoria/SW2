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


// tested by sppde_gs_field

// proc. 0 allocs a large memory area where it receives all the data of the distributed field x
// MJ is the map to access the data gathered
// all procs have MJ, but ONLY PROC. 0 returns the area where the memory is allocated
// if eh==1, external halo is gathered
double *alloc_gather(double *x,map *M, map *MJ,int eh) {

  joinM(M,MJ); // creates MJ from M
  double *y;

  if (quisoc()==0) {
    y=dmem(MJ); // gathered field
  }
  else
    y=NULL;

  gs_field(0,eh,x,M,y,MJ);

  return(y);
}
