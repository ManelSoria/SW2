// sppde: structured parallel pde solver
// Manel Soria 2017 - UPC - ESEIAAT - TUAREG

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
