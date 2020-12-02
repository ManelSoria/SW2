// SW2
// 2018-2020
// Manel Soria, Arnau Prat, Arnau Sabates, Marc Andres, Enrique Garcia-Melendo
// UPC - ESEIAAT - TUAREG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "sppde.h"
#include "sw.h"

void destroy(sw *SW) {

/*

Deallocs all the memory.

- sw *SW (input/output) The Shallow World.

*/

    free(SW->Dxe);
    free(SW->Dxc);
    free(SW->Dxv);
    free(SW->Dyn);
    free(SW->Dyc);
    free(SW->Dyv);
    free(SW->xe);
    free(SW->xc);
    free(SW->xv);
    free(SW->yn);
    free(SW->yc);
    free(SW->yv);


    free(SW->ge);
    free(SW->gn);
    free(SW->gc);
    free(SW->gv);
    free(SW->fe);
    free(SW->fn);
    free(SW->fc);
    free(SW->fv);

    free(SW->tracer);
    free(SW->perturbs);
    free(SW->u);
    free(SW->v);
    free(SW->eta);
    free(SW->hB);
    free(SW->windx);
    free(SW->windy);
    free(SW->Gaussians);
    free(SW->vortices);


    switch (SW->reftype) {
        case 0:
        break;
        case 1: // compare with binary files
            free(SW->reftable);
        break;
        default: CRASH("wrong reftype=%d",SW->reftype);
    }

    pprintf("All memory freed.\n");

}
