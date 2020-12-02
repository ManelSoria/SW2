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


#define LUX_CENT (ac(u_in,i-1,j,SW->M)+ac(u_in,i+1,j,SW->M)-2*ac(u_in,i,j,SW->M))/(POW2(DXE(i,j)))
#define LUY_CENT (ac(u_in,i,j-1,SW->M)+ac(u_in,i,j+1,SW->M)-2*ac(u_in,i,j,SW->M))/(POW2(DYC(i,j)))
#define LVX_CENT (ac(v_in,i-1,j,SW->M)+ac(v_in,i+1,j,SW->M)-2*ac(v_in,i,j,SW->M))/(POW2(DXV(i,j)))
#define LVY_CENT (ac(v_in,i,j-1,SW->M)+ac(v_in,i,j+1,SW->M)-2*ac(v_in,i,j,SW->M))/(POW2(DYN(i,j)))
#define LEX_CENT (ac(e_in,i-1,j,SW->M)+ac(e_in,i+1,j,SW->M)-2*ac(e_in,i,j,SW->M))/(POW2(DXC(i,j)))
#define LEY_CENT (ac(e_in,i,j-1,SW->M)+ac(e_in,i,j+1,SW->M)-2*ac(e_in,i,j,SW->M))/(POW2(DYC(i,j)))

#define LUX_FORW (-5*ac(u_in,i+1,j,SW->M)+4*ac(u_in,i+2,j,SW->M)-ac(u_in,i+3,j,SW->M)+2*ac(u_in,i,j,SW->M))/(POW2((DXE(i+1,j)+DXE(i+2,j)+DXE(i+3,j))/3.0))
#define LUY_FORW (-5*ac(u_in,i,j+1,SW->M)+4*ac(u_in,i,j+2,SW->M)-ac(u_in,i,j+3,SW->M)+2*ac(u_in,i,j,SW->M))/(POW2((DYC(i,j+1)+DYC(i,j+2)+DYC(i,j+3))/3.0))
#define LVX_FORW (-5*ac(v_in,i+1,j,SW->M)+4*ac(v_in,i+2,j,SW->M)-ac(v_in,i+3,j,SW->M)+2*ac(v_in,i,j,SW->M))/(POW2((DXV(i+1,j)+DXV(i+2,j)+DXV(i+3,j))/3.0))
#define LVY_FORW (-5*ac(v_in,i,j+1,SW->M)+4*ac(v_in,i,j+2,SW->M)-ac(v_in,i,j+3,SW->M)+2*ac(v_in,i,j,SW->M))/(POW2((DYN(i,j+1)+DYN(i,j+2)+DYN(i,j+3))/3.0))
#define LEX_FORW (-5*ac(e_in,i+1,j,SW->M)+4*ac(e_in,i+2,j,SW->M)-ac(e_in,i+3,j,SW->M)+2*ac(e_in,i,j,SW->M))/(POW2((DXC(i+1,j)+DXC(i+2,j)+DXC(i+3,j))/3.0))
#define LEY_FORW (-5*ac(e_in,i,j+1,SW->M)+4*ac(e_in,i,j+2,SW->M)-ac(e_in,i,j+3,SW->M)+2*ac(e_in,i,j,SW->M))/(POW2((DYC(i,j+1)+DYC(i,j+2)+DYC(i,j+3))/3.0))

#define LUX_BACK (-5*ac(u_in,i-1,j,SW->M)+4*ac(u_in,i-2,j,SW->M)-ac(u_in,i-3,j,SW->M)+2*ac(u_in,i,j,SW->M))/(POW2((DXE(i,j)+DXE(i-1,j)+DXE(i-2,j))/3.0))
#define LUY_BACK (-5*ac(u_in,i,j-1,SW->M)+4*ac(u_in,i,j-2,SW->M)-ac(u_in,i,j-3,SW->M)+2*ac(u_in,i,j,SW->M))/(POW2((DYC(i,j)+DYC(i,j-1)+DYC(i,j-2))/3.0))
#define LVX_BACK (-5*ac(v_in,i-1,j,SW->M)+4*ac(v_in,i-2,j,SW->M)-ac(v_in,i-3,j,SW->M)+2*ac(v_in,i,j,SW->M))/(POW2((DXV(i,j)+DXV(i-1,j)+DXV(i-2,j))/3.0))
#define LVY_BACK (-5*ac(v_in,i,j-1,SW->M)+4*ac(v_in,i,j-2,SW->M)-ac(v_in,i,j-3,SW->M)+2*ac(v_in,i,j,SW->M))/(POW2((DYN(i,j)+DYN(i,j-1)+DYN(i,j-2))/3.0))
#define LEX_BACK (-5*ac(e_in,i-1,j,SW->M)+4*ac(e_in,i-2,j,SW->M)-ac(e_in,i-3,j,SW->M)+2*ac(e_in,i,j,SW->M))/(POW2((DXC(i,j)+DXC(i-1,j)+DXC(i-2,j))/3.0))
#define LEY_BACK (-5*ac(e_in,i,j-1,SW->M)+4*ac(e_in,i,j-2,SW->M)-ac(e_in,i,j-3,SW->M)+2*ac(e_in,i,j,SW->M))/(POW2((DYC(i,j)+DYC(i,j-1)+DYC(i,j-2))/3.0))

#define ZE(i,j) ac(dseta,i,j,SW->M)
#define DI(i,j) ac(div,i,j,SW->M)
#define VV(i,j) ac(v_in,i,j,SW->M)
#define UU(i,j) ac(u_in,i,j,SW->M)
#define LUU(i,j) ac(u_out,i,j,SW->M)
#define LVV(i,j) ac(v_out,i,j,SW->M)


void laplacian_vels(double *u_in, double *v_in, double  *u_out, double *v_out, sw *SW){
  int i,j;
  double *dseta, *div, a, b;
  dseta = dmem(SW->M);
  div = dmem(SW->M);
  setzero_scaf(dseta,SW->M);
  setzero_scaf(div,SW->M);

  forall(i,j,SW->M){
    a = (1.0/DXV(i+1,j));     b = -(1.0/DYC(i,j+1));
    ZE(i,j) = a*((VV(i+1,j) - VV(i,j)) + b*(DXC(i+1,j+1)*UU(i,j+1)-DXC(i+1,j)*UU(i,j)));

    a = (1.0/DXE(i,j));       b = (1.0/DYN(i,j));
    DI(i,j) = a*((UU(i,j) - UU(i-1,j)) + b*(DXV(i,j)*VV(i,j)-DXV(i,j-1)*VV(i,j-1)));
  }

  halo_update(dseta, SW->M);
  halo_update(div, SW->M);

  if(obN(SW->M)){
    j = SW->M->l1[2];
    for(i = SW->M->l0[1]; i<=SW->M->l1[1]; i++){
      DI(i,j+1) = (UU(i,j) - UU(i-1,j))/DXE(i,j);
    }
  }
  if(obS(SW->M)){
    j = SW->M->l0[2];
    for(i = SW->M->l0[1]; i<=SW->M->l1[1]; i++){
      ZE(i,j-1) = 0;
    }
  }

  forall(i,j,SW->M){
    a = -(1.0/DYN(i,j));      b = (1.0/DXC(i+1,j));
    LUU(i,j) = a*(ZE(i,j)-ZE(i,j-1)) + b*(DI(i+1,j)-DI(i,j));

    a = (1.0/DXV(i,j));       b = (1.0/DYC(i,j+1));
    LVV(i,j) = a*(ZE(i,j)-ZE(i-1,j)) + b*(DI(i,j+1)-DI(i,j));
  }

  halo_update(u_out,SW->M);
  halo_update(v_out,SW->M);

  if(obN(SW->M)){
    j = SW->M->l1[2];
    for(i = SW->M->l0[1]; i<=SW->M->l1[1]; i++){
      LUU(i,j+1) = LUU(i,j);
      LVV(i,j+1) = 0;
    }
  }
  if(obS(SW->M)){
    j = SW->M->l0[2];
    for(i = SW->M->l0[1]; i<=SW->M->l1[1]; i++){
      LUU(i,j-1) = LUU(i,j);
      LVV(i,j-1) = 0;
    }
  }

  free(dseta);
  free(div);
}


void laplacian_eta(double *e_in, double *e_out, sw *SW){
  int i,j;

  double *psi = dmem(SW->M);
  double *xi = dmem(SW->M);

  forall(i,j,SW->M){
    ac(psi,i,j,SW->M) = rM(LATE(i,j))/rZ(LATE(i,j))*(ac(e_in,i+1,j,SW->M)-ac(e_in,i,j,SW->M))/(DXC(i,j));
    ac(xi,i,j,SW->M) = rZ(LATN(i,j))/rM(LATN(i,j))*(ac(e_in,i,j+1,SW->M)-ac(e_in,i,j,SW->M))/(DYC(i,j));
  }

  forall(i,j,SW->M){
    ac(e_out,i,j,SW->M) = 1.0/(rM(LATC(i,j))*rZ(LATC(i,j)))*((ac(psi,i,j,SW->M)-ac(psi,i-1,j,SW->M))/(DXE(i,j))+(ac(xi,i,j,SW->M)-ac(xi,i,j-1,SW->M))/(DYN(i,j)));
  }

  halo_update(e_out, SW->M);

  if(obN(SW->M)){
    j = SW->M->l1[2]+1;
    for(i = SW->M->l0[1]; i <= SW->M->l1[1]; i++){
      ac(e_out,i,j,SW->M) = ac(e_out,i,j-1,SW->M);

      ac(e_out,i,j+1,SW->M) = ac(e_out,i,j,SW->M);
    }
  }
  if(obS(SW->M)){
    j = SW->M->l0[2]-1;
    for(i = SW->M->l0[1]; i <= SW->M->l1[1]; i++){
      ac(e_out,i,j,SW->M) = ac(e_out,i,j+1,SW->M);

      ac(e_out,i,j-1,SW->M) = ac(e_out,i,j-1,SW->M);
    }
  }

  free(psi);
  free(xi);

}
