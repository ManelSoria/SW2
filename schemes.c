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
#include <math.h>
#include "mpi.h"

#include "sppde.h"
#include "sw.h"

// In this file there are the different bodys and schemes of the program


// ........------------- SCHEMES -------------........

#define funPsi_MUSCL(r) (MAX(0,MIN(MIN((1+r)/2.0,2.0*r),2.0)))                  // MUSCL, TVD
#define funPsi_SUPERBEE(r) (MAX(MAX(0,MIN(2.0*r,1.0)),MAX(0,MIN(r,2.0))))       // SUPERBEE, TVD
#define funPsi_MINMOD(r) (MIN(r,1))                                             //MINMOD
#define funPsi_VANLEER(r) ((r+fabs(r))/(1.0+r))                                 // Van Leer flux limiter
#define funPsi_KOREN(r) (MAX(0,MIN(2.0*r,MIN((1+2.0*r)/3.0,2.0))))
#define funPsi_OSPRE(r) ((1.5*(POW2(r)+r))/(POW2(r)+r+1.0))
#define funPsi_VANALBADA1(r) ((POW2(r)+r)/(POW2(r)+1.0))
#define funPsi_UMIST(r) (MAX(0,MIN(MIN(2.0,2.0*r),MIN(0.25+0.75*r,0.75+0.25*r))))
#define funPsi_UPWIND(r) (0)
#define funPsi_LAXWENDROFF(r) (1.0)





// ........------------- BODYS -------------........
#define BODYH(sch)                                                              \
forall(i,j,SW->M) {                                                             \
                                                                                \
    CwP = 0.5 * (UWINDX(i-1,j) + fabs(UWINDX(i-1,j))) * Dt / DXE(i,j);          \
    CwM = 0.5 * (UWINDX(i-1,j) - fabs(UWINDX(i-1,j))) * Dt / DXE(i,j);          \
    CeP = 0.5 * (UWINDX(i,j) + fabs(UWINDX(i,j))) * Dt / DXE(i,j);              \
    CeM = 0.5 * (UWINDX(i,j) - fabs(UWINDX(i,j))) * Dt / DXE(i,j);              \
    CsP = 0.5 * (VWINDY(i,j-1) + fabs(VWINDY(i,j-1))) * Dt / DYN(i,j);          \
    CsM = 0.5 * (VWINDY(i,j-1) - fabs(VWINDY(i,j-1))) * Dt / DYN(i,j);          \
    CnP = 0.5 * (VWINDY(i,j) + fabs(VWINDY(i,j))) * Dt / DYN(i,j);              \
    CnM = 0.5 * (VWINDY(i,j) - fabs(VWINDY(i,j))) * Dt / DYN(i,j);              \
                                                                                \
    rxPprev = r2rP(H(i-2,j), H(i-1,j), H(i,j));                                 \
    rxMprev = r2rM(H(i-1,j), H(i,j), H(i+1,j));                                 \
    rxP = r2rP(H(i-1,j), H(i,j), H(i+1,j));                                     \
    rxM = r2rM(H(i,j), H(i+1,j), H(i+2,j));                                     \
    ryPprev = r2rP(H(i,j-2), H(i,j-1), H(i,j));                                 \
    ryMprev = r2rM(H(i,j-1), H(i,j), H(i,j+1));                                 \
    ryP = r2rP(H(i,j-1), H(i,j), H(i,j+1));                                     \
    ryM = r2rM(H(i,j), H(i,j+1), H(i,j+2));                                     \
                                                                                \
    BwP = B2BPG(sch,H(i-1,j),H(i,j),rxPprev,CwP);                               \
    BwM = B2BMG(sch,H(i-1,j),H(i,j),rxMprev,CwM);                               \
    BeP = B2BPG(sch,H(i,j),H(i+1,j),rxP,CeP);                                   \
    BeM = B2BMG(sch,H(i,j),H(i+1,j),rxM,CeM);                                   \
    BsP = B2BPG(sch,H(i,j-1),H(i,j),ryPprev,CsP);                               \
    BsM = B2BMG(sch,H(i,j-1),H(i,j),ryMprev,CsM);                               \
    BnP = B2BPG(sch,H(i,j),H(i,j+1),ryP,CnP);                                   \
    BnM = B2BMG(sch,H(i,j),H(i,j+1),ryM,CnM);                                   \
                                                                                \
    Deta(i,j)=CwP*BwP+CwM*BwM-CeP*BeP-CeM*BeM+CsP*BsP+CsM*BsM-CnP*BnP-CnM*BnM;  \
                                                                                \
    if (SW->polar!=0) {                                                         \
        dist = get_great_angle(MEAN(SW->lon0,SW->lon1),MEAN(SW->lat0,SW->lat1),LONC(i,j),LATC(i,j))-SW->beta_sponge; \
        if(dist > 0){                                                           \
          sigma_dis = SW->sigma_sponge*(1.0-exp(-SW->tau_sponge*(dist)));       \
        }else{                                                                  \
          sigma_dis = SW->sigma_dis;                                            \
        }                                                                       \
} else {                                                                        \
        sigma_dis = SW->sigma_dis;                                              \
    }                                                                           \
    Deta(i,j) -= 2.0*SW->Dt*(sigma_dis)*ETA(i,j);                               \
}
// End of H

// General body for the advection of Tracer
#define BODYTRACER(sch)                                                         \
forall(i,j,SW->M) {                                                             \
                                                                                \
    CwP = 0.5 * (UWINDX(i-1,j) + fabs(UWINDX(i-1,j))) * Dt / DXE(i,j);          \
    CwM = 0.5 * (UWINDX(i-1,j) - fabs(UWINDX(i-1,j))) * Dt / DXE(i,j);          \
    CeP = 0.5 * (UWINDX(i,j) + fabs(UWINDX(i,j))) * Dt / DXE(i,j);              \
    CeM = 0.5 * (UWINDX(i,j) - fabs(UWINDX(i,j))) * Dt / DXE(i,j);              \
    CsP = 0.5 * (VWINDY(i,j-1) + fabs(VWINDY(i,j-1))) * Dt / DYN(i,j);          \
    CsM = 0.5 * (VWINDY(i,j-1) - fabs(VWINDY(i,j-1))) * Dt / DYN(i,j);          \
    CnP = 0.5 * (VWINDY(i,j) + fabs(VWINDY(i,j))) * Dt / DYN(i,j);              \
    CnM = 0.5 * (VWINDY(i,j) - fabs(VWINDY(i,j))) * Dt / DYN(i,j);              \
                                                                                \
    rxPprev = r2rP(TRACER(i-2,j), TRACER(i-1,j), TRACER(i,j));                  \
    rxMprev = r2rM(TRACER(i-1,j), TRACER(i,j), TRACER(i+1,j));                  \
    rxP = r2rP(TRACER(i-1,j), TRACER(i,j), TRACER(i+1,j));                      \
    rxM = r2rM(TRACER(i,j), TRACER(i+1,j), TRACER(i+2,j));                      \
    ryPprev = r2rP(TRACER(i,j-2), TRACER(i,j-1), TRACER(i,j));                  \
    ryMprev = r2rM(TRACER(i,j-1), TRACER(i,j), TRACER(i,j+1));                  \
    ryP = r2rP(TRACER(i,j-1), TRACER(i,j), TRACER(i,j+1));                      \
    ryM = r2rM(TRACER(i,j), TRACER(i,j+1), TRACER(i,j+2));                      \
                                                                                \
    BwP = B2BPG(sch,TRACER(i-1,j),TRACER(i,j),rxPprev,CwP);                     \
    BwM = B2BMG(sch,TRACER(i-1,j),TRACER(i,j),rxMprev,CwM);                     \
    BeP = B2BPG(sch,TRACER(i,j),TRACER(i+1,j),rxP,CeP);                         \
    BeM = B2BMG(sch,TRACER(i,j),TRACER(i+1,j),rxM,CeM);                         \
    BsP = B2BPG(sch,TRACER(i,j-1),TRACER(i,j),ryPprev,CsP);                     \
    BsM = B2BMG(sch,TRACER(i,j-1),TRACER(i,j),ryMprev,CsM);                     \
    BnP = B2BPG(sch,TRACER(i,j),TRACER(i,j+1),ryP,CnP);                         \
    BnM = B2BMG(sch,TRACER(i,j),TRACER(i,j+1),ryM,CnM);                         \
                                                                                \
    Dtracer(i,j)=CwP*BwP+CwM*BwM-CeP*BeP-CeM*BeM+CsP*BsP+CsM*BsM-CnP*BnP-CnM*BnM; \
                                                                                \
    if (SW->polar!=0) {                                                         \
        dist = get_great_angle(MEAN(SW->lon0,SW->lon1), MEAN(SW->lat0,SW->lat1), LONC(i,j), LATC(i,j))-SW->beta_sponge; \
        if(dist > 0){                                                           \
          sigma_dis = SW->sigma_sponge*(1.0-exp(-SW->tau_sponge*(dist)));       \
        }else{                                                                  \
          sigma_dis = 0;                                                        \
        }                                                                       \
    } else {                                                                    \
        sigma_dis = 0;                                                          \
    }                                                                           \
    Dtracer(i,j) -= 2.0*SW->Dt*(sigma_dis)*TRACER(i,j);                         \
}
// End of tracer


// Generib body for U
#define BODYU(sch)                                                              \
forall(i,j,SW->M) {                                                             \
    CwP = 0.5 * (0.5 * (UWINDX(i-1,j) + fabs(UWINDX(i-1,j))) + 0.5 * (UWINDX(i,j) + fabs(UWINDX(i,j)))) * Dt / DXC(i+1,j); \
    CwM = 0.5 * (0.5 * (UWINDX(i-1,j) - fabs(UWINDX(i-1,j))) + 0.5 * (UWINDX(i,j) - fabs(UWINDX(i,j)))) * Dt / DXC(i+1,j); \
    CeP = 0.5 * (0.5 * (UWINDX(i,j) + fabs(UWINDX(i,j))) + 0.5 * (UWINDX(i+1,j) + fabs(UWINDX(i+1,j)))) * Dt / DXC(i+1,j); \
    CeM = 0.5 * (0.5 * (UWINDX(i,j) - fabs(UWINDX(i,j))) + 0.5 * (UWINDX(i+1,j) - fabs(UWINDX(i+1,j)))) * Dt / DXC(i+1,j); \
    CsP = 0.5 * (0.5 * (VWINDY(i,j-1) + fabs(VWINDY(i,j-1))) + 0.5 * (VWINDY(i+1,j-1) + fabs(VWINDY(i+1,j-1)))) * Dt / DYV(i,j); \
    CsM = 0.5 * (0.5 * (VWINDY(i,j-1) - fabs(VWINDY(i,j-1))) + 0.5 * (VWINDY(i+1,j-1) - fabs(VWINDY(i+1,j-1)))) * Dt / DYV(i,j); \
    CnP = 0.5 * (0.5 * (VWINDY(i,j) + fabs(VWINDY(i,j))) + 0.5 * (VWINDY(i+1,j) + fabs(VWINDY(i+1,j)))) * Dt / DYV(i,j); \
    CnM = 0.5 * (0.5 * (VWINDY(i,j) - fabs(VWINDY(i,j))) + 0.5 * (VWINDY(i+1,j) - fabs(VWINDY(i+1,j)))) * Dt / DYV(i,j); \
                                                                                \
    rxPprev = r2rP(UWINDX(i-2,j), UWINDX(i-1,j), UWINDX(i,j));                  \
    rxMprev = r2rM(UWINDX(i-1,j), UWINDX(i,j), UWINDX(i+1,j));                  \
    rxP = r2rP(UWINDX(i-1,j), UWINDX(i,j), UWINDX(i+1,j));                      \
    rxM = r2rM(UWINDX(i,j), UWINDX(i+1,j), UWINDX(i+2,j));                      \
    ryPprev = r2rP(UWINDX(i,j-2), UWINDX(i,j-1), UWINDX(i,j));                  \
    ryMprev = r2rM(UWINDX(i,j-1), UWINDX(i,j), UWINDX(i,j+1));                  \
    ryP = r2rP(UWINDX(i,j-1), UWINDX(i,j), UWINDX(i,j+1));                      \
    ryM = r2rM(UWINDX(i,j), UWINDX(i,j+1), UWINDX(i,j+2));                      \
                                                                                \
    BwP = B2BPG(sch,UWINDX(i-1,j),UWINDX(i,j),rxPprev,CwP);                          \
    BwM = B2BMG(sch,UWINDX(i-1,j),UWINDX(i,j),rxMprev,CwM);                          \
    BeP = B2BPG(sch,UWINDX(i,j),UWINDX(i+1,j),rxP,CeP);                              \
    BeM = B2BMG(sch,UWINDX(i,j),UWINDX(i+1,j),rxM,CeM);                              \
    BsP = B2BPG(sch,UWINDX(i,j-1),UWINDX(i,j),ryPprev,CsP);                          \
    BsM = B2BMG(sch,UWINDX(i,j-1),UWINDX(i,j),ryMprev,CsM);                          \
    BnP = B2BPG(sch,UWINDX(i,j),UWINDX(i,j+1),ryP,CnP);                              \
    BnM = B2BMG(sch,UWINDX(i,j),UWINDX(i,j+1),ryM,CnM);                              \
                                                                                \
    adv_Pu1=CwP*BwP+CwM*BwM-CeP*BeP-CeM*BeM+CsP*BsP+CsM*BsM-CnP*BnP-CnM*BnM;    \
    adv_Pu2=Dt*(UWINDX(i,j)*(UWINDX(i+1,j)-UWINDX(i-1,j))/(DXE(i,j)+DXE(i+1,j))+UWINDX(i,j)*((VWINDY(i,j)+VWINDY(i+1,j))/2.-(VWINDY(i,j-1)+VWINDY(i+1,j-1))/2.)/DYV(i+1,j)); \
    DU_adv(i,j) = adv_Pu1 + adv_Pu2;                                            \
}
// End of U


// Body for V
#define BODYV(sch)                                                              \
forall(i,j,SW->M) {                                                             \
    CwP = 0.5 * (0.5 * (UWINDX(i-1,j+1) + fabs(UWINDX(i-1,j+1))) + 0.5 * (UWINDX(i-1,j) + fabs(UWINDX(i-1,j)))) * Dt / DXV(i,j); \
    CwM = 0.5 * (0.5 * (UWINDX(i-1,j+1) - fabs(UWINDX(i-1,j+1))) + 0.5 * (UWINDX(i-1,j) - fabs(UWINDX(i-1,j)))) * Dt / DXV(i,j); \
    CeP = 0.5 * (0.5 * (UWINDX(i,j+1) + fabs(UWINDX(i,j+1))) + 0.5 * (UWINDX(i,j) + fabs(UWINDX(i,j)))) * Dt / DXV(i,j); \
    CeM = 0.5 * (0.5 * (UWINDX(i,j+1) - fabs(UWINDX(i,j+1))) + 0.5 * (UWINDX(i,j) - fabs(UWINDX(i,j)))) * Dt / DXV(i,j); \
    CsP = 0.5 * (0.5 * (VWINDY(i,j-1) + fabs(VWINDY(i,j-1))) + 0.5 * (VWINDY(i,j) + fabs(VWINDY(i,j)))) * Dt / DYC(i,j); \
    CsM = 0.5 * (0.5 * (VWINDY(i,j-1) - fabs(VWINDY(i,j-1))) + 0.5 * (VWINDY(i,j) - fabs(VWINDY(i,j)))) * Dt / DYC(i,j); \
    CnP = 0.5 * (0.5 * (VWINDY(i,j) + fabs(VWINDY(i,j))) + 0.5 * (VWINDY(i,j+1) + fabs(VWINDY(i,j+1)))) * Dt / DYC(i,j); \
    CnM = 0.5 * (0.5 * (VWINDY(i,j) - fabs(VWINDY(i,j))) + 0.5 * (VWINDY(i,j+1) - fabs(VWINDY(i,j+1)))) * Dt / DYC(i,j); \
                                                                                \
    rxPprev = r2rP(VWINDY(i-2,j), VWINDY(i-1,j), VWINDY(i,j));                  \
    rxMprev = r2rM(VWINDY(i-1,j), VWINDY(i,j), VWINDY(i+1,j));                  \
    rxP = r2rP(VWINDY(i-1,j), VWINDY(i,j), VWINDY(i+1,j));                      \
    rxM = r2rM(VWINDY(i,j), VWINDY(i+1,j), VWINDY(i+2,j));                      \
    ryPprev = r2rP(VWINDY(i,j-2), VWINDY(i,j-1), VWINDY(i,j));                  \
    ryMprev = r2rM(VWINDY(i,j-1), VWINDY(i,j), VWINDY(i,j+1));                  \
    ryP = r2rP(VWINDY(i,j-1), VWINDY(i,j), VWINDY(i,j+1));                      \
    ryM = r2rM(VWINDY(i,j), VWINDY(i,j+1), VWINDY(i,j+2));                      \
                                                                                \
    BwP = B2BPG(sch,VWINDY(i-1,j),VWINDY(i,j),rxPprev,CwP);                          \
    BwM = B2BMG(sch,VWINDY(i-1,j),VWINDY(i,j),rxMprev,CwM);                          \
    BeP = B2BPG(sch,VWINDY(i,j),VWINDY(i+1,j),rxP,CeP);                              \
    BeM = B2BMG(sch,VWINDY(i,j),VWINDY(i+1,j),rxM,CeM);                              \
    BsP = B2BPG(sch,VWINDY(i,j-1),VWINDY(i,j),ryPprev,CsP);                          \
    BsM = B2BMG(sch,VWINDY(i,j-1),VWINDY(i,j),ryMprev,CsM);                          \
    BnP = B2BPG(sch,VWINDY(i,j),VWINDY(i,j+1),ryP,CnP);                              \
    BnM = B2BMG(sch,VWINDY(i,j),VWINDY(i,j+1),ryM,CnM);                              \
                                                                                \
    adv_Pv1 = CwP * BwP + CwM * BwM - CeP * BeP - CeM * BeM + CsP * BsP + CsM * BsM - CnP * BnP - CnM * BnM; \
    adv_Pv2 = Dt * ( VWINDY(i,j) * ( ( UWINDX(i,j) + UWINDX(i,j+1) ) / 2. - ( UWINDX(i-1,j) + UWINDX(i-1,j+1) ) / 2. ) / DXV(i,j) + VWINDY(i,j) * ( VWINDY(i,j+1) - VWINDY(i,j-1) ) / ( DYN(i,j) + DYN(i,j+1) ) );\
    DV_adv(i,j) = adv_Pv1 + adv_Pv2;                                            \
}
//End of V
