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

// ACCESS MACROS FOR AUXILIARY MATRICES USED IN THE SOLVER

#define Dtracer(i,j) ac(Dtracer_n,i,j,SW->M)
#define Deta(i,j) ac(Deta_n,i,j,SW->M)

// Du at nM2, nM1 and n, and the contributions of the pressure and the advection at the current time step n and the output of the Adams--Bashforth
#define DU_nAB(i,j) ac(Du_nAB,i,j,SW->M)
#define DU_pres(i,j) ac(Du_pres,i,j,SW->M)
#define DU_adv(i,j) ac(Du_adv,i,j,SW->M)

// Dv at nM2, nM1 and n, and the contributions of the pressure and the advection at the current time step n and the output of the Adams--Bashforth
#define DV_nAB(i,j)  ac(Dv_nAB,i,j,SW->M)
#define DV_pres(i,j) ac(Dv_pres,i,j,SW->M)
#define DV_adv(i,j)  ac(Dv_adv,i,j,SW->M)

#define funPsi(r) (MAX(0,MIN(MIN((1+r)/2.0,2.0*r),2.0)))

// ====================================
// efficient multiple schemes in runtime
// TVD scheme related macros
#define B2BP(B1,B2,rP,CP) (B1 + funPsi(rP) * 0.5 * ( 1 - CP ) * ( B2 - B1 ))
#define B2BM(B1,B2,rM,CM) (B2 - funPsi(rM) * 0.5 * ( 1 + CM ) * ( B2 - B1 ))
#define r2rP(Bprev,Bcurr,Bpost) (((Bpost - Bcurr) != 0) ? ((Bcurr - Bprev) / (Bpost - Bcurr)) : 0)
#define r2rM(Bcurr,Bpost,Bnext) (((Bpost - Bcurr) != 0) ? ((Bnext - Bpost) / (Bpost - Bcurr)) : 0)


// generic B2B macros
#define B2BPG(G,B1,B2,rP,CP) (B1 + G(rP) * 0.5 * ( 1 - CP ) * ( B2 - B1 ))
#define B2BMG(G,B1,B2,rM,CM) (B2 - G(rM) * 0.5 * ( 1 + CM ) * ( B2 - B1 ))


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
    tracer_Pv1  = CwP*BwP+CwM*BwM-CeP*BeP-CeM*BeM+CsP*BsP+CsM*BsM-CnP*BnP-CnM*BnM; \
    tracer_Pv2  = Dt * TRACER(i,j) * ( UWINDX(i,j) - UWINDX(i-1,j) ) / DXE(i,j);\
    tracer_Pv2 += Dt * TRACER(i,j) * ( VWINDY(i,j) - VWINDY(i,j-1) ) / DYN(i,j);\
    Dtracer(i,j) = tracer_Pv1 + tracer_Pv2;                                     \
                                                                                \
    if (SW->polar!=0) {                                                         \
        dist = get_great_angle(MEAN(SW->lon0,SW->lon1), MEAN(SW->lat0,SW->lat1), LONC(i,j), LATC(i,j))-SW->beta_sponge; \
        if(dist > 0){                                                           \
          sigma_dis = SW->sigma_sponge*(1.0-exp(-SW->tau_sponge*(dist)));       \
        }else{                                                                  \
          sigma_dis = SW->tracer_dis;                                           \
        }                                                                       \
    } else {                                                                    \
        sigma_dis = SW->tracer_dis;                                             \
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
    BwP = B2BPG(sch,UWINDX(i-1,j),UWINDX(i,j),rxPprev,CwP);                     \
    BwM = B2BMG(sch,UWINDX(i-1,j),UWINDX(i,j),rxMprev,CwM);                     \
    BeP = B2BPG(sch,UWINDX(i,j),UWINDX(i+1,j),rxP,CeP);                         \
    BeM = B2BMG(sch,UWINDX(i,j),UWINDX(i+1,j),rxM,CeM);                         \
    BsP = B2BPG(sch,UWINDX(i,j-1),UWINDX(i,j),ryPprev,CsP);                     \
    BsM = B2BMG(sch,UWINDX(i,j-1),UWINDX(i,j),ryMprev,CsM);                     \
    BnP = B2BPG(sch,UWINDX(i,j),UWINDX(i,j+1),ryP,CnP);                         \
    BnM = B2BMG(sch,UWINDX(i,j),UWINDX(i,j+1),ryM,CnM);                         \
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
    BwP = B2BPG(sch,VWINDY(i-1,j),VWINDY(i,j),rxPprev,CwP);                     \
    BwM = B2BMG(sch,VWINDY(i-1,j),VWINDY(i,j),rxMprev,CwM);                     \
    BeP = B2BPG(sch,VWINDY(i,j),VWINDY(i+1,j),rxP,CeP);                         \
    BeM = B2BMG(sch,VWINDY(i,j),VWINDY(i+1,j),rxM,CeM);                         \
    BsP = B2BPG(sch,VWINDY(i,j-1),VWINDY(i,j),ryPprev,CsP);                     \
    BsM = B2BMG(sch,VWINDY(i,j-1),VWINDY(i,j),ryMprev,CsM);                     \
    BnP = B2BPG(sch,VWINDY(i,j),VWINDY(i,j+1),ryP,CnP);                         \
    BnM = B2BMG(sch,VWINDY(i,j),VWINDY(i,j+1),ryM,CnM);                         \
                                                                                \
    adv_Pv1 = CwP * BwP + CwM * BwM - CeP * BeP - CeM * BeM + CsP * BsP + CsM * BsM - CnP * BnP - CnM * BnM; \
    adv_Pv2 = Dt * ( VWINDY(i,j) * ( ( UWINDX(i,j) + UWINDX(i,j+1) ) / 2. - ( UWINDX(i-1,j) + UWINDX(i-1,j+1) ) / 2. ) / DXV(i,j) + VWINDY(i,j) * ( VWINDY(i,j+1) - VWINDY(i,j-1) ) / ( DYN(i,j) + DYN(i,j+1) ) );\
    DV_adv(i,j) = adv_Pv1 + adv_Pv2;                                            \
}
//End of V


void solve_tracer(double *Dtracer_n, double Dt, sw *SW) {

/*

The tracer variable satisfies Dtracer/Dt = 0, where Dtracer/Dt is the total derivative. The total derivative includes the advection term, which is what is solved here according to the TVD scheme. Dtracer_n is the increment of tracer due to advection.

- double *Dtracer_n (input/output) The increment to be added to the tracer field at the current timestep.
- double Dt (input) The current increment of time.
- sw *SW (input) The Shallow World.

*/


    double CwP, CwM, CeP, CeM, CsP, CsM, CnP, CnM; // Courant numbers
    double rxPprev, rxMprev, rxP, rxM, ryPprev, ryMprev, ryP, ryM;
    double BwP, BwM, BeP, BeM, BsP, BsM, BnP, BnM;
    double tracer_Pv1, tracer_Pv2;
    double sigma_dis, dist;

    int i,j;

    switch(SW->sch) {
        case 0: BODYTRACER(funPsi_MUSCL)      ;  break;
        case 1: BODYTRACER(funPsi_SUPERBEE)   ;  break;
        case 2: BODYTRACER(funPsi_MINMOD)     ;  break;
        case 3: BODYTRACER(funPsi_VANLEER)    ;  break;
        case 4: BODYTRACER(funPsi_KOREN)      ;  break;
        case 5: BODYTRACER(funPsi_OSPRE)      ;  break;
        case 6: BODYTRACER(funPsi_VANALBADA1) ;  break;
        case 7: BODYTRACER(funPsi_UMIST)      ;  break;
        case 8: BODYTRACER(funPsi_UPWIND)     ;  break;
        case 9: BODYTRACER(funPsi_LAXWENDROFF);  break;
        default: CRASH("scheme=%d a la parra, use select_scheme first !!",SW->sch);
    }
}

/*

Computes the increment of the surface perturbation Deta_n according to the equation of the conservation of h and also to the Superbee TVD scheme. Deta_n is the increment of eta due to advection.

See page 21 of https://www.dropbox.com/s/2p0t946nw0h9p88/thesis_Prat_document.pdf?dl=0

- double *Deta_n (input/output) The increment of the surface perturbation or eta.
- double Dt (input) The current increment of time.
- sw *SW (input) The Shallow World.

*/

void solve_cons_h(double *Deta_n, double Dt, sw *SW) {

    int info=1;

    double CwP, CwM, CeP, CeM, CsP, CsM, CnP, CnM; // Courant numbers
    double rxPprev, rxMprev, rxP, rxM, ryPprev, ryMprev, ryP, ryM;
    double BwP, BwM, BeP, BeM, BsP, BsM, BnP, BnM;


    double sigma_dis, dist;

    double *Le = dmem(SW->M);
    double *LLe = dmem(SW->M);
    double *LLLe = dmem(SW->M);


    int i,j;


    switch(SW->sch) {
        case 0: BODYH(funPsi_MUSCL)      ;  break;
        case 1: BODYH(funPsi_SUPERBEE)   ;  break;
        case 2: BODYH(funPsi_MINMOD)     ;  break;
        case 3: BODYH(funPsi_VANLEER)    ;  break;
        case 4: BODYH(funPsi_KOREN)      ;  break;
        case 5: BODYH(funPsi_OSPRE)      ;  break;
        case 6: BODYH(funPsi_VANALBADA1) ;  break;
        case 7: BODYH(funPsi_UMIST)      ;  break;
        case 8: BODYH(funPsi_UPWIND)     ;  break;
        case 9: BODYH(funPsi_LAXWENDROFF);  break;
        default: CRASH("scheme=%d a la parra, use select_scheme first !!",SW->sch);
    }
    // BODYH(funPsi_MUSCL);

    // Application of the Hyperviscosity
    if ( SW->nu2 !=0 || SW->nu4 !=0 || SW->nu6 !=0) {

        setzero_scaf(Le,SW->M);
        setzero_scaf(LLe,SW->M);
        setzero_scaf(LLLe,SW->M);


        if (SW->prn) pprintf("Hypervisosity computing lap2 eta\n");
        laplacian_eta(SW->u,Le,SW);


        if (SW->nu4 !=0 || SW->nu6 !=0) {
            if (SW->prn) pprintf("Hypervisosity computing lap4 eta\n");
            laplacian_eta(Le,LLe,SW);


            if (SW->nu6!=0) {
                if (SW->prn) pprintf("Hypervisosity computing lap6 eta\n");
                laplacian_eta(LLe,LLLe,SW);

            }
        }
        forall(i,j,SW->M){
          ac(Deta_n,i,j,SW->M) += (SW->nu2)*(SW->Dt)*ac(Le,i,j,SW->M);
          ac(Deta_n,i,j,SW->M) += (SW->nu4)*(SW->Dt)*ac(LLe,i,j,SW->M);
          ac(Deta_n,i,j,SW->M) += (SW->nu6)*(SW->Dt)*ac(LLLe,i,j,SW->M);
        }
    }
    free(Le);
    free(LLe);
    free(LLLe);
}

void dissipation_effect(double *Du_n, double *Dv_n, sw *SW){
  /*
    This functions adds the effect of the dissipation for either a sponge boundary or
    the one for the whole domain.

    -double *Du_n (input/output) The increment of the horizontal velocity
    -double *Dv_n (input/output) The increment of the vertical velocity
    - sw *SW (input/output) The Shallow World.
  */


  double distU, sigma_disU, distV, sigma_disV;

  int i,j;

  forall(i,j,SW->M){

    if (SW->polar==0) { // if not polar ...
        sigma_disU=SW->sigma_dis;
        sigma_disV=SW->sigma_dis;

    } else { // for polar ..

        distU = get_great_angle(MEAN(SW->lon0,SW->lon1), MEAN(SW->lat0,SW->lat1), LONE(i,j), LATE(i,j))-SW->beta_sponge;
        distV = get_great_angle(MEAN(SW->lon0,SW->lon1), MEAN(SW->lat0,SW->lat1), LONN(i,j), LATN(i,j))-SW->beta_sponge;

        if(distU > 0){
          sigma_disU = SW->sigma_sponge*(1.0-exp(-SW->tau_sponge*(distU)));
        }else{
          sigma_disU = SW->sigma_dis;
        }

        if(distV > 0){
          sigma_disV = SW->sigma_sponge*(1.0-exp(-SW->tau_sponge*(distV)));
        }else{
          sigma_disV = SW->sigma_dis;
        }

    }
    ac(Du_n,i,j,SW->M) -= U(i,j)*sigma_disU*SW->Dt;
    ac(Dv_n,i,j,SW->M) -= V(i,j)*sigma_disV*SW->Dt;
  }

}


void do_eta_geofactor(double *Deta_nAB, double Dt, sw *SW) {

/*

Modifies the current increment of eta Deta_nAB (after Adams--Bashforth) to take into account the effect of the ellipsoidal shape of the planet.

See transparency 31 of https://www.dropbox.com/s/1d5n1fzsgz7r1e8/thesis_Prat_transparencies.pdf?dl=0

- double *Deta_nAB (input/output) The increment of the surface perturbation after Adams--Bashforth.
- double Dt (input) The current increment of time.
- sw *SW (input/output) The Shallow World.

*/

    int i,j;

    forall(i,j,SW->M) {
        ac(Deta_nAB,i,j,SW->M)-=(Dt*sin(LATC(i,j))/rZ(LATC(i,j))*H(i,j)*0.5*(VWINDY(i,j)+VWINDY(i,j-1)));
    }
}

void solve_adv_u(double *Du_adv, double Dt, sw *SW) {

/*

Computes the increment of the horizontal velocities Du_adv due to advection.

See page 23 of https://www.dropbox.com/s/2p0t946nw0h9p88/thesis_Prat_document.pdf?dl=0

- double *Du_adv (input/output) The increment of the horizontal velocities due to advection.
- double Dt (input) The current increment of time.
- sw *SW (input) The Shallow World.

*/

    double CwP, CwM, CeP, CeM, CsP, CsM, CnP, CnM; // Courant numbers
    double rxPprev, rxMprev, rxP, rxM, ryPprev, ryMprev, ryP, ryM;
    double BwP, BwM, BeP, BeM, BsP, BsM, BnP, BnM;
    double adv_Pu1, adv_Pu2; // Results

    int i,j;

    switch(SW->sch) {
        case 0: BODYU(funPsi_MUSCL)      ;  break;
        case 1: BODYU(funPsi_SUPERBEE)   ;  break;
        case 2: BODYU(funPsi_MINMOD)     ;  break;
        case 3: BODYU(funPsi_VANLEER)    ;  break;
        case 4: BODYU(funPsi_KOREN)      ;  break;
        case 5: BODYU(funPsi_OSPRE)      ;  break;
        case 6: BODYU(funPsi_VANALBADA1) ;  break;
        case 7: BODYU(funPsi_UMIST)      ;  break;
        case 8: BODYU(funPsi_UPWIND)     ;  break;
        case 9: BODYU(funPsi_LAXWENDROFF);  break;
        default: CRASH("scheme=%d a la parra, use select_scheme first !!",SW->sch);
    }


}

void solve_adv_v(double *Dv_adv, double Dt, sw *SW) {

/*

Computes the increment of the vertical velocities Dv_adv due to advection.

See page 26 of https://www.dropbox.com/s/2p0t946nw0h9p88/thesis_Prat_document.pdf?dl=0

- double *Dv_adv (input/output) The increment of the vertical velocities due to advection.
- double Dt (input) The current increment of time.
- sw *SW (input) The Shallow World.

*/

    double CwP, CwM, CeP, CeM, CsP, CsM, CnP, CnM; // Courant numbers
    double rxPprev, rxMprev, rxP, rxM, ryPprev, ryMprev, ryP, ryM;
    double BwP, BwM, BeP, BeM, BsP, BsM, BnP, BnM;
    double adv_Pv1, adv_Pv2; // Results

    int i,j;

    switch(SW->sch) {
        case 0: BODYV(funPsi_MUSCL)      ;  break;
        case 1: BODYV(funPsi_SUPERBEE)   ;  break;
        case 2: BODYV(funPsi_MINMOD)     ;  break;
        case 3: BODYV(funPsi_VANLEER)    ;  break;
        case 4: BODYV(funPsi_KOREN)      ;  break;
        case 5: BODYV(funPsi_OSPRE)      ;  break;
        case 6: BODYV(funPsi_VANALBADA1) ;  break;
        case 7: BODYV(funPsi_UMIST)      ;  break;
        case 8: BODYV(funPsi_UPWIND)     ;  break;
        case 9: BODYV(funPsi_LAXWENDROFF);  break;
        default: CRASH("scheme=%d a la parra, use select_scheme first !!",SW->sch);
    }


}

void solve_pres_u(double *Du_pres, double Dt, sw *SW) {

/*

Computes the increment of the horizontal velocities Du_pres due to the gradient of pressure in the x axis.

See page 28 of https://www.dropbox.com/s/2p0t946nw0h9p88/thesis_Prat_document.pdf?dl=0

- double *Du_pres (output) The increment of the horizontal velocities due to the gradient of pressure.
- double Dt (input) The current increment of time.
- sw *SW (input) The Shallow World.

*/

    int i, j;

    forall(i,j,SW->M)
        DU_pres(i,j) = -GE(i,j)*Dt*(ETA(i+1,j)-ETA(i,j))/DXC(i+1,j);

}

void solve_pres_v(double *Dv_pres, double Dt, sw *SW) {

/*

Computes the increment of the vertical velocities Du_pres due to the gradient of pressure in the y axis.

See page 28 of https://www.dropbox.com/s/2p0t946nw0h9p88/thesis_Prat_document.pdf?dl=0

- double *Du_pres (output) The increment of the vertical velocities due to the gradient of pressure.
- double Dt (input) The current increment of time.
- sw *SW (input) The Shallow World.

*/

    int i, j;

    forall(i,j,SW->M)
        DV_pres(i,j) = -GN(i,j)*Dt*(ETA(i,j+1)-ETA(i,j))/DYC(i,j+1);

}

void solve_Coriolis_u(double *Du_nAB, double Dt, sw *SW) {

/*

Computes the contribution to the horizontal velocities of the Coriolis effect.

See page 29 of https://www.dropbox.com/s/2p0t946nw0h9p88/thesis_Prat_document.pdf?dl=0

- double *Du_nAB (input) The increment of the horizontal velocities after calling do_AdamsBashforth.
- double Dt (input) The current increment of time.
- sw *SW (input/output) The Shallow World.

The term arising from expressing the conservation of momentum in ellipsoidal coordinates is integrated here (only what affects the horizontal velocities).

*/

    int i,j;
    double Alpha, Beta, semiU;

    forall(i,j,SW->M) {
        semiU=0.25*(V(i,j)+V(i+1,j)+V(i,j-1)+V(i+1,j-1));
        Alpha=Dt*FE(i,j) + (Dt*U(i,j)*sin(LATC(i,j)))/rZ(LATE(i,j)); // The second summand arises from the ellipsoidal shape
        Beta=0.25*Alpha*Alpha; // 0.25 is the semi implicit factor squared

        U(i,j)=(U(i,j)+DU_nAB(i,j)-Beta*U(i,j)+Alpha*semiU)/(1.0+Beta);
    }

}

void solve_Coriolis_v(double *Dv_nAB, double Dt, sw *SW) {

/*

Computes the contribution to the vertical velocities of the Coriolis effect.

See page 29 of https://www.dropbox.com/s/2p0t946nw0h9p88/thesis_Prat_document.pdf?dl=0

- double *Dv_nAB (input) The increment of the vertical velocities after calling do_AdamsBashforth.
- double Dt (input) The current increment of time.
- sw *SW (input/output) The Shallow World.

The term arising from expressing the conservation of momentum in ellipsoidal coordinates is integrated here (only what affects the vertical velocities).

*/

    int i,j;
    double Alpha, Beta, semiV;

    forall(i,j,SW->M) {
        semiV=0.25*(U(i-1,j+1)+U(i,j+1)+U(i,j)+U(i-1,j));
        Alpha=Dt*FN(i,j) + (Dt*semiV*sin(LATN(i,j)))/rZ(LATN(i,j)); // The second summand arises from the ellipsoidal shape
        Beta=0.25*Alpha*Alpha; // 0.25 is the semi implicit factor squared

        V(i,j)=(V(i,j)+DV_nAB(i,j)-Beta*V(i,j)-Alpha*semiV)/(1.0+Beta);
    }

}

void do_channel(sw *SW) {

/*

Sets the boundary conditions for a wall at the top and at the bottom of the domain. This function modifies the u, v, windx, windy of SW.

See page 30 of https://www.dropbox.com/s/2p0t946nw0h9p88/thesis_Prat_document.pdf?dl=0

- sw *SW (input/output) The Shallow World.

*/

    int i,j;

    forsouth(i,j,SW->M) { // do j (y axis)

        V(i,j-1) = 0; // this has to be the velocity of the botttom boundary and, therefore, in the halo
        // fill halo for safety and numerical routines
        V(i,j-2) = 0;

        WINDY(i,j-1) = 0;
        WINDY(i,j-2) = 0;

        U(i,j-1) = U(i,j);
        U(i,j-2) = U(i,j-1);

        WINDX(i,j-1) = WINDX(i,j);
        WINDX(i,j-2) = WINDX(i,j-1);

    }

    fornorth(i,j,SW->M) { // do j (y axis)

        V(i,j) = 0; // this has to be the velocity of the top boundary and, therefore, not in the halo
        // fill halo for safety and numerical routines
        V(i,j+1) = 0;
        V(i,j+2) = 0;

        WINDY(i,j) = 0;
        WINDY(i,j+1) = 0;
        WINDY(i,j+2) = 0;

        U(i,j+1) = U(i,j);
        U(i,j+2) = U(i,j+1);

        WINDX(i,j+1) = WINDX(i,j);
        WINDX(i,j+2) = WINDX(i,j+1);
    }

}

void do_AdamsBashforth(double *Dvar_nAB, double *Dvar_n, double *Dvar_nM1, double *Dvar_nM2, int n, map *M) {

/*

Computes Dvar_nAB, which is the increment to be added to the current variable according to Adams--Bashforth expressions. These require the two past results of Dvar_n, which are identified as Dvar_nM1 and Dvar_nM2.

See page 28 of https://www.dropbox.com/s/2p0t946nw0h9p88/thesis_Prat_document.pdf?dl=0

- double *Dvar_nAB (output) The real increment of the variable.
- double *Dvar_n (input) The current increment of the variable.
- double *Dvar_nM1 (input) The past increment of the variable.
- double *Dvar_nM2 (input) The increment of the variable used two time steps before.
- map *M (input) The sppde's map.

*/
    /*
    int i,j;
    if(n<1)
        forall(i,j,M)
            ac(Dvar_nAB,i,j,M) = ac(Dvar_n,i,j,M);
    else if(n<2)
        forall(i,j,M)
            ac(Dvar_nAB,i,j,M) = (3.0 / 2.0) * ac(Dvar_n,i,j,M) - (1.0 / 2.0) * ac(Dvar_nM1,i,j,M);
    else
        forall(i,j,M)
            ac(Dvar_nAB,i,j,M) = (23.0 / 12.0) * ac(Dvar_n,i,j,M) - (4.0 / 3.0) * ac(Dvar_nM1,i,j,M) + (5.0 / 12.0) * ac(Dvar_nM2,i,j,M);

    */

    if(n == 0) {
        // Dvar_nAB=1.0*Dvar_n+0.0*Dvar_n
        // linop_scaf(Dvar_nAB,Dvar_n,Dvar_n,0.0,1.0,M);
        copy_scaf(Dvar_nAB, Dvar_n, 0, M);
    } else if(n == 1) {
        // Dvar_nAB=(3.0 / 2.0)*Dvar_n- (1.0 / 2.0) *Dvar_n
        linop_scaf(Dvar_nAB,Dvar_n,Dvar_nM1,(3.0 / 2.0), - (1.0 / 2.0),M);
    } else {
        // We want to compute Dvar_nAB=(23.0 / 12.0)*Dvar_n - (4.0 / 3.0)*Dvar_nM1 + (5.0 / 12.0)) *Dvar_nM2. We will split the sum into two parts
        // Dvar_nAB=- (4.0 / 3.0)*Dvar_nM1 + (5.0 / 12.0)) *Dvar_nM2
        linop_scaf(Dvar_nAB,Dvar_nM1,Dvar_nM2, - (4.0 / 3.0), + (5.0 / 12.0),M);
            // Update: Dvar_nAB=(23.0 / 12.0)*Dvar_n + 1.0 * Dvar_nAB
        linop_scaf(Dvar_nAB,Dvar_n,Dvar_nAB, (23.0 / 12.0), 1,M); // Asumes matrices can be updated in linop_scaf

        // The above is equivalent to
        // forall(i,j,M)
            // ac(Dvar_nAB,i,j,M) = (23.0 / 12.0) * ac(Dvar_n,i,j,M) - (4.0 / 3.0) * ac(Dvar_nM1,i,j,M) + (5.0 / 12.0) * ac(Dvar_nM2,i,j,M);
    }

}

double compute_Dt(double Courant, sw *SW) {

/*

Returns the required Delta t to be added to the current time t to perform the simulation. From the definition of the Courant number in two dimensions, the function returns the smallest Dt of the domain. For velocity fields of 0, the functions computes Dt according to the maximum possible velocity, given by u_max, v_max = sqrt(gH).

- double Courant (input) The maximum allowed Courant number. Smaller courant numbers are more restrictive.
- sw *SW (input) The Shallow World.

*/

    double max_x = 0, max_y = 0;

    int i,j;

    forall(i,j,SW->M) {
        if(fabs(UWINDX(i,j)/DXE(i,j))>max_x) max_x = fabs(UWINDX(i,j))/DXE(i,j); // For dim 1
        if(fabs(VWINDY(i,j)/DYN(i,j))>max_y) max_y = fabs(VWINDY(i,j))/DYN(i,j); // For dim 2
    }

    if (max_x + max_y == 0) {

        forall(i,j,SW->M) {

            if(fabs(sqrt((SW->depth)*(GE(i,j)))/DXE(i,j))>max_x) max_x = fabs(sqrt((SW->depth)*(GE(i,j))))/DXE(i,j); // For dim 1

            if(fabs(sqrt((SW->depth)*(GN(i,j)))/DYN(i,j))>max_y) max_y = fabs(sqrt((SW->depth)*(GN(i,j))))/DYN(i,j); // For dim 2
        }
    }

    double Dt_CFL;
    double local_Dt_CFL = (Courant/(max_x+max_y));

    checkr(MPI_Allreduce(&(local_Dt_CFL),&Dt_CFL,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD),"Allreduce"); // Get the minimum Dt of all processes

    return Dt_CFL;

}

//  Geostrophic equilibrium


void init_scaf_eta(double *eta, double n, double amp, double loncent, double latcent, double as, double bs, sw *SW){
  int i, j;
  double exponent;
  int info=0;
  forall(i,j,SW->M){
    exponent = pow(POW2((rad2deg(LONC(i,j))-loncent)/as) + POW2((rad2deg(LATC(i,j))-latcent)/bs),n);
    ETA(i,j) += amp*exp(-exponent);
  }
  if(info) print_scaf(eta, "ETA init", 0, SW->M,"%3e ");
}

void init_scaf_u(double *u, double *eta, double n,double loncent, double latcent, double as, double bs, sw *SW){
  int i, j;
  int info=0;
  double psi;
  double R = (SW->rE+SW->rP)/2.0;
  forall(i,j,SW->M){
    psi = POW2((rad2deg(LONE(i,j))-loncent)/as) + POW2((rad2deg(LATE(i,j))-latcent)/bs);
    U(i,j) += 2.0*n*GE(i,j)/(FE(i,j)*POW2(deg2rad(bs))*R)*(LATE(i,j)-deg2rad(latcent))*pow(psi,n-1.0)*ETA(i,j);
  }
  if (info) print_scaf(u, "U init", 0, SW->M,"%3e ");
}

void init_scaf_v(double *v,  double *eta,double n,double loncent, double latcent, double as, double bs, sw *SW){
  int i, j;
  int info=0;
  double psi;
  double R = (SW->rE+SW->rP)/2;
  forall(i,j,SW->M){
    psi = POW2((rad2deg(LONN(i,j))-loncent)/as) + POW2((rad2deg(LATN(i,j))-latcent)/bs);
    V(i,j) += -2.0*(n)*GN(i,j)/(FN(i,j)*R*cos(LATN(i,j))*POW2(deg2rad(as)))*(LONN(i,j)-deg2rad(loncent))*pow(psi,n-1.0)*ETA(i,j);
  }
  if (info) print_scaf(v, "V init", 0, SW->M,"%3e ");
}
