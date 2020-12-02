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
#include <libgen.h>
#include "mpi.h"

#include "sppde.h"
#include "sw.h"

void loop(char *bname, char *output_folder, sw *SW) {

/*

Runs the time iterations, a.k.a. performs the operations to output the solution. In fact this is the core of the program.

- sw *SW (input/output) The Shallow World.

*/

    int info = 0;

    int i, j;

    // Auxiliary matrices
    double *Dtracer_nAB, *Dtracer_n, *Dtracer_nM1, *Dtracer_nM2;
    double *Deta_nAB, *Deta_n, *Deta_nM1, *Deta_nM2;
    double *Du_nAB, *Du_adv, *Du_pres, *Du_n, *Du_nM1, *Du_nM2;
    double *Dv_nAB, *Dv_adv, *Dv_pres, *Dv_n, *Dv_nM1, *Dv_nM2;
    double *Lu, *Lv, *LLu, *LLv, *LLLu, *LLLv; //Hyperviscosity matrices MARC

    // Simulation parameters
    int n;// Number of time steps
    double t0 = SW->t0, t = t0;
    double Dt = SW->Dt;
    double local_Dt_CFL, Dt_CFL;



    // energy
    double ke,ape;

    // Auxiliary matrices allocation
    Dtracer_nAB = dmem(SW->M);
    Dtracer_n = dmem(SW->M);
    Dtracer_nM1 = dmem(SW->M);
    Dtracer_nM2 = dmem(SW->M);
    Deta_nAB = dmem(SW->M);
    Deta_n = dmem(SW->M);
    Deta_nM1 = dmem(SW->M);
    Deta_nM2 = dmem(SW->M);
    Du_nAB = dmem(SW->M);
    Du_adv = dmem(SW->M);
    Du_pres = dmem(SW->M);
    Du_n = dmem(SW->M);
    Du_nM1 = dmem(SW->M);
    Du_nM2 = dmem(SW->M);
    Dv_nAB = dmem(SW->M);
    Dv_adv = dmem(SW->M);
    Dv_pres = dmem(SW->M);
    Dv_n = dmem(SW->M);
    Dv_nM1 = dmem(SW->M);
    Dv_nM2 = dmem(SW->M);
    // hyperv
    Lu = dmem(SW->M);
    Lv = dmem(SW->M);
    LLu = dmem(SW->M);
    LLv = dmem(SW->M);
    LLLu = dmem(SW->M);
    LLLv = dmem(SW->M);



    pprintf("Auxiliary matrices allocated.\n");

    // Update halo for safety
    halo_update(SW->tracer,SW->M);
    halo_update(SW->u,SW->M);
    halo_update(SW->v,SW->M);
    halo_update(SW->eta,SW->M);
    halo_update(SW->hB,SW->M);

    setzero_scaf(SW->perturbs,SW->M);  // Set perturbs field to 0
    add_Gaussians(Dt, t, SW); // The time to be added here is the updated time t0+dt not t0

    halo_update(SW->tracer,SW->M); // Update the halos
    halo_update(SW->eta,SW->M); // Update the halos

    do_channel(SW);

    pprintf("Halos updated.\n");

    cr_start("sim_time",0);

    stats_scaf(SW->tracer,SW->M,"tracer"); // Print field information
    stats_scaf(SW->eta,SW->M,"eta"); // Print field information
    stats_scaf(SW->u,SW->M,"u"); // Print field information
    stats_scaf(SW->v,SW->M,"v"); // Print field information

    int frame = 0; // Start frame indexing from 1
    int nstart;

    if(SW->loadFrom!=0) { // If file number is not 0, it means that a file has to be loaded
        frame = SW->loadFrom;
        int error_load = load_binary_sw(frame,&nstart,&t,bname,output_folder,SW); // All processors have the variables nstart, t updated
        if(error_load) CRASH("File to be loaded not found");
    } else {
        nstart = 0;
        // Save initial conditions as frame 0
        save_binary_sw(frame,n,t,bname, output_folder,SW);
    }

    // calc and print initial energy
    calc_energy_balance(&ke,&ape,SW);
    print_energy(t,ke,ape);

    // Update frame number
    frame++; // Number of files currently saved minus 1 (frame starts at 0)

    // Update nstart number
    nstart++;
    for(n=nstart;n<=SW->N;n++) {

        SW->prn= (n==nstart || n % (SW->IteInfo) == 0);

        //pprintf ("CORE n=%d SW->prn=%d\n",n,SW->prn);
        /*
        Under n=1, the operations carried out are those that happen between the initial conditions (n=0) and the first step(n=1). In other words, when this for starts, n=1.
        */

        cr_start("time_step",0);
        if (SW->prn) pprintf("===SW STARTING TIME STEP %d / %d from t0=%e to t=%e \n",n,SW->N,t,t+Dt);
        t += Dt; // we are calculating instant 't', ie the first we calculate is t0+dt

        /*
        All the operations below are computed for t+Dt, which follow the scheme below:
        t=t0 is step n=0 (initial conditions saved as frame=0) and is done outside the for loop
        t=t+1*Dt is step n=1 (conditions for the first step are saved at the end of the step with these t and n)
        t=t+2*Dt is step n=2 (conditions for the second step are saved at the end of the step with these t and n)
        ...
        */

        cr_start("compute_Dt",0);
        Dt_CFL = compute_Dt(SW->Courant_max, SW); // Depends on the grid, as well as SW->u and SW->v
        if (info || SW->prn) pprintf("Smallest possible timestep is Dt_CFL=%.6e while Dt=%.6e.\n",Dt_CFL,Dt);
        if (Dt>Dt_CFL) CRASH("Smallest possible timestep is Dt_CFL=%.6e while Dt=%.6e.\n",Dt_CFL,Dt);
        cr_end("compute_Dt",0);

        cr_start("solve_tracer",0);
        solve_tracer(Dtracer_n, Dt, SW);
        cr_end("solve_tracer",0);

        cr_start("AB_Dtracer",0);
        do_AdamsBashforth(Dtracer_nAB, Dtracer_n, Dtracer_nM1, Dtracer_nM2, n-nstart, SW->M);
        cr_end("AB_Dtracer",0);

        cr_start("copy_Dtracer",0);
        copy_scaf(Dtracer_nM2,Dtracer_nM1,0,SW->M);
        copy_scaf(Dtracer_nM1,Dtracer_n,0,SW->M);
        cr_end("copy_Dtracer",0);

        cr_start("new_tracer",0);
        forall(i, j, SW->M) {
            TRACER(i,j)=TRACER(i,j)+ac(Dtracer_nAB,i,j,SW->M); // Add increment
        }
        cr_end("new_tracer",0);

        cr_start("halo_update",0);
        halo_update(SW->tracer,SW->M); // Update the halos
        cr_end("halo_update",0);

        cr_start("solve_cons_h",0);
        solve_cons_h(Deta_n, Dt, SW); // Depends on the grid, as well as Dt, u and v
        cr_end("solve_cons_h",0);

        cr_start("AB_Deta",0);
        do_AdamsBashforth(Deta_nAB, Deta_n, Deta_nM1, Deta_nM2, n-nstart, SW->M);
        cr_end("AB_Deta",0);

        cr_start("eta_geof",0);
        do_eta_geofactor(Deta_nAB, Dt, SW);
        cr_end("eta_geof",0);

        // Update/shift Deta
        cr_start("copy_Deta",0);
        copy_scaf(Deta_nM2,Deta_nM1,0,SW->M);
        copy_scaf(Deta_nM1,Deta_n,0,SW->M);
        cr_end("copy_Deta",0);

        cr_start("new_eta",0);
        forall(i, j, SW->M) {
            ETA(i,j)=ETA(i,j)+ac(Deta_nAB,i,j,SW->M); // Add increment
        }
        cr_end("new_eta",0);

        cr_start("halo_update",0);
        halo_update(SW->eta,SW->M); // Update the halos
        cr_end("halo_update",0);

        cr_start("solve_adv",0);
        solve_adv_u(Du_adv, Dt, SW); // Compute Du_adv from the advection in the x axis by means of the Superbee TVD scheme
        solve_adv_v(Dv_adv, Dt, SW); // Compute Dv_adv from the advection in the y axis by means of the Superbee TVD scheme
        cr_end("solve_adv",0);

        cr_start("solve_press",0);
        solve_pres_u(Du_pres, Dt, SW); // Compute Du_pres from the gradient of pressure in the x axis
        solve_pres_v(Dv_pres, Dt, SW); // Compute Dv_pres from the gradient of pressure in the y axis
        cr_end("solve_press",0);

        cr_start("new_vels",0);
        forall(i, j, SW->M) {
            ac(Du_n,i,j,SW->M)=ac(Du_adv,i,j,SW->M)+ac(Du_pres,i,j,SW->M);// Add the contributions of advection and pressure
            ac(Dv_n,i,j,SW->M)=ac(Dv_adv,i,j,SW->M)+ac(Dv_pres,i,j,SW->M);// Add the contributions of advection and pressure
        }
        cr_end("new_vels",0);

        cr_start("dis_vel",0);
        dissipation_effect(Du_n,Dv_n,SW);
        cr_end("dis_vel",0);

        cr_start("hyper_visc",0);
        if ( SW->nu2 !=0 || SW->nu4 !=0 || SW->nu6 !=0) {

            setzero_scaf(Lu,SW->M);
            setzero_scaf(Lv,SW->M);
            setzero_scaf(LLu,SW->M);
            setzero_scaf(LLv,SW->M);
            setzero_scaf(LLLu,SW->M);
            setzero_scaf(LLLv,SW->M);

            if (SW->prn) pprintf("Hypervisosity computing lap2 uv\n");
            laplacian_vels(SW->u, SW->v, Lu, Lv, SW);


            if (SW->nu4 !=0 || SW->nu6 !=0) {
                if (SW->prn) pprintf("Hypervisosity computing lap4 uv\n");
                laplacian_vels(Lu, Lv, LLu, LLv, SW);


                if (SW->nu6!=0) {
                    if (SW->prn) pprintf("Hypervisosity computing lap6 uv\n");
                    laplacian_vels(LLu, LLv, LLLu, LLLv, SW);

                }
            }
            forall(i,j,SW->M){
              ac(Du_n,i,j,SW->M) += (SW->nu2)*(SW->Dt)*ac(Lu,i,j,SW->M);
              ac(Du_n,i,j,SW->M) += (SW->nu4)*(SW->Dt)*ac(LLu,i,j,SW->M);
              ac(Du_n,i,j,SW->M) += (SW->nu6)*(SW->Dt)*ac(LLLu,i,j,SW->M);

              ac(Dv_n,i,j,SW->M) += (SW->nu2)*(SW->Dt)*ac(Lv,i,j,SW->M);
              ac(Dv_n,i,j,SW->M) += (SW->nu4)*(SW->Dt)*ac(LLv,i,j,SW->M);
              ac(Dv_n,i,j,SW->M) += (SW->nu6)*(SW->Dt)*ac(LLLv,i,j,SW->M);
            }
        }
        cr_end("hyper_visc",0);

        cr_start("AB_vels",0);
        do_AdamsBashforth(Du_nAB, Du_n, Du_nM1, Du_nM2, n-nstart, SW->M);
        do_AdamsBashforth(Dv_nAB, Dv_n, Dv_nM1, Dv_nM2, n-nstart, SW->M);
        cr_end("AB_vels",0);

        cr_start("copy_Dvels",0);
        copy_scaf(Du_nM2,Du_nM1,0,SW->M); // Du_nM2 = Du_nM1
        copy_scaf(Du_nM1,Du_n,0,SW->M);
        copy_scaf(Dv_nM2,Dv_nM1,0,SW->M);
        copy_scaf(Dv_nM1,Dv_n,0,SW->M);
        cr_end("copy_Dvels",0);

        cr_start("new_Coriolis",0);
        solve_Coriolis_u(Du_nAB, Dt, SW);
        solve_Coriolis_v(Dv_nAB, Dt, SW);
        cr_end("new_Coriolis",0);

        cr_start("halo_update",0);
        halo_update(SW->u,SW->M); // Update the halos
        halo_update(SW->v,SW->M); // Update the halos
        cr_end("halo_update",0);

        cr_start("Perturbation",0);
        setzero_scaf(SW->perturbs,SW->M);  // Set perturbs field to 0
        add_Gaussians(Dt, t, SW); // The time to be added here is the updated time t0+dt not t0
        add_vortices(Dt, t, SW);
	      cr_end("Perturbation",0);

        cr_start("halo_update",0);
        halo_update(SW->tracer,SW->M); // Update the halos
        halo_update(SW->eta,SW->M); // Update the halos
        cr_end("halo_update",0);

        if (SW->polar==0){
          cr_start("do_channel",0);
          do_channel(SW);
          cr_end("do_channel",0);
        }
        // Print field information and stop simulation if nonsense values are detected
        check_stats_scaf(SW->tracer,SW->M,"tracer",n,-1e10,1e10,SW->prn);
        check_stats_scaf(SW->eta,SW->M,"eta",n,-1e10,1e10,SW->prn);
        check_stats_scaf(SW->u,SW->M,"u",n,-1e10,1e10,SW->prn);
        check_stats_scaf(SW->v,SW->M,"v",n,-1e10,1e10,SW->prn);

        // Save binary files
        cr_start("save_bin",0);
        if(n % SW->saveEvery == 0) {
            save_binary_sw(frame,n,t,bname, output_folder,SW);
            frame++; // Number of files currently saved minus 1 (frame starts at 0)

        }
        cr_end("save_bin",0);

        // Compare with reference fields if needed
        double delta;
        delta=compare_ref(SW,n,t,bname,output_folder);
        if (delta>=0) pprintf("SW compare_ref n=%d t=%e delta=%e\n",n,t,delta);

        if (SW->prn) {
            calc_energy_balance(&ke,&ape,SW);
            print_energy(t,ke,ape);
        }

        if (SW->prn)  // MANEL117d
            pprintf("===TIME STEP END TIME STEP %d wall clock seconds per timestep %.3f = %.3f ms; TOTAL per timestep= %.3f ms\n",
                    n,cr_time("time_step",0),cr_time("time_step",0)*1000,cr_time("time_step",0)*1000*quants());
        cr_end("time_step",0);

    }

    cr_end("sim_time",0);

    pprintf("Core loop ended.\n");

    // Free memory
    free(Dtracer_nAB);
    free(Dtracer_n);
    free(Dtracer_nM1);
    free(Dtracer_nM2);
    free(Deta_nAB);
    free(Deta_n);
    free(Deta_nM1);
    free(Deta_nM2);
    free(Du_nAB);
    free(Du_adv);
    free(Du_pres);
    free(Du_n);
    free(Du_nM1);
    free(Du_nM2);
    free(Dv_nAB);
    free(Dv_adv);
    free(Dv_pres);
    free(Dv_n);
    free(Dv_nM1);
    free(Dv_nM2);
    free(Lu);
    free(Lv);
    free(LLu);
    free(LLv);
    free(LLLu);
    free(LLLv);
    pprintf("Memory for auxiliary variables freed.\n");

}
