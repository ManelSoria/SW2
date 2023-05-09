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


void add_Gaussians(double Dt, double t, sw *SW) {

/*

Assigns to eta the value according to Gaussian functions. The flow rate Q is an input and therefore this function is called at every timestep from t=t0_Gauss to t=t1_Gauss.

- double Dt (input) The increment of time.
- double t (input) The current time.
- sw *SW (input/output) The Shallow World.

*/
    int info=0;
    int j;
    double tau_Gauss,t0_Gauss,t1_Gauss,lat_Gauss,vel_Gauss,lon_Gauss,Q_Gauss,sigma_Gauss,radius_Gauss; // File inputs
    double amp_Gauss; // Amplitude of perturbation computed

    for(j=0; j<=SW->num_Gaussians-1; j++){
    /*


    The positions of the longitudes, latitudes, velocities... in the Gaussians table are known:

    Magnitude:  t0 t1 lon lat vel Q sigma radius
    Offset:     0  1  2   3   4   5 6     7

    max(Offset) + 1, a.k.a. number of columns in the file should match N_COLS_GAUSS

    */

        // TODO_115 ARNAU->115A
        // deg=1: degrees, otherwise meters as in the initial versions
        // int deg=SW->Gaussians[j*N_COLS_GAUSS + posio]

        // 0     1     2      3        4        5             6              7         8              9            10             11
        // t0(s) t1(s) tau(s) lon(deg) lat(deg) isdegvel(0/1) vel(m/s-deg/s) Q(m^3/s)  isdegsize(0/1) sigma(m-deg) radius(m-deg)  inject_tracer(0/1)

        t0_Gauss = SW->Gaussians[j*N_COLS_GAUSS + 0];
        t1_Gauss = SW->Gaussians[j*N_COLS_GAUSS + 1];

        if((t0_Gauss<=t)&&(t1_Gauss>=t)) {

            tau_Gauss = SW->Gaussians[j*N_COLS_GAUSS + 2]; // ARNAU // Intensity factor
            if (tau_Gauss<0) CRASH("Uhhh tau is %e and should be positive",tau_Gauss);

            double mid_LON, mid_LAT;

            lon_Gauss = deg2rad(SW->Gaussians[j*N_COLS_GAUSS + 3]);
            lat_Gauss = deg2rad(SW->Gaussians[j*N_COLS_GAUSS + 4]);

            mid_LON = deg2rad(MEAN(SW->lon0,SW->lon1));
            mid_LAT = deg2rad(MEAN(SW->lat0,SW->lat1));

            int isdegvel = SW->Gaussians[j*N_COLS_GAUSS + 5]; // ARNAU
            vel_Gauss = SW->Gaussians[j*N_COLS_GAUSS + 6]; // ARNAU // Define the horizontal velocity before the longitude as it will be used to define the latter

            if (isdegvel<0 || isdegvel>1) CRASH("Perturbation table error: isdegvel should be 0 or 1");

            double tt=(t-t0_Gauss); // time since perturbation beginning

            double Omega_Gauss;

            if ((isdegvel<0) || (isdegvel>1)) CRASH("Perturbation table error: isdegvel should be 0 or 1");

            if(isdegvel == 0){
                Omega_Gauss = vel_Gauss/rZ(lat_Gauss);
            }
            if(isdegvel == 1){
                Omega_Gauss = deg2rad(vel_Gauss);
            }

            double lon_Gauss0 = lon_Gauss;
            double lat_Gauss0 = lat_Gauss;

            if(SW->polar==0) {
                lon_Gauss = lon_Gauss0 + Omega_Gauss*tt; // Longitude at current time t (lon) MANEL
            }

            // pprintf("lon=%f deg lon=%f rad\n",rad2deg(lon_Gauss), lon_Gauss);

            if((SW->polar==1)||(SW->polar==-1)) {
                double angle0 = 0.0;
                if((lat_Gauss==mid_LAT)&&(lon_Gauss==mid_LON)) CRASH("add_Gaussians Center of coordinates is at the same location where the perturbation has been injected (probably lon,lat=0,0)");
                angle0 = acos((lat_Gauss-mid_LAT)/sqrt(POW2(lat_Gauss-mid_LAT)+POW2(lon_Gauss-mid_LON)));
                if (lon_Gauss-mid_LON>=0.0){
                  angle0 = angle0;
                }else{
                  angle0 = -angle0;
                }

                lon_Gauss = sqrt(POW2(lon_Gauss0)+POW2(lat_Gauss0))*cos(Omega_Gauss*tt-angle0);
                lat_Gauss = -sqrt(POW2(lon_Gauss0)+POW2(lat_Gauss0))*sin(Omega_Gauss*tt-angle0);
            }

            if (SW->prn || info) {
              pprintf("add_Gaussians j=%d lon_Gauss=%f lat_Gauss = %f tt=%e \n",j,rad2deg(lon_Gauss), rad2deg(lat_Gauss),tt);
            }

            Q_Gauss = SW->Gaussians[j*N_COLS_GAUSS + 7]; //  Flow rate

            int isdegsize = SW->Gaussians[j*N_COLS_GAUSS + 8]; // ARNAU
            sigma_Gauss = SW->Gaussians[j*N_COLS_GAUSS + 9]; // ARNAU // Sigma of the normal function
            radius_Gauss = SW->Gaussians[j*N_COLS_GAUSS + 10]; // ARNAU // Distance to the centre of the normal function at which the function will be truncated
            double inject_tracer = SW->Gaussians[j*N_COLS_GAUSS +11];

            if (isdegsize<0 || isdegsize>1) CRASH("Perturbation table error: isdegsize should be 0 or 1");

            if(isdegsize == 1) {
                sigma_Gauss = sigma_Gauss/360.0 * 2*PI*(SW->rE+SW->rP)/2.0; // ARNAU
                radius_Gauss = radius_Gauss/360.0 * 2*PI*(SW->rE+SW->rP)/2.0; // ARNAU
            }

            // if (SW->polar!=0 && vel_Gauss!=0) CRASH("moving perturbations not implemented for polar");

            /* We have only read the Gaussian params up to this point. Now the eta and tracer fields have to be modified. */

            if (info)
                pprintf("Adding perturbation %d/%d at lon=%f (deg) lat=%f (deg) \n",
                    j+1, SW->num_Gaussians,rad2deg(lon_Gauss),rad2deg(lat_Gauss) );

            double Qexp;

            double tauTH = 6*tau_Gauss; // Threshold at 99.7% of max value of step for first order systems

            // if(t1_Gauss - tauTH - (t0_Gauss + tauTH) < 0) CRASH("Perturbation not reaching max value at any point"). // t_intersect is the point where STEP IN and STEP OUT responses cancel out
            double t_intersect = tau_Gauss * log(exp((t1_Gauss - tauTH)/tau_Gauss) + exp(t0_Gauss/tau_Gauss)); // Expression arising from Q*exp(-(t-(t1-tauTH))/T) - Q*(1-exp(-(t-t0)/T)) = 0 and solving for t

            // DEBUG (plot csv)
            // FILE *f = NULL;
            // if(quisoc()==0) {
            //     f = fopen("Q_vs_t.csv", "a");
            //   }

            /*

            See https://www.desmos.com/calculator/36x2czm71c for a temporal plot of the flow rate Q that is being injected.
            To obtain a smooth injection of "perturbation"/flow rate, we have modeled the system as first order, so the
            perturbation reaches the expected/desired value asymptotically.

            */

            // STEP IN: First order system step response
            // tt = t-t0_Gauss
            if((t < tauTH + t0_Gauss) && (t < t_intersect)) {
                Qexp = Q_Gauss*(1.0 - exp(-(t-t0_Gauss)/tau_Gauss));
                if (info) {
                    pprintf("add_Gaussians STEP IN injecting %e of Q=%e\n", Qexp, Q_Gauss);
                    pprintf("add_Gaussians STEP IN: Qexp1=%e (t = %f, [%f %f])\n",Qexp, t, 0.0, tauTH + t0_Gauss);
                }
            }

            // STEP FLAT: Once the value reached is similar to the desired value, we set the flow rate to this value
            if((t0_Gauss+tauTH <= t) && (t <= t1_Gauss-tauTH)) {
                Qexp = Q_Gauss;
                if (info) {
                    pprintf("add_Gaussians STEP FLAT injecting 0.982*Q<=%e of Q=%e or more.\n", Qexp, Q_Gauss);
                    pprintf("add_Gaussians STEP FLAT: Qexp2=%e (t = %f, [%f %f])\n",Qexp, t, t0_Gauss+tauTH, t1_Gauss-tauTH);
                }
            }

            // STEP OUT: First order system response on disappearing step
            if((t > t1_Gauss - tauTH) && (t > t_intersect)) {
                Qexp = Q_Gauss*exp(-(t-(t1_Gauss - tauTH))/tau_Gauss);
                if (info) {
                    pprintf("add_Gaussians STEP OUT injecting %e of Q=%e\n", Qexp, Q_Gauss);
                    pprintf("add_Gaussians STEP OUT: Qexp3=%e (t = %f, [%f inf])\n",Qexp, t, t1_Gauss - tauTH);
                }
            }

            // DEBUG (plot csv)
            // if(quisoc()==0) {
            //     fprintf(f, "%e, %e, %e\n",t, tt, Qexp);
            //     fclose(f);
            // }

            // MPI_Barrier(MPI_COMM_WORLD);

            add_Gaussian_core(lon_Gauss, lat_Gauss, Dt*Qexp, sigma_Gauss, radius_Gauss, inject_tracer, SW); // This creates the distribution. Writes to ETA and TRACER
        }
    }

}


void add_Gaussian_core(double lon, double lat, double volume, double sigma, double radius, double inject_tracer, sw *SW) {

/*

Adds a Gaussian perturbation on the surface of the domain at position (lon,lat) and truncated at r=radius.
The tracer field is updated from 0 to 1 where a perturbation has been injected at get_Gaussian_field(...) function call.
This function evaluates if the injected volume matches the expected volume. The injected volume is rescaled to match the expected volume (volume/Vg).

- double lon (input) The starting longitude in the current instant (maybe outside the domain).
- double lat (input) The starting latitude in the current instant (maybe outside the domain).
- double volume (input) The desired volume to be injected.
- double sigma (input) The standard deviation of the Gaussian perturbation.
- double radius (input) The ending longitude.
- sw *SW (input/output) The Shallow World.

*/

    int i,j;
    int info=0;

    if(IS_PERIODIC_X==0) CRASH("add_Gaussian_core for non periodic domain in X direction is not currently supported.");

    if (SW->prn) pprintf("add_Gaussian_core lon=%e DEG lat=%e DEG volume=%e sigma=%e radius=%e \n",rad2deg(lon),rad2deg(lat),volume,sigma,radius);
    if(SW->polar==0){
        while(lon<SW->lon0){
            lon+=(SW->lon1-SW->lon0);
        }
        while(lon>SW->lon1){
            lon-=(SW->lon1-SW->lon0);
        }
    }

    if((SW->polar==1)||(SW->polar==-1)){
        if((lon<SW->lon0)||(lon>SW->lon1)) CRASH("add_Gaussian_core perturbation out of domain (x direction)");
        if((lat<SW->lat0)||(lat>SW->lat1)) CRASH("add_Gaussian_core perturbation out of domain (y direction)");
    }


        // while(lat<SW->lat0){
        //     lat+=(SW->lat1-SW->lat0);
        // }
        // while(lon>SW->lat1){
        //     lat-=(SW->lat1-SW->lat0);
        // }

    if (info) pprintf("add_Gaussian_core t lon=%.3f DEG lat=%.3f DEG volume=%.3e sigma=%.3e radius=%.3e \n",rad2deg(lon),rad2deg(lat),volume,sigma,radius);

    /*

    If Q=V*Dt and V is the integral of the normal distribution function, its amplitude A can be isolated.
    This expression is derivated in MATLAB as follows:

    >> syms r R sigma A V Dt Q % r: polar variable, R: radius (upper limit of function)
    >> dv = 2*pi*r*A*exp(-r^2/(2*sigma^2)); % Differential volume
    >> V = int(dv,r,0,R) % Volume under normal function from r=0 to r=R
    V =
    -2*A*sigma^2*pi*(exp(-R^2/(2*sigma^2)) - 1)
    >> solve(Dt*Q==V,A) % Find amplitude A
    ans =
    -(Dt*Q)/(2*sigma^2*pi*(exp(-R^2/(2*sigma^2)) - 1)) % Beware!, minus sign at the beginning

    */

    double amp = volume / (2.0*PI*POW2(sigma)*(1.0-exp(-POW2(radius)/(2.0*POW2(sigma))))); // Amplitude = Dt*Q/Volume

    if (info) {
        pprintf("volume=%e\n",volume);
        pprintf("Q=volume/Dt=%e\n",volume/SW->Dt);
        pprintf("PI=%e\n",PI);
        pprintf("sigma=%e\n",sigma);
        pprintf("radius=%e\n",radius);
        pprintf("amp=%e\n",amp);
        pprintf("DT=%e\n",SW->Dt);
        // CRASH("OSTTT!!!");
    }
    if (SW->prn||info) pprintf("add_Gaussian core lon=%.3f lat=%.3f volume=%.3e amplitude=%.3f\n",rad2deg(lon),rad2deg(lat),volume,amp);

    double *perturb_field; // Contains the values to be added to eta
    perturb_field = dmem(SW->M);
    setzero_scaf(perturb_field, SW->M);

    // Calculate the perturbation and store it int perturb_field
    // Because we assume periodic boundary conditions, we call the function three times
    get_Gaussian_field(lon, lat, amp, sigma, radius, perturb_field, inject_tracer, SW); // Perturbation inside the domain (lon,lat are always inside the domain)

    double lonRight=lon+(SW->lon1-SW->lon0);
    double lonLeft=lon-(SW->lon1-SW->lon0);

    if((lon+radius/rZ(lat))>=SW->lon1) get_Gaussian_field(lonLeft, lat, amp, sigma, radius, perturb_field, inject_tracer, SW); // Right periodicity (same perturbation at the right of the domain)
    if((lon-radius/rZ(lat))<=SW->lon0) get_Gaussian_field(lonRight, lat, amp, sigma, radius, perturb_field, inject_tracer, SW); // Left periodicity (same perturbation at the left of the domain)

    double a; // Value of normal function at r (distance to origin)
    double s; // Surface of control volume
    double v; // Injected volume at control volume
    double V = 0.0; // Total injected volume (sum of v accross all the domain) to be compared with input volume
    double Vg; // Global total injected volume

    if(volume == 0.0) {
        // CRASH("Volume to be injected is 0. Set volume to another value."); // ARNAU
        pprintf("Volume to be injected is currently 0.\n"); // ARNAU
    }

    forall(i,j,SW->M){
        a = ac(perturb_field,i,j,SW->M); // Perturbation height at control volume i,j
        // if(a!=0) pprintf("a=%e, s=%e, v=%e\n",a,s,v);
        s = DXC(i,j)*DYC(i,j); // ARNAU // SHOULD BE REVISED (THIS IS USED FOR CHECKING PURPOSES) Surface of control volume i,j
        v = a*s; // Injected volume at control volume i,j
        V += v; // Add to sum (total injected volume)
    }

    checkr(MPI_Allreduce(&V, &Vg, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD),"V");
    // pprintf("V=%e, Vg = %e\n",V,Vg);

    // Up to this point we have computed the total injected volume V, now the perturbation has to be rescaled so Vcheck==V
    double Vcheck = 0.0, Vcheckg;

    if(fabs(volume)>1e-7) { // ignore if too small

      if (Vg==0.0) CRASH("Global total injected volume Vg=%e equals 0",Vg);

      double rf = (volume/Vg); // Reescale factor
      if (SW->prn||info) pprintf("Perturbation reescale factor is %e\n",rf);

      forall(i,j,SW->M){
          a = rf*ac(perturb_field,i,j,SW->M); // Rescale height due to perturbation at control volume i,j
          ETA(i,j)+=a;
          s = DXC(i,j)*DYC(i,j); // Surface of control volume i,j
          v = a*s; // Injected volume at control volume i,j
          Vcheck += v; // Add to sum (total injected volume check)
      }

      checkr(MPI_Allreduce(&Vcheck, &Vcheckg, 1, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD),"Vcheck");
      // pprintf("Vcheck=%e, Vcheckg = %e\n",Vcheck,Vcheckg);

      double epsilon = 1e-1; // Unitary error, a.k.a. maximum assumable relative error between injected volume and desired volume
      if(fabs(Vcheckg-volume)/volume>epsilon) CRASH("Total injected volume Vcheck=%e does not match input volume=%e (difference is=%e)", Vcheckg, volume, fabs(Vcheckg-volume));

    }

    free(perturb_field);
}


void get_Gaussian_field(double lon, double lat, double amp, double sigma, double radius, double *field, double inject_tracer, sw *SW) {

/*
The function calculates the height of the fluid that corresponds to each point of the domain taking into account its distance with respect to the center of the perturbation. It must be taken into account that the Gaussian function is truncated at a radius r = radius.
It also updates the value of the tracer to 1 if a perturbation is addded.

- double lon (input) The starting longitude.
- double lat (input) The starting latitude.
- double amp (input) The amplitude of the Gaussian perturbation.
- double sigma (input) The standard deviation of the Gaussian perturbation.
- double radius (input) The ending longitude.
- sw *SW (input/output) The Shallow World.

When computing the distance between the centre of the perturbation to a point of a region of a spheroid, a formula meant for spheres is used (formula number 80 of Williamson 1992--A Standard Test Set for Numerical Approximations to the Shallow Water Equations in Spherical Geometry). This is what get_great_circle(...) does.

*/

    int i,j;
    double r;
    int fill_mode=1; // Mode to determine the area of influence of the perturbation
    int info=0;

    if (SW->prn||info) pprintf("get_Gaussian_field entering lon=%f DEG lat=%f DEG\n",rad2deg(lon),rad2deg(lat));

    // if(radius<DXC(i,j)||radius<DYC(i,j)) CRASH("Diameter of perturbation is too small or grid is too coarse");
    // r = fabs(get_vincenty(lon, lat, LONC(i,j), LATC(i,j), SW->rE, SW->rP)); // Alternative to get_great_circle(...) which works for ellipses (see description of function)
    switch(fill_mode) {
        case 1: // Or if center of perturbation inside control volume
            forall(i,j,SW->M){
                r = get_great_circle(lon, lat2pclat(lat), LONC(i,j), lat2pclat(LATC(i,j)), (SW->rE+SW->rP)/2.0); // Great circle length
                if(r<radius) {
                    if(info) pprintf("r=%f radius=%f cond=%d\n",r,radius,r<radius);
                    ac(field,i,j,SW->M)+= amp*exp(-POW2(r)/(2.0*POW2(sigma))); // If r<radius, update position i,j with height of perturbation at control volume i,j
                    TRACER(i,j) = (inject_tracer)*exp(-POW2(r)/(2.0*POW2(sigma)));
                    // TRACER(i,j) += (inject_tracer)*exp(-POW2(r)/(2.0*POW2(sigma)));
                    // if (TRACER(i,j) > inject_tracer) TRACER(i,j) = inject_tracer;
                    PERTURBS(i,j) = 1.0;
                }

            }
            break;

        case 2:
            // Define new fill mode
            break;

    }
}

void add_vortices(double Dt, double t, sw *SW){

    int info=1;

    double t0_vortex;
    int j;

    double lonc,latc,as,bs,Umax,n;

    for(j=0; j<=SW->num_vortices-1;j++){

        //  0    1    2  3  4    5 6
        // t0 lonc latc as bs Umax n)

        t0_vortex = SW->vortices[j*N_COLS_VORTICES + 0];

        if((t0_vortex<t)&&(t0_vortex>=(t-Dt))) {

            lonc = SW->vortices[j*N_COLS_VORTICES + 1]; // eg 60 ;    central longitude of vortex [deg]
            latc = SW->vortices[j*N_COLS_VORTICES + 2]; // eg -21;    central (planetographic) latitude of vortex [deg]
            as = SW->vortices[j*N_COLS_VORTICES + 3]; // eg 7.4;    longitude semiaxis [deg]
            bs = SW->vortices[j*N_COLS_VORTICES + 4]; // eg 4.7;      //latitude semiaxis [deg]
            Umax = SW->vortices[j*N_COLS_VORTICES + 5]; // eg -100.0; //maximum tangential velocity [m/s]
            n = SW->vortices[j*N_COLS_VORTICES + 6]; // eg 2.0;       //paramter to control how flat is the Gaussian distribution

            bs = fabs((bs-latc))/pow(1.0-1.0/(2.0*(n)),1.0/(2.0*n));
            as = fabs((as-lonc))/pow(1.0-1.0/(2.0*(n)),1.0/(2.0*n));
            pprintf("as = %e \n",as);
            pprintf("bs = %e \n",bs);
            double amp_eta = 0.5*Umax*rM(deg2rad(latc))*deg2rad(bs)*f(deg2rad(latc))*pow((1.0-1.0/(2.0*n)),((1.0-2.0*n)/(2.0*n)))*exp(1.0-1.0/(2.0*n))/(g_eff(deg2rad(latc))*n);
            pprintf("amp_eta = %e \n",amp_eta);

            if(info) pprintf("add_vortices t0_vortex=%f, lonc=%f DEG, latc=%f DEG, as=%f DEG, BS=%f DEG, Umax=%f, n=%f\n",t0_vortex,lonc,latc,as,bs,Umax,n);

            add_vortex_core(SW->eta, SW->u, SW->v, n, amp_eta, lonc, latc, as, bs, SW);
            halo_update(SW->eta,SW->M);
            halo_update(SW->u,SW->M);
            halo_update(SW->v,SW->M);

        }

    }

} // ARNAUMARC

//  Geostrophic equilibrium

void add_vortex_core(double *eta, double *u, double *v, double n, double amp, double loncent, double latcent, double as, double bs, sw *SW) { // everything in degrees
    int i, j;
    int info=1;
    double exponent;
    double eta_ij;
    double psi_u, psi_v;
    // double R = (SW->rE+SW->rP)/2.0;
    forall(i,j,SW->M){
        exponent = pow(POW2((rad2deg(LONC(i,j))-loncent)/as) + POW2((rad2deg(LATC(i,j))-latcent)/bs),n);
        eta_ij = amp*exp(-exponent);
        ac(eta,i,j,SW->M) += eta_ij;

        psi_u = POW2((rad2deg(LONE(i,j))-loncent)/as) + POW2((rad2deg(LATE(i,j))-latcent)/bs);
        ac(u,i,j,SW->M) += 2.0*(n)*GE(i,j)/(f(deg2rad(latcent))*POW2(deg2rad(bs))*rM(deg2rad(latcent)))*(LATE(i,j)-deg2rad(latcent))*pow(psi_u,n-1.0)*eta_ij;

        psi_v = POW2((rad2deg(LONN(i,j))-loncent)/as) + POW2((rad2deg(LATN(i,j))-latcent)/bs);
        ac(v,i,j,SW->M) += -2.0*(n)*GN(i,j)/(f(deg2rad(latcent))*POW2(deg2rad(as))*rZ(deg2rad(latcent)))*(LONN(i,j)-deg2rad(loncent))*pow(psi_v,n-1.0)*eta_ij;
    }

    if(info) pprintf("add_vortex_core lonc=%f DEG, latc=%f DEG, as=%f DEG, bs=%f DEG\n",loncent,latcent,as,bs);

}


// void vortex_eta(double *eta, double n, double amp, double loncent, double latcent, double as, double bs, sw *SW){
//   int i, j;
//   double exponent;
//   int info=0;
//   forall(i,j,SW->M){
//     exponent = pow(POW2((rad2deg(LONC(i,j))-loncent)/as) + POW2((rad2deg(LATC(i,j))-latcent)/bs),n);
//     ac(eta,i,j,SW->M) += amp*exp(-exponent);
//   }
//   if(info) print_scaf(eta, "ETA init", 0, SW->M,"%3e ");
// }
//
// void vortex_u(double *u, double *eta, double n,double loncent, double latcent, double as, double bs, sw *SW){
//   int i, j;
//   int info=0;
//   double psi;
//   double R = (SW->rE+SW->rP)/2.0;
//   forall(i,j,SW->M){
//     psi = POW2((rad2deg(LONE(i,j))-loncent)/as) + POW2((rad2deg(LATE(i,j))-latcent)/bs);
//     ac(u,i,j,SW->M) += 2.0*n*GE(i,j)/(FE(i,j)*POW2(deg2rad(bs))*R)*(LATE(i,j)-deg2rad(latcent))*pow(psi,n-1.0)*ac(eta,i,j,SW->M);
//   }
//   if (info) print_scaf(u, "U init", 0, SW->M,"%3e ");
// }
//
// void vortex_v(double *v,  double *eta,double n,double loncent, double latcent, double as, double bs, sw *SW){
//   int i, j;
//   int info=0;
//   double psi;
//   double R = (SW->rE+SW->rP)/2;
//   forall(i,j,SW->M){
//     psi = POW2((rad2deg(LONN(i,j))-loncent)/as) + POW2((rad2deg(LATN(i,j))-latcent)/bs);
//     ac(v,i,j,SW->M) += -2.0*(n)*GN(i,j)/(FN(i,j)*R*cos(LATN(i,j))*POW2(deg2rad(as)))*(LONN(i,j)-deg2rad(loncent))*pow(psi,n-1.0)*ETA(i,j);
//   }
//   if (info) print_scaf(v, "V init", 0, SW->M,"%3e ");
// }
