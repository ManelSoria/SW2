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

double linterp_core(double x_i, double *x, double *y, int *i0, int *i1) {

/*

Returns the linearly interpolated value of the function y=y(x) at point x_i.

- double x_i (input) The value at which the interpolation is done (function evaluation point).
- double *x (input) The points whose image of y=y(x) is known.
- double *y (input) The images of x.
- int *i0 (input/output) The search start index of x.
- int *i1 (input/output) The search end index of x.

This function only does an out-of-bonds check. The interpolated value is found using only one while loop.
The search start index is updated with the index at which the interpolated value has been found.

*/

    int info=0;
    int i=(*i0);
    if(x[*i0]>x_i) CRASH("Value out of bounds (smaller than the smallest value of the vector x[first]=%e and x_i=%e)",x[*i0],x_i);
    if(x[*i1]<x_i) CRASH("Value out of bounds (greater than the greatest value of the vector x[last]=%e and x_i=%e)",x[*i1],x_i);
    while((x[i+1]<=x_i)&&(i<(*i1))&&(i>=(*i0))) {
        i++;
        if(info>1) pprintf("NEXT STEP: x[%d]=%e x[%d]=%e (x[i]<=x_i)=%d (i<=(*i1))=%d (i>=(*i0))=%d\n",i,x[i],i-1,x[i-1],(x[i]<=x_i),(i<=(*i1)),(i>=(*i0)));
    }
    if(info) pprintf("Selected position is i=%d (i0=%d, i1=%d)\n",i,(*i0),(*i1));
    if(info>1) pprintf("y[i+1]=%e,y[i]=%e,x_i=%e,x[i]=%e,x[i+1]=%e,x[i]=%e\n",y[i+1],y[i],x_i,x[i],x[i+1],x[i]);
    if(info) pprintf("Returning=%e\n",(y[i+1]-y[i]) * (x_i-x[i]) / (x[i+1]-x[i]) + y[i]);
    (*i0)=i; // Return last valid value
    if(info) pprintf("New index i0=i=%d\n",(*i0));
    return (y[i+1]-y[i]) * (x_i-x[i]) / (x[i+1]-x[i]) + y[i];
}

double linterp(double x_i, double *x, double *y, int n) {

/*

Returns the linearly interpolated value of the function y=y(x) at point x_i.

- double x_i (input) The value at which the interpolation is done (function evaluation point).
- double *x (input) The points whose image of y=y(x) is known.
- double *y (input) The images of x.
- int n (input) The length of vectors x and y.

This is a redefinition of linterp_core for easy use (starting and ending indeces are 0 to length of the array).

*/

    int info=0;
    int i0=0;
    int i1=n-1;
    if(info) pprintf("Interpolation result: %e\n",linterp_core(x_i, x, y, &i0, &i1));
    return linterp_core(x_i, x, y, &i0, &i1);
}

void linterp_vec(double *x_i, double *y_i, int n_i, double *x, double *y, int n) {

/*

Returns y_i, which are the linearly interpolated values of the function y=y(x) at points x_i.

- double *x_i (input) The values at which the interpolation is done (function evaluation points).
- double *y_i (output) The images of $y=y(x)$ at points x_i.
- int n_i (input) The lenght of vectors x_i and y_i.
- double *x (input) The points whose image of y=y(x) is known.
- double *y (input) The images of x.
- int n (input) The length of vectors x and y.

This function assumes that the values are sorted by increasing order. This calls linterp_core as this function returns updated indices (the updated index is that of the last found value of y).

*/

    int info=0;
    int i;
    int i0=0;
    int i1=n-1;
    for(i=0;i<n_i;i++)
        y_i[i]=linterp_core(x_i[i], x, y, &i0, &i1);
    if(info)
        for(i=0;i<n_i;i++) pprintf("%e ",y_i[i]);

}

double get_great_circle(double lon0, double lat0, double lon1, double lat1, double r) {

/*

Computes the great circle distance on a sphere of radius r.

- double lon0 (input) The starting longitude.
- double lat0 (input) The starting latitude.
- double lon1 (input) The ending longitude.
- double lat1 (input) The ending latitude.
- double r (input) The radius of the sphere.

*/

    return (r)*acos(sin(lat0)*sin(lat1)+cos(lat0)*cos(lat1)*cos(lon1-lon0));

}

double get_great_angle(double lon0, double lat0, double lon1, double lat1) {

/*

Computes the great circle angle on a sphere.

- double lon0 (input) The starting longitude.
- double lat0 (input) The starting latitude.
- double lon1 (input) The ending longitude.
- double lat1 (input) The ending latitude.

*/

    return acos(sin(lat0)*sin(lat1)+cos(lat0)*cos(lat1)*cos(lon1-lon0));

}

double get_vincenty(double lon0, double lat0, double lon1, double lat1, double a, double b) {

/*

Computes the geodetic between points (lambda_0,varphi_0) and (lambda_1,varphi_1) according to the Vincenty formulae.

- double lon0 (input) The starting longitude.
- double lat0 (input) The starting [planetographic(?)] latitude.
- double lon1 (input) The ending longitude.
- double lat1 (input) The ending [planetographic(?)] latitude.
- double a (input) The semi-major axis of the ellipse.
- double b (input) The semi-minor axis of the ellipse.

Instead of using the haversine formula on an ellipsoid to determine the geodetic, one shall use the Vincenty algorithm to determine with much more precision the geodetic between two points on a ellipsoid given the angles lambda and phi. The geodetic is used in Shallow Worlds to compute e.g. the area of influence of a vortex.

The Vincenty algorithm may fail to converge for points close to the poles and, while it converges very fast for points close to the equator, the convergence rate decreases for increasing latitudes.
Only use get_vincenty when the distance between relatively close points (such as in the case of boundary to center points of a vortex) is needed.
Also, do not use this formula to compute the distance between antipodal points as it will produce nonsense results.
The GeographicLib implements an algorithm that computes the geodetic even when the Vincenty algorithm fails. It may not be as fast as this function.

*/

    int info=0;
    double s;
    if((lon0==lon1)&&(lat0==lat1)) {
        if(info)pprintf("Same points referenced\n");
        s=0.0;
    }
    double f=(a-b)/a;
    double r_lat0=atan((1.0-f)*tan(lat0));
    double r_lat1=atan((1.0-f)*tan(lat1));
    double cos_r_lat0=cos(r_lat0);
    double cos_r_lat1=cos(r_lat1);
    double sin_r_lat0=sin(r_lat0);
    double sin_r_lat1=sin(r_lat1);
    double Dlon=lon1-lon0;
    double old_lambda;
    double lambda=Dlon;
    double cos_lambda, sin_lambda;
    double sin_sigma, cos_sigma;
    double sigma;
    double sin_alpha;
    double cos2_alpha;
    double cos_2sigmam;
    double C;
    int i=0;
    int niter=1e3;
    do{
        cos_lambda=cos(lambda);
        sin_lambda=sin(lambda);
        old_lambda=lambda;
        if(info>1)pprintf("old_lambda=%e->",old_lambda);
        i++;
        sin_sigma=sqrt(POW2(cos_r_lat1*sin_lambda)+POW2(cos_r_lat0*sin_r_lat1-sin_r_lat0*cos_r_lat1*cos_lambda));
        cos_sigma=sin_r_lat0*sin_r_lat1+cos_r_lat0*cos_r_lat1*cos_lambda;
        sigma=atan2(sin_sigma, cos_sigma);
        sin_alpha=(cos_r_lat0*cos_r_lat1*sin_lambda)/sin_sigma;
        cos2_alpha=1-POW2(sin_alpha);
        if(cos2_alpha==0){
            lambda = Dlon + f*sin_alpha*sigma;
        } else {
            cos_2sigmam=cos_sigma-(2.0*sin_r_lat0*sin_r_lat1)/cos2_alpha;
            C=(f/16.0)*cos2_alpha*(4.0+f*(4.0-3.0*cos2_alpha));
            lambda=Dlon+(1.0-C)*f*sin_alpha*(sigma+C*sin_sigma*(cos_2sigmam+C*cos_sigma*(-1.0+2.0*POW2(cos_2sigmam))));
        }
        if(info>1)pprintf("lambda=%e\n",lambda);
    }while((fabs(old_lambda-lambda)>=1e-10)&&(i<=niter));
    if(info) pprintf("lambda=%e after %d iterations\n", lambda, i);
    if(info&&(i>=niter)) pprintf("Maximum number of iterations achieved %d(>=%d)\n", i, niter);

    double u2=cos2_alpha*(POW2(a)-POW2(b))/POW2(b);
    double A=1.0+(u2/16384.0)*(4096.0+u2*(-768.0+u2*(320.0-175.0*u2)));
    double B=(u2/1024.0)*(256.0+u2*(-128.0+u2*(74.0-47.0*u2)));
    double Dsigma=B*sin(sigma)*(cos_2sigmam+(1.0/4.0)*B*(cos_sigma*(-1.0+2.0*POW2(cos_2sigmam))-(1.0/6.0)*B*cos_2sigmam*(-3.0+4.0*POW2(sin_sigma))*(-3.0+4.0*POW2(cos_2sigmam))));
    s=b*A*(sigma-Dsigma);
    if(info) pprintf("distance=%e\n",s);
    return s;
}

double get_arc_ellipse(double a1, double rP, double epsilon2) {

/*

Returns the arc of ellipse from 0 to a1, in radians.

- double a1 (input) The end angle.
- double rP (input) The semi-minor axis (or b) of the ellipse (polar radius).
- double epsilon2 (input) For an ellipse of radius a and b, varepsilon=a/b. varepsilon is the quotient between the equatorial radius rE and the polar radius rP.

*/

    int N=1e6; // Number of samples
    double da=a1/N,a=0.; // Centric, NOT graphic angles
    double ds=0,s=0;
    int i;
    for(i=0;i<N;i++) {
        a+=da/2.;
        ds=da*rP*sqrt(1-(1-epsilon2)*POW2(sin(a)));
        s+=ds;
    }
    return s;
}

inline double get_arc_circumf(double a1, double r) {

/*

Returns the arc of circumference from 0 to a1, in radians.

- double a1 (input) The end angle.
- double r (input) The radius of the circumference.

*/

    return a1*r;
}

inline double get_linspace_val(int i,double x0,double x1,int nx) {

/*

Returns the value at index i of a vector of n_x values uniformly distributed from x_0 to x_1. For i=0 the value returned is x0.

- int i (input) The element of the vector to be returned.
- double x0 (input) The first value of the vector, found at index i=0.
- double x1 (input) The last value of the vector, found at index i=n_x.
- int nx (input) The number of values in the vector.

At the core of init_coords_ortho, this function is used to distribute the coordinates uniformly. This function is used to create the coorinates of the points of the halo to avoid calls to halo_update, so indices i<0 and i>n_x are also valid.
This function is also used in the definition lone, lonc, lonv, latn, latc, and latv.

*/

    return ( x0+(double)i*(x1-x0)/((double)nx) ); // Set coordinate at division
}

void generate_scaf(double *scaf, double (*funcio)(double,double,double), double x0, double y0, int stgx, int stgy, double t, sw *SW) {

/*

Generates a scalar field according to a mathematical function defined outside of this function.

- int scaf (output) The element of the vector to be returned.
- double (*funcio)(double,double,double) (input) The name of the function f(x,y,t) to call which will be used to compute the values of scaf.
- int stgx(input) If 0, centered coordinates will be used. If 1, staggered coordinates will be used.
- int stgx (input) If 0, centered coordinates will be used. If 1, staggered coordinates will be used.
- double t (input) The current simulation time.
- sw *SW (input) The Shallow World.

*/

    int i,j;

    // identify coordinate field
    if(stgx&&stgy)
        forall(i,j,SW->M) {
            ac(scaf,i,j,SW->M) = funcio( XE(i,j) - x0, YN(i,j) - y0, t);
        }
    else if(stgx)
        forall(i,j,SW->M) {
            ac(scaf,i,j,SW->M) = funcio( XE(i,j) - x0, YC(i,j) - y0, t);
        }
    else if(stgy)
        forall(i,j,SW->M) {
            ac(scaf,i,j,SW->M) = funcio( XC(i,j) - x0, YN(i,j) - y0, t);
        }
    else
        forall(i,j,SW->M) {
            ac(scaf,i,j,SW->M) = funcio( XC(i,j) - x0, YC(i,j) - y0, t);
        }
}

double generate_val(int i, int j, double (*funcio)(double,double,double), double x0, double y0, int stgx, int stgy, double t, sw *SW) {

/*

Generates a scalar field according to a mathematical function defined outside of this function.

- int i (input) The position of the element in the x direction.
- int j (input) The position of the element in the y direction.
- double (*funcio)(double,double,double) (input) The name of the function f(x,y,t) to call which will be used to compute the values of scaf.
- int stgx (input) If 0, centered coordinates will be used. If 1, staggered coordinates will be used.
- int stgx (input) If 0, centered coordinates will be used. If 1, staggered coordinates will be used.
- double t (input) The current simulation time.
- sw *SW (input) The Shallow World.

*/

    // identify coordinate field
    if(stgx&&stgy)
        return funcio( XE(i,j) - x0, YN(i,j) - y0, t);
    else if(stgx)
        return funcio( XE(i,j) - x0, YC(i,j) - y0, t);
    else if(stgy)
        return funcio( XC(i,j) - x0, YN(i,j) - y0, t);
    else
        return funcio( XC(i,j) - x0, YC(i,j) - y0, t);
}
