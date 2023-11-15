#!/bin/env python
import numpy as np, sympy as sy, matplotlib.pyplot as plt

lon0  = np.deg2rad(90)
lat0  = np.deg2rad(-22.5)
a0    = np.deg2rad(22.5)
b0    = np.deg2rad(6.5)
Umax  = -100
nn    = 2

# Define a domain
lon = np.linspace(0,180,100)
lat = np.linspace(-40,-5,100)

llon,llat = np.meshgrid(np.deg2rad(lon),np.deg2rad(lat))

# Jupiter
Omega    = 1.76e-4    # rad/s (Rotation frequency omega = 2*PI/T, where T is the Rotation Period)
rE       = 71492000.0 # m (Equatorial radius)
rP       = 66852000.0 # m (Polar radius)
gE       = 22.58      # m/s2 (Reference gravity)
epsilon2 = (rE/rP)**2


## Derivation of the equations for a geostrophic equilibrium vortex
phi, theta, phi0, theta0, a, b, A, n = sy.symbols('phi theta phi_0 theta_0 a b A n')
rZ, rM, g, f, eta_v, psi_v, umax = sy.symbols('rZ rM g f eta psi u_max')

# Perturbation
psi = ((theta - theta0)/a)**2 + ((phi - phi0)/b)**2
eta = A*sy.exp(-(psi)**n)
print('eta = ',eta)

# Gradients
deta_lon = sy.diff(eta,theta)
deta_lat = sy.diff(eta,phi)
deta_x   = (1/rZ)*deta_lon
deta_y   = (1/rM)*deta_lat

# Velocities
u = (-g/f*deta_y)
v = ( g/f*deta_x)
print('u = ',sy.simplify(u.subs(psi,psi_v).subs(eta,eta_v)))
print('v = ',sy.simplify(v.subs(eta,eta_v).subs(psi,psi_v)))

# Locations of maximum velocity
phi_max   = sy.solve(sy.diff(u.subs(theta,theta0),phi),phi,domain=sy.S.Reals)[1]
print('phi_max = ',sy.simplify(phi_max))
theta_max = sy.solve(sy.diff(v.subs(phi,phi0),theta),theta,domain=sy.S.Reals)[1]
print('theta_max = ',sy.simplify(theta_max))

# Find amplitude for maximum velocity
A_max = sy.solve(umax - u.subs([(theta,theta0),(phi,phi_max)]),A)[0]
print('A = ',sy.simplify(A_max))

# Functions
r     = lambda lat     : rE*rP/np.sqrt((rE*np.sin(lat))**2 + (rP*np.cos(lat))**2) # Planetocentric radius of the spheroid
rZ_f  = lambda lat     : rE/np.sqrt(1.0+(np.tan(lat))**2/epsilon2)
rM_f  = lambda lat     : (rE/epsilon2)*(rZ_f(lat)/(rE*np.cos(lat)))**3
f_f   = lambda lat     : 2.0*Omega*np.sin(lat)                 # Coriolis parameter
g0    = lambda lat     : (gE + Omega**2*rE)*(rE/r(lat))**2     # Reference gravity (potential)
g_eff = lambda lat     : g0(lat) - Omega**2*r(lat)*np.cos(lat) # Effective gravity of the rotating spheroid

#b = np.abs((b - lat0))/np.power(1.0-1.0/(2.0*(n)),1.0/(2.0*n))
#a = np.abs((a - lon0))/np.power(1.0-1.0/(2.0*(n)),1.0/(2.0*n))
#b = np.abs(b)/np.power(1.0-1.0/(2.0*(n)),1.0/(2.0*n))
#a = np.abs(a)/np.power(1.0-1.0/(2.0*(n)),1.0/(2.0*n))
print('a = ',np.rad2deg(a0),'b = ',np.rad2deg(b0))

# Gaussian amplitude
A_fun = sy.lambdify([b,f,rM,umax,n,g],A_max,'numpy')
A_sy  = A_fun(b0,f_f(lat0),rM_f(lat0),Umax,nn,g_eff(lat0))
A_SW  = 0.5*Umax*rM_f(lat0)*b0*f_f(lat0)*np.power((1.0-1.0/(2.0*nn)),((1.0-2.0*nn)/(2.0*nn)))*np.exp(1.0-1.0/(2.0*nn))/(g_eff(lat0)*nn)
print('A = ', A_sy, A_SW)

# eta field
eta_f = sy.lambdify([A,theta,theta0,phi,phi0,a,b,n],eta)
eta   = eta_f(A_SW,llon,lon0,llat,lat0,a0,b0,nn)

# velocity field
u_fun = sy.lambdify([A,theta,theta0,phi,phi0,a,b,n,g,f,rM,rZ],u)
u     = u_fun(A_SW,llon,lon0,llat,lat0,a0,b0,nn,g_eff(lat0),f_f(lat0),rM_f(lat0),rZ_f(lat0))
v_fun = sy.lambdify([A,theta,theta0,phi,phi0,a,b,n,g,f,rM,rZ],v)
v     = v_fun(A_SW,llon,lon0,llat,lat0,a0,b0,nn,g_eff(lat0),f_f(lat0),rM_f(lat0),rZ_f(lat0))

# Plots
plt.figure(figsize=(10,10),dpi=100,facecolor='w',edgecolor='k')
plt.subplot(411)
plt.contourf(lon,lat,eta,100,cmap=plt.cm.jet)
plt.colorbar(label=r'$\eta$')
plt.xlim([lon[0],lon[-1]])
plt.ylim([lat[0],lat[-1]])
plt.subplot(412)
plt.contourf(lon,lat,u,100,cmap=plt.cm.jet)
plt.colorbar(label=r'$u$')
plt.xlim([lon[0],lon[-1]])
plt.ylim([lat[0],lat[-1]])
plt.subplot(413)
plt.contourf(lon,lat,v,100,cmap=plt.cm.jet)
plt.colorbar(label=r'$v$')
plt.xlim([lon[0],lon[-1]])
plt.ylim([lat[0],lat[-1]])
plt.subplot(414)
plt.contourf(lon,lat,np.sqrt(u*u + v*v),100,cmap=plt.cm.jet)
plt.colorbar(label=r'$||u||$')
plt.xlim([lon[0],lon[-1]])
plt.ylim([lat[0],lat[-1]])


plt.show()