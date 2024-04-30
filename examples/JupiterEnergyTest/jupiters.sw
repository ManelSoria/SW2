# Parameters for this run (case.161001)
define(P1,-2.250000e+01)
define(P2,4.520000e+11)
define(P3,3.500000e+06)

# Template (template_61.sw):
# Nature 2019, les dues tempestes polars juntes

define(PI,3.14159265)
define(DAY,86400)
define(DT,50) # time step (s)

sch   1
polar 0

lon0  0.0          # Start longitude (1 deg)
lon1  90.0         # End longitude
lat0 -35.0         # Start latitude (1 deg)
lat1 -10.0         # End latitude

nx 50
ny 50

t0  0.0       # Start time (s)
t1  <DAY*20>  # End time # !!!!!!!!!                  

Dt  DT   # Time step
CFL 0.25 # Max CFL allowed

Tau       1e7
Hypernu   0. 0. 0. 
IteInfo   100
SaveEvery <12*3600/DT>
LoadFrom  0 # Last file

reftype    0 # 0: no comparison

include(JupiterEGM.planet)

depth 1000

geoeq 0

# Table for vortices in geostrophic equilibrium
0    # Rows (0 if empty)
7    # Columns 
# t0     lon     lat       a(deg)      b(deg)     u_max(m/s)      n

# Table for Gaussian perturbations
1 # 
12 # it has to be the number of columns N_COLS_GAUSS !!
#t0   t1(s)    tau   lon   lat   isdegvel vel(m/s) Q(m^3/s) isdegsize sigma(m)  radius(m) inject tracer
0     <7*DAY>  0     45    P1    0        -3.6     P2       0         3.5e6     P3        0 

# Time constant of the dissipation of the tracer
TracerDis 0

# Wind profiles
include(nowind.wind)
