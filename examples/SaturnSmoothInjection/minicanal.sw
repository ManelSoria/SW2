#Parameters for this run (case.123001)
define(P1,6.880000e+01)
define(P2,2.000000e+09)
define(P3,1.500000e+09)

#Template (template_23.sw):
# Nature 2019, tempesta de dalt WS2

define(PI,3.14159265)
define(DAY,86400)
define(DT,60) # time step (s)

polar 0 # no polar

lon0 0.0          # Start longitude (1 deg)
lon1 60.0         # End longitude
lat0 65.0         # Start latitude (1 deg)
lat1 75.0         # End latitude

#petita:
nx 600
ny 100

t0 0.0   # Start time (s)
t1  <DT*100>  # End time # !!!!!!!!!

Dt DT # Time step
Tau 1e9
IteInfo 1 # !!!!!!!!!
SaveEvery 1 # !!!!!!!!!
LoadFrom 0 # Last file

reftype    0 # 0: no comparison

include(SaturnEGM2.planet)

depth 500

geoeq 0

1 #
11 # it has to be the number of columns N_COLS_GAUSS !!
#t0         t1(s)       tau   lon   lat   isdegvel vel(m/s) Q(m^3/s) isdegsize sigma(m)  radius(m)
#<0>  <99*DT>  320   180    67.7 0        59.8     P2       0         121e3     121e3
<0>  <20*DT>  200   30     P1   0        14.2     P3       0         114e3     114e3

#wind profiles
#include(nowind.wind)
include(Saturn_Cassini.wind)
