// SW2
// 2018-2020
// Manel Soria, Arnau Prat, Arnau Sabates, Marc Andres, Enrique Garcia-Melendo
// UPC - ESEIAAT - TUAREG

#define VERSION "SW 0.125"

#define G_universal 6.67408e-11 // Universal gravity constant
#define PI (4.0*atan(1.0)) // Pi value dependant on math.h's atan(...) implementation
// #define PI M_PI // Fixed pi value as defined in math.h

#define MAXOUTLEN 1000


typedef struct {

    /*
    Preferably, avoid changing the variables of this struct on every loop iteration (including *Gaussians), except for u, v, eta.
    */
  int sch;

  int polar;

    int N; // Number of timesteps
    int loadFrom; // Load from frame
    int saveEvery; // Period of data saving
    int IteInfo; // Print information every timesteps
    int prn; // print information on this time step (specified in the main loop by core)
    double t0, t1, Dt; // Time info
    double Courant_max;  // Maximum Courant
    int nx, ny, nz; // Number of cells in the centered mesh (there are nx+1, ny+1,... division lines)
    double lon0, lat0; // Start longitude and latitude
    double lon1, lat1; // End longitude and latitude
    map m; // Intrinsic properties of the domain
    map *M; // Pointer to m, to ease notation

    /* Scalar values */
    double Omega; // Angular velocity
    double rE,rP; // Equatorial and polar radii of the spheroid
    double epsilon2; // Radius ratio of the spheroid squared
    double gE; // Reference gravity at equator
    // double mass; // Spheroid average mass
    double depth; // Reference depth of the fluid layer
    double sigma_dis; //Magnitude of the dissipation
    double tracer_dis; // Tracer dissipation

    double beta_sponge;
    double sigma_sponge;
    double tau_sponge;

    /* Scalar fields */

    double *windx; // Winds in x direction
    double *windy; // Winds in y direction
    double *u,*v,*eta,*hB; // Horizontal and vertical velocities, surface perturbation and terrain roughness
    double *tracer, *perturbs; // Passive tracer
    double *xe; // Coordinates of the easterners/staggered x/east points
    double *xc; // Coordinates of the centered points
    double *xv; // Coordinates of the vertex points
    double *yn; // Coordinates of the northerners/staggered y/north points
    double *yc; // Coordinates of the centered points
    double *yv; // Coordinates of the vertex points
    double *Dxe,*Dxc,*Dxv,*Dyn,*Dyc,*Dyv; // Distances between centered, staggered and vertex points, respectively for x and y
    double *ge,*gn,*gc,*gv,*fe,*fn,*fc,*fv; //  Gravity and Coriolis parameter


    double nu2, nu4, nu6;


    /* Reference fields */
    int reftype;
    char refname[MAXOUTLEN];
    double *reftable;
    int reftable_nr,reftable_nc;

    /* Geostrophic equilibrium */
    int geoeq; // =1 : we begin from geostrophic eq

    /* Perturbations */
    int num_Gaussians;
    double *Gaussians;

    int num_vortices;
    double *vortices;

    /* Other parameters */
    // Temporal resolution (i.e. save time) to be defined here if needed
} sw;

#define N_COLS_GAUSS 12 // Number of columns in Gaussian perturbations files
#define N_COLS_VORTICES 7 // Number of columns in vortices files
#define IS_PERIODIC_X 1 // Periodicity in x direction flag

// ACCESS MACROS FOR THE SHALLOW WORLD VARIABLES

// Commonly used macros
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define POW2(a) ((a)*(a))
#define POW3(a) ((a)*(a)*(a))
#define deg2rad(a)  ((a)*(PI/180.0))
#define rad2deg(b)  ((b)*(180.0/PI))
#define MEAN(a,b) (((a)+(b))/2.0)

// Spheroidal coordinates
#define lone(i) (get_linspace_val(i,SW->lon0,SW->lon1,SW->nx)) // Longitude of the staggered x/east point. Longitudes can only be staggered in x
#define lonc(i) ((lone(i)+lone(i-1))/2.) // Longitude of the centered point. The centered point falls in the middle of staggered points
#define lonv(i) (lone(i))
#define latn(j) (get_linspace_val(j,SW->lat0,SW->lat1,SW->ny)) // Planetographic latitude of the staggered y/north point. Latitudes can only be staggered in y
#define latc(j) ((latn(j)+latn(j-1))/2.) // Planetographic latitude of the centered point. The centered point falls in the middle of staggered points
#define latv(j) (latn(j))

// Latitude and longitude matrices (correlation table for a structured grid)
#define LONE(i,j) lone(i) // Longitude of a staggered x/east point
#define LONN(i,j) lonc(i) // Longitude of a staggered y/north point
#define LONC(i,j) lonc(i) // Longitude of a centered point
#define LONV(i,j) lonv(i) // Longitude of a vertex point
#define LATE(i,j) latc(j) // Planetographic latitude of a  staggered x/east point
#define LATN(i,j) latn(j) // Planetographic latitude of a staggered y/north point
#define LATC(i,j) latc(j) // Planetographic latitude of a centered point
#define LATV(i,j) latv(j) // Planetographic latitude of a vertex point

#define POLAR_LATC(i,j) (SW->polar)*(PI/2.0-get_great_angle(LONC(i,j),LATC(i,j),MEAN(SW->lon0,SW->lon1),MEAN(SW->lat0,SW->lat1))) //Latitude with respect to the pole C
#define POLAR_LATE(i,j) (SW->polar)*(PI/2.0-get_great_angle(LONE(i,j),LATE(i,j),MEAN(SW->lon0,SW->lon1),MEAN(SW->lat0,SW->lat1))) //Latitude with respect to the pole E
#define POLAR_LATN(i,j) (SW->polar)*(PI/2.0-get_great_angle(LONN(i,j),LATN(i,j),MEAN(SW->lon0,SW->lon1),MEAN(SW->lat0,SW->lat1))) //Latitude with respect to the pole N
#define POLAR_LATV(i,j) (SW->polar)*(PI/2.0-get_great_angle(LONV(i,j),LATV(i,j),MEAN(SW->lon0,SW->lon1),MEAN(SW->lat0,SW->lat1))) //Latitude with respect to the pole V


// #define DLONE(i,j) (LONE(i,j)-LONE(i-1,j))
// #define DLONN(i,j) (LONN(i,j)-LONN(i-1,j))
// #define DLONC(i,j) (LONC(i,j)-LONC(i-1,j))
// #define DLONV(i,j) (LONV(i,j)-LONV(i-1,j))
// #define DLATE(i,j) (LATE(i,j)-LATE(i,j-1))
// #define DLATN(i,j) (LATN(i,j)-LATN(i,j-1))
// #define DLATC(i,j) (LATC(i,j)-LATC(i,j-1))
// #define DLATV(i,j) (LATV(i,j)-LATV(i,j-1))

// Coordinates of the centered points and staggered points
// Even though a plane or a spheroid do not require to store the coordinates of each point in-memory (a vector for each axis would suffice), because the core of the solver has been developed to be as generalized as possible and to avoid expensive operations (in terms of computation) such as for an ellipsoid, a provisional decision to store the information of every point in-memory was made. A great deal of memory could be freed if the solver will only solve revolution objects
#define DXE(i,j) ac(SW->Dxe,i,j,SW->M) // Distance on x between staggered x/east points
#define DXC(i,j) ac(SW->Dxc,i,j,SW->M) // Distance on x between centered points
#define DXV(i,j) ac(SW->Dxv,i,j,SW->M) // Distance on x between vertex points
#define DYN(i,j) ac(SW->Dyn,i,j,SW->M) // Distance on y between staggered y/north points
#define DYC(i,j) ac(SW->Dyc,i,j,SW->M) // Distance on y between centered points
#define DYV(i,j) ac(SW->Dyv,i,j,SW->M) // Distance on y between vertex points
#define XE(i,j) ac(SW->xe,i,j,SW->M)   // x coordinates of a staggered x/east point
#define XC(i,j) ac(SW->xc,i,j,SW->M)   // x coordinates of a centered point
#define XV(i,j) ac(SW->xv,i,j,SW->M)   // x coordinates of a vertex point
#define YN(i,j) ac(SW->yn,i,j,SW->M)   // y coordinates of a staggered y/north point
#define YC(i,j) ac(SW->yc,i,j,SW->M)   // y coordinates of a centered point
#define YV(i,j) ac(SW->yv,i,j,SW->M)   // y coordinates of a vertex point

// Spatial 3d (not on surface but in a Cartesian 3d grid) coordinates of centered and staggered points
#define XC3D(i,j) r(lat2pclat(LATC(i,j)))*cos(LONC(i,j))*cos(LATC(i,j)) // Spatial x coordinate of centered point in Cartesian 3d grid
#define XV3D(i,j) r(lat2pclat(LATV(i,j)))*cos(LONV(i,j))*cos(LATV(i,j)) // Spatial y coordinate of centered point in Cartesian 3d grid
#define YC3D(i,j) r(lat2pclat(LATC(i,j)))*sin(LONC(i,j))*cos(LATC(i,j)) // Spatial z coordinate of centered point in Cartesian 3d grid
#define YV3D(i,j) r(lat2pclat(LATV(i,j)))*sin(LONV(i,j))*cos(LATV(i,j)) // Spatial x coordinate of vertex point (nodes) in Cartesian 3d grid
#define ZC3D(i,j) r(lat2pclat(LATC(i,j)))*sin(LATC(i,j))                // Spatial y coordinate of vertex point (nodes) in Cartesian 3d grid
#define ZV3D(i,j) r(lat2pclat(LATV(i,j)))*sin(LATV(i,j))                // Spatial z coordinate of vertex point (nodes) in Cartesian 3d grid

// #define ds_dpclat(pclat) ((SW->rP)*sqrt(1-(1-(SW->epsilon2))*POW2(sin(pclat)))) // Same as rM
#define pclat2lat(pclat) atan((SW->epsilon2)*tan(pclat)) // Planetocentric latitude to planetographic latitude operator
#define lat2pclat(lat) atan(1.0/(SW->epsilon2)*tan(lat)) // Planetographic latitude to planetocentric latitude operator
#define fEP (((SW->rE)-(SW->rP))/(SW->rP)) // Second flattening of the ellipse
#define r(pclat) ((SW->rE)*(SW->rP)/sqrt(POW2(SW->rE*sin(pclat))+POW2(SW->rP*cos(pclat)))) // Planetocentric radius of the spheroid
// #define r(pclat) ((SW->rE)*(1.0-fEP*POW2(sin(pclat)))) // An approximation of \mymacro{r}
#define rZ(lat) ((SW->rE)/sqrt(1.0+POW2(tan(lat))/(SW->epsilon2))) // Zonal radius of the spheroid
#define rM(lat) ((SW->rE)/(SW->epsilon2)*POW3(rZ(lat)/((SW->rE)*cos(lat)))) // Meridional radius of the spheroid

// Gravity definition
#define g0(pclat) ((SW->gE + POW2((SW->Omega))*(SW->rE))*POW2((SW->rE)/(r(pclat)))) // Reference gravity (potential)
#define g_eff(pclat) (g0(pclat)-POW2((SW->Omega))*r(pclat)*cos(pclat)) // Effective gravity of the rotating spheroid

// Gravitiy matrices
#define GEmacro(i,j) g_eff(lat2pclat(LATE(i,j))) // Effective gravity at a staggered x/east point
#define GNmacro(i,j) g_eff(lat2pclat(LATN(i,j))) // Effective gravity at a staggered y/north point
#define GCmacro(i,j) g_eff(lat2pclat(LATC(i,j))) // Effective gravity at a centered point
#define GVmacro(i,j) g_eff(lat2pclat(LATV(i,j))) // Effective gravity at a vertex point

#define GE(i,j) ac(SW->ge,i,j,SW->M) // Effective gravity at a staggered x/east point
#define GN(i,j) ac(SW->gn,i,j,SW->M) // Effective gravity at a staggered y/north point
#define GC(i,j) ac(SW->gc,i,j,SW->M) // Effective gravity at a centered point
#define GV(i,j) ac(SW->gv,i,j,SW->M) // Effective gravity at a vertex point

// Definition of the Coriolis parameter
#define f(lat) (2.0*(SW->Omega)*sin(lat)) // Coriolis parameter

// Matrices with the Coriolis parameter
#define FEmacro(i,j) f(LATE(i,j)) // Coriolis parameter at a staggered x/east point
#define FNmacro(i,j) f(LATN(i,j)) // Coriolis parameter at a staggered y/north point
#define FCmacro(i,j) f(LATC(i,j)) // Coriolis parameter at a centered point
#define FVmacro(i,j) f(LATV(i,j)) // Coriolis parameter at a staggered y/north point

#define FE(i,j) ac(SW->fe,i,j,SW->M) // Coriolis parameter at a staggered x/east point
#define FN(i,j) ac(SW->fn,i,j,SW->M) // Coriolis parameter at a staggered y/north point
#define FC(i,j) ac(SW->fc,i,j,SW->M) // Coriolis parameter at a centered point
#define FV(i,j) ac(SW->fv,i,j,SW->M) // Coriolis parameter at a staggered y/north point

// Main matrices
#define TRACER(i,j) ac(SW->tracer,i,j,SW->M) // Tracer at centered point
#define PERTURBS(i,j) ac(SW->perturbs,i,j,SW->M) // Perturbation (0:no, 1: yes) at centered point
#define ETA(i,j) ac(SW->eta,i,j,SW->M) // The surface perturbation at a centered point
#define U(i,j) ac(SW->u,i,j,SW->M) // The horizontal velocity at a staggered x/east point
#define V(i,j) ac(SW->v,i,j,SW->M) // The vertical velocity at a staggered y/north point
#define HB(i,j) ac(SW->hB,i,j,SW->M) // The surface relief at a centered point
#define H(i,j) ((ETA(i,j))+(SW->depth)-(HB(i,j))) // READ ONLY. The layer depth at a centered point
#define WINDX(i,j) ac(SW->windx,i,j,SW->M) // Wind in x direction at a staggered x/east point.
#define WINDY(i,j) ac(SW->windy,i,j,SW->M) // Wind in y direction at a staggered y/north point.
#define UWINDX(i,j) (U(i,j)+WINDX(i,j)) // READ ONLY. The sum of the WINDX and the horizontal velocity U at a staggered x/east point
#define VWINDY(i,j) (V(i,j)+WINDY(i,j)) // READ ONLY. The sum of the WINDY and the vertical velocity V at a staggered y/north point

// PROTOTYPES
int wrapper(sw *SW);

int initialize(int npx,int npy,int force1proc,char *bname,int copyinput,sw *SW);
void loop(char *bname, char *output_folder, sw *SW);
void destroy(sw *SW);

void init_coords_ortho(sw *SW);
void init_coords_spheroid(sw *SW);

// numeric integration functions
void solve_tracer(double *Dtracer_n, double Dt, sw *SW);
void get_pvort(sw *SW);
void solve_cons_h(double *Deta_n, double Dt, sw *SW);
void do_eta_geofactor(double *Deta_nAB, double Dt, sw *SW);
void solve_adv_u(double *Du_adv, double Dt, sw *SW);
void solve_adv_v(double *Dv_adv, double Dt, sw *SW);
void solve_pres_u(double *Du_pres, double Dt, sw *SW);
void solve_pres_v(double *Dv_pres, double Dt, sw *SW);
void solve_Coriolis_u(double *Du_nAB, double Dt, sw *SW);
void solve_Coriolis_v(double *Dv_nAB, double Dt, sw *SW);
void do_channel(sw *SW);
void do_AdamsBashforth(double *Dvar_nAB, double *Dvar_n, double *Dvar_nM1, double *Dvar_nM2, int n, map *M);
double compute_Dt(double Courant, sw *SW);
void add_Gaussians(double Dt, double t, sw *SW);
void add_Gaussian_core(double lon, double lat, double volume, double sigma, double radius, double inject_tracer, sw *SW);
void get_Gaussian_field(double lon, double lat, double amp, double sigma, double radius, double *field, double inject_tracer, sw *SW);
void add_vortices(double Dt, double t, sw *SW);
void add_vortex_core(double *eta, double *u, double *v, double n, double amp, double loncent, double latcent, double as, double bs, sw *SW);

double linterp_core(double x_i, double *x, double *y, int *i0, int *i1);
double linterp(double x_i, double *x, double *y, int n);
void linterp_vec(double *x_i, double *y_i, int n_i, double *x, double *y, int n);
double get_great_circle(double lon0, double lat0, double lon1, double lat1, double r);
double get_vincenty(double lon0, double lat0, double lon1, double lat1, double a, double b);
double get_arc_ellipse(double a1, double rP, double epsilon2);
double get_arc_circumf(double a1, double r);
double get_linspace_val(int i,double x0,double x1,int nx);
void generate_scaf(double *scaf, double (*funcio)(double,double,double), double x0, double y0, int stgx, int stgy, double t, sw *SW);
double generate_val(int i, int j, double (*funcio)(double,double,double), double x0, double y0, int stgx, int stgy, double t, sw *SW);

// Geostrophic eq
void init_scaf_eta(double *eta, double n, double amp, double loncent, double latcent, double as, double bs, sw *SW);
void init_scaf_u(double *u, double *eta, double n, double loncent, double latcent, double as, double bs, sw *SW);
void init_scaf_v(double *v,  double *eta, double n, double loncent, double latcent, double as, double bs, sw *SW);

// Sponge and polar regions
void dissipation_effect(double *Du_n, double *Dv_n, sw *SW);
double get_great_angle(double lon0, double lat0, double lon1, double lat1);
void init_polar_region(sw *SW);
void init_channel_region(sw *SW);


void read_wind(FILE *fin, sw *SW);
void read_Gaussians(FILE *fin, sw *SW);
void read_vortices(FILE *fin, sw *SW);


void calc_energy_balance(double *ke,double *ape,sw *SW);
void print_energy(double t,double ke,double ape); // just print

// binary IO
void save_binary_sw(int frame, int nstep,double timedouble,char *bname, char *output_folder,sw *SW);
int load_binary_sw(int frame, int *nstep,double *timedouble,char *bname, char *binaries_folder,sw *SW);
int load_binary_header_sw(int frame,int *nstep,double *timedouble, char *bname,char *binaries_folder);

// reference fields
double compare_ref(sw *SW,int n,double t, char *bname,char *output_folder);



//Hyperviscosity
void laplacian_vels(double *u_in, double *v_in, double  *u_out, double *v_out, sw *SW);
void laplacian_eta(double *e_in, double *e_out, sw *SW);


// double funPsi(double r);
void select_scheme(int s);
int query_scheme(int s);

/*
  RECORD OF CHANGES:
  - 113dbg3:
      o Change macros of F and G for scalar fields
  - 116:
      o Addition of polar regions



  TESTS:
  - 113dbg3:
      o Minitest canal with winds. Results:
          -> With CHECK activated, no errors are displayed
          -> When compared to 113dbg2:
             # Small mesh, short simualtion: no difference observed in results, yet yes in speed (x2.5 faster)
             # Dense mesh, short simulation: no difference observed, yet yes in speed (x2 faster)


*/
