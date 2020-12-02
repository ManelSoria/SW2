// sppde: structured parallel pde solver
// Manel Soria 2017 - UPC - ESEIAAT - TUAREG


// Use #define CHECK to activate matrix range check. This makes the code slow but safe
// Don't use it for production runs
// In case of problems, recompile and test the code with CHECK activated
//#define CHECK // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


// MPI interfase =========================================
void checkr(int r, char *txt); // to check if an MPI call has been sucessful
int quisoc(void); // my processor number, from 0 to quants()
int quants(void); // number of processors

// basic parallel IO and error handling
#define CRASH(...) crash(__LINE__,__FILE__,__VA_ARGS__) // halts nicely the program

char* poutp(); // returns path of the output files (used by pprintf)
int pprintf(char *fmt, ...); // subtitutes printf, redirecting output to a file if needed
void init_pprintf(int mode, int flushall, char *path); // inits io with pprintf
// mode:
// 0: all processors to file
// 1: processor 0 to stdout, the rest to the file
// 2: processor 0 to stdout, the rest ignore output

// flushall: files are flushed after every pprintf (very slow)
// path: path to create stdout files (if needed)

void end_pprintf(void); // ends output
void crash(int line, char *file, char *fmt, ...); // (internal) auxiliary of CRASH

// memory management
// memory map =====================
#define MAXP 1920 // max number of procs, change to increase
// axis ranges from 1 to 2 (3 if 3d)
typedef struct { // never use this directly except in sppde
	int sh; // halo size (>=1)
	int nd; // number of dimensions2 or 3
	int l0[4]; // start of my owned area, NOT HALO, directions [1..3] (0 is not used)
 	int l1[4]; // end
  	int gl0[4];  // IDEM, global domain
  	int gl1[4];
  	// distribution of data for other processors
  	int al0[MAXP][4]; // al0[p][a] start area proc=p, axis=a
  	int al1[MAXP][4]; // al1[p][a] end   area proc=p, axis=a
  	int pco[MAXP][4]; // processor corrdinates [p][a]

	int s[4]; // size of my area WITH HALO, at each axis
	int S; // total size of my area
 	int gpc[4]; // number of processors at each axis
	int pc[4]; // my coordinates in the processor grid  (starting at 1,1,1)
	int nb_e,nb_w,nb_n,nb_s; // neighbour processor at each direction (-1 if any)
	int nb_t,nb_b; // top, botom (3d)
	int bc0[4]; // Do I own at least a part of the BC in initial semiaxis [a] ? (a=1,2,3)
	int bc1[4]; // Do I own at least a part of the BC in final semiaxis [a] ? (a=1,2,3)
	int per[4]; // Periodicity in each axis
} map;

int crax(int i, int a, int line, char *file, map *M); // internal

void *mallocc(int n); // malloc or crash
void *callocc(int n, int s) ; // calloc or crash

#define check(q,a,M)  ( ( (q)>=(M->l0[a])-M->sh && (q)<=(M->l1[a])+M->sh ) ? (q) : crax(q,a,__LINE__,__FILE__,M) )



// ==================== macros to query the map --- NEVER use directly the structure !

#define HS(M) ((M)->sh) // halo size

// limits of owned area (NOT halos)
#define lsx(M) ((M)->l0[1]) // local start , x axis
#define lsy(M) ((M)->l0[2]) //             , y axis
#define lsz(M) ((M)->l0[3]) //             , z axis

#define lex(M) ((M)->l1[1]) // local end   , x axis
#define ley(M) ((M)->l1[2]) //             , y axis
#define lez(M) ((M)->l1[3]) //             , z axis

// limits of global domain (NOT external halos)
#define gsx(M) ((M)->gl0[1]) // global start, x axis
#define gsy(M) ((M)->gl0[2]) //             , y axis
#define gsz(M) ((M)->gl0[3]) //             , z axis

#define gex(M) ((M)->gl1[1]) // global end  , x axis
#define gey(M) ((M)->gl1[2]) //             , y axis
#define gez(M) ((M)->gl1[3]) //             , z axis

// absolute limites of the universe (WITH external halos)
#define asx(M) ((M)->gl0[1]-M->sh) // global start, x axis
#define asy(M) ((M)->gl0[2]-M->sh) //             , y axis
#define asz(M) ((M)->gl0[3]-M->sh) //             , z axis

#define aex(M) ((M)->gl1[1]+M->sh) // global end  , x axis
#define aey(M) ((M)->gl1[2]+M->sh) //             , y axis
#define aez(M) ((M)->gl1[3]+M->sh) //             , z axis

#define NPX(M) ((M)->gpc[1]) // size processor grid, x axis
#define NPY(M) ((M)->gpc[2]) //                      y axis
#define NPZ(M) ((M)->gpc[3]) //                      z axis

#define PGX(M) ((M)->gpc[1]) // my processor grid coordinate, x axis (1..SPGX)
#define PGY(M) ((M)->gpc[2]) //                               y axis
#define PGZ(M) ((M)->gpc[3]) //                               z axis

#define pex(M) ((M)->per[1]) // is the domain periodic, x axis ?
#define pey(M) ((M)->per[2]) //             , y axis
#define pez(M) ((M)->per[3]) //             , z axis

// manel
#define obN(M) ((M)->bc1[2]) // do I own a part of North BC ? (returns 0 if periodic in Y)
#define obS(M) ((M)->bc0[2]) // south ?
#define obE(M) ((M)->bc1[1]) // east ?
#define obW(M) ((M)->bc0[1]) // west ?

// =========================================== 2d ==================
// =========================================== 2d ==================
// 2d access macros

#ifdef CHECK

// evaluate the position in the vector (internal)
#define LP(i,j,M) (    check(i,1,(M))-((M)->l0[1])+(M)->sh + \
                      (check(j,2,(M))-((M)->l0[2])+(M)->sh )*((M)->s[1])   )
#else
#define LP(i,j,M) ( (i)- ((M)->l0[1])+(M)->sh + \
                      ((j)-((M)->l0[2])+(M)->sh )*((M)->s[1])  )
#endif

// acces to a generic 2d scalar field
#define ac(v,i,j,M) *((v)+LP(i,j,M))

// macros for 1d looping  (internal) ======================================

// only inner nodes:
#define for1(i,a,M) for (i=(M)->l0[a]; i<=(M)->l1[a] ; i++)

// inner nodes plus halo
#define for1h(i,a,M) for (i=(M)->l0[a]-(M)->sh; i<=(M)->l1[a]+(M)->sh ; i++)

// macros for 2d looping (this is the standard way to operate with fields) =========

// loop interior, NO HALOS, 2d (this is be the most frequent loop)
#define forall(i,j,M) for1(j,2,M) for1(i,1,M)

// loop incloent halos, 2d
#define forallh(i,j,M) for1h(j,2,M) for1h(i,1,M)

// loop my inner nodes (not halos), EXCEPT external boundaries
// periodic boundary conditions are not looped
#define foralli(i,j,M) \
 for (j=(M)->l0[2]+((M)->bc0[2] == 1); j<=(M)->l1[2]-((M)->bc1[2] == 1); j++) \
 for (i=(M)->l0[1]+((M)->bc0[1] == 1); i<=(M)->l1[1]-((M)->bc1[1] == 1); i++)

// loop inner nodes + south + west halos (even if they are in external halos).. not very common
#define for1Xnodes(i,a,M) for (i=(M)->l0[a]-1; i<=(M)->l1[a] ; i++)
#define forallX(i,j,M) for1Xnodes(j,2,M) for1Xnodes(i,1,M)

// access to a position only if owned
#define ifowned(i,j,M) if ((i)>=(M)->l0[1] && (i)<=(M)->l1[1] && (j)>=(M)->l0[2] && (j)<=(M)->l1[2])

// access to a position only if owned or in my halos
#define ifownedh(i,j,M) if ((i)>=(M)->l0[1]-(M)->sh && (i)<=(M)->l1[1]+(M)->sh && (j)>=(M)->l0[2]-(M)->sh && (j)<=(M)->l1[2]+M->sh)

// sweep only certain boundaries

// exclude vertices if non-periodic
#define fornorth(i,j,M)  j=(M)->l1[2]; if ((M)->bc1[2]) for (i=(M)->l0[1]+(M)->bc0[1];i<=(M)->l1[1]-(M)->bc1[1];i++)
#define foreast(i,j,M)   i=(M)->l1[1]; if ((M)->bc1[1]) for (j=(M)->l0[2]+(M)->bc0[2];j<=(M)->l1[2]-(M)->bc1[2];j++)
#define forsouth(i,j,M)  j=(M)->l0[2]; if ((M)->bc0[2]) for (i=(M)->l0[1]+(M)->bc0[1];i<=(M)->l1[1]-(M)->bc1[1];i++)
#define forwest(i,j,M)   i=(M)->l0[1]; if ((M)->bc0[1]) for (j=(M)->l0[2]+(M)->bc0[2];j<=(M)->l1[2]-(M)->bc1[2];j++)

// just vertices
#define forne(i,j,M) i=(M)->gl1[1]; j=(M)->gl1[2]; if ((M)->bc1[1] && (M)->bc1[2]) ifowned(i,j,M)
#define fornw(i,j,M) i=(M)->gl0[1]; j=(M)->gl1[2]; if ((M)->bc0[1] && (M)->bc1[2]) ifowned(i,j,M)
#define forse(i,j,M) i=(M)->gl1[1]; j=(M)->gl0[2]; if ((M)->bc1[1] && (M)->bc0[2]) ifowned(i,j,M)
#define forsw(i,j,M) i=(M)->gl0[1]; j=(M)->gl0[2]; if ((M)->bc0[1] && (M)->bc0[2]) ifowned(i,j,M)

// end macros for 2d looping ===================================


// some macros to access generic fields X Y Z
#define X(i,j,M) *((x)+LP(i,j,M))
#define Y(i,j,M) *((y)+LP(i,j,M))
#define Z(i,j,M) *((z)+LP(i,j,M))


// =========================================== 3d ==================
// =========================================== 3d ==================

// 3d access macros
#ifdef CHECK
#define LP3(i,j,k,M) ((check(i,1,M)-(M->l0[1])+M->sh)+\
                      (check(j,2,M)-(M->l0[2])+M->sh)*(M->s[1])+\
                      (check(k,3,M)-(M->l0[3])+M->sh)*(M->s[1])*(M->s[2])\
                      )
#else
#define LP3(i,j,k,M) (((i)-(M->l0[1])+M->sh)+\
                      ((j)-(M->l0[2])+M->sh)*(M->s[1])+\
                      ((k)-(M->l0[3])+M->sh)*(M->s[1])*(M->s[2])\
				  ) // QRAT
#endif

// acces 3d a un camp escalar generic v
#define ac3(v,i,j,k,M) *((v)+LP3(i,j,k,M))


#define X3(i,j,k,M) *((x)+LP3(i,j,k,M))
#define Y3(i,j,k,M) *((y)+LP3(i,j,k,M))
#define Z3(i,j,k,M) *((z)+LP3(i,j,k,M))

#define Ifowned(i,j,k,M) if ( \
	(i)>=M->l0[1] && (i)<=M->l1[1] && \
	(j)>=M->l0[2] && (j)<=M->l1[2] && \
	(k)<=M->l0[3] && (k)<=M->l1[3] )

// loop interior, NO HALOS, 3d
#define Forall(i,j,k,M) for1(k,3,M) for1(j,2,M) for1(i,1,M)

// loop incloent halos, 3d (no es el mes usual)
#define Forallh(i,j,k,M) for1h(k,3,M) for1h(j,2,M) for1h(i,1,M)

// loop my inner nodes (not halos), EXCEPT external boundaries
#define Foralli(i,j,k,M) \
 for (k=M->l0[3]+(M->bc0[3] == 1); k<=M->l1[3]-(M->bc1[3] == 1); k++) \
 for (j=M->l0[2]+(M->bc0[2] == 1); j<=M->l1[2]-(M->bc1[2] == 1); j++) \
 for (i=M->l0[1]+(M->bc0[1] == 1); i<=M->l1[1]-(M->bc1[1] == 1); i++)

// these macros sweep all the lateral faces including the edges and vertices
#define Foreast(i,j,k,M)   i=M->l1[1]; if (M->bc1[1]) for1(k,3,M) for1(j,2,M)
#define Fornorth(i,j,k,M)  j=M->l1[2]; if (M->bc1[2]) for1(k,3,M) for1(i,1,M)
#define Fortop(i,j,k,M)    k=M->l1[3]; if (M->bc1[3]) for1(j,2,M) for1(i,1,M)

#define Forwest(i,j,k,M)   i=M->l0[1]; if (M->bc0[1]) for1(k,3,M) for1(j,2,M)
#define Forsouth(i,j,k,M)  j=M->l0[2]; if (M->bc0[2]) for1(k,3,M) for1(i,1,M)
#define Forbottom(i,j,k,M) k=M->l0[3]; if (M->bc0[3]) for1(j,2,M) for1(i,1,M)

// these macros sweep the lateral faces but not the edges ============================= 	 TO BE TESTED

#define for1I(i,a,M) for (i=M->l0[a]+M->bc0[a]; i<=M->l1[a]-M->bc1[a] ; i++)

#define ForE(i,j,k,M) i=M->l1[1]; if (M->bc1[1]) for1I(k,3,M) for1I(j,2,M)
#define ForW(i,j,k,M) i=M->l0[1]; if (M->bc0[1]) for1I(k,3,M) for1I(j,2,M)

#define ForN(i,j,k,M) j=M->l1[2]; if (M->bc1[2]) for1I(k,3,M) for1I(i,1,M)
#define ForS(i,j,k,M) j=M->l0[2]; if (M->bc0[2]) for1I(k,3,M) for1I(i,1,M)

#define ForT(i,j,k,M) k=M->l1[3]; if (M->bc1[3]) for1I(j,2,M) for1I(i,1,M)
#define ForB(i,j,k,M) k=M->l0[3]; if (M->bc0[3]) for1I(j,2,M) for1I(i,1,M)

// these macros sweep the edges but not the vertices ============================= 	 TO BE TESTED
// axis 1 and 2

// Work for AnBo: all these loops have to be simplified with the new for1I, just like ForEN
// then test the sweeps with the functions in test_core
//#define ForEN(i,j,k,M) i=M->l1[1]; j=M->l1[2]; if (M->bc1[1] && M->bc1[2]) for (k=M->l0[3]+M->bc0[3]; k<=M->l1[3]-M->bc1[3] ; k++)

#define ForEN(i,j,k,M) i=M->l1[1]; j=M->l1[2]; if (M->bc1[1] && M->bc1[2]) for1I(k,3,M)

#define ForES(i,j,k,M) i=M->l1[1]; j=M->l0[2]; if (M->bc1[1] && M->bc0[2]) for (k=M->l0[3]+M->bc0[3]; k<=M->l1[3]-M->bc1[3] ; k++)
#define ForWN(i,j,k,M) i=M->l0[1]; j=M->l1[2]; if (M->bc0[1] && M->bc1[2]) for (k=M->l0[3]+M->bc0[3]; k<=M->l1[3]-M->bc1[3] ; k++)
#define ForWS(i,j,k,M) i=M->l0[1]; j=M->l0[2]; if (M->bc0[1] && M->bc0[2]) for (k=M->l0[3]+M->bc0[3]; k<=M->l1[3]-M->bc1[3] ; k++)

// axis 2 and 3
#define ForNT(i,j,k,M) j=M->l1[2]; k=M->l1[3]; if (M->bc1[2] && M->bc1[3]) for (i=M->l0[1]+M->bc0[1]; i<=M->l1[1]-M->bc1[1] ; i++)
#define ForNB(i,j,k,M) j=M->l1[2]; k=M->l0[3]; if (M->bc1[2] && M->bc0[3]) for (i=M->l0[1]+M->bc0[1]; i<=M->l1[1]-M->bc1[1] ; i++)
#define ForST(i,j,k,M) j=M->l0[2]; k=M->l1[3]; if (M->bc0[2] && M->bc1[3]) for (i=M->l0[1]+M->bc0[1]; i<=M->l1[1]-M->bc1[1] ; i++)
#define ForSB(i,j,k,M) j=M->l0[2]; k=M->l0[3]; if (M->bc0[2] && M->bc0[3]) for (i=M->l0[1]+M->bc0[1]; i<=M->l1[1]-M->bc1[1] ; i++)

// axis 1 and 3
#define ForET(i,j,k,M) i=M->l1[1]; k=M->l1[3]; if (M->bc1[1] && M->bc1[3]) for (j=M->l0[2]+M->bc0[2]; j<=M->l1[2]-M->bc1[2] ; j++)
#define ForEB(i,j,k,M) i=M->l1[1]; k=M->l0[3]; if (M->bc1[1] && M->bc0[3]) for (j=M->l0[2]+M->bc0[2]; j<=M->l1[2]-M->bc1[2] ; j++)
#define ForWT(i,j,k,M) i=M->l0[1]; k=M->l1[3]; if (M->bc0[1] && M->bc1[3]) for (j=M->l0[2]+M->bc0[2]; j<=M->l1[2]-M->bc1[2] ; j++)
#define ForWB(i,j,k,M) i=M->l0[1]; k=M->l0[3]; if (M->bc0[1] && M->bc0[3]) for (j=M->l0[2]+M->bc0[2]; j<=M->l1[2]-M->bc1[2] ; j++)

// these macros "sweep" the vertices === TO BE TESTED
#define ForWSB(i,j,k,M) i=M->l0[1]; j=M->l0[2]; k=M->l0[3]; if (M->bc0[1] && M->bc0[2] && M->bc0[3])
#define ForWST(i,j,k,M) i=M->l0[1]; j=M->l0[2]; k=M->l1[3]; if (M->bc0[1] && M->bc0[2] && M->bc1[3])
#define ForWNB(i,j,k,M) i=M->l0[1]; j=M->l1[2]; k=M->l0[3]; if (M->bc0[1] && M->bc1[2] && M->bc0[3])
#define ForWNT(i,j,k,M) i=M->l0[1]; j=M->l1[2]; k=M->l1[3]; if (M->bc0[1] && M->bc1[2] && M->bc1[3])
#define ForESB(i,j,k,M) i=M->l1[1]; j=M->l0[2]; k=M->l0[3]; if (M->bc1[1] && M->bc0[2] && M->bc0[3])
#define ForEST(i,j,k,M) i=M->l1[1]; j=M->l0[2]; k=M->l1[3]; if (M->bc1[1] && M->bc0[2] && M->bc1[3])
#define ForENB(i,j,k,M) i=M->l1[1]; j=M->l1[2]; k=M->l0[3]; if (M->bc1[1] && M->bc1[2] && M->bc0[3])
#define ForENT(i,j,k,M) i=M->l1[1]; j=M->l1[2]; k=M->l1[3]; if (M->bc1[1] && M->bc1[2] && M->bc1[3])

// =========================================== 3d ==================
// =========================================== 3d ==================


// ========== macros for 2d/3d.. k is ignored where 2d

#define gac(v,i,j,k,M) *( (v) + ( (M->nd==3) ? (LP3(i,j,k,M)) : (LP(i,j,M)) ) )

#define gOwned(i,j,k,M) \
 	(  (i)>=M->l0[1] && (i)<=M->l1[1] && \
	   (j)>=M->l0[2] && (j)<=M->l1[2] && \
	 (((k)>=M->l0[3] && (k)<=M->l1[3]) || M->nd==2) )

double *dmem(map *M); // allocates scalar field
void printmap(char *name,map *M); // prints map
void createM(int nd,int *gl0,int *gl1,int *np,int sh,int *per,map *M); // creates map
void createM_2d(int npx,int npy,int gsx,int gex,int gsy,int gey,int sh,int perx,int pery,map *M); // idem, 2d, "easy"

void joinM(map *M, map *MJ);



void halo_update(double *x,map *M);


double *alloc_gather(double *x,map *M, map *MJ,int eh); // allocates global map MJ and gathers all the data


// this will replace gather_field
// given two fields x and y
//   x is distributed x according to map M
//   y is global, according to MJ (consistent with M), only available to proc 0
// m==0: gather  proc 0 receives all the field x and stores in the preallocated global field y
// m==1: scatter proc 0 sends global field y to the procs, that store it in their prealloc areas of x

// eh controls the behaviour of the halos
// if eh==0, halos are not gathered nor scattered
// if eh==1, external halos are gathered to 0, and **all** the halos are scattered
void gs_field(int m,int eh,double *x,map *M,double *y,map *MJ);


// field operators ========================================
void print_scaf(double *x, char *label, int h, map *M, char *fmt);
void setzero_scaf(double *r, map *M);
void linop_scaf(double *r, double *a, double *b, double k1, double k2, map *M);
double norm_scaf(double *r,int l,map *M); // norm of r, l=0: Max(Abs(r)), otherwise Sum(r^2)
double norm_dif_scaf(double *a,double *b,int l,map *M); // norm of a-b
void copy_scaf(double *matrixCopy, double *matrixInput, int copy_halo, map *M); // Copies the matrix matrixInput to matrixCopy

// field statistics (for debugging mainly)
void get_stats_scaf(double *r,map *M,char *fname,double *gmin_,double *gmax_,double *gavg_); // eval and return (not print)
void stats_scaf(double *r,map *M,char *fname); // eval and print
void check_stats_scaf(double *r,map *M,char *fname,int step,double MIN, double MAX,int prn); // eval, print (if prn==1) and crash if out of range

// time control functions

void cr_reset(void);
/* modes de comptar el temps d'execucio */
/*  0 -> Temps fisic. / VALOR PER DEFECTE /
    1 -> Temps d'usuari + sistema
    2 -> Temps d'usuari                 */
void cr_setmode (int mode);
void cr_info(void);
void cr_start(char *n,int suff);
void cr_end(char *n,int suff);
double cr_time(char *n,int suff);


// binary IO

void fwrite_scaf(double *x, FILE *f,map *M);
void fread_scaf(double *x, FILE *f,map *M);
void text_write_scaf(double *x, char *fname,map *M); // text write actually
void text_read_scaf(double *x, char *fname, map *M);


// generates an auxiliary field for testing
void set_field(double *x,int f,map *M);
