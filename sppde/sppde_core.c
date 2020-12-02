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
#include <sys/time.h>
#include <sys/times.h>
#include <stdio.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "mpi.h"
#include "sppde.h"

static void split1d(int ini, int fi, int NP, int *i, int *f, int minsize); // internal use only

// mpi =========================

void checkr(int r,char *txt) {
  if (r!=MPI_SUCCESS) {
	int a,b;
  	a=MPI_Comm_rank(MPI_COMM_WORLD,&b);
    crash(__LINE__,__FILE__,"MPI error proc=%d error code=%s \n",a,txt);
    exit(-1);
  }
}

int quisoc() {
  int a,b;
  a=MPI_Comm_rank(MPI_COMM_WORLD,&b);
  checkr(a,"quisoc");
  return(b);
}

int quants() {
  int a,b;
  a=MPI_Comm_size(MPI_COMM_WORLD,&b);
  checkr(a,"quants");
  return(b);
}

void *mallocc(int n) {
  void *x;
  x=malloc(n);
  if (x==NULL) CRASH("mallocc n=%d",n);
  return(x);
}

void *callocc(int n, int s) {
  void *x;
  x=calloc(n,s);
  if (x==NULL) CRASH("callocc n=%d",n);
  return(x);
}

FILE *myout=NULL;
int flushall=0;
int initout=0;

int pprintf(char *fmt,...) {
  if (initout==0) { fprintf(stderr,"initout=%d call init_pprintf before pprintf !!\n",initout); MPI_Finalize(); exit(1); }
  int r;
  va_list ap;
  va_start(ap,fmt);
  r=vfprintf(myout,fmt,ap);
  va_end(ap);
  if (flushall==1)
    fflush(myout);
  return(r);
}

// 0: all processors to file
// 1: processor 0 to stdout, the rest to the file
// 2: processor 0 to stdout, the rest ignore output

#define MAXL 1000

char outpath_[MAXL];

char* poutp() { // returns output path
  return outpath_;
}

void init_pprintf(int mode,int flushall_, char *path) {
  char q[MAXL];
  int openfile=0,tonull=0;
  int info=0;

  if (info) printf("init_pprintf mode=%d \n",mode);
  flushall=flushall_;

  if (quisoc()==0) {
    if (mode==0) openfile=1;
  }
  else {
    if (mode==0 || mode==1) openfile=1;
  }

  if (mode==2 && quisoc()!=0) tonull=1;


  if (info) printf("init_pprintf openfile=%d\n",openfile);

  // create output folder, even if proc 0 won't print to it
  if (quisoc()==0) {
      struct stat st = {0};
      if (info) fprintf(stdout,"init_pprintf stat=%d \n",stat(path, &st) );
      if (stat(path, &st) == -1) {
        sprintf(q,"mkdir -p %s",path);
        fprintf(stdout,"Processor 0 creating folder '%s' for output\n",path);
        int r=system(q);
        if (r!=0) { fprintf(stderr,"initout can't run command <%s> !!\n",q); MPI_Finalize(); exit(1); }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);


  if (openfile==1) {
    MPI_Barrier(MPI_COMM_WORLD);
    sprintf(q,"%s/stdout%02d.txt",path,quisoc());
    if (info) printf("q=%s\n",q);
    if (tonull!=1) myout=fopen(q,"w");
  } else
    myout=stdout;

  if (tonull==1) myout=fopen("/dev/null","w");

  if (myout==NULL) CRASH("Can't open <%s> to write output\n",q);

  initout=1;
  if (info) printf("initout=%d\n",initout);

  if (flushall==1) pprintf("Warning: fflush is called after each pprintf\n");
#ifdef CHECK
  pprintf("Warning: Matrix range is checked\n");
#else
  pprintf("Warning: Matrix range is NOT checked\n");
#endif

  strcpy(outpath_,path);
}

void end_pprintf() {
  pprintf("End of output\n");
  fflush(NULL);
  if (myout!=stdout) fclose(myout);
  initout=0;
}

// el macro posa line i file automaticament
void crash(int line,char *file,char *fmt,...) {
	va_list ap;
  FILE *errout;

  errout=stderr;

  //Minimalistic CRASH
  //MPI_Finalize(); exit(1);

  fflush(NULL);

  int flag;
  MPI_Initialized(&flag);

  if (quisoc()==0 || flag==0) {
    fprintf(errout,"CRASH: proc=%d FILE=<%s> LINE=%d\n",quisoc(),file,line);
   	va_start(ap,fmt);
    vfprintf(errout,fmt,ap);
   	va_end(ap);
    fprintf(errout,"\n");
   	fflush(NULL);
    MPI_Abort(MPI_COMM_WORLD,-1);
   }

  MPI_Abort(MPI_COMM_WORLD,-1);
  //MPI_Finalize();
	 exit(0);
}


double *dmem(map *M) { // alloc memory for a scaf 2d or 3d
	double *r;
	int info=0;
	if (info) pprintf("dmem S=%d \n",M->S);
	r=(double *)malloc(sizeof(double)*(M->S));
	if (r==NULL) {
		CRASH("oops out of memory S=%d",M->S);
		return(0);
	}
	return(r);
}


/////
int crax(int i,int a,int line, char *file,map *M) {
	crash(line,file,"Wrong index axis=%d v=%d l0=%d l1=%d PROCESSOR=%d\n",a,i,M->l0[1],M->l1[1], quisoc());
	return(0);
}


// h can be 0 or 1 (1: prints ALL the halos)
// 2d or 3d
void print_scaf(double *x, char *label, int h, map *M,char *fmt) {

	int i,j,k;

	if (M->nd==2) {
		pprintf("scalar field 2d %s global=%d:%d x %d:%d printhalo=%d sh=%d\n",label,M->gl0[1],M->gl0[2],M->gl1[1],M->gl1[2],h,M->sh);
		for (j=M->l1[2]+h*M->sh ; j>=M->l0[2]-h*M->sh ; j--) {
			pprintf("%s j=%03d; i=%03d:%03d ::: ",label,j,M->l0[1]-h*M->sh,M->l1[1]+h*M->sh);
			for (i=M->l0[1]-h*M->sh ; i<=M->l1[1]+h*M->sh ; i++ )
				pprintf(fmt,X(i,j,M) );
			pprintf("\n");
		}
	} else {
		pprintf("scalar field 3d <%s> global=%d:%d:%d x %d:%d:%d  printhalo=%d sh=%d\n",
			label,M->gl0[1],M->gl0[2],M->gl0[3],
				  M->gl1[1],M->gl1[2],M->gl1[3],h,M->sh);
		for (k=M->l0[3]-h*M->sh; k <= M->l1[3]+h*M->sh ; k++) {
			pprintf("print_scaf %s k=%d === \n",label,k);
			for (j=M->l1[2]+h*M->sh ; j>=M->l0[2]-h*M->sh ; j--) {
				pprintf("j=%03d; i=%03d:%03d ::: ",j,M->l0[1]-h*M->sh,M->l1[1]+h*M->sh);
				for (i=M->l0[1]-h*M->sh ; i<=M->l1[1]+h*M->sh ; i++ )
					pprintf(fmt,X3(i,j,k,M) );
				pprintf("\n");
			}
		}
	}


}

// basic functions:

void copy_scaf(double *matrixCopy, double *matrixInput, int copy_halo, map *M) {

/*

Copies the matrix matrixInput to matrixCopy, both characterized by M. The contents of the halo are passed to matrixCopy if copy_halo=1.
- double *matrixCopy (output) The copied matrix with the values of matrixInput.
- double *matrixInput (input) The origin of the values that are passed to matrixCopy.
- double *copy_halo (input) If copy_halo=0, only the domain is copied. If copy_halo=1, the values of the halo are also copied.
- map *M (input) The sppde's map.

This function can be used for passing values between time steps, as the function copies the current value of a given variable matrixInput to matrixCopy, which is the value of the previous time step.
*/

// USE sppde's linop_scaf?

    int i, j;

    if(copy_halo==0) {
        forall(i,j,M) {
            ac(matrixCopy,i,j,M) = ac(matrixInput,i,j,M);
        }
    } else {
        forallh(i,j,M) {
            ac(matrixCopy,i,j,M) = ac(matrixInput,i,j,M);
        }
    }
}


void setzero_scaf(double *r,map *M) { // 2d or 3d, including halos !!
 	int i,j,k;

	Forallh(i,j,k,M) {
    	ac3(r,i,j,k,M)=0;
  	}
}

// R=k1.A+k2.Bxw scaf only inner regions, 2d or 3d
void linop_scaf(double *r,double *a,double *b,double k1,double k2,map *M) {
	int i,j,k;
	if (M->nd==2) {
		forall(i,j,M) {
    		ac(r,i,j,M)=k1*ac(a,i,j,M)+k2*ac(b,i,j,M);
  		}
  	}
	else {
        CRASH("KRAT why 2d forall no halo and 3d Forallh with halo "); // REMOVE THIS AFTER FINDING THE ERROR
		Forallh(i,j,k,M) {
			ac3(r,i,j,k,M)=k1*ac3(a,i,j,k,M)+k2*ac3(b,i,j,k,M);
		}
	}
}

// norm of a scalar field r; if l==2 it is the sqrt of sum(ri^2), otherwise the max(abs(ri))
// uses global communications
// 2d or 3d
#define MAXa(a,b) ( fabs((a)>fabs((b))?fabs((a)):fabs((b)) ) )

double norm_scaf(double *r,int l,map *M) {
	int i,j,k;
	double s=0,gs=0;
	if (M->nd==2) {
		forall(i,j,M) {
			if (l==2) s+=ac(r,i,j,M)*ac(r,i,j,M); else s=MAXa( s , ac(r,i,j,M) );
  		}
  	}
	else {
        CRASH("KRAT why 2d forall no halo and 3d Forallh with halo "); // REMOVE THIS AFTER FINDING THE ERROR
		Forallh(i,j,k,M) {
			if (l==2) s+=ac3(r,i,j,k,M)*ac3(r,i,j,k,M); else s=MAXa( s , ac3(r,i,j,k,M) );
		}
	}
	if (l==2) {
	    MPI_Allreduce(&s, &gs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // add all
	    gs=sqrt(gs);
	}
	else {
	    MPI_Allreduce(&s, &gs, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // max
	}

	return(gs);
}


double norm_dif_scaf(double *a,double *b,int l,map *M) {
  double *c;
  double v;
  int info=0;
  int i,j;
  c=dmem(M);
  forall(i,j,M)
    ac(c,i,j,M)=ac(a,i,j,M)-ac(b,i,j,M);
  if (info>1) print_scaf(c,"delta",0,M,"%9.3e ");
  if (info) stats_scaf(c,M,"delta");
  v=norm_scaf(c,l,M);
  free(c);
  return(v);
}

void get_stats_scaf(double *r,map *M,char *fname,double *gmin_,double *gmax_,double *gavg_) { // eval and return (not print)
  int i,j,k;
  double s=0,gs=0;
  double n=0,gn=0;
  double min = 0,gmin=0;
  double max = 0,gmax=0;
  int info=0;

  if (M->nd==2) {
    forall(i,j,M) {
      n=n+1;
      s+=ac(r,i,j,M);
      if (ac(r,i,j,M)<min ||n==1.0) min=ac(r,i,j,M);
      if (ac(r,i,j,M)>max ||n==1.0) max=ac(r,i,j,M);
      }
    }
  else {
    CRASH("OOps only implemented for 2d ... please finish it (pecador)");
  }

  MPI_Allreduce(&s  , &gs  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&n  , &gn  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&min, &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


  if (info) pprintf("get_stats_scaf fname=<%s> min=%e max=%e avg=%e \n",fname,gmin,gmax,gs/gn);

  *gmin_=gmin;
  *gmax_=gmax;
  *gavg_=gs/gn;
}

// molt lenta i mal implementada amb 4 all-to-all !!! nomes per diagnostics !!!
void stats_scaf(double *r,map *M,char *fname) {
  double gmin,gmax,gavg;
  get_stats_scaf(r,M,fname,&gmin,&gmax,&gavg);
  pprintf("stats_scaf fname=<%s> min=%e max=%e avg=%e \n",fname,gmin,gmax,gavg);
}

void check_stats_scaf(double *r,map *M,char *fname,int step,double MIN, double MAX,int prn) {
  double gmin,gmax,gavg;
  get_stats_scaf(r,M,fname,&gmin,&gmax,&gavg);
  if (prn) pprintf("stats_scaf fname=<%s> step=%d min=%e max=%e avg=%e \n",fname,step,gmin,gmax,gavg);
  if (gmin<MIN || gmax>MAX || isfinite(gmin)==0 || isfinite(gmax)==0 ) {
    CRASH("OOOPS nonsense value detected min=%e max=%e isfinite(min)=%d isfinite(max)=%d range acceptable = %e to %e \n",
      gmin,gmax,isfinite(gmin),isfinite(gmax),MIN,MAX);
  }

}

void printmap(char *name,map *M) {
	pprintf("printmap name=%s\n",name);
	pprintf("printmap M->nd  = %d dimensions (2 o 3)\n",M->nd);
	pprintf("printmap M->sh  = %d halo size\n",M->sh);
	pprintf("printmap M->gl0 = %d %d %d domain global limits inicial \n",M->gl0[1],M->gl0[2],M->gl0[3]);
	pprintf("printmap M->gl1 = %d %d %d domain global limits final   \n",M->gl1[1],M->gl1[2],M->gl1[3]);
	pprintf("printmap M->l0  = %d %d %d local limit inicial EXCLUDING halo \n",M->l0[1],M->l0[2],M->l0[3]);
	pprintf("printmap M->l1  = %d %d %d local limit final   EXCLUDING halo\n",M->l1[1],M->l1[2],M->l1[3]);
	pprintf("printmap M->s   = %d %d %d local size in each axis (for internal use) \n",M->s[1],M->s[2],M->s[3]);
	pprintf("printmap M->S   = %d total size (for internal use) \n",M->S);
	pprintf("printmap gpc    = %d x %d x %d global number of procs in each axis\n",M->gpc[1],M->gpc[2],M->gpc[3]);
	pprintf("printmap pc     = %d x %d x %d my local coordinates in the proc grid \n",M->pc[1],M->pc[2],M->pc[3]);
	pprintf("printmap M->nb_ e=%d w=%d n=%d s=%d t=%d b=%d my neighbouring processors (-1: no nb.)\n",
				M->nb_e,M->nb_w,M->nb_n,M->nb_s,M->nb_t,M->nb_b);
	pprintf("printmap bc0    = %d %d %d Do I own a part of the BC ?? initial semiaxis \n",M->bc0[1],M->bc0[2],M->bc0[3]);
	pprintf("printmap bc1    = %d %d %d Do I own a part of the BC ?? final semiaxis \n"  ,M->bc1[1],M->bc1[2],M->bc1[3]);
	pprintf("printmap periodicity in each axis = %d %d %d \n",M->per[1],M->per[2],M->per[3]);
	for (int p=0;p<=M->gpc[1]*M->gpc[2]*(M->nd==3? M->gpc[3]:1)-1; p++ ) {
		pprintf("printmap proc=%d ",p);
		for (int a=1;a<=M->nd;a++) {
			pprintf("al0[%2d]=%3d al1[%2d]=%3d ; ",a,M->al0[p][a],a,M->al1[p][a] );
		}
		pprintf("\n");
	}

  for (int p=0;p<=M->gpc[1]*M->gpc[2]*(M->nd==3? M->gpc[3]:1)-1; p++ ) {
    pprintf("printmap pco[%2d][1:3]=%2d %2d %2d \n",p,M->pco[p][1],M->pco[p][2],M->pco[p][3]);
  }
}


/////////
#define MAXP1 256 // maxim nombre de procs en 1 eix

static void split1d(int ini,int fi,int NP,int *i,int *f,int minsize) {
	int n[MAXP1+1];
	int l=fi-ini+1;
	int p;
	int info=0;
	if (NP>MAXP1) CRASH("split1d NP=%d and should be less than %d",NP,MAXP1);
	for (p=1;p<=NP;p++) {
		n[p]=l/NP;
    if (n[p]<minsize) CRASH("One of the partitions is too small, reduce number of processors");
	}
	int d=l-(l/NP)*NP;
	for (p=1;p<=NP;p++) {
		if (d>0) {
			n[p]++;
			d--;
		}
	}
	i[1]=ini;
	for (p=1;p<=NP;p++) {
		if (p>1) i[p]=f[p-1]+1;
		f[p]=i[p]+n[p]-1;
	}
	if (info) {
		pprintf("ini=%d fi=%d NP=%d \n",ini,fi,NP);
		for (p=1;p<=NP;p++) {
			pprintf("p=%d i=%d f=%d n[p]=%d \n",p,i[p],f[p],n[p]);
		}
	}
}

void createM(int nd,int *gl0,int *gl1,int *np,int sh,int *per,map *M) {

 	int npx,npy,npz,i1,f1,i2,f2,i3,f3;

  	npx=np[1]; npy=np[2]; npz=np[3];
	if (npx*npy*npz > MAXP) CRASH("Fariseus pecadors !!! Too many processors !! P=%d and should be less than %d ",npx*npy*npz,MAXP);
  	i1=gl0[1]; f1=gl1[1];
  	i2=gl0[2]; f2=gl1[2];
  	i3=gl0[3]; f3=gl1[3];

  	if (nd<2 || nd>3) CRASH("nd=%d has to be 2 or 3",nd);

  	M->nd=nd;

	if (M->nd == 2 && npx*npy != quants() )
		CRASH("Inconsistent 2d division npx=%d npy=%d np=%d ",npx,npy,quants());

	if (M->nd == 3 && npx*npy*npz != quants() )
		CRASH("Inconsistent 3d division npx=%d npy=%d npz=%d np=%d ",npx,npy,npz,quants());

	if (nd==3 && npz!=1)
		CRASH("For 3d maps, only decomposition in x and y is supported so far");

	int iv1[MAXP1+1],fv1[MAXP1+1];
	int iv2[MAXP1+1],fv2[MAXP1+1];

	split1d(i1,f1,npx,iv1,fv1,3);
	split1d(i2,f2,npy,iv2,fv2,3);

  	M->sh=sh;
  	M->gl0[0]=0;
  	M->gl0[1]=i1;
  	M->gl0[2]=i2;
  	M->gl0[3]=(nd==3)? i3: 1;
  	M->gl1[0]=0;
  	M->gl1[1]=f1;
  	M->gl1[2]=f2;
  	M->gl1[3]=(nd==3)? f3: 1;

	int mypx,mypy,px,py,p;
  	int done=0,pp=0;

	// store domain decomposition vectors
	p=0;
	for (py=1;py<=npy;py++)
	for (px=1;px<=npx;px++) {
		M->al0[p][1]=iv1[px]; M->al1[p][1]=fv1[px];
		M->al0[p][2]=iv2[py]; M->al1[p][2]=fv2[py];
		M->al0[p][3]=i3;       M->al1[p][3]=f3;	 // ============= careful with 3d dec =================
		p++;
	}

	// grid of processors
	p=0;
	for (py=1;py<=npy;py++)
	for (px=1;px<=npx;px++) {
		if (done== 0 && p==quisoc()) { mypx=px; mypy=py; done=1; pp=p; }
    M->pco[p][1]=px; // ===================================<<<< careful when z dec is available !!
    M->pco[p][2]=py;
    M->pco[p][3]=1;
		p++;
	}
  	M->gpc[0]=0; M->gpc[1]=npx;  M->gpc[2]=npy;
  	M->gpc[3]=1; // ===================================<<<< careful when z dec is available !!
  	M->pc[0]=0;  M->pc[1]=mypx; M->pc[2]=mypy; M->pc[3]=1;

  	M->l0[0]=0;
  	M->l0[1]=iv1[mypx]; // ini x
  	M->l0[2]=iv2[mypy]; // ini y
  	M->l0[3]=i3; // ===================================<<<< careful when z dec is available !!

  	M->l1[0]=0;
  	M->l1[1]=fv1[mypx];  // fi x
  	M->l1[2]=fv2[mypy]; // fi y
  	M->l1[3]=f3;  // ===================================<<<< careful when z dec is available !!

  	M->s[0]=0;
  	M->s[1]=M->l1[1]-M->l0[1]+2*M->sh+1;
  	M->s[2]=M->l1[2]-M->l0[2]+2*M->sh+1;
  	M->s[3]=M->l1[3]-M->l0[3]+2*M->sh+1;

  	M->S=(M->s[1])*(M->s[2])*(M->s[3]);

  	// check ..
  	for (int a=1; a<=M->nd;a++) {
	  	if (M->al0[quisoc()][a] != M->l0[a])
	  		CRASH("global map/0 %d != %d uhh??",M->al0[quisoc()][a],M->l0[a]);
	  	if (M->al1[quisoc()][a] != M->l1[a])
	  		CRASH("global map/1 %d != %d uhh??",M->al1[quisoc()][a],M->l1[a]);
  	}

  	// NEIGHBOURING PROCESSORS
  	if (mypx<npx) M->nb_e=pp+1;   else M->nb_e=-1;
  	if (mypx>1  ) M->nb_w=pp-1;   else M->nb_w=-1;
  	if (mypy<npy) M->nb_n=pp+npx; else M->nb_n=-1;
  	if (mypy>1  ) M->nb_s=pp-npx; else M->nb_s=-1;
  	M->nb_t=-1; M->nb_b=-1;

  	// include periodic domains
  	if (per[3]==1) CRASH("periodic in z not implemented yet, sorry !!");

  	// if (per[1]==1 && npx%2 ==1 && mypx>1) CRASH("pecador!!! posa un nombre parell de procs als eixos periodics (x)"); // QRAT
  	if (per[1]==1) { // periodic in x
    	if (mypx==npx) // if I'm the last
        M->nb_e=pp-npx+1;
        if (mypx==1)
     	M->nb_w=pp+npx-1;
  	}

    // if (per[2]==1 && npy%2 ==1 && mypy>1) CRASH("pecador!!! posa un nombre parell de procs als eixos periodics (y)"); // QRAT
    if (per[2]==1) { // periodic in y
        if (mypy==npy) // if I'm the last
        M->nb_n=pp-(npy-1)*npx;
        if (mypy==1)
        M->nb_s=pp+(npy-1)*npx;
    }


  	// Flags to signal that I own part of each boundary
  	M->bc1[3] = 1; // TOP ========================== careful with zdec !!!!!!
  	M->bc1[2] = (M->pc[2] == M->gpc[2] );  // NORTH BC: 1==I own at least part of top boundary condition
  	M->bc1[1] = (M->pc[1] == M->gpc[1] );  // EAST
  	M->bc0[3] = 1; // TOP ========================== careful with zdec !!!!!!
  	M->bc0[2] = (M->pc[2] == 1 );          // SOUTH
  	M->bc0[1] = (M->pc[1] == 1 );          // WEST

  	// in periodic domains there are no boundary conditions
  	if (per[1]==1) { M->bc1[1]=0; M->bc0[1]=0; }
    if (per[2]==1) { M->bc1[2]=0; M->bc0[2]=0; }

  	M->per[1]=per[1];
  	M->per[2]=per[2];
  	M->per[3]=per[3];
}


void createM_2d(int npx,int npy,int gsx,int gex,int gsy,int gey,int sh,int perx,int pery,map *M) {
    int gl0[4]={0,gsx,gsy,0}; // global limits (start)
    int gl1[4]={0,gex,gey,0}; // global limits (end)
    int np[4] ={0,npx,npy,1}; // processor grid
    int per[4]={0,perx,pery,0}; //  periodicity
    createM(2,gl0,gl1,np,sh,per,M);
}


// from map M, creates a map MJ to contain all the data shared among the different processors
void joinM(map *M, map *MJ) {
    int a;
	if (M->nd==3) CRASH("joinM 3d not implemented yet sorry.. ");

	MJ->nd=M->nd;
    MJ->sh=M->sh; // MANEL preserve halo size
    // global size
    for (a=1;a<=3;a++) MJ->gl0[a]=M->gl0[a];
    for (a=1;a<=3;a++) MJ->gl1[a]=M->gl1[a];
    MJ->gl0[3]=1; MJ->gl1[3]=1; // MANEL 2d
    // MJ local size is equal to global size
    for (a=1;a<=3;a++) MJ->l0[a]=MJ->gl0[a];
    for (a=1;a<=3;a++) MJ->l1[a]=MJ->gl1[a];
    // s vector
    for (a=1;a<=3;a++) MJ->s[a]=MJ->l1[a]-MJ->l0[a]+2*MJ->sh+1;
    MJ->s[3]=1; // MANEL 2d
    MJ->S=(MJ->s[1])*(MJ->s[2])*(MJ->s[3]);
    // number of procs
    for (a=1;a<=3;a++) MJ->gpc[a]=1;
    for (a=1;a<=3;a++) MJ->pc[a]=1;
    // No neighbours
    MJ->nb_e=-1;
    MJ->nb_w=-1;
    MJ->nb_n=-1;
    MJ->nb_s=-1;
    MJ->nb_t=-1;
    MJ->nb_b=-1;
    for (a=1;a<=3;a++) MJ->bc0[a]=1;
    for (a=1;a<=3;a++) MJ->bc1[a]=1;
    // distribution ...
   	for (a=1;a<=3;a++) {
		MJ->al0[0][a]=MJ->gl0[a]; MJ->al1[0][a]=MJ->gl1[a];
   	}
    // processor coordinates
    MJ->pco[0][1]=1; MJ->pco[0][2]=1;  MJ->pco[0][3]=1;
    // inherit periodicity
    for (a=1;a<=3;a++) MJ->per[a]=M->per[a];
    if (MJ->per[1]) { MJ->bc0[1]=0;  MJ->bc1[1]=0; }
    if (MJ->per[2]) { MJ->bc0[2]=0;  MJ->bc1[2]=0; }
    if (MJ->per[3]) CRASH("uuups can't be periodic in z");
}
