// sppde: structured parallel pde solver
// Manel Soria 2017 - UPC - ESEIAAT - TUAREG

#include <stdio.h>
#include <sys/time.h>
#include <sys/times.h>
#include <string.h>
#include "mpi.h"
#include "sppde.h"

/* ================================================================================= */
/* Modul de control del temps ====================================================== */
/* ================================================================================= */


/* no son per a utilitzar fora d'aqui : */
#define NCH 80
#define MAXLCHNAME 14

static char chname[NCH][MAXLCHNAME+1]; /* noms dels canals */
static double tmax[NCH]; /* temps maxims de cada canal */
static double tmin[NCH]; /* minims */
static double tsum[NCH]; /* totals */
static int    nop[NCH];  /* nombre d'operacions */
static double tini[NCH]; /* instant d'inici; si es == 0 vol dir que el canal no esta cronometrant  */
static int nch=0; /* nombre de canals ocupats */
static int cr_mode=0;

void cr_reset(void) { /* elimina tots els canals i torna a comen�ar */
  int i;
  for (i=0;i<=NCH-1;i++) {
     strcpy(chname[i],"");
     tmax[i]=0;
     tmin[i]=0;
     tsum[i]=0;
     nop[i]=0;
     tini[i]=0;
  }
  nch=0;
  cr_mode=0;
}

/* crea un canal nou, torna el numero */
static int newch(char *n) {
  if (nch>=NCH) CRASH("cr_newch: obrint <%s>. Ja hi ha %d canals oberts i aixo es el maxim",n,NCH);
  if (strlen(n)>MAXLCHNAME) CRASH("cr_newch: string massa llarg <%s>",n);
  strcpy(chname[nch],n);
  nch++;
  return(nch-1);
}

/* busquem el canal (de forma ineficient pero facileta) */
static int findch(char *n) {
  int i,r=-1;
  for (i=0;i<=nch-1;i++) {
    if (strcmp(chname[i],n)==0) {
      r=i;
      break;
    }
  }
  return(r);
}

/* busquem el canal i si no el trobem peta */
static int findch_crash(char *n) {
  int ch;
  ch=findch(n);
  if (ch==-1) {
    pprintf("findch_crash petant. Estat actual:\n");
    cr_info();
    CRASH("cr: el canal <%s> no existeix",n);
  }
  return(ch);
}

/* busquem el canal i si no el trobem el creem */
static int findch_create(char *n) {
  int ch;
  ch=findch(n);
  if (ch==-1) ch=newch(n);
  return(ch);
}

/* torna el nombre de segons que ha passat des d'un instant arbitrari pero fixat
   el valor retornat sera sempre > 0 */
static double gettime() {
 /*
 struct timeval	now;
 struct tms buf;
 double sec=0,fakt=(1./CLOCKS_PER_SEC);
 switch (cr_mode) {
   case 0:
     if (gettimeofday(&now, (struct timezone *) 0)) { printf("cr:gettimeofday no funciona\n"); }
     sec = (double) now.tv_sec + (double) now.tv_usec / 1000000;
     break;
   case 1:
     times (&buf);
     sec = fakt*( buf.tms_utime + buf.tms_stime + buf.tms_cutime + buf.tms_cstime);
     break;
   case 2:
     times (&buf);
     sec = fakt*( buf.tms_utime + buf.tms_cutime);
     break;
   default: CRASH("gettime: invalid cr_mode=%d",cr_mode);
   }
 return(sec+1.0);
 */
 return(MPI_Wtime());
}



static void addsuf(char *in,int suff,char *out) {
  if (strlen(in)<1 || strlen(in)+2>MAXLCHNAME) CRASH("cr:El nom de canal <%s> es massa llarg o massa curt",in);
  if (suff<0)
    sprintf(out,"%s",in);
  else
    sprintf(out,"%s%02d",in,suff);
}


/* imprimir informacio - ordenem de major a menor per la suma*/
void cr_info(void) {
  int i,j,pj;
  int chp[NCH];
  int sort=1;
  double max;

  pprintf("cr_info\tnom\to\tN\tmode\ttmin\ttmax\ttsum\ttavg(s)\n");

  for (j=0;j<=nch-1;j++)
    chp[j]=0;

  for (i=0;i<=nch-1;i++) {
    if (sort) {
      max=-1.0;
      pj=0;
      for (j=0;j<=nch-1;j++) {
	if (chp[j]) continue;
	if (tsum[j]>max) { max=tsum[j]; pj=j; }
      }
      chp[pj]=1;
    } else {
      pj=i;
    }
    pprintf("cr_info:\t%13s\t%2d\t%5d\t%2d\t%4d\t%e\t%e\t%e\t%e",
	    chname[pj],
	    pj,
	    (tini[pj]!=0),
	    nop[pj],
	    cr_mode,
	    tmin[pj],
	    tmax[pj],
	    tsum[pj] ,
	    tsum[pj]/(1.* nop[pj]) );
    pprintf("\n");
  }
}

/* iniciar crono d'un canal */
void cr_start(char *n,int suff) {
  int ch;
  char tmp[MAXLCHNAME+1];
  addsuf(n,suff,tmp);
  ch=findch_create(tmp);
  if (tini[ch]!=0) CRASH("cr:El canal <%s> ja estava cronometrant (2 cr_start seguits)",n);
  tini[ch]=gettime(); /* lo ultim es engegar el crono */
}

/* parar crono d'un canal */
void cr_end(char *n,int suff) {
  int ch;
  char tmp[MAXLCHNAME+1];
  double end,time;
  end=gettime(); /* abans que res parem el crono ! */
  addsuf(n,suff,tmp);
  ch=findch_crash(tmp);
  time=end-tini[ch];

  nop[ch]=nop[ch]+1;
  if (time > tmax[ch] || nop[ch]==1 ) tmax[ch]=time;
  if (time < tmin[ch] || nop[ch]==1 ) tmin[ch]=time;
  tsum[ch] += time;

  tini[ch]=0;
}

/* et dona el temps d'un crono que estigui contant; el crono segueix contant  */
double cr_time(char *n,int suff) {
  int ch;
  char tmp[MAXLCHNAME+1];
  double end,time;
  end=gettime(); /* abans que res parem el crono ! */
  addsuf(n,suff,tmp);
  ch=findch_crash(tmp);
  time=end-tini[ch];
  return(time);
}

/* modes de comptar el temps d'execucio */
/*  0 -> Temps fisic.
    1 -> Temps d'usuari + sistema
    2 -> Temps d'usuari                 */
void cr_setmode (int mode) {
  if (nch) CRASH("cr_setmode: There's no sense in changing the cr_mode with opened channels");
  if ( !(mode >=0 && mode <=2) ) CRASH("cr_setmode: mode = %d does not exist.",mode);
  cr_mode=mode;
}
