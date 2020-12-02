// sppde: structured parallel pde solver
// Manel Soria 2017 - UPC - ESEIAAT - TUAREG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h> // needed for crash
#include <mpi.h>

#include "sppde.h"
#include "sppde_tinyexpr.h"
#include "sppde_parser.h"



// **********************************************************************

#define MAXL 1000 // buffer
char qq[MAXL];
int debug; // if 1, prints all input after preprocessing


int isthischar(char c,char *sep) {
//	pprintf("isthischar c=%c sep='%s'\n",c,sep);
	while (*sep!='\0') {
		if (*sep==c) return(1);
		sep++;
	}
	return(0);
}

char *get_first_token(char *l,char *first) { // writes on first, returns pointer to rest
	char sep[]=" \t";
	int q;
//	pprintf("getonlyfirsttoken0='%s' \n",l);
	while (isthischar(*l,sep)==1) { l++; } // skip to first non-sep
//	pprintf("getonlyfirsttoken1='%s' \n",l);
	q=0;
	do { // copy first token
//		pprintf("loopget *l=%c q=%d \n",*l,q);
		first[q]=*l;
		l++;
		q++;
	} while (isthischar(*l,sep)==0);
	first[q]='\0';
//	pprintf("getonly btok1='%s' l='%s'\n",first,l);
	while (isthischar(*l,sep)==1) { l++; } // skip to first non-sep
//	pprintf("get_first_token btok1='%s' l='%s'\n",first,l);
	return(l);
}



int cleanstr(char *qq) { // removes comments; returns 1 if the line is non-empty
	int text=0;
	int info=0;
	if (info) pprintf("cleanstr entering qq=<%s>\n",qq);
	while (*qq!='\0') {
		if (info) pprintf("cleanstr *qq=%c\n",*qq);
		if (*qq=='\n' || *qq=='#') { *qq='\0'; break;}
		if (*qq!=' ' && *qq!='\t') text=1;
		qq++;
	}
	return(text);
}

// if one or more "< expression>" are found in qq, calls the expression evaluator and replaces expression by its value
// e.g. < 2+3^2 > will be replaced by 11
// stops if there are errors, after printing message
// '<' and '>' that are not to be evaluated should be preceded by '\'
// ie: "\<this is not \evaluated"  will produce: "<this is not \evaluated"
// returns the number of expressions evaluated

int evalstr(char *input) {
	char *qq=input;
	char *start=NULL;
	char result[MAXL];
	char numval[MAXL];
	char expr[MAXL];
	int copy=0,p,pr=0;
	int nch;
	int error;
	int nexpr=0;
	int backs=0;
	int skip=0;
	while (*qq!='\0') {

		if (backs==1) {
			if (*qq=='<' || *qq=='>') {skip=1; backs=0; } // skip expression, dont produce '\' symbol
		}

		if (*qq=='\\') { // oops.. wait to next char
			backs=1;
			qq++;
			continue;
		}

		if (*qq=='>' && skip==0 ) { // end of expression detected
			copy=0;
			expr[p]='\0';   // end expr string
			double val=te_interp(expr,&error);
			if (error) CRASH("Error. Expression='%s' can not be evaluated",expr);
			nch=sprintf(numval,"%.15e",val);
			for (int i=0;i<=nch-1;i++) { // insert result in output string
				result[pr]=numval[i];
				pr++;
			}
			qq++;
			nexpr++;
			continue;
		}
		if (*qq=='<' && skip==0 ) { // beginning of expression detected
			copy=1; p=0;
			qq++;
			continue;
		}
		if (copy==1) { // copy to expression
			expr[p]=*qq;
			p++;
		} else { // copy to output string
			if (backs==1) { result[pr]='\\'; pr++; backs=0;} // recover the missing backslash
			result[pr]=*qq;
			pr++;
		}
		qq++;
	}
	if (copy==1) CRASH("oops missing '>' while evaluating '%s'",input);
	result[pr]='\0';
	for (int i=0;i<=strlen(result);i++) // include '\0'
		input[i]=result[i];
	return(nexpr);
	//pprintf("bye: result='%s'\n",result);
}

char *p_nextline(FILE *fin) {
	int info=0;
	if (info) pprintf("p_nexline quisoc=%d\n",quisoc());
	while (fgets(qq,MAXL,fin) != NULL) {
		if (cleanstr(qq)) {
			if (debug>1) pprintf("INPUT='%s'\n",qq);
			int ne=evalstr(qq);
			if (debug>1 && ne) pprintf("EVALUATED='%s'\n",qq);
			return(qq);
		}
	}
	CRASH("Unexpected end of file");
	return((char *)NULL);
}


char *getstrchecktag(FILE *f,char *tag) {
	char *l;
	char first[MAXL];
	char *rest;
	l=p_nextline(f);
//	pprintf("get_str_check_tag l='%s' tag='%s'\n",l,tag);
	if (strlen(tag)>0) {
		rest=get_first_token(l,first);
//		pprintf("get_str_check_tag first='%s' rest='%s'\n",first,rest);
		if (strcmp(first,tag)!=0) CRASH("expecting '%s' got '%s'",tag,first);
		return(rest);
	} else
		return(l);
}


FILE *p_startinput(char *fname, int copyinput) {
	FILE *fin;
	int info=0;
	int r=0,rcopy=0,rg;

	char command[MAXL];
	if (strlen(fname)>MAXL-20) CRASH("oops input file name=<%s> is too long",fname);

	if (quisoc()==0) { // input file is pre-processed with m4, the result is stored in standard output folder
		if (copyinput==1) {
			sprintf(command,"cp %s %s",fname,poutp()); // copy input file to output folder, for later use
			pprintf("sppde: command=<%s>\n",command);
			rcopy=system(command);
			if (info) pprintf("rcopy=%d\n",rcopy);
		}
		sprintf(command,"m4 %s > %s/%s.prem4",fname,poutp(),fname);
		pprintf("sppde: command=<%s>\n",command);
		r=system(command);
		if (info) pprintf("m4 returns r=%d\n",r);
		r=r+rcopy;
	}


	checkr(MPI_Allreduce(&r, &rg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD),"Allreduce_startinput");

	if (rg!=0) CRASH("Error could not process '%s'\nCheck data file '%s' and if present, then preprocessor commands such as define and include",command,fname);

	sprintf(command,"%s/%s.prem4",poutp(),fname);
	pprintf("sppde: pre-processed input file = <%s>\n",fname);
	fin=fopen(command,"r");
	if (fin==NULL) CRASH("oops could not open '%s'",command);

	return(fin);
}

void p_setdebug(int dval) {
	debug=dval;
}

void p_stopdebug() {
	debug=0;
}

void p_endinput(FILE *f) {
	pclose(f);
}

char* p_getstring(FILE *f,char *tag) {
	char *l;
	l=getstrchecktag(f,tag);
	if (debug) {
		pprintf("read string ");
		if (strlen(tag)) pprintf("%s=",tag);
		pprintf("'%s'",l);
		pprintf("\n");
	}
	return(l);
}

int p_getint(FILE *f, char *tag) {
	int nd,nextint;
	double aux,iptr;
	char *nl=getstrchecktag(f,tag);
	nd=sscanf(nl,"%le",&aux);
	if (nd!=1) CRASH("error, expecting valid integer and found '%s' ",nl);
	if (modf(aux,&iptr) !=0) CRASH("error, expecting an integer and found '%s'",nl);
	if (debug) {
		pprintf("read int ");
		if (strlen(tag)) pprintf("%s=",tag);
		pprintf("%d",(int)aux);
		pprintf("\n");
	}
	return((int)aux);
}

double p_getdouble(FILE *f,char *tag) {
	int nd;
	double nextdouble;
	char *nl=getstrchecktag(f,tag);
	nd=sscanf(nl,"%le",&nextdouble);
	if (nd!=1) CRASH("error, reading tag '%s' expected valid double and found '%s'",tag,nl);
	if (debug) {
		pprintf("read double ");
		if (strlen(tag)) pprintf("%s=",tag);
		pprintf("%e",nextdouble);
		pprintf("\n");
	}
	return(nextdouble);
}

void p_getndoubles(FILE *f, int nd, double *d,char *tag) {

	char *l=getstrchecktag(f,tag);

	char *token;
	int ok;
	for (int i=1;i<=nd;i++) {
		if (i==1)
			token=strtok(l," \t");
		else
			token=strtok(NULL," \t");
		if (token==NULL) CRASH("error, expecting more doubles nd=%d i=%d",nd,i);
		// pprintf("----> l=<%s>\n",token);
		ok=sscanf(token,"%lf",d+i-1);
		if (ok!=1)  CRASH("error, expecting one double and found=<%s> ",token);
	}
	if (debug) {
		pprintf("read n doubles ");
		if (strlen(tag)) pprintf("%s=",tag);
		for (int i=1;i<=nd;i++) {
			pprintf("%e ",*(d+i-1));
		}
		pprintf("\n");
	}
}

double *p_alloc_gettable(FILE *f,int *nr,int *nc, int transpose) {
	double *table;
	double onerow[100];
	if (debug) pprintf("Reading table\n");
	*nr=p_getint(f,"");
	if (debug) pprintf("number of rows=%d\n",*nr);
	*nc=p_getint(f,"");
	if (debug) pprintf("number of columns=%d\n",*nc);
	if (*nc>100) { pprintf("uhhh? columns=%d ... too many \n",*nc); exit(0); }

	table=(double *)malloc(sizeof(double)*(*nr)*(*nc));
	for (int r=1;r<=*nr;r++) {
		p_getndoubles(f,*nc,onerow,"");
		for (int c=1;c<=*nc;c++) {
			if (!transpose) { // PRAT Added by Arnau Prat Gasull
				table[ (r-1)*(*nc) + c-1 ] = onerow[c-1]; // Rows contiguous in memory
			} else {
				table[ r-1 + (c-1)*(*nr) ] = onerow[c-1]; // Columns contiguous in memory
			}
		}
	}

	if (debug) pprintf("Ended reading table.\n");
	return(table);
}


/*
int main(int argc,char **argv) {

	checkr(MPI_Init(&argc,&argv),"init"); // MPI should be started before calling p_startinput

	FILE *fin;
	fin=p_startinput("inputexample.sw");
	p_setdebug(1); // turn on debug, extra level
	double radius=p_getdouble(fin,"planet_radius");
	double startl=p_getdouble(fin,"start_latitude");
	double endl=p_getdouble(fin,"end_latitude");
	double pars[2];
	p_getndoubles(fin,2,pars,"parameters");
	char *outputfile=p_getstring(fin,"output"); // if string has to be used, it must be copied to another location

	double *table;
	int nr,nc;
	table=p_alloc_gettable(fin,&nr,&nc,0);
	for (int i=0;i<=nr*nc-1;i++)
		printf("i=%d table=%e \n",i,table[i]);
	free(table);
	pclose(fin); // end input


	MPI_Finalize();
}

*/
