/*
This file is part of SW2 code

2018-2020
Manel Soria, Arnau Prat, Arnau Sabates, Marc Andres-Carcasona, Arnau Miro, Enrique Garcia-Melendo
UPC - ESEIAAT - TUAREG

(c) Manel Soria, Enrique Garcia-Melendo 2018-2020

LICENSED UNDER: Attribution 4.0 International

*/


/* parser for text input files
Example of file:

#: comment
define(R, 33.34) # macro definition
define(o, 0.12) # macro definition
4 <R*cos(o)> <R*sin(o)> # expressions within <> are evaluated
include(tutu) # include a file
# blank lines are ignored

  this_is_a_tag 4.0 # tags are the first token after spaces or tabs (if any)
                    # functions can check the tag (if tag string is non-empty)
#functions that use tags are:
one_double <2+2>
ont_integer -34
multiple_doubles_in_a_line 4.5 <2+2> 8.9
string filename
*/


// opens file for reading. At the end, it should be closed with p_endinput() and NOT fclose
// if copyinput==1, it copies the file to the working folder

FILE *p_startinput(char *fname, int copyinput);
void p_endinput(FILE *f);

void p_setdebug(int dval); // sets debug level, 1 or 2
void p_stopdebug(); // stops debugging (level set to 0)


char *p_nextline(FILE *fin); // gets a pointer to next text line. If is to be reused, data has to be copied to another string

int p_getint(FILE *f, char *tag); // checks tag (if not empty) and gets an integer (rest of the line is ignored)

double p_getdouble(FILE *f, char *tag); // checks tag (if not empty) and gets a double (rest of the line is ignored)

void p_getndoubles(FILE *f, int nd, double *d, char *tag); // checks tag (if not empty) and  gets n doubles and stores in a preallocated vector d

char* p_getstring(FILE *f,char *tag); // checks tag (if not empty) and gets a string

double *p_alloc_gettable(FILE *f,int *nr,int *nc,int transpose);
// allocates and reads a table from f, with the following format
// r rows
// c columns
// 1rst row..
// ...
// r row
// (the table may contain expressions to be evaluated)
// if transpose==0, the table is read line by line
// ie:
// 1 2 3
// 4 5 6
// is stored in memory as 1 2 3 4 5 6
// else, it is stored as 1 4 2 5 3 6
