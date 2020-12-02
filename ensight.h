// SW2
// 2018-2020
// Manel Soria, Arnau Prat, Arnau Sabates, Marc Andres, Enrique Garcia-Melendo
// UPC - ESEIAAT - TUAREG

#define N_OUT_VARIABLES 7 // Number of variables of which output files are wanted

void save_ens_case(char *bname, char **var_names, int nstart, int nsteps, double *timevals);
void save_ens_geo_sw(char *bname, int mode, sw *SW);
void save_ens_geo(char *bname, double **coords, int *nelems, map *M);
void save_ens_step(char *bname, char **var_names, double **var_fields, double (*var_gridpos)[2], int frame, sw *SW);
void save_ens_var(char *bname, char *scafname, double *scaf, int step, int cx, int cy, map *M);
