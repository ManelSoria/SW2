// sppde: structured parallel pde solver
// Manel Soria 2017 - UPC - ESEIAAT - TUAREG

// internal functions

int pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,map *M); // 2d
int un_pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,map *M);

int g_pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,int k1,int k2,map *M);
int g_un_pack(double *b,int q0,double *x,int i1,int i2,int j1,int j2,int k1,int k2,map *M);


void easy_s(int nb, double *b, int ndata); // Defined in sppde.h
void easy_r(int nb,double *b,int ndata); // Defined in sppde.h
void easy_sr(int nb,double *bs,double *br,int ndata);
void halo_update_y(double *x,map *M);
void halo_update_x(double *x,map *M);
void halo_update_single_map_y(double *x,map *M);
void halo_update_single_map_x(double *x,map *M);



// test functions ========================================= TO BE REMOVED
void testmemo1(map *M);
void test_distributed_arrays();
void test_halo_update_3d(map *M);

void test_distributed_arrays_2d(map *M);
void benchmark_halo_update_2d(int NX,int NY,int NPX,int NPY);
