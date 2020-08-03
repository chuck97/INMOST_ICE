#ifndef   __ANI_C2F_H
#define   __ANI_C2F_H


#ifdef   __cplusplus
extern "C" {
#endif


typedef int      f77i;
typedef double   f77r;
typedef char     f77c;
typedef int      f77b;

static const int ANI2D_CONTROL_MaxSkipE = 0;
static const int ANI2D_CONTROL_MaxQItr  = 1;
static const int ANI2D_CONTROL_status   = 2;
static const int ANI2D_CONTROL_flagAuto = 3;
static const int ANI2D_CONTROL_iPrint   = 4;
static const int ANI2D_CONTROL_iErrMesg = 5;
 

typedef struct{
/* variables for main calls */
  f77i      nv, nvmax, nvfix;
  f77r      *vrt;
  f77i      *labelV, *fixedV;
  
  f77i      nb, nc, nbmax, nbfix;
  f77i      *bnd;
  f77r      *crv;
  f77i      *labelB, *fixedB, *labelC;

  f77i      nt, ntmax, ntfix;
  f77i      *tri, *labelT, *fixedT;

  f77r      Quality;
  f77i      *control, nEStar;

  f77i      MaxWr, MaxWi;
} ani2D;


/* simple useful routines */
#if defined(__WINDOWS__)
__inline int C2Fmax(int a, int b);
__inline int C2Fmin(int a, int b);
#else
inline int C2Fmax(int a, int b);
inline int C2Fmin(int a, int b);
#endif


/* C INTERFACES TO FORTRAN ROUTINES */
void loadmani_(
     f77i*, f77i*, f77i*, f77r*, f77i*, f77i*, 
     f77i*, f77i*, f77i*, f77i*, f77i*, f77i*, 
     f77i*,               f77r*, f77i*,
     f77i*, f77i*, f77i*, f77i*, f77i*, f77i*, 
     const f77c*);

void loadmaniheader_(
     f77i*, f77i*, f77i*, f77i*, 
     f77i*, f77i*, f77i*, const f77c*);

void savemani_(
     f77i*, f77i*, f77r*, f77i*, f77i*, 
     f77i*, f77i*, f77i*, f77i*, f77i*, 
     f77i*,        f77r*, f77i*,
     f77i*, f77i*, f77i*, f77i*, f77i*, 
     const f77c*);

void savemgmv_(
      f77i*, f77i*, f77r*, f77i*, const f77c*);

void mbaanalytic_(
/* M */ f77i*, f77i*, f77i*, f77r*, f77i*, f77i*, 
        f77i*, f77i*, f77i*, f77i*, f77i*, f77i*,
        f77i*,               f77r*, f77i*, void*, 
        f77i*, f77i*, f77i*, f77i*, f77i*, f77i*,
/* CONTROL */  
        f77i*, f77r*, f77i*, void*, 
/* W */ f77i*, f77i*, f77r*, f77i*, f77i*);

void mbanodal_(
/* M */ f77i*, f77i*, f77i*, f77r*, f77i*, f77i*,
        f77i*, f77i*, f77i*, f77i*, f77i*, f77i*,
        f77i*,               f77r*, f77i*, void*,
        f77i*, f77i*, f77i*, f77i*, f77i*, f77i*,
/* CONTROL */
        f77i*, f77r*, f77i*, f77r*,
/* W */ f77i*, f77i*, f77r*, f77i*, f77i*);


void draw_(f77i*, f77i*, f77i*, f77r*, f77i*, f77i*, 
           f77i*, f77i*, const f77c*);


/* BASIC ROUTINES */
int ani2D_INIT(ani2D* ani, f77i nEStar);
int ani2D_KILL(ani2D* ani);

int ani2D_load_mesh(ani2D* ani, const char* fname);
int ani2D_save_mesh(ani2D* ani, const char* fname);
int ani2D_save_mesh_gmv(ani2D* ani, const char* fname);

int ani2D_analytic(ani2D* ani, void* metric, void* crvfunc);
int ani2D_nodal(ani2D* ani, double* metric, void* crvfunc);

int ani2D_draw_mesh(ani2D* ani, const char* fname);


/* ELEMENTAL ROUTINES */
int ani2D_number_points(ani2D* ani);
int ani2D_number_edges(ani2D* ani);
int ani2D_number_curved_edges(ani2D* ani);
int ani2D_number_elements(ani2D* ani);

void ani2D_set_max_points(ani2D* ani, int nvmax);
void ani2D_set_max_edges(ani2D* ani, int nbmax);
void ani2D_set_max_elements(ani2D* ani, int ntmax);

void ani2D_get_point(ani2D* ani, int i, double* xy);
void ani2D_set_point(ani2D* ani, int i, double* xy);
void ani2D_fix_point(ani2D* ani, int i);

void ani2D_get_edge(ani2D* ani, int i, int* edge, int* icrv, int* label);
void ani2D_set_edge(ani2D* ani, int i, int* edge, int  icrv, int  label);
void ani2D_fix_edge(ani2D* ani, int i);

void ani2D_get_crv(ani2D* ani, int i, double* par, int* iFnc);
void ani2D_set_crv(ani2D* ani, int i, double* par, int  iFnc);

void ani2D_get_element(ani2D* ani, int i, int* tri, int* label);
void ani2D_set_element(ani2D* ani, int i, int* tri, int  label);
void ani2D_fix_element(ani2D* ani, int i);

void ani2D_get_quality(ani2D* ani, double* rQ);
void ani2D_set_quality(ani2D* ani, double   Q);

void ani2D_get_status(ani2D* ani, int* status);
void ani2D_set_status(ani2D* ani, int  status);

void ani2D_get_max_iters(ani2D* ani, int* max_iters);
void ani2D_set_max_iters(ani2D* ani, int  max_iters);

void ani2D_get_max_basket(ani2D* ani, int* max_basket);
void ani2D_set_max_basket(ani2D* ani, int  max_basket);

void ani2D_get_print_level(ani2D* ani, int* iprint);
void ani2D_set_print_level(ani2D* ani, int  iprint);


/* TRIANGLE ROUTINES */
double triangle_area(ani2D* ani, int i);
void randr_(f77r*, f77r*, f77r*, f77r*, f77r*);


#ifdef   __cplusplus
           }
#endif

#endif

