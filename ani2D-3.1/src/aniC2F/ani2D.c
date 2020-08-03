#include <stdlib.h>

#include "ani2D.h"


int ani2D_INIT(ani2D* ani, f77i nEStar)
{
  int mem;

  /* default parameters */
  (*ani).nv = 0;
  (*ani).nb = 0;
  (*ani).nt = 0;

  (*ani).nvfix = 0;
  (*ani).nvfix = 0;
  (*ani).ntfix = 0;

  /* default control parameters */
  (*ani).nEStar = nEStar;
  (*ani).Quality = 0.7;

  ani->control = malloc(6 * sizeof(f77i));

  (*ani).control[ANI2D_CONTROL_flagAuto] = 1;
  (*ani).control[ANI2D_CONTROL_status]   = 0;
  (*ani).control[ANI2D_CONTROL_MaxSkipE] = C2Fmax(200,  nEStar / 10);
  (*ani).control[ANI2D_CONTROL_MaxQItr]  = C2Fmax(1000, nEStar * 10);
  (*ani).control[ANI2D_CONTROL_iPrint]   = 1;
  (*ani).control[ANI2D_CONTROL_iErrMesg] = 1;

  /* default mesh memory */
  double mem_factor = 1.5;
  (*ani).ntmax = nEStar * mem_factor;
  (*ani).nbmax = nEStar * mem_factor / 4;
  (*ani).nvmax = nEStar * mem_factor / 2;
 
  return 0;
}


int ani2D_KILL(ani2D* ani)
{
  free(ani->vrt);
  free(ani->labelV);

  free(ani->bnd);
  free(ani->labelB);
  free(ani->crv);
  free(ani->labelC);

  free(ani->tri);
  free(ani->labelT);

  free(ani->control);

  return 0;
}


/* C INTERFACE TO FORTRAN ROUTINES */
int ani2D_load_mesh(ani2D* ani, const char* fname)
{
  int mem;

  loadmaniheader_(
      &ani->nv, &ani->nb, &ani->nc, &ani->nt, 
      &ani->nvfix, &ani->nbfix, &ani->ntfix, fname);

  mem = ani->nvmax = C2Fmax(ani->nvmax, ani->nv);
  ani->vrt    = malloc(2*mem * sizeof(f77r));
  ani->labelV = malloc(  mem * sizeof(f77i));
  ani->fixedV = malloc(  mem * sizeof(f77i));

  mem = ani->nbmax = C2Fmax(ani->nbmax, ani->nb);
  ani->crv    = malloc(2*mem * sizeof(f77r));
  ani->bnd    = malloc(4*mem * sizeof(f77i));
  ani->labelB = malloc(  mem * sizeof(f77i));
  ani->fixedB = malloc(  mem * sizeof(f77i));
  ani->labelC = malloc(  mem * sizeof(f77i));

  mem = ani->ntmax = C2Fmax(ani->ntmax, ani->nt);
  ani->tri    = malloc(3*mem * sizeof(f77i));
  ani->labelT = malloc(  mem * sizeof(f77i));
  ani->fixedT = malloc(  mem * sizeof(f77i));

  loadmani_(
      &ani->nv, &ani->nvfix, &ani->nvmax, ani->vrt, ani->labelV, ani->fixedV,
      &ani->nb, &ani->nbfix, &ani->nbmax, ani->bnd, ani->labelB, ani->fixedB,
      &ani->nc,                           ani->crv, ani->labelC, 
      &ani->nt, &ani->ntfix, &ani->ntmax, ani->tri, ani->labelT, ani->fixedT,
      fname);
  
  return 0;
}

void ani2D_load_mesh_header(ani2D* ani, const char* fname)
{
  loadmaniheader_(
      &ani->nv, &ani->nb, &ani->nc, &ani->nt, 
      &ani->nvfix, &ani->nbfix, &ani->ntfix, fname);
}

int ani2D_save_mesh(ani2D* ani, const char* fname)
{
  savemani_(
      &ani->nv, &ani->nvfix, ani->vrt, ani->labelV, ani->fixedV,
      &ani->nb, &ani->nbfix, ani->bnd, ani->labelB, ani->fixedB,
      &ani->nc,              ani->crv, ani->labelC, 
      &ani->nt, &ani->ntfix, ani->tri, ani->labelT, ani->fixedT,
      fname);

  return 0;
}

int ani2D_save_mesh_gmv(ani2D* ani, const char* fname)
{
  savemgmv_(&ani->nv, &ani->nt, ani->vrt, ani->tri, fname);
  return 0;
}


int ani2D_analytic(ani2D* ani, int (*metricfunc)(double*, double*, double*), 
                   void (*crvfunc)(double*, double*, int*))
{
  int     mem, *iW;
  double  *rW;
  f77i    iERR;

  /* working memory */
  mem =  6 * ani->nvmax + 10 * ani->nv 
      + 19 * ani->nbmax 
      + 11 * ani->ntmax + 12 * ani->nt + 10000;

  ani->MaxWi = mem;
  iW = malloc(mem * sizeof(f77i));

  mem =  4 * ani->nvmax + 10 * ani->nv 
      +      ani->nbmax 
      +      ani->ntmax + 10000;

  ani->MaxWr = mem;
  rW = malloc(mem * sizeof(f77r));

  /* run the code */
  mbaanalytic_(
/* M */ 
      &ani->nv, &ani->nvfix, &ani->nvmax, ani->vrt, ani->labelV, ani->fixedV,
      &ani->nb, &ani->nbfix, &ani->nbmax, ani->bnd, ani->labelB, ani->fixedB,
      &ani->nc,                           ani->crv, ani->labelC, crvfunc,
      &ani->nt, &ani->ntfix, &ani->ntmax, ani->tri, ani->labelT, ani->fixedT,
/* CONTROL */ 
      &ani->nEStar, &ani->Quality, ani->control, metricfunc, 
/* W */ 
      &ani->MaxWr, &ani->MaxWi, rW, iW, &iERR);


  /* free working memory */
  free(iW);
  free(rW);

  return iERR;
}


int ani2D_nodal(ani2D* ani, double* metric, void (*crvfunc)(double*, double*, int*))
{
  int i, mem, *iW;
  double *rW;
  f77i iERR;

  /* working memory */
  mem =  6 * ani->nvmax + 10 * ani->nv 
      +  9 * ani->nbmax 
      + 10 * ani->ntmax +  6 * ani->nt + 100000;

  ani->MaxWi = mem;
  iW = malloc(mem * sizeof(f77i));

  mem = 4 * ani->nvmax + 10 * ani->nv 
      +     ani->nbmax 
      +     ani->ntmax + 10000;

  ani->MaxWr = mem;
  rW = malloc(mem * sizeof(f77r));


  /* run the code */
  mbanodal_(
/* M */ 
      &ani->nv, &ani->nvfix, &ani->nvmax, ani->vrt, ani->labelV, ani->fixedV,
      &ani->nb, &ani->nbfix, &ani->nbmax, ani->bnd, ani->labelB, ani->fixedB,
      &ani->nc,                           ani->crv, ani->labelC, crvfunc,
      &ani->nt, &ani->ntfix, &ani->ntmax, ani->tri, ani->labelT, ani->fixedT,
/* CONTROL */ 
      &ani->nEStar, &ani->Quality, ani->control, metric, 
/* W */ 
      &ani->MaxWr, &ani->MaxWi, rW, iW, &iERR);

  /* free working memory */
  free(iW);
  free(rW);

  return iERR;
}


int  ani2D_draw_mesh(ani2D* ani, const char* fname)
{
  int  mem, *iW;

  mem = ani->nv;
  iW = calloc(mem, sizeof(f77i));

  draw_(
      &ani->nv, &ani->nb, &ani->nt, 
      ani->vrt, iW, ani->bnd, ani->labelB, ani->tri, 
      fname);

  free(iW);
  return 0;
}


/* ELEMENTAL ROUTINES */
inline int ani2D_number_points(ani2D* ani) { return ani->nv; }
inline int ani2D_number_edges(ani2D* ani) { return ani->nb; }
inline int ani2D_number_curved_edges(ani2D* ani) { return ani->nc; }
inline int ani2D_number_elements(ani2D* ani) { return ani->nt; }

int C2Fmax(int a, int b) { return (a>b)? a : b; }
int C2Fmin(int a, int b) { return (a>b)? b : a; }


/* set maximal number of points, edges and elements */
void ani2D_set_max_points(ani2D* ani, int nvmax) 
{ 
  ani->nvmax = nvmax; 
  ani->vrt = realloc(ani->vrt, 2*nvmax * sizeof(f77r));
}

void ani2D_set_max_edges(ani2D* ani, int nbmax) 
{ 
  ani->nbmax = nbmax; 
  ani->labelC = realloc(ani->labelC,   nbmax * sizeof(f77i));
  ani->crv    = realloc(ani->crv,    2*nbmax * sizeof(f77r));
  ani->bnd    = realloc(ani->bnd,    4*nbmax * sizeof(f77i));
}

void ani2D_set_max_elements(ani2D* ani, int ntmax) 
{ 
  ani->ntmax = ntmax; 
  ani->labelT = realloc(ani->labelT,   ntmax * sizeof(f77i));
  ani->tri    = realloc(ani->tri,    3*ntmax * sizeof(f77i));
}


/* operations with points */
inline void ani2D_get_point(ani2D* ani, int i, double* xy)
{
  int k;
  for(k=0; k<2; k++) xy[k] = ani->vrt[2*i+k];   
}

void ani2D_set_point(ani2D* ani, int i, double* xy)
{
  int k;
  for(k=0; k<2; k++) ani->vrt[2*i+k] = xy[k];   
}

void ani2D_fix_point(ani2D* ani, int i)
{
  int k, nvfix, flag, mem;

  flag = 0;
  nvfix = ani->nvfix;

  for(k=0; k<nvfix; k++) 
    if(ani->fixedV[k] == i) { flag = 1; break; }

  if(!flag) { 
    nvfix++; 
    ani->fixedV[nvfix] = i; 
  };
}


/* operations with boundary edges */
void ani2D_get_edge(ani2D* ani, int i, int* edge, int* icrv, int* label)
{
  int k;
  for(k=0; k<2; k++) edge[k] = ani->bnd[2*i+k];
  *icrv  = ani->labelC[i];
  *label = ani->labelB[i];
}

void ani2D_set_edge(ani2D* ani, int i, int* edge, int icrv, int label)
{
  int k;
  for(k=0; k<2; k++) ani->bnd[2*i+k] = edge[k];
  ani->labelC[i] = icrv;
  ani->labelB[i] = label;
}

void ani2D_fix_edge(ani2D* ani, int i)
{
  int k, nbfix, flag, mem;

  flag = 0;
  nbfix = ani->nbfix;

  for(k=0; k<nbfix; k++) 
    if(ani->fixedB[k] == i) { flag = 1; break; }

  if(!flag) { 
    nbfix++; 
    ani->fixedB[nbfix] = i; 
  };
}


/* operations with curved edges */
void ani2D_get_crv(ani2D* ani, int i, double* par, int* iFnc)
{
  int k;
  for(k=0; k<2; k++) par[k] = ani->crv[2*i+k];
  *iFnc = ani->labelC[i];
}

void ani2D_set_crv(ani2D* ani, int i, double* par, int iFnc)
{
  int k;
  for(k=0; k<2; k++) ani->crv[2*i+k] = par[k];
  ani->labelC[i] = iFnc;
}


/* operations with elements */
void ani2D_get_element(ani2D* ani, int i, int* tri, int* label)
{
  int k;
  for(k=0; k<3; k++) tri[k] = ani->tri[3*i+k] - 1;
  *label = ani->labelT[i];
}

void ani2D_set_element(ani2D* ani, int i, int* tri, int label)
{
  int k;
  for(k=0; k<3; k++) ani->tri[3*i+k] = tri[k];
  ani->labelT[i] = label;
}

void ani2D_fix_element(ani2D* ani, int i)
{
  int k, ntfix, flag, mem;

  flag = 0;
  ntfix = ani->ntfix;

  for(k=0; k<ntfix; k++) 
    if(ani->fixedT[k] == i) { flag = 1; break; }

  if(!flag) { 
    ntfix++; 
    ani->fixedT[ntfix] = i; 
  };
}


/* control parameers */
inline void ani2D_get_quality(ani2D* ani, double* Q) { *Q = ani->Quality; }
inline void ani2D_set_quality(ani2D* ani, double  Q) { ani->Quality = Q; }

inline void ani2D_get_status(ani2D* ani, int* status) { *status = ani->control[ANI2D_CONTROL_status]; }
inline void ani2D_set_status(ani2D* ani, int  status) { ani->control[ANI2D_CONTROL_status] = status; }

inline void ani2D_get_max_iters(ani2D* ani, int* max_iters) { *max_iters = ani->control[ANI2D_CONTROL_MaxQItr]; }
inline void ani2D_set_max_iters(ani2D* ani, int  max_iters) { ani->control[ANI2D_CONTROL_MaxQItr] = max_iters; }

inline void ani2D_get_max_basket(ani2D* ani, int* max_basket) { *max_basket = ani->control[ANI2D_CONTROL_MaxSkipE]; }
inline void ani2D_set_max_basket(ani2D* ani, int  max_basket) { ani->control[ANI2D_CONTROL_MaxSkipE] = max_basket; }

inline void ani2D_get_print_level(ani2D* ani, int* iprint) { *iprint = ani->control[ANI2D_CONTROL_iPrint]; }
inline void ani2D_set_print_level(ani2D* ani, int  iprint) { ani->control[ANI2D_CONTROL_iPrint] = iprint; }


/* triangle operations */
double triangle_area(ani2D* ani, int i) {
  int tri[3], label;
  double xy1[2], xy2[2], xy3[2], area;

  ani2D_get_element(ani, i, tri, &label);
  ani2D_get_point(ani, (int)tri[0], xy1);
  ani2D_get_point(ani, (int)tri[1], xy2);
  ani2D_get_point(ani, (int)tri[2], xy3);

  area = ((xy1[0] - xy3[0]) * (xy2[1] - xy3[1]) -
          (xy1[1] - xy3[1]) * (xy2[0] - xy3[0])) / 2;

  return area;
}
