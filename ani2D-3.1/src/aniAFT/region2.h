#ifndef H_REGION2_MESH2D
#define H_REGION2_MESH2D

#include <stdio.h>
#include <string.h>
#include "struct2.h"

#define MAX_POINT 50000
#define EMPTY 0


double sizeFace(double x,double y);
void  makeAdvancedFront(int *boundVert, int *n_bound, int *cut, int *curve, double *par_t, char *outFileName);
void  boundary0( double t, double *x, double *y );
void  boundary( int n, double t, double *x, double *y );
void  initRegion( void );
void  initAddRegion( int num );
double sizeFace( double x, double y );


#endif
