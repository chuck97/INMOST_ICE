#ifndef H_USER2_MESH2D
#define H_USER2_MESH2D

#include"struct2.h"

/* exported  function  function */
void initUserParameters(char *pc);
void outParameters(char *outFilename);
double userSizeFace(double *xy);
void userOutMesh(void);
void userBoundary(int *i, double *t, double *x, double *y);

#endif








