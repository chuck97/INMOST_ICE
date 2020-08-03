//to call floating-point check from the Fortran routine
#include <math.h>

#ifndef isinf
#include <ieeefp.h>
int isinf(double x) { return !finite(x) && x==x; }
#endif


void fpcheck_(double *a, int *flag) 
{ 
   *flag = isnan(*a) && !isinf(*a); 
}

