#include<ctype.h>
#include<stdlib.h>
#include<math.h>


void errorExit2(int group,char *number);


void boundaryCirc(double *t, double *x, double *y)
{
   double PI = 4*atan(1.0);  // 3.14159265358979323846
   x[0] = cos(PI * (*t));  // 0 <= t <= 1.5 
   y[0] = sin(PI * (*t));

   return;
}  // boundaryCirc



void userboundary_(int *i, double *t, double *x, double *y)
{
  switch (*i) {
    case 1:
      boundaryCirc(t, x, y);
      break;
    default:
      errorExit2(3," i is wrong   in  boundary ");
  }

  return;
}  // userboundary

