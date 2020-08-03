#include<ctype.h>
#include<stdlib.h>
#include<math.h>


void errorExit2(int group, char *number);

static double X1=0.30,X2=0.70,XE=1.,Y2=0.5;
static double C1=0.17735, C2=-0.075597, C3=-0.212836, C4=0.17363, C5=-0.062547;
static double wingAngle=0; 


void boundaryWing(double *t, double *x, double *y)
{
  int sign;
  double tt, xx, yy;
  double CX2;

  CX2 = 0.5 * (X1+X2);

  tt = 2 * (*t);
  if(tt > 1.0) {
    tt -= 1.0;
    sign = -1;
  }
  else
  {
    tt = 1.0 - tt;
    sign = 1;
  }
  xx = X1 + tt * (X2-X1);
  yy = Y2 + sign * (X2-X1) * (C1*sqrt(tt) + C2*tt + C3*tt*tt + C4*tt*tt*tt + C5*tt*tt*tt*tt);

  x[0] = (CX2 + (xx-CX2) * cos(wingAngle) + (yy-Y2) * sin(wingAngle));
  y[0] = (Y2  - (xx-CX2) * sin(wingAngle) + (yy-Y2) * cos(wingAngle));

  return;
}  // boundaryWing



void boundarySlit(double *t, double *x, double *y)
{
  double w, XB;

  w = 0.5 * (X2-X1);
  XB = X2 - w * (1.0 - cos(wingAngle));

  x[0] = XB + (XE-XB) * (*t);
  y[0] = Y2 + w * sin(-wingAngle);

  return;
}  // boundarySlit



void userboundary_(int *i, double *t, double *x, double *y)
{
  switch (*i) {
    case 1:
      boundaryWing(t, x, y);
      break;
    case 2:
      boundarySlit(t, x, y);
      break;
    default:
       errorExit2(3," i is wrong in boundary");
   }

   return;
}  // userboundary

