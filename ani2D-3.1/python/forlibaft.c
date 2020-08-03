#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>


void errorExit2(int group, char *number);


void boundaryTopHose(double t, double *x, double *y)
{
   x[0] = -2.0 + 4.0 * t;   /* 0 <= t <= 1.0 */
   y[0] =  1.0 + 0.5 * atan(8*t-4) / atan(4.0);

   return;
}


void boundaryBottomHose(double t, double *x, double *y)
{
   x[0] = -2.0 + 4.0 * t;  /* 0 <= t <= 1.0 */
   y[0] = -1.0 - 0.5 * atan(8*t-4) / atan(4.0);


   return;
}


void boundaryTopCircle(double t, double *x, double *y)
{
   double PI = 4*atan(1.0); /*3.14159265358979323846;*/
   double R  = 0.5;
 
   x[0] = R * cos( PI*t ); /* 0 <= t <= 1.0 */
   y[0] = R * sin( PI*t );

   return;
}


void boundaryBottomCircle(double t, double *x, double *y)
{
   double PI = 4*atan(1.0); /*3.14159265358979323846;*/
   double R  = 0.5;
 
   x[0] = -R * cos( PI*t ); /* 0 <= t <= 1.0 */
   y[0] = -R * sin( PI*t );

   return;
}


void userboundary_(int *i, double *t, double *x, double *y)
{
   switch (*i) {
      case 1:  
         boundaryTopHose(*t, x, y); break;
      case 2:  
         boundaryBottomHose(*t, x, y); break;
      case 3:  
         boundaryTopCircle(*t, x, y); break;
      case 4:  
         boundaryBottomCircle(*t, x, y); break; 
      default:
	 errorExit2(3, "Wrong id for the curvilinear boundary\n");
   }

   return;
}/*userboundary*/

