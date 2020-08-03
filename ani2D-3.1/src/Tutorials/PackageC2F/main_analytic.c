#include  <stdio.h>
#include  <stdlib.h>

#include  "ani2D.h"


int metric_user(double*, double*, double*);
void crvfunc_user(double*, double*, int*);


int main() 
{
/* local variables */
  int    i;
  int    nv, nb, nt;
  int    edge[2], tri[3], icrv, ifnc, label, max_iters;
  double xy[2], Q, Lp;

/* create a pointer to the main structure */
  ani2D  ani;

/* mesh control parameters */
  int nEStar = 1000;
  ani2D_INIT(&ani, nEStar);


/* load a mesh from a file */
  ani2D_load_mesh(&ani, "../data/simple.ani");


/* example of accessing mesh points */
  printf("C2F: usage of FORTRAN library animba2D-3.x.a in a C-code\n");

  nv = ani2D_number_points(&ani); 
  printf("     Number of points is%5d  (we print two)\n", nv);

  for(i=0; i<2; i++) {
     ani2D_get_point(&ani, i, xy);
     printf("%8d %6.3f %6.3f\n", i, xy[0], xy[1]);
  }


/* example of accessing boundary edges */
  nb = ani2D_number_edges(&ani); 
  printf("     Number of boundary edges is%5d  (we print two)\n", nb);

  for(i=0; i<2; i++) {
     ani2D_get_edge(&ani, i, edge, &icrv, &label);

     printf("%8d  pts: %3d %3d   label =%2d", i, edge[0], edge[1], label);
     if(icrv > 0) printf("(curliniear edge %4d)", icrv); 
     printf("\n");
  }


/* example of accessing elements */
  nt = ani2D_number_elements(&ani); 
  printf("     Number of triangles is%5d  (we print two)\n", nt);
  
  for(i=0; i<2; i++) {
     ani2D_get_element(&ani, i, tri, &label);
     printf("%8d  pts: %3d %3d %3d   label =%2d\n", i, tri[0], tri[1], tri[2], label);
  }


/* control parameter */
  Q = 0.7;
  ani2D_set_quality(&ani, Q);
  i = 5000;
  ani2D_set_max_iters(&ani, i);


/* draw initial mesh.  The name must terminate with .ps*/
  ani2D_draw_mesh(&ani, "initial_mesh.ps");


/* adapt the mesh */
  ani2D_analytic(&ani, &metric_user, &crvfunc_user);


/* save the mesh in a file */
  ani2D_save_mesh(&ani, "../data/save.ani");


/* draw final mesh. The name must terminate with .ps */
  ani2D_draw_mesh(&ani, "final_mesh.ps");


/* kill the mesh */
  ani2D_KILL(&ani);

  exit(0);
}


/* USER WRITTEN ROUTINES */
/*
  Three Fortran routines below create a metric field which
  is 2x2 variable positive definite symmetric tensor,
                | f_(x,y)  h_(x,y) |
      M(x, y) = |                  |
                | h_(x,y)  g_(x,y) |

*/
int metric_user(double* x, double* y, double* M) 
{ 
  double f = 2, g = 1, h = 1;

  M[0] = f;
  M[1] = h;
  M[2] = g;
  M[3] = h;

  return 0; 
}

/*
  The routine computes the Cartesian coordinates of point
  xyc from its parametric coordinate tc.

  tc    - the given parametric coordinate of point
  xy(2) - the Cartesian coordinate of the same point
  iFnc  - the function number for computing

  On input :  tc, iFnc
  On output:  xyc(2)
*/

void crvfunc_user(double* tc, double* xy, int* iFnc) {};


