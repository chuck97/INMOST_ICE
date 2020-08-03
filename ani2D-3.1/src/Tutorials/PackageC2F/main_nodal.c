#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>

#include  "ani2D.h"

void crvfunc_user( double*, double*, int* );


int main() 
{
/* local variables */
  int     i, n;
  int     nv, nb, nt;
  int     edge[2], tri[3], icrv, ifnc, label, max_iters;
  double  xy[2], Q, Lp;

/* create a pointer to the main structure */
  ani2D   ani;

/* mesh control parameters */
  int     nEStar;
  double  *metric;


/* print information about this program */
  printf("\nThe program illustrates use of FORTRAN library animba2D-3.x.a in a C-code\n\n");


/* initialize a mesh with nEStar triangles */
  nEStar = 4000;
  ani2D_INIT(&ani, nEStar);


/* load a mesh from a file */
  ani2D_load_mesh(&ani, "../data/square.ani");


/* allocate more memory than it is set by default */
  ani2D_set_max_points(&ani, 7000);
  //ani2D_set_max_elements(&ani, 12000);


/* example of accessing mesh points */
  nv = ani2D_number_points(&ani); 
  printf("\nNumber of points is%5d  (we print two)\n", nv);

  for( i=0; i<2; i++ ) {
     ani2D_get_point(&ani, i, xy);
     printf("%4d %6.3f %6.3f\n", i, xy[0], xy[1]);
  }


/* example of accessing boundary edges */
  nb = ani2D_number_edges(&ani); 
  printf("\nNumber of boundary edges is%5d  (we print two)\n", nb);

  for( i=0; i<2; i++ ) {
     ani2D_get_edge(&ani, i, edge, &icrv, &label);

     printf("%4d  pts: %4d %4d   label =%2d", i, edge[0], edge[1], label);
     if( icrv > 0 ) printf(" (curvilinear edge #%2d)", icrv); 
     printf("\n");
  }


/* example of accessing elements */
  nt = ani2D_number_elements(&ani); 
  printf("\nNumber of triangles is%5d  (we print two)\n", nt);
  
  for( i=0; i<2; i++ ) {
     ani2D_get_element(&ani, i, tri, &label);
     printf("%4d  pts: %4d %4d %4d   label =%2d\n", i, tri[0], tri[1], tri[2], label);
  }


/* control parameter */
  Q = 0.7;
  ani2D_set_quality(&ani, Q);
  ani2D_get_max_iters(&ani, &max_iters);


/* draw ninitial mesh. The name must terminate with .ps */
  ani2D_draw_mesh(&ani, "initial_mesh.ps");


/* create an artificial diagonal metric */
  nv = ani2D_number_points(&ani);
  metric = malloc(3 * nv * sizeof(double));

  for( n=0,i=0; i<nv; i++) {
    ani2D_get_point(&ani, i, xy);
    metric[n++] = 1 + xy[1] * xy[1];   /* M_11 = 1 + y^2 */ 
    metric[n++] = 2 + xy[0] * xy[0];   /* M_22 = 2 + x^2 */ 
    metric[n++] = 0;                   /* M_12 = M_21 = 0 */ 
  }
  

/* adapt the mesh to solution sol */
  ani2D_nodal(&ani, metric, &crvfunc_user);


/* draw final mesh. The name must terminate with .ps */
  ani2D_draw_mesh(&ani, "final_mesh.ps");


/* save the mesh in a file */
  ani2D_save_mesh(&ani, "../data/save.ani");


/* kill the mesh */
  ani2D_KILL(&ani);
  exit(0);
}



/* USER WRITTEN ROUTINES */
/*
  The routine computes the Cartesian coordinates of point
  xyc from its parametric coordinate tc.

  tc    - the given parametric coordinate of point
  xy(2) - the Cartesian coordinate of the same point
  iFnc  - the function number for computing

  On input :  tc, iFnc
  On output:  xyc(2)
*/
void crvfunc_user(double* tc, double* xy, int* iFnc) 
{
  double t = *tc;

  if(*iFnc == 100) {
    xy[0] = 0.2 - 2 * t * (0.3 - t);
    xy[1] = t;

  } else if(*iFnc == 200) {
    xy[0] = 1 - 2 * (1 - t) * (t - 0.7);
    xy[1] = t;
  
  } else {
    printf("Undefined function = %d\n", *iFnc);
    exit(1);
  }
}


