#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>

#include  "ani2D.h"


int main(int argc, char *argv[]) 
{
/* local variables */
  int i, n;
  int nv, nb, nc, nt;
  int edge[2], tri[3], icrv, ifnc, label;
  double xy1[2], xy2[2], xy3[2], xy[2], xymin[2], xymax[2];
  double area, area_min, area_max, area_avg, tmp;
  double rR, rIn, rOut, rRmin, rRmax, rRavg;
  double len, len_min, len_max, len_avg;
  const char *FILE_NAME;

/* create a pointer to the main structure */
  ani2D ani;
  int nEStar;

  if (argc < 2) {
    //FILE_NAME = "../data/simple.ani";
    FILE_NAME = "../data/country.ani";
  }
  else {
    FILE_NAME = argv[1];
  }


/* load and draw a mesh from a file */
  nEStar = 1000;
  ani2D_INIT(&ani, nEStar);

  ani2D_load_mesh(&ani, FILE_NAME);
  ani2D_draw_mesh(&ani, "initial_mesh.ps");


/* general information about the mesh */
  nv = ani2D_number_points(&ani); 
  nt = ani2D_number_elements(&ani); 
  nb = ani2D_number_edges(&ani); 
  nc = ani2D_number_curved_edges(&ani); 
  
  printf("Mesh audit...\n");
  printf("  Number of points:         %5d \n", nv);
  printf("  Number of triangles:      %5d \n", nt);
  printf("  Number of boundary edges: %5d \n", nb);
  printf("  Number of curved edges:   %5d \n", nc);


/* domain size */
  for(i=0; i<nv; i++) {
    ani2D_get_point(&ani, i, xy);
    if(i==0) {
      xymin[0] = xymax[0] = xy[0];
      xymin[1] = xymax[1] = xy[1];
    }
    else {
      if(xymin[0] > xy[0]) xymin[0] = xy[0];
      if(xymin[1] > xy[1]) xymin[1] = xy[1];

      if(xymax[0] < xy[0]) xymax[0] = xy[0];
      if(xymax[1] < xy[1]) xymax[1] = xy[1];
    }
  }
  printf("\n  Bounding box [%-lg, %-lg] x [%-lg, %-lg]\n", xymin[0], xymax[0], xymin[1], xymax[1]);
  

/* total area */
  area = 0.0;
  area_min = 1e+24;
  area_max = 0.0;
  for(i=0; i<nt; i++) {
    tmp = fabs(triangle_area(&ani, i));
    area += tmp;
    if (area_min > tmp)  area_min = tmp;
    if (area_max < tmp)  area_max = tmp;
  }
  area_avg = area / nt;

  printf("  Total area = %-lg\n", area);
  printf("    min area = %-lg\n", area_min);
  printf("    max area = %-lg\n", area_max);
  printf("    avg area = %-lg\n", area_avg);


/* finite element regularity */
  rRavg = 0.0;
  rRmin = 1e+24;
  rRmax = 0.0;
  for (i=0; i<nt; i++) {
    ani2D_get_element(&ani, i, tri, &label);

    ani2D_get_point(&ani, tri[0], xy1);
    ani2D_get_point(&ani, tri[1], xy2);
    ani2D_get_point(&ani, tri[2], xy3);

    randr_(xy1, xy2, xy3, &rOut, &rIn);

    rR = rOut / rIn;
    if (rRmax < rR)  rRmax = rR;
    if (rRmin > rR)  rRmin = rR;
    rRavg += rR;
  }
  rRavg = rRavg / nt;

  printf("\n  Triangles statistics\n");
  printf("    min aspect ratio = %-lg\n", rRmin);
  printf("    max aspect ratio = %-lg\n", rRmax);
  printf("    avg aspect ratio = %-lg\n", rRavg);

 
/* edges */
  len_min = 1e+24;
  len_avg = 0;
  for (i=0; i<nt; i++) {
    ani2D_get_element(&ani, i, tri, &label);

    ani2D_get_point(&ani, tri[0], xy1);
    ani2D_get_point(&ani, tri[1], xy2);
    ani2D_get_point(&ani, tri[2], xy3);
 
    len = sqrt(pow(xy2[0] - xy1[0], 2.0) + pow(xy2[1] - xy1[1], 2.0)); 
    len_avg += len;
    if( len_max < len) len_max = len;
    if( len_min > len) len_min = len;
 
    len = sqrt(pow(xy3[0] - xy1[0], 2.0) + pow(xy3[1] - xy1[1], 2.0)); 
    len_avg += len;
    if( len_max < len) len_max = len;
    if( len_min > len) len_min = len;
 
    len = sqrt(pow(xy3[0] - xy2[0], 2.0) + pow(xy3[1] - xy2[1], 2.0)); 
    len_avg += len;
    if( len_max < len) len_max = len;
    if( len_min > len) len_min = len;
  }
  len_avg /= (3*nt);

  printf("\n  Edge distribution\n");
  printf("    min edge length = %-lg\n", len_min);
  printf("    max edge length = %-lg\n", len_max);
  printf("    avg edge length = %-lg\n", len_avg);


/* kill the mesh */
  ani2D_KILL(&ani);
  exit(0);
}



