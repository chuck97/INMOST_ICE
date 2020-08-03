#ifndef H_TREE2_MESH2D
 #define H_TREE2_MESH2D


#include"struct2.h"


typedef struct{
   PStrucNode2d root;            /*  root  of  the  quadtree  */
   PStrucFace2  *face;           /*  array of  the  faces  in  advanced  front */
   int          nFace;
   long         maxFace;         /*  number  &  max  ...  of  faces  in  ... */

   PStrucFace2  *vicinityFace;   /*  array of  faces  the  vicinity  */
   int          nVicinityFace,maxVicinityFace;   /*  number  &  max  ...  of  faces  in  ... */
   double       xVicinity,yVicinity;   /*  center  of  vicinity */
   double       sVicinity;       /*  size  of  vicinity */

   int          fill,empty;      /*  global  for  recursive  remove  function  */
   double       xc,yc,side;      /*  global  for  recursive  remove  &  insert  function  */
   double       x,y;             /*  global  for  recursive  remove  &  insert  function  */

}  StrucTree2;


#define  EMPTY  0
#define  FACE   1
#define  TREE   2


/* exported  function */
int  direction(double x,double y,double xc,double yc);
void center(double *x,double *y,double side,int d);
double  distance(double x,double y,double xc,double yc);
double distanceS(double x,double y,double xc,double yc);
double distance9(double x,double y,double xc,double yc);
PStrucFace2 addFace(int v1, int v2, int twin );
void remFace( PStrucFace2  face );
double nearest2( int *vert, double x, double y, double size );
void vicinityFaces( double x, double y, double size );


#endif



