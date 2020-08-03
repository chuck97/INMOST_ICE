#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "struct2.h"
#include "tree2.h"
#include "memory2.h"
#include "user2.h"
#include "region2.h"

StrucMesh2  mesh2;
extern StrucTree2  tree2;
extern int boolFAF;
int boolSize=0,boolNav=0,boolRegul=0;

extern int boolFAFglobal, mesh2nSmoothglobal;

void initParameter( void)
{
   char  *dummy=0;

   boolFAF = boolFAFglobal;
   boolSize = 1;
   mesh2.nSmooth = mesh2nSmoothglobal;
   initUserParameters(dummy);
   return;
}  /*initParameter*/


void init( void )
{
/* printf(" Generator  of  unstructured  2-D  meshes.\n");
   printf(" (c) Copyright  1993-95.\n");
   printf(" LabNumMath ( INM ).\n");
   printf(" Version July, 18, 1995.\n");
*/
   mesh2.nPoint = 0;
   mesh2.nTria = 0;
   mesh2.nRegion = 1;
   mesh2.nRLine = 0;
   mesh2.iRLine = 0;
   mesh2.nSmooth = 5;
   strcpy(mesh2.outFileName,"");
   strcpy(mesh2.inFileName,"");
   strcpy(mesh2.debFileName,"debug");

   mesh2.debug = 1;
   tree2.nFace = 0;

   initMemory();
   initParameter();

//   tree2.root->flag = EMPTY;
	tree2.root->entrycount = 0;
	tree2.root->parent = 0;
	tree2.root->nodelist[0] = 0;
	tree2.root->nodelist[1] = 0;
	tree2.root->nodelist[2] = 0;
	tree2.root->nodelist[3] = 0;
	tree2.root->firstentry = 0;
   return;
}/*init*/


void addPoint( double x, double y )
{
   if( mesh2.nPoint >= mesh2.maxPoint ) {
       fprintf(stderr, "aniAFT: mesh2.nPoint >= mesh2.maxPoint\n");
       return;
   }
   mesh2.neigbor[0][mesh2.nPoint] = 0;
   mesh2.x[mesh2.nPoint] = x;
   mesh2.y[mesh2.nPoint] = y;
   mesh2.nPoint++;

   return;
}/*addPoint*/


void addTria( int v1, int v2, int v3, int lab )
{
   if( mesh2.nTria >= mesh2.maxTria ) {
       fprintf(stderr, "aniAFT: mesh2.nTria >= mesh2.maxTria\n");
       return;
   }
   mesh2.v1[mesh2.nTria] = v1;
   mesh2.v2[mesh2.nTria] = v2;
   mesh2.v3[mesh2.nTria] = v3;
   mesh2.label[mesh2.nTria] = lab;
   mesh2.nTria++;

   return;
}/*addTria*/


extern int    nVRTglobal;
extern double *vrtglobal;
extern int    nTRIglobal, *triglobal, *labtriglobal;

void outMesh( void )
{
   int vv1,vv2,vv3,i;
   double x,y/*,s*/;
//   FILE *f;
//   char name[128];
//   char ext_p[]=".ps",ext_t[]=".tri",ext_v[]=".vrt";

   nVRTglobal = mesh2.nPoint;
   vrtglobal = (double*) malloc(sizeof(double)*(2*nVRTglobal));

   for(i=0;i<mesh2.nPoint;i++){
          x=*(mesh2.x+i);   y=*(mesh2.y+i);
          vrtglobal[2*i] = x;
          vrtglobal[2*i+1] = y;
   }

   nTRIglobal = mesh2.nTria;
   triglobal = (int*) malloc(sizeof(int)*(3*nTRIglobal));
   labtriglobal = (int*) malloc(sizeof(int)*(nTRIglobal));
   for(i=0;i<mesh2.nTria;i++){
          vv1=1+*(mesh2.v1+i);   vv2=1+*(mesh2.v2+i);   vv3=1+*(mesh2.v3+i);
          triglobal[3*i] = vv1;
          triglobal[3*i+1] = vv2;
          triglobal[3*i+2] = vv3;
          labtriglobal[i] = mesh2.label[i];
   }

   if( boolNav )
     userOutMesh();

   return;
}/*outMesh*/





