//#include <malloc.h>  // deprecated on Mac OSX ??
#include <stdlib.h>
#include <stdio.h>

#include "memory2.h"
#include "tree2.h"


/* extern  variables */
extern  StrucMesh2  mesh2;
extern  StrucTree2  tree2;


/* static  function */
static char Huge *hugeAccess( char Huge * p, unsigned long n );
static void initFarMemory( void );
static void initNearMemory( void );


static char Huge *hugeAccess( char Huge * p, unsigned long n )
{
   p += n;
   return( p );
} /*hugeAccess*/


void Near *nearAlloc( unsigned int n )
{
   void Near *p;

   p = malloc(n);
   if (p == NULL) perror("nearAlloc");

   return( p );
} /*nearAlloc*/


char Huge *    ppMemory;
PStrucFace2  *ptree2face;

static void initFarMemory( void )
{
   int            i,allSize=0,size[100];
   unsigned long  maxMemory=0,allMemory=0,sMemory[100];
   char Huge *    pMemory=NULL;


   mesh2.neigbor = farmalloc( (MAX_NEIGBOR2+1)*sizeof(int Huge *) );
   if( mesh2.neigbor == NULL )  perror("farMemory");

   maxMemory = 48000000;

   if( maxMemory == 0 )  perror("farMemory");

/* init  distribution  of  memory */
   size[0] = 1 * ( 2*sizeof(double) + (MAX_NEIGBOR2+1)*sizeof(int) );
   size[1] = 3 *   4*sizeof(int);

   for(i=0;i<2;i++)  allSize += size[i];
   for(i=0;i<2;i++){
      sMemory[i] = maxMemory*size[i]/allSize;
      allMemory += sMemory[i];
   }
/* end  distribution  of  memory */

/* init  memory  ptr */
   pMemory = farmalloc( allMemory );
   ppMemory = pMemory;
   if( pMemory == NULL )  perror("farMemory");

   mesh2.maxPoint = sMemory[0]/( 2*sizeof(double) + (MAX_NEIGBOR2+1)*sizeof(int) );
   mesh2.x = (double Huge *)pMemory;
   pMemory = hugeAccess( pMemory, mesh2.maxPoint*sizeof(double) );
   mesh2.y = (double Huge *)pMemory;
   pMemory = hugeAccess( pMemory, mesh2.maxPoint*sizeof(double) );
   for(i=0;i<=MAX_NEIGBOR2;i++){
      mesh2.neigbor[i] =(int Huge *) pMemory;
      pMemory = hugeAccess( pMemory, mesh2.maxPoint*sizeof(int) );
   }

   mesh2.maxTria = sMemory[1]/( 4*sizeof(int) );
   mesh2.v1 = (int Huge *)pMemory;
   pMemory  = hugeAccess( pMemory, mesh2.maxTria*sizeof(int) );
   mesh2.v2 = (int Huge *)pMemory;
   pMemory  = hugeAccess( pMemory, mesh2.maxTria*sizeof(int) );
   mesh2.v3 = (int Huge *)pMemory;
   pMemory  = hugeAccess( pMemory, mesh2.maxTria*sizeof(int) );
   mesh2.label = (int Huge *)pMemory;
   pMemory  = hugeAccess( pMemory, mesh2.maxTria*sizeof(int) );
/* end  init  memory  ptr */

   return;
}/*initFarMemory*/


static void initNearMemory( void )
{
   tree2.root = nearAlloc( S_StrucNode2d );
   tree2.maxFace = mesh2.maxPoint/5;
   tree2.face = nearAlloc( tree2.maxFace*sizeof(PStrucFace2) );
   ptree2face = tree2.face;
   tree2.maxVicinityFace = 200;
   tree2.vicinityFace = nearAlloc( tree2.maxVicinityFace*sizeof(StrucFace2) );

   return;
}/*initNearMemory*/


void initMemory( void )
{
   initFarMemory();
   initNearMemory();

   return;
}/*initMemory*/

void freeMemory(void) {
    int i;
    free(mesh2.neigbor),  mesh2.neigbor = NULL;
    free(tree2.vicinityFace),  tree2.vicinityFace = NULL;
    free(mesh2.nRTria),  mesh2.nRTria = NULL;
    free(mesh2.nRPoint),  mesh2.nRPoint = NULL;
    free(mesh2.bCut),  mesh2.bCut = NULL;
    free(mesh2.region),  mesh2.region = NULL;
    free(mesh2.nVert),  mesh2.nVert = NULL;
    for (i=0; i<mesh2.iRLine; i++)  free(mesh2.boundVert[i]);
    free(mesh2.boundVert),  mesh2.boundVert = NULL;
    while (tree2.nFace > 0)  remFace(tree2.face[0]);
    free(tree2.root),  tree2.root = NULL;
}
