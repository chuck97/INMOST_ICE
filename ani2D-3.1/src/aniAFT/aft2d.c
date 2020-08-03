#include <stdlib.h>
#include "struct2.h"
#include "user2.h"
#include "region2.h"
#include "tria2.h"
#include "memory2.h"

extern  StrucMesh2  mesh2;
unsigned int _stklen=24000u;

double ReferenceCrd[2];
double ScalingFactor;

double *bvglobal;
double *bltailglobal;
int    *blglobal;
int    nVertglobal, nLineglobal;

int    nTrglobal, *trglobal;
int    nBrglobal, *brglobal;
int    nVrglobal;
double *vrglobal, *vrbrglobal;

int    boolFAFglobal, mesh2nSmoothglobal;
double SIglobal, S1global, SMglobal;

int    nVRTglobal;
double *vrtglobal;
int    nTRIglobal, *triglobal, *labtriglobal, nBNDglobal, *bndglobal;
int    nCRVglobal;
double *crvglobal;
int    *iFNCglobal;
int    StopAfterinitRegion;

void  userBoundaryDummy(int *i, double *t, double *x, double *y)
{
   (void) i,  (void) t,  (void) x,  (void) y;
   return;
}/*userBoundaryDummy*/

typedef void (*userfn_t) (int *, double *, double *, double *);
userfn_t userfn = userBoundaryDummy;

void registeruserfn_(void * p)
{
   userfn = (userfn_t) p ;
   return;
}/*registeruserfn*/

typedef double (*sizefn_t) (double *);
sizefn_t sizefn = userSizeFace;

void registersizefn_(void * p)
{
   sizefn = (sizefn_t) p ;
   return;
}/*registersizefn*/

int aft2dfront_( 
             int *pnBr, int *br, int *pnVr, double *vrbr, 
             int *pnVRT, double *vrt, 
             int *pnTRI, int *tri, int *labtri, 
             int *pnBND, int *bnd, int *labbnd ) 
{
   int i,j, last,k, err;
   double xmax,ymax;

/*  Input to makeTria()  */


   nBrglobal = *pnBr;
   nVrglobal = *pnVr;
 if (nBrglobal) {
	   brglobal = (int*) malloc(sizeof(int)*(2*nBrglobal));
	   vrbrglobal = (double*) malloc(sizeof(double)*(2*nVrglobal));

	   for( i = 0; i < nBrglobal; i++) 
	      for( j = 0; j < 2; j++) brglobal[i*2+j] = br[i*2+j];
	   for( i = 0; i < nVrglobal; i++) 
	      for( j = 0; j < 2; j++) vrbrglobal[i*2+j] = vrbr[i*2+j];

 } else {
	   nBrglobal = nVrglobal;
	   brglobal = (int*) malloc(sizeof(int)*(2*nBrglobal));
	   vrbrglobal = (double*) malloc(sizeof(double)*(2*nVrglobal));
	   last=-1; k=0;
	   for( i = 0; i < nBrglobal; i++) {
		   brglobal[k*2] = i+1;
		   brglobal[k*2+1] = i+2;
		   for( j = 0; j < 2; j++) vrbrglobal[i*2+j] = vrbr[i*2+j];
		   if (last>=0 && vrbr[last*2]==vrbr[i*2] && 
				   vrbr[last*2+1]==vrbr[i*2+1]) {
			   brglobal[(k-1)*2+1] = last+1;
			   last = -1;
		   } else {
			   k++;
			   if (last<0) last = i;
		   }
	   }
	   nBrglobal = k;
   }

/* scaling to the unit square */
   ReferenceCrd[0] = vrbrglobal[0];
   ReferenceCrd[1] = vrbrglobal[1];
   xmax            = vrbrglobal[0];
   ymax            = vrbrglobal[1];
   for( i = 1; i < nVrglobal; i++) {
       if ( ReferenceCrd[0] > vrbrglobal[i*2+0] ) 
            ReferenceCrd[0] = vrbrglobal[i*2+0];                    
       if ( ReferenceCrd[1] > vrbrglobal[i*2+1] ) 
            ReferenceCrd[1] = vrbrglobal[i*2+1];                    
       if ( xmax < vrbrglobal[i*2+0] ) 
            xmax = vrbrglobal[i*2+0];                    
       if ( ymax < vrbrglobal[i*2+1] ) 
            ymax = vrbrglobal[i*2+1];                    
   }
   if ( xmax-ReferenceCrd[0] < ymax - ReferenceCrd[1] ) {
      ScalingFactor = 1.0 / ( ymax - ReferenceCrd[1] );
   } else {
      ScalingFactor = 1.0 / ( xmax - ReferenceCrd[0] );
   }
   for( i = 0; i < nVrglobal; i++) {
      for( j = 0; j < 2; j++) vrbrglobal[i*2+j] = (vrbrglobal[i*2+j]-ReferenceCrd[j]) * ScalingFactor;
   }

   boolFAFglobal   = 1;
   S1global        = 0.1 ;
   SIglobal        = 0.2;
   SMglobal        = 0.35;
   mesh2nSmoothglobal  = 5;
   StopAfterinitRegion = 0;


   err = makeTria();

   *pnVRT = nVRTglobal;
   *pnTRI = nTRIglobal;
   *pnBND = nBNDglobal;
/* scaling back */
   for(i=0;i<nVRTglobal;i++){
          vrt[2*i]   = vrtglobal[2*i]   / ScalingFactor + ReferenceCrd[0];
          vrt[2*i+1] = vrtglobal[2*i+1] / ScalingFactor + ReferenceCrd[1];
   }
   for(i=0;i<nTRIglobal;i++){
          tri[3*i]   = triglobal[3*i] ;
          tri[3*i+1] = triglobal[3*i+1];
          tri[3*i+2] = triglobal[3*i+2];
          labtri[i]  = labtriglobal[i];
   }
   for(i=0;i<nBNDglobal;i++){
          bnd[2*i]   = bndglobal[5*i] ;
          bnd[2*i+1] = bndglobal[5*i+1];
        /*bnd[4*i+2] = bndglobal[5*i+2];*/
          labbnd[i]  = bndglobal[5*i+4];
   }
   
   freeMemory();

   free(brglobal);
   free(vrbrglobal);

   free(vrtglobal);
   free(bndglobal);
   free(crvglobal);
   free(iFNCglobal);
   if ( StopAfterinitRegion == 0 ) free(triglobal);
   if ( StopAfterinitRegion == 0 ) free(labtriglobal);

   return err;
}/*aft2dfront*/



int aft2dboundary_( int *pnVert, double *bv,  
             int *pnLine, int *bl, double *bltail, double *hsze,
             int *pnVRT, double *vrt, 
             int *pnTRI, int *tri, int *labtri, 
             int *pnBND, int *bnd, int *labbnd,
             int *pnCRV, double *crv, int *iFNC )  
{ 
   int i,j, err;
   double xmax,ymax;

/*  Input to makeTria()  */

   nVertglobal = *pnVert;
   bvglobal = (double*) malloc(sizeof(double)*(2*nVertglobal));

   for( i = 0; i < nVertglobal; i++) {
      bvglobal[2*i  ] = bv[2*i  ];
      bvglobal[2*i+1] = bv[2*i+1];
   }

   nLineglobal = *pnLine;
   blglobal = (int*) malloc(sizeof(int)*(7*nLineglobal));
   bltailglobal = (double*) malloc(sizeof(double)*(2*nLineglobal));

   for( i = 0; i < nLineglobal; i++) {
      for( j = 0; j < 7; j++) blglobal[i*7+j] = bl[i*7+j];
      for( j = 0; j < 2; j++) bltailglobal[i*2+j] = bltail[i*2+j];
   }

/* scaling to the unit square */
   ReferenceCrd[0] = bvglobal[0];
   ReferenceCrd[1] = bvglobal[1];
   xmax            = bvglobal[0];
   ymax            = bvglobal[1];
   for( i = 1; i < nVertglobal; i++) {
       if ( ReferenceCrd[0] > bvglobal[i*2+0] ) 
            ReferenceCrd[0] = bvglobal[i*2+0];                    
       if ( ReferenceCrd[1] > bvglobal[i*2+1] ) 
            ReferenceCrd[1] = bvglobal[i*2+1];                    
       if ( xmax < bvglobal[i*2+0] ) 
            xmax = bvglobal[i*2+0];                    
       if ( ymax < bvglobal[i*2+1] ) 
            ymax = bvglobal[i*2+1];                    
   }
   if ( xmax-ReferenceCrd[0] < ymax - ReferenceCrd[1] ) {
      ScalingFactor = 1.0 / ( ymax - ReferenceCrd[1] );
   } else {
      ScalingFactor = 1.0 / ( xmax - ReferenceCrd[0] );
   }
   for( i = 0; i < nVertglobal; i++) {
      for( j = 0; j < 2; j++) bvglobal[i*2+j] = (bvglobal[i*2+j]-ReferenceCrd[j]) * ScalingFactor;
   }


   boolFAFglobal = 0;
   S1global        = hsze[0];
   mesh2nSmoothglobal = 5;
   StopAfterinitRegion = 0;


   err = makeTria();

   *pnVRT = nVRTglobal;
   *pnTRI = nTRIglobal;
   *pnBND = nBNDglobal;
   *pnCRV = nCRVglobal;
/* scaling back */
   for(i=0;i<nVRTglobal;i++){
          vrt[2*i]   = vrtglobal[2*i]   / ScalingFactor + ReferenceCrd[0];
          vrt[2*i+1] = vrtglobal[2*i+1] / ScalingFactor + ReferenceCrd[1];
   }
   for(i=0;i<nTRIglobal;i++){
          tri[3*i]   = triglobal[3*i] ;
          tri[3*i+1] = triglobal[3*i+1];
          tri[3*i+2] = triglobal[3*i+2];
          labtri[i]  = labtriglobal[i];
   }
   for(i=0;i<nBNDglobal;i++){
          bnd[2*i]   = bndglobal[5*i] ;
          bnd[2*i+1] = bndglobal[5*i+1];
       /* bnd[2*i+2] = bndglobal[5*i+2]; */
          labbnd[i]  = bndglobal[5*i+4];

          if( (j = bndglobal[5*i+2]) > 0 ) {
                 j--;
                 crv[2*i]   = crvglobal[2*j] ;
                 crv[2*i+1] = crvglobal[2*j+1];
                 iFNC[i]    = iFNCglobal[j];
          } else {
                 iFNC[i] = 0;
          }
   }

   /*
   for(i=0;i<nCRVglobal;i++){
          crv[2*i]   = crvglobal[2*i] ;
          crv[2*i+1] = crvglobal[2*i+1];
          iFNC[i]    = iFNCglobal[i];
   }
   */

   freeMemory();

   free(bvglobal);
   free(blglobal);
   free(bltailglobal);

   free(trglobal);
   free(vrglobal);
   free(brglobal);
   free(vrbrglobal);

   free(vrtglobal);
   free(bndglobal);
   free(crvglobal);
   free(iFNCglobal);
   if ( StopAfterinitRegion == 0 ) free(triglobal);
   if ( StopAfterinitRegion == 0 ) free(labtriglobal);

   return err;
}/*aft2dboundary*/


/* Added for backward compatibility */
#define  GE_Memory     1
#define  GE_Critical   2
#define  GE_User       3
#define  GE_Temp       4

void errorExit2(int group, const char *number)
{
    switch(group){
	case  GE_Memory:
	    printf("\nOut of memory!\n");
	    break;
	case  GE_Critical:
	    printf("\nCritical error in program!\n");
	    break;
	case  GE_User:
	    printf("\nError in user function!\n");
	    break;
	case  GE_Temp:
	    printf("\nTempoparal error!\n");
	    break;
    }
    printf("User message: %s\n",number);
    exit(1);
    return;
}/* error */


