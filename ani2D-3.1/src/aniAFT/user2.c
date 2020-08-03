#include <ctype.h>
#include <stdlib.h>
#include <math.h>


extern double  S1global;

void initUserParameters( char *pc )
{
   (void) pc;
   return;
} /* initUserParameters */

void outParameters( char *outFilename )
{
   (void) outFilename;
   return;
} /* outParameters */


/***************************************************************/
/********************    userSizeFace      *********************/
/***************************************************************/
double userSizeFace( double *xy )
{
   (void) xy;
   return  S1global;


}/*userSizeFace*/


void  userOutMesh( void )
{
   return;
}/*userOutMesh*/
