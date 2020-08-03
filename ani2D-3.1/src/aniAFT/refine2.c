#include <stdlib.h>
#include "tree2.h"
#include "region2.h"
#include "refine2.h"

#include "tree2.h"
#include "region2.h"
#include "refine2.h"


extern  StrucMesh2  mesh2;
extern  int  boolRegul;
extern  int  boolFAF;

double angle( int v1, int v2, int v )
{  double p,x1,x2,y1,y2,xx,yy;
   x1=*(mesh2.x+v1);   y1=*(mesh2.y+v1);   x2=*(mesh2.x+v2);   y2=*(mesh2.y+v2);
   xx=*(mesh2.x+v);    yy=*(mesh2.y+v);    x2-=xx;  x1-=xx; y1-=yy;  y2-=yy;
   p=1./(double)sqrt( (x1*x1+y1*y1) * (x2*x2+y2*y2) );
   p*=(x1*x2+y1*y2);
   return( acos(p) );
} /* angle */

void test_quality()
{
   int i,j,p1,p2,p3,n[20];
   double pp,x1,x2,x3,y1,y2,y3,xc,yc,a1,a2,a3,da=M_PI/6,size,min=10.0,max=0.0;
   double min_angle=3.0,max_angle=0.0;
/*   FILE *f;
*/

   for(i=0;i<20;i++) n[i]=0;
   for(i=0;i<mesh2.nTria;i++){
          p1=*(mesh2.v1+i);   p2=*(mesh2.v2+i);   p3=*(mesh2.v3+i);
	  a1=angle(p2,p3,p1);  a2=angle(p3,p1,p2);  a3=angle(p1,p2,p3);
	  for(j=0;j<6;j++){
		 if( (j*da<a1) && (a1<=j*da+da) )   n[j]++ ;
		 if( (j*da<a2) && (a2<=j*da+da) )   n[j]++ ;
		 if( (j*da<a3) && (a3<=j*da+da) )   n[j]++ ;
		 if( a1>max_angle ) max_angle=a1;if( a1<min_angle ) min_angle=a1;
		 if( a2>max_angle ) max_angle=a2;if( a2<min_angle ) min_angle=a2;
		 if( a3>max_angle ) max_angle=a3;if( a3<min_angle ) min_angle=a3;
	  }
	  x1=*(mesh2.x+p1);   y1=*(mesh2.y+p1);      x2=*(mesh2.x+p2);   y2=*(mesh2.y+p2);
	  x3=*(mesh2.x+p3);   y3=*(mesh2.y+p3);
	  xc=0.3333333*(x1+x2+x3);     yc=0.3333333*(y1+y2+y3);
	  if (!boolFAF)
		  size=1.0/sizeFace(xc,yc);
	  else
		  size=3.0/(distance(x1,y1,x2,y2)+distance(x2,y2,x3,y3)+distance(x3,y3,x1,y1));
	  pp=size*distance(x1,y1,x2,y2);
	  if(pp<min) min=pp;     if(pp>max) max=pp;
	  pp=size*distance(x3,y3,x2,y2);
	  if(pp<min) min=pp;     if(pp>max) max=pp;
	  pp=size*distance(x1,y1,x3,y3);
	  if(pp<min) min=pp;     if(pp>max) max=pp;
   }/* for */

   for(i=0;i<20;i++)
     n[i] = 0;
   for(i=0;i<mesh2.nPoint;i++){
      j = mesh2.neigbor[0][i];
      if( j == -1 )
        j = 0;
      n[j]++;
   }
/* printf("Neigbor number for nP = %7d nT = %7d\n",mesh2.nPoint,mesh2.nTria);
   for(i=0;i<16;i++)
     printf("neig =%6d  nPoint =%6d  perc = %lf  perc = %lf\n",i,n[i],(double)n[i]/mesh2.nPoint,(double)n[i]/(mesh2.nPoint-n[0]));
*/
   return;
} /* test_quality */


void smoothing( void )
{
   int     i,j,n;
   double  xx,yy,x,y;

   for(i=0;i<mesh2.nPoint;i++){
      xx=0.;  yy=0.;
      n = mesh2.neigbor[0][i];
      for(j=0;j<n;j++){
	 x = mesh2.x[ mesh2.neigbor[j+1][i] ];
	 y = mesh2.y[ mesh2.neigbor[j+1][i] ];
	 xx += x;  yy += y;
      }
      if( n > 0 ){
	 xx /= n;  yy /= n;
	 mesh2.x[i] = xx;
	 mesh2.y[i] = yy;
      }
   }/* for i */

   return;
}/*smoothing*/


void delPoint( int n )
{
   mesh2.bPoint[n] = 0;

   return;
}/*delPoint*/


void changePoint( int v1, int v2 )
{
   mesh2.x[v1] = 0.5*(mesh2.x[v1] + mesh2.x[v2]);
   mesh2.y[v1] = 0.5*(mesh2.y[v1] + mesh2.y[v2]);

   return;
}/*changePoint*/


void delTria( int n )
{
   mesh2.bTria[n] = 0;

   return;
}/*delTria*/


void changeTria( int iTria, int v1, int v2, int v3 )
{
   mesh2.v1[iTria] = v1;
   mesh2.v2[iTria] = v2;
   mesh2.v3[iTria] = v3;

   return;
}/*changeTria*/


void sortNeigbor( int iVert )
{
   int     i,j,k,n,neig[MAX_NEIGBOR2];
   double  x,y,x0,y0,p,angle[MAX_NEIGBOR2];

   n = mesh2.neigbor[0][iVert];
   for(i=0;i<n;i++)
     neig[i] = mesh2.neigbor[i+1][iVert];
   x0 = mesh2.x[iVert];
   y0 = mesh2.y[iVert];
   for(i=0;i<n;i++){
      x = mesh2.x[neig[i]] - x0;
      y = mesh2.y[neig[i]] - y0;
      p = sqrt(x*x+y*y);
      if( p == 0. )  continue;
      x /= p;
      y /= p;
      if( x < 0. )
        angle[i] = 2. - y;
      else if( y >= 0. )
        angle[i] = y;
      else
        angle[i] = 4. + y;
   }

   for(i=0;i<n;i++)
     for(j=i+1;j<n;j++)
       if( angle[i] > angle[j] ){
          p = angle[i];
          k = neig[i];
          angle[i] = angle[j];
          neig[i] = neig[j];
          angle[j] = p;
          neig[j] = k;
       }

   for(i=0;i<n;i++)
     mesh2.neigbor[i+1][iVert] = neig[i];

   return;
}/*sortNeigbor*/


void calcNeigTria( void )
{
   int     i,j,v1;

   for(i=0;i<mesh2.nPoint;i++)
     mesh2.neigTria[0][i] = 0;
   for(i=0;i<mesh2.nTria;i++){
      v1 = mesh2.v1[i];
      j = mesh2.neigTria[0][v1];
      if( j >= MAX_NEIGBOR2 ) {
	  fprintf(stderr, "aniAFT: Smoothing MAX_NEIGBOR2 exceeded\n");
	  continue;
      }
      mesh2.neigTria[0][v1]++;
      mesh2.neigTria[j+1][v1] = i;

      v1 = mesh2.v2[i];
      j = mesh2.neigTria[0][v1];
      if( j >= MAX_NEIGBOR2 ) {
	  fprintf(stderr, "aniAFT: Smoothing MAX_NEIGBOR2 exceeded\n");
	  continue;
      }
      mesh2.neigTria[0][v1]++;
      mesh2.neigTria[j+1][v1] = i;

      v1 = mesh2.v3[i];
      j = mesh2.neigTria[0][v1];
      if( j >= MAX_NEIGBOR2 ) {
	  fprintf(stderr, "aniAFT: Smoothing MAX_NEIGBOR2 exceeded\n");
	  continue;
      }
      mesh2.neigTria[0][v1]++;
      mesh2.neigTria[j+1][v1] = i;
   }

   return;
}/*calcNeigTria*/


void calcNeigbor( void )
{
   int     i,j,k,n,iTria,vert[3*MAX_NEIGBOR2];

   for(i=0;i<mesh2.nPoint;i++){
      if( mesh2.neigbor[0][i] == -1 )
        continue;
      mesh2.neigbor[0][i] = 0;
   }
   for(i=0;i<mesh2.nPoint;i++){
      if( mesh2.neigbor[0][i] == -1 )
        continue;
      n = 0;
      for(j=0;j<mesh2.neigTria[0][i];j++){
         iTria = mesh2.neigTria[j+1][i];
         vert[3*j+0] = mesh2.v1[iTria];
         vert[3*j+1] = mesh2.v2[iTria];
         vert[3*j+2] = mesh2.v3[iTria];
      }
      for(j=0;j<3*mesh2.neigTria[0][i];j++){
         if( vert[j] == i )
           continue;
         iTria = 0;
         for(k=0;k<n;k++)
           if( vert[j] == mesh2.neigbor[k+1][i] )
             iTria = 1;
         if( iTria )
           continue;
         n++;
         mesh2.neigbor[n][i] = vert[j];
      }
      mesh2.neigbor[0][i] = n;
   }

   return;
}/*calcNeigbor*/


void calcEdge( void )
{
   int     i,iTria,vert,j,k,n,bool;
   unsigned long p;

   p = (unsigned long)3*mesh2.maxPoint*sizeof(int);
   mesh2.vb = (int Huge *)farmalloc(p);
   mesh2.ve = (int Huge *)farmalloc(p);
   mesh2.tria1 = (int Huge *)farmalloc(p);
   mesh2.tria2 = (int Huge *)farmalloc(p);
   if( mesh2.vb == NULL || mesh2.vb == NULL || mesh2.tria1 == NULL || mesh2.tria2 == NULL ) {
       fprintf(stderr, "aniAFT: out of memory\n");
       return;
   }

   mesh2.nEdge = 0;
   for(i=0;i<mesh2.nPoint;i++){
      n = mesh2.neigbor[0][i];
      if( n == -1 )
        continue;
      for(j=0;j<n;j++){
         vert = mesh2.neigbor[j+1][i];
         if( vert <= i )
           continue;
         if( mesh2.neigbor[0][vert] == -1 )
           continue;

         if( mesh2.nEdge >= 3*mesh2.maxPoint ) {
	     fprintf(stderr, "Smootinh: nT >= 3 * nV\n");
	 }
         mesh2.vb[mesh2.nEdge] = i;
         mesh2.ve[mesh2.nEdge] = vert;
         bool = 1;
         for(k=0;k<n;k++){
            iTria = mesh2.neigTria[k+1][i];
            if( mesh2.v1[iTria] == vert || mesh2.v2[iTria] == vert || mesh2.v3[iTria] == vert ){
               if( bool ){
                  mesh2.tria1[mesh2.nEdge] = iTria;
                  bool = 0;
               }
               else
                  mesh2.tria2[mesh2.nEdge] = iTria;
            }
         }
         bool = 1;
         iTria = mesh2.tria1[mesh2.nEdge];
         if( mesh2.v1[iTria] == i && mesh2.v2[iTria] == vert )
           bool = 0;
         else if( mesh2.v2[iTria] == i && mesh2.v3[iTria] == vert )
           bool = 0;
         else if( mesh2.v3[iTria] == i && mesh2.v1[iTria] == vert )
           bool = 0;
         if( bool ){
            mesh2.tria1[mesh2.nEdge] = mesh2.tria2[mesh2.nEdge];
            mesh2.tria2[mesh2.nEdge] = iTria;
         } 
         mesh2.nEdge++;
      }
   }

   return;
}/*calcEdge*/


void initRegularity( void )
{
   int  i;
   unsigned long p;

   p = (unsigned long)mesh2.maxPoint*sizeof(char);
   mesh2.bPoint = (char Huge *)farmalloc(p);
   if( mesh2.bPoint == NULL ) {
       fprintf(stderr, "aniAFT: out of memory\n");
       return;
   }
   p = (unsigned long)mesh2.maxTria*sizeof(char);
   mesh2.bTria = (char Huge *)farmalloc(p);
   if( mesh2.bTria == NULL ) {
       fprintf(stderr, "aniAFT: out of memory\n");
       return;
   }

   for(i=0;i<mesh2.nPoint;i++)
     mesh2.bPoint[i] = 1;
   for(i=mesh2.nPoint;i<mesh2.maxPoint;i++)
     mesh2.bPoint[i] = 0;
   for(i=0;i<mesh2.nTria;i++)
     mesh2.bTria[i] = 1;
   for(i=mesh2.nTria;i<mesh2.maxTria;i++)
     mesh2.bTria[i] = 0;

   mesh2.neigTria = farmalloc( (MAX_NEIGBOR2+1)*sizeof(int Huge *) );
   if( mesh2.neigTria == NULL ) {
       fprintf(stderr, "aniAFT: out of memory\n");
       return;
   }
   for(i=0;i<=MAX_NEIGBOR2;i++){
      mesh2.neigTria[i] = (int Huge *) farmalloc(mesh2.maxPoint*sizeof(int));
      if( mesh2.neigTria[i] == NULL ) {
	  fprintf(stderr, "aniAFT: out of memory\n");
	  return;
      }
   }
   calcNeigTria();

   return;
}/*initRegularity*/


void pack( void )
{
   int  i,j;

   for(i=0;i<mesh2.nTria;i++){
      if( !mesh2.bTria[i] ){
         mesh2.nTria--;
         mesh2.bTria[i] = mesh2.bTria[mesh2.nTria];
         mesh2.v1[i] = mesh2.v1[mesh2.nTria];
         mesh2.v2[i] = mesh2.v2[mesh2.nTria];
         mesh2.v3[i] = mesh2.v3[mesh2.nTria];
         i--;
      }
   }
   for(i=0;i<mesh2.nPoint;i++){
      if( !mesh2.bPoint[i] ){
         mesh2.nPoint--;
         mesh2.bPoint[i] = mesh2.bPoint[mesh2.nPoint];
         mesh2.x[i] = mesh2.x[mesh2.nPoint];
         mesh2.y[i] = mesh2.y[mesh2.nPoint];
         for(j=0;j<mesh2.nTria;j++){
            if( mesh2.v1[j] == mesh2.nPoint )
              mesh2.v1[j] = i;
            if( mesh2.v2[j] == mesh2.nPoint )
              mesh2.v2[j] = i;
            if( mesh2.v3[j] == mesh2.nPoint )
              mesh2.v3[j] = i;
         }
         i--;
      }
   }

   return;
}/*pack*/


void deletePoint( void )
{
   int  i,j,n1,n2,n3,n4,sum,swap;

   for(i=0;i<mesh2.nPoint;i++){
      j = mesh2.neigbor[0][i];
      if( j == 3 ){
         delPoint(i);
         sortNeigbor(i);
         changeTria(mesh2.neigTria[1][i],
                    mesh2.neigbor[2][i],mesh2.neigbor[1][i],mesh2.neigbor[3][i]);
         delTria(mesh2.neigTria[2][i]);
         delTria(mesh2.neigTria[3][i]);
      }
      else if( j == 4 ){
         delPoint(i);
         sortNeigbor(i);
         n1 = mesh2.neigbor[0][ mesh2.neigbor[1][i] ] - 7;
         n2 = mesh2.neigbor[0][ mesh2.neigbor[2][i] ] - 7;
         n3 = mesh2.neigbor[0][ mesh2.neigbor[3][i] ] - 7;
         n4 = mesh2.neigbor[0][ mesh2.neigbor[4][i] ] - 7;
         sum = (n1+1)*(n1+1) + (n3+1)*(n3+1) + n2*n2 + n4*n4;
         swap = (n2+1)*(n2+1) + (n4+1)*(n4+1) + n1*n1 + n3*n3;
/*printf("swap %4d%4d n(%4d%4d%4d%4d) \n",sum,swap,n1,n3,n2,n4);*/
         if( swap < sum ){
            changeTria(mesh2.neigTria[1][i],
                       mesh2.neigbor[2][i],mesh2.neigbor[1][i],mesh2.neigbor[4][i]);
            changeTria(mesh2.neigTria[2][i],
                       mesh2.neigbor[4][i],mesh2.neigbor[3][i],mesh2.neigbor[2][i]);
         }
         else{
            changeTria(mesh2.neigTria[1][i],
                       mesh2.neigbor[2][i],mesh2.neigbor[1][i],mesh2.neigbor[3][i]);
            changeTria(mesh2.neigTria[2][i],
                       mesh2.neigbor[4][i],mesh2.neigbor[3][i],mesh2.neigbor[1][i]);
         }
         delTria(mesh2.neigTria[3][i]);
         delTria(mesh2.neigTria[4][i]);
      }
   }
   pack();
   calcNeigTria();
   calcNeigbor();

   return;
}/*deletePoint*/


void mergePoint( void )
{
   int  i,iTria,j,v1,v2,n1,n2;

   calcEdge();
   for(i=0;i<mesh2.nEdge;i++){
      v1 = mesh2.vb[i];
      v2 = mesh2.ve[i];
      n1 = mesh2.neigbor[0][v1];
      n2 = mesh2.neigbor[0][v2];
      if( n1 == 5 && n2 == 5 ){
         mesh2.neigbor[0][v1] = 6;
         mesh2.neigbor[0][v2] = 6;
         changePoint(v1,v2);
         delPoint(v2);
         for(j=0;j<mesh2.neigTria[0][v2];j++){
            iTria = mesh2.neigTria[j+1][v2];
            if( mesh2.v1[iTria] == v2 )
              mesh2.v1[iTria] = v1;
            if( mesh2.v2[iTria] == v2 )
              mesh2.v2[iTria] = v1;
            if( mesh2.v3[iTria] == v2 )
              mesh2.v3[iTria] = v1;
         }
         delTria(mesh2.tria1[i]);
         delTria(mesh2.tria2[i]);
      }
   }
   free(mesh2.vb);
   free(mesh2.ve);
   free(mesh2.tria1);
   free(mesh2.tria2);

   pack();
   calcNeigTria();
   calcNeigbor();

   return;
}/*mergePoint*/


void splitPoint( void )
{
   int     i,j,vert,n1,n2,lab,it;
   double  size;

   for(i=0;i<mesh2.nPoint;i++){
      n1 = mesh2.neigbor[0][i];
      if( n1 <= 7 )
        continue;
      n2 = n1/2;
      n1 = n1 - n2;
      sortNeigbor(i);

	  if (!boolFAF)
		  size = sizeFace(mesh2.x[i],mesh2.y[i]);
	  else
		  size = distance(mesh2.x[i],mesh2.y[i],
				  mesh2.x[mesh2.neigbor[n2][i]],mesh2.y[mesh2.neigbor[n2][i]]);
      vert = mesh2.nPoint;
      mesh2.neigbor[0][vert] = 1;
      addPoint(mesh2.x[i]+0.1*size,mesh2.y[i]);
      for(j=1;j<=n1+n2;j++){
         if( mesh2.neigbor[0][ mesh2.neigbor[j][i] ] != -1 )
           mesh2.neigbor[0][ mesh2.neigbor[j][i] ] = 1;
      }
lab=-1;
      for(j=1;j<n1;j++){
it=mesh2.neigTria[j][i]; if (lab==-1) lab=mesh2.label[it];
if (lab != mesh2.label[it])  fprintf(stderr, "aniAFT: different materials inside region\n");
         changeTria(mesh2.neigTria[j][i],
                    mesh2.neigbor[j+1][i],mesh2.neigbor[j][i],i);
      }
      for(j=n1+1;j<n1+n2;j++){
it=mesh2.neigTria[j][i]; if (lab==-1) lab=mesh2.label[it]; 
if (lab != mesh2.label[it])  fprintf(stderr, "aniAFT: different materials inside region\n");
         changeTria(mesh2.neigTria[j][i],
                    mesh2.neigbor[j+1][i],mesh2.neigbor[j][i],vert);
      }
it=mesh2.neigTria[n1][i]; if (lab==-1) lab=mesh2.label[it]; 
if (lab != mesh2.label[it])  fprintf(stderr, "aniAFT: different materials inside region\n");
      changeTria(mesh2.neigTria[n1][i],
                 mesh2.neigbor[n1][i],i,vert);
it=mesh2.neigTria[n1+n2][i]; if (lab==-1) lab=mesh2.label[it]; 
if (lab != mesh2.label[it])  fprintf(stderr, "aniAFT: different materials inside region\n");
      changeTria(mesh2.neigTria[n1+n2][i],
                 mesh2.neigbor[n1+1][i],mesh2.neigbor[n1][i],vert);
      addTria(mesh2.neigbor[n1+n2][i],vert,i,lab);
      addTria(mesh2.neigbor[n1+n2][i],i,mesh2.neigbor[1][i],lab);
   }
/*   pack();*/
   calcNeigTria();
   calcNeigbor();
   for(i=0;i<5;i++)
     smoothing();

   return;
}/*splitPoint*/


void swapEdge( void )
{
   int  i,vb,ve,v1,v2,nb,ne,n1,n2,tria1,tria2,sum,swap;

   calcEdge();
   for(i=0;i<mesh2.nEdge;i++){
      vb = mesh2.vb[i];
      ve = mesh2.ve[i];
      tria1 = mesh2.tria1[i];
      tria2 = mesh2.tria2[i];
      v1 = mesh2.v1[tria1];
      if( v1 == vb || v1 == ve ){
         v1 = mesh2.v2[tria1];
         if( v1 == vb || v1 == ve ){
            v1 = mesh2.v3[tria1];
         }
      }
      v2 = mesh2.v1[tria2];
      if( v2 == vb || v2 == ve ){
         v2 = mesh2.v2[tria2];
         if( v2 == vb || v2 == ve ){
            v2 = mesh2.v3[tria2];
         }
      }
      if( mesh2.neigbor[0][v1] == -1 || mesh2.neigbor[0][v2] == -1 )
        continue;
      if( mesh2.neigTria[0][vb] == -1 && mesh2.neigTria[0][ve] == -1 )
        continue;
      nb = mesh2.neigbor[0][vb] - 6;
      ne = mesh2.neigbor[0][ve] - 6;
      n1 = mesh2.neigbor[0][v1] - 6;
      n2 = mesh2.neigbor[0][v2] - 6;
      if( nb < 1 && ne < 1 && n1 < 1 && n2 < 1 )
        continue;
      sum = nb*nb + ne*ne + n1*n1 + n2*n2;
      swap = (nb-1)*(nb-1) + (ne-1)*(ne-1) + (n1+1)*(n1+1) + (n2+1)*(n2+1);
/*printf("swap %4d%4d v(%4d%4d%4d%4d) n(%4d%4d%4d%4d) \n",sum,swap,vb,ve,v1,v2,nb,ne,n1,n2);*/
      if( swap < sum ){
         mesh2.neigbor[0][vb] = nb + 5;
         mesh2.neigbor[0][ve] = ne + 5;
         mesh2.neigbor[0][v1] = n1 + 7;
         mesh2.neigbor[0][v2] = n2 + 7;
         mesh2.neigTria[0][vb] = -1;
         mesh2.neigTria[0][ve] = -1;
         mesh2.neigTria[0][v1] = -1;
         mesh2.neigTria[0][v2] = -1;
         changeTria(tria1,v2,v1,vb);
         changeTria(tria2,v1,v2,ve);
      }
   }
   free(mesh2.vb);
   free(mesh2.ve);
   free(mesh2.tria1);
   free(mesh2.tria2);

   calcNeigTria();
   calcNeigbor();

   return;
}/*swapEdge*/


void regularity( void )
{
   int  i;

if( boolRegul ){

   initRegularity();
   for(i=0;i<5;i++)
     smoothing();
/*test_quality();
getch();*/

   deletePoint();
   for(i=0;i<5;i++)
     smoothing();
/*test_quality();
getch();*/

   mergePoint();
   mergePoint();
   mergePoint();
   for(i=0;i<5;i++)
     smoothing();
/*test_quality();
getch();*/

   splitPoint();
   splitPoint();
   splitPoint();
/*test_quality();
printf("SWAP:\n");
getch();*/

   swapEdge();
   for(i=0;i<5;i++)
     smoothing();
   swapEdge();
   for(i=0;i<5;i++)
     smoothing();
   swapEdge();
   for(i=0;i<5;i++)
     smoothing();
   swapEdge();

/*test_quality();
getch();*/

}/*if( boolRegul )*/

   for(i=0;i<mesh2.nSmooth;i++)
     smoothing();

   return;
}/*regularity*/
