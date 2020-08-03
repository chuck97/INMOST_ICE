#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include"memory2.h"
#include"tree2.h"


/* extern  variables */
extern StrucMesh2  mesh2;
extern int boolFAF;

/* global  variables */
StrucTree2  tree2;


/* static  function */
static void addFaceArray( PStrucFace2  face );
static void remFaceArray( int  face );
static int  sectQuad( void );
static void vicinityFacesRec( PStrucNode2d  node );


int direction(double x,double y,double xc,double yc)
{
   int  dx=0,dy=0;
   if( x > xc )  dx = 1;   if( y > yc )  dy = 2;
   return( dx+dy );
}


void center(double *x,double *y,double side,int d)
{
   int  dx=0,dy=0;

   if( (d<0) || (d>3) )  d = 0;
   dy = d/2;  dx = d%2;

   if( dx )  x[0] += side;   else  x[0] -= side;
   if( dy )  y[0] += side;   else  y[0] -= side;

   return;
}


double distance(double x,double y,double xc,double yc)
{  return(sqrt( (x-xc)*(x-xc)+(y-yc)*(y-yc) ) );
}/*distance*/

double distanceS(double x,double y,double xc,double yc)
{  return( (x-xc)*(x-xc)+(y-yc)*(y-yc) );
}/*distanceS*/

double distance9(double x,double y,double xc,double yc)
{
   double  p1,p2;

   p1 = fabs(x-xc);
   p2 = fabs(y-yc);
   if( p1 < p2 )
     p1 = p2;

   return  p1;
} /*distance9*/


static void addFaceArray( PStrucFace2  face )
{
   int son,fath,i;
   PStrucFace2  cface,*tface=tree2.face;
   if( tree2.nFace >= tree2.maxFace ) {
       fprintf(stderr, "aniAFT: no memory for new face in front\n");
       return;
   }

   tface[tree2.nFace] = face;  face->f = tree2.nFace;
   son=++tree2.nFace;
aaa:fath=son/2;
   if( (fath>0) && (tface[son-1]->s < tface[fath-1]->s)  ){
      cface=tface[son-1];
      tface[son-1]=tface[fath-1];
      tface[fath-1]=cface;
      i=cface->f;  cface->f=tface[son-1]->f;  tface[son-1]->f=i;
      son=fath; if(son>1) goto aaa;
   }

   return;
} /*addFaceArray*/


static void remFaceArray( int  fath )
{
   int son1,son2,i;
   PStrucFace2 face,*tface;
   tface = tree2.face;
   if( (fath >= tree2.nFace) || (fath<0) )  {
       fprintf(stderr, "aniAFT: delete face: tree is empty\n");
       return;
   }

   tree2.nFace--;
   tface[fath]=tface[tree2.nFace];
   tface[fath]->f=fath;
aaa:i = 0;
   if( 2*fath+1 < (tree2.nFace-1) ){
      son1=2*fath+1;  son2=son1++;
      if((tface[fath]->s<=tface[son1]->s)&&(tface[fath]->s<=tface[son2]->s)) i=0;
      if((tface[fath]->s>=tface[son1]->s)&&(tface[son1]->s>=tface[son2]->s)) i=son2;
      if((tface[fath]->s>=tface[son2]->s)&&(tface[son2]->s>tface[son1]->s)) i=son1;
      if((tface[son1]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son2]->s)) i=son2;
      if((tface[son2]->s>=tface[fath]->s)&&(tface[fath]->s>tface[son1]->s)) i=son1;
   }
   if( 2*fath+1 == (tree2.nFace-1) ){
      son1=2*fath+1;
      if((tface[fath]->s>tface[son1]->s)) i=son1;else i=0;
   }

   if(i)
   {  face=tface[i];  tface[i]=tface[fath];   tface[fath]=face;
      son1=face->f;  face->f=tface[i]->f;  tface[i]->f=son1;
      fath=i;   goto aaa;
   }

   return;
} /*remFaceArray*/


PStrucFace2 addFace(int v1, int v2, int twin )
{
	int         i,d;
	double      x,y,x1=-1.0,y1=-1.0;
	double      xc=0.5,yc=0.5,side=0.5;
	PStrucFace2 face;
	PStrucNode2d old_node,node;
	entrychain  *faceentry, *new_faceentry;

	node=tree2.root;
	face = nearAlloc( S_StrucFace2 );
	face->v1 = v1;   face->v2 = v2;
	if( twin == 0 ){
		x = (mesh2.x[v1]+mesh2.x[v2])/2.;
		y = (mesh2.y[v1]+mesh2.y[v2])/2.;
	} else {
		x = 0.49*mesh2.x[v1]+0.51*mesh2.x[v2];
		y = 0.49*mesh2.y[v1]+0.51*mesh2.y[v2];
	}
	face->x=x;
	face->y=y;
	face->s = distance(mesh2.x[v1],mesh2.y[v1],mesh2.x[v2],mesh2.y[v2]);
	if( (x<0.)||(x>1.)||(y<0.)||(y>1) ) {
		printf( "x= %lf,  y= %lf ", x, y );
		fprintf(stderr, "aniAFT: out of bouding box (wrong orientation in input?)\n");
	}
	while (side>face->s) {
		d = direction(x,y,xc,yc);
		side *= 0.5;
		center(&xc,&yc,side,d);
		old_node = node;
		node = node->nodelist[d];
		if (!node) {
			node = nearAlloc( S_StrucNode2d );
			for (i=0;i<4;i++)
				node->nodelist[i] = 0;
			node->parent = old_node;
			node->entrycount = 0;
			node->firstentry = 0;
			old_node->nodelist[d] = node;
		}
	}
	faceentry = node->firstentry;
	for (i=0;i<node->entrycount;i++) {
		x1 = faceentry->face->x;
		y1 = faceentry->face->y;
		if ( (x1==x)&&(y1==y)&&(v1==faceentry->face->v1)&&(v2==faceentry->face->v2) ) {
		    	printf("x,  y:  %lf, %lf\nx1, y1: %lf, %lf\nv1:   %d, v2:   %d\n v1_f: %d, v2_f: %d\n", x, y, x1, y1, v1, v2, faceentry->face->v1, faceentry->face->v2);
			fprintf(stderr, "aniAFT: dups in front (corrupted front?)\n");
		}
		faceentry = faceentry->next;
	}
	faceentry = node->firstentry;
	new_faceentry = nearAlloc( S_entrychain );
	node->firstentry = new_faceentry;
	new_faceentry->face = face;
	new_faceentry->next = faceentry;
	node->entrycount++;
	addFaceArray( face );
	
	return face;
} /*addFace*/

void remFace( PStrucFace2  face )
{
	int         i,d;
	double      x,y;
	double      xc=0.5,yc=0.5,side=0.5;
	PStrucNode2d node, old_node;
	entrychain  *faceentry, *old_faceentry;

	remFaceArray(face->f);
	x = face->x;
	y = face->y;
	node = tree2.root;
	while (side>face->s) {
		d = direction(x,y,xc,yc);
		side *= 0.5;
		center(&xc,&yc,side,d);
		node = node->nodelist[d];
		if (!node) {
		    fprintf(stderr, "aniAFT: delete face: no such element\n");
		    return;
		}
	}
	faceentry = node->firstentry;
	if (faceentry->face == face) {
		node->firstentry = faceentry->next;
		free(face);
		free(faceentry);
		node->entrycount--;
	} else {
		do {
			old_faceentry = faceentry;
			faceentry = faceentry->next;
		} while ( (faceentry)&&(faceentry->face != face) );
		if (faceentry) {
			old_faceentry->next = faceentry->next;
			free(face);
			free(faceentry);
			node->entrycount--;
		} else {
		    fprintf(stderr, "aniAFT: delete face: no such element\n");
		    return;
		}
	}
	while ( (node->entrycount==0)&&(!node->nodelist[0]&&!node->nodelist[1]&&
				!node->nodelist[2]&&!node->nodelist[3])&&(node!=tree2.root) ) {
		old_node = node;
		node = node->parent;
		for (i=0;i<4;i++)
			if (node->nodelist[i]==old_node) break;
		node->nodelist[i] = 0;
		free(old_node);
	}
	return;
} /* remFace */

double nearest2( int *vert, double x, double y, double size )
{
   int          i,j,vn=0,v[2];
   double       p=0,dist=2.0;
   PStrucFace2  face;

   vicinityFaces(x,y,size);
   for(i=0;i<tree2.nVicinityFace;i++){
      face = tree2.vicinityFace[i];
      v[0] = face->v1;  v[1] = face->v2;
      for(j=0;j<2;j++){
	 p = distance(mesh2.x[v[j]],mesh2.y[v[j]],x,y);
	 if(p < dist){
	    dist = p;
	    vn = v[j];
	 }
      }/* for  2  vert  of  face */
   }/* for  i < tree2.nVicinityFace */
   vert[0] = vn;

   return( dist );
} /*nearest2*/


static int sectQuad( void )
/*  need  Quad     : tree2.xc;  tree2.yc;
		     tree2.side;
	  Vicinity : tree2.xVicinity;  tree2.yVicinity;
		     tree2.sVicinity;
*/
{
   if( tree2.xc+3*tree2.side < tree2.xVicinity-tree2.sVicinity )  return(0);
   if( tree2.xc-3*tree2.side > tree2.xVicinity+tree2.sVicinity )  return(0);
   if( tree2.yc+3*tree2.side < tree2.yVicinity-tree2.sVicinity )  return(0);
   if( tree2.yc-3*tree2.side > tree2.yVicinity+tree2.sVicinity )  return(0);

   return(1);

} /*sectQuad*/


static void vicinityFacesRec( PStrucNode2d  node )
/*  NEED    tree2.xc=0.5;tree2.yc=0.5;tree2.side=0.5;
            tree2.nVicinityFace=0;
            tree2.sVicinity , tree2.xVicinity , tree2.yVicinity ,
*/
{
	int         i;
	double      x,y,s;
	entrychain  *faceentry;
	
	if (!node) return;
	
	for (faceentry = node->firstentry;faceentry;faceentry=faceentry->next) {
		x = faceentry->face->x;
		y = faceentry->face->y;
		s = distance(x,y,tree2.xVicinity,tree2.yVicinity)-0.5*faceentry->face->s;
		if (s<=tree2.sVicinity) {
			tree2.vicinityFace[tree2.nVicinityFace] = faceentry->face;
			tree2.nVicinityFace++;
			if (tree2.nVicinityFace > tree2.maxVicinityFace) {
			    fprintf(stderr, "aniAFT: partial result in tree search\n");
			    return;
			}
		}
	}
	tree2.side *= 0.5;
	x = tree2.xc;
	y = tree2.yc;
	s = tree2.side;
	if (node->nodelist[0]) {
	    center(&tree2.xc, &tree2.yc, tree2.side, 0);
	    if (sectQuad()) vicinityFacesRec( node->nodelist[0] );
	}
	for (i=1;i<4;i++) {
		if (!node->nodelist[i]) continue;
		tree2.side = s;
		tree2.xc = x;
		tree2.yc = y;
		center(&tree2.xc, &tree2.yc, tree2.side, i);
		if (sectQuad()) vicinityFacesRec( node->nodelist[i] );
	}
	return;
} /* vicinityFacesRec */

void vicinityFaces( double x, double y, double size )
{
	tree2.xc=0.5; tree2.yc=0.5; tree2.side=0.5;
	tree2.nVicinityFace = 0;
	tree2.sVicinity = size;
	tree2.xVicinity = x; tree2.yVicinity = y;
	vicinityFacesRec( tree2.root );

	return;
} /* vicinityFaces */

