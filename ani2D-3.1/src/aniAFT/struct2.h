#ifndef H_STRUCT2_MESH2D
#define H_STRUCT2_MESH2D

#define getch       getchar
#define farfree     free
#define farmalloc   malloc
#define Huge
#define Far
#define Near

#define  GRAF     0


typedef struct{
   int     nPoint,nTria,nSmooth;
   int     nRegion,nRLine,iRLine; /**/
   int     **boundVert;
   int     *nVert,*region,*bCut,*nRPoint,*nRTria;
   long    maxPoint,maxTria;             /**/
   double  Huge *x;
   double  Huge *y;
   int     Huge *v1;
   int     Huge *v2;
   int     Huge *v3;
   int     Huge *label;
   int Huge *( Huge * neigbor ); /**/

   char    Huge *bPoint;
   char    Huge *bTria;
   int Huge *( Huge * neigTria ); /**/
   int     Huge *vb;
   int     Huge *ve;
   int     Huge *tria1;
   int     Huge *tria2;
   int     nEdge;

   double  size;
   char    outFileName[128];
   char    inFileName[128];
   char    parFileName[128];
   char    debFileName[128];
   char    debug;
}  StrucMesh2;


typedef struct{
   int     v1,v2;     /**/
   int     f;          /* number  of  face       */
   double  x,y,s;    /**/
}  StrucFace2;
typedef  StrucFace2  *PStrucFace2;


typedef struct _MtEC{
    PStrucFace2  face;
    struct _MtEC *next;
} entrychain;

typedef struct _MtSN2{
    int           entrycount;
    entrychain    *firstentry;
    struct _MtSN2 *parent;
    struct _MtSN2 *nodelist[4];
} StrucNode2d;
typedef  StrucNode2d  *PStrucNode2d;

typedef struct {
    int           flag;
    union {
	PStrucFace2 faceNode;
	int v;
    } u;
} StrucNode2;
typedef  StrucNode2  *PStrucNode2;


typedef  PStrucNode2  StrucList2[4];
typedef  StrucList2   *PStrucList2;

#define  MAX_NEIGBOR2   16
#define  S_StrucList2  sizeof(StrucList2)
#define  S_StrucNode2  sizeof(StrucNode2)
#define  S_StrucFace2  sizeof(StrucFace2)

#define  S_StrucNode2d  sizeof(StrucNode2d)
#define  S_entrychain  sizeof(entrychain)


/* exported  function  function */
void init( void );
void addPoint( double x, double y );
void addTria( int v1, int v2, int v3, int lab );
void outMesh( void );


#endif
