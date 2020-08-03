#include<math.h>
#include<stdlib.h>

#include"region2.h"
#include"refine2.h"
#include"tree2.h"
#include"tria2.h"
#include"user2.h"

#define SHOWPROGRESS 0

extern  StrucMesh2  mesh2;
extern  StrucTree2  tree2;

extern  int boolFAF;

double COARSEFACTOR = 1.5;

static int     new_vert;

static double minrho,alpha,beta;

static PStrucFace2 intedge;
static int ntodelete;
static PStrucFace2 todelete[3];

typedef int edge[2];
typedef double vertex[3];
typedef struct {
    int *ia;
    int *ja;
} neigh_edges;

typedef struct {
    int *ia;
    edge *ja;
} neigh_trias;

static neigh_edges eadj = {NULL, NULL};
static neigh_trias tadj = {NULL, NULL};


static double dist(int  a, int b) {
    return  sqrt((mesh2.x[a]-mesh2.x[b])*(mesh2.x[a]-mesh2.x[b]) + (mesh2.y[a]-mesh2.y[b])*(mesh2.y[a]-mesh2.y[b]));
}

static double func_q(int k1, int k2, int k) {
    double L, S;
    S = (mesh2.x[k2]*mesh2.y[k1] - mesh2.x[k1]*mesh2.y[k2] + mesh2.x[k]*(mesh2.y[k2]-mesh2.y[k1]) + mesh2.y[k]*(mesh2.x[k1]-mesh2.x[k2]))/2.0;
    L = 2.0*mesh2.x[k]*mesh2.x[k] + 2.0*mesh2.x[k1]*mesh2.x[k1] + 2.0*mesh2.x[k2]*mesh2.x[k2] - 2.0*mesh2.x[k]*mesh2.x[k1] - 2.0*mesh2.x[k]*mesh2.x[k2] - 2.0*mesh2.x[k1]*mesh2.x[k2] +
	2.0*mesh2.y[k]*mesh2.y[k] + 2.0*mesh2.y[k1]*mesh2.y[k1] + 2.0*mesh2.y[k2]*mesh2.y[k2] - 2.0*mesh2.y[k]*mesh2.y[k1] - 2.0*mesh2.y[k]*mesh2.y[k2] - 2.0*mesh2.y[k1]*mesh2.y[k2];
    return  4.0*sqrt(3.0)*S/L;
}
static int func_xy(int k, int k1, int k2, double dx[2], double delta, int n) {
    double L, S, Lx, Ly, H, Hx, Hy, f, r;
    n = 8;
    S = (mesh2.x[k2]*mesh2.y[k1] - mesh2.x[k1]*mesh2.y[k2] + mesh2.x[k]*(mesh2.y[k2]-mesh2.y[k1]) + mesh2.y[k]*(mesh2.x[k1]-mesh2.x[k2]))/2.0;
    H = S + sqrt(S*S + 4.0*delta*delta);
    //    H = 2.0*S;
    L = 2.0*mesh2.x[k]*mesh2.x[k] + 2.0*mesh2.x[k1]*mesh2.x[k1] + 2.0*mesh2.x[k2]*mesh2.x[k2] - 2.0*mesh2.x[k]*mesh2.x[k1] - 2.0*mesh2.x[k]*mesh2.x[k2] - 2.0*mesh2.x[k1]*mesh2.x[k2] +
    2.0*mesh2.y[k]*mesh2.y[k] + 2.0*mesh2.y[k1]*mesh2.y[k1] + 2.0*mesh2.y[k2]*mesh2.y[k2] - 2.0*mesh2.y[k]*mesh2.y[k1] - 2.0*mesh2.y[k]*mesh2.y[k2] - 2.0*mesh2.y[k1]*mesh2.y[k2];
    Lx = 4.0*mesh2.x[k] - 2.0*(mesh2.x[k1]+mesh2.x[k2]);
    Ly = 4.0*mesh2.y[k] - 2.0*(mesh2.y[k1]+mesh2.y[k2]);
    Hx = (mesh2.y[k2]-mesh2.y[k1])/2.0 + (S*(mesh2.y[k2]-mesh2.y[k1]))/2.0/sqrt(S*S + 4.0*delta*delta);
    Hy = (mesh2.x[k1]-mesh2.x[k2])/2.0 + (S*(mesh2.x[k1]-mesh2.x[k2]))/2.0/sqrt(S*S + 4.0*delta*delta);
    //    Hx = (mesh2.y[k2]-mesh2.y[k1]);
    //    Hy = (mesh2.x[k1]-mesh2.x[k2]);
    f = 2.0*L/H/4.0/sqrt(3.0);
    dx[0] = (2.0*Lx*H - 2.0*Hx*L)/H/H/4.0/sqrt(3.0);
    dx[1] = (2.0*Ly*H - 2.0*Hy*L)/H/H/4.0/sqrt(3.0);
    while (n>1) {
	dx[0] *= f;
	dx[1] *= f;
	n--;
    }
    /* 
      r = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
      if (r>1.0) {
      r = (1.0 + log(r))/r;
      dx[0] *= r,  dx[1] *= r;
      }
    */
    return 0;
}

static long double lldet2(long double a, long double b,
                          long double c, long double d) {
    return  a*d - b*c;
}
static double orient2d(double a1, double a2, double b1, double b2, double c1, double c2) {
    long double a = a1 - c1,  b = b1 - c1;
    long double c = a2 - c2,  d = b2 - c2;
    return lldet2(a, b,  c, d);
}
static double det2i3(int v1, int v2, int v3) {
    double r;
    r = orient2d(mesh2.x[v1], mesh2.y[v1], mesh2.x[v2], mesh2.y[v2], mesh2.x[v3], mesh2.y[v3]);
    return  r;
}
static int idet2i3(int v1, int v2, int v3) {
    double d = det2i3(v2, v1, v3);
    if (d > 1e-16) return  +1;
    else if (d < -1e-16) return  -1;
    else return  0;
}


static int intsect(int a, int b, int c, int u, int v) {
    int uv, dup = 0;

    if ((u == a) || (u == b) || (u == c)) dup++;
    if ((v == a) || (v == b) || (v == c)) dup++;
    if (dup == 0) {
        uv = idet2i3(u, v, a) + idet2i3(u, v, b) + idet2i3(u, v, c);
        if ((uv == 3) || (uv == -3)) return 0;
        if (idet2i3(b, c, u) + idet2i3(b, c, v) == -2) return  0;
        if (idet2i3(c, a, u) + idet2i3(c, a, v) == -2) return  0;
        if (idet2i3(a, b, u) + idet2i3(a, b, v) == -2) return  0;
    } else if (dup == 1) {
        if (idet2i3(b, c, u) + idet2i3(b, c, v) == -1) return  0;
        if (idet2i3(c, a, u) + idet2i3(c, a, v) == -1) return  0;
        if (idet2i3(a, b, u) + idet2i3(a, b, v) == -1) return  0;
    } else return  0;
    return  1;
}
static int check(PStrucFace2 e, int pn) {
    int v1, v2, i, p1, p2;

    intedge = NULL;
    v1 = e->v1,  v2 = e->v2;
    if (idet2i3(v1, v2, pn) != 1) {
        //fprintf(stderr, "\ninverted? %d (%d %d %d)\n", idet2i3(v1, v2, pn), v1, v2, pn);
        return 1;
    }
    for (i = 0; i < tree2.nVicinityFace; i++) {
        p1 = tree2.vicinityFace[i]->v1,  p2 = tree2.vicinityFace[i]->v2;
        if (intsect(v1, v2, pn, p1, p2)) {
            intedge = tree2.vicinityFace[i];
            return  1;
        }
    }
    return  0;
}

static double height(PStrucFace2 e) {
    int v0, v1;
    double x, y, x0, y0, x1, y1, s, r, c;
    v0 = e->v1,  v1 = e->v2;
    x = mesh2.x[new_vert],  y = mesh2.y[new_vert];
    x0 = mesh2.x[v0],  y0 = mesh2.y[v0];
    x1 = mesh2.x[v1],  y1 = mesh2.y[v1];
    c = ((x-x0)*(x1-x0) + (y-y0)*(y1-y0)) / ((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    if (c<=0.0)  return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
    else if (c>=1.0)  return sqrt((x1-x)*(x1-x) + (y1-y)*(y1-y));
    s = (x0-x)*(y1-y) - (y0-y)*(x1-x);
    r = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
    return fabs(s)/r;
}

#define CND_MAX 1000
static int new_point(PStrucFace2 e) {
    double p, a, b, dist, r, q;
    double x, y, x1, y1, x2, y2;
    int i, j, k, m, pn;
    int v1, v2, nearVert;
    int ncnd, cnd[CND_MAX];
    double cnd_q[CND_MAX];
    int nchk, chk[CND_MAX];
    int neari;
    double radius, rmin, rv;
    double hc, h;
    int dirty;

    v1 = e->v1,  v2 = e->v2;
    x1 = mesh2.x[v1],  y1 = mesh2.y[v1],  x2 = mesh2.x[v2],  y2 = mesh2.y[v2];
    b = x1 - x2,  a = y2 - y1,  p = sqrt(a * a + b * b);
    if (p == 0.0)  return -1;
    a /= p,  b /= p;
    r = p/2.0;
    if (!boolFAF) {
	x = 0.5 * (x1 + x2) + a * 0.3 * p,  y = 0.5 * (y1 + y2) + b * 0.3 * p;
	p = sizeFace(x, y);
	if (p*p-r*r < p*p*3.0/4.0)  p = sqrt(r*r + p*p*3.0/4.0);
    } else {
	p = COARSEFACTOR * sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    }
    if (p<1.1*r)  p = 1.1*r;
    x = 0.5 * (x1 + x2) + a * sqrt(p*p-r*r),  y = 0.5 * (y1 + y2) + b * sqrt(p*p-r*r);
    addPoint(x, y);
    mesh2.nPoint--;
    new_vert = mesh2.nPoint;
    /*____________________ TEST __________ TEST ____________________________*/
    radius = distance(x, y, x1, y1);
    dist = nearest2(&nearVert, x, y, radius);
    r = radius * 1.0001220703125; 
    rmin = (beta*r + minrho + sqrt((beta*r - minrho)*(beta*r - minrho) + alpha))/2.0;
    r = radius * 2.0;

    vicinityFaces(x, y, r);
    
    cnd[0] = new_vert,  cnd_q[0] = func_q(v1, v2, new_vert),  ncnd = 1,  nchk = 0;
    neari = -1,  rv = rmin;
    hc = height(e),  dirty = 0;
    for (i = 0; i < tree2.nVicinityFace; i++) {
	h = height(tree2.vicinityFace[i]);
	if (h < 0.5*hc)  dirty++;//,  printf("dirty: %lf, %lf\n", hc, h);
	for (j = 0, pn = tree2.vicinityFace[i]->v1; j < 2; j++, pn = tree2.vicinityFace[i]->v2) {
	    if ((pn==v1) || (pn==v2))  continue;
	    if (idet2i3(v1, v2, pn) != 1)  continue;
	    r = distance(mesh2.x[pn], mesh2.y[pn], x, y);
	    if (r < rv)  rv = r,  neari = pn;
	    q = func_q(v1, v2, pn);
	    for (k=0; k<ncnd; k++) {
		if (pn == cnd[k])  break;
		if (q > cnd_q[k])  break;
	    }
	    if (k>=CND_MAX)  continue;
	    if (k==ncnd) {
		cnd[ncnd] = pn,  cnd_q[ncnd] = q,  ncnd++;
	    } else {
		if (pn == cnd[k])  continue;
		for (m=ncnd; m>k; m--)  cnd[m] = cnd[m-1],  cnd_q[m] = cnd_q[m-1];
		cnd[k] = pn,  cnd_q[k] = q,  ncnd++;
	    }
	}
    }
    if (dirty && ncnd > 1)  {
	chk[nchk++] = new_vert;
	new_vert = cnd[1];
    }
    if (neari >= 0)  new_vert = neari;

    for ( ; nchk < CND_MAX; nchk++) {
	chk[nchk] = new_vert;
	if (check(e, new_vert) == 0)  return 0;
	for (k=0; k < ncnd; k++) {
	    new_vert = cnd[k];
	    for (m=0; m <= nchk; m++) {
		if (new_vert == chk[m])  break;
	    }
	    if (m>nchk)  break;
	}
	if (k>=ncnd) {
	    return -2;
	}
    }
    return  1;
} /* new_point */ 

static int chknadd(int v1, int v2) {
    int i, p1, p2;
    
    for (i = 0; i < tree2.nVicinityFace; i++) {
	p1 = tree2.vicinityFace[i]->v1,  p2 = tree2.vicinityFace[i]->v2;
	if ((p1==v2) && (p2==v1)) {
	    todelete[ntodelete++] = tree2.vicinityFace[i];
	    return 1;
	}
    }
    addFace(v1, v2, 0);
    return 0;
}

static int newTria(int lab) {
    int nn;
    int v1, v2;
    PStrucFace2 e1;

    e1 = tree2.face[0];
    v1 = e1->v1,  v2 = e1->v2;
    nn = new_point(e1);
    if (nn != 0)  return nn;
    
    if (new_vert == mesh2.nPoint)  mesh2.nPoint++;

    todelete[0] = e1,  ntodelete = 1;
    addTria(v1, v2, new_vert, lab);
    chknadd(new_vert, v2);
    chknadd(v1, new_vert);
    for (nn = 0; nn < ntodelete; nn++)  remFace(todelete[nn]);
    return 0;
} /*newTria*/

static int fill_eadj(void);
static int fill_tadj(void);
static int opt_func(int nfixed);
static void savedump(void);
extern  char Huge * ppMemory;
extern  PStrucFace2  *ptree2face;
extern  int StopAfterinitRegion;
extern int    nVRTglobal;
extern int    nTRIglobal;


/* error codes:
 *  0 - success
 * -1 - zero sized edge (error in user data)
 * -2 - internal search failed (most likely self intersection in front)
 ****************************************************************************/
int makeTria(void) {
    int i = 0, j, err = 0, smooth = 0;

    init();
    initRegion();

    if ( StopAfterinitRegion != 0 ) {
        free(ppMemory);
        free(ptree2face);
        return 0;
    }

    for (i = 1; i <= mesh2.nRegion; i++) {
        mesh2.nRPoint[i-1] = mesh2.nPoint;
        mesh2.nRTria[i-1] = mesh2.nTria;
	minrho = dist(tree2.face[0]->v1, tree2.face[0]->v2);

	for (j=0; j<tree2.nFace; j++) {
	    if (minrho > dist(tree2.face[j]->v1, tree2.face[j]->v2))
		minrho = dist(tree2.face[j]->v1, tree2.face[j]->v2);
	}

	if (boolFAF)  minrho *= 0.95;
	else minrho = 0.0;
	beta = 0.5;
	minrho *= beta;
	alpha = 2.0*minrho*(1.0-beta)/beta;
	alpha *= alpha;
        while (tree2.nFace > 0) {
	    err = newTria(i);
	    if (err)  break;
            if (SHOWPROGRESS && (mesh2.nTria % 100 == 0))
		printf("\r Number  of  Point  =  %u    Number  of  Tria  =  %u    ", mesh2.nPoint,mesh2.nTria),  fflush(stdout);
            if (0)  savedump();
        }
	if (err)  break;
	eadj.ia   = (int*   )realloc(eadj.ia, sizeof(int   )*(mesh2.nPoint + 1));
	eadj.ja   = (int*   )realloc(eadj.ja, sizeof(int   )*(6*mesh2.nTria));
	tadj.ia   = (int*   )realloc(tadj.ia, sizeof(int   )*(mesh2.nPoint + 1));
	tadj.ja   = (edge*  )realloc(tadj.ja, sizeof(edge  )*(3*mesh2.nTria));
	if (SHOWPROGRESS)  printf("\n");
	fill_eadj();
	fill_tadj();
	if (0) savedump();
	smooth += opt_func(mesh2.nRPoint[i-1]);
	
	if (i != mesh2.nRegion) {
	    initAddRegion(i + 1);
	}
	mesh2.nRPoint[i-1] = mesh2.nPoint - mesh2.nRPoint[i-1];
        mesh2.nRTria[i-1] = mesh2.nTria - mesh2.nRTria[i-1];
    }


    /*	printf("\nRESULT :  %5u     %5u    \n",mesh2.nPoint,mesh2.nTria);*/
    if (!smooth)  regularity();

    //test_quality();

    outMesh();
    if (!err) {
	free(ppMemory);
	free(ptree2face);
	free(eadj.ia),  free(eadj.ja),  free(tadj.ia),  free(tadj.ja);
	eadj.ia = NULL,  eadj.ja = NULL,  tadj.ia = NULL,  tadj.ja = NULL;
    } else {
	nVRTglobal = mesh2.nPoint;
	nTRIglobal = mesh2.nTria;
    }

    return err;
} /*makeTria*/


typedef struct {
    int p;
    int n;
} plist;

static int add_glist(int a, int b, int *png, plist *glist, int *s) {
    int c = s[a];
    while (c>=0) {
	if (glist[c].p == b) return 0;
	c = glist[c].n;
    }
    glist[*png].p = b;
    glist[*png].n = s[a];
    s[a] = *png;
    (*png)++;
    return 1;
}

static int fill_eadj(void) {
    int *s;
    int ng;
    plist *glist;
    int i, j, n, l;
    int a, b, c;
    
    glist = (plist*)malloc(sizeof(plist)*6*mesh2.nTria);
    ng = 0;
    s = (int*)malloc(sizeof(int)*mesh2.nPoint);
    for (j=0; j<mesh2.nPoint; j++)  s[j] = -1;
    
    for (i=0; i<mesh2.nTria; i++) {
	a = mesh2.v1[i],  b = mesh2.v2[i],  c = mesh2.v3[i];
	add_glist(a, b, &ng, glist, s);
	add_glist(a, c, &ng, glist, s);
	add_glist(b, c, &ng, glist, s);
	add_glist(b, a, &ng, glist, s);
	add_glist(c, a, &ng, glist, s);
	add_glist(c, b, &ng, glist, s);
    }
    n = 0;
    for (j=0; j<mesh2.nPoint; j++) {
	eadj.ia[j] = n;
	l = s[j];
	while (l>=0) {
	    eadj.ja[n++] = glist[l].p;
	    l = glist[l].n;
	}
    }
    eadj.ia[mesh2.nPoint] = n;
    
    free(s),  free(glist);
    return n;
}


typedef struct {
    edge e;
    int n;
} elist;

static int add_hlist(int a, edge b, int *pnh, elist *hlist, int *s) {
    int c = s[a];
    while (c>=0) {
	if (hlist[c].e[0] == b[0] && hlist[c].e[1] == b[1]) return 0;
	c = hlist[c].n;
    }
    hlist[*pnh].e[0] = b[0];
    hlist[*pnh].e[1] = b[1];
    hlist[*pnh].n = s[a];
    s[a] = *pnh;
    (*pnh)++;
    return 1;
}
static int fill_tadj(void) {
    int *s;
    int nh;
    elist *hlist;
    int i, j, n, l;
    int a;
    edge b;
    
    hlist = (elist*)malloc(sizeof(elist)*3*mesh2.nTria);
    nh = 0;
    s = (int*)malloc(sizeof(int)*mesh2.nPoint);
    for (j=0; j<mesh2.nPoint; j++)  s[j] = -1;
    
    for (i=0; i<mesh2.nTria; i++) {
	a = mesh2.v1[i],  b[0] = mesh2.v2[i],  b[1] = mesh2.v3[i],  add_hlist(a, b, &nh, hlist, s);
	a = mesh2.v2[i],  b[0] = mesh2.v3[i],  b[1] = mesh2.v1[i],  add_hlist(a, b, &nh, hlist, s);
	a = mesh2.v3[i],  b[0] = mesh2.v1[i],  b[1] = mesh2.v2[i],  add_hlist(a, b, &nh, hlist, s);
    }
    n = 0;
    for (j=0; j<mesh2.nPoint; j++) {
	tadj.ia[j] = n;
	l = s[j];
	while (l>=0) {
	    tadj.ja[n][0] = hlist[l].e[0];
	    tadj.ja[n][1] = hlist[l].e[1];
	    n++;
	    l = hlist[l].n;
	}
    }
    tadj.ia[mesh2.nPoint] = n;
    
    free(s),  free(hlist);
    return n;
}

static int opt_func(int nfixed) {
    double d, delta, r, ds, x[2], z[2], rs, as;
    int j, l, s;
    int p, b, c;
    int pn, n_iters;
    int *ps;
    vertex *dx, *bkp;
    
    
    ps = (int*)malloc(sizeof(int)*(mesh2.nPoint - nfixed));
    pn = 0;
    ds = 0.0;
    for (j=nfixed; j<mesh2.nPoint; j++) {
	p = j;
	ps[pn++] = p;
	d = 0.0,  s = 0;
	for (l=eadj.ia[p]; l<eadj.ia[p+1]; l++) {
	    c = eadj.ja[l];
	    d += dist(p, c);
	    s++;
	}
	if (s > 0)  d /= s;
	ds += d;
    }
    if (pn == 0)  {
	free(ps);
	return 0;
    }
    
    ds /= pn;
    //printf("ds=%lf\n", ds);
    
    dx = (vertex*)malloc(sizeof(vertex)*pn);
    bkp = (vertex*)malloc(sizeof(vertex)*pn);
    for (j=0; j<pn; j++) {
	p = ps[j];
	bkp[j][0] = mesh2.x[p];
	bkp[j][1] = mesh2.y[p];
    }
    
    delta = 1.0;
    n_iters = 100; // 4000
    d = 0.25 * ds * ds / n_iters;  // 1.0
    while (1) {
	for (s=0; s<n_iters; s++) {
	    delta = 0.01 * ds * ds * (1.0 - 0.9*s/n_iters);
	    
	    //	delta = -as / n_iters;
	    rs = 0.0;
	    for (j=0; j<pn; j++) {
		p = ps[j];
		x[0] = 0.0,  x[1] = 0.0;
		for (l=tadj.ia[p]; l<tadj.ia[p+1]; l++) {
		    b = tadj.ja[l][0];
		    c = tadj.ja[l][1];
		    func_xy(p, b, c, z, delta, 8);
		    x[0] -= z[0],  x[1] -= z[1];
		}
		r = sqrt(x[0]*x[0] + x[1]*x[1]);
		if (r>1.0/ds) {
		    r = 1.0/ds*(1.0 + log(r*ds))/r;
		    x[0] *= r,  x[1] *= r;
		}
		dx[j][0] = x[0],  dx[j][1] = x[1];
		r = sqrt(x[0]*x[0] + x[1]*x[1]);
		if (rs < r)  rs = r;
	    }
	    for (j=0; j<pn; j++) {
		p = ps[j];
		mesh2.x[p] += d*dx[j][0];
		mesh2.y[p] += d*dx[j][1];
	    }
	    //printf("%3d: delta = %12.8lf, d = %12.8lf\n", s, delta, d*rs/ds);
	    //savedump();
	    //scanf("%d", &j);
	    //	if (d*rs/ds < 0.01)  break;
	    if (SHOWPROGRESS)  printf(" Smoothing %d%%\r", s+1),  fflush(stdout);
	}
	if (SHOWPROGRESS)  printf(" Smoothing done\n");
	break;
    }
    s = 0;
    for (j=0; j<mesh2.nTria; j++) {
	if (func_q(mesh2.v1[j], mesh2.v2[j], mesh2.v3[j]) <= 0.0)  s++;
    }
    if (s) {
	fprintf(stderr,"aniAFT: Quality improvement failed, falling back to simple smoothing\n");
	for (j=0; j<pn; j++) {
	    p = ps[j];
	    mesh2.x[p] = bkp[j][0];
	    mesh2.y[p] = bkp[j][1];
	}
    }
    
    free(ps);
    free(dx);
    free(bkp);

    return s;
}


static void savedump(void) {
    FILE *f;
    int i;
    
    f=fopen("dump.txt", "w");
    if (!f)  return;
    fprintf(f, "%d %d %d\n", mesh2.nPoint, mesh2.nTria, tree2.nFace);
    for (i=0; i<mesh2.nPoint; i++)  fprintf(f, "%lf %lf\n", mesh2.x[i], mesh2.y[i]);
    for (i=0; i<mesh2.nTria; i++)  fprintf(f, "%d %d %d\n", mesh2.v3[i]+1, mesh2.v2[i]+1, mesh2.v1[i]+1);
    for (i=0; i<tree2.nFace; i++)  fprintf(f, "%d %d\n", tree2.face[i]->v1+1, tree2.face[i]->v2+1);
    fclose(f);
}
