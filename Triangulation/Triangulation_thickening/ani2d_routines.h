#pragma once

namespace ani
{
int aft2dboundary_(int *pnVert, double *bv,  
				  int *pnLine, int *bl, double *bltail, double *hsze,
				  int *pnVRT, double *vrt, 
				  int *pnTRI, int *tri, int *labtri, 
				  int *pnBND, int *bnd, int *labbnd,
				  int *pnCRV, double *crv, int *iFNC);
				  
void graph_(int *nv, double *vrt, int *nt, int *tri, char fName[100]);

void registersizefn_(void* p);
}
