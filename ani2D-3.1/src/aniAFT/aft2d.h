void registeruserfn_( void * p );
void registersizefn_( void * p );

int aft2dfront_( 
             int *pnBr, int *br, int *pnVr, double *vrbr, 
             int *pnVRT, double *vrt, 
             int *pnTRI, int *tri, int *labtri, 
             int *pnBND, int *bnd, int *labbnd);

int aft2dboundary_( int *pnVert, double *bv,  
             int *pnLine, int *bl, double *bltail, double *hsze,
             int *pnVRT, double *vrt, 
             int *pnTRI, int *tri, int *labtri, 
             int *pnBND, int *bnd, int *labbnd,
             int *pnCRV, double *crv, int *iFNC );
