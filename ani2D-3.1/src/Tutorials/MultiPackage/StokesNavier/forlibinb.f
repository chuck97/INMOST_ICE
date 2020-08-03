c ======================================================================
      Subroutine prevec(iprevec, dummy, x, y, iwork, dwork)
c ======================================================================
c This user routine evaluates the preconditioner-vector product for the 
c unknown jacobian. The precondtioner here is the factorized Stokes system.
c
c INPUT:
c    iprevec - control parameters for the routine (not used)
c    dummy   - different dummy allow to use different preconditioners
c    x       - input vector
c    iwork   - integer data arrays for the preconditioner
c    dwork   - real*8 data arrays for the preconditioner
c
c OUTPUT:
c    y       - preconditioner-vector product
c ======================================================================
      Integer    iprevec(*), dummy, iwork(*)
      Real*8     x(*),y(*), dwork(*), luinfo(90)

      Integer    sys
c ======================================================================
      sys = 0

      Call umf4sol(sys, y, x, iwork, dwork, luinfo)

      If(luinfo(1).LT.0) Then
         Write(*,*) 'Error occurred in umf4sol in prevec: ', luinfo(1)
         Stop
      End if

      Return
      End



c ======================================================================
      Subroutine fnlin(n, xcur, fcur, rpar, ipar, itrmf)
c ======================================================================
c This user routine evaluates the nonlinear residual for the FEM problem.
c
c  INPUT:
c     n       dimension of vectors
c     xcur    current vector
c     rpar    double precision user-supplied parameters
c     ipar    integer user-supplied parameters
c
c OUTPUT:
c     fcur    nonlinear residual vector (zero for the solution)
c     itrmf   flag for successful termination of the routine
c ======================================================================
      implicit none
      include 'assemble.fd'

      Integer  itrmf, n, ipar(*)
      Real*8   xcur(n), fcur(n), rpar(*)

      Integer  nv,nt,nb,MaxWi,na,nedges,nRow,nCol
      Integer  iux,iuy,ip,iuxe,iuye
      Integer  ipvrt,ipRHS,ipA,iplabelP,ipbnd,iplabelB,iptri
      Integer  iplabelT,ipIA,ipJA,ipiW, i,i1,i2,i3, MaxWr

      Real*8   rW(1), dotProd
      Integer  iDATA(1), Dbc, controlFEM(3)
      EXTERNAL Dbc,  FEM2Dext

c ======================================================================
c  Set local variables from user-supplied parameters.
      Call unpackMesh(nv, ipvrt, iplabelP, nb, ipbnd, iplabelB, 
     &                nt, iptri, iplabelT, ipar, rpar)

      na    = ipar(6) 
      MaxWi = ipar(7)

      nedges = (n - 3*nv) / 2  

c  structure of solution vector xcur
      iux = 1
      iuy = iux + nv
      ip  = iuy + nv
      iuxe= ip  + nv
      iuye= iuxe+ nedges

c  structure of rpar
      ipRHS = ipvrt + 2*nv
      ipA   = ipRHS + n

c  structure of ipar
      ipIA     = iplabelT + nt
      ipJA     = ipIA + n+1
      ipiW     = ipJA + na

      MaxWr = 1

c general sparse matrix in the AMG format 
      controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_AMG)
      controlFEM(2) = 0

      Call BilinearFormTemplate(
     &     nv, nb, nt, rpar(ipvrt), ipar(iplabelP), 
     &                 ipar(ipbnd), ipar(iplabelB),
     &                 ipar(iptri), ipar(iplabelT),
     &     FEM2Dext, xcur, iDATA, controlFEM,
     &     n+1,na,ipar(ipIA),ipar(ipJA),rpar(ipA),rpar(ipRHS),nRow,nCol,
     &     MaxWi, MaxWr, ipar(ipiW), rW)

c  Evaluate nonlinear residual.
      Call mulAcsr(nRow, ipar(ipIA),ipar(ipJA),rpar(ipA), xcur,fcur)
      Call daxpy(nRow,-1d0,rpar(ipRHS),1,fcur,1)

c  Set  termination flag for success.
      itrmf = 0

      Return
      End



