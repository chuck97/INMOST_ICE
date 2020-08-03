C ======================================================================
c The program demonstrates metric generation from edge or cell based
c data representing the local error. In the program, the interpolation 
C error is calculated using a given function Func. The error calculation 
C may be replaced by any a posteriori error estimates.
C ======================================================================


C ======================================================================
      Real*8 Function Func(xy)
C ======================================================================
      implicit none
      Real*8 xy(2), a

      a = 2.00D0
      Func = exp(2*(dabs(xy(1))**a + dabs(xy(2))**a))

      Return
      End



C ======================================================================
      Program Est2MetricRecovery
C ======================================================================
      implicit none
      integer nvmax,ntmax,nbmax
c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
      parameter(nvmax = 150 000, ntmax = 2*nvmax, nbmax = 10 000)

c working memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr =  1 000 000, MaxWi = 5 000 000)

C ======================================================================
C group (M)
      Integer  nv, nvfix, labelV(nvmax), fixedV(nvmax)
      Real*8   vrt(2, nvmax)

      Integer  nb, nbfix, bnd(2, nbmax), labelB(nbmax), fixedB(nbmax)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2, nbmax)

      Integer  nt, ntfix, tri(3, ntmax), labelT(ntmax), fixedT(ntmax)

C group (CONTROL)
      Real*8   Func, Metric(3, nvmax)
      EXTERNAL Func, ANI_CrvFunction

C group (W)
      Real*8   rW(MaxWr)
      Integer  iW(MaxWi)


C LOCAL VARIABLES
      Real*8   error(3*ntmax), U(nvmax)
      Integer  i,n, ipIRE,ipWork, MaxWiWork, nx
      

c === load the initial mesh
      Call loadMani(nv, nvfix, nvmax, vrt, labelV, fixedV,
     &              nb, nbfix, nbmax, bnd, labelB, fixedB,
     &              nc,               crv, labelC, 
     &              nt, ntfix, ntmax, tri, labelT, fixedT,
     &              "../data/simple.ani")

c ... refine uniformly the mesh 3 times. 
c ... ANI_CrvFunction, crv, labelC, and rW are dummy here. 
c ... iW(ipIRE) is a work array of size 3*nt; MaxWiWork >= 3*nt + nv
      Do i = 1, 3
         ipIRE = 1
         ipWork = ipIRE + 3 * nt
         MaxWiWork = MaxWi - 3 * nt

         Call uniformRefinement(
     &        nv, nvmax, nb, nbmax, nt, ntmax,
     &        vrt, bnd, labelB, tri, labelT, 
     &        ANI_CrvFunction, crv, labelC, iW(ipIRE), 
     &        rW, 1, iW(ipWork), MaxWiWork)            
      End do

      Write(*,'(A,I6,A)') 'The refined mesh containes ',nt,' triangles'


C ======================================================================
c Step 2. Compute Errors
C ======================================================================
         Do i = 1, nv
           U(i) = Func(vrt(1,i))
         End do


         nx = 10
c ... element-based error estimates
c        Call  triIntL8(U, nv, vrt, nt, tri, nx, error)

c ... edge-based error estimates 
         Call edgeIntL8(U, nv, vrt, nt, tri, nx, error)


C ======================================================================
c Step 3. Compute Metric
C ======================================================================
c ... generate nodal isotropic metric from element-based errors
c ... it requires calling of triIntL8()

c        Call CellEst2MetricZZ(nv, nt, vrt, tri,  
c    &                         error, Metric, 
c    &                         MaxWr, MaxWi, rW, iW)


c ... generate nodal tensor metric from edge-based errors by LeastSquares
c ... it requires calling of edgeIntL8()

         Call EdgeEst2MetricLS(nv, nt, vrt, tri, error, Metric, 
     &                         MaxWr, MaxWi, rW, iW)


c ... generate nodal tensor metric from edge-based errors through cell-based metric
c ... it requires calling of edgeIntL8()

c        Call EdgeEst2MetricMAX(error, nv, nt, vrt, tri, 
c    &                          Metric, MaxWr, rW)


c ... save the metric
      Open(10, file='metric')
        Do n = 1, nv
           Write(10, *) (Metric(i, n), i=1,3)
        End do
      Close(10)
      Write(*,'(A,/)') 'Recovered metric is stored in file bin/metric'


c ... testing the results
      Do n = 1, nv
         Do i = 1, 2
            If(Metric(i, n).LE.0D0) Stop 911
         End do
      End do

      Stop 
      End


