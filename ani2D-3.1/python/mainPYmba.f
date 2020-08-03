C ======================================================================
      Program Main
C ======================================================================
      implicit none

c Maximum number of nodes, elements and boundary edges
      Integer   nvmax,  ntmax,  nbmax
      Parameter(nvmax = 20000, nbmax = 10000, ntmax = 2 * nvmax)

c Available memory 
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 2 000 000)

C =====================================================================
C group (M)
      Integer  nv, nvfix, labelV(nvmax), fixedV(nvmax)
      Real*8   vrt(2, nvmax)

      Integer  nb, nbfix, bnd(2, nbmax), labelB(nbmax), fixedB(nbmax)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2, nbmax)
      EXTERNAL CrvFunction_disk

      Integer  nt, ntfix, tri(3, ntmax), labelT(ntmax), fixedT(ntmax)

C group (CONTROL)
      Integer  control(6)
      Real*8   Quality
c The routine which defines the nodal metric
      Integer  MetricFunction_user
      External MetricFunction_user

C group (W)
      Real*8  rW(MaxWr)
      Integer iW(MaxWi)


C LOCAL VARIABLES
      Integer  i, nEStar, iERR

      Real*8   varF1, varG1, varH1
      COMMON  /metricPY/ varF1, varG1, varH1

C ======================================================================
c ... load the initial mesh. The extension must be .ani
      Call loadMani(nv, nvfix, nvmax, vrt, labelV, fixedV,
     &              nb, nbfix, nbmax, bnd, labelB, fixedB,
     &              nC,              Crv, labelC, 
     &              nt, ntfix, ntmax, tri, labelT, fixedT,
     &      "../data/python.ani")

c ... load the parameters
      Open(10, file="aniPY_mba.txt")
      Read(10, *) varF1
      Read(10, *) varG1
      Read(10, *) varH1
      Close(10)


c === draw the initial mesh
c the name must terminate with .ps
      Do i = 1, nv
         iW(i) = 0
      End do
      Call draw(nv,nb,nt, vrt,labelV, bnd,labelB, tri,'mesh_initial.ps')


c === generate adaptive mesh
      control(1) = 500     !  MaxSkipE
      control(2) = 25000   !  MaxQItr
      control(3) = 1       !  status
      control(4) = 1       !  flagAuto
      control(5) = 2       !  iPrint:   average level of output information
      control(6) = 0       !  iErrMesgt: only critical termination allowed

      Quality = 0.8D0      !  request shape-regular triangles in metric
      nEStar  = 1000       !  desired number of triangles

      Call mbaAnalytic(
c group (M)
     &      nv, nvfix, nvmax, vrt, labelV, fixedV,
     &      nb, nbfix, nbmax, bnd, labelB, fixedB,
     &      nc,               crv, labelC, CrvFunction_disk,
     &      nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &      nEStar, Quality, control, MetricFunction_user,
c group (W)
     &      MaxWr, MaxWi, rW, iW, iERR)


c === draw the final mesh 
c     The name must terminate with .ps
      Do i = 1, nv
         iW(i) = 0
      End do
      Call graph(nv,vrt, nt,tri,'mesh_final.ps')
      Call graph_py(nv,vrt, nt,tri,'mesh_final.py')

      Stop
      End

 

C ======================================================================
      Integer Function MetricFunction_user(x, y, Metric)
C ======================================================================
C  This routine creates a metric at the given point (x,y). The
C  metric is a 2x2 positive definite symmetric tensor:
C                M11   M12
C      Metric =  
C                M12   M22
C
C  Only the upper triangular part of array Metric must be defined.
C ======================================================================
      implicit none
      Real*8  x, y, Metric(2, 2)

      Real*8  varF1, varG1, varH1
      COMMON /metricPY/varF1, varG1, varH1

      Metric(1,1) = varF1
      Metric(2,2) = varG1
      Metric(1,2) = varH1

      MetricFunction_user = 0

      Return
      End



C ======================================================================
      Subroutine CrvFunction_disk(tc, xyc, iFnc)
C ======================================================================
C  The routine computes the Cartesian coordinates of point
C  xyc from its parametric coordinate tc.
C
C  tc     - the given parametric coordinate of point
C  xyc(2) - the Cartesian coordinate of the same point
C  iFnc   - the function number for computing
C
C  On input :  tc, iFnc
C  On output:  xyc(2)
C
C  *** Remarks:
C         1. There is 1 parametrization of a circle.
C
C ======================================================================
      Real*8  tc, xyc(2)
      Real*8  PI

      If(iFnc.EQ.1) Then
         PI=4*datan(1d0)
         xyc(1) = dcos( PI*tc )
         xyc(2) = dsin( PI*tc )
      Else
         Write(*, '(A,I5)') 'Undefined function =', iFnc
         Stop
      End if

      Return
      End






