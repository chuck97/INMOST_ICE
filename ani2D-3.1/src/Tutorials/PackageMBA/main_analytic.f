C =====================================================================
      Program Main
C =====================================================================
      implicit none

c Maximum number of nodes, elements and boundary edges
      Integer   nvmax, ntmax, nbmax
      Parameter(nvmax = 20000, nbmax = 10000, ntmax = 2 * nvmax)

c Available memory 
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 1 000 000)

C =====================================================================
C group (M)
      Integer  nv, nvfix, labelV(nvmax), fixedV(nvmax)
      Real*8   vrt(2, nvmax)

      Integer  nb, nbfix, bnd(2, nbmax), labelB(nbmax), fixedB(nbmax)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2, nbmax)
      EXTERNAL CrvFunction

      Integer  nt, ntfix, tri(3, ntmax), labelT(ntmax), fixedT(ntmax)

C group (CONTROL)
      Integer  control(6)
      Real*8   Quality
c The routine which defines the nodal metric
      Integer  MetricFunction
      External MetricFunction

C group (W)
      Real*8  rW(MaxWr)
      Integer iW(MaxWi)


C LOCAL VARIABLES
      Integer  nEStar, iERR

C =====================================================================
c === load the initial mesh. The extension must be .ani
      Call loadMani(nv, nvfix, nvmax, vrt, labelV, fixedV,
     &              nb, nbfix, nbmax, bnd, labelB, fixedB,
     &              nC,               crv, labelC, 
     &              nt, ntfix, ntmax, tri, labelT, fixedT,
     &             "../data/wing.ani")


c === draw the initial mesh (no labels of mesh points is available)
c     The name must have extensione with .ps
c     Demo graphics has been activated
c     Call draw(nv,nb,nt, vrt,labelV, bnd,labelB, tri,'mesh_initial.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_initial.ps', 
     &               'Initial weakly anisotropic mesh')


c === generate adaptive mesh
      control(1) = 500     !  MaxSkipE
      control(2) = 50000   !  MaxQItr
      control(3) = 1       !  status
      control(4) = 1       !  flagAuto
      control(5) = 1       !  iPrint:   average level of output information
      control(6) = 0       !  iErrMesgt: only critical termination allowed

      Quality = 0.8D0      !  request shape-regular triangles in metric
      nEStar  = 3000       !  desired number of triangles

      Call mbaAnalytic(
c group (M)
     &      nv, nvfix, nvmax, vrt, labelV, fixedV,
     &      nb, nbfix, nbmax, bnd, labelB, fixedB,
     &      nc,               crv, labelC, CrvFunction,
     &      nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &      nEStar, Quality, control, MetricFunction,
c group (W)
     &      MaxWr, MaxWi, rW, iW, iERR)


c === draw final mesh (mesh routine saved labels of points in iW(1:nv))
c the file name must have extension .ps
c demo graphics has been activated
c     Call draw(nv, nb, nt, vrt, iW, bnd, labelB, tri, 'mesh_final.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps', 
     &               'Final mesh quasi-uniform in isotropic metric')


c ... save the final mesh. The extension must be .ani
      Call saveMani(nv, nvfix, vrt, labelV, fixedV,
     &              nb, nbfix, bnd, labelB, fixedB,
     &              nc,        crv, labelC, 
     &              nt, ntfix, tri, labelT, fixedT,
     &             "save.ani")


c ... testing the results
      If(Quality.LT.0.65) Stop 911

      Stop
      End



C =====================================================================
      Integer Function MetricFunction(x, y, Metric)
C =====================================================================
C  This routine creates a metric at the given point (x,y). The
C  metric is a 2x2 positive definite symmetric tensor:
C                M11   M12
C      Metric =  
C                M12   M22
C
C  Only the upper triangular part of 2x2 array Metric may be defined.
C  
C  In this example, the metric is constant and isotropic.
C =====================================================================
      Real*8  a, x, y, Metric(2, 2)

      a = 0.8
      Metric(1,1) = max(1D0/((x-0.4)**2 + (y-0.5)**2 + 0.0001) ** a,
     &                  1D0/((x-0.6)**2 + (y-0.5)**2 + 0.0001) ** a)
      Metric(2,2) = Metric(1,1)
      Metric(1,2) = 0D0

      MetricFunction = 0

      Return
      End


C =====================================================================
      Subroutine CrvFunction(tc, xyc, iFnc)
C =====================================================================
C  The routine computes the Cartesian coordinates of point
C  xyc from its parametric coordinate tc.
C
C  tc     - [input]   the given parametric coordinate of point
C  xyc(2) - [output]  the Cartesian coordinate of the same point
C  iFnc   - [input]   the function number for computing
C
C  On input :  tc, iFnc
C  On output:  xyc(2)
C
C  *** Remarks:
C         1. In this example there are 2 different parametrizations 
C            identifed by iFnc = 2 and iFnc = 5.
C
C =====================================================================
      Real*8  tc, xyc(2)

c local variables
      Real*8  CX2, tt, xx, yy, w, XB
      Real*8  X1, X2, XE, Y2, C1, C2, C3, C4, C5, wingAngle

C =====================================================================
c .. parametrization of the wing
      If(iFnc.EQ.2) Then
         C1 = 0.17735D0
         C2 =-0.075597D0
         C3 =-0.212836D0
         C4 = 0.17363D0
         C5 =-0.062547D0

         wingAngle = 0D0
         
         X1 = 0.45D0
         X2 = 0.55D0
         XE = 1D0
         Y2 = 5D-1

         X1 = 0.40D0
         X2 = 0.60D0

         CX2 = 5D-1 * (X1 + X2)

         tt = 2 * tc
         if(tt.GT.1D0) Then
            tt = tt - 1D0
            isgn = -1
         Else
            tt = 1D0 - tt
            isgn = 1
         End if

         xx = X1 + tt * (X2 - X1)
         yy = Y2 + isgn * (X2 - X1) * (C1 * dsqrt(tt) + C2 * tt +
     &        C3 * tt * tt + C4 * tt * tt * tt +
     &        C5 * tt * tt * tt * tt)

         xyc(1) = CX2 + (xx - CX2) * dcos(wingAngle) +
     &                  (yy - Y2)  * dsin(wingAngle)

         xyc(2) = Y2 - (xx - CX2) * dsin(wingAngle) +
     &                 (yy - Y2)  * dcos(wingAngle)

c ... parametrization of a straight tail line behind the wing
      Else If(iFnc.EQ.5) Then
         wingAngle = 0D0

         X1 = 0.45D0
         X2 = 0.55D0
         XE = 1D0
         Y2 = 5D-1

         X1 = 0.40D0
         X2 = 0.60D0

         w = 5D-1 * (X2 - X1)
         XB = X2 - w * (1D0 - dcos(wingAngle))

         xyc(1) = XB + (XE - XB) * tc
         xyc(2) = Y2 + w * dsin(-wingAngle)
      Else
         Write(*, '(A,I5)') 'Undefined function =', iFnc
         Stop
      End if

      Return
      End






