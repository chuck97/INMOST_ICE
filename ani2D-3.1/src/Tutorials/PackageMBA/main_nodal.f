C =====================================================================
      Program Main
C =====================================================================
      implicit none

c Maximum number of nodes, elements and boundary edges
      Integer   MaxP,  MaxE,  MaxF
      Parameter(MaxP = 20000, MaxF = 10000, MaxE = 2 * MaxP)

c Available memory 
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 2 000 000, MaxWi = 3 000 000)

C =====================================================================
C group (M)
      Integer  nP, nPfix, lbP(MaxP), fixP(MaxP)
      Real*8   XYP(2, MaxP)

      Integer  nF, nFfix, IPF(2, MaxF), lbF(MaxF), fixF(MaxF)

      Integer  nC, lbC(MaxF)
      Real*8   Crv(2, MaxF)
      EXTERNAL CrvFunction

      Integer  nE, nEfix, IPE(3, MaxE), lbE(MaxE), fixE(MaxE)

C group (CONTROL)
      Integer  control(6)
      Real*8   Metric_h(3, MaxP), Quality

C group (W)
      Real*8  rW(MaxWr)
      Integer iW(MaxWi)


C LOCAL VARIABLES
      Integer       i, nEStar, iERR
      Real*8        x, y

C =====================================================================
c ... load the initial mesh. The extension must be .ani
      Call loadMani(nP, nPfix, MaxP, XYP, lbP, fixP,
     &              nF, nFfix, MaxF, IPF, lbF, fixF,
     &              nC,              Crv, lbC, 
     &              nE, nEfix, MaxE, IPE, lbE, fixE,
     &              "../data/wing.ani")

c ... draw initial mesh. The name must have extension .ps
c     Call draw(nP, nF, nE, XYP, lbP, IPF, lbF, IPE, 'mesh_initial.ps')
      Call graph_demo(nP,XYP, nE,IPE, 'mesh_initial.ps',
     &               'Initial weakly anisotropic mesh')


c ... define a piecewise-linear metric 
c           20.5  -19.5  
c          -19.5   20.5
c
      Do i = 1, nP
         x = XYP(1, i)
         y = XYP(2, i)

         Metric_h(1, i) = 20.5
         Metric_h(2, i) = 20.5
         Metric_h(3, i) =-19.5
      End do


c ... generate adaptive mesh
      control(1) = 200     !  MaxSkipE
      control(2) = 15000   !  MaxQItr
      control(3) = 1       !  status
      control(4) = 1       !  flagAuto
      control(5) = 1       !  iPrint:   minimal level of output information
      control(6) = 0       !  iErrMesgt: only critical termination allowed

      Quality = 0.8D0      !  request shape-regular triangles in metric
      nEStar  = 2000       !  desired number of triangles

      Call mbaNodal(
c group (M)
     &      nP, nPfix, MaxP, XYP, lbP, fixP,
     &      nF, nFfix, MaxF, IPF, lbF, fixF,
     &      nC,              Crv, lbC, CrvFunction,
     &      nE, nEfix, MaxE, IPE, lbE, fixE,
c group (CONTROL)
     &      nEStar, Quality, control, Metric_h,
c group (W)
     &      MaxWr, MaxWi, rW, iW, iERR)


c ... draw the final mesh (we use colors saved in iW(1:nP)
c     the name must have extension .ps
c     Demo graphics has been activated
c     Call draw(nP, nF, nE, XYP, iW, IPF, lbF, IPE, 'mesh_final.ps')
      Call graph_demo(nP,XYP, nE,IPE, 'mesh_final.ps',
     &               'Final mesh quasi-unifom in analytic metric')


c ... save the final mesh. The extension must be .ani
      Call saveMani(nP, nPfix, XYP, lbP, fixP,
     &              nF, nFfix, IPF, lbF, fixF,
     &              nC,        Crv, lbC, 
     &              nE, nEfix, IPE, lbE, fixE,
     &              "save.ani")


c ... testing the result
      If(Quality.LT.0.3) Stop 911

      Stop 
      End



C =====================================================================
      Subroutine CrvFunction(tc, xyc, iFnc)
C =====================================================================
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
C         1. There are 2 different parametrizations with 
C            positive identificators (called function 
C            numbers) 2 and 5.
C
C =====================================================================
      Real*8  tc, xyc(2)

C group (Local variables)
      Real*8  CX2, tt, xx, yy, w, XB
      Real*8  X1, X2, XE, Y2, C1, C2, C3, C4, C5, wingAngle

C =====================================================================
c .. closed profile of the wing model 
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

c ... tail line after the wing
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






