C ======================================================================
      Program Main
C ======================================================================
c Maximum number of nodes, elements and boundary edges
      Integer   MaxP,  MaxE,  MaxF
      Parameter(MaxP = 20000, MaxF = 10000, MaxE = 2 * MaxP)

c Available memory 
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 2 000 000)


C ======================================================================
C group (M)
      Integer  nP, nPfix, lbP(MaxP), fixP(MaxP)
      Real*8   XYP(2, MaxP)

      Integer  nF, nFfix, IPF(2, MaxF), lbF(MaxF), fixF(MaxF)

      Integer  nC, lbC(MaxF)
      Real*8   Crv(2, MaxF)
      EXTERNAL ANI_CrvFunction

      Integer  nE, nEfix, IPE(3, MaxE), lbE(MaxE), fixE(MaxE)

C group (CONTROL)
      Integer  control(6)
      Real*8   Quality
c The routine which defines the nodal metric
      Integer  ANI_MetricFunction
      External ANI_MetricFunction

C group (W)
      Real*8  rW(MaxWr)
      Integer iW(MaxWi)


C LOCAL VARIABLES
      Integer i, nEStar, iERR

C ======================================================================
c ... load the initial mesh. The extension must be .ani
      Call loadMani(nP, nPfix, MaxP, XYP, lbP, fixP,
     &              nF, nFfix, MaxF, IPF, lbF, fixF,
     &              nC,              Crv, lbC, 
     &              nE, nEfix, MaxE, IPE, lbE, fixE,
     &              "../data/texas.ani")


c ... shift  node 53 in order to get bad triangles
      XYP(1,53) = XYP(1,53) - 1d0


c ... draw the initial mesh (no labels of mesh points is available)
c     the name must have extension .ps
c     Call draw(nP, nF, nE, XYP, lbP, IPF, lbF, IPE, 'mesh_initial.ps')
      Call graph_demo(nP,XYP, nE,IPE, 'mesh_initial.ps',
     &               'Initial distorted mesh')


c ... generate adaptive mesh
      control(1) = 300     !  MaxSkipE
      control(2) = 15000   !  MaxQItr
      control(3) = 32      !  status
      control(4) = 1       !  flagAuto
      control(5) = 1       !  iPrint:   minimal level of output information
      control(6) = 0       !  iErrMesgt: only critical termination allowed

      Quality = 0.8D0      !  request shape-regular triangles in metric
      nEStar  = 400        !  desired number of triangles


      Call mbaFixShape(
c group (M)
     &     nP, nPfix, MaxP, XYP, lbP, fixP,
     &     nF, nFfix, MaxF, IPF, lbF, fixF,
     &     nC,              Crv, lbC, ANI_CrvFunction,
     &     nE, nEfix, MaxE, IPE, lbE, fixE,
c group (CONTROL)
     &     nEStar, Quality, control, ANI_MetricFunction,
c group (W)
     &     MaxWr, MaxWi, rW, iW, iERR)


c ... draw final mesh (aniMBA saved labels of mesh points in iW(1:nP))
c     the name must have extension .ps
c     demo graphics has been activated
c     Call draw(nP, nF, nE, XYP, iW, IPF, lbF, IPE, 'mesh_final.ps')
      Call graph_demo(nP,XYP, nE,IPE, 'mesh_final.ps',
     &               'Final improved mesh')


c ... save the final mesh. The extension must be .ani
      Call saveMani(nP, nPfix, XYP, lbP, fixP,
     &              nF, nFfix, IPF, lbF, fixF,
     &              nC,        Crv, lbC, 
     &              nE, nEfix, IPE, lbE, fixE,
     &     "save.ani")


c ... testing the results
      If(Quality.LT.0.7) Stop 911

      Stop
      End




