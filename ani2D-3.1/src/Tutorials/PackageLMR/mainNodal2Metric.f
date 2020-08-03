C ======================================================================
C The program demonstrates the local metric recovery from discrete
C solution defined at nodes of the mesh. The metric depends on what Lp 
c norm of the error the user wants to minimize.
C ======================================================================
      Program mainNodal2Metric
C ======================================================================
      implicit none

c Maximum number of nodes, elements and boundary edges
      Integer   MaxP,  MaxE,  MaxF
      Parameter(MaxP = 20000, MaxF = 10000, MaxE = 2 * MaxP)

c Available memory 
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 2 000 000)

C ======================================================================
C group (M)
      Integer  nP, nPfix, lbP(MaxP), fixP(MaxP)
      Real*8  XYP(2, MaxP)

      Integer  nF, nFfix, IPF(2, MaxF), lbF(MaxF), fixF(MaxF)

      Integer  nC, lbC(MaxF)
      Real*8   Crv(2, MaxF)

      Integer  nE, nEfix, IPE(3, MaxE), lbE(MaxE), fixE(MaxE)

C group (CONTROL)
      Real*8   Metric(3, MaxP), Sol(MaxP), Quality, Lp

C group (W)
      Real*8  rW(MaxWr)
      Integer iW(MaxWi)


C LOCAL VARIABLES
      Integer  n, i
      Real*8  x, y

C ==========================================================
c ... load the initial mesh
      Call loadMani(nP, nPfix, MaxP, XYP, lbP, fixP,
     &              nF, nFfix, MaxF, IPF, lbF, fixF,
     &              nC,              Crv, lbC, 
     &              nE, nEfix, MaxE, IPE, lbE, fixE,
     &      "../data/wing.ani")


c ... define a discrete solution at mesh points
      Do n = 1, nP
         x = XYP(1, n)
         y = XYP(2, n)

         Sol(n) = dexp( -(x-0.5)**2 - (y-0.5)**2) 
      End do


c ... generate metric (from SOL) optimal for the L_p norm
      Lp = 0             ! maximum norm
c     Lp = 1             ! L_1 norm

c choose the recovery method (ZZ or variational)
      Call Nodal2MetricVAR(Sol,
     &                     XYP, nP, IPE, nE, IPF, nF, Metric,
     &                     MaxWr, rW, MaxWi, iW)

      If(Lp.GT.0) Call Lp_norm(nP, Lp, Metric)


c ... save the metric
      Open(10, file='metric')
        Do n = 1, nP
           Write(10, *) (Metric(i, n), i=1,3)
        End do
      Close(10)
      write(*,'(A,/)') 'Recovered metric is stored in file bin/metric'


c ... testing the results
      Do n = 1, nP
         Do i = 1, 2
            If(Metric(i, n).LE.0D0) Stop 911
         End do
      End do

      Stop
      End


