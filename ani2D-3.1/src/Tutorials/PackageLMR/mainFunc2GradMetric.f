C ======================================================================
C The program demonstrates building the optimal metric for the gradient
C of the P1 interpolation error. The metric depends on Lp norm of error
C the user wants to minimize. The interpolation error is calculated
c using a given function Func. The interpolation error may be replaced
C with any error estimates.
C ======================================================================


C ======================================================================
      Real*8 Function Func(xy)
C ======================================================================
      implicit none
      Real*8  xy(2)

      Func = xy(1)**2 * xy(2)
     &     + xy(2)**3 + dtanh(10*(dsin(5*xy(2)) - 2*xy(1)))

      Return
      End



C ======================================================================
      Subroutine FuncGrd(xy, fGrd)
C ======================================================================
      implicit none
      Real*8   xy(2), fGrd(2), u

      u = (dcosh(10*(dsin(5*xy(2)) - 2*xy(1))))**2
      fGrd(1) = 2*xy(1)*xy(2) - 20d0/u
      fGrd(2) = xy(1)**2 + 3*xy(2)**2 + 50*dcos(5*xy(2)) / u

      Return
      End



C ======================================================================
      Program mainFunc2Metric
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
      Real*8  XYP(2, MaxP)

      Integer  nF, nFfix, IPF(2, MaxF), lbF(MaxF), fixF(MaxF)

      Integer  nC, lbC(MaxF)
      Real*8   Crv(2, MaxF)

      Integer  nE, nEfix, IPE(3, MaxE), lbE(MaxE), fixE(MaxE)

C group (CONTROL)
      Real*8   Metric(3, MaxP), Lp

C group (W)
      Real*8  rW(MaxWr)
      Integer  iW(MaxWi)


C LOCAL VARIABLES
      Real*8   Func
      EXTERNAL Func

C ======================================================================
c ... load the initial mesh
      Call loadMani(nP, nPfix, MaxP, XYP, lbP, fixP,
     &              nF, nFfix, MaxF, IPF, lbF, fixF,
     &              nC,              Crv, lbC, 
     &              nE, nEfix, MaxE, IPE, lbE, fixE,
     &              "../data/square.ani")


c ... generate metric (from SOL) optimal for the L_p norm
      Lp = 0             ! maximum norm
c     Lp = 1             ! L_1 norm

c the recovery of continuous metric based on methods of shifts for
c minimizing gradient of the interpolation error
      Call Func2GradMetricMAX(Func, nP, nE, XYP, IPE, 
     &                        Metric, MaxWr, rW)

      If(Lp.GT.0) Call Lp_gradnorm(nP, Lp, Metric)


c ... save the metric
      Open(10, file='metric')
        Do n = 1, nP
           Write(10, *) (Metric(i, n), i=1,3)
        End do
      Close(10)
      Write(*,'(A,/)') 'Recovered metric is stored in file bin/metric'


c ... testing the results
      Do n = 1, nP
         Do i = 1, 2
            If(Metric(i, n).LE.0D0) Stop 911
         End do
      End do

      Stop
      End


