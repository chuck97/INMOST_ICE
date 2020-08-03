C ================================================================
c Routine generates isotropic metric from element-based error 
c estimates. The nodal metric is recovered with the ZZ 
c interpolation algorithm. 
C ================================================================
      Subroutine CellEst2MetricZZ(nP, nE, XYP, IPE, 
     &                            Error, Metric,       
     &                            MaxWr, MaxWi, rW, iW) 
C ================================================================
c Input:
      Integer nP, nE       ! numbers of nodes and elements
      Real*8  XYP(2, *)    ! coordinates of mesh nodes
      Integer IPE(3, *)    ! connectivity table for elements
      Real*8  Error(*)     ! element-based error estimates

c Output:
      Real*8  Metric(3, *) ! node-based metric

c Working arrays:
      Integer MaxWr, MaxWi
      Integer iW(MaxWi)
      Real*8  rW(MaxWr)

C ================================================================
c ... memory test
      iEnd = 3 * nP + 3 * nE
      If(iEnd.GT.MaxWi) Call errMesLMR(1001, 'CellEst2MetricZZ', 
     &                      'increase size of work memory MaxWi')

c ... ZZ-interpolation of element-based metric into nodes
c     the allocated dmemory checked above is ssufficient
      Call P02P1(nP, nE, XYP, IPE, Error, rW, MaxWi, iW, iERR)

c ... create the metric
      Do i = 1, nP
         Metric(1, i) = rW(i)
         Metric(2, i) = Metric(1, i)
         Metric(3, i) = 0
      End do

      Return
      End

