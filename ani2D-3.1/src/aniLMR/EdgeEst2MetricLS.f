C ======================================================================
C Routine generates nodal tensor metric from edge-based error estimates
C using the Least Square solution of the following system
C
C    (M(a_i) e_k,e_k) = \eta_k
C
C for all edges k incident to a mesh node a_i.
C ======================================================================
      Subroutine EdgeEst2MetricLS(nP, nE, XYP, IPE,
     &                            Error, Metric,       
     &                            MaxWr, MaxWi, rW, iW) 
C ======================================================================
      include 'makS.fd'
      include 'magic.fd'
C ======================================================================
c Input:
      Integer nP, nE       ! number of nodes and elements
      Real*8  Error(3, *)  ! edge error estimates data
      Real*8  XYP(2, *)    ! coordinates of mesh nodes
      Integer IPE(3, *)    ! connectivity table

c Output:
      Real*8  Metric(3, *) ! node-based metric 

c Working arrays:
      Integer MaxWr, MaxWi
      Integer iW(MaxWi)
      Real*8  rW(MaxWr)

c Local variables
      Real*8  ebuf(MaxS), vbuf(2, MaxS), det

C ======================================================================
c ... form list edge -> elements
      iIRE = 1
      inEP = iIRE + 3 * nE
      iIEP = inEP + nP
      iEnd = iIEP + 3 * nE
      If(iEnd.GT.MaxWi) Call errMesLMR(1001, 
     &                      'EdgeEst2MetricLS', 'Please increase MaxWi')

      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))


c ... collect errors on edges
      Do i = 1, nE
         Do j = 1, 3
            iR = iW(iIRE + (i-1)*3 + j-1)
            rW(iR) = error(j, i)
         End do
      End do


c ... form list edges -> points
      MaxR = nR
      iIPR = 1
      inEP = iIPR + 2 * MaxR
      iIEP = inEP + nP
      iEnd = iIEP + 3 * nE
      If(iEnd.GT.MaxWi) Call errMesLMR(1001, 
     &                      'EdgeEst2MetricLS', 'Please increase MaxWi')

      Call listR2P(nP, nR, nE, MaxR, IPE, iW(iIPR), iW(inEP), iW(iIEP))


c ... form list points -> edges 
      inRP = inEP
      iIRP = inRP + nP
      iEnd = iIRP + 2 * nR
      If(iEnd.GT.MaxWi) Call errMesLMR(1001, 
     &                      'EdgeEst2MetricLS', 'Please increase MaxWi')

      Call backReferences(nP, nR, 2,2, iW(iIPR), iW(inRP), iW(iIRP))


c ... compute nodal values of metric
      i2 = 0
      Do i = 1, nP
         i1 = i2 + 1
         i2 = iW(inRP + i - 1)

         m = 0
         Do j = i1, i2
            iR  = iW(iIRP + j - 1)

            iP1 = iW(iIPR + (iR-1) * 2)
            iP2 = iW(iIPR + (iR-1) * 2 + 1)

            m = m + 1
            ebuf(m)    = rW(iR)
            vbuf(1, m) = XYP(1, iP1) - XYP(1, iP2)
            vbuf(2, m) = XYP(2, iP1) - XYP(2, iP2)
         End do    

         Call NodalMetric(m, ebuf, vbuf, Metric(1, i))

c ... take spectral modulo of metric, |M|
         Call SpectralModule(Metric(1, i), det)
      End do

      Return
      End
      


C ================================================================
      Subroutine NodalMetric(k, values, xy, metric)
C ================================================================
C  This routine uses least square linear solution to the system
C  (Metric xy_i, xy_i) = values_i, i=1,\dots,k
C ================================================================
      implicit none

      Integer  k
      Real*8   xy(2,*), values(*), metric(3)

c (local variables)
      Real*8   A(3, 3), S(3), work(30)
      Integer  i, j, ipiv(3), info

C ================================================================
      Do i = 1, 3
         Do j = 1, 3
            A(i, j) = 0D0
         End do
         S(i) = 0D0
      End do

c ... generate the least squares matrix
      Do i = 1, k
         A(1,1) = A(1,1) + xy(1, i)**4 
         A(1,2) = A(1,2) + 2 * xy(1, i)**3 * xy(2, i)
         A(1,3) = A(1,3) + xy(1, i)**2 * xy(2, i)**2

         A(2,2) = A(2,2) + 4 * xy(1, i)**2 * xy(2, i)**2
         A(2,3) = A(2,3) + 2 * xy(1, i) * xy(2, i)**3
         A(3,3) = A(3,3) + xy(2, i)**4 
      End do

c ... generate the RHS
      Do i = 1, k
         S(1) = S(1) + xy(1, i)**2 * values(i)
         S(2) = S(2) + 2 * xy(1, i) * xy(2, i) * values(i)
         S(3) = S(3) + xy(2, i)**2 * values(i)
      End do

c ... fix the matrix (a bug???)
      A(2,1) = A(1,2)
      A(3,1) = A(1,3)
      A(3,2) = A(2,3)
      Do i = 1, 3
         A(i, i) = A(i, i) + 1e-16
      End do

      Call dsysv('U', 3, 1, A, 3, ipiv, S, 3, work, 30, info)

      If(info.NE.0) Call errMesLMR(1001, 
     &                  'NodalMetric', 'Lapack routine dsysv')

c ... write the metric
      metric(1) = S(1) ! m_xx
      metric(3) = S(2) ! m_xy
      metric(2) = S(3) ! m_yy

      Return
      End

