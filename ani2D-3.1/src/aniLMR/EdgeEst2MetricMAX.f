C ======================================================================
      Subroutine EdgeEst2MetricMAX(Error, nP, nE, XYP, IPE, 
     &                             Metric, MaxWr, rW)
C ======================================================================
C  Routine uses method of shifts to calculate the nodal metric 
C  corresponding to edge-based error estimates in array Error.
C ======================================================================
C  Input:
C     Error(3, nE) - error estimates on edges
C
C     nP - the number of nodes
C     nE - the number of triangles
C
C     XYP(2, nP) - Cartesian coordinates of the nodes
C     IPE(3, nE) - the connecticity table
C
C  Output: 
C     Metric(3, nP) - nodal metric 
C
C  Work arrays: 
C     rW(MaxWr) - Real*8  of length MaxWr, MaxWr > nP
C ======================================================================
      implicit none

C Mesh
      Integer  nP, nE
      Real*8   XYP(2, *)
      Integer  IPE(3, *)

C Function and metric 
      Real*8   Error(3, *), Metric(3, *)

C Work arrays
      Integer  MaxWr
      Real*8   rW(*)

C Local variables
      Integer  i, j, n, iP1, iP2, iP3
      Real*8   det, H(3)

C ======================================================================
c ... memory check
      If(nP.GT.MaxWr) Call errMesLMR(1002, 'EdgeEst2MetricMAX', 
     &                     'Not enough memory for Real*8 arrays')

      Do n = 1, nP
         rW(n) = 0D0
      End do


c ... loop over elements
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call TriMetric(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), 
     &                  Error(1, n), H)

         det = H(1)*H(2) - H(3)**2


c  ...  update the maximum metric
         Do i = 1, 3
            iP1 = IPE(i, n)

            If(det.GT.rW(iP1)) Then
                rW(iP1) = det
                Do j = 1, 3
                   Metric(j, iP1) = H(j)
                End do
            End if
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine EdgeEst2MetricCell(Error, nE, XYP, IPE, 
     &                              Metric)
C ======================================================================
C  Input:
C     nE - the number of triangles
C
C     Error(3, nE) - given function on edges
C
C     XYP(2,  *) - Cartesian coordinates of the nodes
C     IPE(3, nE) - the connecticity table
C
C  Output: 
C     Metric(3, nE) - cell metric 
C
C ======================================================================
      implicit none

C Mesh
      Integer  nE, IPE(3, *)
      Real*8   XYP(2, *)

C Function and metric 
      Real*8   Error(3, *)
      Real*8   Metric(3, *)

C Local variables
      Integer  n, iP1, iP2, iP3

C ======================================================================

      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call TriMetric(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), 
     &                  Error(1, n), Metric(1, n))

      End do


      Return
      End






