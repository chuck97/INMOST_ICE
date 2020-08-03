C ======================================================================
      Subroutine EdgeEst2GradMetricMAX(Error, nP, nE, XYP, IPE, 
     &                                 Metric, 
     &                                 MaxWr, rW)
C ======================================================================
C  Routines calculates a nodal metric using the method of shifts. 
C  The input data are error estimates prescribed to mesh edges.
C ======================================================================
C  Input:
C     nP  - the number of nodes
C     nE  - the number of triangles
C
C     Error(3, nE) - error estimates given on edges
C
C     XYP(2, nP) - Cartesian coordinates of the nodes
C     IPE(3, nE) - the connecticity table
C
C  Output: 
C     Metric(3, nE) - nodal tensor metric 
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
      Real*8   Error(3, *),  Metric(3, *)

C Work arrays
      Integer  MaxWr
      Real*8   rW(*)

C Local variables
      Integer  i, j, n, iP1, iP2, iP3
      Real*8   det, H(3), errL2loc, errH1loc

C ======================================================================
c ... memory check
      If(nP.gt.MaxWr) Call errMesLMR(1002, 'EdgeEst2GradMetricMAX', 
     &                     'Not enough memory for Real*8 arrays')

      Do n = 1, nP
         rW(n) = 0D0
      End do


c ... loop over elements
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call TriGradMetric(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), 
     &                      Error(1, n), H, errL2loc, errH1loc)

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
      Subroutine EdgeEst2GradMetricMAXVector(Error1, Error2,  
     &                                       nP, nE, XYP, IPE, 
     &                                       Metric, MaxWr, rW)
C ======================================================================
C  Routines calculates a nodal metric using the method of shifts. 
C  The input data are error estimates for two velocity components
C  prescribed to mesh edges.
C ======================================================================
C  Input:
C     nP  - the number of nodes
C     nE  - the number of triangles
C
C     Error1(3, nE) - error estimates for the first  component given on edges
C     Error2(3, nE) - error estimates for the second component given on edges
C
C     XYP(2, nP) - Cartesian coordinates of the nodes
C     IPE(3, nE) - the connecticity table
C
C  Output: 
C     Metric(3, nE) - nodal tensor metric 
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
      Real*8   Error1(3, *), Error2(3,*),  Metric(3, *)

C Work arrays
      Integer  MaxWr
      Real*8   rW(*)

C Local variables
      Integer  i, j, n, iP1, iP2, iP3
      Real*8   det, H(3)

C ======================================================================
c ... memory check
      If(nP.gt.MaxWr) Call errMesLMR(1002, 'EdgeEst2GradMetricMAX', 
     &                     'Not enough memory for Real*8 arrays')

      Do n = 1, nP
         rW(n) = 0D0
      End do

c ... loop over elements
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call TriGradMetricVector(XYP(1,iP1), XYP(1,iP2), XYP(1,iP3),
     &                            Error1(1, n), Error2(1, n), H)

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



C=======================================================================
      Subroutine EdgeEst2GradMetricCell(ErrorEst, nE, XYP, IPE, 
     &                                  Metric)
C=======================================================================
C  Routine calculates the cell-based metric using the edge-based error
C  estimates ErrorEst(:,1:nE). 
C=======================================================================
C  Input:
C     nE  - the number of triangles
C
C     ErrorEst(3, nE)  - given error estimates on edges
C
C     XYP(2,  *) - Cartesian coordinates of the nodes
C     IPE(3, nE) - the connectivity table
C
C  Output: 
C     Metric(3, nE) - cell metric 
C=======================================================================
      implicit none

C Mesh
      Integer  nE
      Real*8   XYP(2, *)
      Integer  IPE(3, *)

C Function and metric 
      Real*8   ErrorEst(3, *), Metric(3, *)

C Local variables
      Integer  i, n, iP1, iP2, iP3
      Real*8   H(3), errL2loc, errH1loc

C=======================================================================

c ... loop over elements
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call TriGradMetric(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), 
     &                      ErrorEst(1, n), H, errL2loc, errH1loc)

         Do i = 1, 3
            Metric(i, n) = H(i)
         End do

      End do

      Return
      End





