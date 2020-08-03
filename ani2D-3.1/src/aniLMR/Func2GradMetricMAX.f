C=======================================================================
      Subroutine Func2GradMetricMAX(Func, nP, nE, XYP, IPE, 
     &                              Metric, 
     &                              MaxWr, rW)
C=======================================================================
C  Input:
C     Func       - Real*8 Function f(xy), where xy(2)
C
C     nP         - the number of nodes
C     nE         - the number of triangles
C
C     XYP(2, nP) - Cartesian coordinates of the nodes
C     IPE(3, nE) - the connecticity table
C
C  Output: 
C     Metric(3, nE) - nodal tensor metric 
C
C  Work arrays: 
C     rW(MaxWr) - Real*8  of length MaxWr
C=======================================================================
      implicit none

C Mesh
      Integer  nP, nE
      Real*8   XYP(2, *)
      Integer  IPE(3, *)

C Function and metric 
      Real*8   Func, Metric(3, *)
      EXTERNAL Func

C Work arrays
      Integer  MaxWr
      Real*8   rW(*)

C Local variables
      Integer  i, j, n, i1, i2, iP1, iP2, iP3, iref(4)
      Real*8   det, Up1, Up2, Ur1, xyt(2), H(3), Error(3)
      Real*8   errL2loc, errH1loc

      DATA     iref/1, 2, 3, 1/

C=======================================================================
c ... memory check
      If(nP.gt.MaxWr) Call errMesLMR(1002, 'Func2GradMetricMAX', 
     &                     'Not enough memory for Real*8 arrays')

      Do n = 1, nP
         rW(n) = 0D0
      End do


c ... loop over elements
      Do n = 1, nE
         Do i1 = 1, 3
            i2 = iref(i1 + 1)
 
            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            Up1 = Func(XYP(1, iP1))
            Up2 = Func(XYP(1, iP2))

            Do i = 1, 2
               xyt(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do

            Ur1 = Func(xyt)

            Error(i1) = dabs(4 * (Ur1 - (Up1 + Up2) / 2))
         End do

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call TriGradMetric(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), 
     &                      Error, H, errL2loc, errH1loc)

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
      Subroutine Func2GradMetricCell(Func, nE, XYP, IPE, 
     &                               Metric)
C=======================================================================
C  Input:
C     Func       - Real*8 Function f(xy), where xy(2)
C
C     nE         - the number of triangles
C
C     XYP(2, nP) - Cartesian coordinates of the nodes
C     IPE(3, nE) - the connecticity table
C
C  Output: 
C     Metric(3, nE) - cell metric 
C=======================================================================
      implicit none

C Mesh
      Integer  nE, IPE(3, *)
      Real*8   XYP(2, *)

C Function and metric 
      Real*8   Func, Metric(3, *)
      EXTERNAL Func

C Local variables
      Integer  i, n, i1, i2, iP1, iP2, iP3, iref(4)
      Real*8   Up1, Up2, Ur1, xyt(2), H(3), Error(3)
      Real*8   errL2loc, errH1loc

      DATA     iref/1, 2, 3, 1/

C=======================================================================

c ... loop over elements
      Do n = 1, nE
         Do i1 = 1, 3
            i2 = iref(i1 + 1)
 
            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            Up1 = Func(XYP(1, iP1))
            Up2 = Func(XYP(1, iP2))

            Do i = 1, 2
               xyt(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do

            Ur1 = Func(xyt)

            Error(i1) = dabs( 4 * (Ur1 - (Up1 + Up2) / 2) )
         End do

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call TriGradMetric(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3),
     &                      Error, H, errL2loc, errH1loc)

         Do i = 1, 3
            Metric(i, n) = H(i)
         End do

      End do


      Return
      End


