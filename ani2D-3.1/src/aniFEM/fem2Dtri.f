C ======================================================================
      Subroutine fem2Dtri(XY1, XY2, XY3, 
     &                    operatorA, FEMtypeA, operatorB, FEMtypeB, 
     &                    label, D, dDATA, iDATA, iSYS, order,  
     &                    LDA, A, nRow, nCol)
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'
C ======================================================================
C  *** ITRODUCTION ***
C  The routine computes elemental matrix for a bilinear form 
C
C  (1)            <D OpA(u), OpB(v)>               
C
C  where D is a tensor, OpA and OpB are linear operators, and u and v 
C  are finite element functions, u in "A", and v in "B". Note that 
C  solution of a non-linear problem may involve a Newton-like iterative 
C  loop. In this case D may depend on a discrete function (e.g. 
C  approximation from the previous iterative step). If this is the case, 
C  evaluation of D maybe quite complex and require additional data. At 
C  the moment only one plug for additional data has been reserved for 
C  the user (see formula (3)). 
C 
C  In order to compute the right hand side, we can use the following
C  trick:
C
C  (2)             f(v) = < D v, FEM_P0 > 
C
C  where action of D is given by function f, read also comments about
C  formula (3) below.
C
C   
C ======================================================================
C  *** BILINER FORMS ***
C ======================================================================
C  The constant basis function is E0.
C  The vertex-based unknowns are L1, L2, L3 (follow the order 
C  of vertices). The edge-based unknwons are B1, B2, B3 (correspond 
C  to edges 12, 23 and 31). The element-based unknown is E1.
C  The possible choices for these spaces are:
C 
C    FEM_P0 - piecewise constant. 
C             Order of unknowns: E1
C    FEM_P1 - continuous piecewise linear. 
C             Order of unknowns: L1,L2,L3  
C    FEM_P2 - continuous piecewise quadratic. 
C             Order of unknowns: L1,L2,L3, B1,B2,B3
C    FEM_P3 - continuous piecewise cubic. 
C             Order of unknowns: L1,L2,L3, B1,B2,B3, B1,B2,B3, E1
C    FEM_P4 - continuous piecewise quartic. 
C             Order of unknowns: L1,L2,L3, B1,B2,B3, B1,B2,B3, B1,B2,B3, E1,E2,E3
C    FEM_P1vector - vector continuous piecewise linear. 
C             Order of unknowns: L1,L2,L3 (x-component), L1,L2,L3 (y-component) 
C    FEM_P2vector  - vector continuous piecewise quadratic. 
C             Order of unknowns: L1,L2,L3, B1,B2,B3 (x-component),
C                                L1,L2,L3, B1,B2,B3 (y-component) 
C    FEM_P2reduced - Bernardi-Fortin-Raugel finite element: vector 
C             continuous piecewise linear functions enriched by edge bubbles. 
C             Order of unknowns: L1,L2,L3 (x-component), B1,B2,B3 (bubbles),
C                                L1,L2,L3 (y-component) 
C    FEM_MINI - vector continuous piecewise linear functions enriched 
C             by a central bubble. 
C             Order of unknowns: L1,L2,L3 (x-component), E1 (bubble),
C                                L1,L2,L3 (y-component), E0 (bubble) 
C    FEM_RT0 - lowest order Raviart-Thomas finite elements
C             Order of unknowns: B1,B2,B3
C    FEM_BDM1 - lowest order Brezzi-Douglas-Marini finite elements
C             Order of unknowns: B1,B2,B3 (closer to vertices P1,P2,P3), 
C                                B1,B2,B3 (closer to vertices P2,P3,P1)
C    FEM_CR1 - Crouzeix-Raviart finite elements
C             Order of unknowns: B1,B2,B3 
C    FEM_CR1vector - Crouzeix-Raviart vector finite elements
C             Order of unknowns: B1,B2,B3 (x-component), 
C                                B1,B2,B3 (y-component) 
C
C  The available operators are:
C    IDEN  - identity operator 
C    GRAD  - gradient operator
C    DIV   - divergence operator
C    CURL  - rotor operator
C    DUDX  - partial derivative d/dx
C    DUDY  - partial derivative d/dy
C
C  A quadrature formula is chosen as follows:
C    order = 1  -  quadrature formula with one center point
C    order = 2  -  quadrature formula with  3 points inside triangle
C    order = 5  -  quadrature formula with  7 points inside triangle
C    order = 6  -  quadrature formula with 12 points inside triangle
C    order = 9  -  quadrature formula with 19 points inside triangle
C    order =13  -  quadrature formula with 37 points inside triangle
C
C  The matrix A is defined by bilinear form (1). The following rules 
C  are applied for numbering basis functions (see details above):
C      A) First, basis function associated with vertices are enumerated
C         in the same order as the vertices XYi, i = 1,2 and 3
C
C      B) Second, basis function associated with edges are enumerated
C         in the same order as egdes 12,23 and 31, i.e. counter clockwise. 
C
C         If there are more than one degree of freedom per edge, they
C         are ordered by the distance from the lowest end-point.
C
C      C) The vector basis functions with N degrees of freedom per a
C         mesh object (vertex, edge, element) are enumerated in loops
C         by the corresponding mesh objects, e.g. P1,P2,P3,P1,P2,P3,...
C 
C      D) First, basis functions corresponding to x-component are enumerated.
C         Second, basis functions corresponding to y-component are enumerated.
C
C      D) High-order basis functions are enumerate by their momentum order.
C         No such finite elements ewre implemented yet.
C  
C     LDA   -  leading dimention of matrix A(LDA, LDA)
C     nRow  -  the number of rows of A
C     nCol  -  the number of columns of A
C
C
C ======================================================================
C *** DESCRIPTION OF THE TENSOR ***  
C ======================================================================
C  The external function D has the following STANDARD format
C
C  (3)     INTEGER FUNCTION D(x, y, label, dDATA, iDATA, iSYS, Diff)
C
C    The function returns type of the tendor Diff which is defined
C    in fem2Dtri.fd.
C
C    (x, y)   - Real*8 Cartesian coordinates of a 2D point where
C               tensor Diff should be evaluated
C
C    label    - identificator of a mesh element
C
C    dDATA    - Real*8  user given data (a number or an array)
C    iDATA    - Integer user given data (a number or an array)
C
C    iSYS(12) - system buffer for information exchange:
C           On Input:
C                iSYS(3) <- triangle number
C                iSYS(4) <- 1st vertex number
C                iSYS(5) <- 2nd vertex number
C                iSYS(6) <- 3rd vertex number
C
C                iSYS(7) <- 1st edge number, connects vertices 1 & 2
C                iSYS(8) <- 2nd edge number, connects vertices 2 & 3
C                iSYS(9) <- 3rd edge number, connects vertices 3 & 1
C
C                iSYS(10) <- total number of points,   nP
C                iSYS(11) <- total number of edges,    nR
C                iSYS(12) <- total number of elements, nE
C
C                iSYS(7:9) and iSYS(11) may be empty if the edge
C                degrees of freedom are not used. Two routines in
C                isys.f will help to encode and to decode iSYS().
C
C           On Output:
C                iD = iSYS(1) -> number of rows in tensor Diff 
C                jD = iSYS(2) -> number of columns in tensor Diff
C
C
C    Diff(iD, jD) - Real*8 matrix with leading dimension 4. In the case 
C           of vector finite elements U = (u, v) the following 
C           ordering should be used: u_x, u_y, v_x, v_y.
C                  
C           Examples. 
C              A) isotropic diffusion problem:  iD = jD = 1
C                 Diff = diffusion value at the point (x,y)
C                 the user returns TENSOR_SCALAR
C
C              B) anisotropic diffusion problem:  iD = jD = 2
C                 Diff(i, j) = diffusion value at the point (x,y)
C                 the user returns TENSOR_SYMMETRIC
C
C              C) convection problem:  iD = 1, jD = 2
C                 Diff(1, j) = velocity value at the point (x,y)
C                 the user returns TENSOR_GENERAL
C
C
C ======================================================================
C *** ADMISSIBLE (Yes), IMPOSSIBLE (NO), and TESTED (T) OPERATIONS ***
C ======================================================================
C
C                    IDEN  GRAD  DIV   CURL  DUDX  DUDY  DUDN
C                  -------------------------------------------
C     FEM_P0          Yes   Yes   NO    Yes   Yes   Yes   Yes
C     FEM_P1          Yes   Yes   NO    Yes   Yes   Yes   Yes
C     FEM_P2          Yes   Yes   NO    Yes   Yes   Yes   Yes
C     FEM_P3          Yes   Yes   NO    Yes    T     T    Yes
C     FEM_P4          Yes   Yes   NO    Yes    T     T    Yes
C     FEM_P1vector    Yes   Yes   Yes   Yes    T     T    Yes
C     FEM_P2vector    Yes   Yes   Yes    T    Yes   Yes    T
C     FEM_P2reduced   Yes   Yes   Yes    T     T     T     T
C     FEM_MINI        Yes   Yes   Yes   Yes    T     T     T
C     FEM_RT0         Yes   Yes   Yes   Yes    T     T     T
C     FEM_RT1          T     T     T     T     T     T     T
C     FEM_BDM1        Yes   Yes   Yes   Yes    T     T     T
C     FEM_CR1         Yes   Yes   NO    Yes   Yes   Yes    T
C     FEM_CR1vector   Yes   Yes   Yes    T     T     T     T
C     HERMIT          Yes   Yes   NO     T     T     T     T
C     ARGYRIS          T     T    NO     T     T     T     T
C 
C ======================================================================
      Real*8   XY1(*), XY2(*), XY3(*)

      Integer  FEMtypeA, FEMtypeB, operatorA, operatorB
      Integer  label, order, LDA, D, nRow, nCol

      Real*8   dDATA(*)
      Integer  iDATA(*), iSYS(*)
      EXTERNAL D

      Real*8   A(LDA, *)

C ======================================================================
c Local variables
      Real*8  XYP(2, 3), XYN(2, 3)
      Real*8  det, vol, s
   
      Real*8  PSI(2, 2)
      Real*8  U(4, MaxDOFs, MaxPointGauss), V(4, MaxDOFs, MaxPointGauss)
      Real*8  Diff(4, 4, MaxPointGauss),   DU(4, MaxDOFs, MaxPointGauss)

      Real*8  w(MaxPointGauss), XYG(2, MaxPointGauss)
      Real*8  XYL(3, AllPointGauss)

      Integer i,j,k,n, nfa,nfb, idim,jdim, iD,jD, tensor, iGauss, iL
      Logical FEMtype, operator, ifXtensor

C ======================================================================
      DATA XYL/3 * T1A, 
c ... 3 points (order 2)
     &  T2A,T2B,T2B,  T2B,T2A,T2B,  T2B,T2B,T2A,
c ... 7 points (order 5)
     &  3 * T5A,      T5B,T5C,T5C,  T5C,T5B,T5C,  T5C,T5C,T5B,
     &                T5D,T5E,T5E,  T5E,T5D,T5E,  T5E,T5E,T5D,
c ... 12 points (order 6)
     &  T6A,T6B,T6B,  T6B,T6A,T6B,  T6B,T6B,T6A,  
     &  T6C,T6D,T6D,  T6D,T6C,T6D,  T6D,T6D,T6C,  
     &  T6E,T6F,T6G,  T6G,T6E,T6F,  T6F,T6G,T6E,  
     &  T6E,T6G,T6F,  T6G,T6F,T6E,  T6F,T6E,T6G,
c ... 19 points (order 9)
     &  3 * T9A,      T9B,T9C,T9C,  T9C,T9B,T9C,  T9C,T9C,T9B,  
     &  T9D,T9E,T9E,  T9E,T9D,T9E,  T9E,T9E,T9D,
     &  T9F,T9G,T9G,  T9G,T9F,T9G,  T9G,T9G,T9F,
     &  T9H,T9I,T9I,  T9I,T9H,T9I,  T9I,T9I,T9H,
     &  T9J,T9K,T9L,  T9L,T9J,T9K,  T9K,T9L,T9J,  T9J,T9L,T9K,  
     &                              T9L,T9K,T9J,  T9K,T9J,T9L,
c ... 37 points (order 13)
     &  3 * TXA,      TXB,TXC,TXC,  TXC,TXB,TXC,  TXC,TXC,TXB,
     &  TXD,TXE,TXE,  TXE,TXD,TXE,  TXE,TXE,TXD,  
     &  TXF,TXG,TXG,  TXG,TXF,TXG,  TXG,TXG,TXF,  
     &  TXH,TXI,TXI,  TXI,TXH,TXI,  TXI,TXI,TXH,
     &  TXJ,TXK,TXK,  TXK,TXJ,TXK,  TXK,TXK,TXJ,  
     &  TXL,TXM,TXM,  TXM,TXL,TXM,  TXM,TXM,TXL,
     &  TXN,TXO,TXP,  TXP,TXN,TXO,  TXO,TXP,TXN,  TXN,TXP,TXO,  
     &                              TXP,TXO,TXN,  TXO,TXN,TXP,
     &  TXR,TXS,TXT,  TXT,TXR,TXS,  TXS,TXT,TXR,  TXR,TXT,TXS,
     &                              TXT,TXS,TXR,  TXS,TXR,TXT,
     &  TXU,TXV,TXW,  TXW,TXU,TXV,  TXV,TXW,TXU,  TXU,TXW,TXV,
     &                              TXW,TXV,TXU,  TXV,TXU,TXW/

C ======================================================================
      If(order.LE.0 .OR. order.GT.13) 
     &  Call errMesFEM(2001, 'fem2Dtri', 'quadrature order is wrong')

      FEMtype = FEMtypeA.EQ.FEMtypeB
      operator = operatorA.EQ.operatorB

c ... transformation of variables y = PSI * (x - x_0)
      Do i = 1, 2
         XYP(i, 1) = 0D0
         XYP(i, 2) = XY2(i) - XY1(i)
         XYP(i, 3) = XY3(i) - XY1(i)
      End do

      Call solve2x2(XYP(1, 2), XYP(2, 2), PSI(1, 1),
     &              XYP(1, 3), XYP(2, 3), PSI(1, 2), det)

      Call solve2x2(XYP(1, 3), XYP(2, 3), PSI(2, 1),
     &              XYP(1, 2), XYP(2, 2), PSI(2, 2), det)

c ... weights and points
      vol = dabs(det) / 2
      Call WeightsPoints(XY1, XY2, XY3, vol, order, XYG, w, iGauss, iL)

c ... exterior normal vectors
      If(FEMtypeA.EQ.FEM_P2reduced .OR. 
     &   FEMtypeB.EQ.FEM_P2reduced .OR.
     &   FEMtypeA.EQ.FEM_ARGYRIS) Then
         Call calNormalExt(XY1, XY2, XY3, XYN(1, 1))
         Call calNormalExt(XY2, XY3, XY1, XYN(1, 2))
         Call calNormalExt(XY3, XY1, XY2, XYN(1, 3))
      End if

c ... compute operatorA * FEMtypeA
      If(operatorA.EQ.GRAD) Then
         Call applyGRAD(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, XYN)

      Else If(operatorA.EQ.DIV) Then
         Call applyDIV( iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, XYN, det)

      Else If(operatorA.EQ.IDEN) Then
         Call applyIDEN(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, XYN, det)

      Else If(operatorA.EQ.CURL) Then
         Call applyCURL(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, det)

      Else If(operatorA.EQ.DUDX) Then
         Call applyDUDX(iGauss, XYL(1, iL), PSI, FEMtypeA, nfa, idim, U)

      Else If(operatorA.EQ.DUDY) Then
         Call applyDUDY(iGauss, XYL(1, iL), PSI, FEMtypeA, nfa, idim, U)

      Else
         Call errMesFEM(2001, 'fem2Dtri', 'operatorA is not supported')
      End if

      If(nfa.GT.LDA) Call errMesFEM(2001, 'fem2Dtri',
     &     'the local matrix leading dimension, LDA, is too small')


c ... compute operatorB * FEMtypeB
      nfb = nfa
      jdim = idim
      if(operator .AND. FEMtype) Goto 100

      If(operatorB.EQ.GRAD) Then
         Call applyGRAD(iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                  nfb, jdim, V, XYP, XYN)

      Else If(operatorB.EQ.DIV) Then
         Call applyDIV( iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                  nfb, jdim, V, XYP, XYN, det)

      Else If(operatorB.EQ.IDEN) Then
         Call applyIDEN(iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                  nfb, jdim, V, XYP, XYN, det)

      Else If(operatorB.EQ.CURL) Then
         Call applyCURL(iGauss, XYL(1, iL), PSI, FEMtypeB, 
     &                  nfb, jdim, V, XYP, det)

      Else If(operatorB.EQ.DUDX) Then
         Call applyDUDX(iGauss, XYL(1, iL), PSI, FEMtypeB, nfb, jdim, V)

      Else If(operatorB.EQ.DUDY) Then
         Call applyDUDY(iGauss, XYL(1, iL), PSI, FEMtypeB, nfb, jdim, V)

      Else
         Call errMesFEM(2001, 'fem2Dtri', 'operatorB is not supported')
      End if

c     If(nfb.GT.LDA) Call errMesFEM(2001, 'fem2Dtri',
c    &     'the local matrix second dimension, LDA, is too small')
 

c ... compute D * U
 100  iD = jdim
      jD = idim

      Do n = 1, iGauss
         tensor = D(XYG(1,n), XYG(2,n), label, 
     &              dDATA, iDATA, iSYS, Diff(1, 1, n))
         If(ifXtensor(tensor, TENSOR_NULL)) Goto 200
      End do 
      If(ifXtensor(tensor, BC_DIRICHLET) .OR. 
     &   ifXtensor(tensor, BC_NEUMANN) .OR. 
     &   ifXtensor(tensor, BC_ROBIN) .OR. 
     &   ifXtensor(tensor, BC_ROBIN_COEF) .OR. 
     &   ifXtensor(tensor, BC_NULL)) tensor = TENSOR_GENERAL


 200  Continue
      If(ifXtensor(tensor, TENSOR_NULL)) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 
     &        'fem2Dtri', 'Operators A and B are not compatible')

         Do n = 1, iGauss
            Do i = 1, idim
               Do k = 1, nfa
                  DU(i, k, n) = U(i, k, n) * w(n)
               End do
            End do
         End do
      Else If(ifXtensor(tensor, TENSOR_SCALAR)) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 'fem2Dtri', 
     &        'Operators A and B are not compatible')

         Do n = 1, iGauss
            s = Diff(1, 1, n) * w(n) 
            Do i = 1, idim
               Do k = 1, nfa
                  DU(i, k, n) = U(i, k, n) * s
               End do
            End do
         End do
      Else If(ifXtensor(tensor, TENSOR_SYMMETRIC) .OR. 
     &        ifXtensor(tensor, TENSOR_GENERAL)) Then
         iD = iSYS(1)
         jD = iSYS(2)

         If(jD.NE.idim .OR. iD.NE.jdim) Call errMesFEM(2001, 
     &        'fem2Dtri', 'the operators A and B are not compatible')

         Do n = 1, iGauss
            Do i = 1, iD
               Do k = 1, nfa
                  s = 0D0
                  Do j = 1, jD
                     s = s + Diff(i, j, n) * U(j, k, n)
                  End do
                  DU(i, k, n) = s * w(n)
               End do
            End do
         End do
      Else
         Call errMesFEM(2001, 'fem2Dtri', 
     &        'the given tensor is not supported') 
      End if


c ... compute <D U, V>
      If(operator .AND. 
     &   FEMtype .AND. 
     &   .NOT. ifXtensor(tensor, TENSOR_GENERAL)) Then
         Do i = 1, nfa
            Do j = 1, i - 1
               A(i, j) = A(j, i)
            End do

            Do j = i, nfa
               s = 0D0
               Do k = 1, iD
                  Do n = 1, iGauss
                     s = s + DU(k, i, n) * U(k, j, n)
                  End do
               End do
               A(i, j) = s
            End do
         End do
      Else If(operator .AND. FEMtype) Then
         Do i = 1, nfa
            Do j = 1, nfb
               s = 0D0
               Do k = 1, iD
                  Do n = 1, iGauss
                     s = s + DU(k, i, n) * U(k, j, n)
                  End do
               End do
               A(i, j) = s
            End do
         End do
      Else 
         Do i = 1, nfa
            Do j = 1, nfb
               s = 0D0
               Do k = 1, iD
                  Do n = 1, iGauss
                     s = s + DU(k, i, n) * V(k, j, n)
                  End do
               End do
               A(i, j) = s
            End do
         End do
      End if

      nRow = nfa
      nCol = nfb

      Return
      End



c ======================================================================
      Subroutine WeightsPoints(XY1, XY2, XY3, vol, order, 
     &                         XYG, w, iGauss, iLref)
c ======================================================================
      implicit none
      include 'fem2Dtri.fd'
c ======================================================================
C The procedure is used for effective computing points for numerical
C integration, b/c of a symmetry.
C
C Remark: A 1-to-1 coresspondance between XYL and XYG should be hold.
c ======================================================================
      Real*8  XY1(2), XY2(2), XY3(2), vol
      Integer order

      Real*8  XYG(2, *), w(*)
      Integer iGauss, iLref, i

c ======================================================================
      If(order.EQ.1) Then
         iGauss = LDG1
         iLref = 1
         w(1) = W1A * vol

         Do i = 1, 2
            XYG(i, 1) = T1A * (XY1(i) + XY2(i) + XY3(i))
         End do

      Else If(order.EQ.2) Then
         iGauss = LDG2
         iLref = LDG1 + 1
         Do i = 1, iGauss
            w(i) = W2A * vol
         End do

         Do i = 1, 2
            XYG(i, 1) = T2B * (XY2(i) + XY3(i)) 
            XYG(i, 2) = T2B * (XY1(i) + XY3(i)) 
            XYG(i, 3) = T2B * (XY1(i) + XY2(i)) 
         End do

      Else If(order.LE.5) Then
         iGauss = LDG5
         iLref = LDG1 + LDG2 + 1

         w(1) = W5A * vol
         Do i = 1, 2
            XYG(i, 1) = T5A * (XY1(i) + XY2(i) + XY3(i))
         End do

         Do i = 2, 4
            w(i) = W5B * vol
         End do
         Do i = 1, 2
            XYG(i, 2) = T5B * XY1(i) + T5C * (XY2(i) + XY3(i)) 
            XYG(i, 3) = T5B * XY2(i) + T5C * (XY1(i) + XY3(i)) 
            XYG(i, 4) = T5B * XY3(i) + T5C * (XY1(i) + XY2(i)) 
         End do

         Do i = 5, 7
            w(i) = W5C * vol
         End do
         Do i = 1, 2
            XYG(i, 5) = T5D * XY1(i) + T5E * (XY2(i) + XY3(i)) 
            XYG(i, 6) = T5D * XY2(i) + T5E * (XY1(i) + XY3(i)) 
            XYG(i, 7) = T5D * XY3(i) + T5E * (XY1(i) + XY2(i)) 
         End do

      Else If(order.LE.6) Then
         iGauss = LDG6
         iLref = LDG1 + LDG2 + LDG5 + 1

         Do i = 1, 3
            w(i) = W6A * vol
         End do
         Do i = 1, 2
            XYG(i, 1) = T6A * XY1(i) + T6B * (XY2(i) + XY3(i)) 
            XYG(i, 2) = T6A * XY2(i) + T6B * (XY1(i) + XY3(i)) 
            XYG(i, 3) = T6A * XY3(i) + T6B * (XY1(i) + XY2(i)) 
         End do

         Do i = 4, 6
            w(i) = W6B * vol
         End do
         Do i = 1, 2
            XYG(i, 4) = T6C * XY1(i) + T6D * (XY2(i) + XY3(i)) 
            XYG(i, 5) = T6C * XY2(i) + T6D * (XY1(i) + XY3(i)) 
            XYG(i, 6) = T6C * XY3(i) + T6D * (XY1(i) + XY2(i)) 
         End do

         Do i = 7, 12
            w(i) = W6C * vol
         End do
         Do i = 1, 2
            XYG(i, 7) = T6E * XY1(i) + T6F * XY2(i) + T6G * XY3(i)
            XYG(i, 8) = T6G * XY1(i) + T6E * XY2(i) + T6F * XY3(i)
            XYG(i, 9) = T6F * XY1(i) + T6G * XY2(i) + T6E * XY3(i)

            XYG(i,10) = T6E * XY1(i) + T6G * XY2(i) + T6F * XY3(i)
            XYG(i,11) = T6G * XY1(i) + T6F * XY2(i) + T6E * XY3(i)
            XYG(i,12) = T6F * XY1(i) + T6E * XY2(i) + T6G * XY3(i)
         End do

      Else If(order.LE.9) Then
         iGauss = LDG9
         iLref = LDG1 + LDG2 + LDG5 + LDG6 + 1

         w(1) = W9A * vol
         Do i = 1, 2
            XYG(i, 1) = T9A * (XY1(i) + XY2(i) + XY3(i))
         End do

         Do i = 2, 4
            w(i) = W9B * vol
         End do
         Do i = 1, 2
            XYG(i, 2) = T9B * XY1(i) + T9C * (XY2(i) + XY3(i)) 
            XYG(i, 3) = T9B * XY2(i) + T9C * (XY1(i) + XY3(i)) 
            XYG(i, 4) = T9B * XY3(i) + T9C * (XY1(i) + XY2(i)) 
         End do

         Do i = 5, 7
            w(i) = W9C * vol
         End do
         Do i = 1, 2
            XYG(i, 5) = T9D * XY1(i) + T9E * (XY2(i) + XY3(i)) 
            XYG(i, 6) = T9D * XY2(i) + T9E * (XY1(i) + XY3(i)) 
            XYG(i, 7) = T9D * XY3(i) + T9E * (XY1(i) + XY2(i)) 
         End do

         Do i = 8, 10 
            w(i) = W9D * vol
         End do
         Do i = 1, 2
            XYG(i,  8) = T9F * XY1(i) + T9G * (XY2(i) + XY3(i)) 
            XYG(i,  9) = T9F * XY2(i) + T9G * (XY1(i) + XY3(i)) 
            XYG(i, 10) = T9F * XY3(i) + T9G * (XY1(i) + XY2(i)) 
         End do
 
         Do i = 11, 13 
            w(i) = W9E * vol
         End do
         Do i = 1, 2
            XYG(i, 11) = T9H * XY1(i) + T9I * (XY2(i) + XY3(i)) 
            XYG(i, 12) = T9H * XY2(i) + T9I * (XY1(i) + XY3(i)) 
            XYG(i, 13) = T9H * XY3(i) + T9I * (XY1(i) + XY2(i)) 
         End do
 
         Do i = 14, 19
            w(i) = W9F * vol
         End do
         Do i = 1, 2
            XYG(i, 14) = T9J * XY1(i) + T9K * XY2(i) + T9L * XY3(i)
            XYG(i, 15) = T9L * XY1(i) + T9J * XY2(i) + T9K * XY3(i)
            XYG(i, 16) = T9K * XY1(i) + T9L * XY2(i) + T9J * XY3(i)

            XYG(i, 17) = T9J * XY1(i) + T9L * XY2(i) + T9K * XY3(i)
            XYG(i, 18) = T9L * XY1(i) + T9K * XY2(i) + T9J * XY3(i)
            XYG(i, 19) = T9K * XY1(i) + T9J * XY2(i) + T9L * XY3(i)
         End do

      Else If(order.LE.13) Then
         iGauss = LDGX
         iLref = LDG1 + LDG2 + LDG5 + LDG6 + LDG9 + 1

         w(1) = WXA * vol
         Do i = 1, 2
            XYG(i, 1) = TXA * (XY1(i) + XY2(i) + XY3(i))
         End do

         Do i = 2, 4
            w(i) = WXB * vol
         End do
         Do i = 1, 2
            XYG(i, 2) = TXB * XY1(i) + TXC * (XY2(i) + XY3(i)) 
            XYG(i, 3) = TXB * XY2(i) + TXC * (XY1(i) + XY3(i)) 
            XYG(i, 4) = TXB * XY3(i) + TXC * (XY1(i) + XY2(i)) 
         End do

         Do i = 5, 7
            w(i) = WXC * vol
         End do
         Do i = 1, 2
            XYG(i, 5) = TXD * XY1(i) + TXE * (XY2(i) + XY3(i)) 
            XYG(i, 6) = TXD * XY2(i) + TXE * (XY1(i) + XY3(i)) 
            XYG(i, 7) = TXD * XY3(i) + TXE * (XY1(i) + XY2(i)) 
         End do

         Do i = 8, 10 
            w(i) = WXD * vol
         End do
         Do i = 1, 2
            XYG(i,  8) = TXF * XY1(i) + TXG * (XY2(i) + XY3(i)) 
            XYG(i,  9) = TXF * XY2(i) + TXG * (XY1(i) + XY3(i)) 
            XYG(i, 10) = TXF * XY3(i) + TXG * (XY1(i) + XY2(i)) 
         End do
 
         Do i = 11, 13 
            w(i) = WXE * vol
         End do
         Do i = 1, 2
            XYG(i, 11) = TXH * XY1(i) + TXI * (XY2(i) + XY3(i)) 
            XYG(i, 12) = TXH * XY2(i) + TXI * (XY1(i) + XY3(i)) 
            XYG(i, 13) = TXH * XY3(i) + TXI * (XY1(i) + XY2(i)) 
         End do
 
         Do i = 14, 16 
            w(i) = WXF * vol
         End do
         Do i = 1, 2
            XYG(i, 14) = TXJ * XY1(i) + TXK * (XY2(i) + XY3(i)) 
            XYG(i, 15) = TXJ * XY2(i) + TXK * (XY1(i) + XY3(i)) 
            XYG(i, 16) = TXJ * XY3(i) + TXK * (XY1(i) + XY2(i)) 
         End do
 
         Do i = 17, 19 
            w(i) = WXG * vol
         End do
         Do i = 1, 2
            XYG(i, 17) = TXL * XY1(i) + TXM * (XY2(i) + XY3(i)) 
            XYG(i, 18) = TXL * XY2(i) + TXM * (XY1(i) + XY3(i)) 
            XYG(i, 19) = TXL * XY3(i) + TXM * (XY1(i) + XY2(i)) 
         End do
 
         Do i = 20, 25
            w(i) = WXH * vol
         End do
         Do i = 1, 2
            XYG(i, 20) = TXN * XY1(i) + TXO * XY2(i) + TXP * XY3(i)
            XYG(i, 21) = TXP * XY1(i) + TXN * XY2(i) + TXO * XY3(i)
            XYG(i, 22) = TXO * XY1(i) + TXP * XY2(i) + TXN * XY3(i)

            XYG(i, 23) = TXN * XY1(i) + TXP * XY2(i) + TXO * XY3(i)
            XYG(i, 24) = TXP * XY1(i) + TXO * XY2(i) + TXN * XY3(i)
            XYG(i, 25) = TXO * XY1(i) + TXN * XY2(i) + TXP * XY3(i)
         End do
 
         Do i = 26,31 
            w(i) = WXI * vol
         End do
         Do i = 1, 2
            XYG(i, 26) = TXR * XY1(i) + TXS * XY2(i) + TXT * XY3(i)
            XYG(i, 27) = TXT * XY1(i) + TXR * XY2(i) + TXS * XY3(i)
            XYG(i, 28) = TXS * XY1(i) + TXT * XY2(i) + TXR * XY3(i)

            XYG(i, 29) = TXR * XY1(i) + TXT * XY2(i) + TXS * XY3(i)
            XYG(i, 30) = TXT * XY1(i) + TXS * XY2(i) + TXR * XY3(i)
            XYG(i, 31) = TXS * XY1(i) + TXR * XY2(i) + TXT * XY3(i)
         End do
 
         Do i = 32,37 
            w(i) = WXJ * vol
         End do
         Do i = 1, 2
            XYG(i, 32) = TXU * XY1(i) + TXV * XY2(i) + TXW * XY3(i)
            XYG(i, 33) = TXW * XY1(i) + TXU * XY2(i) + TXV * XY3(i)
            XYG(i, 34) = TXV * XY1(i) + TXW * XY2(i) + TXU * XY3(i)

            XYG(i, 35) = TXU * XY1(i) + TXW * XY2(i) + TXV * XY3(i)
            XYG(i, 36) = TXW * XY1(i) + TXV * XY2(i) + TXU * XY3(i)
            XYG(i, 37) = TXV * XY1(i) + TXU * XY2(i) + TXW * XY3(i)
         End do
      End if
      Return
      End



c ======================================================================
      Real*8 Function Lfun(i, x, y)
c ======================================================================
      implicit none
      Integer  i
      Real*8   x, y

      If(i.EQ.1) Then
         Lfun = 1D0 - x - y
      Else If(i.EQ.2) Then
         Lfun = x
      Else If(i.EQ.3) Then
         Lfun = y
      End if
      Return
      End



c ======================================================================
      Subroutine solve2x2(a11, a12, a,
     &                    a21, a22, b, det)
c ======================================================================
      implicit none
      Real*8  a11, a12, a
      Real*8  a21, a22, b, det

c Local variables
      Real*8  s

c ======================================================================
      det = a11 * a22 - a21 * a12
      s = 1D0 / det

      a = a22 * s
      b =-a21 * s
      Return
      End



c ======================================================================
      Subroutine invert3x3(A, B, det)
c ======================================================================
      implicit none
      Real*8  A(3, 3), B(3, 3), det

c Local variables
      Real*8  s11, s12, s13

c ======================================================================
      s11 = A(2,2) * A(3,3) - A(3,2) * A(2,3)
      s12 = A(2,1) * A(3,3) - A(3,1) * A(2,3) 
      s13 = A(2,1) * A(3,2) - A(3,1) * A(2,2)

      det = A(1,1) * s11 - A(1,2) * s12 + A(1,3) * s13

      B(1,1) = s11 / det
      B(2,1) =-s12 / det
      B(3,1) = s13 / det

      B(1,2) =-(A(1,2) * A(3,3) - A(3,2) * A(1,3)) / det
      B(2,2) = (A(1,1) * A(3,3) - A(3,1) * A(1,3)) / det
      B(3,2) =-(A(1,1) * A(3,2) - A(3,1) * A(1,2)) / det

      B(1,3) = (A(1,2) * A(2,3) - A(2,2) * A(1,3)) / det
      B(2,3) =-(A(1,1) * A(2,3) - A(2,1) * A(1,3)) / det
      B(3,3) = (A(1,1) * A(2,2) - A(2,1) * A(1,2)) / det

      Return
      End



C ======================================================================
      Subroutine calNormalVec(xy1, xy2, xyn)
C ======================================================================
C Routine computes a normal vector to the oriented edge xy1 -> xy2
C ======================================================================
      implicit none
      Real*8 xy1(2), xy2(2), xyn(2)
      Real*8 ax, ay

      ax = xy2(1) - xy1(1)
      ay = xy2(2) - xy1(2)

      xyn(1) =  ay
      xyn(2) = -ax
     
      Return
      End



C ======================================================================
      Subroutine calNormalExt(xy1, xy2, xy3, xyn)
C ======================================================================
C Routines compute external normal vectors to the edge {xy1, xy2}
C of triangle {xy1, xy2, xy3}. This is the twin sister of a similar
C routine from package aniMBA.
C ======================================================================
      implicit none
      Real*8 xy1(2), xy2(2), xy3(2), xyn(2)
      Real*8 x, y, d

      x = xy2(1) - xy1(1)
      y = xy2(2) - xy1(2)

      d = dsqrt(x * x + y * y)
 
      xyn(1) = -y / d
      xyn(2) =  x / d

c ... orientation
      x = xy3(1) - xy1(1)
      y = xy3(2) - xy1(2)

      If(x*xyn(1) + y*xyn(2).GT.0D0) Then
         xyn(1) = -xyn(1)
         xyn(2) = -xyn(2)
      End if 

      Return
      End



C ======================================================================
      Subroutine copyGauss(iGauss, nfa, idim, U)
C ======================================================================
      implicit none
      include "fem2Dtri.fd"

      Integer iGauss, nfa, idim
      Real*8  U(4, MaxDOFs, *)

      Integer i, k, n

      Do n = 2, iGauss
         Do i = 1, nfa
            Do k = 1, idim
               U(k, i, n) = U(k, i, 1)
            End do
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine clearU(iGauss, nfa, idim, U)
C ======================================================================
      implicit none
      include "fem2Dtri.fd"

      Integer iGauss, nfa, idim
      Real*8  U(4, MaxDOFs, *)

      Integer i, k, n

      Do n = 1, iGauss
         Do i = 1, nfa
            Do k = 1, idim
               U(k, i, n) = 0D0
            End do
         End do
      End do

      Return
      End

