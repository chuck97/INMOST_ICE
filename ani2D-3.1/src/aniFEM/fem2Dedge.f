C ======================================================================
      Subroutine fem2Dedge(XY1, XY2, XY3,
     &                     operatorA, FEMtypeA, operatorB, FEMtypeB, 
     &                     label, D, dDATA, iDATA, iSYS, order, 
     &                     LDA, A, nRow, nCol)
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'
C ======================================================================
C Routine implements a bilinear form on an edge {XY1,XY2} of triangle
C {XY1,XY2,XY3}. 
C The notations and order of unknowns are described in fem2Dtri.f
C ======================================================================
      Real*8   XY1(2), XY2(2), XY3(2)

      Real*8   dDATA(*)
      Integer  iDATA(*), label, D, iSYS(*)
      EXTERNAL D

      Integer  FEMtypeA, FEMtypeB, operatorA, operatorB
      Integer  order, LDA, nRow, nCol

      Real*8   A(LDA, *)

C ======================================================================
C Local variables
      Real*8  XYP(2, 3), XYN(2)

      Real*8  PSI(2, 2)
      Real*8  U(2, MaxDOFs, MaxPnt2DGauss), V(2, MaxDOFs, MaxPnt2DGauss)
      Real*8  Diff(2, 2, MaxPnt2DGauss),   DU(2, MaxDOFs, MaxPnt2DGauss)

      Real*8  w(MaxPnt2DGauss), XYG(2, MaxPnt2DGauss)
      Real*8  XYL(2, AllPnt2DGauss)

      Integer i,j,k,n, nfa,nfb, idim,jdim, iD,jD, iGauss, iL, tensor
      Real*8  det, s
      Logical FEMtype, operator, ifXtensor

C ======================================================================
      DATA XYL/0.5D0, 0.5D0,
c ... 3 points (order 5)
     &         0.112701665379258D0, 0.887298334620742D0, 0.5D0, 0.5D0, 
     &         0.887298334620742D0, 0.112701665379258D0,
c ... 5 points (order 9)
     &         0.046910077030668D0, 0.953089922969332D0, 
     &         0.230765344947158D0, 0.769234655052841D0, 0.5D0, 0.5D0,
     &         0.769234655052841D0, 0.230765344947158D0, 
     &         0.953089922969332D0, 0.046910077030668D0/

C ======================================================================
      If(order.LE.0 .OR. order.GT.9) 
     &  Call errMesFEM(2001, 'fem2Dedge', 'quadrature order is wrong')

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
      Call Weights1D(XY1, XY2, order, XYL, XYG, w, iGauss, iL)


c ... compute operatorA * FEMtypeA
      If(operatorA.EQ.IDEN) Then
         Call applyTRACE(iGauss, XYL(1,iL), FEMtypeA, nfa, idim, U)
      Else If(operatorA.EQ.DUDN) Then
         Call calNormalExt(XY1, XY2, XY3, XYN)
         Call applyDUDN(iGauss, XYL(1,iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYN)
      Else
         Call errMesFEM(2001, 'fem2Dedge', 'operatorA is not supported')
      End if
      If(nfa.GT.LDA) Call errMesFEM(2001, 'fem2Dedge',
     &     'local matrix leading dimension, LDA, is too small')


c ... compute operatorB * FEMtypeB
      nfb = nfa
      jdim = idim
      if(operator .AND. FEMtype) Goto 100

      If(operatorB.EQ.IDEN) Then
         Call applyTRACE(iGauss, XYL(1,iL), FEMtypeB, nfb, jdim, V)
      Else If(operatorB.EQ.DUDN) Then
         Call calNormalExt(XY1, XY2, XY3, XYN)
         Call applyDUDN(iGauss, XYL(1,iL), PSI, FEMtypeB, 
     &                  nfb, jdim, V, XYN)
      Else
         Call errMesFEM(2001, 'fem2Dedge', 'operatorB is not supported')
      End if
 

c ... compute D * U
 100  iD = jdim
      jD = idim

      Do n = 1, iGauss
         tensor = D(XYG(1, n), XYG(2, n), label, 
     &              dDATA, iDATA, iSYS, Diff(1,1,n))
         If(ifXtensor(tensor, TENSOR_NULL)) Goto 200
      End do


 200  Continue
      If(ifXtensor(tensor, TENSOR_NULL)) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 
     &        'fem2Dedge', 'Operators A and B are not compatible')

         Do n = 1, iGauss
            Do i = 1, idim
               Do k = 1, nfa
                  DU(i, k, n) = U(i, k, n) * w(n)
               End do
            End do
         End do

      Else If(ifXtensor(tensor, TENSOR_SCALAR)) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 'fem2Dedge', 
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
     &        'fem2Dedge', 'the operators A and B are not compatible')

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
      End if


c ... compute <D U, V>
      If(operator .AND. FEMtype) Then
         Do i = 1, nfa
            Do j = 1, nfb
               s = 0D0
               Do k = 1, iD
                  Do n = 1, iGauss
                     s = s + DU(k, i, n) * U(k, j, n)
                  End do
                  A(i, j) = s
               End do
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
                  A(i, j) = s
               End do
            End do
         End do
      End if

      nRow = nfa
      nCol = nfb

      Return
      End



c ======================================================================
      Subroutine Weights1D(XY1, XY2, order, XYL, XYG, w, iGauss, iLref)
c ======================================================================
      implicit none
      include 'fem2Dtri.fd'
c ======================================================================
      Real*8  XY1(2), XY2(2), XYL(2, *), XYG(2, *), w(*)
      Integer order, iGauss, iLref

      Real*8 length
      Integer i, n, m
c ======================================================================
      length = dsqrt((XY2(1) - XY1(1)) ** 2 + (XY2(2) - XY1(2)) ** 2) 

      If(order.EQ.1) Then
         iGauss = 1
         iLref = 1
         w(1) = length

         Do i = 1, 2
            XYG(i, 1) = (XY1(i) + XY2(i)) / 2
         End do

      Else If(order.LE.5) Then
         iGauss = 3
         iLref = 2

         w(1) = length * 5D0 / 18
         w(2) = length * 8D0 / 18
         w(3) = length * 5D0 / 18

         Do n = 1, iGauss 
            m = iLref + n - 1
            Do i = 1, 2
               XYG(i, n) = XYL(1, m) * XY1(i) + XYL(2, m) * XY2(i)
            End do
         End do

      Else If(order.LE.9) Then
         iGauss = 5
         iLref = 5

         w(1) = length * 0.118463442528095D0
         w(2) = length * 0.239314335249683D0
         w(3) = length * 0.284444444444444D0
         w(4) = w(2) 
         w(5) = w(1) 

         Do n = 1, iGauss 
            m = iLref + n - 1
            Do i = 1, 2
               XYG(i, n) = XYL(1, m) * XY1(i) + XYL(2, m) * XY2(i)
            End do
         End do
      End if

      Return
      End



