C ======================================================================
      Subroutine fem2Dsub(XY1, XY2, XY3,  P1, P2, P3,  Q1, Q2, Q3,
     &                    operatorA, FEMtypeA, operatorB, FEMtypeB, 
     &                    label, D, dDATA, iDATA, iSYS, order,  
     &                    LDA, A, nRow, nCol)
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'
C ======================================================================
C  *** ITRODUCTION ***
C  The routine computes elemental matrix for the bilinear form 
C
C  (1)            <OpA(u), OpB(v)>               
C
C  where OpA and OpB are linear operators, and u and v are finete
C  element functions, u in "A", and v in "B". The integration is made
C  over a sub-triangle XY1,XY2,XY3 which is a part of triangles P1,P2,P3
C  and Q1,Q2,Q3.
C
C  Remarks.
C    1. Update this routine together with the family
C 
C ======================================================================
      Real*8   XY1(*), XY2(*), XY3(*)
      Real*8   P1(*), P2(*), P3(*),  Q1(*), Q2(*), Q3(*)

      Integer  FEMtypeA, FEMtypeB, operatorA, operatorB
      Integer  label, order, LDA, D, nRow, nCol

      Real*8   dDATA(*)
      Integer  iDATA(*), iSYS(*)
      EXTERNAL D

      Real*8   A(LDA, *)

C ======================================================================
c Local variables
      Real*8  XYP(2, 3), XYQ(2, 3), XYN(2, 3), XYM(2, 3)
      Real*8  detp, detq, vol, s
   
      Real*8  PSI(2, 2), QSI(2, 2)
      Real*8  U(4, MaxDOFs, MaxPointGauss), V(4, MaxDOFs, MaxPointGauss)
      Real*8  Diff(4, 4, MaxPointGauss), DU(4, MaxDOFs, MaxPointGauss)

      Real*8  w(MaxPointGauss), XYG(2, MaxPointGauss)
      Real*8  XYL(3, AllPointGauss)
      Real*8  XPL(3, MaxPointGauss), XQL(3, MaxPointGauss)

      Integer i,j,k,n, nfa,nfb, idim,jdim, iD,jD, tensor, iGauss, iL
      Logical ifXtensor

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

C ================================================================
      If(order.LE.0 .OR. order.GT.13) 
     &  Call errMesFEM(2001, 'fem2Dsub', 'quadrature order is wrong')

c ... transformation of variables y = PSI * (x - x_0)
      Do i = 1, 2
         XYP(i, 1) = 0D0
         XYP(i, 2) = P2(i) - P1(i)
         XYP(i, 3) = P3(i) - P1(i)
      End do

      Call solve2x2(XYP(1, 2), XYP(2, 2), PSI(1, 1),
     &              XYP(1, 3), XYP(2, 3), PSI(1, 2), detp)

      Call solve2x2(XYP(1, 3), XYP(2, 3), PSI(2, 1),
     &              XYP(1, 2), XYP(2, 2), PSI(2, 2), detp)


c ... transformation of variables y = QSI * (x - x_0)
      Do i = 1, 2
         XYQ(i, 1) = 0D0
         XYQ(i, 2) = Q2(i) - Q1(i)
         XYQ(i, 3) = Q3(i) - Q1(i)
      End do

      Call solve2x2(XYQ(1, 2), XYQ(2, 2), QSI(1, 1),
     &              XYQ(1, 3), XYQ(2, 3), QSI(1, 2), detq)

      Call solve2x2(XYQ(1, 3), XYQ(2, 3), QSI(2, 1),
     &              XYQ(1, 2), XYQ(2, 2), QSI(2, 2), detq)


c ... weights and points for subtriangle
      vol = dabs((XY2(1) - XY1(1)) * (XY3(2) - XY1(2)) 
     &         - (XY2(2) - XY1(2)) * (XY3(1) - XY1(1))) / 2

      Call WeightsPoints(XY1, XY2, XY3, vol, order, XYG, w, iGauss, iL)


c ... exterior normal vectors
      If(FEMtypeA.EQ.FEM_P2reduced) Then
         Call calNormalExt(P1, P2, P3, XYN(1, 1))
         Call calNormalExt(P2, P3, P1, XYN(1, 2))
         Call calNormalExt(P3, P1, P2, XYN(1, 3))
      End if

      If(FEMtypeB.EQ.FEM_P2reduced) Then
         Call calNormalExt(Q1, Q2, Q3, XYM(1, 1))
         Call calNormalExt(Q2, Q3, Q1, XYM(1, 2))
         Call calNormalExt(Q3, Q1, Q2, XYM(1, 3))
      End if


c ... calculate barycentric coordinates for parents
      Call barycentric(XY1, XY2, XY3,  P1, XYP, 
     &                 iGauss, XYL(1,iL), XPL, detp)
      Call barycentric(XY1, XY2, XY3,  Q1, XYQ, 
     &                 iGauss, XYL(1,iL), XQL, detq)


c ... compute operatorA * FEMtypeA
      If(operatorA.EQ.GRAD) Then
         Call applyGRAD(iGauss, XPL, PSI, FEMtypeA, nfa, idim, U, XYN)

      Else If(operatorA.EQ.DIV) Then
         Call applyDIV( iGauss, XPL, PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, XYN, detp)

      Else If(operatorA.EQ.IDEN) Then
         Call applyIDEN(iGauss, XPL, PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, XYN, detp)

      Else If(operatorA.EQ.CURL) Then
         Call applyCURL(iGauss, XPL, PSI, FEMtypeA, 
     &                  nfa, idim, U, XYP, detp)

      Else If(operatorA.EQ.DUDX) Then
         Call applyDUDX(iGauss, XPL, PSI, FEMtypeA, nfa, idim, U)

      Else If(operatorA.EQ.DUDY) Then
         Call applyDUDY(iGauss, XPL, PSI, FEMtypeA, nfa, idim, U)

      Else
         Call errMesFEM(2001, 'fem2Dsub', 'operatorA is not supported')
      End if

      If(nfa.GT.LDA) Call errMesFEM(2001, 'fem2Dsub',
     &     'the local matrix leading dimension, LDA, is too small')


c ... compute operatorB * FEMtypeB
      If(operatorB.EQ.GRAD) Then
         Call applyGRAD(iGauss, XQL, QSI, FEMtypeB, nfb, jdim, V, XYM)

      Else If(operatorB.EQ.DIV) Then
         Call applyDIV( iGauss, XQL, QSI, FEMtypeB, 
     &                  nfb, jdim, V, XYQ, XYM, detq)

      Else If(operatorB.EQ.IDEN) Then
         Call applyIDEN(iGauss, XQL, QSI, FEMtypeB, 
     &                  nfb, jdim, V, XYQ, XYM, detq)

      Else If(operatorB.EQ.CURL) Then
         Call applyCURL(iGauss, XQL, QSI, FEMtypeB, 
     &                  nfb, jdim, V, XYQ, detq)

      Else If(operatorB.EQ.DUDX) Then
         Call applyDUDX(iGauss, XQL, QSI, FEMtypeB, nfb, jdim, V)

      Else If(operatorB.EQ.DUDY) Then
         Call applyDUDY(iGauss, XQL, QSI, FEMtypeB, nfb, jdim, V)

      Else
         Call errMesFEM(2001, 'fem2Dsub', 'operatorB is not supported')
      End if

      If(nfb.GT.LDA) Call errMesFEM(2001, 'fem2Dsub',
     &     'the local matrix second dimension, LDA, is too small')
 

c ... compute D * U
      iD = idim
      jD = jdim

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
     &        'fem2Dsub', 'Operators A and B are not compatible')

         Do n = 1, iGauss
            Do i = 1, idim
               Do k = 1, nfa
                  DU(i, k, n) = U(i, k, n) * w(n)
               End do
            End do
         End do
      Else If(ifXtensor(tensor, TENSOR_SCALAR)) Then
         If(idim.NE.jdim) Call errMesFEM(2001, 'fem2Dsub', 
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
     &        'fem2Dsub', 'the operators A and B are not compatible')

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
         Call errMesFEM(2001, 'fem2Dsub', 
     &        'the given tensor is not supported') 
      End if


c ... compute <D U, V>
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

      nRow = nfa
      nCol = nfb

      Return
      End



C ======================================================================
      Subroutine barycentric(XY1, XY2, XY3, P1, XYP, 
     &                       iGauss, XYL, XPL, det)
C ======================================================================
C Routine transforms local (subtriangle XY1,XY2,XY3) barycentric 
C coordinates XYL to global (triangle P1,P2,P3) barycentric 
C coordinates XPL
C ======================================================================
      implicit none
      Integer  iGauss
      Real*8   XY1(2), XY2(2), XY3(2), P1(2), XYP(2, 3)
      Real*8   XYL(3, *), XPL(3, *), det

      Integer  i, j, n
      Real*8   XYsub(2, 3), L(3, 3), v, s
C ======================================================================
      Do i = 1, 2
         XYsub(i, 1) = XY1(i) - P1(i)
         XYsub(i, 2) = XY2(i) - P1(i)
         XYsub(i, 3) = XY3(i) - P1(i)
      End do

      Do i = 1, 3
         v = XYsub(1, i) * XYP(2, 2) - XYsub(2, i) * XYP(1, 2)
         L(3, i) = dabs(v / det)

         v = XYsub(1, i) * XYP(2, 3) - XYsub(2, i) * XYP(1, 3)
         L(2, i) = dabs(v / det)

         L(1, i) = 1 - L(2, i) - L(3, i)
      End do

      Do n = 1, iGauss
         Do i = 1, 3
            s = 0D0
            Do j = 1, 3
               s = s + L(i, j) * XYL(j, n)
            End do
            XPL(i, n) = s
         End do 
      End do

      Return
      End



