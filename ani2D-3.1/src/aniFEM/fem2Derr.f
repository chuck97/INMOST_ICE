C ======================================================================
      Subroutine fem2Derr(XY1, XY2, XY3, Lp,
     &                    operatorA, FEMtypeA, Uh, Fu, dDATAFU, iDATAFU,
     &                    label, D, dDATA, iDATA, iSYS, order, ERR)
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'
C ======================================================================
C  *** ITRODUCTION ***
C  The routine computes the weighted norm of error:  
C
C        Integral (D (OpA(u) - Fu) (OpA(u) - Fu))^{Lp/2}
C
C  where D is a tensor, OpA is a linear operator, u is the discrete
C  solution and Fu is the target (analytic or exact) function.
C
C  Details of input the parameters are described in file fem2Dtri.f. 
C  Here we focus on additional parameters, Lp, Uh and Fu, only.
C 
C  Lp - norm for which the error has to be calculated:
C             Lp > 0  means the L^p norm
C             Lp = 0  means the maximum norm (L^infinity)
C
C  Uh - Real*8 vector of discrete solution over triangle. It must have 
C       the same size as the matrix of bilinear form <OpA(u), OpA(u)>
C
C  Fu - Integer function for the exact solution. The standard format:
C       Fu(x, y, label, dDATAFU, iDATAFU, iSYS, Diff)
C
C       The function returns type of the tendor Diff which is the
C       value of this function. See fem2Dtri.f for more details.
C
C ======================================================================
      Real*8   XY1(*), XY2(*), XY3(*), Lp
      Integer  operatorA, FEMtypeA, label, order, D, Fu

      Real*8   dDATA(*), dDATAFU(*), Uh(*), ERR
      Integer  iDATA(*), iDATAFU(*), iSYS(*)

      EXTERNAL D, Fu

C ======================================================================
c Local variables
      Real*8  XYP(2, 3), XYN(2, 3)
      Real*8  det, vol, s, p
   
      Real*8  PSI(2, 2), U(4, MaxDOFs, MaxPointGauss)
      Real*8  Diff(4, 4, MaxPointGauss), UUh(4, MaxPointGauss)
      Real*8  V(4, 4, MaxPointGauss),   DUUh(4, MaxPointGauss)

      Real*8  w(MaxPointGauss), XYG(2, MaxPointGauss)
      Real*8  XYL(3, AllPointGauss)

      Integer i,j,k,n, iD,jD, nfa, idim, iL, iGauss, tensor
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
     &  Call errMesFEM(2001, 'fem2Derr', 'quadrature order is wrong')

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
      If(FEMtypeA.EQ.FEM_P2reduced) Then
         Call calNormalExt(XY1, XY2, XY3, XYN(1, 1))
         Call calNormalExt(XY2, XY3, XY1, XYN(1, 2))
         Call calNormalExt(XY3, XY1, XY2, XYN(1, 3))
      End if

c ... compute operatorA * FEMtypeA
      If(operatorA.EQ.GRAD) Then
         Call applyGRAD(iGauss, XYL(1, iL), PSI, FEMtypeA, 
     &                  nfa, idim, U, XYN)

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
         Call errMesFEM(2001, 'fem2Derr', 'operatorA is not supported')
      End if


c ... compute U * Uh and save in UUh
      Do n = 1, iGauss
         Do i = 1, idim
            UUh(i, n) = 0D0
            Do k = 1, nfa
               UUh(i, n) = UUh(i, n) + U(i, k, n) * Uh(k)
            End do
         End do
      End do


c ... compute values of Fu at quadrature points
      Do n = 1, iGauss
         k = Fu(XYG(1,n), XYG(2,n), label, 
     &          dDATAFU, iDATAFU, iSYS, V(1,1,n))

         If(iSYS(1).ne.idim) Call errMesFEM(2001, 'fem2Derr',
     &      'Analytic function must return a column vector, iSYS(2)=1')

         Do i = 1, idim
            UUh(i, n) = UUh(i, n) - V(i, 1, n)
         End do
      End do


c ... compute values of D * UUh and save it in DUUh
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
         Do n = 1, iGauss
            Do i = 1, idim
               DUUh(i, n) = UUh(i, n)
            End do
         End do

      Else If(ifXtensor(tensor, TENSOR_SCALAR)) Then
         Do n = 1, iGauss
            Do i = 1, idim
               DUUh(i, n) = UUh(i, n) * Diff(1,1,n)
            End do
         End do

      Else If(ifXtensor(tensor, TENSOR_SYMMETRIC) .OR. 
     &        ifXtensor(tensor, TENSOR_GENERAL)) Then
         iD = iSYS(1)
         jD = iSYS(2)

         Do n = 1, iGauss
            Do i = 1, iD
               s = 0D0
               Do j = 1, jD
                  s = s + Diff(i, j, n) * UUh(j, n)
               End do
               DUUh(i, n) = s
            End do
         End do
      Else
         Call errMesFEM(2001, 'fem2Derr', 
     &        'the given tensor is not supported') 
      End if


c ... compute |D (OpA(u) - Fu) (OpA(u) - Fu)|^p
c ... for maximum norm (p=0), we calculate norm of 2-D vector
      p = Lp / 2

      ERR = 0D0
      Do n = 1, iGauss
         s = 0D0
         Do k = 1, idim
            s = s + DUUh(k, n) * UUh(k, n)
         End do

         If(p.EQ.0D0) Then
            ERR = max(ERR, dsqrt(s))
         Else
            ERR = ERR + (s**p) * w(n)
         End if
      End do

      Return
      End
