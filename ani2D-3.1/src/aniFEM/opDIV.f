c ======================================================================
      Subroutine applyDIV(iGauss, XYL, PSI, FEMtype, 
     &                    nfa, dim, U, XYP, XYN, det)
c ======================================================================
      Implicit None
      include 'fem2Dtri.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(2, 2), XYL(3, *), U(4, MaxDOFs, *)
      Real*8  XYP(2, 3), XYN(2, 3), det

c Local variables 
      Real*8  V(4, MaxDOFs, MaxPointGauss)
      Real*8  GRAD_P1(2, 3), GRAD_P3r(2)
      Integer IPF(3, 3)

      Integer i, k, n, iP1, iP2, nfc
      Real*8  vol, edge, edge_length, Lfun, x,y, s1,s2,s3

      DATA    GRAD_P1/-1D0,-1D0, 1D0,0D0, 0D0,1D0/
      DATA    IPF/1,2,3, 2,3,1, 1,3,2/ 

c ======================================================================
      If(FEMtype.EQ.FEM_RT0 .OR. FEMtype.EQ.FEM_BDM1) Then
         nfa = 3
         dim = 1
         vol = dabs(det)

         Do i = 1, nfa
            iP1 = IPF(1, i)
            iP2 = IPF(2, i)
            
            edge = edge_length(XYP(1, iP1), XYP(1, iP2)) 
           
            Do n = 1, iGauss
               U(1, i, n) = 2 * edge / vol 
            End do
         End do

         If(FEMtype.EQ.FEM_BDM1) Then
            nfa = 6
            Do i = 4, nfa 
               Do n = 1, iGauss
                  U(1, i, n) = 0D0
               End do
            End do
         End if

c ... next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         Call applyGRAD(1, XYL, PSI, FEM_P1, nfc, dim, V, XYN)

         nfa = 6
         dim = 1 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            U(1, i,     1) = V(1, i, 1)
            U(1, i + 3, 1) = V(2, i, 1)
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1vector) Then
         Call applyGRAD(1, XYL, PSI, FEM_CR1, nfc, dim, V, XYN)

         nfa = 6
         dim = 1 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            U(1, i,     1) = V(1, i, 1)
            U(1, i + 3, 1) = V(2, i, 1)
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_MINI) Then
         Call applyGRAD(1, XYL, PSI, FEM_P1, nfc, dim, V, XYN)

         nfa = 8
         dim = 1 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            U(1, i,     1) = V(1, i, 1)
            U(1, i + 4, 1) = V(2, i, 1)
         End do
         Call copyGauss(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            s1 = Lfun(1, x, y)
            s2 = Lfun(2, x, y)
            s3 = Lfun(3, x, y)

            Do k = 1, 2
               GRAD_P3r(k) = 27 * (GRAD_P1(k, 1) * s2 * s3 +
     &                             GRAD_P1(k, 2) * s1 * s3 +
     &                             GRAD_P1(k, 3) * s1 * s2)
            End do

            Do k = 1, 2
               U(1, 4*k, n) = PSI(1, k) * GRAD_P3r(1) 
     &                      + PSI(2, k) * GRAD_P3r(2)
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P2reduced) Then
         Call applyGRAD(1, XYL, PSI, FEM_P1, nfc, k, V, XYN)

         nfa = 9
         dim = 1 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            U(1, i,     1) = V(1, i, 1) 
            U(1, i + 6, 1) = V(2, i, 1) 
         End do
         Call copyGauss(iGauss, nfa, dim, U)

         Call applyGRAD(iGauss, XYL, PSI, FEM_P2, nfc, k, V, XYN)

         Do n = 1, iGauss
            Do i = 1, 3
               U(1, i + 3, n) = XYN(1, i) * V(1, i + 3, n)
     &                        + XYN(2, i) * V(2, i + 3, n)
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P2vector) Then
         Call applyGRAD(iGauss, XYL, PSI, FEM_P2, nfc, k, V, XYN)

         nfa = 12
         dim = 1
         Call clearU(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            Do i = 1, 6
               U(1, i,     n) = V(1, i, n)
               U(1, i + 6, n) = V(2, i, n)
            End do
         End do

      Else
         nfa = 0
         dim = -FEMtype
         Call errMesFEM(2001, 'fem2Dtri', 
     &        'unsupported operation for the given element type')
      End if
      Return
      End
