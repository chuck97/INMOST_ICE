c ======================================================================
      Subroutine applyCURL(iGauss, XYL, PSI, FEMtype, 
     &                     nfa, dim, U, XYP, det)
c ======================================================================
      Implicit None
      include 'fem2Dtri.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(2, 2), XYL(3, *), U(4, MaxDOFs, *)
      Real*8  XYP(2, *), det

      Integer i, j, k, l, n, mfa, i1,i2,i3
      Real*8  edge_length, edge, Lfun, x, y, s1, s2, s3, grad_xy(2)

      Real*8  GRAD_P1(2, 3), GRAD_P2(2,  6, MaxPointGauss)
      Real*8                 GRAD_P3(2, 10, MaxPointGauss), GRAD_P3r(2)

      Integer iref(4)
      DATA    iref /1,2,3,1/
      DATA    GRAD_P1/-1D0,-1D0, 1D0,0D0, 0D0,1D0/

c ======================================================================
      If(FEMtype.EQ.FEM_P0) Then
         nfa = 1
         dim = 1
         Do n = 1, iGauss
            U(1, 1, n) = 0D0
         End do

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_P1) Then
         nfa = 3
         dim = 2 
         Do i = 1, nfa
            Do k = 1, dim
               U(k, i, 1) = PSI(1, 3-k) * GRAD_P1(1, i) 
     &                    + PSI(2, 3-k) * GRAD_P1(2, i)
            End do
            U(2, i, 1) = -U(2, i, 1)
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_CR1) Then
         nfa = 3
         dim = 2 
         Do i = 1, nfa
            l = iref(i + 1)
            Do k = 1, dim
               U(k, l, 1) = PSI(1, 3-k) * GRAD_P1(1, i)
     &                    + PSI(2, 3-k) * GRAD_P1(2, i)
            End do
            U(2, l, 1) = -U(2, l, 1)
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_RT0 .OR. FEMtype.EQ.FEM_BDM1) Then
         nfa = 3
         dim = 1

         Do i1 = 1, nfa
            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)
            
            edge = edge_length(XYP(1, i1), XYP(1, i2)) 
            edge = edge / (det**2)

            Do n = 1, iGauss
               U(1, i1, n) = 0D0

               If(FEMtype.EQ.FEM_BDM1) Then
                  U(1, i1+nfa, n) = edge * 12 *
     &               ((XYP(1,i1) - XYP(1,i3)) * (XYP(1,i2) - XYP(1,i3))
     &               +(XYP(2,i1) - XYP(2,i3)) * (XYP(2,i2) - XYP(2,i3)))
               End if
            End do
         End do

         If(FEMtype.EQ.FEM_BDM1) nfa = 6

c ... next FEM
      Else If(FEMtype.EQ.FEM_P2) Then
         nfa = 6
         dim = 2

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               Do k = 1, dim
                  GRAD_P2(k, i, n) = GRAD_P1(k, i) *
     &                               (4 * Lfun(i, x, y) - 1D0) 
               End do
            End do

            mfa = 3
            Do i = 1, 3
               j = iref(i + 1)
               mfa = mfa + 1
               Do k = 1, dim
                  GRAD_P2(k, mfa, n) = 
     &               4 * (Lfun(i, x, y) * GRAD_P1(k, j) +
     &                    Lfun(j, x, y) * GRAD_P1(k, i))
               End do
            End do
         End do

         Do n = 1, iGauss
            Do i = 1, nfa
               Do k = 1, dim
                  U(k, i, n) = 0D0
                  Do j = 1, 2
                     U(k, i, n) = U(k, i, n) 
     &                          + PSI(j, 3-k) * GRAD_P2(j, i, n)
                  End do
               End do

               U(2, i, n) = -U(2, i, n)
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P3) Then
         nfa = 10
         dim = 2

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)

               s1 = Lfun(i, x, y)
               s2 = Lfun(j, x, y)

               Do k = 1, dim
                  GRAD_P3(k, i, n) = 
     &                GRAD_P1(k, i) * (13.5D0 * s1 * (s1 - 1D0) + 1D0)

                  GRAD_P3(k, i + 3, n) = 
     &                GRAD_P1(k, i) * (6D0 * s1 - 1D0) * s2 +
     &                GRAD_P1(k, j) * (3D0 * s1 - 1D0) * s1

                  GRAD_P3(k, i + 6, n) = 
     &                GRAD_P1(k, j) * (6D0 * s2 - 1D0) * s1 +
     &                GRAD_P1(k, i) * (3D0 * s2 - 1D0) * s2
               End do
            End do

            s1 = Lfun(1, x, y)
            s2 = Lfun(2, x, y)
            s3 = Lfun(3, x, y)

            Do k = 1, dim
               GRAD_P3(k, 10, n) = 27 * 
     &            (GRAD_P1(k, 1) * s2 * s3 +
     &             GRAD_P1(k, 2) * s1 * s3 +
     &             GRAD_P1(k, 3) * s1 * s2)
            End do
         End do

         Do n = 1, iGauss
            Do i = 1, nfa
               Do k = 1, dim
                  U(k, i, n) = 0D0
                  Do j = 1, 2
                     U(k, i, n) = U(k, i, n) 
     &                          + PSI(j, 3-k) * GRAD_P3(j, i, n)
                  End do
               End do

               U(2, i, n) = -U(2, i, n)
            End do
         End do
 
c ... next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         nfa = 6
         dim = 1
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            Do k = 1, 2
               grad_xy(k) = 0D0
               Do j = 1, 2
                  grad_xy(k) = grad_xy(k) + PSI(j, k) * GRAD_P1(j, i)
               End do

               U(1, i,     1) =  grad_xy(2)
               U(1, i + 3, 1) = -grad_xy(1)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_MINI) Then
         nfa = 8
         dim = 1
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            Do k = 1, 2
               grad_xy(k) = 0D0
               Do j = 1, 2
                  grad_xy(k) = grad_xy(k) + PSI(j, k) * GRAD_P1(j, i)
               End do

               U(1, i,     1) =  grad_xy(2)
               U(1, i + 4, 1) = -grad_xy(1)
            End do
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
               grad_xy(k) = PSI(1, k) * GRAD_P3r(1) 
     &                    + PSI(2, k) * GRAD_P3r(2)
            End do

            U(1, 4, n) =  grad_xy(2)
            U(1, 8, n) = -grad_xy(1)
         End do

      Else
         nfa = 0
         dim = -FEMtype
         Call errMesFEM(2001, 'fem2Dtri', 
     &        'unsupported operation for the given element type')
      End if
      Return
      End
