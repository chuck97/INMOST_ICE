c ======================================================================
      Subroutine applyGRAD(iGauss, XYL, PSI, FEMtype, 
     &                     nfa, dim, U, XYP, XYN)
c ======================================================================
      implicit none
      include 'fem2Dtri.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(2, 2), U(4, MaxDOFs, *), XYP(2, 3), XYN(2, 3)

c Data for the reference triangle
      Real*8  XYL(3, *)
      Real*8  GRAD_P1(2, 3), GRAD_PX(2, MaxDOFs, MaxPointGauss)
      Real*8  GRAD_P2r(2), GRAD_P3r(2)

      Integer i, j, k, n, l, mfa
      Real*8  x, y, Lfun, s1, s2, s3, gx, gy, dx, dy, ax, ay, a, b, c

      Integer iref(4)

      DATA    GRAD_P1/-1D0,-1D0, 1D0,0D0, 0D0,1D0/
      DATA    iref /1,2,3,1/

c ======================================================================
      If(FEMtype.EQ.FEM_P0) Then
         nfa = 1
         dim = 1
         Do n = 1, iGauss
            U(1, 1, n) = 0D0
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_P1) Then
         nfa = 3
         dim = 2 
         Do i = 1, nfa
            Do k = 1, dim
               U(k, i, 1) = 0D0
               Do j = 1, 2
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c === next FEM
      Else If(FEMtype.EQ.FEM_CR1) Then
c  ...  use formula psi_l = 1 - 2 phi_i
         nfa = 3
         dim = 2 
         Do i = 1, nfa
            l = iref(i + 1)
            Do k = 1, dim
               U(k, l, 1) = 0D0
               Do j = 1, 2
                  U(k, l, 1) = U(k, l, 1) - 2*PSI(j, k) * GRAD_P1(j, i)
               End do
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c === next FEM
      Else If(FEMtype.EQ.FEM_CR1vector) Then
         nfa = 6
         dim = 4
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            l = iref(i + 1)
            Do k = 1, 2
               U(k, l, 1) = 0D0
               Do j = 1, 2
                  U(k, l, 1) = U(k, l, 1) - 2*PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 2, l + 3, 1) = U(k, l, 1)
            End do
         End do

         Call copyGauss(iGauss, nfa, dim, U)

c === next three FEMs
      Else If(FEMtype.EQ.FEM_P2 .OR. FEMtype.EQ.FEM_P2vector) Then
         nfa = 6
         dim = 2
         If(FEMtype.EQ.FEM_P2vector)  Call clearU(iGauss, 12, 4, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            Do i = 1, 3
               Do k = 1, dim
                  GRAD_PX(k, i, n) = GRAD_P1(k, i) *
     &                               (4 * Lfun(i, x, y) - 1D0) 
               End do
            End do

            mfa = 3
            Do i = 1, 3
               j = iref(i + 1)
               mfa = mfa + 1
               Do k = 1, dim
                  GRAD_PX(k, mfa, n) = 
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
     &                          + PSI(j, k) * GRAD_PX(j, i, n)
                  End do
               End do
            End do
         End do

         If(FEMtype.EQ.FEM_P2vector) Then
            nfa = 12
            dim = 4

            Do n = 1, iGauss
               Do i = 1, 6
                  Do k = 1, 2
                     U(k + 2, i + 6, n) = U(k, i, n)
                  End do
               End do
            End do
         End if

c === next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         nfa = 6
         dim = 4 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            Do k = 1, 2
               Do j = 1, 2
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 2, i + 3, 1) = U(k, i, 1)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c === next FEM
      Else If(FEMtype.EQ.FEM_MINI) Then
         nfa = 8
         dim = 4 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            Do k = 1, 2
               Do j = 1, 2
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 2, i + 4, 1) = U(k, i, 1)
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
               U(k, 4, n) = PSI(1, k) * GRAD_P3r(1) 
     &                    + PSI(2, k) * GRAD_P3r(2)
               U(k + 2, 8, n) = U(k, 4, n)
            End do
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_P2reduced) Then
         nfa = 9
         dim = 4 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            Do k = 1, 2
               Do j = 1, 2
                  U(k, i, 1) = U(k, i, 1) + PSI(j, k) * GRAD_P1(j, i)
               End do
               U(k + 2, i + 6, 1) = U(k, i, 1)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)

               Do k = 1, 2
                  GRAD_P2r(k) = 4 * (Lfun(i, x, y) * GRAD_P1(k, j) +
     &                               Lfun(j, x, y) * GRAD_P1(k, i))
               End do

               Do k = 1, 2
                  s1 = PSI(1, k) * GRAD_P2r(1) + PSI(2, k) * GRAD_P2r(2)
                  U(k,     i + 3, n) = XYN(1, i) * s1
                  U(k + 2, i + 3, n) = XYN(2, i) * s1
               End do
            End do
         End do

c === next FEM
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
                  GRAD_PX(k, i, n) = 
     &                GRAD_P1(k, i) * (4.5D0 * s1 * (3*s1-2) + 1D0)

                  GRAD_PX(k, i + 3, n) = 4.5D0 *
     &               (GRAD_P1(k, i) * (6*s1-1) * s2 +
     &                GRAD_P1(k, j) * (3*s1-1) * s1)

                  GRAD_PX(k, i + 6, n) = 4.5D0 *
     &               (GRAD_P1(k, j) * (6*s2-1) * s1 +
     &                GRAD_P1(k, i) * (3*s2-1) * s2)
               End do
            End do

            s1 = Lfun(1, x, y)
            s2 = Lfun(2, x, y)
            s3 = Lfun(3, x, y)

            Do k = 1, dim
               GRAD_PX(k, 10, n) = 27 * 
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
     &                          + PSI(j, k) * GRAD_PX(j, i, n)
                  End do
               End do
            End do
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_P4) Then
         nfa = 15
         dim = 2

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)
               l = iref(j + 1)

               s1 = Lfun(i, x, y)
               s2 = Lfun(j, x, y)
               s3 = Lfun(l, x, y)

               Do k = 1, dim
                  GRAD_PX(k, i, n) = 
     &                GRAD_P1(k, i) * (s1*s1*(128*s1-144) + 44*s1-3) / 3

                  GRAD_PX(k, i + 3, n) = 8D0 / 3 *
     &               (GRAD_P1(k, i) * (s1*(48*s1-24) + 2) * s2 +
     &                GRAD_P1(k, j) *  s1*(4*s1-1)*(4*s1-2))

                  GRAD_PX(k, i + 6, n) = 4 *
     &               (GRAD_P1(k, i) * (8*s1-1)*s2*(4*s2-1) +
     &                GRAD_P1(k, j) * (8*s2-1)*s1*(4*s1-1))

                  GRAD_PX(k, i + 9, n) = 8D0 / 3 *
     &               (GRAD_P1(k, j) * (s2*(48*s2-24) + 2) * s1 +
     &                GRAD_P1(k, i) *  s2*(4*s2-1)*(4*s2-2))

                  GRAD_PX(k, i + 12, n) = 32 * 
     &               (GRAD_P1(k, i) * (8*s1-1) * s2 * s3 +
     &                GRAD_P1(k, j) * s1 * (4*s1-1) * s3 +
     &                GRAD_P1(k, l) * s1 * (4*s1-1) * s2)
               End do
            End do
         End do

         Do n = 1, iGauss
            Do i = 1, nfa
               Do k = 1, dim
                  U(k, i, n) = 0D0
                  Do j = 1, 2
                     U(k, i, n) = U(k, i, n) 
     &                          + PSI(j, k) * GRAD_PX(j, i, n)
                  End do
               End do
            End do
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_HERMIT) Then
         nfa = 10 
         dim = 2

         Do i = 1, 3
            Do k = 1, 2
               GRAD_PX(k, i, 1) = PSI(1, k) * GRAD_P1(1, i) 
     &                          + PSI(2, k) * GRAD_P1(2, i)
            End do
         End do

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)
               k = iref(j + 1)

               gx = GRAD_PX(1, i, 1)
               gy = GRAD_PX(2, i, 1)

               s1 = Lfun(i, x, y)
               s2 = Lfun(j, x, y)
               s3 = Lfun(k, x, y)

               c = s1 * s1
               a = s1 * s2 * s3
               ax = gx * s2 * s3 + GRAD_PX(1, j, 1) * s1 * s3 
     &                           + GRAD_PX(1, k, 1) * s1 * s2
               ay = gy * s2 * s3 + GRAD_PX(2, j, 1) * s1 * s3 
     &                           + GRAD_PX(2, k, 1) * s1 * s2

               dx = x * XYP(1, 2) + y * XYP(1, 3) - XYP(1, i)
               dy = x * XYP(2, 2) + y * XYP(2, 3) - XYP(2, i)

               b = 1 - 2 * gx * dx - 2 * gy * dy - s1
               U(1, i, n) = 2 * s1 * gx * b - 7 * ax
               U(2, i, n) = 2 * s1 * gy * b - 7 * ay

               U(1, i+3, n) = 2 * s1 * gx * dx + s1 * s1 - 3 * dx * ax
               U(2, i+3, n) = 2 * s1 * gy * dy + s1 * s1 - 3 * dx * ay

               U(1, i+6, n) = 2 * s1 * gx * dy + s1 * s1 - 3 * dy * ax
               U(2, i+6, n) = 2 * s1 * gy * dy + s1 * s1 - 3 * dy * ay
            End do

            U(1, 10, n) = 27 * ax
            U(2, 10, n) = 27 * ay
         End do

      Else
         nfa = 0
         dim = -FEMtype
         Call errMesFEM(2001, 'fem2Dtri', 
     &        'unsupported operation for this element type')
      End if
      Return
      End
