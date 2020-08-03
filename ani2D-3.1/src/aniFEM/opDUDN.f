c ======================================================================
      Subroutine applyDUDN(iGauss, XYL, PSI, FEMtype, nfa, dim, U, XYN)
c ======================================================================
c Calculates normal derivative of basis functions at points XYL(2,*) 
c on the 1st edge of considered triangle.
c ======================================================================
      implicit none
      include 'fem2Dtri.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(2, 2), U(2, MaxDOFs, *), XYN(2)

c Data for the reference triangle
      Real*8  XYL(2, *)
      Real*8  GRAD_P1(2, 3), GRAD_PX(2, MaxDOFs, MaxPnt2DGauss)

      Integer iref(4), i, j, k, l, n, mfa
      Real*8  Lfun, s1, s2, s3, x, s

      DATA    GRAD_P1/-1D0,-1D0, 1D0,0D0, 0D0,1D0/
      DATA    iref /1,2,3,1/

c ======================================================================
      If(FEMtype.EQ.FEM_P0) Then
         nfa = 1
         dim = 1
         Do n = 1, iGauss
            U(1, 1, n) = 0D0
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P1) Then
         nfa = 3
         dim = 1 
         Do i = 1, nfa
            U(1, i, 1) = 0D0
            Do k = 1, 2
               s = PSI(1, k) * GRAD_P1(1, i) 
     &           + PSI(2, k) * GRAD_P1(2, i)
               U(1, i, 1) = U(1, i, 1) + s * XYN(k)   
            End do
         End do

         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         nfa = 6
         dim = 2 
         Call clearU(1, nfa, dim, U)

         Do i = 1, 3
            U(1, i, 1) = 0D0
            Do k = 1, 2
               s = PSI(1, k) * GRAD_P1(1, i)
     &           + PSI(2, k) * GRAD_P1(2, i)
               U(1, i, 1) = U(1, i, 1) + s * XYN(k)
            End do
            U(2, i + 3, 1) = U(1, i, 1)
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_P2) Then
         nfa = 6
         dim = 1

         Do n = 1, iGauss
            x = XYL(2, n)
            Do i = 1, 3
               Do k = 1, 2
                  GRAD_PX(k, i, n) = GRAD_P1(k, i) *
     &                               (4 * Lfun(i, x, 0D0) - 1D0) 
               End do
            End do

            mfa = 3
            Do i = 1, 3
               j = iref(i + 1)
               mfa = mfa + 1
               Do k = 1, 2
                  GRAD_PX(k, mfa, n) = 
     &               4 * (Lfun(i, x, 0D0) * GRAD_P1(k, j) +
     &                    Lfun(j, x, 0D0) * GRAD_P1(k, i))
               End do
            End do
         End do

         Do n = 1, iGauss
            Do i = 1, nfa
               U(1, i, n) = 0D0
               Do k = 1, 2
                  s = PSI(1, k) * GRAD_PX(1, i, n) 
     &              + PSI(2, k) * GRAD_PX(2, i, n)
                  U(1, i, n) = U(1, i, n) + s * XYN(k)   
               End do
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P3) Then
         nfa = 10
         dim = 1

         Do n = 1, iGauss
            x = XYL(2, n)

            Do i = 1, 3
               j = iref(i + 1)

               s1 = Lfun(i, x, 0D0)
               s2 = Lfun(j, x, 0D0)

               Do k = 1, 2
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

            s1 = Lfun(1, x, 0D0)
            s2 = Lfun(2, x, 0D0)
            s3 = Lfun(3, x, 0D0)

            Do k = 1, 2
               GRAD_PX(k, 10, n) = 27 * 
     &            (GRAD_P1(k, 1) * s2 * s3 +
     &             GRAD_P1(k, 2) * s1 * s3 +
     &             GRAD_P1(k, 3) * s1 * s2)
            End do
         End do

         Do n = 1, iGauss
            Do i = 1, nfa
               U(1, i, n) = 0D0
               Do k = 1, 2
                  s = PSI(1, k) * GRAD_PX(1, i, n) 
     &              + PSI(2, k) * GRAD_PX(2, i, n)
                  U(1, i, n) = U(1, i, n) + s * XYN(k)   
               End do
            End do
         End do

c ... next FEM
      Else If(FEMtype.EQ.FEM_P4) Then
         nfa = 15
         dim = 1

         Do n = 1, iGauss
            x = XYL(2, n)

            Do i = 1, 3
               j = iref(i + 1)
               l = iref(j + 1)

               s1 = Lfun(i, x, 0D0)
               s2 = Lfun(j, x, 0D0)
               s3 = Lfun(l, x, 0D0)

               Do k = 1, 2
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
               U(1, i, n) = 0D0
               Do k = 1, 2
                  s = PSI(1, k) * GRAD_PX(1, i, n) 
     &              + PSI(2, k) * GRAD_PX(2, i, n)
                  U(1, i, n) = U(1, i, n) + s * XYN(k)   
               End do
            End do
         End do
      End if

      Return
      End

