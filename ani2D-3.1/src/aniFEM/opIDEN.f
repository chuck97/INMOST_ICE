c ======================================================================
      Subroutine applyIDEN(iGauss, XYL, PSI, FEMtype, 
     &                     nfa, dim, U, XYP, XYN, det)
c ======================================================================
      implicit none
      include 'fem2Dtri.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(2, 2), XYL(3, *), U(4, MaxDOFs, *)
      Real*8  XYP(2, *), XYN(2, *), det

c Local variables 
      Integer i, j, k, l, n, i1, k1, mfa, iP1, iP2, iP3
      Real*8  vol, edge, edge_length
      Real*8  x, y, Lfun, s, s1, s2, s3
      Real*8  GRAD_LX(2, 3), DN(3, 3) 
      Real*8  gx, gy, dx, dy, dx2, dxy, dy2, ex, ey, exy, a, b, c

c Data for the reference triangle
      Real*8  GRAD_P1(2, 3)
      Integer IPF(3, 3), iref(4)

      DATA    GRAD_P1/-1D0,-1D0, 1D0,0D0, 0D0,1D0/
      DATA    IPF/1,2,3, 2,3,1, 1,3,2/ 
      DATA    iref /1,2,3,1/

c ======================================================================
      If(FEMtype.EQ.FEM_P0) Then
         nfa = 1
         dim = 1
         Do n = 1, iGauss
            U(1, 1, n) = 1D0
         End do

c === next two FEMs
      Else If(FEMtype.EQ.FEM_P1) Then
         nfa = 3
         dim = 1
         Do n = 1, iGauss
            Do i = 1, nfa
               U(1, i, n) = XYL(i, n)
            End do
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_CR1) Then
c ... use formula psi_l = 1 - 2 phi_i
         nfa = 3
         dim = 1
         Do n = 1, iGauss
            Do i = 1, nfa
               l = iref(i + 1)
               U(1, l, n) = 1 - 2 * XYL(i, n)  
            End do
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_CR1vector) Then
         nfa = 6
         dim = 2
         Call clearU(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            Do i = 1, nfa
               l = iref(i + 1)
               U(1, l,     n) = 1 - 2 * XYL(i, n)  
               U(2, l + 3, n) = 1 - 2 * XYL(i, n)  
            End do
         End do

c === next two FEMs
      Else If(FEMtype.EQ.FEM_P2 .OR. FEMtype.EQ.FEM_P2vector) Then
         nfa = 6
         dim = 1
         If(FEMtype.EQ.FEM_P2vector) Call clearU(iGauss, 12, 2, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            Do i = 1, 3
               s = Lfun(i, x, y)
               U(1, i, n) = (2 * s - 1D0) * s
            End do

            mfa = 3
            Do i = 1, 3
               j = iref(i + 1)
               mfa = mfa + 1
               U(1, mfa, n) = 4 * Lfun(i, x, y) * Lfun(j, x, y)
            End do
         End do

         If(FEMtype.EQ.FEM_P2vector) Then
            nfa = 12
            dim = 2

            Do n = 1, iGauss
               Do i = 1, 6
                  U(2, i + 6, n) = U(1, i, n)
               End do
            End do
         End if

c === next FEM
      Else If(FEMtype.EQ.FEM_P3) Then
         nfa = 10 
         dim = 1

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)

               s1 = Lfun(i, x, y)
               s2 = Lfun(j, x, y)

               U(1, i,     n) = s1 * (3*s1 - 1) * (3*s1 - 2) / 2
               U(1, i + 3, n) = s1 * (3*s1 - 1) * s2 * 4.5D0
               U(1, i + 6, n) = s1 * (3*s2 - 1) * s2 * 4.5D0
            End do

            U(1, 10, n) = 27 * Lfun(1, x, y) 
     &                       * Lfun(2, x, y) * Lfun(3, x, y)
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_P4) Then
         nfa = 15 
         dim = 1

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)
               k = iref(j + 1)

               s1 = Lfun(i, x, y)
               s2 = Lfun(j, x, y)
               s3 = Lfun(k, x, y)

               U(1, i,      n) = s1 * (4*s1-1) * (4*s1-2) * (4*s1-3) / 6
               U(1, i +  3, n) = s1 * (4*s1-1) * (4*s1-2) * s2 * 8 / 3D0
               U(1, i +  6, n) = s1 * (4*s1-1) * s2 * (4*s2-1) * 4 
               U(1, i +  9, n) = s1 * s2 * (4*s2-1) * (4*s2-2) * 8 / 3D0
               U(1, i + 12, n) = s1 * (4*s1-1) * s2 * s3 * 32
            End do
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_P1vector) Then
         nfa = 6
         dim = 2

         Call clearU(iGauss, nfa, dim, U)
         Do n = 1, iGauss
            Do i = 1, 3
               U(1, i,     n) = XYL(i, n)
               U(2, i + 3, n) = XYL(i, n)
            End do
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_MINI) Then
         nfa = 8
         dim = 2

         Call clearU(iGauss, nfa, dim, U)
         Do n = 1, iGauss
            Do i = 1, 3
               U(1, i,     n) = XYL(i, n)
               U(2, i + 4, n) = XYL(i, n)
            End do
         End do

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            U(1, 4, n) = 27 * Lfun(1, x, y) 
     &                      * Lfun(2, x, y) * Lfun(3, x, y)
            U(2, 8, n) = U(1, 4, n) 
         End do

c === next FEM
      Else If(FEMtype.EQ.FEM_P2reduced) Then
         nfa = 9
         dim = 2

         Call clearU(iGauss, nfa, dim, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               U(1, i,     n) = XYL(i, n)
               U(2, i + 6, n) = XYL(i, n)
            End do

            Do i = 1, 3
               j = iref(i + 1)

               s1 = 4 * Lfun(i, x, y) * Lfun(j, x, y)
               U(1, i + 3, n) = s1 * XYN(1, i)
               U(2, i + 3, n) = s1 * XYN(2, i)
            End do
         End do

c === next two FEMs
      Else If(FEMtype.EQ.FEM_RT0 .OR. FEMtype.EQ.FEM_BDM1) Then
         nfa = 3
         dim = 2
         vol = dabs(det)

         Do i = 1, nfa
            iP1 = IPF(1, i)
            iP2 = IPF(2, i)
            iP3 = IPF(3, i)
            
            edge = edge_length(XYP(1, iP1), XYP(1, iP2)) 
            edge = edge / vol
            Do n = 1, iGauss
               Do k = 1, dim
                  U(k, i, n) = edge *
     &                 (XYL(iP1, n) * (XYP(k, iP1) - XYP(k, iP3))
     &                + XYL(iP2, n) * (XYP(k, iP2) - XYP(k, iP3)))
               End do

               If(FEMtype.EQ.FEM_BDM1) Then
                  Do k = 1, dim
                     U(k, i+nfa, n) = edge * 6 *
     &                   (- XYL(iP1, n) * (XYP(k, iP1) - XYP(k, iP3))
     &                    + XYL(iP2, n) * (XYP(k, iP2) - XYP(k, iP3)))
                  End do
               End if
            End do
         End do

         If(FEMtype.EQ.FEM_BDM1) nfa = 6

c === next FEM
      Else If(FEMtype.EQ.FEM_HERMIT) Then
         nfa = 10 
         dim = 1

         Do i = 1, 3
            Do k = 1, 2
               GRAD_LX(k, i) = PSI(1, k) * GRAD_P1(1, i) 
     &                       + PSI(2, k) * GRAD_P1(2, i)
            End do
         End do

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)
               k = iref(j + 1)

               gx = GRAD_LX(1, i)
               gy = GRAD_LX(2, i)

               s1 = Lfun(i, x, y)
               s2 = Lfun(j, x, y)
               s3 = Lfun(k, x, y)

               c = s1 * s1
               a = s1 * s2 * s3

               dx = x * XYP(1, 2) + y * XYP(1, 3) - XYP(1, i)
               dy = x * XYP(2, 2) + y * XYP(2, 3) - XYP(2, i)
               U(1, i, n) = c * (1 - 2 * gx * dx - 2 * gy * dy) - 7 * a
               U(1, i+3, n) = c * dx - 3 * dx * a
               U(1, i+6, n) = c * dy - 3 * dy * a 
            End do

            U(1, 10, n) = 27 * a
         End do


c === next FEM
      Else If(FEMtype.EQ.FEM_ARGYRIS) Then
         nfa = 21
         dim = 1

         Do i = 1, 3
            Do k = 1, 2
               GRAD_LX(k, i) = PSI(1, k) * GRAD_P1(1, i) 
     &                       + PSI(2, k) * GRAD_P1(2, i)
            End do
         End do

c ... calculate normal derivatives on edges (j=edge, i=vertex) 
         Do i = 1, 3
            Do j = 1, 3
               DN(j, i) = GRAD_LX(1, i) * XYN(1, j)  
     &                  + GRAD_LX(2, i) * XYN(2, j) 
            End do
         End do

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)

            Do i = 1, 3
               j = iref(i + 1)
               k = iref(j + 1)

               gx = GRAD_LX(1, i)
               gy = GRAD_LX(2, i)

               dx = x * XYP(1, 2) + y * XYP(1, 3) - XYP(1, i)
               dy = x * XYP(2, 2) + y * XYP(2, 3) - XYP(2, i)
               dx2 = dx * dx
               dxy = dx * dy
               dy2 = dy * dy

               s1 = Lfun(i, x, y)
               s2 = Lfun(j, x, y)
               s3 = Lfun(k, x, y)

               a = gx*gx * dx2 + 2 * gx*gy * dxy + gy*gy * dy2
               b = gx * dx + gy * dy
               c = s1 * s1 * s1
               U(1, i, n) = c * (6 * a - 3 * b + 1)

               U(1, i+ 3, n) = c * (-3 * (gx * dx2 + gy * dxy) + dx)
               U(1, i+ 6, n) = c * (-3 * (gx * dxy + gy * dy2) + dy)

               U(1, i+ 9, n) = c * dx2 / 2
               U(1, i+12, n) = c * dy2 / 2
               U(1, i+15, n) = c * dxy 

               U(1, i+18, n) = 16 * s1 * s1 * s2 * s2 * s3 / DN(i, k)
            End do
         End do

c ... normalize normal derivatives on edges
         Do i = 1, 3
            j = iref(i + 1)
            k = iref(j + 1)

            i1 = i + 18
            k1 = k + 18

            a = DN(i, i) / 8
            b = DN(k, i) / 8
            Do n = 1, iGauss
               U(1, i, n) = U(1, i, n) + a * U(1,i1,n) + b * U(1,k1,n)
            End do

            dx = XYP(1, j) - XYP(1, i)
            ex = XYP(1, k) - XYP(1, i)
            a = (24 * dx * DN(i, i) + 3 * XYN(1, i)) / 16
            b = (24 * ex * DN(k, i) + 3 * XYN(1, k)) / 16
            Do n = 1, iGauss
               U(1,i+3,n) = U(1,i+3,n) - a * U(1,i1,n) - b * U(1,k1,n)
            End do

            dy = XYP(2, j) - XYP(2, i)
            ey = XYP(2, k) - XYP(2, i)
            a = (24 * dy * DN(i, i) + 3 * XYN(2, i)) / 16
            b = (24 * ey * DN(k, i) + 3 * XYN(2, k)) / 16
            Do n = 1, iGauss
               U(1,i+6,n) = U(1,i+6,n) - a * U(1,i1,n) - b * U(1,k1,n)
            End do

            a = (3 * DN(i, i) * dx + 2 * XYN(1, i)) * dx / 32 
            b = (3 * DN(k, i) * ex + 2 * XYN(1, k)) * ex / 32
            Do n = 1, iGauss
               U(1,i+9,n) = U(1,i+9,n) - a * U(1,i1,n) - b * U(1,k1,n)
            End do

            a = (3 * DN(i, i) * dy + 2 * XYN(2, i)) * dy / 32
            b = (3 * DN(k, i) * ey + 2 * XYN(2, k)) * ey / 32
            Do n = 1, iGauss
               U(1,i+12,n) = U(1,i+12,n) - a * U(1,i1,n) - b * U(1,k1,n)
            End do

            dxy = dx * dy
            exy = ex * dy
            a = XYN(1, i) * dy + XYN(2, i) * dx
            b = XYN(1, k) * ey + XYN(2, k) * ex
            a = (3 * DN(i, i) * dxy + 2 * a) / 16
            b = (3 * DN(k, i) * exy + 2 * b) / 16
            Do n = 1, iGauss
               U(1,i+15,n) = U(1,i+15,n) - a * U(1,i1,n) - b * U(1,k1,n)
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

