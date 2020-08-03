c ======================================================================
      Subroutine applyDUDX(iGauss, XYL, PSI, FEMtype, nfa, dim, U)
c ======================================================================
      Implicit None
      include 'fem2Dtri.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(2, 2), U(4, MaxDOFs, *)

c Data for the reference triangle
      Real*8  XYL(3, *)
      Real*8  GRAD_P1(2, 3), GRAD_PX(2, MaxDOFs, MaxPointGauss)

      Integer iref(4), i,j,k,l,n, ixy, mfa
      Real*8  x, y, Lfun

      DATA    GRAD_P1/-1D0,-1D0, 1D0,0D0, 0D0,1D0/
      DATA    iref /1,2,3,1/

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
         dim = 1 

c  ...   choose ixy = 1 for DuDx and ixy = 2 for DuDy
         ixy = 1
         Do i = 1, nfa
            U(1, i, 1) = 0D0
            Do j = 1, 2
               U(1, i, 1) = U(1, i, 1) + PSI(j, ixy) * GRAD_P1(j, i)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1) Then
         nfa = 3
         dim = 1 

c  ...   choose ixy = 1 for DuDx and ixy = 2 for DuDy
         ixy = 1
         Do i = 1, nfa
            l = iref(i + 1)
            U(1, l, 1) = 0D0
            Do j = 1, 2
               U(1, l, 1) = U(1, l, 1) - 2*PSI(j, ixy) * GRAD_P1(j, i)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_P2 .OR. FEMtype.EQ.FEM_P2vector) Then
         nfa = 6
         dim = 1

         If(FEMtype.EQ.FEM_P2vector) Call clearU(iGauss, 12, 2, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            Do i = 1, 3
               Do k = 1, 2
                  GRAD_PX(k, i, n) = GRAD_P1(k, i) *
     &                               (4 * Lfun(i, x, y) - 1D0) 
               End do
            End do

            mfa = 3
            Do i = 1, 3
               j = iref(i + 1)
               mfa = mfa + 1
               Do k = 1, 2
                  GRAD_PX(k, mfa, n) = 
     &               4 * (Lfun(i, x, y) * GRAD_P1(k, j) +
     &                    Lfun(j, x, y) * GRAD_P1(k, i))
               End do
            End do
         End do

c  ...   choose ixy = 1 for DuDx and ixy = 2 for DuDy
         ixy = 1
         Do n = 1, iGauss
            Do i = 1, nfa
               U(1, i, n) = 0D0
               Do j = 1, 2
                  U(1, i, n) = U(1, i, n) 
     &                       + PSI(j, ixy) * GRAD_PX(j, i, n)
               End do
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

      Else
         nfa = 0
         dim = -FEMtype
         Call errMesFEM(2001, 'fem2Dtri', 
     &        'unsupported operation for the given element type')
      End if

      Return
      End



c ======================================================================
      Subroutine applyDUDY(iGauss, XYL, PSI, FEMtype, nfa, dim, U)
c ======================================================================
      include 'fem2Dtri.fd'
c ======================================================================
      Integer iGauss, FEMtype, nfa, dim
      Real*8  PSI(2, 2), U(4, MaxDOFs, *)

c Data for the reference triangle
      Real*8  XYL(3, *)
      Real*8  GRAD_P1(2, 3), GRAD_PX(2, MaxDOFs, MaxPointGauss)
      Real*8  x, y, Lfun

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

c ... next two FEMs
      Else If(FEMtype.EQ.FEM_P1) Then
         nfa = 3
         dim = 1 

c  ...   choose ixy = 1 for DuDx and ixy = 2 for DuDy
         ixy = 2
         Do i = 1, nfa
            U(1, i, 1) = 0D0
            Do j = 1, 2
               U(1, i, 1) = U(1, i, 1) + PSI(j, ixy) * GRAD_P1(j, i)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next FEM
      Else If(FEMtype.EQ.FEM_CR1) Then
         nfa = 3
         dim = 1 

c  ...   choose ixy = 1 for DuDx and ixy = 2 for DuDy
         ixy = 2
         Do i = 1, nfa
            l = iref(i + 1)
            U(1, l, 1) = 0D0
            Do j = 1, 2
               U(1, l, 1) = U(1, l, 1) - 2*PSI(j, ixy) * GRAD_P1(j, i)
            End do
         End do
         Call copyGauss(iGauss, nfa, dim, U)

c ... next three FEMs
      Else If(FEMtype.EQ.FEM_P2 .OR. FEMtype.EQ.FEM_P2vector) Then
         nfa = 6
         dim = 1

         If(FEMtype.EQ.FEM_P2vector)  Call clearU(iGauss, 12, 2, U)

         Do n = 1, iGauss
            x = XYL(2, n)
            y = XYL(3, n)
            Do i = 1, 3
               Do k = 1, 2
                  GRAD_PX(k, i, n) = GRAD_P1(k, i) *
     &                               (4 * Lfun(i, x, y) - 1D0) 
               End do
            End do

            mfa = 3
            Do i = 1, 3
               j = iref(i + 1)
               mfa = mfa + 1
               Do k = 1, 2
                  GRAD_PX(k, mfa, n) = 
     &               4 * (Lfun(i, x, y) * GRAD_P1(k, j) +
     &                    Lfun(j, x, y) * GRAD_P1(k, i))
               End do
            End do
         End do

c  ...   choose ixy = 1 for DuDx and ixy = 2 for DuDy
         ixy = 2
         Do n = 1, iGauss
            Do i = 1, nfa
               U(1, i, n) = 0D0
               Do j = 1, 2
                  U(1, i, n) = U(1, i, n) 
     &                       + PSI(j, ixy) * GRAD_PX(j, i, n)
               End do
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


      Else
         nfa = 0
         dim = -FEMtype
         Call errMesFEM(2001, 'fem2Dtri', 
     &        'unsupported operation for the given element type')
      End if

      Return
      End

