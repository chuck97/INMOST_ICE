c ======================================================================
      Subroutine applyTRACE(iGauss, XYL, FEMtype, nfa, dim, U)
c ======================================================================
c Routine calculates trace of a finite element function at edge 
c barycentric points given by XYL(2).
c ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer iGauss, FEMtype, nfa, dim
      Real*8  XYL(2, *), U(2, MaxDOFs, *)

      Integer i, n
      Real*8  x

c ======================================================================
      If(FEMtype.EQ.FEM_P0) Then
         nfa = 1
         dim = 1
         Do n = 1, iGauss
            U(1, 1, n) = 1D0
         End do

      Else If(FEMtype.EQ.FEM_P1) Then
         nfa = 3
         dim = 1
         Do n = 1, iGauss
            x = XYL(2, n) 
            U(1, 1, n) = 1 - x
            U(1, 2, n) = x
            U(1, 3, n) = 0D0
         End do

      Else If(FEMtype.EQ.FEM_P2) Then
         nfa = 6
         dim = 1
         Do n = 1, iGauss
            x = XYL(2, n) 
            U(1, 1, n) = (1 - x) * (1 - 2*x) 
            U(1, 2, n) = x * (2*x - 1) 
            U(1, 3, n) = 0D0
            U(1, 4, n) = 4 * x * (1 - x) 
            U(1, 5, n) = 0D0
            U(1, 6, n) = 0D0
         End do

      Else If(FEMtype.EQ.FEM_P3) Then
         nfa = 10
         dim = 1
         Do n = 1, iGauss
            Do i = 3, nfa
               U(1, i, n) = 0D0
            End do

            x = XYL(2, n) 
            U(1, 1, n) = (1 - x) * (1 - 3*x) * (2 - 3*x) / 2
            U(1, 2, n) = x * (3*x - 1) * (3*x - 2) / 2
            U(1, 4, n) = 9 * x * (2 - 3*x) * (1 - x) / 2 
            U(1, 7, n) = 9 * x * (3*x - 1) * (1 - x) / 2 
         End do

      Else
         nfa = 0
         dim = -FEMtype
         Call errMesFEM(2001, 'fem2Dedge', 
     &        'unsupported operation for the given element type')
      End if

      Return
      End

