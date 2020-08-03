c ======================================================================
      Program  mainMatrix
c ======================================================================
c This program creates and plots a simple three-diagonal matrix. 
c ======================================================================
      implicit none
      integer nfmax
c ... nfmax - maximum number of matrix columns and rows
      parameter(nfmax = 20)

c ... matrix
      Integer IA(nfmax+1), JA(3*nfmax)
      
C LOCAL VARIABLEs      
      integer i, nz

c ======================================================================
      nz = 0
      IA(1) = 1

      Do i = 1, nfmax
         nz = nz + 1
         JA(nz) = i

         If(i.GT.1) Then
            nz = nz + 1
            JA(nz) = i - 1
         End if

         If(i.LT.nfmax) Then
            nz = nz + 1
            JA(nz) = i + 1
         End if

         IA(i + 1) = nz + 1
      End do

c Make a PostScript file of the matrix. The name must terminate with .ps
      Call draw_matrix(nfmax, IA, JA, 'matrix.ps')

      Stop 
      End


