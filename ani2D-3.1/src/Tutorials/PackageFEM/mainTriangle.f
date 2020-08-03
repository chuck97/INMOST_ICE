c =====================================================================
      Program  mainTriangle
c =====================================================================
c Generation of elemental mass and stiffness matrices.
c =====================================================================
      implicit none

c ... standard mesh arrays (see doc/user_guide.pdf for more detail)
c ... number of points, triangles, boundary edges, and curved edges
      Integer  nv, nt, nb

c ... coordinates of mesh points 
      Real*8   vrt(2,3)

c ... connectivity table for triangles and triangle labels
      Integer  tri(3,1), material(1)

c ... connectivity table for boundary edges, and edge labels
      Integer  bnd(4,3)


c ... library FEM
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Real*8   dDATA(1), A(Maxsize, MaxSize)
      Integer  iDATA(1), nRow, nCol, order, LDA, iSYS(MAXiSYS)

      EXTERNAL Ddiff
      Integer  Ddiff


c ... local variables
      Integer  i, j, label

c =====================================================================
c ... Define a simple mesh consisting of 1 triangle 
c     used convetion: 'v'=vertex, 't'=triangle, 'b'=boundary edge
      nv = 3
      vrt(1,1) = 0.0d0
      vrt(2,1) = 0.0d0
      vrt(1,2) = 1.0d0
      vrt(2,2) = 0.0d0
      vrt(1,3) = 0.0d0
      vrt(2,3) = 1.0d0

      nt = 1
      tri(1,1) = 1
      tri(2,1) = 2
      tri(3,1) = 3
      material(1) = 1 

      Write(*,'(A)') 'Triangle with the following vertices:'
      Do i = 1, 3
         Write(*,'(I3,3F8.3)') i, vrt(1, i), vrt(2, i)
      End do

c     no boundary edges
      nb = 0

c     no extra data is provided for the user subroutine Ddiff
      dDATA(1) = 0D0
      iDATA(1) = 0

c     order of the quadrature rule 
      order = 2

c     leading dimension of the local matrix A
      LDA = MaxSize

c ... L2-product of BDM1 functions. No need to populate iSYS
      label = material(1)
      Call fem2Dtri(vrt(1,1), vrt(1,2), vrt(1,3),   
     &              IDEN, FEM_BDM1, IDEN, FEM_BDM1, 
     &              label, Ddiff, dDATA, iDATA, iSYS, order,  
     &              LDA, A, nRow, nCol)

      Write(*,'(/,A)') 'TEST 1: Bilinear form <BDM1, BDM1>'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol

      Do i = 1, nRow
         Write(*, '(100E14.6)') (A(i, j), j = 1, nCol)
         If(A(i, i).LE.0D0) Stop 911
      End do


c ... integral of div(BDM1), No need to populate iSYS
      label = material(1)
      Call fem2Dtri(vrt(1,1), vrt(1,2), vrt(1,3), 
     &              DIV, FEM_BDM1, IDEN, FEM_P0, 
     &              label, Ddiff, dDATA, iDATA, iSYS, order,  
     &              LDA, A, nRow, nCol)

      Write(*,'(/,A)') 'TEST 2: Linear form <div(BDM1), 1>'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol

      Do i = 1, nRow
         Write(*, '(100E14.6)') (A(i, j), j = 1, nCol)
      End do


c ... L2-product of CURL(P1) functions. No need to populate iSYS
      label = material(1)
      Call fem2Dtri(vrt(1,1), vrt(1,2), vrt(1,3), 
     &              CURL, FEM_P1, CURL, FEM_P1, 
     &              label, Ddiff, dDATA, iDATA, iSYS, order,  
     &              LDA, A, nRow, nCol)

      Write(*,'(/,A)') 'TEST 3: Bilinear form <curl(P1), curl(P1)>'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol

      Do i = 1, nRow
         Write(*, '(100E14.6)') (A(i, j), j = 1, nCol)
         If(A(i, i).LE.0D0) Stop 911
      End do


c ... L2-product of grad(P3) functions. No need to populate iSYS
      label = material(1)
      order = 9
      Call fem2Dtri(vrt(1,1), vrt(1,2), vrt(1,3), 
     &              GRAD, FEM_P3, GRAD, FEM_P3, 
     &              label, Ddiff, dDATA, iDATA, iSYS, order,  
     &              LDA, A, nRow, nCol)

      Write(*,'(/,A)') 'TEST 4: Bilinear form <grad(P3), grad(P3)>'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol

      Do i = 1, nRow
         Write(*, '(100E12.4)') (A(i, j), j = 1, nCol)
         If(A(i, i).LE.0D0) Stop 911
      End do


c ... L2-product of P1 functions on edge 12. No need to populate iSYS
      label = material(1)
      order = 5
      Call fem2Dedge(vrt(1,1), vrt(1,2), vrt(1,3), 
     &               DUDN, FEM_P2, IDEN, FEM_P1, 
     &               label, Ddiff, dDATA, iDATA, iSYS, order,  
     &               LDA, A, nRow, nCol)

      Write(*,'(/,A)') 'TEST 5: Bilinear form <grad(P1).n, P1> on edge'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol

      Do i = 1, nRow
         Write(*, '(100E12.4)') (A(i, j), j = 1, nCol)
      End do


c ... testing...
      label = material(1)
      order = 9
      Call fem2Dtri(vrt(1,1), vrt(1,2), vrt(1,3), 
     &              GRAD, FEM_HERMIT, GRAD, FEM_HERMIT, 
     &              label, Ddiff, dDATA, iDATA, iSYS, order,  
     &              LDA, A, nRow, nCol)
 
      Write(*,'(/,A)') 'TEST 6: Bilinear form <HERMIT, HERMIT>'
      Write(*,'(A,2I4)') 'Elemental matrix size:', nRow, nCol
 
      Do i = 1, nRow
         Write(*, '(100E14.6)') (A(i, j), j = 1, nCol)
      End do


c ... testing the results

      Stop 
      End



C ======================================================================
C  The user defined routines required above
C ======================================================================
c ... 2x2 diffusion tensor Diff
      Integer Function Ddiff(x, y, label, dDATA, iDATA, iSYS, Coef)
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = 1D0
      Ddiff = TENSOR_SCALAR

      Return
      End



