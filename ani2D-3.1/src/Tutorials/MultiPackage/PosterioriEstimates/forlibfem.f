C ======================================================================
      Subroutine FEM2DextLinear(XY1, XY2, XY3, 
     &                          lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &                          LDA, A, F, nRow, nCol,
     &                          templateR, templateC)
C ======================================================================
C Generates a 3x3 stiffness matrix for the P_1 finite elements and the 
C right-hand side. dDATA() and iDATA() are not used.
C ======================================================================
      implicit none
      Include 'fem2Dtri.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*)
      
      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Integer LDA, nRow, nCol
      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Drhs, Dbc
      External Ddiff, Drhs, Dbc

      Integer  i, k, ir, ic
      Real*8   x, y, eBC(1)

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up DOFs (we have 1 group, skipping additional bits)
      Do i = 1, 3
         templateR(i) = Vdof ! + Group1
         templateC(i) = Vdof ! + Group1
      End do

c ... compute the stiffness matrix 
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              lbE, Ddiff, dDATA, iDATA, iSYS, 2,
     &              LDA, A, ir, ic)

c ... compute right hand side
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 2,
     &              1, F, ir, ic)

c ... impose Dirichlet boundary conditions (assume nRow = nCol)
      Do k = 1, 3
         If(lbP(k).GT.0) Then
            If(k.EQ.1) Then
               x = XY1(1)
               y = XY1(2)
            ElseIf(k.EQ.2) Then
               x = XY2(1)
               y = XY2(2)
            ElseIf(k.EQ.3) Then
               x = XY3(1)
               y = XY3(2)
            End if

            i = Dbc(x, y, lbP(k), dDATA, iDATA, iSYS, eBC)

            Do i = 1, nRow
               F(i) = F(i) - A(i, k) * eBC(1)
               A(k, i) = 0D0
               A(i, k) = 0D0
            End do 

            A(k, k) = 1D0
            F(k) = eBC(1)
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine FEM2DextBubble( 
C ======================================================================
     &           XY1, XY2, XY3, 
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
C Generates a 3x3 stiffness matrix for bubbles and the right-hand side
C vector using the residual of the solution at vertices. dDATA() is the
C solution vector at mesh vertices. iDATA() is not used.
C ======================================================================
      implicit none
      Include 'fem2Dtri.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*)
      
      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Integer LDA, nRow, nCol
      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Drhs, Dbc
      External Ddiff, Drhs, Dbc

      Integer  i, j, k, iP1, iP2, iP3, ir, ic
      Real*8   B(6, 6), G(6), Uv(3), Ur(3)

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templated degrees of freedom, all three on triangle edges
      Do i = 1, 3
         templateR(i) = Rdof
         templateC(i) = Rdof
      End do

c ... compute the stiffness matrix B
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P2, GRAD, FEM_P2,
     &              lbE, Ddiff, dDATA, iDATA, iSYS, 2,
     &              6, B, ir, ic)

c ... compute right hand side
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P2,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 2,
     &              1, G, ir, ic)

c ... bubble block to the stiffness matrix
      Do i = 1, 3
         Do j = 1, 3
            A(i, j) = B(i+3, j+3)
         End do
         F(i) = G(i+3)
      End do

c ... elliminate the solution at nodes
      iP1 = iSYS(4)
      iP2 = iSYS(5)
      iP3 = iSYS(6)

      Uv(1) = dDATA(iP1)
      Uv(2) = dDATA(iP2)
      Uv(3) = dDATA(iP3)

      Ur(1) = (dDATA(iP1) + dDATA(iP2)) / 2
      Ur(2) = (dDATA(iP2) + dDATA(iP3)) / 2
      Ur(3) = (dDATA(iP3) + dDATA(iP1)) / 2

      Do j = 1, 3
         Do i = 1, 3
            F(i) = F(i) - B(i+3, j) * Uv(j) - B(i+3, j+3) * Ur(j) 
         End do
      End do

c ... impose Dirichlet boundary conditions (assume nRow = nCol)
      Do k = 1, 3
         If(lbF(k).GT.0) Then
            Do i = 1, nRow
               A(k, i) = 0
               A(i, k) = 0
            End do 

            A(k, k) = 1
            F(k) = 0D0
         End if
      End do

      Return
      End



C ======================================================================
C Identity tensor
C ======================================================================
      Integer Function Ddiff(x, y, label, dDATA, iDATA, iSYS, Coef)
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 2

      Do i = 1, 2
         Do j = 1, 2
            Coef(i, j) = 0D0
         End do
      End do

      Coef(1, 1) = 1D0
      Coef(2, 2) = 1D3

      Call rotateTensor(4, Coef, 250*(x+y))

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C Essential boundary conditions
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, eBC(*)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      eBC(1) = 0D0 
      Dbc = BC_DIRICHLET

      Return
      End



C ======================================================================
C Sourse data
C ======================================================================
      Integer Function Drhs(x, y, label, dDATA, iDATA, iSYS, F)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, F(*)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      F(1) = 1D0
      Drhs = TENSOR_SCALAR

      Return
      End


