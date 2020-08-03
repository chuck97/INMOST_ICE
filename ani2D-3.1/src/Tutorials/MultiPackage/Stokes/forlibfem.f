C ======================================================================
      Subroutine FEM2Dext(XY1, XY2, XY3, 
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      implicit none
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*)
      
      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Integer templateR(*), templateC(*), LDA, nRow, nCol
      Real*8  A(LDA, *), F(*)

C Local variables
      Integer  Ddiff, Drhs, Dbc, ANI_Dnull
      External Ddiff, Drhs, Dbc, ANI_Dnull

      Integer  i,j,k, label, ir,ic, ibc
      Real*8   B(12, 12), x, y, eBC(2)
      Logical  ifXbc

C ======================================================================
      nRow = 15
      nCol = 15

c Populate degrees of freedom for rows and columns: 9 at vertices and 6 on edges.
c We use convention 'V'=vertex d.o.f. and 'R'=edge d.o.f.
c We breake unknowns into groups  to obtain a 5x5 block structure of the global matrix.
c The first (nodal u_x), second (nodal u_y) and third (nodal p)  diagonal blocks will have the same size equal to the number of mesh points.
c The remaining diagonal blocks (edge u_x, edge u_y) will have size equal to the number of mesh edges.
c For the ordering of unknowns inside these blocks please refer to section "Assembling utilities" of the user_guide
      Do i = 1, 3
         templateR(i)     = Vdof + Group1  !u_x
         templateR(i + 3) = Rdof + Group4

         templateR(i + 6) = Vdof + Group2  !u_y
         templateR(i + 9) = Rdof + Group5

         templateR(i + 12)= Vdof + Group3 !p
      End do

      Do k = 1, nRow
         templateC(k) = templateR(k)
      End do

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:12,1:12) is elemental vector Laplacian matrix for velocity components
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P2vector, GRAD, FEM_P2vector,
     &              label, Ddiff, dDATA, iDATA, iSYS, 2,
     &              LDA, A, ir, ic)


c     B(1:3,1:12) is elemental divergence matrix for velocity
      Call fem2Dtri(XY1, XY2, XY3,
     &              DIV, FEM_P2vector, IDEN, FEM_P1,
     &              label, ANI_Dnull, dDATA, iDATA, iSYS, 4,
     &              12, B, ir, ic)

      Do i = 1, 12
         Do j = 1, 3
            A(i, j + 12) = B(i, j)
            A(j + 12, i) = B(i, j)
         End do
      End do

      Do i = 1, 3
         Do j = 1, 3
            A(i + 12, j + 12) = 0D0
         End do
      End do

c ... compute right hand side vector
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P2vector, 
     &              lbE, Drhs, dDATA, iDATA, iSYS, 2,
     &              1, F, ir, ic)

      Do i = 1, 3
         F(i + 12) = 0
      End do

c ... impose boundary conditions (assume nRow = nCol) at nodes of triangle
      Do k = 1, 3
         If(lbP(k).ne.0) Then
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

            ibc = Dbc(x, y, lbP(k), dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_DIRICHLET)) Then
               Call applyDIR(LDA, nRow, A, F, k,   eBC(1))
               Call applyDIR(LDA, nRow, A, F, k+6, eBC(2))
            End if
         End if
      End do

c ... impose boundary conditions (assume nRow = nCol) at middles of edges
      Do k = 1, 3
         If(lbF(k).GT.0) Then
            If(k.EQ.1) Then
               x = (XY1(1) + XY2(1)) / 2
               y = (XY1(2) + XY2(2)) / 2
            ElseIf(k.EQ.2) Then
               x = (XY2(1) + XY3(1)) / 2
               y = (XY2(2) + XY3(2)) / 2
            ElseIf(k.EQ.3) Then
               x = (XY3(1) + XY1(1)) / 2
               y = (XY3(2) + XY1(2)) / 2
            End if

            ibc = Dbc(x, y, lbF(k), dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_DIRICHLET)) Then
               Call applyDIR(LDA, nRow, A, F, k+3, eBC(1))
               Call applyDIR(LDA, nRow, A, F, k+9, eBC(2))
            End if
         End if
      End do

      Return
      End



C ======================================================================
C Identity tensor
C ======================================================================
      Integer Function Ddiff(x, y, label, dDATA, iDATA, iSYS, Coef)
      implicit none
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, 4)
      Integer iDATA(*), iSYS(*), label, i, j

      iSYS(1) = 4
      iSYS(2) = 4

      Do i = 1, 4
         Do j = 1, 4
            Coef(i, j) = 0
         End do
         Coef(i, i) = 1
      End do

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C Essential boundary conditions
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
      implicit none
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, eBC(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*), label

      iSYS(1) = 2
      iSYS(2) = 1

      If(label.eq.3) Then
         eBC(1,1) = 1 - 4*y*y
         eBC(2,1) = 0
         Dbc = BC_DIRICHLET

      Else If(label.eq.1) Then
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = BC_NEUMANN

      Else
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = BC_DIRICHLET
      End if

      Return
      End



C ======================================================================
C Right hand side
C ======================================================================
      Integer Function Drhs(x, y, label, dDATA, iDATA, iSYS, F)
      implicit none
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, F(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*), label

      iSYS(1) = 2
      iSYS(2) = 1
     
      F(1,1) = 0D0
      F(2,1) = 0D0

      Drhs = TENSOR_GENERAL

      Return
      End 


