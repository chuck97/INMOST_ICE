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
      Integer  Ddiff, Dconv, Drhs, Dbc
      External Ddiff, Dconv, Drhs, Dbc

      Integer  i, j, k, label, ir, ic, ibc
      Real*8   B(3, 3), x, y, eBC(1)
      Logical  ifXbc

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templates (we have 1 group, skipping additional bits)

      Do i = 1, nRow
         templateR(i) = Vdof ! + Group1
      End do

      Do k = 1, nCol
         templateC(k) = templateR(k)
      End do

c ... compute the stiffness matrix 
      label = lbE

      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, dDATA, iDATA, iSYS, 1,
     &              LDA, A, ir, ic)

      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P1, GRAD, FEM_P1,
     &              label, Dconv, dDATA, iDATA, iSYS, 1,
     &              3, B, ir, ic)

      Do i = 1, 3
         Do j = 1, 3
            A(i, j) = A(i, j) + B(i, j)
         End do
      End do

c ... compute right hand side
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 1,
     &              1, F, ir, ic)

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

c  ...  Dbc may change iSYS
            ibc = Dbc(x, y, lbP(k), dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_DIRICHLET)) Then
               Call applyDIR(LDA, nRow, A, F, k, eBC(1))
            End if
         End if
      End do

      Return
      End



C ======================================================================
C  Diffusion tensor             
C ======================================================================
      Integer Function Ddiff(x, y, label, dDATA, iDATA, iSYS, Coef)
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)


      iSYS(1) = 2
      iSYS(2) = 2

      Do i = 1, 2
         Do j = 1, 2
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = 1D-2
      End do

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C  Convection tensor             
C ======================================================================
      Integer Function Dconv(x, y, label, dDATA, iDATA, iSYS, Conv)
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Conv(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 1

      Conv(1,1) = dDATA(1)
      Conv(2,1) = dDATA(2)

      Dconv = TENSOR_GENERAL

      Return
      End



C ======================================================================
C Boundary condition: whole boundary is Dirichlet
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, eBC(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      If(label.EQ.1) Then ! homogeneous Dirichlet on edges with label 1
         Dbc = BC_DIRICHLET
         eBC(1, 1)  = 0D0
      Else
         Write(*,*) 'Dbc: wrong label=', label
         Stop
      End if

      Return
      End



C ======================================================================
C Right hand side
C ======================================================================
      Integer Function Drhs(x, y, label, dDATA, iDATA, iSYS, F)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, F(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      F(1, 1) = 1D0
      Drhs = TENSOR_SCALAR

      Return
      End


