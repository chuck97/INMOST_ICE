C ======================================================================
      Subroutine FEM2DextLinear(XY1, XY2, XY3, 
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
C Generates a 3x3 stiffness matrix for P_1 finite elements and the 
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

c ... set up templates (we have 1 group, skipping additional bits)
      Do i = 1, 3
         templateR(i) = Vdof ! + Group1
      End do
      Do i = 1, 3
         templateC(i) = templateR(i)
      End do


c ... compute the stiffness matrix 
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              lbE, Ddiff, dDATA, iDATA, iSYS, 2,
     &              LDA, A, ir, ic)

c ... compute right hand side
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 5,
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
      Subroutine FEM2DextBubble(XY1, XY2, XY3, 
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
C Generates a 3x3 stiffness matrix for bubbles and the right-hand 
C side using the residual of the solution at vertices. dDATA(*) is the
C solution vector at mesh vertices.
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
      Real*8   B(6, 6), C(6, 6), G(6), Uv(3), eBC(1), Ur, x, y

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templates 
      Do i = 1, 3
         templateR(i) = Rdof
         templateC(i) = Rdof
      End do

c ... compute the stiffness matrix B
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P2, GRAD, FEM_P2,
     &              lbE, Ddiff, dDATA, iDATA, iSYS, 5,
     &              6, B, ir, ic)

c ... compute right hand side
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P2,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 5,
     &              1, G, ir, ic)

c ... bubble block to the stiffness matrix
      Do i = 1, 3
         Do j = 1, 3
            A(i, j) = B(i+3, j+3)
         End do
         F(i) = G(i+3)
      End do

c ... elliminate the solution at nodes
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P2, GRAD, FEM_P1,
     &              lbE, Ddiff, dDATA, iDATA, iSYS, 2,
     &              6, C, ir, ic)

      iP1 = iSYS(4)
      iP2 = iSYS(5)
      iP3 = iSYS(6)

      Uv(1) = dDATA(iP1)
      Uv(2) = dDATA(iP2)
      Uv(3) = dDATA(iP3)

      Do j = 1, 3
         Do i = 1, 3
            F(i) = F(i) - C(i+3, j) * Uv(j)
         End do
      End do

c ... impose non-homogeneous Dirichlet boundary conditions (assume nRow = nCol)
      Do k = 1, 3
         If(lbF(k).GT.0 .AND. lbF(k).LE.5) Then
            If(k.EQ.1) Then
               x = (XY1(1) + XY2(1)) / 2
               y = (XY1(2) + XY2(2)) / 2
               Ur = (Uv(1) + Uv(2)) / 2
            Else If(k.EQ.2) Then
               x = (XY2(1) + XY3(1)) / 2
               y = (XY2(2) + XY3(2)) / 2
               Ur = (Uv(2) + Uv(3)) / 2
            Else If(k.EQ.3) Then
               x = (XY3(1) + XY1(1)) / 2
               y = (XY3(2) + XY1(2)) / 2
               Ur = (Uv(3) + Uv(1)) / 2
            End if

            i = Dbc(x, y, lbF(k), dDATA, iDATA, iSYS, eBC)
            eBC(1) = eBC(1) - Ur

            Do i = 1, nRow
               F(i) = F(i) - A(i, k) * eBC(1)
               A(k, i) = 0
               A(i, k) = 0
            End do 

            A(k, k) = 1
            F(k) = eBC(1)
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine FEM2DextFull(XY1, XY2, XY3, 
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
C Generates a 6x6 stiffness matrix for bubbles and the right-hand.
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
      Real*8   B(6, 6), C(6, 6), Uv(3), eBC(1), Ur, x, y

C ======================================================================
      nRow = 6
      nCol = 6

c ... set up templates 
      Do i = 1, 3
         templateR(i) = Vdof
         templateC(i) = Vdof
      End do

      Do i = 4, 6
         templateR(i) = Rdof
         templateC(i) = Rdof
      End do

c ... compute the stiffness matrix A
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              lbE, Ddiff, dDATA, iDATA, iSYS, 2,
     &              6, B, ir, ic)

      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P2, GRAD, FEM_P1,
     &              lbE, Ddiff, dDATA, iDATA, iSYS, 2,
     &              6, C, ir, ic)

      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P2, GRAD, FEM_P2,
     &              lbE, Ddiff, dDATA, iDATA, iSYS, 2,
     &              LDA, A, ir, ic)

      Do i = 1, 3
         Do j = 1, 3
            A(i,   j) = B(i,   j)
            A(i+3, j) = C(i+3, j)
            A(i, j+3) = C(i+3, j)
         End do
      End do

c ... compute right hand side
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P2,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 5,
     &              1, F, ir, ic)

c ... elliminate the solution at nodes
      iP1 = iSYS(4)
      iP2 = iSYS(5)
      iP3 = iSYS(6)

      Uv(1) = dDATA(iP1)
      Uv(2) = dDATA(iP2)
      Uv(3) = dDATA(iP3)

      Do j = 1, 3
         Do i = 1, 3
            F(i)   = 0
            F(i+3) = F(i+3) - A(i+3, j) * Uv(j)
         End do
      End do

c ... impose homogeneous Dirichlet boundary conditions (assume nRow = nCol)
      Do k = 1, 3
         If(lbP(k).GT.0) Then
            Do i = 1, nRow
               A(k, i) = 0
               A(i, k) = 0
            End do 

            A(k, k) = 1
            F(k) = 0D0
         End if
      End do

c ... impose non-homogeneous Dirichlet boundary conditions (assume nRow = nCol)
      Do k = 1, 3
         If(lbF(k).GT.0 .AND. lbF(k).LE.5) Then
            If(k.EQ.1) Then
               x = (XY1(1) + XY2(1)) / 2
               y = (XY1(2) + XY2(2)) / 2
               Ur = (Uv(1) + Uv(2)) / 2
            ElseIf(k.EQ.2) Then
               x = (XY2(1) + XY3(1)) / 2
               y = (XY2(2) + XY3(2)) / 2
               Ur = (Uv(2) + Uv(3)) / 2
            ElseIf(k.EQ.3) Then
               x = (XY3(1) + XY1(1)) / 2
               y = (XY3(2) + XY1(2)) / 2
               Ur = (Uv(3) + Uv(1)) / 2
            End if

            i = Dbc(x, y, lbF(k), dDATA, iDATA, iSYS, eBC)
            eBC(1) = eBC(1) - Ur

            Do i = 1, nRow
               F(i+3) = F(i+3) - A(i+3, k) * eBC(1)
               A(k+3, i) = 0
               A(i, k+3) = 0
            End do 

            A(k+3, k+3) = 1
            F(k+3) = eBC(1)
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
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = 1D0

      Ddiff = TENSOR_NULL

      Return
      End



C ======================================================================
C Essential boundary conditions
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
C ======================================================================
      implicit none
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, eBC(MaxTensorSize, *), r, phi, pi
      Integer iDATA(*), iSYS(*), label

      DATA    pi/3.14159265358979D0/

      iSYS(1) = 1
      iSYS(2) = 1

      If(label.EQ.6) Then !homogeneous Neumann on edges with label 7
         Dbc = BC_NEUMANN
         eBC(1,1) = 0D0

      Else If(label.EQ.5) Then !homogeneous Dirichlet on edges with labels 6
         Dbc = BC_DIRICHLET
         eBC(1,1) = 0D0

      Else If(label.GE.1 .AND. label.LE.4) Then !nonhomogeneous Dirichlet on edges with labels 1-4
         Dbc = BC_DIRICHLET
         If(x.GT.1d-10) Then
            phi = datan(y/x)
            If(phi.lt.0) phi = phi + 2*pi
         Else If(x.lt.-1d-10) Then
            phi = datan(y/dabs(x))
            phi = pi - phi
         Else
            If(y.GT.0) phi = pi/2
            If(y.LT.0) phi = 3*pi/2
         End if
         eBC(1,1) = dsin(phi/4)

      Else
         Write(*,*) 'Dbc: wrong label=', label
         Stop
      End if

      Return
      End



C ======================================================================
C Source data
C ======================================================================
      Integer Function Drhs(x, y, label, dDATA, iDATA, iSYS, Coef)
C ======================================================================
      implicit none
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1,1) = 0D0
      Drhs = TENSOR_SCALAR

      Return
      End



C ======================================================================
C Exact solution
C ======================================================================
      Integer Function DexactU(x, y, label, dDATA, iDATA, iSYS, Coef)
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *), r, phi, pi
      Integer iDATA(*), label, iSYS(*)

      DATA    pi/3.14159265358979D0/

      iSYS(1) = 1
      iSYS(2) = 1

      If(x.GT.1d-10) Then
         phi = datan(y/x)
         If(phi.lt.0) phi = phi + 2*pi
      Else If(x.lt.-1d-10) Then
         phi = datan(y/dabs(x))
         phi = pi - phi
      Else
         If(y.GT.0) phi = pi/2
         If(y.LT.0) phi = 3*pi/2
      End if

      r = (x*x + y*y + 1D-24) ** 0.125D0
      Coef(1, 1) = r * dsin(phi/4)

      DexactU = TENSOR_SCALAR

      Return 
      End



C ======================================================================
c Gradient of the exact solution
C ======================================================================
      Integer Function DgradU(x, y, label, dDATA, iDATA, iSYS, Coef)
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *), r, phi, pi
      Integer iDATA(*), label, iSYS(*)

      DATA    pi/3.14159265358979D0/

      iSYS(1) = 2
      iSYS(2) = 1

      If(x.GT.1d-10) Then
         phi = datan(y/x)
         If(phi.lt.0) phi = phi + 2*pi

      Else If(x.lt.-1d-10) Then
         phi = datan(y/dabs(x))
         phi = pi - phi

      Else
         If(y.GT.0) phi = pi/2
         If(y.LT.0) phi = 3*pi/2
      End if

      r = x*x + y*y + 1D-16
      r = 4 * (r ** 0.875D0)

      Coef(1, 1) = (x * dsin(phi/4) - y * dcos(phi/4)) / r
      Coef(2, 1) = (y * dsin(phi/4) + x * dcos(phi/4)) / r

      DgradU = TENSOR_GENERAL

      Return
      End

