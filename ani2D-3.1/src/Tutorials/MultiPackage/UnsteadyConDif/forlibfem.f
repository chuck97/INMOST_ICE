C ======================================================================
      Subroutine FEM2DextStif(XY1, XY2, XY3,
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*)

      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Dconv, Drhs, Dbc
      External Ddiff, Dconv, Drhs, Dbc

      Integer  ANI_Dnull
      Real*8   tri_area, edge_length, tri_diameter
      EXTERNAL tri_area, edge_length, tri_diameter, ANI_Dnull

      Integer  idum, iT
      Real*8   B(3, 3), C(3, 3), FF, area, DiffCoef(4,4), ConvCoef(4,4)
      Real*8   x, y, eBC(1), Dbb, bb, bnrm, h, h1, h2, h3
      Real*8   Pe, delta, dT, ConcAvr

      Logical  ifXbc

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templates (we have 1 group, skipping additional bits)
      Do i = 1, 3
         templateR(i) = Vdof  ! + Group1
         templateC(i) = Vdof  ! + Group1
      End do



c ... extract data from exchange arrays
      iT      = iSYS(3)     ! index of triangle
      dT      = dDATA(3)     ! time step
      ConcAvr = dDATA(3+iT)  ! mean value of 2u^{n-1}-u^{n}/2


c ... compute the stiffness matrix 
      label = lbE

      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, dDATA, iDATA, iSYS, 1,
     &              LDA, A, ir, ic)

c ... compute the convection matrix 
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P1, GRAD, FEM_P1,
     &              label, Dconv, dDATA, iDATA, iSYS, 2,
     &              3, B, ir, ic)

c ... compute the mass matrix 
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P1, IDEN, FEM_P1,
     &              label, ANI_Dnull, dDATA, iDATA, iSYS, 2,
     &              3, C, ir, ic)

      Do i = 1, 3
         Do j = 1, 3
            A(i, j) = A(i, j) + B(i, j) + C(i, j) * 1.5d0 / dT
         End do
      End do

c ... compute right hand side
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 1,
     &              1, F, ir, ic)


c ... SUPG contribution
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P0,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 1,
     &              1, FF, ir, ic)

      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, GRAD, FEM_P1,
     &              label, Dconv, dDATA, iDATA, iSYS, 1,
     &              3, C, ir, ic)

      area = dabs(tri_area(XY1, XY2, XY3))
      Call tri_center(XY1, XY2, XY3, x, y)

      idum = Ddiff(x, y, label, dDATA, iDATA, iSYS, DiffCoef)
      idum = Dconv(x, y, label, dDATA, iDATA, iSYS, ConvCoef)
      
      Dbb = DiffCoef(1,1) * ConvCoef(1,1)**2 +
     &      DiffCoef(2,2) * ConvCoef(2,1)**2 + 
     &      DiffCoef(1,2) * ConvCoef(1,1)*ConvCoef(2,1)*2

      bb  = ConvCoef(1,1)**2 + ConvCoef(2,1)**2 
      bnrm  = dsqrt( bb )

      h = tri_diameter(XY1, XY2, XY3)

      If(bb.NE.0D0) Then
          Dbb = Dbb / bb
          Pe = h * bnrm / Dbb
      Else
          Pe = 0D0
      End if
     
      delta = 0D0
      If(Pe.GE.1) Then
         delta = (h / (2*bnrm)) * (1D0 - 1D0/Pe) 
      End if 

      Do i = 1, 3
         Do j = 1, 3
            A(i, j) = A(i, j) + delta * C(1,i) * 
     &               (C(1, j) / area + 0.5d0 / dT)
         End do
      End do

      Do i = 1, 3
         F(i) = F(i) + delta * C(1,i) * (FF / area + ConcAvr / dT)
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
               Call applyDIR(LDA, nRow, A, F, k, eBC(1))
            End if
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine FEM2DextMass(XY1, XY2, XY3,
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*)
      Integer lbE, lbF(3), lbP(3)

      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  ANI_D1x1_scalar
      External ANI_D1x1_scalar

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templates
      Do i = 1, nRow
         templateR(i) = Vdof   
         templateC(i) = Vdof   
      End do

c ... compute the mass matrix 
      label = lbE

      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P1, IDEN, FEM_P1,
     &              label, ANI_D1x1_scalar, dDATA, iDATA, iSYS, 2,
     &              LDA, A, ir, ic)

       Return
       End



C ======================================================================
C  Diffusion tensor             
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
         Coef(i, i) = 1D-4
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
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 1

      Conv(1,1) = dDATA(1)
      Conv(2,1) = dDATA(2)

      Dconv = TENSOR_GENERAL

      Return
      End



C ======================================================================
C Boundary condition
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, eBC(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      If(label.eq.1) Then      ! homogeneous Dirichlet on edges labeled 1
         Dbc = BC_DIRICHLET 
         eBC(1,1) = 0D0

      Else If(label.eq.2) Then ! homogeneous Dirichlet on edges labeled 2
         Dbc = BC_DIRICHLET
         eBC(1,1) = 1D0

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

      F(1,1) = 0D0
      Drhs = TENSOR_SCALAR

      Return
      End


