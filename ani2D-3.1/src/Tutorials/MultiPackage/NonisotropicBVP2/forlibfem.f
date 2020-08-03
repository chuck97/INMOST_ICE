C ======================================================================
      Subroutine FEM2Dext(XY1, XY2, XY3, 
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      Implicit None
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*)
      
      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*), LDA, nRow, nCol

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Drhs, Dbc, ANI_Dnull
      External Ddiff, Drhs, Dbc, ANI_Dnull

      Integer  i,j,k, ir,ic, label, ibc
      Real*8   B(3, 3), x, y, eBC(1)
      Logical  ifXbc

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templates (we have 1 group, skipping additional bits)
      Do i = 1, 3
         templateR(i) = Vdof  ! + Group1
         templateC(i) = Vdof  ! + Group1
      End do

c ... compute the stiffness matrix A
      label = lbE
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, dDATA, iDATA, iSYS, 1,
     &              LDA, A, ir, ic)

      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P1, DUDX, FEM_P1,
     &              label, ANI_Dnull, dDATA, iDATA, iSYS, 2,
     &              3, B, ir, ic)

       Do i = 1, 3
          Do j = 1, 3
             A(i, j) = A(i, j) + B(i, j)
          End do
       End do

c ... compute right hand side F
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P1,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 2,
     &              1, F, ir, ic)


c ... impose boundary conditions (assume nRow = nCol)
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

            ibc = Dbc(x, y, lbP(k), dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_DIRICHLET)) Then
               Call applyDIR(LDA, nRow, A, F, k, eBC)
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

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, 4)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 2

      Do i = 1, 2
         Do j = 1, 2
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = 1D-3
      End do

      Ddiff = TENSOR_SYMMETRIC

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

      If(label.EQ.1) Then ! homogeneous Dirichlet on edges with label 1
         Dbc = BC_DIRICHLET
         eBC(1,1) = 0D0

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

      F(1,1) = 1D0
      Drhs = TENSOR_SCALAR

      Return
      End


