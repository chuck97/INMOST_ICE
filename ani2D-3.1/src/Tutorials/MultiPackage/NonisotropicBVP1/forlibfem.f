C ======================================================================
      Subroutine FEM2Dext(XY1, XY2, XY3, 
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      Implicit none
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
      Integer  Ddiff, Drhs, Dbc
      External Ddiff, Drhs, Dbc

      Integer  i,k, ir,ic, label, ibc
      Real*8   x, y, eBC(1)
      Logical  ifXbc

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

c ... compute the stiffness matrix A
      label = lbE
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, dDATA, iDATA, iSYS, 1,
     &              LDA, A, ir, ic)

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
      End do

      If(label.eq.1) Then  !isotropic in domain with material 1
         Coef(1, 1) = 1D1
         Coef(2, 2) = 1D1

      Else If(label.eq.2) Then !anisotropic in domain with material 2
         Coef(1, 1) = 1D0
         Coef(2, 2) = 1D2

      Else 
         Write(*,*) 'Ddiff: wrong label=', label
         Stop
      End if

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

      If (label.eq.5) Then !homogeneous Neumann on edges with label 5
         Dbc = BC_NEUMANN
         eBC(1,1) = 0D0

      Else If (label.eq.6) Then !interior boundary edges with label 6 are ignored
         Dbc = BC_NULL   
         eBC(1,1) = 0D0

      Else If ( label.ge.1 .and. label.le.4 ) Then !homogeneous Dirichlet on edges with labels 1-4
         Dbc = BC_DIRICHLET
         eBC(1,1) = 0D0

      Else
         write(*,*) 'Dbc: wrong label=', label
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


