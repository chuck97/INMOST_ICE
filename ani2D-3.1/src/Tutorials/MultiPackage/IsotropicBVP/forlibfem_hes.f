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

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1,1) = 1D0

      Ddiff = TENSOR_NULL

      Return
      End



C ======================================================================
C Boundary condition
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, eBC(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      Real*8  pi, phi
      Data    pi/3.14159265358979D0/

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
C Right hand side
C ======================================================================
      Integer Function Drhs(x, y, label, dDATA, iDATA, iSYS, F)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, F(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      F(1, 1) = 0D0
      Drhs = TENSOR_SCALAR

      Return
      End



c ======================================================================
c  Exact (analytical) solution
c ======================================================================
      Integer Function DexactU(x, y, label, dDATA, iDATA, iSYS, U)
      Include 'fem2Dtri.fd'

      Real*8   dDATA(*), x, y, U(MaxTensorSize, *)
      Integer  iDATA(*), label, iSYS(*)

      Real*8   phi, r, pi

c ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

      pi = 4d0*datan(1d0)

      If(x.gt.1d-10) Then
         phi = datan(y/x)
         If(phi.lt.0) phi = phi+2*pi

      Else If(x.lt.-1d-10) Then
         phi = datan(y/dabs(x))
         phi = pi - phi

      Else
         If(y.gt.0) phi = pi/2
         If(y.le.0) phi = 3*pi/2
      End if

      r = dsqrt(x**2 + y**2)
      U(1,1) = r**0.25d0 * dsin(phi/4)

      DexactU = TENSOR_SCALAR

      Return
      End



c ======================================================================
c  Exact (analytical) gradient of the solution
c ======================================================================
      Integer Function DgradU(x, y, label, dDATA, iDATA, iSYS, U)
      implicit none
      Include 'fem2Dtri.fd'

      Real*8    dDATA(*), x, y, U(MaxTensorSize, *)
      Integer   iDATA(*), label, iSYS(*)
        
      Integer   DexactU, i
      EXTERNAL  DExactU

      Real*8    delta, xD, yD, u0(1), uD(1)
      PARAMETER(delta = 1d-9)

c ==========================================================
      xD = x + delta

      i = DexactU(x,  y, label, dDATA, iDATA, iSYS, u0)
      i = DexactU(xD, y, label, dDATA, iDATA, iSYS, uD)

      U(1,1) = (uD(1) - u0(1)) / delta

      If(y.gt.0D0) Then
         yD = y + delta
         i  = DexactU(x, yD, label, dDATA, iDATA, iSYS, uD)
         U(2,1) = (uD(1) - u0(1)) / delta
      Else
         yD = y - delta
         i  = DexactU(x, yD, label, dDATA, iDATA, iSYS, uD)
         U(2,1) = -(uD(1) - u0(1)) / delta
      End if

      iSYS(1) = 2
      iSYS(2) = 1
      DgradU = TENSOR_GENERAL

      Return
      End





