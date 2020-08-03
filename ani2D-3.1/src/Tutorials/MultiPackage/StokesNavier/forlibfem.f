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

      Integer LDA, nRow, nCol, templateR(*), templateC(*)
      Real*8  A(LDA, *), F(*)

C Local variables
      Integer  Ddiff, Dconv, Drhs, Dbc, ANI_Dnull
      External Ddiff, Dconv, Drhs, Dbc, ANI_Dnull

      Integer  i,j,k, ibc, ir,ic, label, iE1,iE2,iE3
      Integer  iDATANUL(1), iE, iP1,iP2,iP3, iR1,iR2,iR3, nP,nR,nE
      Real*8   dDATANUL(1), dDATACONV(18), B(12, 12), x,y, eBC(2)
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

         templateR(i + 12) = Vdof + Group3 !p
      End do

c     same degrees of freedom for columns (symmetric sparsity structure)
      Do k = 1, nRow
         templateC(k) = templateR(k)
      End do


c ... preprocess data arrays
      Call decodeISYS(iE, iP1,iP2,iP3, iR1,iR2,iR3, 
     &                iE1,iE2,iE3, nP,nR,nE, iSYS) 

      dDATACONV( 1) = XY1(1)
      dDATACONV( 2) = XY2(1)
      dDATACONV( 3) = XY3(1)
      dDATACONV( 4) = XY1(2)
      dDATACONV( 5) = XY2(2)
      dDATACONV( 6) = XY3(2)

      dDATACONV( 7) = dDATA(iP1)
      dDATACONV( 8) = dDATA(iP2)
      dDATACONV( 9) = dDATA(iP3)
      dDATACONV(10) = dDATA(3*nP + iR1)
      dDATACONV(11) = dDATA(3*nP + iR2)
      dDATACONV(12) = dDATA(3*nP + iR3)

      dDATACONV(13) = dDATA(nP + iP1)
      dDATACONV(14) = dDATA(nP + iP2)
      dDATACONV(15) = dDATA(nP + iP3)
      dDATACONV(16) = dDATA(3*nP + nR + iR1)
      dDATACONV(17) = dDATA(3*nP + nR + iR2)
      dDATACONV(18) = dDATA(3*nP + nR + iR3)
     

c ... compute the convection-diffusion matrix (M)
      label = lbE

c     A(1:12,1:12) is elemental vector stiffness matrix for velocity components
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P2vector, GRAD, FEM_P2vector,
     &              label, Ddiff, dDATANUL, iDATANUL, iSYS, 2,
     &              LDA, A, ir, ic)

c     B(1:12,1:12) is elemental vector convection matrix for velocity components
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P2vector, GRAD, FEM_P2vector,
     &              label, Dconv, dDATACONV, iDATANUL, iSYS, 2,
     &              12, B, ir, ic)

c     summing diffusion and convection matrices
      Do i = 1, 12
         Do j = 1, 12
            A(i, j) = A(i, j) + B(i, j)
         End do
      End do


c     B(1:3,1:12) is elemental divergence matrix for velocity
      Call fem2Dtri(XY1, XY2, XY3,
     &              DIV, FEM_P2vector, IDEN, FEM_P1,
     &              label, ANI_Dnull, dDATANUL, iDATANUL, iSYS, 2,
     &              12, B, ir, ic)

      Do i = 1, 12
         Do j = 1, 3
            A(i, j + 12) = B(i, j) 
            A(j + 12, i) = B(i, j)
         End do
      End do

      Do i = 1, 3
         Do j = 1, 3
            A(i + 12, j + 12) = 0
         End do
      End do

c ... compute right hand side vector
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P2vector, 
     &              lbE, Drhs, dDATA, iDATA, iSYS, 2,
     &              1, F, ir, ic)

      Do i = 1, 3
         F(i + 12) = 0D0
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
C Diffusion tensor = viscosity * Identity tensor
C ======================================================================
      Integer Function Ddiff(x, y, label, dDATA, iDATA, iSYS, Coef)
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)
 
      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1,1) = 3D-4

      Ddiff = TENSOR_SCALAR

      Return
      End



C ======================================================================
C Dirichlet boundary conditions
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)

      Include 'fem2Dtri.fd'
      Real*8  dDATA(*), x, y, eBC(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 1

      If(label.eq.1) Then
         eBC(1,1) = 1
         eBC(2,1) = 0
         Dbc = BC_DIRICHLET

      Else
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = BC_DIRICHLET
      End if

      Return
      End



C ======================================================================
C  Convection tensor
C ======================================================================
      Integer Function Dconv(x, y, label, dDATA, iDATA, iSYS, Conv)
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Conv(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      Real*8  s1,s2,s3, s12,s13,s23, x1,x2,x3, y1,y2,y3
      Real*8  d1,d2,d3,d, lam1,lam2,lam3, g1,g2,g3
C =====================================================================
      iSYS(1) = 4
      iSYS(2) = 2

      Do i = 1, 4
         Do j = 1, 2
            Conv(i, j) = 0
         End do
      End do

      x1 = dDATA(1)
      x2 = dDATA(2)
      x3 = dDATA(3)
      y1 = dDATA(4)
      y2 = dDATA(5)
      y3 = dDATA(6)

      d1 = (x2-x)*(y3-y)-(y2-y)*(x3-x)
      d2 = (x3-x)*(y1-y)-(y3-y)*(x1-x)
      d3 = (x1-x)*(y2-y)-(y1-y)*(x2-x)
      d = d1 + d2 + d3

c     barycentric coordinates of x,y
      lam1 = d1 / d
      lam2 = d2 / d
      lam3 = d3 / d

c quadratic interpolation of ux
      s1 = dDATA( 7)
      s2 = dDATA( 8)
      s3 = dDATA( 9)
      s12= dDATA(10)
      s23= dDATA(11)
      s13= dDATA(12)

      Conv(1,1) = s1*lam1*(2*lam1-1d0) + s2*lam2*(2*lam2-1d0)
     &          + s3*lam3*(2*lam3-1d0)
     &          + s12*4*lam1*lam2 + s13*4*lam1*lam3 + s23*4*lam2*lam3
      Conv(3,2) = Conv(1,1)

c quadratic interpolation of uy
      s1 = dDATA(13)
      s2 = dDATA(14)
      s3 = dDATA(15)
      s12= dDATA(16)
      s23= dDATA(17)
      s13= dDATA(18)

      Conv(2,1) = s1*lam1*(2*lam1-1d0) + s2*lam2*(2*lam2-1d0)
     &          + s3*lam3*(2*lam3-1d0)
     &          + s12*4*lam1*lam2 + s13*4*lam1*lam3 + s23*4*lam2*lam3
      Conv(4,2) = Conv(2,1)

      Dconv = TENSOR_GENERAL

      Return
      End



C ======================================================================
C Right hand side
C ======================================================================
      Integer Function Drhs(x, y, label, dDATA, iDATA, iSYS, F)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, F(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 1
     
      F(1,1) = 0D0
      F(2,1) = 0D0

      Drhs = TENSOR_GENERAL

      Return
      End 



C ======================================================================
C Templated routine for the P1 mass elemental matrix. It calls standard 
C bilinear forms and imposes no boundary conditions.
C ======================================================================
      Subroutine FEM2DPrj1ext(XY1, XY2, XY3,
     &                    lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &                    LDA, A, F, nRow, nCol,
     &                    templateR, templateC)
C ======================================================================
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'

      Real*8  XY1(*), XY2(*), XY3(*)

      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  ANI_Dnull
      External ANI_Dnull
C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templated degrees of freedom for rows and columns.
      Do i = 1, 3
         templateR(i) = Vdof
         templateC(i) = Vdof
      End do

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:3,1:3) is elemental identity operator;
c     in other words, for the bilinear form <iden(P1), iden(P1)>
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P1, IDEN, FEM_P1,
     &              label, ANI_Dnull, dDATA, iDATA, iSYS, 2,
     &              LDA, A, ir, ic)

      Return
      End



C ======================================================================
C Templated routine for the P2 mass elemental matrix. It calls standard 
C bilinear forms and imposes no boundary conditions.
C ======================================================================
      Subroutine FEM2DPrj2ext(XY1, XY2, XY3,
     &                    lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &                    LDA, A, F, nRow, nCol,
     &                    templateR, templateC)
C ======================================================================
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'

      Real*8  XY1(*), XY2(*), XY3(*)

      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  ANI_Dnull
      External ANI_Dnull
C ======================================================================
      nRow = 6
      nCol = 6

c ... set up templated degrees of freedom for rows and columns.
      Do i = 1, 3
         templateR(i)     = Vdof
         templateR(i + 3) = Rdof
      End do

      Do i = 1, 6
          templateC(i) = templateR(i)
      End do 

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:6,1:6) is elemental identity operator;
c     in other words, for the bilinear form <iden(P2), iden(P2)>
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P2, IDEN, FEM_P2,
     &              label, ANI_Dnull, dDATA, iDATA, iSYS, 4,
     &              LDA, A, ir, ic)

      Return
      End

