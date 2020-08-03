C ======================================================================
      Subroutine FEM2Dext(XY1, XY2, XY3, 
     &           lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      Include 'fem2Dtri.fd'
C ======================================================================
      Real*8  XY1(*), XY2(*), XY3(*)
      
      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Drhs, Dbc
      External Ddiff, Drhs, Dbc

      Real*8   edge_length
      Real*8   M(3, 3), Q(3, 3), B(3), L(3), G(1), det, alpha, s
      Real*8   x, y, eBC(1)

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templates (we have 1 group, skipping additional bits)

      Do i = 1, 3
         templateR(i) = Rdof ! + Group1
      End do

      Do i = 1, 3
         templateC(i) = templateR(i)
      End do

c ... compute the stiffness matrix (M)
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_RT0, IDEN, FEM_RT0,
     &              label, Ddiff, dDATA, iDATA, iSYS, 2,
     &              3, M, ir, ic)

c ... compute right hand side
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P0,
     &              lbE, Drhs, dDATA, iDATA, iSYS, 2,
     &              1, G, ir, ic)

c ... compute the constraint matrix
      B(1) = edge_length(XY1, XY2)
      B(2) = edge_length(XY2, XY3)
      B(3) = edge_length(XY3, XY1)

c ... elliminate fluxes
      Call invert3x3(M, Q, det)

      Do i = 1, 3
         Do j = 1, 3
            Q(i, j) = Q(i, j) * B(j)
         End do
      End do

      Do i = 1, 3
         L(i) = Q(1, i) * B(1) + Q(2, i) * B(2) + Q(3, i) * B(3) 
      End do

      alpha = L(1) + L(2) + L(3)

      Do i = 1, 3
         Do j = i, 3
            A(i, j) = B(i) * Q(i, j) - L(i) * L(j) / alpha
            A(j, i) = A(i, j)
         End do
         F(i) = L(i) * G(1) / alpha
      End do

c ... impose boundary conditions (assume nRow = nCol)
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

            i = Dbc(x, y, lbF(k), dDATA, iDATA, iSYS, eBC)

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
C Identity tensor
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
         Coef(i, i) = 1D0
      End do

c     The hybrid method requires inverse of the diffusion tensor.
c     Here we are lucky, the inverse of identity is identity:
C     Call invertSymmetricMatrix(2, 4, Coef) 

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C Essential boundary conditions
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, eBC(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      eBC(1, 1) = 0D0 
      Dbc = BC_DIRICHLET

      Return
      End



C ======================================================================
C Sourse data
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


