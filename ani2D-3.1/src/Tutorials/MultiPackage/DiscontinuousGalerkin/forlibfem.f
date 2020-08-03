C ======================================================================
      Subroutine FEM2DextDG(XY1, XY2, XY3, 
     &                      lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &                      LDA, A, F, nRow, nCol,
     &                      templateR, templateC)
C ======================================================================
      Implicit none
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'
C ======================================================================
      Real*8  XY1(2,2), XY2(2,2), XY3(2,2)
      
      Integer lbE, lbF(3,3), lbP(3,2)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*), LDA, nRow, nCol

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  ANI_Dnull, ANI_D1x1_scalar, Drhs, Dbc
      External ANI_Dnull, ANI_D1x1_scalar, Drhs, Dbc

      Integer  i,j,k,m,n, i1,i2,i3, ir,ic, label, iE, iEt
      Integer  FEMtypeA, FEMtypeB(3), IEE(3), ks(3)

      Real*8   B(40,40), XYZ(2,3,2), dJUMP(1), eBC(1)
      Real*8   x1,y1, x2,y2, x3,y3, x,y, d

      Integer  iref(4)
      DATA     iref/1,2,3,1/

C ======================================================================
c ... required call to pre-process iSYS
      Call DG_init(XY1, XY2, XY3, iSYS, XYZ, iE, IEE) 

c ... set up templates
      FEMtypeA = iDATA(iE)
      Do i = 1, 3
         If(IEE(i).GT.0) FEMtypeB(i) = iDATA(IEE(i))
      End do
      Call DG_template(iE,IEE, FEMtypeA,FEMtypeB, nRow, templateR, ks) 

      nCol = nRow
      Do i = 1, nRow
         templateC(i) = templateR(i)
      End do
      If(iSYS(1).LT.0) Return


c ... compute the stiffness matrix A
c ... add jumps of edges
      dJUMP(1) = 1D2
      Call DGjump_AJ(XYZ, iE, IEE, ks,
     &               FEMtypeA, FEMtypeB, 
     &               ANI_D1x1_scalar, dJUMP, iDATA, iSYS, 9,
     &               LDA, A, nRow, nCol)

c     Call DGjump_SIP(XYZ, iE, IEE, ks,
      Call DGjump_NIP(XYZ, iE, IEE, ks,
     &                FEMtypeA, FEMtypeB, 
     &                ANI_Dnull, dDATA, iDATA, iSYS, 9,
     &                40, B, nRow, nCol)

      Do i = 1, nRow
         Do j = 1, nCol
            A(i, j) = A(i, j) + B(i, j)
         End do
      End do


c ... add stiffness matrix for the central element
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEMtypeA, GRAD, FEMtypeA,
     &              label, ANI_Dnull, dDATA, iDATA, iSYS, 9,
     &              40, B, ir, ic)

      Do i = 1, ir
         Do j = 1, ic
            A(i, j) = A(i, j) + B(i, j)
         End do
      End do


c ... compute the right-hand side
      Do i = 1, nRow
         F(i) = 0D0
      End do

      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEMtypeA, 
     &              label, Drhs, dDATA, iDATA, iSYS, 9,
     &              1, F, ir, ic)


c ... apply boundary conditions on central triangle
c     n=0 for P1, n=1 for P2, n=2 for P3
      Do i1 = 1, 3
         x1 = XYZ(1,i1,1)
         y1 = XYZ(2,i1,1)

         label = lbP(i1, 1)
         If(label.GT.0) Then
            k = Dbc(x1, y1, label, dDATA, iDATA, iSYS, eBC)
            Call applyDIR(LDA, nRow, A, F, i1, eBC(1))
         End if

         n = 0
         If(iDATA(iE).EQ.FEM_P2) n = 1
         If(iDATA(iE).EQ.FEM_P3) n = 2

         label = lbF(i1, 1)
         If(label.GT.0) Then
            i2 = iref(i1 + 1)

            x2 = XYZ(1,i2,1)
            y2 = XYZ(2,i2,1)

            Do i = 1, n
               d = dble(i) / (n + 1)
               m = i1 + 3*i

               x = x1 * (1D0-d) + x2 * d
               y = y1 * (1D0-d) + y2 * d

               k = Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
               Call applyDIR(LDA, nRow, A, F, m, eBC(1))
            End do
         End if
      End do

c ... apply boundary conditions on neigboring triangles
      Do i1 = 1, 3
         iEt = IEE(i1)

         If(iEt.GT.0) Then
            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            label = lbP(i1, 1)
            If(label.GT.0) Then
               x1 = XYZ(1,i1,1)
               y1 = XYZ(2,i1,1)

               k = Dbc(x1, y1, label, dDATA, iDATA, iSYS, eBC)
               Call applyDIR(LDA, nRow, A, F, ks(i1)+1, eBC(1))
            End if

            label = lbP(i2, 1)
            If(label.GT.0) Then
               x2 = XYZ(1,i2,1)
               y2 = XYZ(2,i2,1)

               k = Dbc(x2, y2, label, dDATA, iDATA, iSYS, eBC)
               Call applyDIR(LDA, nRow, A, F, ks(i1)+2, eBC(1))
            End if

            label = lbP(i3, 2)
            If(label.GT.0) Then
               x3 = XYZ(1,i3,2)
               y3 = XYZ(2,i3,2)

               k = Dbc(x3, y3, label, dDATA, iDATA, iSYS, eBC)
               Call applyDIR(LDA, nRow, A, F, ks(i1)+3, eBC(1))
            End if

c           xi,yi are already defined for the Dirichlet problem
            n = 0
            If(iDATA(iEt).EQ.FEM_P2) n = 1
            If(iDATA(iEt).EQ.FEM_P3) n = 2

            label = lbF(i1, 2)
            If(label.GT.0) Then
               Do i = 1, n
                  d = dble(i) / (n + 1)
                  m = ks(i1) + 3*(i-1) + 5 

                  x = x2 * (1D0-d) + x3 * d
                  y = y2 * (1D0-d) + y3 * d

                  k = Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
                  Call applyDIR(LDA, nRow, A, F, m, eBC(1))
               End do
            End if

            label = lbF(i1, 3)
            If(label.GT.0) Then
               Do i = 1, n
                  d = dble(i) / (n + 1)
                  m = ks(i1) + 3*(i-1) + 6

                  x = x3 * (1D0-d) + x1 * d
                  y = y3 * (1D0-d) + y1 * d

                  k = Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
                  Call applyDIR(LDA, nRow, A, F, m, eBC(1))
               End do
            End if
         End if
      End do

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

      F(1, 1) = -2D0
      Drhs = TENSOR_SCALAR

      Return
      End



C ======================================================================
C Essential boundary conditions
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, eBC(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*), label

      iSYS(1) = 1
      iSYS(2) = 1

      eBC(1,1) = x*x
      Dbc = BC_DIRICHLET

      Return
      End



C ======================================================================
C Exact solution
C ======================================================================
      Integer Function Dexact(x, y, label, dDATA, iDATA, iSYS, U)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, U(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      U(1, 1) = x*x
      Dexact = TENSOR_SCALAR

      Return
      End





