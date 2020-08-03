c ======================================================================
      Program  mainBC
c ======================================================================
c This program generates a finite element system for the diffusion problem
c
c  -div K grad u + A u = F     in  Omega
c                    u = U_0   on  Gamma_D
c        K du/dn       = G_0   on  Gamma_N
c        K du/dn + S u = G_1   on  Gamma_R
c
c where Omega = [0,1]^2. The user-defined coefficients are  
c    K(x)   - positive definite tensor
c    A(x)   - non-negative reaction 
c    F(x)   - right-hand side
c    U_0(x) - essential (Dirichlet) boundary condition
c    G_0(x) - Neumann boundary condition
c    G_1(x) - Robin boundary condition
c    S(x)   - Robin boundary coefficient
c
c These coefficients are implemented in the following routines: 
c    K->Ddiff,  A->Dreact,  F->Drhs,  {U_0,G_0,G_1}->Dbc,  S->DbcRobCoef
c
c ======================================================================
      implicit none

      integer nvmax, ntmax, nbmax
      integer namax
c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
c namax - maximum number of non-zero matrix entries
      parameter(nvmax = 150 000, ntmax = 2*nvmax, nbmax = 10 000)
      parameter(namax = 90 000)

c working memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 5 000 000)

      Real*8   rW(MaxWr)
      Integer  iW(MaxWi)


C ======================================================================
c For library aniMBA
c ======================================================================
C group (M)
      Integer  nv, nvfix, labelV(nvmax), fixedV(nvmax)
      Real*8   vrt(2, nvmax)

      Integer  nb, nbfix, bnd(2, nbmax), labelB(nbmax), fixedB(nbmax)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2, nbmax)

      Integer  nt, ntfix, tri(3, ntmax), labelT(ntmax), fixedT(ntmax)

C group (CONTROL)
      Real*8   Metric(3, nvmax)


c ======================================================================
c For library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax), nRow, nCol
      Real*8    A(namax), RHS(nvmax)

      Integer  iDATAFEM(1), controlFEM(3)
      Real*8   dDATAFEM(1)

      EXTERNAL FEM2Dext


C LOCAL VARIABLES
      Integer  i

c =====================================================================
c ... load the initial mesh with 4 boundary edges. 
      Call loadMani(nv, nvfix, nvmax, vrt, labelV, fixedV,
     &              nb, nbfix, nbmax, bnd, labelB, fixedB,
     &              nc,               crv, labelC, 
     &              nt, ntfix, ntmax, tri, labelT, fixedT,
     &      "../data/simple.ani")


c ... Assemble the stifness matrix
c     no extra data is provided for the user subroutine Ddiff
      dDATAFEM(1) = 0D0
      iDATAFEM(1) = 0

c     mark the Dirichlet points on boundary edges 1 & 2
      labelV(1) = 1
      labelV(2) = 1
      labelV(4) = 1

      labelV(3) = 0
      labelV(5) = 0


c     general sparse matrix in the AMG format 
      controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_AMG)
      controlFEM(2) = 1

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &     FEM2Dext, dDATAFEM, iDATAFEM, controlFEM,
     &     nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &     MaxWi, MaxWr, iW, rW)


c ... save the matrix
      Open(1,file='CSRsystem')
        Write(1,*)  nRow, IA(nRow+1)-1
        Write(1,*) (IA(i), i=1,nRow+1)
        Write(1,*) (JA(i), i=1,IA(nRow+1)-1)
        Write(1,*) ( A(i), i=1,IA(nRow+1)-1)
        Write(1,*) (RHS(i),i=1,nRow)
      Close(1)


c ... draw the sparcity structure of A
      Call draw_matrix(nRow, IA, JA, "matrix_fin.ps")


c ... testing the results
      Do i = 1, nRow
         If(A(IA(i)).LE.0) Stop 911
      End do

      Stop
      End



C ======================================================================
C  The user defined routines required above
C ======================================================================
c Templated routine for the elemental matrix. It calls standard bilinear
c forms and imposes the boundary conditions using the provided labels.
C ======================================================================
      Subroutine FEM2Dext(XY1, XY2, XY3, 
     &                    lbT, lbB, lbV, dDATA, iDATA, iSYS,
     &                    LDA, A, F, nRow, nCol,
     &                    templateR, templateC)
C ======================================================================
      implicit none
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'

      Real*8  XY1(*), XY2(*), XY3(*)
      
      Integer lbT, lbB(3), lbV(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*), LDA

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Dreact, Drhs, Dbc, DbcRobCoef
      External Ddiff, Dreact, Drhs, Dbc, DbcRobCoef

      Integer  i,j,k,l,m, ir,ic, label, ibc, nRow, nCol
      Real*8   B(3, 3), C(3, 3), G(3), XYP(2, 3)
      Real*8   x, y, eBC(1)
      Logical  ifXbc
 
      Integer  iref(4)
      DATA     iref/1, 2, 3, 1/

C ======================================================================
      nRow = 3
      nCol = 3

c ... set up templated degrees of freedom for rows and columns. 
c     used convention 'V'=vertex d.o.f. and 'R'=edge d.o.f.
      Do i = 1, 3
         templateR(i) = Vdof
         templateC(i) = Vdof
      End do

c ... compute the stiffness matrix (M)
      label = lbT

c     A(1:3,1:3) is elemental vector elliptic operator;
c     in other words, for the bilinear form <grad(P1), grad(P1)>
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P1, GRAD, FEM_P1,
     &              label, Ddiff, dDATA, iDATA, iSYS, 2,
     &              LDA, A, ir, ic)

c     B(1:3,1:3) is elemental mass matrix;
c     in other words, for the bilinear form <P1, P1>
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P1, IDEN, FEM_P1,
     &              label, Dreact, dDATA, iDATA, iSYS, 4,
     &              3, B, ir, ic)

      Do i = 1, 3
         Do j = 1, 3
            A(i, j) = A(i, j) + B(i, j) 
         End do
      End do


c ... compute the right hand side vector using external function Drhs
c     in other words, the linear form <Drhs(x), P1> 
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P1, 
     &              lbT, Drhs, dDATA, iDATA, iSYS, 4,
     &              1, F, ir, ic)


c ... impose the Neumann and Robin boundary conditions
      Do i = 1, 2
         XYP(i, 1) = XY1(i)
         XYP(i, 2) = XY2(i)
         XYP(i, 3) = XY3(i)
      End do

      Do k = 1, 3
         If(lbB(k).GT.0) Then
            l = iref(k + 1)
            m = iref(l + 1)

            x = (XYP(1, k) + XYP(1, l)) / 2
            y = (XYP(2, k) + XYP(2, l)) / 2

            label = lbB(k)
            ibc = Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_NEUMANN) .OR. 
     &         ifXbc(ibc, BC_ROBIN)) Then
               label = lbB(k)

               Call fem2Dedge(XYP(1, k), XYP(1, l), XYP(1, m),
     &                        IDEN, FEM_P0, IDEN, FEM_P1, 
     &                        label, Dbc, dDATA, iDATA, iSYS, 4, 
     &                        1, G, ir, ic)

               F(k) = F(k) + G(1)
               F(l) = F(l) + G(2)
            End if

            If(ifXbc(ibc, BC_ROBIN)) Then
               label = lbB(k)

               Call fem2Dedge(XYP(1, k), XYP(1, l), XYP(1, m),
     &                        IDEN, FEM_P1, IDEN, FEM_P1,
     &                        label, DbcRobCoef, dDATA, iDATA, iSYS, 4,
     &                        3, C, ir, ic)

               A(k, k) = A(k, k) + C(1, 1)
               A(l, l) = A(l, l) + C(2, 2)
               A(k, l) = A(k, l) + C(1, 2)
               A(l, k) = A(l, k) + C(2, 1)
            End if

         End if
      End do

c ... impose Dirichlet boundary conditions at triangle nodes
c     this condition should go the last one
      Do k = 1, 3
         If(lbV(k).NE.0) Then
            x = XYP(1, k)
            y = XYP(2, k)

            label = lbV(k)
            ibc = Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_DIRICHLET)) Then
               Call applyDIR(LDA, nRow, A, F, k, eBC(1))
            End if
         End if
      End do

      Return
      End



C ======================================================================
c 2x2 diffusion tensor K
C ======================================================================
      Integer Function Ddiff(x, y, label, dDATA, iDATA, iSYS, Coef)
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 2

      Coef(1, 1) = 1D0
      Coef(2, 2) = 1D1
      Coef(1, 2) =-1D0
      Coef(2, 1) =-1D0

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
c Reaction coefficient A
C ======================================================================
      Integer Function Dreact(x, y, label, dDATA, iDATA, iSYS, Coef)
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = 1D0
      Dreact = TENSOR_SCALAR

      Return
      End



C ======================================================================
c Boundary conditions 
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, Coef)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      If(label.LE.2) Then        ! Dirichlet on sides labeled 1 & 2
         Coef(1, 1) = x + y      
         Dbc = BC_DIRICHLET 
      Else If(label.EQ.3) Then   ! Neumann on side labeled 3
         Coef(1, 1) = x - y
         Dbc = BC_NEUMANN    
      Else If(label.EQ.4) Then
         Coef(1, 1) = x
         Dbc = BC_ROBIN
      End if

      Return
      End



C ======================================================================
c Coefficient in Robin boundary condition
C ======================================================================
      Integer Function DbcRobCoef(x, y, label, dDATA, iDATA, iSYS, Coef)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      If(label.LE.3) Then
         DbcRobCoef = BC_NULL
      Else If(label.EQ.4) Then
         Coef(1, 1) = 1D0
         DbcRobCoef = BC_ROBIN_COEF
      End if

      Return
      End



C ======================================================================
c Right hand side F
C ======================================================================
      Integer Function Drhs(x, y, label, dDATA, iDATA, iSYS, Coef)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = 1D0
      Drhs = TENSOR_SCALAR

      Return
      End


