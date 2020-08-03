c ======================================================================
      Program  mainTemplate
c ======================================================================
c This program generates a finite element system for the Stokes problem
c
c  -div grad u  + grad p = 0   in Omega
c        div u           = 0   in Omega
c
c                      u = u_0 on dOmega_1
c                      u = 0   on dOmega_2
c              du/dn - p = 0   on dOmega_3
c
c where Omega is the unit square,
c dOmega_1 is the side at x=0, dOmega_3 is the side at x=1, 
c and dOmega_2 is the rest of the boundary. The non-homogeneous
c boundary condition is 
c
c    u_0 = { 4*y*(1-y), 0 }.
c
c We use the P2 finite elements for the velocity u and P1 finite 
c elements for the pressure p, which are known to be the stable 
c pair. The discretization method results in a symmetric indefinite 
c matrix.
c
c Step 1: we generate a quasi-uniform mesh using library libmba2D-2.x.a
c step 2: we generate  an algebraic system using library libfem2D-2.x.a
c ======================================================================
      implicit none

      integer nvmax,ntmax,nbmax,niamax,namax
c ... nvmax - maximum number of mesh nodes (v)
c ... ntmax - maximum number of mesh triangles (t)
c ... nbmax - maximum number of boundary edges (b)
      parameter(nvmax = 1000, ntmax = 2*nvmax, nbmax = 1000)

c ... niamax - maximum number of matrix rows
c ... namax  - maximum number of non-zero matrix entries
      parameter(niamax = nvmax*4 + nvmax*4 + nvmax)
      parameter(namax = 30*niamax )

c ... work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 100 000, MaxWi = 100 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)

c ======================================================================
c mesh description
c ======================================================================
C group (M)
      Integer  nv, nvfix, labelV(nvmax), fixedV(nvmax)
      Real*8   vrt(2, nvmax)

      Integer  nb, nbfix, bnd(2, nbmax), labelB(nbmax), fixedB(nbmax)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2, nbmax)

      Integer  nt, ntfix, tri(3, ntmax), labelT(ntmax), fixedT(ntmax)


c ======================================================================
c For library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(niamax), JA(namax), nRow, nCol
      Real*8    A(namax), RHS(niamax)

      Integer  iDATAFEM(1), iSYS(MAXiSYS), controlFEM(3), Dbc
      Real*8   dDATAFEM(1)

      EXTERNAL FEM2Dext, Dbc
  
      Logical  ifXbc


C LOCAL VARIABLES
      Integer  i, iERR
      Real*8   center, radius

c ==========================================================
c ... load the initial mesh with 4 boundary edges. 
      Call loadMani(nv, nvfix, nvmax, vrt, labelV, fixedV,
     &              nb, nbfix, nbmax, bnd, labelB, fixedB,
     &              nc,               crv, labelC, 
     &              nt, ntfix, ntmax, tri, labelT, fixedT,
     &              "../data/simple.ani")


      labelB(1) = 1  ! nonhomogeneous Dirichlet edge
      labelB(2) = 2  ! homogeneous Dirichlet edge
      labelB(3) = 3  ! homogeneous Neumann edge
      labelB(4) = 2  ! homogeneous Dirichlet edge


c ... generate the finite element system with P2-P1 elements
c     no data is provided for the user subroutines FEM2Dext and Dbc
      dDATAFEM(1) = 0D0
      iDATAFEM(1) = 0

c ... mark the Dirichlet points using labelV. The assembling routine
c     will pass these labels to elemental routines where user will 
c     set the Dirichlet boundary conditions.
      Call markDIR(nv,vrt, labelV, nb,bnd, labelB, 
     &             Dbc, dDATAFEM, iDATAFEM, iSYS)


c ... general sparce matrix in the AMG format 
      controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_AMG)
c ... verbosity level, from 0 to 9
      controlFEM(2) = 9

c ... we can call either of the two assembling routines below
      Call BilinearFormTemplate(  ! keeps some of zero entries
c     Call AssembleFromTemplate(  ! removes zero entries
     &     nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &     FEM2Dext, dDATAFEM, iDATAFEM, controlFEM,
     &     niamax, namax, IA, JA, A, RHS, nRow, nCol,
     &     MaxWi, MaxWr, iW, rW)


c ... save the matrix
      Open(1,file='CSRsystem',status='UNKNOWN')
         Write(1,*) nRow, IA(nRow+1)-1
         Write(1,*) (IA(i),i=1,nRow+1)  
         Write(1,*) (JA(i),i=IA(1),IA(nRow+1)-1)  
         Write(1,*) ( A(i),i=IA(1),IA(nRow+1)-1)  
         Write(1,*) (RHS(i),i=1,nRow)  
      Close(1)


c ... print information about the matrix
      Call gershgorin(nRow, IA, JA, A, center, radius)

c ... draw the sparsity structure of A
      Call draw_matrix(nRow, IA, JA, "matrix_fin.ps")
      Write(*,'(2(A,E10.3))') 
     &    '     the biggest Gershgorin circle: center =', 
     &    center, ', radius =', radius



      Stop 
      End



C ======================================================================
C  Here come user defined routines
C ======================================================================
c Templated routine for the elemental matrix. It calls standard bilinear
c forms and imposes the boundary conditions using the provided labels.
C ======================================================================
      Subroutine FEM2Dext( 
     &           XY1, XY2, XY3, 
     &           lbT, lbB, lbV, dDATA, iDATA, iSYS,
     &           LDA, A, F, nRow, nCol,
     &           templateR, templateC)
C ======================================================================
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'

      Real*8  XY1(*), XY2(*), XY3(*)
      
      Integer lbT, lbB(3), lbV(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  Ddiff, Drhs, Dbc, ANI_Dnull
      External Ddiff, Drhs, Dbc, ANI_Dnull

      Real*8   B(12, 12), x, y, eBC(2)
      Logical  ifXbc

C ======================================================================
      nRow = 15
      nCol = 15

c Populate degrees of freedom for rows and columns. 
c We use convention 'V'=vertex d.o.f. and 'R'=edge d.o.f.
c We breake unknowns into groups  to obtain a 5x5 block structure of the global matrix.
c The first, third and fifth diagonal blocks will have the same size equal to the number of mesh points.
c The remaining diagonal blocks will have size equal to the number of mesh edges.
c For the ordering of unknowns inside these blocks please refer to section "Assembling utilities" of the user_guide

      Do i = 1, 3
         templateR(i)     = Vdof + Group1  !u_x
         templateR(i + 3) = Rdof + Group2

         templateR(i + 6) = Vdof + Group3  !u_y
         templateR(i + 9) = Rdof + Group4

         templateR(i + 12)= Vdof + Group5  !p
      End do

      Do k = 1, nRow
         templateC(k) = templateR(k)
      End do

c ... compute the stiffness matrix (M)
      label = lbT

c     A(1:12,1:12) is elemental vector Laplacian matrix for velocity components
c     in other words, the bilinear form <grad(P2^2), grad(P2^2)>
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, FEM_P2vector, GRAD, FEM_P2vector,
     &              label, Ddiff, dDATA, iDATA, iSYS, 2,
     &              LDA, A, ir, ic)

c     B(1:3,1:12) is elemental divergence matrix for velocity
c     in other words, the bilinear form <div(P2^2), P1>
      Call fem2Dtri(XY1, XY2, XY3,
     &              DIV, FEM_P2vector, IDEN, FEM_P1,
     &              label, ANI_Dnull, dDATA, iDATA, iSYS, 2,
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

c ... compute the right hand side vector using external function Drhs
c     in other words, the linear form <Drhs(x), P2^2> 
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P2vector, 
     &              lbT, Drhs, dDATA, iDATA, iSYS, 2,
     &              1, F, ir, ic)

      Do i = 1, 3
         F(i + 12) = 0
      End do

c ... impose boundary conditions (assume nRow = nCol) for d.o.f. at triangle nodes
      Do k = 1, 3
         If(lbV(k).ne.0) Then
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

            label = lbV(k)
            ibc = Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_DIRICHLET)) Then
               Call applyDIR(LDA, nRow, A, F, k,   eBC(1))
               Call applyDIR(LDA, nRow, A, F, k+6, eBC(2))
            End if
         End if
      End do

c ... impose boundary conditions (assume nRow = nCol) at d.o.f. on triangle edges
      Do k = 1, 3
         If(lbB(k).GT.0) Then
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

            label = lbB(k)
            ibc = Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_DIRICHLET)) Then
               Call applyDIR(LDA, nRow, A, F, k+3, eBC(1))
               Call applyDIR(LDA, nRow, A, F, k+9, eBC(2))
            End if
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

      iSYS(1) = 4
      iSYS(2) = 4

      Do i = 1, 4
         Do j = 1, 4
            Coef(i, j) = 0
         End do
         Coef(i, i) = 1
      End do

      Ddiff = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C boundary conditions
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, Coef)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 1

      If(label.eq.1) Then
         Coef(1, 1) = y * (1 - y) * 4
         Coef(2, 1) = 0
         Dbc = BC_DIRICHLET

      Else If(label.eq.3) Then
         Coef(1, 1) = 0
         Coef(2, 1) = 0
         Dbc = BC_NEUMANN

      Else
         Coef(1, 1) = 0
         Coef(2, 1) = 0
         Dbc = BC_DIRICHLET
      End if

      Return
      End



C ======================================================================
C Right hand side
C ======================================================================
      Integer Function Drhs(x, y, label, dDATA, iDATA, iSYS, Coef)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 1

      Coef(1,1) = 0D0
      Coef(2,1) = 0D0

      Drhs = TENSOR_GENERAL

      Return
      End 






