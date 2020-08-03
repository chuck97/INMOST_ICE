c ======================================================================
      Program Stokes
c ======================================================================
c The adaptive solution for the following bvp:
c
c  -div grad u  + grad p = 0   in Omega
c        div u           = 0   in Omega
c
c                      u = u_0 on dOmega_1
c                      u = 0   on dOmega_2
c                  du/dn = 0   on dOmega_3
c
c where Omega is a domain with a circle hole,  dOmega_1 is the side 
c at x=-2, dOmega_3 is the side at x=2, and dOmega_2 is the rest of 
c the boundary. The non-homogeneous boundary condition is 
c
c    u_0 = { 1-4y^2), 0 }.
c
c We use the P2 finite elements for the velocity u and P1 finite elements 
c for the pressure p, which are known to be stable. The discretization 
c method results in a symmetric indefinite matrix.
c
c ======================================================================
      implicit none

c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
c namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,nfmax,namax
      parameter(nvmax =  50 000, ntmax = 2*nvmax, nbmax = 5 000)
      parameter(nfmax = 200 000, namax = 5 000 000)

c work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 3 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


c ======================================================================
c Mesh definition
c ======================================================================
      Integer  nv, nvfix, labelV(nvmax), fixedV(6)
      Real*8   vrt(2,nvmax)

      Integer  nb, nbfix, bnd(2,nbmax), labelB(nbmax), fixedB(1)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2,nbmax)

      Integer  nt, ntfix, tri(3,ntmax), labelT(ntmax), fixedT(1)

      DATA     nvfix/6/, fixedV/1,2,3,4,5,6/,  nbfix/0/, ntfix/0/


c ======================================================================
c for library aniAFT 
c ======================================================================
      Real*8   bv(2,6), bltail(2,6)
      Integer  Nbv, Nbl, bl(7,6)

      Data     Nbv/6/, Nbl/6/
      Data     bv/2.0,-1.5,  2.0,1.5,  -2.0,0.5,  -2.0,-0.5, 
     &            0.5,0.0,  -0.5,0.0/
      Data     bl/2,1,0,-1,1,1,0,    1,4,2,-1,2,1,0,      ! four exterior boundary edges
     &            4,3,0,-1,3,1,0,    3,2,1,-1,4,1,0,
     &            5,6,3,-1,5,1,0,    6,5,4,-1,6,1,0/      ! two boundary edges for circle
      Data     bltail/0.0,0.0, 1.0,0.0, 0.0,0.0, 0.0,1.0, ! data for curved boundary edges
     &                0.0,1.0, 0.0,1.0/

      Integer   aft2dboundary, userboundary
      EXTERNAL  aft2dboundary, userboundary


c ======================================================================
c for library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nfmax), JA(namax)
      Real*8    A(namax), RHS(nfmax), SOL(nfmax), RES(nfmax)

      Integer  iDATA(30), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATA(1)

      Integer  Dbc
      EXTERNAL Dbc, FEM2Dext


c ======================================================================
c for library aniLU
c ======================================================================
      Integer  symbolic(2), numeric(2), sys
      Real*8   lucontrol(20), luinfo(90)


c ======================================================================
c for library aniLMR
c ======================================================================
      Real*8   Lp
      Real*8   Metric(3,nvmax)


c ======================================================================
c for library aniMBA
c ======================================================================
      Integer  control(6), nEStar
      Real*8   Quality
      EXTERNAL CrvFunction_user


c LOCAL VARIABLEs
      Integer  i, iLoop, nLOOPs, nRow, nCol, iERR
      Integer  nr, iux,iuy, physics, ip, ibc, iv1,iv2
      Real*8   rmax, h, x, y, eBC(2)
      Character*30 file_mesh, file_velocity, file_pressure

c ======================================================================
c ... load the parameters
      Open(10, file="aniPY_stokes.txt")
      Read(10, *) nLOOPs  ! number of adaptive loops
      Read(10, *) nEStar
      Read(10, *) Quality
      Read(10, *) physics
      Do i = 1, 7
         Read(10, *) iDATA(i)
      End do
      Close(10)



c ... pass the user function userBoundary (forlibaft.c) 
      Call registeruserfn( userboundary ) ! register the name in the


c === generate a quasi-uniform mesh with mesh step h
      h = 0.196                          
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary(
     &       Nbv, bv, Nbl, bl, bltail,  h,  ! input geometric data
     &       nv, vrt,                       ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)

      If(iERR.NE.0) Stop 'Error in function aft2dboundary'
      If(nv.GT.nvmax) stop ' too many nodes, increase nvmax'
      If(nt.GT.ntmax) stop ' too many triangles, increase ntmax'
      If(nb.GT.nbmax) stop ' too many boundary edges, increase nbmax'

      Write(*,'(A,3(I6,A))') '     mesh with', nt, ' triangles', 
     &                       nv, ' vertices ', nb, ' boundary edges' 
          

c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c === integer data are provided for the user subroutines 
         dDATA(1) = 0D0

c mark the Dirichlet points with the maximal edge color
         Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &                Dbc, dDATA, iDATA, iSYS)


c === general sparse matrix in a 0-based CSC format used in UMFPACK
         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
         controlFEM(2) = 9

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2Dext, dDATA, iDATA, controlFEM,
     &        nfmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c === call the driver for LU factorization and solution
         Call CSC2CSC0(nCol, IA, JA)

c set up default control parameters & print only error messages
         Call umf4def(lucontrol)
         lucontrol(1) = 1

c pre-order and symbolic analysis
         Call umf4sym(nCol, nCol, IA,JA,A, symbolic,lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5001

c numeric factorization
         Call umf4num(IA,JA,A, symbolic,numeric,lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5002

c free the symbolic analysis data
         Call umf4fsym(symbolic)

c solve Ax=b, without iterative refinement
         sys = 0
         Call umf4sol(sys, SOL, RHS, numeric, lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5003

c free the numeric factorization data
         Call umf4fnum (numeric)


c check the residual
         Call mulAcsc0(nRow, IA, JA, A, SOL, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(A,E12.6)') 'LU:  maximal norm of residual: ', rmax


c === unknowns in SOL: nodal x-component of velocity starts at SOL(iux),
c                      nodal y-component of velocity starts at SOL(iuy),
c                      nodal pressure starts at SOL(ip)
c                      mid-edge x-component of velocity starts at SOL(iuxe)
c                      mid-edge y-component of velocity starts at SOL(iuye)
c        iux = 1
c        iuy = iux + nv
c        ip  = iuy + nv
c        iuxe= ip  + nv
c        iuye= iuxe+ nr


c === draw mesh and solution isolines as ps-figures (mid-edges are ignored)
         nr  = (nRow - 3*nv) / 2

         iux = 1
         iuy = iux + nv
         ip  = iuy + nv

c  PostScript file names must have extension .ps
c  demo graphics has been activated
c  Isolines of velocity components can be done with the following calls:
c  Call isolines(SOL(iux), nv,vrt, nt,tri, nb,bnd,'velocity_x.ps',20)
c  Call isolines(SOL(iuy), nv,vrt, nt,tri, nb,bnd,'velocity_y.ps',20)
         If(iLoop.EQ.1) Then
            file_mesh     = 'mesh_initial.ps'
            file_velocity = 'streamlines_ini.ps'
            file_pressure = 'pressure_ini.ps'
         Else
            file_mesh     = 'mesh_final.ps'
            file_velocity = 'streamlines_fin.ps'
            file_pressure = 'pressure_fin.ps'
         End if

         Call stream_function(nv, nr, nb, nt, vrt, bnd, tri,
     &                        SOL(iux), rW, MaxWi, iW)

c        Call isolines(rW,     nv,vrt, nt,tri, nb,bnd,file_velocity, 20)
c        Call isolines(SOL(ip),nv,vrt, nt,tri, nb,bnd,file_pressure, 40)
c        Call graph(           nv,vrt, nt,tri,        file_mesh)
         If(iLoop.EQ.1) Then
            Call isolines_demo(rW, nv,vrt, nt,tri, nb,bnd,file_velocity,
     &           20, 'Velocity streamlines on the initial mesh')
            Call isolines_demo(SOL(ip),nv,vrt, nt,tri, nb,bnd,
     &           file_pressure,
     &           40, 'Pressure isolines on the initial mesh')
            Call graph_demo(nv,vrt, nt,tri, file_mesh, 
     &           'Initial quasi-uniform mesh')
            Call draw_matrix_tol(nRow, IA, JA, A, 0D0, "matrix_ini.ps")
         Else
            Call isolines_demo(rW, nv,vrt, nt,tri, nb,bnd,file_velocity,
     &           20, 'Velocity streamlines on the final mesh')
            Call isolines_demo(SOL(ip),nv,vrt, nt,tri, nb,bnd,
     &           file_pressure,
     &           40, 'Pressure isolines on the final mesh')
            Call graph_demo(nv,vrt, nt,tri, file_mesh, 
     &           'Final mesh adapted to x-component of velocity')
            Call draw_matrix_tol(nRow, IA, JA, A, 0D0, "matrix_fin.ps")
         End if

         If(iLoop.eq.nLOOPs) Goto 500


c === generate metric (from SOL) optimal for the L_p norm
c        Lp = 0             ! maximum norm
         Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(SOL(physics * nv + 1), 
     &                        vrt, nv, tri, nt, bnd, nb, Metric,
     &                        MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c === generate the adaptive mesh to u_x
         control(1) = 100     ! MaxSkipE
         control(2) = 10000   ! MaxQItr
         control(3) = 1       ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1       ! flagAuto
         control(5) = 1       ! iPrint = minimal level of output information
         control(6) = 0       ! iErrMesgt: only critical termination allowed

         Call mbaNodal(
c group (M)
     &        nv, nvfix, nvmax, vrt, labelV, fixedV,
     &        nb, nbfix, nbmax, bnd, labelB, fixedB,
     &        nc,               crv, labelC, CrvFunction_user,
     &        nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &        nEStar, Quality, control, Metric,
c group (W)
     &        MaxWr, MaxWi, rW, iW, iERR)
 
         If(iERR.GT.1000) Call errMesMBA(iERR, 'main',
     &                        'unspecified error if mbaNodal')

 500     Continue
      End do

      Stop 


c ... Error messages
 5001 Continue
      Write(*,*) 'Error occurred in umf4sym: ', luinfo(1)
      Stop 

 5002 Continue
      Write(*,*) 'Error occurred in umf4num: ', luinfo(1)
      Stop

 5003 Continue
      Write(*,*) 'Error occurred in umf4sol: ', luinfo(1)

      Stop
      End



C ======================================================================
      Subroutine FEM2Dext(XY1, XY2, XY3, 
     &                    lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &                    LDA, A, F, nRow, nCol,
     &                    templateR, templateC)
C ======================================================================
      implicit none
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
      Integer  Ddiff, Drhs, Dbc, ANI_Dnull, femu, femp 
      External Ddiff, Drhs, Dbc, ANI_Dnull

      Real*8   B(12, 12), XYP(2, 3), XYN(2), Atmp(11, 11)
      Real*8   x,y, x1,y1, x2,y2, s, eBC(2), eBC1(2), eBC2(2)

      Integer  i,j,k,l,m, nu,np, iuy, label, ibc, ir,ic
      Integer  iref(4), map(11)
      DATA     iref/1, 2, 3, 1/

      Logical  ifXbc

C ======================================================================
      femu = iDATA(6)
      femp = iDATA(7)  

      If(femu.EQ.FEM_P2vector) Then
         nu  = 12
         iuy = 6
         Do i = 1, 3
            templateR(i)     = Vdof + Group1 !u_x
            templateR(i + 3) = Rdof + Group5
            templateR(i + 6) = Vdof + Group2 !u_y
            templateR(i + 9) = Rdof + Group6
         End do

      Else If(femu.EQ.FEM_P2reduced) Then
         nu  = 9
         iuy = 6
         Do i = 1, 3
            templateR(i)     = Vdof + Group1       !u_x
            templateR(i + 3) = RdofOrient + Group5 !bubbles
            templateR(i + 6) = Vdof + Group2       !u_y
         End do

      Else If(femu.EQ.FEM_MINI) Then
         Do i = 1, 11
            map(i) = i
         End do

         nu  = 8
         iuy = 4
         Do i = 1, 3
            templateR(i)     = Vdof + Group1 !u_x
            templateR(i + 4) = Vdof + Group2 !u_y
            map(i + 3) = i + 4
         End do
         templateR(4) = Edof + Group4  !u_x
         templateR(8) = Edof + Group4  !u_y
         map(7) = 4
      End If

      If(femp.EQ.FEM_P1) Then
         np = 3
         Do i = 1, 3
            templateR(i + nu) = Vdof + Group3
         End do

      Else If(femp.EQ.FEM_P0) Then
         np = 1
         templateR(nu + 1) = Edof + Group4
      End if
 
      nRow = nu + np
      nCol = nRow

      Do k = 1, nRow
         templateC(k) = templateR(k)
      End do

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:12,1:12) is elemental vector Laplacian matrix for velocity components
      Call fem2Dtri(XY1, XY2, XY3,
     &              GRAD, femu, GRAD, femu,
     &              label, Ddiff, dDATA, iDATA, iSYS, 4,
     &              LDA, A, ir, ic)


c     B(1:3,1:12) is elemental divergence matrix for velocity
      Call fem2Dtri(XY1, XY2, XY3,
     &              DIV, femu, IDEN, femp,
     &              label, ANI_Dnull, dDATA, iDATA, iSYS, 4,
     &              12, B, ir, ic)

      Do i = 1, nu
         Do j = 1, np
            A(i, j + nu) = B(i, j)
            A(j + nu, i) = B(i, j) 
         End do
      End do

      Do i = 1, np
         Do j = 1, np
            A(i + nu, j + nu) = 0D0
         End do
      End do

c ... compute right hand side vector
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, femu, 
     &              lbE, Drhs, dDATA, iDATA, iSYS, 2,
     &              1, F, ir, ic)

      Do i = 1, np
         F(i + nu) = 0
      End do

c ... impose boundary conditions at nodes of triangle
      Do i = 1, 2
         XYP(i, 1) = XY1(i)
         XYP(i, 2) = XY2(i)
         XYP(i, 3) = XY3(i)
      End do

      Do k = 1, 3
         If(lbP(k).ne.0) Then
            x = XYP(1, k)
            y = XYP(2, k)

            ibc = Dbc(x, y, lbP(k), dDATA, iDATA, iSYS, eBC)

            If(ifXbc(ibc, BC_DIRICHLET)) Then
               Call applyDIR(LDA, nRow, A, F, k,     eBC(1))
               Call applyDIR(LDA, nRow, A, F, k+iuy, eBC(2))
            End if
         End if
      End do

c ... impose boundary conditions at middles of edges
      If(femu.EQ.FEM_P2vector .OR. femu.EQ.FEM_P2reduced) Then
         Do k = 1, 3
            If(lbF(k).GT.0) Then
               l = iref(k + 1)
               m = iref(l + 1)

               x1 = XYP(1, k)
               y1 = XYP(2, k)

               x2 = XYP(1, l)
               y2 = XYP(2, l)

               x = (x1 + x2) / 2
               y = (y1 + y2) / 2
               ibc = Dbc(x, y, lbF(k), dDATA, iDATA, iSYS, eBC)

               If(ifXbc(ibc, BC_DIRICHLET) .AND. 
     &            femu.EQ.FEM_P2reduced) Then
                  Call calNormalExt(XYP(1,k), XYP(1,l), XYP(1,m), XYN)
   
                  ibc = Dbc(x1, y1, lbF(k), dDATA, iDATA, iSYS, eBC1)
                  ibc = Dbc(x2, y2, lbF(k), dDATA, iDATA, iSYS, eBC2)

                  eBC(1) = eBC(1) - (eBC1(1) + eBC2(1)) / 2
                  eBC(2) = eBC(2) - (eBC1(2) + eBC2(2)) / 2

                  s = eBC(1) * XYN(1) + eBC(2) * XYN(2)
                  Call applyDIR(LDA, nRow, A, F, k+3, s)

               Else If(ifXbc(ibc, BC_DIRICHLET) .AND. 
     &                 femu.EQ.FEM_P2vector) Then
                  Call applyDIR(LDA, nRow, A, F, k+3, eBC(1))
                  Call applyDIR(LDA, nRow, A, F, k+9, eBC(2))
               End if
            End if
         End do
      End if

c     reorder the first matrix block to comply with the ordering rules
      If(femu.EQ.FEM_MINI) Then
         Do i = 1, nRow
            Do j = 1, nRow
               Atmp(i, j) = A(i, j)
            End do
         End do

         Do i = 1, nRow
            Do j = 1, nRow
               A(i, j) = Atmp(map(i), map(j))
            End do
         End do

         Do i = 1, nRow
            Atmp(i, 1) = F(i)
         End do

         Do i = 1, nRow
            F(i) = Atmp(map(i), 1) 
            templateC(i) = templateR(map(i)) 
         End do
         Do i = 1, nRow
            templateR(i) = templateC(i) 
         End do
      End if

      Return
      End



C ======================================================================
C Identity tensor
C ======================================================================
      Integer Function Ddiff(x, y, label, dDATA, iDATA, iSYS, Coef)
      include 'fem2Dtri.fd'

      Real*8  x, y, dDATA(*), Coef(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)

      iSYS(1) = 4
      iSYS(2) = 4

      Do i = 1, 4
         Do j = 1, 4
            Coef(i, j) = 0
         End do
         Coef(i, i) = 1
      End do

      Ddiff = TENSOR_SYMMETRIC + BC_NULL

      Return
      End



C ======================================================================
C Essential boundary conditions
C ======================================================================
      Integer Function Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
      Include 'fem2Dtri.fd'

      Real*8  x, y, dDATA(*), eBC(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)
      Logical left, right, top, bottom

      iSYS(1) = 2
      iSYS(2) = 1
   
      left = iDATA(3).EQ.1
      right = iDATA(4).EQ.1
      bottom = iDATA(2).EQ.1
      top = iDATA(1).EQ.1

      If(label.eq.3 .AND. left) Then
         eBC(1,1) = 1 - 4*y*y
         eBC(2,1) = 0
         Dbc = TENSOR_GENERAL + BC_DIRICHLET

      Else If(label.eq.3 .AND. .NOT.left) Then
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = TENSOR_GENERAL + BC_NEUMANN

      Else If(label.eq.1 .AND. right) Then
         eBC(1,1) = 1 - 4*y*y
         eBC(2,1) = 0
         Dbc = TENSOR_GENERAL + BC_DIRICHLET

      Else If(label.eq.1 .AND. .NOT.right) Then
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = TENSOR_GENERAL + BC_NEUMANN

      Else If(label.eq.2 .AND. bottom) Then
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = TENSOR_GENERAL + BC_DIRICHLET

      Else If(label.eq.2 .AND. .NOT.bottom) Then
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = TENSOR_GENERAL + BC_NEUMANN

      Else If(label.eq.4 .AND. top) Then
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = TENSOR_GENERAL + BC_DIRICHLET

      Else If(label.eq.4 .AND. .NOT.top) Then
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = TENSOR_GENERAL + BC_NEUMANN

      Else
         eBC(1,1) = 0
         eBC(2,1) = 0
         Dbc = TENSOR_GENERAL + BC_DIRICHLET
      End if

      Return
      End



C ======================================================================
C Right hand side
C ======================================================================
      Integer Function Drhs(x, y, label, dDATA, iDATA, iSYS, F)
      Include 'fem2Dtri.fd'

      Real*8  x, y, dDATA(*), F(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 1
     
      If(iDATA(5).EQ.0) Then
         F(1,1) = 0D0
      Else
         F(1,1) = x
      End if
      F(2,1) = 0D0

      Drhs = TENSOR_GENERAL + BC_NULL

      Return
      End 



C ======================================================================
C Parameterization of the unit circle centered at the origin
C ======================================================================
      Subroutine CrvFunction_user(tc, xyc, iFnc)
C ======================================================================
C  The routine computes the Cartesian coordinates of point
C  xyc from its parametric coordinate tc.
C
C  tc     - the given parametric coordinate of point
C  xyc(2) - the Cartesian coordinate of the same point
C  iFnc   - the function number for computing
C
C  On input :  tc, iFnc
C  On output:  xyc(2)
C
C *** Remarks:
C        1. This is the default routine to link library libani2D-2.x.a
C ======================================================================
      Real*8  tc, xyc(2)
      Real*8  PI, R

C ======================================================================
      PI = 4*datan(1D0)
      R  = 0.5D0

      If(iFnc.EQ.1) Then
         xyc(1) = -2.0 + 4.0 * tc
         xyc(2) =  1.0 + 0.5 * datan(8*tc-4) / datan(4D0)

      Else If(iFnc.EQ.2) Then
         xyc(1) = -2.0 + 4.0 * tc
         xyc(2) = -1.0 - 0.5 * datan(8*tc-4) / datan(4D0)

      Else If(iFnc.EQ.3) Then
         xyc(1) = R * dcos( PI*tc )
         xyc(2) = R * dsin( PI*tc )

      Else If(iFnc.EQ.4) Then
         xyc(1) = -R * dcos( PI*tc )
         xyc(2) = -R * dsin( PI*tc )

      Else
         Write(*, '(A,I5)') 'Undefined function =', iFnc
         Stop
      End if

      Return
      End



