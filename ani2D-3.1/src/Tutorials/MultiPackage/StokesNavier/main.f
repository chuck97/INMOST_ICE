c ======================================================================
      Program NavierStokes
c ======================================================================
c The adaptive solution for the lid driven cavity boundary value problem:
c
c -div nu grad u  + (u.grad) u + grad p = 0   in Omega
c          div u                        = 0   in Omega
c
c                                     u = u_0 on dOmega_1
c                                     u = 0   on dOmega_2
c
c where Omega is the unit square, dOmega_1 is the side at y=1, 
c and dOmega_2 is the rest of the boundary. The non-homogeneous
c boundary condition is 
c
c    u_0 = { 1, 0 }.
c
c We use the stable pair of finite elements: P2 for velocity u and P1 
c for the pressure p. The convection term is NOT stabilized by SUPG. 
c The discretization method results in a system of nonlinear algebraic 
c equations with the saddle point structure. This system is solved by 
c the Newton-Krylov method. The initial guess and the precoditioner 
c (LU-factorized Oseen problem) is formed after  Picard  iteration.
c Once the mesh is adapted, the solution is L2-projected onto the new mesh
c using the remapping package aniPRJ.
c
c Stabilization of convective terms is provided by mesh adaptation.
c Since the viscosity is small (nu = 3e-4) and the initial mesh is not adapted,
c the nonlinear vector-function is not smooth and Newton-Krylov does not converge.
c Therefore we solve the first nonlinear problem by limited number (15) of Picard iterations,
c followed by Newton-Krylov iterations. The nonlinear problems on adapted meshes
c have already good initial guess due to the L2-projection of the solution from the previous mesh
c and the adapted mesh stabilizes the solution. Therefore we apply only one Picard iteration
c in order to compute the preconditioner and the initial guess which
c provides fast convergence of Newton-Krylov.
c ======================================================================
      implicit none

c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
c nfmax - maximum order of system + 1
c namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,nfmax,namax
      parameter(nvmax =  50 000, ntmax = 2*nvmax, nbmax = 5 000)
      parameter(nfmax = 250 000, namax = 5 000 000)

c work memory
      Integer   MaxWr, MaxWi
c     Parameter(MaxWr = 3 000 000, MaxWi = 5 000 000)
      Parameter(MaxWr = 6 000 000, MaxWi = 10 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


c ======================================================================
c Mesh definition: vertices (v), boundary edges (b), and triangles (t)
c ======================================================================
      Integer  nv, nvfix, labelV(nvmax), fixedV(4)
      Real*8   vrt(2,nvmax)

      Integer  nb, nbfix, bnd(2,nbmax), labelB(nbmax), fixedB(1)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2,nbmax)

      Integer  nt, ntfix, tri(3,ntmax), labelT(ntmax), fixedT(1)

      DATA     nvfix/4/, fixedV/1,2,3,4/,  nbfix/0/, ntfix/0/


c ======================================================================
c for library aniPRJ: second mesh and meta-mesh
c ======================================================================
      Integer   nvMetaMax, ntMetaMax
      parameter(nvMetaMax = 5*nvmax, ntMetaMax = 5*ntmax)

      Integer  nv0, labelV0(nvmax)
      Integer  nb0, bnd0(2,nbmax), labelB0(nbmax)
      Integer  nt0, tri0(3,ntmax), labelT0(ntmax)
      Real*8   vrt0(2,nvmax)

      Integer  nvMeta,ntMeta, triMeta(3,ntMetaMax), parents(2,ntMetaMax)
      Real*8   vrtMeta(2,nvMetaMax)


c ======================================================================
c for library aniAFT  (we define domain Omega)
c ======================================================================
      Real*8   bv(2,4),bltail(2,4)
      Integer  Nbv,Nbl,bl(7,4)

      Data     Nbv/4/,Nbl/4/
      Data     bv/0.0,1.0, 1.0,1.0, 1.0,0.0, 0.0,0.0/  ! boundary vertices
      Data     bl/1,2,0,-1,1,1,0,     2,3,0,-1,2,1,0,  ! boundary edges
     &            3,4,0,-1,3,1,0,     4,1,0,-1,4,1,0/
      Data     bltail/8 * 0.0/  ! all edges are straight

      Integer  aft2dboundary
      EXTERNAL aft2dboundary

c user defined meshSize function for aniAFT (override quasi-uniform mesh)
      Real*8   meshSize
      EXTERNAL meshSize 


c ======================================================================
c for library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nfmax), JA(namax)
      Real*8    A(namax), RHS(nfmax), SOL(nfmax), RES(nfmax)
      Real*8   SOL0(nfmax)

      Integer  iDATAFEM(1), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATAFEM(1)

      Integer  Dbc
      EXTERNAL Dbc, FEM2Dext, FEM2DPrj1ext, FEM2DPrj2ext


c ======================================================================
c for library aniLU
c ======================================================================
      Integer  symbolic(2), numeric(2), sys
      Real*8   lucontrol(20), luinfo(90)


c ======================================================================
c for library aniINB
c ======================================================================
      Integer  IPREVEC, infoINB(5)
      Real*8   RESID, STPTOL
      EXTERNAL fnlin, prevec


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
      EXTERNAL ANI_CrvFunction

c ======================================================================
c for library blas
c ======================================================================
      Real*8   dnrm2
      EXTERNAL dnrm2, daxpy 

c LOCAL VARIABLEs
      Integer  i,j,k,n, iLoop, nLOOPs, iPicard, nPICARDs
      Integer  nRow, nCol, iERR, iux,iuy,ip, iuxe,iuye, ne0
      Integer  LenWork, nnz, nr, ne, ipRMesh, ipIMesh
      Integer  iIRE, inEP, iIEP, iEnd
      Real*8   rmax, h, x,y, eBC(2), TOLERANCE, Viscosity 
      Character*30 file_mesh, file_velx, file_vely, file_pres, file_strm

c ======================================================================
      Call registersizefn(meshSize) ! register function in library aniAFT

c number of adaptive loops
      nLOOPs = 5

c tolerance for eucledian norm of nonlinear residual
      TOLERANCE  = 1d-6

c initial guess for the INB is provided by Picard iterations, maximal number of iterations 
      nPICARDs = 15 

c === generate a shape regular mesh with mesh size controlled by meshSize.
c parameter h in this case is dummy since meshSize will take over.

      iERR = aft2dboundary(
     &       Nbv, bv, Nbl, bl, bltail, h,  ! input geometric data
     &       nv, vrt,                      ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)

      If(iERR.NE.0) Then
         Write(*,*) 'Error in function aft2dboundary'
         Stop 911
      End if

      Write(*,'(A,3I7)') 'Mesh nodes/bnd edges/triangles:', nv, nb, nt
      Write(*,'(A,3I7)') 'Max  nodes/bnd edges/triangles:', 
     &                   nvmax, nbmax, ntmax
      If(nv.GT.nvmax .OR. nb.GT.nbmax .OR. nt.GT.ntmax) Then
         Write(*,*) 'Too many nodes/boundary edges/triangles'
         Stop 911
      End if


c === Initial guess for the first adaptive step. 
c It is passed to the user subroutine Dconv 
      Do i = 1, nfmax
         SOL(i) = 0D0
      End Do

c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c mark the Dirichlet points with maximal edge color
c cornerpoints (0,1) and (1,1) get color >1 (homogeneous Dirichlet)
         Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &                Dbc, dDATAFEM, iDATAFEM, iSYS)


c if the initial guess for Newton-Krylov is good, it is enough to apply one Picard iteration
         If (iLoop.gt.1) nPICARDs = 1

c ======================================================================
c Begin PICARD iterations
         Write(*,'(A)')   '===> Picards iterations'
         Do iPicard = 1, nPICARDs
c === generate matrix using library aniFEM
            controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
            controlFEM(2) = 0

            Call BilinearFormTemplate(
     &           nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &           FEM2Dext, SOL, iDATAFEM, controlFEM,
     &           nfmax, namax, IA, JA, A, RHS, nRow, nCol,
     &           MaxWi, MaxWr, iW, rW)

c === check the residual, if converged then exit from the loop
            Call mulAcsc(nRow, IA, JA, A, SOL, RES)
            Call daxpy(nRow, -1D0, RHS, 1, RES, 1)
            RESID = dnrm2(nRow, RES, 1)

            Write(*,'(A,I2,A,E12.5)') 
     &         'itr=', iPicard, ',  ||F(SOL)|| = ', RESID

            If(RESID.LT.TOLERANCE) GOTO 777

c === solve the problem using library aniLU
            Call CSC2CSC0(nCol, IA,JA,A)  ! 0-based format

            If(iPicard.EQ.1) Then
               Call umf4def(lucontrol)  ! default control parameters
               lucontrol(1) = 1

c pre-order and symbolic analysis
               Call umf4sym(nCol,nCol,IA,JA,A,symbolic,lucontrol,luinfo)
               If(luinfo(1).LT.0) Goto 5001
            End If

c numeric factorization
            Call umf4num(IA,JA,A, symbolic,numeric, lucontrol,luinfo)
            If(luinfo(1).lt.0) Goto 5002

c solve Ax=b, without iterative refinement
            sys = 0
            Call umf4sol(sys, SOL,RHS, numeric, lucontrol,luinfo)
            If(luinfo(1).LT.0) Goto 5003

c free the numerical data
            If(iPicard.LT.nPICARDs) Call umf4fnum(numeric)
         End do
c free the symbolic analysis data
         Call umf4fsym(symbolic)

c END of the PICARD iterations
c SOL is the approximation to the solution of the Oseen problem
 777     Continue
c ======================================================================

c === Solve non-linear problem using library AniINB
c pack mesh for the internal function fnlin 
         nnz = IA(nRow + 1) - 1
         iW(6) = nnz + 2  !reserve some space
         iW(7) = MaxWi - 10 - 4*nt - nv - 3*nb - (nRow+1) - (nnz+2)

         If(iW(7).LE.0) Then
            Write(*,*) 'Size of the integer array iW is small'
            Stop 911
         End if
         LenWork = 11*nRow  ! work memory for INB

         If(LenWork + 2*nv + nRow + nnz .GT. MaxWr) Then
            Write(*,*) 'Size of the Real*8 array rW is small'
            Stop 911
         End if
        
         ipIMesh = 1
         ipRMesh = LenWork + 1

         Call packMesh(nv,vrt,labelV, nb,bnd,labelB, nt,tri,labelT, 
     &                 iW(ipIMesh), rW(ipRMesh))

c backup the original mesh to mesh0
         nv0 = nv
         Do n = 1, nv
            vrt0(1, n) = vrt(1, n)
            vrt0(2, n) = vrt(2, n)
            labelV0(n) = labelV(n)
         End do

         nb0 = nb
         Do n = 1, nb
            bnd0(1, n) = bnd(1, n)
            bnd0(2, n) = bnd(2, n)
            labelB0(n) = labelB(n)
         End do
         
         nt0 = nt
         Do n = 1, nt
            tri0(1, n) = tri(1, n)
            tri0(2, n) = tri(2, n)
            tri0(3, n) = tri(3, n)
            labelT0(n) = labelT(n)
         End do

         ne = (nRow - 3*nv) / 2  ! number of edges
         ne0 = ne


c call driver for the inexact Newton-Krylov backtracking method
         Write(*,'(A)') '===> Start Library LIBINB:'

         RESID = TOLERANCE 
         STPTOL = 1d-8
         IPREVEC = nRow
         infoINB(1) = 0   ! initializing successful termination flag for Newton
         infoINB(2) =  5  ! maximal number of linear iterations
         infoINB(3) = 10  ! maximal number of nonlinear iterations
         infoINB(4) = 10  ! maximal number of backtracks
         infoINB(5) = 1   ! print level (0 mum, 1 nonlinear residual, 2 linear residual)
        
         Call slInexactNewton(prevec, IPREVEC, numeric, lucontrol,
     &                        fnlin, rW(ipRMesh), iW(ipIMesh), 
     &                        nRow, SOL,
     &                        RESID, STPTOL,
     &                        rW, LenWork,
     &                        infoINB)
        
         If(infoINB(1).ne.0) then
            Write(*,*) 'failed to solve the problem with',
     &                 ' desired accuracy, info =', infoINB(1)
         Else
            Write(*,*)'||F(SOL)|| = ', RESID
            Write(*,*)'Number of linear iterations:   ', infoINB(2)
            Write(*,*)'Number of nonlinear iterations:', infoINB(3)
            Write(*,*)'Number of backtracks:          ', infoINB(4)
            Write(*,*)'Number of function evaluations:', infoINB(5)
         End if

c free the numeric factorization data
         Call umf4fnum (numeric)

c backup the solution SOL to SOL0
         Call dcopy(nRow, SOL, 1, SOL0, 1)

c partition of SOL: nodal x-component of velocity starts at SOL(iux),
c                   nodal y-component of velocity starts at SOL(iuy),
c                   nodal pressure starts at SOL(ip)
c                   mid-edge x-component of velocity starts at SOL(iuxe)
c                   mid-edge y-component of velocity starts at SOL(iuye)
c        ne  = (nRow - 3*nv) / 2
c        iux = 1
c        iuy = iux + nv
c        ip  = iuy + nv
c        iuxe= ip  + nv
c        iuye= iuxe+ ne


c === draw mesh and solution isolines in PostScript files using library aniVIEW.
         iux = 1
         iuy = iux + nv
         ip  = iuy + nv
         Call stream_function(nv, nr, nb, nt, vrt, bnd, tri,
     &                        SOL(iux), rW, MaxWi, iW)
c rW(1) contains streamlines of velocities
c The names must terminate with .ps. 
         If(iLoop.EQ.1) Then
            file_mesh = 'mesh_initial.ps'
            file_velx = 'velx_ini.ps'
            file_vely = 'vely_ini.ps'
            file_pres = 'pres_ini.ps'
            file_strm = 'streamlines_ini.ps'
         Else
            file_mesh = 'mesh_final.ps'
            file_velx = 'velx_fin.ps'
            file_vely = 'vely_fin.ps'
            file_pres = 'pres_fin.ps'
            file_strm = 'streamlines_fin.ps'
         End if

c        Call isolines(SOL(iux), nv,vrt, nt,tri, nb,bnd, file_velx, 20)
c        Call isolines(SOL(iuy), nv,vrt, nt,tri, nb,bnd, file_vely, 20)
c        Call isolines(SOL(ip),  nv,vrt, nt,tri, nb,bnd, file_pres, 20)
c        Call isolines(rW,       nv,vrt, nt,tri, nb,bnd, file_strm, 20)
c        Call graph(             nv,vrt, nt,tri,         file_mesh)
         If(iLoop.EQ.1) Then
            Call isolines_demo(SOL(iux), nv,vrt,nt,tri,nb,bnd,file_velx,
     &           20, 'Velocity isolines (u_x) on the initial mesh')
            Call isolines_demo(SOL(iuy), nv,vrt,nt,tri,nb,bnd,file_vely,
     &           20, 'Velocity isolines (u_y) on the initial mesh')
c           Call isolines_demo(SOL(ip),  nv,vrt,nt,tri,nb,bnd,file_pres,
c    &           20, 'Pressure isolines on the initial mesh')
            Call isolines_demo(rW, nv,vrt, nt,tri, nb,bnd,file_strm,
     &           50, 'Velocity streamlines on the initial mesh')
            Call graph_demo(nv,vrt, nt,tri, file_mesh, 
     &           'Initial shape regular mesh')
            Call draw_matrix(nRow, IA, JA, 'matrix_ini.ps')
         Else
            Call isolines_demo(SOL(iux), nv,vrt,nt,tri,nb,bnd,file_velx,
     &           20, 'Velocity isolines (u_x) on the final mesh')
            Call isolines_demo(SOL(iuy), nv,vrt,nt,tri,nb,bnd,file_vely,
     &           20, 'Velocity isolines (u_y) on the final mesh')
c           Call isolines_demo(SOL(ip),  nv,vrt,nt,tri,nb,bnd,file_pres,
c    &           20, 'Pressure isolines on the final mesh')
            Call isolines_demo(rW, nv,vrt, nt,tri, nb,bnd,file_strm,
     &           50, 'Velocity streamlines on the final mesh')
            Call graph_demo(nv,vrt, nt,tri, file_mesh, 
     &           'Final mesh adapted to velocity module')
            Call draw_matrix(nRow, IA, JA, 'matrix_fin.ps')
         End if

         If(iLoop.eq.nLOOPs) goto 500


c === Denerate metric from SOL using library aniLMR.
c This metric is optimal for the L_p norm of |u|.
c        Lp = 0  ! maximum norm
         Lp = 2  ! L_2 norm

         Do i = 0, nv-1
            RES(i + 1) = dsqrt(SOL(iux + i)**2 + SOL(iuy + i)**2)
         End Do

c        Call Nodal2MetricVAR(SOL(iux),  ! mesh is adapted to u_x
c        Call Nodal2MetricVAR(SOL(iuy),  ! mesh is adapted to u_y
         Call Nodal2MetricVAR(RES,       ! mesh is adapted to velocity module
     &                        vrt, nv, tri, nt, bnd, nb, Metric,
     &                        MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c === generate mesh adapted to solution using library AniMBA.
         control(1) = 100    ! MaxSkipE
         control(2) = 10000  ! MaxQItr
         control(3) = 1      ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1      ! flagAuto
         control(5) = 1      ! iPrint = minimal level of output information
         control(6) = 0      ! iErrMesgt: only critical termination allowed

         nEStar = 10000
         Quality = 0.7

         Call mbaNodal(
     &        nv, nvfix, nvmax, vrt, labelV, fixedV,  ! Mesh
     &        nb, nbfix, nbmax, bnd, labelB, fixedB,
     &        nc,               crv, labelC, ANI_CrvFunction,
     &        nt, ntfix, ntmax, tri, labelT, fixedT,
     &        nEStar, Quality, control, Metric,  ! Control parameters
     &        MaxWr, MaxWi, rW, iW, iERR)  ! Work arrays

         If(iERR.GT.1000) Call errMesMBA(iERR, 'main',
     &                        'unspecified error if mbaNodal')


c ======================================================================
c === project solution onto the new mesh using four libraries aniMBA, 
c aniPRJ, aniFEM, and aniLU. We do it in three steps, two for the 
c velocity components and one for the pressure
c === find number of edges using library aniMBA
         inEP = 1
         iIEP = inEP + nv
         iIRE = iIEP + 3 * nt
         iEnd = iIRE + 3 * nt
         If(iEnd.GT.MaxWi) Goto 500

         Call listE2R(nv, ne, nt, tri, iW(iIRE), iW(inEP), iW(iIEP))

c === create a meta-mesh using library aniPRJ
         Call dcopy(nRow, SOL, 1, SOL0, 1)  ! Backup solution

         Call MetaMesh(nv, vrt, nt, tri, nv0, vrt0, nt0, tri0,
     &                 nvMeta, nvMetaMax, vrtMeta, 
     &                 ntMeta, ntMetaMax, triMeta, parents,
     &                 MaxWi, MaxWr, iW, rW, iERR)
 

c === generate a mass matrix on meta-mesh with the user subroutine FEM2DPrj2ext()
         dDATAFEM(1) = 0D0
         iDATAFEM(1) = 0

         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
         controlFEM(2) = 0
         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2DPrj2ext, dDATAFEM, iDATAFEM, controlFEM,
     &        nfmax, namax, IA, JA, A, RES, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)

c === solve the projection problem on using library aniLU
         Call CSC2CSC0(nCol, IA,JA,A)  ! converting to 0-based format

         Call umf4def(lucontrol)  ! default control parameters
         lucontrol(1) = 1

c pre-order and symbolic analysis
         Call umf4sym(nCol,nCol,IA,JA,A,symbolic,lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5001

c numeric factorization
         Call umf4num(IA,JA,A, symbolic,numeric, lucontrol,luinfo)
         If(luinfo(1).lt.0) Goto 5002

c project x-component of velocity
         Call dcopy(nv0, SOL0, 1, RES, 1)
         Call dcopy(ne0, SOL0(3*nv0+1), 1, RES(nv0+1), 1)

         Call assemble_rhs(nv, vrt, nt, tri, nv0, vrt0, nt0, tri0,
     &                     nvMeta, vrtMeta, ntMeta, triMeta, parents,
     &                     IDEN, FEM_P2, IDEN, FEM_P2,
     &                     RHS, RES, MaxWi, iW)

c solve Ax=b, without iterative refinement
         sys = 0
         Call umf4sol(sys, RES, RHS, numeric, lucontrol, luinfo)
         If(luinfo(1).LT.0) Goto 5003

c copy the projected x-component to the primary vector
         Call dcopy(nv, RES, 1, SOL, 1)
         Call dcopy(ne, RES(nv+1), 1, SOL(3*nv+1), 1)

c project y-component of velocity
         Call dcopy(nv0, SOL0(nv0+1), 1, RES, 1)
         Call dcopy(ne0, SOL0(3*nv0+ne0+1), 1, RES(nv0+1), 1)

         Call assemble_rhs(nv, vrt, nt, tri, nv0, vrt0, nt0, tri0,
     &                     nvMeta, vrtMeta, ntMeta, triMeta, parents,
     &                     IDEN, FEM_P2, IDEN, FEM_P2,
     &                     RHS, RES, MaxWi, iW)

c solve Ax=b, without iterative refinement
         sys = 0
         Call umf4sol(sys, RES,RHS, numeric, lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5003

c copy the projected y-component to the primary vector
         Call dcopy(nv, RES, 1, SOL(nv+1), 1)
         Call dcopy(ne, RES(nv+1), 1, SOL(3*nv+ne+1), 1)

         Call umf4fnum(numeric)  ! Free numerical data
         Call umf4fsym(symbolic)  ! Free symbolic analysis data

c === generate a mass matrix on meta-mesh with the user subroutine FEM2DPrj1ext().
         dDATAFEM(1) = 0D0
         iDATAFEM(1) = 0

         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
         controlFEM(2) = 0
         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2DPrj1ext, dDATAFEM, iDATAFEM, controlFEM,
     &        nfmax, namax, IA, JA, A, RES, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)

c === solve the projection problem on using library aniLU
         Call CSC2CSC0(nCol, IA,JA,A)  ! Convert to 0-based format
         Call umf4def(lucontrol)  ! Default parameters
         lucontrol(1) = 1

c pre-order and symbolic analysis
         Call umf4sym(nCol,nCol,IA,JA,A,symbolic,lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5001

c numeric factorization
         Call umf4num(IA,JA,A, symbolic,numeric, lucontrol,luinfo)
         If(luinfo(1).lt.0) Goto 5002

c project pressure
         Call dcopy(nv0, SOL0(2*nv0+1), 1, RES, 1)

         Call assemble_rhs(nv, vrt, nt, tri, nv0, vrt0, nt0, tri0,
     &                     nvMeta, vrtMeta, ntMeta, triMeta, parents,
     &                     IDEN, FEM_P1, IDEN, FEM_P1,
     &                     RHS, RES, MaxWi, iW)

c solve Ax=b, without iterative refinement
         sys = 0
         Call umf4sol(sys, RES, RHS, numeric, lucontrol, luinfo)
         If(luinfo(1).LT.0) Goto 5003

c copy the projected pressure to the primary vector
         Call dcopy(nv, RES, 1, SOL(2*nv+1), 1)

         Call umf4fnum(numeric)  ! Free numerical data
         Call umf4fsym(symbolic)  ! Free symbolic data

c End PROJECTION of solution from mesh 0 to the updted mesh.
c ======================================================================
c        If(iLoop.GT.1) Then
c           Call isolines_demo(
c    &           SOL(1), nv,vrt, nt,tri, nb,bnd, 'velx_tst.ps',
c    &           20, 'Test of projection velocity u_x')
c        End if

 500     Continue
      End do


c === verifying the results for CTest
      If(nLOOPs.gt.1 .AND. Quality.LT.0.3) Stop 911
      Stop 


c error messages
 5001 Continue
      Write(*,*) 'Error occurred in umf4sym: ', luinfo(1)
      Stop 911

 5002 Continue
      Write(*,*) 'Error occurred in umf4num: ', luinfo(1)
      Stop 911

 5003 Continue
      Write(*,*) 'Error occurred in umf4sol: ', luinfo(1)
      Stop 911

      End



