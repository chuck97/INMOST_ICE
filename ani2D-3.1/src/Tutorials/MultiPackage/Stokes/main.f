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
      parameter(nvmax = 50 000, ntmax = 2*nvmax, nbmax = 5 000)
      parameter(nfmax = 200 000, namax = 2 000 000)

c work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 4 000 000, MaxWi = 6 000 000)

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

      Integer   aft2dboundary
      EXTERNAL  aft2dboundary, userboundary


c ======================================================================
c for library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nfmax), JA(namax)
      Real*8    A(namax), RHS(nfmax), SOL(nfmax), RES(nfmax)

      Integer  iDATAFEM(1), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATAFEM(1)

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
      Integer  nr, iux,iuy, ip, ibc, iv1,iv2
      Real*8   rmax, h, x, y, eBC(2)
      Character*30 file_mesh, file_velocity, file_pressure

c ======================================================================
c number of adaptive loops
      nLOOPs = 5

c ... pass the user function userBoundary (forlibaft.c) 
      Call registeruserfn(userboundary) ! register the name in the


c === generate a quasi-uniform mesh with mesh step h
      h = 0.096
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary(
     &       Nbv, bv, Nbl, bl, bltail, h,    ! input geometric data
     &       nv, vrt,                        ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)

      If(iERR.NE.0)   Stop 'Error in function aft2dboundary'
      If(nv.GT.nvmax) Stop 'too many nodes, increase nvmax'
      If(nt.GT.ntmax) Stop 'too many triangles, increase ntmax'
      If(nb.GT.nbmax) Stop 'too many boundary edges, increase nbmax'

      Write(*,'(A,3(I6,A))') '     mesh with', nt, ' triangles', 
     &                       nv, ' vertices ', nb, ' boundary edges' 


c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c === no data is provided for the user subroutine Ddiff
         dDATAFEM(1) = 0D0
         iDATAFEM(1) = 0

c mark the Dirichlet points with the maximal edge color
         Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &                Dbc, dDATAFEM, iDATAFEM, iSYS)


c === general sparse matrix in a 0-based CSC format used in UMFPACK
         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
         controlFEM(2) = 1

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2Dext, dDATAFEM, iDATAFEM, controlFEM,
     &        nfmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c  ===   call the driver for LU factorization and solution
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
c                       nodal y-component of velocity starts at SOL(iuy),
c                       nodal pressure starts at SOL(ip)
c                       mid-edge x-component of velocity starts at SOL(iuxe)
c                       mid-edge y-component of velocity starts at SOL(iuye)
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
            Call draw_matrix(nRow, IA, JA, "matrix_ini.ps")
         Else
            Call isolines_demo(rW, nv,vrt, nt,tri, nb,bnd,file_velocity,
     &           20, 'Velocity streamlines on the final mesh')
            Call isolines_demo(SOL(ip),nv,vrt, nt,tri, nb,bnd,
     &           file_pressure,
     &           40, 'Pressure isolines on the final mesh')
            Call graph_demo(nv,vrt, nt,tri, file_mesh, 
     &           'Final mesh adapted to x-component of velocity')
            Call draw_matrix(nRow, IA, JA, "matrix_fin.ps")
         End if

         If(iLoop.eq.nLOOPs) Goto 500


c  ===   generate metric (from SOL) optimal for the L_p norm
c        Lp = 0             ! maximum norm
         Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(SOL(iux),           ! mesh is adapted to u_x
c        Call Nodal2MetricVAR(SOL(iuy),           ! mesh is adapted to u_y
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

         nEStar = 7000
         Quality = 0.7

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


c === testing the results
      If(Quality.LT.0.5) Stop 911
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

