c ======================================================================
      Program ConvectionDiffusion
c ======================================================================
c The adaptive iterative solution of the following boudary value problem:
c
c -div D grad u + vx * du/dx + vy * du/dy = 1 in Omega
c                                       u = 0 on dOmega
c
c where Omega has a backstep-like shape,
c     D  = diag{0.01,0.01} 
c     vx = 3
c     vy = 2
c
c The SUPG stabilization is NOT used here, since the adaptive mesh 
c itself stabilizes the approximation.
c
c The program generates optimal mesh with nEStar elements.
c ======================================================================
      implicit none

c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
c namax - maximum number of non-zero matrix entries
      integer nvmax,ntmax,nbmax,namax
      parameter(nvmax = 150 000, ntmax = 2*nvmax, nbmax = 10 000)
      parameter(namax = 900 000)

c work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 1 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


c ======================================================================
c Mesh definition
c ======================================================================
      Integer  nv, nvfix, labelV(nvmax), fixedV(4)
      Real*8   vrt(2,nvmax)

      Integer  nb, nbfix, bnd(2,nbmax), labelB(nbmax), fixedB(1)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2,nbmax)

      Integer  nt, ntfix, tri(3,ntmax), labelT(ntmax), fixedT(1)

      DATA     nvfix/4/, fixedV/1,2,3,4/,  nbfix/0/, ntfix/0/


c ======================================================================
c for library aniAFT 
c ======================================================================
      Real*8   bv(2,6), bltail(2,6)
      Integer  Nbv, Nbl, bl(7,6)

      Data     Nbv/6/, Nbl/6/
      Data     bv/0.0,0.5, 0.0,1.0, 1.0,1.0,        ! boundary nodes
     &            1.0,0.0, 0.5,0.0, 0.5,0.5/
      Data     bl/1,2,0,-1,1,1,0, 2,3,0,-1,1,1,0,   ! outer boundary edge
     &            3,4,0,-1,1,1,0, 4,5,0,-1,1,1,0,
     &            5,6,0,-1,1,1,0, 6,1,0,-1,1,1,0 /  
      Data     bltail/0,0,0,0,0,0,0,0,0,0,0,0/      ! curved  boundary edges

      Integer   aft2dboundary
      EXTERNAL  aft2dboundary


c ======================================================================
c for library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax)
      Real*8    A(namax), RHS(nvmax), SOL(nvmax), RES(nvmax)

      Real*8   dDATAFEM(2)
      Integer  iDATAFEM(1)

      Integer  Dbc, nRow, nCol, iSYS(MAXiSYS), controlFEM(3)
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
      Real*8    Quality
      EXTERNAL  ANI_CrvFunction


c LOCAL VARIABLEs
      Integer  i, iLoop, iERR, nLOOPs
      Real*8   h, rmax

      
c ======================================================================
c number of adaptive loops
      nLOOPs = 9


c === generate a quasi-uniform mesh with mesh step h
      h = 0.026
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary(
     &     Nbv, bv, Nbl, bl, bltail,  h,  ! input geometric data
     &     nv, vrt,                       ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)
          
      If(iERR.NE.0)   Stop 'Error in function aft2dboundary'
      If(nv.GT.nvmax) Stop 'Too many nodes'
      If(nb.GT.nbmax) Stop 'Too many boundary edges'
      If(nt.GT.ntmax) Stop 'Too many triangles'

      Write(*,'(A,2I6)') '     number of triangles/vertices ',nt,nv 


c     plot the initial mesh. The name must terminate with .ps
c     Demo graphics has been activated
c     Call graph(nv,vrt, nt,tri, 'mesh_initial.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_initial.ps',
     &               'Initial quasi-uniform mesh')


c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c === mark end points of the Dirichlet edges 
         dDATAFEM(1) = 0
         Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &                Dbc, dDATAFEM, iDATAFEM, iSYS)

c data provided for the user subroutine Dconv employed in FEM2Dext
         dDATAFEM(1) = 3  ! v_x
         dDATAFEM(2) = 2  ! v_y
         iDATAFEM(1) = 0

c general sparce matrix in the complessed column format
         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
         controlFEM(2) = 1

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2Dext, dDATAFEM, iDATAFEM, controlFEM,
     &        nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c  ===   call the driver for LU factorization and solution
         Call CSC2CSC0(nCol, IA, JA)

c set up default control parameters & print only error messages
         Call umf4def(lucontrol)
         lucontrol(1) = 1

c        pre-order and symbolic analysis
         Call umf4sym(nCol, nCol, IA,JA,A, symbolic, lucontrol, luinfo)
         If(luinfo(1).LT.0) Goto 5001

c        numeric factorization
         Call umf4num(IA,JA,A, symbolic, numeric, lucontrol, luinfo)
         If(luinfo(1).lt.0) Goto 5002

c        free the symbolic analysis data
         Call umf4fsym(symbolic)

c solve Ax=b, without iterative refinement, x=SOL, b=RHS
         sys = 0
         Call umf4sol(sys, SOL, RHS, numeric, lucontrol, luinfo)
         If(luinfo(1).LT.0) Goto 5003

c        free the numeric factorization data
         Call umf4fnum (numeric)

c check the residual
         Call mulAcsc0(nRow, IA, JA, A, SOL, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(A,E12.6)') 'LU:  maximal norm of residual: ', rmax


c === plot 20 isolines. The name of PostScript fila must terminate with .ps
c        demo raphics has been activated
         If (iLoop.eq.1) then
c           Call isolines(SOL, nv,vrt, nt,tri, nb,bnd, 'iso_ini.ps',20)
            Call isolines_demo(SOL, nv,vrt, nt,tri, nb,bnd,'iso_ini.ps',
     &                    20, 'Solution isolines on the initial mesh')
            Call draw_matrix(nRow, IA, JA, 'matrix_ini.ps')
         Else
c           Call isolines(SOL, nv,vrt, nt,tri, nb,bnd, 'iso_fin.ps',20)
            Call isolines_demo(SOL, nv,vrt, nt,tri, nb,bnd,'iso_fin.ps',
     &                    20, 'Solution isolines on the final mesh')
            Call draw_matrix(nRow, IA, JA, 'matrix_fin.ps')
         End if

c  ===   generate metric (from SOL) optimal for the L_p norm
         Lp = 0             ! maximum norm
c        Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(SOL,
     &                        vrt, nv, tri, nt, bnd, nb, Metric,
     &                        MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c === build a mesh quasi-uniform in metric Metric
         control(1) = 100     ! MaxSkipE
         control(2) = 50000   ! MaxQItr
         control(3) = 1       ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1       ! flagAuto
         control(5) = 1       ! iPrint = minimal level of output information
         control(6) = 0       ! iErrMesgt: only critical termination allowed

         nEStar = 4000
         Quality = 0.7

         Call mbaNodal(
c group (M)
     &        nv, nvfix, nvmax, vrt, labelV, fixedV,
     &        nb, nbfix, nbmax, bnd, labelB, fixedB,
     &        nc,               crv, labelC, ANI_CrvFunction,
     &        nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &        nEStar, Quality, control, Metric,
c group (W)
     &        MaxWr, MaxWi, rW, iW, iERR)

         If(iERR.GT.1000) Call errMesMBA(iERR, 'main',
     &                        'unspecified error in mbaNodal')


c === plot the adapted mesh. The name must terminate with .ps
c        Call graph(nv,vrt, nt,tri, 'mesh_final.ps')
         Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps',
     &                  'Final mesh adapted to the solution')
      End do


c === testing the results
      If(Quality.LT.0.4) Stop 911
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





