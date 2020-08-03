c ======================================================================
      Program AnisotropicAdaptation2
c ======================================================================
c The adaptive iterative solution of the following boudary value problem:
c
c -div D grad u + du/dx = 1 in Omega
c                     u = 0 on dOmega
c
c where Omega = [0;1]x[0;1] - [0;1/3]x[0;0.5] (double rump example)
c           D = diag{0.001,0.001} 
c
c The program generates optimal mesh with nEStar elements
c
c The SUPG stabilization is not used here, since the adaptive mesh itself 
c stabilizes the approximation. We exploit here the fact that the 
c convection operator is simple and coinsides with the existing operator
c DUDX. In case of a general convection field refer to 
c TutorialAdvanced/ConvectionDiffusion/main.f
c ======================================================================
      implicit none

c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
c namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,namax
      parameter(nvmax = 150 000, ntmax = 2*nvmax, nbmax = 10 000)
      parameter(namax = 900 000)

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
      integer  Nbv, Nbl, bl(7,6)

      Data     Nbv/6/, Nbl/6/
      Data     bv/0,0.5, 0,1, 1,1, 1,0, 0.33333,0, 0.33333,0.5 / ! boundary nodes
      Data     bl/1,2,0,-1,1,1,0, 2,3,0,-1,1,1,0,                ! outer boundary edge
     &            3,4,0,-1,1,1,0, 4,5,0,-1,1,1,0,                ! outer boundary edge
     &            5,6,0,-1,1,1,0, 6,1,0,-1,1,1,0 /               ! outer boundary edge
      Data     bltail/0,0,0,0,0,0,0,0,0,0,0,0/                   ! curved  boundary edges

      Integer   aft2dboundary
      EXTERNAL  aft2dboundary


c ======================================================================
c for library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax)
      Real*8    A(namax), RHS(nvmax), SOL(nvmax), RES(nvmax)

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
      EXTERNAL ANI_CrvFunction


c LOCAL VARIABLEs
      Integer  i, iLoop, nLOOPs, nRow, nCol, iERR
      Real*8   rmax, h
      
c ======================================================================
c number of adaptive loops
      nLOOPs = 10


c === generate a quasi-uniform mesh with mesh step h
      h = 0.02                           
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary(
     &     Nbv, bv, Nbl, bl, bltail,  h,  ! input geometric data
     &     nv, vrt,                       ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)

      If(iERR.NE.0)   Stop 'Error in function aft2dboundary'
      If(nv.GT.nvmax) Stop 'too many nodes, increase nvmax'
      If(nt.GT.ntmax) Stop 'too many triangles, increase ntmax'
      If(nb.GT.nbmax) Stop 'too many boundary edges, increase nbmax'

      Write(*,'(A,3(I6,A))') '     mesh with', nt, ' triangles', 
     &                       nv, ' vertices ', nb, ' boundary edges' 
          

c === draw the initial mesh. The name must terminate with .ps
c demo has been activated.
c     Call graph(nv,vrt, nt,tri, 'mesh_initial.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_initial.ps',
     &               'Initial quasi-uniform mesh')


c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c === no data is provided for the user subroutine Ddiff
         dDATAFEM(1) = 0D0
         iDATAFEM(1) = 0

c mark the Dirichlet points with the maximal edge color
         Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &                Dbc, dDATAFEM, iDATAFEM, iSYS)

c general sparce matrix in the compressed column format
         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
         controlFEM(2) = 1

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2Dext, dDATAFEM, iDATAFEM, controlFEM,
     &        nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c === call the driver for LU factorization and solution
         Call CSC2CSC0(nCol, IA,JA,A)

c set up default control parameters & print only error messages
         Call umf4def(lucontrol)
         lucontrol(1) = 1

c pre-order and symbolic analysis
         Call umf4sym(nCol, nCol, IA,JA,A, symbolic, lucontrol,luinfo)
         If(luinfo(1).LT. 0) Goto 5001

c numeric factorization
         Call umf4num(IA,JA,A, symbolic,numeric, lucontrol,luinfo)
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


c === plot solution isolines. The name must terminate with .ps
c        demo graphics has been activated
         If(iLOOP.eq.1) then
c           Call isolines(SOL, nv,vrt, nt,tri, nb,bnd, 'iso.ps', 20)
            Call isolines_demo(SOL, nv,vrt, nt,tri, nb,bnd, 'iso.ps',
     &                    20, 'Solution isolines on the initial mesh') 
         Else
c           Call isolines(SOL, nv,vrt, nt,tri, nb,bnd, 'his.ps', 20)
            Call isolines_demo(SOL, nv,vrt, nt,tri, nb,bnd, 'his.ps',
     &                    20, 'Solution isolines on the adapted mesh') 
         End if

         If(iLoop.eq.nLOOPs) goto 500

c  ===   generate metric (from SOL) optimal for the L_p norm
         Lp = 0             ! maximum norm
c        Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(SOL,
     &                        vrt, nv, tri, nt, bnd, nb, Metric,
     &                        MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c === generate the adaptive mesh
         control(1) = 100     ! MaxSkipE
         control(2) = 50000   ! MaxQItr
         control(3) = 1       ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1       ! flagAuto
         control(5) = 1       ! iPrint = minimal level of output information
         control(6) = 0       ! iErrMesgt: only critical termination allowed

         nEStar = 10000
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
     &                        'unspecified error if mbaNodal')


c === draw the adapted mesh. The name must terminate with .ps
c demo graphics has been activated.
c        Call graph(nv,vrt, nt,tri, 'mesh_final.ps')
         Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps',
     &                  'Final mesh adapted to solution')

 500     Continue
      End do


c === testing the results
      If(Quality.LT.0.4) Stop 911
      Stop 


c === Error messages
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





