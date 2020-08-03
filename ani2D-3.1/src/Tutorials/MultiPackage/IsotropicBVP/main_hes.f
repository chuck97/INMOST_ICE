c ======================================================================
      Program IsotropicAdaptation
c ======================================================================
c The adaptive solution for the following boudary value problem (bvp):
c
c   -div   grad u = 0   in  Omega
c               u = g   on  dOmega_D
c          du/dn  = 0   on  dOmega_N
c
c where Omega = { x^2 + y^2 < 1 \ {y=0, 0<x<1} }
c    dOmega_N = { y=0, 0<x<1 } from  bottom
c    dOmega_D = { x^2 + y^2 = 1} .or. {y=0, 0<x<1} 
c           g = r^0.25 sin phi/4
c
c Exact solution of this problem is u=r^0.25 sin phi/4  (in polar coordinates)
c
c The program generates optimal mesh with nEStar elements.
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
c for library aniAFT  (we define domain Omega)
c ======================================================================
      Real*8    bv(2,6), bltail(2,6)
      Integer   Nbv, Nbl, bl(7,6)

      DATA      Nbv/6/, Nbl/6/
      DATA      bv/0,-1, -1,0, 0,1, 1,0, 0,0, 1,0/    ! boundary nodes (last coord.pare is fictitious)
      DATA      bl/1,2,1,-1,1,1,0, 2,3,1,-1,2,1,0,    ! outer boundary edges
     &             3,4,1,-1,3,1,0, 4,5,0,-1,5,1,0,    ! outer boundary edges
     &             5,6,0,-1,6,1,0, 6,1,1,-1,4,1,0/    ! outer boundary edges
      DATA      bltail/1.5,1, 1,0.5, 0.5,0, 0,0,0,0,1.999999,1.5/   ! data for each curved boundary edge

      Integer   aft2dboundary
      EXTERNAL  aft2dboundary
      EXTERNAL  userboundary      


c ======================================================================
c for library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax)
      Real*8    A(namax), RHS(nvmax), SOL(nvmax), RES(nvmax)

      Integer  iDATA(1), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATA(1)

      Integer  Dbc
      EXTERNAL Dbc,  FEM2Dext


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
      Integer  i,n, iLoop, iv1,iv2,iv3, label, nRow, nCol, nLOOPs, iERR
      Real*8   trierr, errLp, errLpGrad, edge_length, calNorm
      Real*8   SOL_loc(3), xyt(2), rmax, hmax, h

      Integer  DexactU, DgradU, ANI_Dnull
      EXTERNAL DexactU, DgradU, ANI_Dnull, edge_length, calNorm

c ======================================================================
c number of adaptive loops
      nLOOPs = 10

c === call library AniAFT to generate a mesh
c     pass the user function userBoundary (forlibaft.c)
      Call registeruserfn(userboundary)  ! register the name to be used in the library


c     generate a quasi-uniform mesh with mesh step h
      h = 0.04  
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary( 
     &       Nbv, bv, Nbl, bl, bltail, h,  ! input geometric data
     &       nv, vrt,                      ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)

      If(iERR.NE.0) Stop 'Error in function aft2dboundary'
      If(nv.GT.nvmax) stop ' too many nodes, increase nvmax'
      If(nt.GT.ntmax) stop ' too many triangles, increase ntmax'
      If(nb.GT.nbmax) stop ' too many boundary edges, increase nbmax'
          
      Write(*,'(A,3(I6,A))') '     mesh with', nt, ' triangles', 
     &                       nv, ' vertices ', nb, ' boundary edges' 

          
c === draw Postscript file of the initial mesh, name must have extension .ps
c     Demo graphics as been activated
c     Call graph(nv,vrt, nt,tri, 'mesh_initial.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_initial.ps',
     &               'Initial quasi-uniform mesh')


c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c === no data is provided for the user subroutine FEM2Dext
         dDATA(1) = 0D0
         iDATA(1) = 0

c mark the Dirichlet points with the maximal edge color
         Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &                Dbc, dDATA, iDATA, iSYS)

c general sparce matrix in the compressed row format
         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
         controlFEM(2) = 1

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2Dext, dDATA, iDATA, controlFEM,
     &        nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c === call library aniLU
c converting to the 0-based format
         Call CSC2CSC0(nCol, IA,JA,A)

c set up default control parameters & print only error messages
         Call umf4def(lucontrol)
         lucontrol(1) = 1

c pre-order and symbolic analysis
         Call umf4sym(nCol, nCol, IA,JA,A, symbolic, lucontrol,luinfo)
         If(luinfo(1).LT.0) Goto 5001

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


c === make PostScript file with solution isolines. The name must have extension .ps
c demo graphics has been activated
         If(iLoop.EQ.1) then
c           Call isolines(SOL, nv,vrt, nt,tri, nb,bnd, 'iso.ps',20)
            Call isolines_demo(SOL, nv,vrt, nt,tri, nb,bnd, 'iso.ps',
     &                    20, 'Solution isolines on the initial mesh')
         Else
c           Call isolines(SOL, nv,vrt, nt,tri, nb,bnd, 'his.ps',20)
            Call isolines_demo(SOL, nv,vrt, nt,tri, nb,bnd, 'his.ps',
     &                    20, 'Solution isolines on the final mesh')
         End if


c === call library aniFEM to compute Lp-error
         Lp = 0
         errLp = 0D0
         errLpGrad = 0D0

         Do n = 1, nt
            label = labelT(n)

            Do i = 1, 3
               iv3 = tri(i, n)
               SOL_loc(i) = SOL(iv3)
            End do

            iv1 = tri(1, n)
            iv2 = tri(2, n)
            iv3 = tri(3, n)

            Call fem2Derr(vrt(1,iv1), vrt(1,iv2), vrt(1,iv3), Lp,
     &                    IDEN, FEM_P1, SOL_loc, DexactU, dDATA, iDATA,
     &                    label, ANI_Dnull, dDATA,iDATA,iSYS, 5, trierr)

            If (Lp.eq.0) then
               errLp = max(errLp, trierr)
            Else
               errLp = errLp + trierr
            End if

            Call fem2Derr(vrt(1,iv1), vrt(1,iv2), vrt(1,iv3), Lp,
     &                    GRAD, FEM_P1, SOL_loc, DgradU, dDATA, iDATA,
     &                    label, ANI_Dnull, dDATA,iDATA,iSYS, 5, trierr)
 
            If(Lp.EQ.0) Then
               errLpGrad = max(errLpGrad, trierr)
            Else
               errLpGrad = errLpGrad + trierr
            End if
         End do

         If(Lp.GT.0) Then
            errLp     = errLp     ** (1D0 / Lp)
            errLpGrad = errLpGrad ** (1D0 / Lp)
         End if

         Write(*,'(2(A,E11.5))') 'FEM: Interpolation error = ', errLp, 
     &                           ',  gradient error = ', errLpGrad

         If (iLoop.eq.nLOOPs) goto 500


c === call library aniLMR to generate a metric from solution SOL
         Lp = 0             ! maximum norm
c        Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(SOL,
     &                        vrt, nv, tri, nt, bnd, nb, Metric,
     &                        MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c === call library aniMBA to build a new mesh
         control(1) = 100     ! MaxSkipE
         control(2) = 10000   ! MaxQItr
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
     &        nc,               crv, labelC, CrvFunction_user,
     &        nt, ntfix, ntmax, tri, labelT, fixedT,
c group (CONTROL)
     &        nEStar, Quality, control, Metric,
c group (W)
     &        MaxWr, MaxWi, rW, iW, iERR)

         If(iERR.GT.1000) Call errMesMBA(iERR, 'main',
     &                        'unspecified error in mbaNodal')


c === draw Postscript file of the adapted mesh
c        The name must have extension with .ps
c        Call graph(nv,vrt, nt,tri, 'mesh_final.ps')
         Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps',
     &                  'Final mesh adapted to the solution')


c === save distribution of triangle sizes w.r.t. distance to the origin
         Open(10, file = 'gnu')
         Do n = 1, nt
            iv1 = tri(1, n) 
            iv2 = tri(2, n) 
            iv3 = tri(3, n) 

            hmax = max(edge_length(vrt(1,iv1), vrt(1,iv2)),
     &                 edge_length(vrt(1,iv1), vrt(1,iv3)),
     &                 edge_length(vrt(1,iv2), vrt(1,iv3)))

            Do i = 1, 2
               xyt(i) = (vrt(i,iv1) + vrt(i,iv2) + vrt(i,iv3)) / 3
            End do
            rmax = calNorm(xyt)

            Write(10, *) hmax, rmax
         End do
         Close(10)

 500     Continue
      End do


c === testing the results
      If(Quality.LT.0.6) Stop 911
      Stop 


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

