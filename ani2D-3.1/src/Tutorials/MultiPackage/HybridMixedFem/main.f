c ======================================================================
      Program HybridMixedFEM
c ======================================================================
c The adaptive solution for the boundary value problem:
c
c
c    div u = 1       in  Omega    (mass conservation eqn)
c        u =-grad p  in  Omega    (constitutive equation)
c        p = 0       on  dOmega   (boundary condition)
c
c
c Here Omega is the convex pentagon. 
c We use the mixed-hybrid finite element method with the lowest order 
c Raviart-Thomas elements. The method results in a problem with a 
c symmetric positive definite stifness matrix for the Lagrange 
c multipliers.
c ======================================================================
      implicit none

c  nvmax - maximum number of mesh nodes
c  ntmax - maximum number of mesh triangles
c  nbmax - maximum number of boundary edges
c  namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,namax
      parameter(nvmax = 100 000, ntmax = 2*nvmax, nbmax = 10 000)
      parameter(namax = 900 000)

c  work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 3 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


c ======================================================================
c Mesh definition: vertices (v), boundary edges (b), triangles (t)
c ======================================================================
      Integer  nv, nvfix, labelV(nvmax), fixedV(5) 
      Real*8   vrt(2,nvmax)

      Integer  nb, nbfix, bnd(2,nbmax), labelB(nbmax), fixedB(1)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2,nbmax)

      Integer  nt, ntfix, tri(3,ntmax), labelT(ntmax), fixedT(1)

      DATA     nvfix/5/, fixedV/1,2,3,4,5/,  nbfix/0/, ntfix/0/


c ======================================================================
c for library aniAFT: convex pentagon
c ======================================================================
      Real*8   bv(2,5), bltail(2,5)
      Integer  Nbv, Nbl, bl(7,5)

      DATA     Nbv/5/, Nbl/5/
      DATA     bv/0.0,0.0, 0.0,1.0, 0.7,1.0, 1.0,0.7, 1.0,0.0/  ! 5 nodes
      DATA     bl/1,2,0,-1,1,1,0, 2,3,0,-1,2,1,0,  ! boundary edges
     &            3,4,0,-1,3,1,0, 4,5,0,-1,4,1,0,
     &            5,1,0,-1,5,1,0/    
      DATA     bltail/0,0, 0,0, 0,0, 0,0, 0,0/  ! data for curved edges

      Integer  aft2dboundary
      EXTERNAL aft2dboundary


c ==========================================================
c for library aniFEM
c ==========================================================
c  include constants from Package src/aniFEM
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax)
      Real*8    A(namax), RHS(nvmax), SOL(nvmax), RES(nvmax)

      Integer  iDATAFEM(1), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATAFEM(1)

      Integer  Ddiff, Drhs, Dbc, nRow, nCol
      EXTERNAL FEM2Dext, Ddiff, Drhs, Dbc


c ======================================================================
c for library aniILU
c ======================================================================
      EXTERNAL matvec, prevec2
      Integer  imatvec(1), iprevec(1)

      Integer  iter, info, nunit, verb, UsedWr, UsedWi
      Real*8   resid, tau1, tau2, partlur, partlurout 


c ======================================================================
c for library aniLMR  (allocate data for metric)
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
      Integer  i,k,n, ip1,ip2,ip3, ir,ic, ir1,ir2,ir3, nr, label, iERR
      Integer  iIRE, inEP, iIEP, iiW, ipBCG, MaxWj, MaxWd, iLoop, nLOOPs
      Real*8   h, localA(4, 4), localF(4), rmax

      Real*8   edge_length
      External edge_length
   
c ======================================================================
c number of adaptive loops
      nLOOPs = 10


c === Generate a quasi-uniform mesh with mesh step h using library AniAFT
      h = 0.1
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary( 
     &       Nbv, bv, Nbl, bl, bltail, h,    ! input geometric data
     &       nv, vrt,                        ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)
 
      If(iERR.NE.0) Then
         Write(*,*) 'Error in function aft2dboundary'
         Stop 911
      End if

      Write(*,'(A,3I6)') 'Mesh nodes/bnd edges/triangles:', nv, nb, nt
      Write(*,'(A,3I6)') 'Max  nodes/bnd edges/triangles:', 
     &                   nvmax, nbmax, ntmax
      If(nv.GT.nvmax .OR. nb.GT.nbmax .OR. nt.GT.ntmax) Then
         Write(*,*) 'Too many nodes/boundary edges/triangles'
         Stop 911
      End if

          
c === Draw the initial mesh using library AniVIEW. 
c The name must terminate with .ps
c     Call graph(nv,vrt,nt,tri,'mesh_initial.ps')
      Call graph_demo(nv,vrt,nt,tri,'mesh_initial.ps',
     &               'Initial quasi-uniform mesh')


c Begin adaptive iterative loop.
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c === Generate matrix for Lagrange multiplies usign library AniFEM.
c No data is provided for the user subroutine Ddiff.
         dDATAFEM(1) = 0
         iDATAFEM(1) = 0

c Fill labels of mesh points with zeros
         Do i = 1, nv
            labelV(i) = 0
         End do
 
         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_AMG)
         controlFEM(2) = 1

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2Dext, dDATAFEM, iDATAFEM, controlFEM,
     &        nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c === Solve the problem with matrix A using library AniILU.
c Initialization of the incomplete LU preconditioner
         verb    = 0     ! verbose no
         tau1    = 1d-2  ! absolute threshold for L,U
         tau2    = 1d-3  ! absolute threshold for T,R
         partlur = 0.5   ! even partition of memory between LU and R
         iERR    = 0     ! error code

         Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &              rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &              UsedWr, UsedWi, iERR)

         If(iERR.NE.0) Then
            Write(*,*) 'Initialization(1) of iluoo failed, iERR=', iERR
            Stop 911
         End if

c solve the problem
         If(UsedWr + 8*nRow.GT.MaxWr) Then
            Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr + 8*nRow
            Stop 911
         End if
         ipBCG = UsedWr + 1

         iter  = 10000   ! max number of iterations
         resid = 1d-12   ! final residual
         info  = 0       ! no troubles on input
         nunit = 6       ! output to display 
         Do i = 1, nRow  ! initial guess
            SOL(i) = 0d0
         End do
 
         iprevec(1) = nRow
         imatvec(1) = nRow
         Call slpbcgs(prevec2, iprevec, iW, rW,
     &                matvec,  imatvec, IA, JA, A,
     &                rW(ipBCG), nRow, 8,
     &                nRow, RHS, SOL,
     &                iter, resid, info, nunit)

         If(info.NE.0) Stop 'BiCGStab had failed'

c Verify the residual
         Call mulAcsr(nRow, IA, JA, A, SOL, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(5X,A,E12.6)') 'maximal norm of residual: ', rmax


c === Interpolate solution at mesh vertices using library AniMBA.
         iIRE = 1
         inEP = iIRE + 3 * nt
         iIEP = inEP + nv
         iiW  = iIEP + 3 * nt

         If(iiW .gt.MaxWi) Stop 'Please increase MaxWi'
         If(nCol.gt.MaxWr) Stop 'Please increase MaxWr'

         Call listE2R(nv, nr, nt, tri, iW(iIRE), iW(inEP), iW(iIEP))

         Do n = 1, nCol
            rW(n) = SOL(n)
         End do

         MaxWj = MaxWi - inEP 
         MaxWd = MaxWr - nCol

c map a discontinuous piece-wise linear function with d.o.f on edges
c to a continuous piece-wise linear function with d.o.f. at vertices
         Call DG2P1(nv,nt, vrt, tri, iW(iIRE), rW, 
     &              SOL, MaxWj, iW(inEP), iERR)
         If(iERR.NE.0) Stop 'Please increase MaxWi'

c restore boundary conditions
         Do n = 1, nb 
            label = labelB(n)

            Do i = 1, 2
               ip1 = bnd(i, n)
               k = Dbc(vrt(1,ip1), vrt(2,ip1), 
     &                 label, dDATAFEM, iDATAFEM, iSYS, SOL(ip1))
            End do
         End do


c === Draw solution isolines in a PostScript file using library aniVIEW.
         If(iLoop.EQ.1) Then
c           Call isolines(SOL, nv,vrt, nt,tri, nb,bnd, 'iso.ps',20)
            Call isolines_demo(SOL, nv,vrt, nt,tri, nb,bnd, 'iso.ps',
     &           20, 'Solution isolines on the initial mesh')
         Else
c           Call isolines(SOL, nv,vrt, nt,tri, nb,bnd, 'his.ps',20)
            Call isolines_demo(SOL, nv,vrt, nt,tri, nb,bnd, 'his.ps',
     &           20, 'Solution isolines on the adapted mesh')
         End if


c === Generate a metric from solution SOL using library AniLMR.
c The metric is optimal for the L_p norm of error.
         Lp = 0             ! maximum norm
c        Lp = 1             ! L_1 norm
         Call Nodal2MetricVAR(SOL,
     &                        vrt, nv, tri, nt, bnd, nb, Metric,
     &                        MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c === Generate a mesh quasi-uniform in a metric using library AniMBA.
         If(iLoop.EQ.nLOOPs) Goto 500  ! Skip last adaptation

         control(1) = 100    ! MaxSkipE
         control(2) = 10000  ! MaxQItr
         control(3) = 1      ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1      ! flagAuto
         control(5) = 1      ! iPrint = minimal level of output information
         control(6) = 0      ! iErrMesgt: only critical termination allowed

         nEStar = 5000
         Quality = 0.7

         Call mbaNodal(
     &        nv, nvfix, nvmax, vrt, labelV, fixedV,  ! mesh
     &        nb, nbfix, nbmax, bnd, labelB, fixedB,
     &        nc,               crv, labelC, ANI_CrvFunction,
     &        nt, ntfix, ntmax, tri, labelT, fixedT,
     &        nEStar, Quality, control, Metric,  ! control parameters
     &        MaxWr, MaxWi, rW, iW, iERR)  ! work arrays

         If(iERR.GT.1000) Call errMesMBA(iERR, 'main',
     &                        'unspecified error in mbaNodal')


c === Draw the adapted mesh using library AniVIEW.
c        Call graph(nv,vrt,nt,tri,'mesh_final.ps')
         Call graph_demo(nv,vrt,nt,tri,'mesh_final.ps',
     &                  'Final mesh adapted to solution')

 500     Continue
      End do


c === Compute velocities u on mesh edges and cell-centered pressure p 
c using library AniFEM. This is provided only for illustration purpose.
      Call listE2R(nv, nr, nt, tri, iW(iIRE), iW(inEP), iW(iIEP))

      Do n = 1, nt
         ip1 = tri(1, n)
         ip2 = tri(2, n)
         ip3 = tri(3, n)

         ir1 = iW(iIRE + 3*(n-1))
         ir2 = iW(iIRE + 3*(n-1) + 1)
         ir3 = iW(iIRE + 3*(n-1) + 2)

c compute mass matrix using material properties
c we use the user-specified routine Ddiff()
         Call encodeISYS(n, ip1,ip2,ip2, ir1,ir2,ir3, 
     &                   0,0,0, nv,nr,nt, iSYS)

         Call fem2Dtri(vrt(1, ip1), vrt(1, ip2), vrt(1, ip3),
     &                 IDEN, FEM_RT0, IDEN, FEM_RT0,
     &                 labelT(n), Ddiff, dDATAFEM, iDATAFEM, iSYS, 2,
     &                 4, localA, ir, ic)

c add mass conservation equation
c we use the divergent theorem here as in forlibfem.f
         localA(1, 4) = edge_length(vrt(1, ip1), vrt(1, ip2))
         localA(2, 4) = edge_length(vrt(1, ip2), vrt(1, ip3))
         localA(3, 4) = edge_length(vrt(1, ip3), vrt(1, ip1))

c compute right hand side (Lagrange muptipliers are in rW)
         Do i = 1, 3
            ir1 = iW(iIRE + 3*(n-1) + i-1)
            localF(i) = localA(i, 4) * rW(ir1)
         End do

c we use the user-specified routine Drhs()
         Call fem2Dtri(vrt(1, ip1), vrt(1, ip2), vrt(1, ip3),
     &                 IDEN, FEM_P0, IDEN, FEM_P0,
     &                 labelT(n), Drhs, dDATAFEM, iDATAFEM, iSYS, 2,
     &                 1, localF(4), ir, ic)

c Solve the 4x4 sytem localA x = localF using LAPACK. The first three 
c entries in the solution x will be velocities on triangle edges. 
c The 4-th entry is the solution p at the mass center of the triangle. 
         Call dsysv('U', 4, 1, localA, 4, iW, localF, 4, rW, 40, info)
      End do


c === Testing the results for CTest.
      If(Quality.LT.0.6) Stop 911
      Stop 
      End

