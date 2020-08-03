c ======================================================================
      Program APosterioriEstimates
c ======================================================================
c The iterative adaptive solution for the following boudary value problem (bvp):
c
c   -div grad u = 0   in  Omega
c             u = g   on  dOmega_D
c        du/dn  = 0   on  dOmega_N
c
c where Omega = { x^2 + y^2 < 1 \ {y=0, 0<x<1} }
c    dOmega_N = { y=0, 0<x<1 } from  bottom
c    dOmega_D = { x^2 + y^2 = 1} .or. {y=0, 0<x<1}
c           g = r^0.25 sin phi/4
c
c Exact solution of this problem is u=r^0.25 sin phi/4  (in polar coordinates)
c
c The program generates optimal mesh with nEStar elements.
c We use the P_1 finite element method with edge-based a posteriori error 
c estimates.
c ======================================================================
      implicit none

c  nvmax - maximum number of mesh nodes
c  ntmax - maximum number of mesh triangles
c  nbmax - maximum number of boundary edges
c  namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,namax
      parameter(nvmax = 200 000, ntmax = 2*nvmax, nbmax = 10 000)
      parameter(namax = 3 000 000)

c  work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 3 000 000, MaxWi = 9 000 000)

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
c  define computational domain Omega
      Real*8    bv(2,6),bltail(2,6)
      Integer   Nbv,Nbl,bl(7,6)

      Data      Nbv/6/,Nbl/6/
      Data      bv/0,-1, -1,0, 0,1, 1,0, 0,0, 1,0/    ! boundary nodes (last coord.pare is fictitious)
      Data      bl/1,2,1,-1,1,1,0, 2,3,1,-1,2,1,0,    ! outer boundary edges
     &             3,4,1,-1,3,1,0, 4,5,0,-1,5,1,0,    ! outer boundary edges
     &             5,6,0,-1,6,1,0, 6,1,1,-1,4,1,0/    ! outer boundary edges
      Data      bltail/1.5,1, 1,0.5, 0.5,0, 0,0,0,0,1.999999,1.5/   ! data for each curved boundary edge

      Integer   aft2dboundary
      EXTERNAL  aft2dboundary
      EXTERNAL  userboundary      


c ======================================================================
c for ibrary aniFEM
c ======================================================================
c  include constants from Package src/aniFEM
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax)
      Real*8    A(namax), RHS(namax), SOL(namax), RES(namax)
      Real*8   SOL_L(namax)

      Integer  iDATA(1), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATA(1)

      Integer  Ddiff, Drhs, Dbc
      EXTERNAL Ddiff, Drhs, Dbc
      EXTERNAL FEM2DextLinear, FEM2DextBubble, FEM2DextFull

      Integer  DgradU
      EXTERNAL DgradU


c ======================================================================
c for library aniILU
c ======================================================================
      EXTERNAL matvec, prevec2
      Integer  imatvec(1), iprevec(1)

      Integer  iter, info, nunit, verb, UsedWr, UsedWi
      Real*8   resid, tau1,tau2,partlur,partlurout 


c ======================================================================
c for library aniLMR  (we allocate data for metric)
c ======================================================================
      Real*8   Lp
      Real*8   Error(3,ntmax), Metric(3,nvmax)


c ======================================================================
c for library aniMBA  (we define adaptation strategy here)
c ======================================================================
      Integer   control(6), nEStar
      Real*8    Quality
      EXTERNAL  CrvFunction_user

c LOCAL VARIABLEs
      Integer  i,j,k,n, iv1,iv2,iv3, iLoop, ipBCG, label, iERR, nLOOPs
      Integer  nRow, nCol
      Real*8   rmax, enemin, SOLloc(6), trierr, h1global
      Real*8   h

c ======================================================================
c number of adaptive loops
      nLOOPs = 12

c === pass the user function userBoundary (forlibaft.c)
      Call registeruserfn(userboundary)  ! register the name to be used in the library

c generate a quasi-uniform mesh with mesh step h
      h = 0.035
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary( 
     &       Nbv, bv, Nbl, bl, bltail, h,    ! input geometric data
     &       nv, vrt,                        ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)

      If(iERR.NE.0) Stop 'Error in function aft2dboundary'
      If(nv.GT.nvmax) stop ' too many nodes, increase nvmax'
      If(nt.Gt.ntmax) stop ' too many triangles, increase ntmax'
      If(nb.GT.nbmax) stop ' too many boundary edges, increase nbmax'
 
      Write(*,'(A,3(I6,A))') '     mesh with', nt, ' triangles', 
     &                       nv, ' vertices ', nb, ' boundary edges' 

          
c === draw the initial mesh. The name must terminate with .ps
      Call graph(nv,vrt, nt,tri, 'mesh_initial.ps')


c begin adaptive iterative loop
      enemin = 1D+24
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c === no data is provided for the user subroutine Ddiff
         dDATA(1) = 0D0
         iDATA(1) = 0

c        mark Dirichlet boundary points with maximal edge color
         Do i = 1, nv
            labelV(i) = 0
         End do

         Do i = 1, nb
            If(labelB(i).LE.5) Then
               labelV(bnd(1,i)) = max(labelV(bnd(1,i)), labelB(i))
               labelV(bnd(2,i)) = max(labelV(bnd(2,i)), labelB(i))
            End if
         End do
 
c        general sparce matrix in the AMG format 
         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_AMG)
         controlFEM(2) = 1

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2DextLinear, dDATA, iDATA, controlFEM,
     &        nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c === ILU SOLVER 
c        initialization of the preconditioner
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

c        iterative solution
         If(UsedWr + 8*nRow.GT.MaxWr) Then
            Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr + 8*nRow
            Stop 911
         End if

         ipBCG = UsedWr + 1

         iter  = 10000  ! max number of iterations
         resid = 1d-12  ! final residual
         info  = 0      ! no troubles on input
         nunit = 6      ! output to display 
         Do i = 1, nRow ! initial guess
            SOL_L(i) = 0d0
         End do
 
         iprevec(1) = nRow
         imatvec(1) = nRow
         Call slpbcgs(prevec2, iprevec, iW, rW,
     &                matvec,  imatvec, IA, JA, A,
     &                rW(ipBCG), nRow, 8,
     &                nRow, RHS, SOL_L,
     &                iter, resid, info, nunit)

         If(info.NE.0) Stop 'BiCGStab had failed'

c        check the residual
         Call mulAcsr(nRow, IA, JA, A, SOL_L, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(5X,A,E12.6)') 'maximal norm of residual: ', rmax


c === plot solution isolines 
         If(iLoop.EQ.1) then
c           Call isolines(SOL_L, nv,vrt, nt,tri, nb,bnd, 'iso.ps',20)
            Call isolines_demo(SOL_L, nv,vrt, nt,tri, nb,bnd, 'iso.ps',
     &                    20, 'Solution isolines on the initial mesh')
         Else
c           Call isolines(SOL_L, nv,vrt, nt,tri, nb,bnd, 'his.ps',20)
            Call isolines_demo(SOL_L, nv,vrt, nt,tri, nb,bnd, 'his.ps',
     &                    20, 'Solution isolines on the final mesh')
         End if



c === error calculation with library aniFEM
         Lp = 2
         h1global = 0D0

         Do n = 1, nt
            label = labelT(n)

            Do i = 1, 3
               iv1 = tri(i, n)
               SOLloc(i) = SOL_L(iv1)
            End do

c           compute H1 error for Uh 
            iv1 = tri(1, n)
            iv2 = tri(2, n)
            iv3 = tri(3, n)
            Call fem2Derr(vrt(1,iv1), vrt(1,iv2), vrt(1,iv3), Lp,
     &                    GRAD, FEM_P1, SOLloc, DgradU, dDATA, iDATA,
     &                    label, Ddiff, dDATA, iDATA, iSYS, 5, trierr)
 
            h1global = h1global + trierr
         End do
         Write(*,'(A,E12.4)') 'FEM: true H1-error =', dsqrt(h1global)


c === A POSTERIORI ERROR ESTIMATES (GENERATE A SYSTEMS FOR BUBBLES)
c        general sparce matrix in the AMG format (old control is used)
         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2DextBubble, SOL_L, iDATA, controlFEM,
     &        nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)

         Write(*,'(/,A,I8)') 
     &     'FEM: number of edges (= matrix size):', nCol


c === A POSTERIORI ERROR ESTIMATES (SOLVE A SYSTEMS FOR BUBBLES)
c        initialization of the preconditioner (old parameters are used)
         Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &              rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &              UsedWr, UsedWi, iERR)

         If(iERR.NE.0) Then
            Write(*,*) 'Initialization(2) of iluoo failed, iERR=', iERR
            Stop
         End if

c        iterative solution
         If(UsedWr + 8*nRow.GT.MaxWr) Then
            Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr + 8*nRow
            Stop
         End if

         ipBCG = UsedWr + 1

         iter  = 10000  ! max number of iterations
         resid = 1d-8  ! final residual
         info  = 0      ! no troubles on input
         nunit = 6      ! output to display 
         Do i = 1, nRow ! initial guess
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

c        check the residual
         Call mulAcsr(nRow, IA, JA, A, SOL, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(5X,A,E12.6)') 'maximal norm of residual: ', rmax

c === compute error estimate
         rmax = 0
         Do i = 1, nRow
            rmax = rmax + SOL(i) * RHS(i)
         End do
         rmax = dsqrt(rmax)

         If(rmax.LT.enemin) Then
            enemin = rmax
         End if

         Write(*,'(5X,2(A,E11.5))') 
     &      'estimated energy norm of error = ', rmax, 
     &              ',  minimal in the loop = ', enemin


c === GENERATE A METRIC (from SOL) optimal for the L_p norm 
         Lp = 2

c        copy data from edges to triangles for analysis  
         Do i = 1, nRow
            SOL(i) = dabs(SOL(i))
         End do

         Call copyEdge2Triangle(nv, nt, tri, SOL, Error, iW)

         Call EdgeEst2GradMetricMAX(Error, nv, nt, vrt, tri, 
     &                              Metric, MaxWr, rW)

         If(Lp.GT.0) Call Lp_gradnorm(nv, Lp, Metric)


c === A POSTERIORI ERROR ESTIMATES (GENERATE A FULL SYSTEMS)
         Goto 1000

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2DextFull, SOL_L, iDATA, controlFEM,
     &        nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c === A POSTERIORI ERROR ESTIMATES (SOLVE THE FULL SYSTEMS)
c        initialization of the preconditioner 
         verb    = 0     ! verbose no
         tau1    = 1d-2  ! absolute threshold for L,U
         tau2    = 1d-3  ! absolute threshold for T,R
         partlur = 0.5   ! even partition of memory between LU and R
         iERR    = 0     ! error code

         Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &              rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &              UsedWr, UsedWi, iERR)

         If(iERR.NE.0) Then
            Write(*,*) 'Initialization(3) of iluoo failed, iERR=', iERR
            Stop
         End if

c        iterative solution
         If(UsedWr + 8*nRow.GT.MaxWr) Then
            Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr + 8*nRow
            Stop
         End if

         ipBCG = UsedWr + 1

         iter  = 2000   ! max number of iterations
         resid = 1d-6   ! final residual
         info  = 0      ! no troubles on input
         nunit = 6      ! output to display 
         Do i = 1, nRow ! initial guess
            SOL(i) = 0d0
         End do
 
         iprevec(1) = nRow
         imatvec(1) = nRow
         Call slpbcgs(prevec2, iprevec, iW, rW,
     &                matvec,  imatvec, IA, JA, A,
     &                rW(ipBCG), nRow, 8,
     &                nRow, RHS, SOL,
     &                iter, resid, info, nunit)

         If(info.NE.0) Write(*,*) '  BiCGStab had failed'

c        check the residual
         Call mulAgen(nRow, IA, JA, A, SOL, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(5X,A,E12.6)') 'maximal norm of residual: ', rmax

c === estimate the error
         rmax = 0
         Do k = nv + 1, nRow
            Do i = IA(k), IA(k + 1) - 1
               j = JA(i)
               If(j.GT.nv) rmax = rmax + A(i) * SOL(j) * SOL(k)
            End do
         End do

         Write(*,'(5X,A,E11.5)') 
     &      'estimated energy norm of error = ', dsqrt(rmax)


c === GENERATE ADAPTIVE MESH using Metric
 1000    Continue
         If(iLoop.EQ.nLOOPs) Goto 500  !Skip last adaptation

         control(1) = 400     ! MaxSkipE
         control(2) = 20000   ! MaxQItr
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


c === draw the adapted mesh
c        Call graph(nv,vrt, nt,tri, 'mesh_final.ps')
         Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps',
     &                  'Final mesh adapted to the solution')

 500     Continue
      End do


c === testing the results
      If(Quality.LT.0.3) Stop 911

      Stop 
      End


