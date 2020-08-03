c ======================================================================
      Program APosterioriEstimates
c ======================================================================
C The program demonstrates the iterative adaptive solution of
C the boundary value problem:
C
C   -div K grad p = f    in  Omega   
C               p = g    on  dOmega
C
C where Omega is a bounded domain. We shall refer to p as the pressure. 
c We use the P_1 finite element method with edge-based a posteriori error 
c estimates.
c ======================================================================
      implicit none

c  nvmax - maximum number of mesh nodes
c  ntmax - maximum number of mesh triangles
c  nbmax - maximum number of boundary edges
c  namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,namax
      parameter(nvmax = 150 000, ntmax = 2*nvmax, nbmax = 50 000)
      parameter(namax = 2 000 000)

c  work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 9 000 000, MaxWi = 9 000 000)

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
c for Library AFT  (we define computational domain Omega here)
c ======================================================================
      Real*8   bv(2,4), bltail(2,4)
      Integer  Nbv, Nbl, bl(7,4)

      data     Nbv/4/, Nbl/4/
      data     bv/0,0, 0,1, 1,1, 1,0/                ! boundary nodes
      data     bl/1,2,0,-1,1,1,0, 2,3,0,-1,2,1,0,    ! boundary edges
     &            3,4,0,-1,3,1,0, 4,1,0,-1,4,1,0/    
      data     bltail/0,0, 0,0, 0,0, 0,0/            ! curved data for boundary edge

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

      EXTERNAL FEM2DextLinear, FEM2DextBubble, Ddiff, Drhs, Dbc
      Integer                                  Ddiff, Drhs, Dbc


c ======================================================================
c for Library ILU
c ======================================================================
      External matvec, prevec2
      Integer  imatvec(1), iprevec(1)

      Integer  iter, info, nunit
      Real*8   resid

      Real*8   tau1,tau2,partlur,partlurout
      Integer  verb, UsedWr, UsedWi


c ======================================================================
c for library aniLMR
c ======================================================================
      Real*8   Lp
      Real*8   Error(3,ntmax), Metric(3,nvmax)


c ======================================================================
c for library aniMBA
c ======================================================================
      Integer  control(6), nEStar
      Real*8   Quality
      EXTERNAL ANI_CrvFunction


c LOCAL VARIABLEs
      Integer  i, iLoop, nLOOPs, nRow, nCol, iERR, ipBCG
      Real*8   rmax, h, enemin

c ======================================================================
c number of adaptive loops
      nLOOPs = 7


c === generate a quasi-uniform mesh with mesh step h
      h = 0.009
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
          
          
c === draw the initial mesh. The name must terminate with .ps
c     Call graph(nv,vrt, nt,tri, 'mesh_initial.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_initial.ps',
     &               'Initial quasi-uniform mesh')


c begin adaptive iterative loop
      enemin = 1D+24
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c ======================================================================
c SOLVE A LINEAR SYSTEM
c ======================================================================
c === no data is provided for the user subroutine Ddiff
         dDATAFEM(1) = 0D0
         iDATAFEM(1) = 0

c mark the Dirichlet points with the maximal edge color
         Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &                Dbc, dDATAFEM, iDATAFEM, iSYS)

c general sparce matrix in the AMG format
         controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_AMG)
         controlFEM(2) = 1

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2DextLinear, dDATAFEM, iDATAFEM, controlFEM,
     &        nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c === initialization of the ILU preconditioner
         verb    = 0     ! verbose no
         tau1    = 1d-2  ! absolute threshold for L,U
         tau2    = 1d-3  ! absolute threshold for T,R
         partlur = 0.5   ! even partition of memory between LU and R
         iERR    = 0     ! error code

         Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &              rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &              UsedWr, UsedWi, iERR)

         If(iERR.NE.0) Then
            Write(*,*) 'Initialization of iluoo has failed, iERR=', iERR
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
            SOL(i) = 0d0
         End do
 
         iprevec(1) = nRow
         imatvec(1) = nRow
         Call slpbcgs(prevec2, iprevec, iW, rW,
     &                matvec,  imatvec, IA, JA, A,
     &                rW(ipBCG), nRow, 8,
     &                nRow, RHS, SOL,
     &                iter, resid, info, nunit)

         If(info.NE.0) Then
            Write(*,*) 'BiCGStab had failed'
            Stop 911
         End if

c        check the residual
         Call mulAcsr(nRow, IA, JA, A, SOL, RES)
         rmax = 0
         Do i = 1, nRow
            rmax = max(rmax, RES(i) - RHS(i))
         End do
         Write(*,'(A,E12.6)') '     maximal norm of residual: ', rmax

c === draw solution isolines: demo graphics has been activated
c        Call isolines(SOL, nv,vrt, nt,tri, nb,bnd, 'iso.ps',20)
         Call isolines_demo(SOL, nv,vrt, nt,tri, nb,bnd, 'iso.ps',
     &                 20, 'Solution isolines on the adapted mesh')


c ======================================================================
c SOLVE A SYSTEM FOR BUBBLES (A POSTERIORI ERROR ESTIMATES) 
c ======================================================================
c === general sparce matrix in the AMG format (old control is used)
         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2DextBubble, SOL, iDATAFEM, controlFEM,
     &        nvmax, namax, IA, JA, A, RHS, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)

c === initialization of the preconditioner (old parameters are used)
         Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &              rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &              UsedWr, UsedWi, iERR)

         If(iERR.NE.0) Then
            Write(*,*) 'Initialization of iluoo has failed, iERR=', iERR
            Stop 911
         End if

c        iterative solution
         If(UsedWr + 8*nRow.GT.MaxWr) Then
            Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr + 8*nRow
            Stop 911
         End if

         ipBCG = UsedWr + 1

         iter  = 2000   ! max number of iterations
         resid = 1d-4   ! final residual
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

         If(info.NE.0) Then
            Write(*,*) 'BiCGStab had failed'
            Stop 911
         End if

c        compute error estimate
         rmax = 0
         Do i = 1, nRow
            rmax = rmax + SOL(i) * RHS(i)
         End do
         rmax = dsqrt(rmax)

         enemin = min(enemin, rmax)  
         Write(*,'(2(A,E11.5))') 
     &      '     estimated energy norm of error = ', rmax, 
     &      ',   minimal in the loop = ', enemin

         If(iLoop.EQ.nLOOPs) Goto 500  ! Skip last adaptation


c === generate metric (from SOL) optimal for the L_p norm
c        Lp = 0             ! maximum norm
         Lp = 2             ! L_2 norm

c        copy data from edges to triangles for analysis  
         Do i = 1, nRow
            SOL(i) = dabs(SOL(i))
         End do

         Call copyEdge2Triangle(nv, nt, tri, SOL, Error, iW)

         Call EdgeEst2GradMetricMAX(Error, nv, nt, vrt, tri, 
     &                              Metric, MaxWr, rW)

         If(Lp.GT.0) Call Lp_gradnorm(nv, Lp, Metric)


c === generate the adaptive mesh
         control(1) = 400     ! MaxSkipE
         control(2) = 20000   ! MaxQItr
         control(3) = 1       ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1       ! flagAuto
         control(5) = 1       ! iPrint = minimal level of output information
         control(6) = 0       ! iErrMesgt: only critical termination allowed

         nEStar = 24000
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


c === draw the adapted mesh. The name must terminate with .ps
c        Call graph(nv,vrt, nt,tri, 'mesh_final.ps')
         Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps',
     &                  'Final mesh adapted to the solution')

 500     Continue
      End do


c === testing the results
      If(Quality.LT.0.4) Stop 911
      Stop 

      End

