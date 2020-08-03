c ======================================================================
      Program DiscontinuousGalerkin
c ======================================================================
c The  solution for the following boudary value problem (bvp):
c
c   -div  grad u = -2  in  Omega
c              u = g   on  dOmega_D
c
c where Omega = { x^2 + y^2 < 1 } and  g = x*x
c
c Exact solution of this problem is u=x*x. To solve this problem we use
c a mixture of P2 and P3 finite element. For triangles with centers 
c belonging to the disk or radius 0.5, we use the quadratic finite 
c elements. For the remaining triangles, we use the cubic finite elements. 
c These elements represent the solution exactly, therefore the expected
c error should be zero (see the output). 
c ======================================================================
      implicit none

c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
c namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,numax,namax
      parameter(nvmax = 150 000, ntmax = 2*nvmax, nbmax = 10 000)
      parameter(numax = 3*ntmax, namax = 5 000 000)

c work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 3 000 000, MaxWi = 5 000 000)

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
c for library aniAFT  (we define domain Omega)
c ======================================================================
      Real*8    bv(2,4), bltail(2,4)
      Integer   Nbv, Nbl, bl(7,4)

      DATA      Nbv/4/, Nbl/4/
      DATA      bv/0,-1, -1,0, 0,1, 1,0/    ! boundary nodes (last coord.pare is fictitious)
      DATA      bl/1,2,1,-1,1,1,0, 2,3,1,-1,2,1,0,    ! outer boundary edges
     &             3,4,1,-1,3,1,0, 4,1,1,-1,4,1,0/    ! outer boundary edges
      DATA      bltail/1.5,1, 1,0.5, 0.5,0, 1.99999999,1.5/   ! data for each curved boundary edge

      Integer   aft2dboundary
      EXTERNAL  aft2dboundary
      EXTERNAL  userboundary      


c ======================================================================
c for library aniFEM
c ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(numax), JA(namax)
      Real*8    A(namax), RHS(numax), SOL(numax), RES(numax)
      Real*8   SOL_V(numax)

      Integer  iDATA(ntmax), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATA(1)

      EXTERNAL FEM2DextDG


c ======================================================================
c for library aniILU
c ======================================================================
      EXTERNAL matvec, prevec2
      Integer  imatvec(1), iprevec(1)

      Integer  iter, info, nunit, verb, UsedWr, UsedWi
      Real*8   resid, tau1,tau2,partlur,partlurout 


c LOCAL VARIABLEs
      Integer  nRow, nCol, iERR
      Integer  i,k,n, iLoop, iv1,iv2,iv3, label, ipBCG, iEnd
      Real*8   trierr, errLp, rmax, h, x,y, r, U(4)

      Integer  ANI_Dnull, Dbc, Dexact
      EXTERNAL ANI_Dnull, Dbc, Dexact

c ======================================================================
c === call library AniAFT to generate a mesh
c     pass the user function userBoundary (forlibaft.c)
      Call registeruserfn(userboundary)  ! register the name to be used in the library


c     generate a quasi-uniform mesh with mesh step h
      h = 0.1
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary( 
     &       Nbv, bv, Nbl, bl, bltail, h,    ! input geometric data
     &       nv, vrt,                        ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)

      If(iERR.NE.0) Stop 'Error in function aft2dboundary'
      If(nv.GT.nvmax) stop ' too many nodes, increase nvmax'
      If(nt.GT.ntmax) stop ' too many triangles, increase ntmax'
      If(nb.GT.nbmax) stop ' too many boundary edges, increase nbmax'
          
c     Call loadMani(nv, nvfix, nvmax, vrt, labelV, fixedV,
c    &              nb, nbfix, nbmax, bnd, labelB, fixedB,
c    &              nc,               crv, labelC, 
c    &              nt, ntfix, ntmax, tri, labelT, fixedT,
c    &              "../data/simple.ani")

      Write(*,'(A,3(I6,A))') '     mesh with', nt, ' triangles', 
     &                       nv, ' vertices ', nb, ' boundary edges' 

          
c === draw Postscript file of the initial mesh, name must have extension .ps
c     Demo graphics as been activated
c     Call graph(nv,vrt, nt,tri, 'mesh.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh.ps',
     &               'Quasi-uniform mesh')


c === no real*8 data is provided for the user subroutine FEM2DextDG
      dDATA(1) = 0D0
c     polynomial order is provided to user routines
      Do n = 1, nt
         iv1 = tri(1, n)
         iv2 = tri(2, n)
         iv3 = tri(3, n)
         Call tri_center(vrt(1,iv1), vrt(1,iv2), vrt(1,iv3), x, y)

         r = dsqrt(x*x + y*y)
         iDATA(n) = FEM_P3
         If(r.LE.0.5D0) iDATA(n) = FEM_P2
      End do

c mark the Dirichlet points with the maximal edge color
      Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &             Dbc, dDATA, iDATA, iSYS)

c general sparce matrix in the compressed row format
      controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_AMG)
      controlFEM(2) = 1

      Call AssembleFromTemplate(
     &     nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &     FEM2DextDG, dDATA, iDATA, controlFEM,
     &     numax, namax, IA, JA, A, RHS, nRow, nCol,
     &     MaxWi, MaxWr, iW, rW)
         
      Call draw_matrix(nRow, IA, JA, "matrix.ps")


c === ILU SOLVER 
c initialization of the preconditioner
      verb    = 0     ! verbose no
      tau1    = 1d-2  ! absolute threshold for L,U
      tau2    = 1d-3  ! absolute threshold for T,R
      partlur = 0.5   ! even partition of memory between LU and R
      iERR    = 0     ! error code

      Call iluoo(nRow, IA, JA, A, tau1, tau2, verb,
     &           rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &           UsedWr, UsedWi, iERR)

      If(iERR.NE.0) Goto 5001

c iterative solution
      iEnd = UsedWr + 8*nRow
      If(iEnd.GT.MaxWr) Goto 5002

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
     &             matvec,  imatvec, IA, JA, A,
     &             rW(ipBCG), nRow, 8,
     &             nRow, RHS, SOL,
     &             iter, resid, info, nunit)

      If(info.NE.0) Stop 'BiCGStab had failed'

c check the residual
      Call mulAcsr(nRow, IA, JA, A, SOL, RES)
      rmax = 0
      Do i = 1, nRow
         rmax = max(rmax, RES(i) - RHS(i))
      End do
      Write(*,'(A,E12.6)') 'ILU: maximal norm of residual: ', rmax


c === extract solution at points 
      Call DG_nodal(nv, nt, tri, iDATA, SOL, SOL_V)

      errLp = 0D0
      Do n = 1, nt
         Do i = 1, 3
            iv1 = tri(i, n)
            x = vrt(1, iv1)
            y = vrt(2, iv1)
            k = Dexact(x, y, 1, dDATA, iDATA, iSYS, U)

            errLp = errLp + (SOL_V(iv1) - U(1)) ** 2
         End do
      End do
      Write(*,'(A,E12.5)') 
     &    'FEM: DEBUG: L2 norm of error:', sqrt(errLp / nt)

c === make PostScript file with solution isolines.
c The name must have extension .ps
c     Call isolines(SOL_V, nv,vrt, nt,tri, nb,bnd, 'iso.ps',20)
      Call isolines_demo(SOL_V, nv,vrt, nt,tri, nb,bnd,'iso.ps',
     &                   20, 'Solution isolines')


c === testing the results
      Stop 


 5001 Continue
      Write(*,*) 'Initialization(1) of iluoo failed, iERR=', iERR
      Stop 911

 5002 Continue
      Write(*,'(A,I7)') 'Increase MaxWr to ', iEnd
      Stop 911

      End




