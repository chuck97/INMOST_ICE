c =====================================================================
      Program UnsteadyConvectionDiffusion
c =====================================================================
c The finite element solution of the boundary value problem:
c
c  du/dt - div (D grad u) + vx * du/dx + vy * du/dy = 0  in  Omega
c                                                u = g  on  dOmega
c
c where Omega  = [0;1]x[0;1] 
c           D  = diag{0.0001, 0.0001} 
c           vx = 1
c           vy = 0
c
c SUPG stabilization is used in this example.
c =====================================================================
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
      Parameter(MaxWr = 3 000 000, MaxWi = 1 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


c ======================================================================
c Mesh definition
c ======================================================================
      Integer  nv, nvfix, labelV(nvmax), fixedV(1)
      Real*8   vrt(2,nvmax)

      Integer  nb, nbfix, bnd(2,nbmax), labelB(nbmax), fixedB(1)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2,nbmax)

      Integer  nt, ntfix, tri(3,ntmax), labelT(ntmax), fixedT(1)

      DATA     nvfix/0/,  nbfix/0/, ntfix/0/


c =====================================================================
c for library aniAFT  (we define domain Omega)
c =====================================================================
      Real*8   bv(2,6), bltail(2,6)
      integer  Nbv, Nbl, bl(7,6)

      Data     Nbv/6/, Nbl/6/
      Data     bv/0,0, 0,0.25, 0,0.75, 0,1, 1,1, 1,0 /! boundary nodes
      Data     bl/1,2,0,-1,1,1,0, 2,3,0,-1,2,1,0,   ! outer boundary edge
     &            3,4,0,-1,1,1,0, 4,5,0,-1,1,1,0,   ! outer boundary edge
     &            5,6,0,-1,1,1,0, 6,1,0,-1,1,1,0 /  ! outer boundary edge
      Data     bltail/0,0,0,0,0,0,0,0,0,0,0,0/      ! curved  boundary edges

      Integer   aft2dboundary
      EXTERNAL  aft2dboundary


c =====================================================================
c for library aniFEM
c =====================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax)
      Real*8    A(namax), F(nvmax), G(nvmax)

      Integer  IB(nvmax), JB(namax)
      Real*8    B(namax), U(nvmax), Uprev(nvmax)

      Integer  iDATAFEM(1), iDATAMASS(1), iSYS(MAXiSYS), controlFEM(3)
      Real*8   dDATAFEM(3+ntmax), dDATAMASS(1)

      Real*8   edge_length
      Integer  Dbc 
      EXTERNAL Dbc, FEM2DextStif, FEM2DextMass, edge_length


c ======================================================================
c for library aniLU
c ======================================================================
      Integer  symbolic(2), numeric(2), sys
      Real*8   controlLU(20), infoLU(90)


c LOCAL VARIABLEs
      Integer   i, i1,i2,i3, ibc, iLoop, nLOOPs, iERR, iTime, nRow,nCol
      Real*8    h, x,y, eBC(2), Time, DeltaT, ta,tb

c ======================================================================
c time integration step
      DeltaT = 0.005D0

c === generate a quasi-uniform mesh with mesh step h
      h = 0.015
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary(
     &       Nbv, bv, Nbl, bl, bltail,  h,  ! input geometric data
     &       nv, vrt,                       ! output mesh data
     &       nt, tri, labelT,
     &       nb, bnd, labelB,
     &       nc, crv, labelC)

      If(iERR.NE.0)   Stop 'Error in function aft2dboundary'
      If(nv.GT.nvmax) Stop 'Too many nodes, increase nvmax'
      If(nt.GT.ntmax) Stop 'Too many triangles, increase ntmax'
      If(nb.GT.nbmax) Stop 'Too many boundary edges, increase nbmax'
          
      Write(*,'(A,3(I6,A))') '     mesh with', nt, ' triangles', 
     &                       nv, ' vertices ', nb, ' boundary edges' 

c === draw the initial mesh. The name must terminate with .ps
c     demo graphics has been activated
c     Call graph(nv, vrt, nt, tri, 'mesh_initial.ps')
      Call graph_demo(nv, vrt, nt, tri, 'mesh_initial.ps',
     &               'Quasi-uniform mesh')


c  ===  assemble the stiffness matrix
c mark the Dirichlet points via with maximal edge color
      Call markDIR(nv, vrt, labelV, nb, bnd, labelB, 
     &             Dbc, dDATAFEM, iDATAFEM, iSYS)

c data provided for the user subroutine Dconv used in FEM2DextStif
      dDATAFEM(1) = 1  ! v_x
      dDATAFEM(2) = 0  ! v_y
      dDATAFEM(3) = DeltaT  

      Do i = 1, nt
         dDATAFEM(3+i) = 0D0
      End do

c general sparce matrix in the AMG format
      controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_CSC)
      controlFEM(2) = 1

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &     FEM2DextStif, dDATAFEM, iDATAFEM, controlFEM,
     &     nvmax, namax, IB, JB, B, F, nRow, nCol,
     &     MaxWi, MaxWr, iW, rW)


c === call the driver for LU factorization and solution
c converting the AMG row format to the CSC 0-based column format
      Call CSC2CSC0(nCol,IB,JB)

c set up default control parameters & print only error messages
      Call umf4def(controlLU)
      controlLU(1) = 1

c pre-order and symbolic analysis
      Call umf4sym(nCol, nCol, IB,JB,B, symbolic,controlLU,infoLU)
      If(infoLU(1).lt.0) Goto 5001

c numeric factorization
      Call umf4num(IB,JB,B, symbolic,numeric, controlLU,infoLU)
      If(infoLU(1).lt.0) Goto 5002

c free the symbolic analysis data
      Call umf4fsym(symbolic)


c === make the mass matrix scaled by 1.5/DeltaT
c generate bilinear form using build-in function ANI_D1x1_scalar
c general sparce matrix in the AMG format (old control is used)
      dDATAMASS(1) = 1.5d0 / DeltaT

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &     FEM2DextMass, dDATAMASS, iDATAMASS, controlFEM,
     &     nvmax, namax, IB, JB, B, U, nRow, nCol,
     &     MaxWi, MaxWr, iW, rW)


c === set the initial state U = 0
      Call dcopy(nv, 0d0, 0, U, 1)

      Do i = 1, nv
         If(labelV(i).EQ.2 .OR. labelV(i).EQ.1) Then
            x = vrt(1, i)
            y = vrt(2, i)

            ibc = Dbc(x, y, labelV(i), dDATAFEM, iDATAFEM, iSYS, eBC)

            U(i)= eBC(1) 
         End if
      End do
         
c copy U to Uprev
      Call dcopy(nv, U, 1, Uprev, 1)


c ======================================================================
c Time iteration loop
c ======================================================================
      Time  = 0
      iTime = 0

 2    Continue
         iTime = iTime + 1
         Time  = Time + DeltaT

         write(*,'(A,I3,A,F5.2)') 'Time step ', iTime, ' time=', Time

c === regenerate RHS due to SUPG
         ta = 2D0 / 3
         tb = 1D0 / 6
         Do i = 1, nt
            i1 = tri(1,i)
            i2 = tri(2,i)
            i3 = tri(3,i)

            dDATAFEM(3 + i) = ta * (U(i1) + U(i2) + U(i3))  
     &                      - tb * (Uprev(i1) + Uprev(i2) + Uprev(i3))
         End do

         controlFEM(2) = 0  ! no output

         Call BilinearFormTemplate(
     &        nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &        FEM2DextStif, dDATAFEM, iDATAFEM, controlFEM,
     &        nvmax, namax, IA, JA, A, F, nRow, nCol,
     &        MaxWi, MaxWr, iW, rW)


c === BFD time step: (1.5 u_i - 2 u_{i-1} + 0.5 u_{i-2}) / DeltaT
c       note: mass matrix is scaled with 1.5/DeltaT
         Call daxpy(nRow, -4d0, U,1, Uprev,1)
         Call dscal(nRow, -1d0/3d0, Uprev,1)
         Call mulAcsr(nRow, IB, JB, B, Uprev, G)

         Call daxpy(nRow, 1d0, F,1, G,1)
         Call dcopy(nRow, U,1, Uprev,1)


c solve Ax=b, without iterative refinement
         sys = 0
         Call umf4sol(sys, U, G, numeric, controlLU,infoLU)
         If(infoLU(1).lt.0) Goto 5003

c assign the Dirichlet boundary condition
         Do i = 1, nv
            If(labelV(i).EQ.2 .OR. labelV(i).eq.1) then
               x =  vrt(1, i)
               y =  vrt(2, i)

               ibc = Dbc(x, y, labelV(i), dDATAFEM, iDATAFEM, iSYS, eBC)
               U(i) = eBC(1) 
            End if
         End do

         If(iTime.EQ.1) 
c           Call isolines(U,nv,vrt,nt,tri,nb,bnd,'isolines_initial.ps',10)
     &      Call isolines_demo(U, nv,vrt, nt,tri, nb,bnd,
     &          'isolines_initial.ps',10,
     &          'Solution isolines at time T=0')

      If(Time.LT.0.5) Goto 2
c END OF TIME INTEGRATION

c free the numeric factorization data
      Call umf4fnum (numeric)


c === draw 10 isolines of solution. The name must terminate with .ps
c     umin = 0 and umax = 1 explicitly there!
c     Call isolines(U,nv,vrt,nt,tri,nb,bnd,'isolines_final.ps',10)
      Call isolines_demo(U, nv,vrt, nt,tri, nb,bnd,
     &                  'isolines_final.ps', 10,
     &                  'Solution isolines at time T=0.5')


c === testing the results
      Stop 


c error messages
 5001 Continue
      Write(*,*) 'Error occurred in umf4sym: ', infoLU(1)
      Stop 911

 5002 Continue
      Write(*,*) 'Error occurred in umf4num: ', infoLU(1)
      Stop 911

 5003 Continue
      Write(*,*) 'Error occurred in umf4sol: ', infoLU(1)
      Stop 911

      End





