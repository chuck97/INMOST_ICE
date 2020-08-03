C ======================================================================
c The program illustrates organization of an adaptive loop for building
c a mesh minimizing the gradient of the interpolation error for a given 
c function (Func). The function is taken from D'Azevedo's paper on 
c optimal triangulations. 
C ======================================================================
      program GradInterpolation
C ======================================================================
      implicit none

c nvmax   - maximum number of mesh nodes
c ntmax   - maximum number of mesh triangles
c nbmax   - maximum number of boundary edges
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
      Integer  nv, nvfix, labelV(nvmax), fixedV(4)
      Real*8   vrt(2,nvmax)

      Integer  nb, nbfix, bnd(2,nbmax), labelB(nbmax), fixedB(1)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2,nbmax)

      Integer  nt, ntfix, tri(3,ntmax), labelT(ntmax), fixedT(1)

      DATA     nvfix/4/, fixedV/1,2,3,4/,  nbfix/0/, ntfix/0/


c ======================================================================
c for library aniAFT (advanced front technique)
c ======================================================================
      Real*8   bv(2,4),bltail(2,4)
      Integer  Nbv,Nbl,bl(7,4)

      DATA     Nbv/4/,Nbl/4/
      DATA     bv/-1,-1, -1,1, 1,1, 1,-1/           ! boundary nodes
c     DATA     bv/0,0, 0,1, 1,1, 1,0/               ! boundary nodes

      DATA     bl/1,2,0,-1,1,1,0, 2,3,0,-1,2,1,0,   ! outer boundary edges
     &            3,4,0,-1,3,1,0, 4,1,0,-1,4,1,0/   ! outer boundary edges
      DATA     bltail/0,0, 0,0, 0,0, 0,0/           ! curved DATA for outer boundary edges

      Integer   aft2dboundary
      EXTERNAL  aft2dboundary


c ======================================================================
c for library aniFEM (finite element matrices)
c ======================================================================
      include  'fem2Dtri.fd'
      include  'assemble.fd'

      Real*8    dDATA(1), U(nvmax), Metric(3,nvmax)
      Integer   iDATA(1), iSYS(MAXiSYS)

      Integer   ANI_Dnull, DexactU, DgradU
      EXTERNAL  ANI_Dnull, DexactU, DgradU


c ======================================================================
c for library aniMBA
c ======================================================================
      Integer   control(6), nEStar
      Real*8    Quality
      EXTERNAL  ANI_CrvFunction


c LOCAL VARIABLEs
      Real*8    Func
      EXTERNAL  Func

      Integer   i,n, iv1,iv2,iv3, iLoop, label, nLOOPs, iERR
      Real*8    errLp,    errLpGrad, Uloc(3), trierr, x,y, Lp, h
      Real*8    errLpMin, errLpGradMin

c ======================================================================
c number of adaptive loops
      nLOOPs = 7


c === generate a quasi-uniform mesh with mesh step h
      h = 0.047
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary(Nbv, bv, Nbl, bl, bltail, h,   ! geomtric DATA
     &                   nv, vrt,                       ! mesh DATA on output
     &                     nt, tri, labelT,
     &                     nb, bnd, labelB,
     &                     nc, crv, labelC)
      If(iERR.NE.0)   Stop 'Error in function aft2dboundary'
      If(nv.GT.nvmax) Stop 'Too many nodes'
      If(nt.GT.ntmax) Stop 'Too many triangles'
      If(nb.GT.nbmax) Stop 'Too many boundary edges'

      Write(*,'(A,3(I6,A))') '     mesh with', nt, ' triangles', 
     &                       nv, ' vertices ', nb, ' boundary edges' 


c === PostScript figure of the initial mesh. The name must terminate with .ps
c Demo version is activated
c     Call graph(nv,vrt, nt,tri, 'mesh_initial.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_initial.ps',
     &                'Initial quasi-uniform mesh')


c begin adaptive iterative loop
      errLpMin     = 1D+24
      errLpGradMin = 1D+24

      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c define the function at the nodes and at mid-edges
         Do i = 1, nv
            x = vrt(1, i)
            y = vrt(2, i)

            n = DexactU(x, y, label, dDATA, iDATA, iSYS, U(i))
         End do


c === make PosctScript file with 20 isolines of interpolated mesh function U
c The name must terminate with .ps
c Demo version has been activated
         If(iLoop.eq.1) Then
c           Call isolines(U, nv,vrt, nt,tri, nb,bnd, 'iso.ps',20)
            Call isolines_demo(U, nv,vrt, nt,tri, nb,bnd, 'iso.ps',20,
     &                         'Solution isolines on the initial mesh')
         Else
c           Call isolines(U, nv,vrt, nt,tri, nb,bnd, 'his.ps',20)
            Call isolines_demo(U, nv,vrt, nt,tri, nb,bnd, 'his.ps',20,
     &                        'Solution isolines on the final mesh')
         End if


c === compute Lp-error
c        Lp = 0             ! maximum norm
         Lp = 2             ! L_2 norm

         errLp     = 0D0
         errLpGrad = 0D0
         Do n = 1, nt
            label = labelT(n)

            Do i = 1, 3
               iv1 = tri(i, n)
               Uloc(i) = U(iv1)
            End do

            iv1 = tri(1, n)
            iv2 = tri(2, n)
            iv3 = tri(3, n)
            Call fem2Derr(vrt(1,iv1), vrt(1,iv2), vrt(1,iv3), Lp,
     &                    IDEN, FEM_P1, Uloc, DexactU, dDATA, iDATA,
     &                    label, ANI_Dnull, dDATA,iDATA,iSYS, 5, trierr)

            errLp = errLp + trierr

            Call fem2Derr(vrt(1,iv1), vrt(1,iv2), vrt(1,iv3), Lp,
     &                    GRAD, FEM_P1, Uloc, DgradU, dDATA, iDATA,
     &                    label, ANI_Dnull, dDATA,iDATA,iSYS, 5, trierr)

            errLpGrad = errLpGrad + trierr
         End do
         errLp = dsqrt(errLp)
         errLpGrad = dsqrt(errLpGrad)

         errLpMin     = min(errLpMin,     errLp)
         errLpGradMin = min(errLpGradMin, errLpGrad)

         write(*,'(2(2(A,E9.3),/))') 
     &      'FEM: Lp norm of gradient = ', errLpGrad, 
     &                 ',  minimum in loop = ', errLpGradMin,
     &      '     Lp norm of error    = ', errLp,     
     &                 ',  minimum in loop = ', errLpMin


c === generate metric (from U) optimal for the L_p norm of gradient
         Call Func2GradMetricMAX(Func, nv, nt, vrt, tri,
     &                           Metric, MaxWr, rW)

c        Call Nodal2MetricVAR(U, vrt, nv, tri, nt, bnd, nb,
c    &                        Metric, MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_gradnorm(nv, Lp, Metric)


c === adapt the mesh using metric Metric
         control(1) = 400     ! MaxSkipE
         control(2) = 150000  ! MaxQItr
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
     &                        'Unspecified error in mbaNodal')


c === make PostScript figure of the adaptive mesh. The name must terminate with .ps
c        Call graph(nv,vrt, nt,tri, 'mesh_final.ps')
         Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps',
     &       'Mesh minimizing gradient of interpolation error')
      End do


c === testing the results
      If(Quality.LT.0.4) Stop 911
      Stop 

      End



C ======================================================================
      Real*8  Function Func(xy)
C ======================================================================
      implicit none
      include  'assemble.fd'

      Real*8   dDATA(1), Coef(4, 4), xy(2), x, y
      Integer  iDATA(1), DexactU, iSYS(MAXiSYS), i, label
      EXTERNAL DexactU

C ======================================================================
      x = xy(1)
      y = xy(2)

      i = DexactU(x, y, label, dDATA, iDATA, iSYS, Coef)

      Func = Coef(1,1)

      Return
      End


