c ======================================================================
c The program illustrates organization of an adaptive loop for building 
c a mesh minimizing the interpolation error for a given function DexactU. 
c The function in [0,1]^2 is taken from D'Azevedo's paper on optimal 
c triangulations. The second example (commented out) is a snake function 
c in [-1,1]^2. 
c ======================================================================
      program Interpolation
c ======================================================================
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
      Integer  nv, nvfix, labelV(nvmax), fixedV(1)
      Real*8   vrt(2,nvmax)

      Integer  nb, nbfix, bnd(2,nbmax), labelB(nbmax), fixedB(1)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2,nbmax)

      Integer  nt, ntfix, tri(3,ntmax), labelT(ntmax), fixedT(1)

      DATA     nvfix/0/, fixedV/0/,  nbfix/0/, ntfix/0/


c ======================================================================
c for library aniAFT (advanced front technique)
c ======================================================================
c Specify initial front with 4 segments
      Integer   Nbv, Nbl, bl(7,4)
      DATA      Nbv/4/, Nbl/4/

      Real*8    bv(2,4),  bltail(2,4)
      DATA      bv/-1,-1, -1,1, 1,1, 1,-1/            ! boundary nodes for snake function
c     DATA      bv/ 0, 0,  0,1, 1,1, 1, 0/            ! boundary nodes for D'Azevedo function

      DATA      bl/1,2,0,-1,1,1,0,  2,3,0,-1,2,1,0,   ! outer boundary edges
     &             3,4,0,-1,3,1,0,  4,1,0,-1,4,1,0/   ! outer boundary edges
      DATA      bltail/0,0, 0,0, 0,0, 0,0/            ! curved DATA for each outer boundary edge

      Integer   aft2dboundary
      EXTERNAL  aft2dboundary


c ======================================================================
c for library aniFEM (finite element matrices)
c ======================================================================
      include  'fem2Dtri.fd'
      include  'assemble.fd'

      Real*8    U(nvmax), Metric(3,nvmax), dDATA(1)
      Integer   iDATA(1), iSYS(MAXiSYS)

      Integer   ANI_Dnull, DexactU
      EXTERNAL  ANI_Dnull, DexactU


c ======================================================================
c for library aniMBA
c ======================================================================
      Integer  control(6), nEStar
      Real*8   Quality
      EXTERNAL ANI_CrvFunction


c LOCAL VARIABLEs
      Integer   i,n, iLoop, iv1, iv2, iv3, label, iERR, nLOOPs
      Real*8    Lp, h, Uloc(3), trierr, lpglobal, x,y

c ======================================================================
c number of adaptive loops
      nLOOPs = 5


c === generate a quasi-uniform mesh with mesh step h
      h = 0.043
      Write(*,'(A,F7.3)') 'AFT: h = ', h 

      iERR = aft2dboundary(Nbv, bv, Nbl, bl, bltail, h,    ! boundary mesh
     &                     nv, vrt,                        ! output mesh
     &                     nt, tri, labelT,
     &                     nb, bnd, labelB,
     &                     nc, crv, labelC)

      If(iERR.NE.0)   Stop 'Error in function aft2dboundary'
      If(nv.gt.nvmax) Stop 'too many nodes, increase nvmax'
      If(nt.gt.ntmax) Stop 'too many triangles, increase ntmax'
      If(nb.gt.nbmax) Stop 'too many boundary edges, increase nbmax'
          
      Write(*,'(A,3(I6,A))') '     mesh with', nt, ' triangles', 
     &                       nv, ' vertices ', nb, ' boundary edges' 


c === make PostScript figure of the initial mesh.  The name must terminate 
c with .ps. Demo version is activated here.
c     Call graph(nv,vrt, nt,tri, 'mesh_initial.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_initial.ps', 
     &               'Initial quasi-uniform mesh')


c begin adaptive iterative loop
      Do iLoop = 1, nLOOPs
         Write(*,'(/,A,I2)') '===> LOOP: ', iLoop

c define the function at mesh nodes
         label = 1
         Do i = 1, nv
            x = vrt(1, i)
            y = vrt(2, i)

            n = DexactU(x, y, label, dDATA, iDATA, iSYS, U(i))
         End do


c === make PosctScript file with 20 isolines of mesh function U
c The name must terminate with .ps
c Demo vertion is activated here.
         If(iLoop.eq.1) Then
c           Call isolines(U, nv,vrt, nt,tri, nb,bnd, 'iso.ps',20)
            Call isolines_demo(U, nv,vrt, nt,tri, nb,bnd, 'iso.ps',20,
     &                         'Solution isolines on the initial mesh')
         Else
c           Call isolines(U, nv,vrt, nt,tri, nb,bnd, 'his.ps',20)
            Call isolines_demo(U, nv,vrt, nt,tri, nb,bnd, 'his.ps',20,
     &                         'Solution isolines on the final mesh')
         End if


c === generate a metric from nodal values of U that is optimal 
c for the L_p norm of interpolation error
         Lp = 0             ! maximum norm
c        Lp = 2             ! L_2 norm
         Call Nodal2MetricVAR(U,
     &                        vrt, nv, tri, nt, bnd, nb, Metric,
     &                        MaxWr, rW, MaxWi, iW)

         If(Lp.GT.0) Call Lp_norm(nv, Lp, Metric)


c === compute Lp-error
         lpglobal = 0D0
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

            If(Lp.EQ.0) Then
               lpglobal = max(lpglobal, trierr)
            Else
               lpglobal = lpglobal + trierr
            End if
         End do

         If(Lp.GT.0) lpglobal = lpglobal ** (1 / Lp)
         Write(*,'(2(A,E12.4))') 'FEM: Lp =', Lp, 
     &                           '  Lp-error =', lpglobal


c === adapt the mesh using metric Metric
         control(1) = 100     ! MaxSkipE
         control(2) = 10000   ! MaxQItr
         control(3) = 1       ! status = forbid boundary triangles (see aniMBA/status.fd)
         control(4) = 1       ! flagAuto
         control(5) = 1       ! iPrint = minimal level of output information
         control(6) = 0       ! iErrMesgt: only critical termination allowed

         nEStar = 5000
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


c === PostScript figure of the adapted mesh. The name must terminate with .ps
c Demo version is activated.
c        Call graph(nv,vrt, nt,tri, 'mesh_final.ps')
         Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps', 
     &                   'Final mesh minimizing interpolation error')
      End do


c ... testing the results
      If(Quality.LT.0.6) Stop 911
      Stop 

      End


