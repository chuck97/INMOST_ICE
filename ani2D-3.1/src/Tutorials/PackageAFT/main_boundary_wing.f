C ======================================================================
C Example of using mesh generator with an analytic geometry description
C ======================================================================
      program main  

      implicit none
      integer nvmax,ntmax,nbmax
c     nvmax   - maximum number of mesh nodes
c     ntmax   - maximum number of mesh triangles
c     nbmax   - maximum number of boundary edges
      parameter(nvmax=150000,ntmax=2*nvmax,nbmax=10000)

c mesh generator data specifying domain analytically

c rectangle
c     double precision bv(2,4),bltail(2,4) 
c     integer          Nbv,Nbl,bl(7,4)
c     data             Nbv/4/,Nbl/4/
c     data             bv/0.0,0.25, 0.0,0.5, 1.0,0.5, 1.0,0.25/  ! boundary nodes
c     data             bl/1,2,0,0,1,1,0, 3,4,0,0,2,1,0,          ! outer boundary edges
c    &                    2,3,0,0,3,1,0, 4,1,0,0,4,1,0/          ! outer boundary edges
c     data             bltail/0,0, 0,0, 0,0, 0,0/                ! curved data for each outer boundary edge

c complement of a wing NACA0012 to the unit square
      double precision bv(2,7),bltail(2,8)
      integer          Nbv,Nbl,bl(7,8)
      data             Nbv/7/,Nbl/8/
      data             bv/0,0, 0,1, 1,1, 1,0, .4,.5, .6,.5, 1,.5/       ! boundary nodes
      data             bl/1,2,0,-1,-1,1,0, 4,1,0,-1,-1,1,0,             ! outer boundary edges
     &                    2,3,0,-1,1,1,0,  7,4,0,-1,1,1,0,              ! outer boundary edges
     &                    3,7,0,-1,1,1,0,                               ! outer boundary edges
     &                    6,7,2,0,11,1,1,                               ! slit  boundary edges
     &                    6,5,1,-1,2,1,0,  5,6,1,-1,2,1,0/              ! wing  boundary edges
      data             bltail/0,0, 0,0, 0,0, 0,0, 0,0, 0,1, 0,.5, .5,1/ ! curved data for each outer boundary edge

      integer nv,nt,nb,nc
      double precision crv(2,nbmax), vrt(2,nvmax)
      integer          iFNC(nbmax), labelT(ntmax),
     &                 tri(3,ntmax), bnd(2,nbmax), labelB(nbmax)
      double precision h
      integer ierr

c AFT2D library function
      Integer  aft2dboundary
      EXTERNAL aft2dboundary

c the name of the user written function (see file crv_model.c)
      EXTERNAL userboundary
C ======================================================================
      call registeruserfn(userboundary)  ! register the name in the library


      h = 0.02  ! mesh step of the quasi-uniform mesh
C Generate quasiuniform mesh with meshstep h
      ierr = aft2dboundary( Nbv, bv, Nbl, bl, bltail, h,  ! geometric data
     &                      nv, vrt,                      ! mesh data on output
     &                      nt, tri, labelT,
     &                      nb, bnd, labelB,
     &                      nc, crv, iFNC)
      If(ierr.ne.0) stop ' error in function aft2dboundary'
      write(*,*) 'mesh: number of triangles/vertices ',nt,nv 

c  The name must terminate with .ps
c  Demo graphics has been activated
c     call graph(nv,vrt, nt,tri, 'mesh_final.ps')
      call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps', 
     &               'Mesh built from boundary given analytically')

      stop  
      end

