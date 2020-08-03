C ======================================================================
c The program is the simplest example of using library libaft2D.a for
c generating a quasi-uniform mesh in the unit square. This program uses 
c two libraries libaft2D.a and libview.a
C ======================================================================
      Program main  
C ======================================================================
      implicit none
      integer nvmax,ntmax,nbmax,namax
c nvmax   - maximum number of mesh nodes
c ntmax   - maximum number of mesh triangles
c nbmax   - maximum number of boundary edges
      parameter(nvmax=150000,ntmax=2*nvmax,nbmax=10000)

C ======================================================================
c We define boundary edges for the library libaft2D.a. 
c Note: crv_model.c is not needed here because the computational domain
c is the unit square and has no curved edges. For more details see the 
c PDF document doc/user_guide.pdf.
C ======================================================================

c Define the boundary
      double precision bv(2,4),bltail(2,4)
      integer          Nbv,Nbl,bl(7,4)

c ... Number of boundary nodes and edges
      data             Nbv/4/,Nbl/4/

c ... Coordinates of boundary nodes
      data             bv/0,0, 0,1, 1,1, 1,0/
 
c ... Boundary edges (7 parameters per edge)
      data             bl/1,2,0,-1,1,1,0, 2,3,0,-1,2,1,0,
     &                    3,4,0,-1,3,1,0, 4,1,0,-1,4,1,0/

c ... bltail is dummy here since edges are not curved
      data             bltail/0,0, 0,0, 0,0, 0,0/


c Standard desription of the mesh
c ... Number of vertices, triangles, boundary edges, and curved edges
      integer          nv,nt,nb,nc

c ... Cartesian coordinates of mesh nodes
      double precision vrt(2,nvmax)

c ... connectivity table for triangles, their labels
      integer          tri(3,ntmax), labelT(ntmax)

c ... conenctivity table for edges
      integer          bnd(2,nbmax), labelB(nbmax)

c ... parameterezation of curved edges, ID of functions describing these curves
      double precision crv(2,nbmax)
      integer          iFNC(nbmax)

c ... AFT2D library function
      Integer   aft2dboundary
      EXTERNAL  aft2dboundary


c ... local variables
      double precision h
      integer ierr 

c override quasi-uniform mesh 
      external meshSize               ! function at the bottom of the file
      Call registersizefn( meshSize ) ! register function in the library
c Note: parameter h in this case is dummy


      ierr = aft2dboundary(Nbv, bv, Nbl, bl, bltail, h,  ! geometry
     &                                                   ! mesh data on output
     &                     nv, vrt,                      ! coordinates of nodes
     &                     nt, tri, labelT,              ! triangles and their material
     &                     nb, bnd, labelB,              ! boundary edges and their labels
     &                     nc, crv, iFNC)                ! curved edges and parametrization (dummy)
      If (ierr.ne.0) stop ' error in function aft2dboundary'


      Write(*,*) 'mesh: number of triangles/vertices ', nt, nv 
      If(nv.gt.nvmax) stop ' too many nodes'
      If(nt.gt.ntmax) stop ' too many triangles'
      If(nb.gt.nbmax) stop ' too many boundary edges'
          

c make PostScript figure of the final mesh; the name must terminate with .ps
c demo graphics has been activated
c     Call graph(nv,vrt, nt,tri, 'mesh_final.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps', 
     &               'Mesh built from boundary given analytically')

      Stop 
      End



C ======================================================================
c Example of user defined mesh size function
C ======================================================================
      double precision  Function meshSize(XY)
C ======================================================================
c ... Coordinates of the point
      double precision  XY(2)
c ... local variables
      double precision  R, m0, c0

c distance to point (0,0)
      R = dsqrt(XY(1) ** 2 + XY(2) ** 2)
c minimal mesh size
      m0 = 0.01
c coarsening speed
      c0 = 0.1

      meshSize = dsqrt(m0 ** 2 + (R*c0) ** 2)
      Return
      End

