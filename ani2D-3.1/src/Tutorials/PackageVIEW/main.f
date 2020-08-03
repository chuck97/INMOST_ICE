c ======================================================================
      Program  mainView
c ======================================================================
c This program creates a simple mesh, defines a mesh solution,
c and visualize both mesh and isolines of the solution.
c ======================================================================
      implicit none
      integer nvmax,ntmax,nbmax
c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh triangles
c ... nbmax - maximum number of boundary edges
      parameter(nvmax = 1 000, ntmax = 2*nvmax, nbmax = 1 000)


c ... standard mesh arrays (see description in doc/aft_guide.pdf)
      Integer  nv, nt, nb
      Real*8   vrt(2,nvmax), sol(nvmax)
      Integer  tri(3,ntmax), bnd(2,nbmax), labelB(nbmax)

      
C LOCAL VARIABLEs      
      integer i

c ======================================================================
c Define a simple mesh in [0,1]^2 with 2 triangles
      nv = 4
      vrt(1,1) = 0d0
      vrt(2,1) = 0d0
      vrt(1,2) = 0d0
      vrt(2,2) = 1d0
      vrt(1,3) = 1d0
      vrt(2,3) = 1d0
      vrt(1,4) = 1d0
      vrt(2,4) = 0d0

      nt = 2
      tri(1,1) = 1
      tri(2,1) = 2
      tri(3,1) = 3
      tri(1,2) = 3
      tri(2,2) = 4
      tri(3,2) = 1

      nb = 5
      bnd(1,1) = 1
      bnd(2,1) = 2
      labelB(1) = 1    ! outer face - label 1
      bnd(1,2) = 2
      bnd(2,2) = 3
      labelB(2) = 2    ! outer face - label 2
      bnd(1,3) = 3
      bnd(2,3) = 4
      labelB(3) = 3    ! outer face - label 3
      bnd(1,4) = 4
      bnd(2,4) = 1
      labelB(4) = 4    ! outer face - label 4
      bnd(1,5) = 1
      bnd(2,5) = 3
      labelB(5) = 5    ! Interior face - label 5

c Define the linear solution
      Do i = 1, nv
         sol(i) = vrt(1,i) + vrt(2,i)     ! sol=x+y
      End do

c Make a ps-file with the mesh. The name must terminate with .ps
      Call graph(nv,vrt, nt,tri, 'mesh.ps')

c Make a ps-file with 20 isolines of sol (defined at nodes) and boundary 
c edges of the mesh.
c ATTENTION: since there exists an interior boundary edge, it will be 
c shown as well! The name must terminate with .ps
      Call isolines(sol, nv,vrt, nt,tri, nb,bnd,'isolines.ps',20)


c When either mbaMetric() or mbaSolution() from library libmba2D.a is used, 
c one could plot more mesh information. Besides the standard mesh data 
c (nv,nb,nt, vrt,bnd,tri), mbaMetric or mbaSolution fills the first nv 
c entries of the input/output work array iW() with vertex colors.
c These colors may be visualized by calling
c     Call draw(nv, nb, nt, vrt, iW, bnd, labelB, tri, 'mesh.ps')  

      Stop 
      End


