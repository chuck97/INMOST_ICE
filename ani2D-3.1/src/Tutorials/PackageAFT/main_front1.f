      program main ! example of calls of mesh generator with segment data
      implicit none
      integer nvmax,ntmax,nbmax
c     nvmax   - maximum number of mesh nodes
c     ntmax   - maximum number of mesh triangles
c     nbmax   - maximum number of boundary edges
      parameter(nvmax=150000,ntmax=2*nvmax,nbmax=10000)

c mesh generator data specifying domain via in the segment format
      double precision vbr(2,nbmax)
      integer          Nbr,Nvr,brd(2,nbmax)


      integer          nv,nt,nb
      double precision vrt(2,nvmax)
      integer          labelB(nbmax),labelT(ntmax)
      integer          tri(3,ntmax),bnd(2,nbmax)

c ... AFT2D library function
      Integer  aft2dfront   
      EXTERNAL aft2dfront   

      Integer  i,j,k, ierr

C Read input file that contains coordinates of boundary points
      Open(10, file='../data/front/front_lake.txt') 
        Read(10, *) Nvr, Nbr               
        Do i = 1, Nvr
            Read(10, *) (vbr(j,i), j=1,2)
        End do
        Do i = 1, Nbr
           Read(10, *) (brd(3-j,i), j=1,2)
        End do
      Close(1)

c draw initial front
      nv = 0
      nt = 0
      Call graph_front(Nvr, vbr, Nbr, brd, 'mesh_initial.ps')

C Generate a mesh  starting  from boundary mesh
      ierr=aft2dfront(
     &           Nbr, brd, Nbr, vbr,                         ! segment data
     &           nv, vrt,                                    ! mesh data on output
     &           nt, tri, labelT,
     &           nb, bnd, labelB)
      If (ierr.ne.0) stop ' error in function aft2dfront'   
      Write(*,5000) Nbr, nt, nv

c demo graphics has been activated
c     Call graph(nv,vrt, nt,tri, 'mesh_final.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps',
     &               'Mesh built from the discrete boundary')

 5000 Format(
     & 'Program generates a mesh using the advanced front technique.',/,
     & '  The initial front has', I5, ' edges.',/,
     & '  The final mesh has', I5, ' triangles and', I5, ' vertices.',/,
     & 'Create Postscript figure $(ANIHOME)/bin/mesh_final.ps.',/)

      Stop
      End



