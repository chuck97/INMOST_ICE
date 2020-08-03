      program main ! example of calls of mesh generator with segment data
      implicit none
      integer nvmax,ntmax,nbmax
c     nvmax   - maximum number of mesh nodes
c     ntmax   - maximum number of mesh triangles
c     nbmax   - maximum number of boundary edges
      parameter(nvmax=150000,ntmax=2*nvmax,nbmax=10000)

c mesh generator data specifying domain via in the segment format
      double precision vbr(2,nbmax)
      integer          Nbr

      integer          nv,nt,nb
      double precision vrt(2,nvmax)
      integer          labelB(nbmax), labelT(ntmax)
      integer          tri(3,ntmax), bnd(2,nbmax)
c ... AFT2D library function
      Integer   aft2dfront   
      EXTERNAL  aft2dfront   


      integer          i,j,dummy,ierr


C Read input file that contains coordinates of boundary points
      Open(1,file='../data/front/front2.txt') 
         Read(1,*) Nbr                   
         Do i = 1, Nbr
            Read(1,*) (vbr(j,i),j=1,2)
         End do
      Close(1)

C Generate a mesh  starting  from boundary mesh
      ierr=aft2dfront(
     &           0, dummy, Nbr, vbr,                         ! segment data
     &           nv, vrt,                                    ! mesh data on output
     &           nt, tri, labelT,
     &           nb, bnd, labelB)
      If (ierr.ne.0) stop ' error in function aft2dfront'   
      Write(*,5000) Nbr, nt, nv

c demo graphics has been activated
c     Call graph(nv,vrt,nt,tri,'mesh_final.ps')
      Call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps',
     &               'Mesh built from the discrete boundary')

 5000 Format(
     & 'Program generates a mesh using the advanced front technique.',/,
     & '  The initial front has', I4, ' edges.',/,
     & '  The final mesh has', I4, ' triangles and', I4, ' vertices.',/,
     & 'Created Postscript figure $(ANIHOME)/bin/mesh_final.ps.',/)

      Stop
      End



