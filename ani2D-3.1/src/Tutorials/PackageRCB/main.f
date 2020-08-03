c ======================================================================
      Program  mainRCB
c ======================================================================
c The program demonstrates the use of local refinement/coarsening tool 
c librcb2D-2.x.a based on the Marked Edge Bisection. It starts from a 
c simple coarse triangulation, refines it nlevel times according to the 
c rule given in RefineRule, and then coarses the refined mesh nlevel 
c times according to the rule given in CoarseRule. Example of RefineRule
c and CoarseRule are at the file end. Initial bisection of each element 
c is defined in Initialize.
c
c The program uses routine graph() from libview2D-2.x.a. Therefore, this
c visualization library has to be installed prior linking:
c
c   cd src/aniVIEW; make lib
c
c On output 3 PostScript figures of the initial mesh, refined mesh and 
c coarsed mesh are saved in files.
c
c Caution! Even if the coarsed mesh looks like the initial, orientation 
c of triangles may be different!
c ======================================================================
      implicit none

c ... user defined procedures 
      external  RefineRule,  CoarseRule 

c ... nvmax - maximum number of mesh nodes
c ... ntmax - maximum number of mesh triangles
c ... nbmax - maximum number of boundary edges
      Integer nvmax,ntmax,nbmax
      Parameter(nvmax = 100 000, ntmax = 2*nvmax, nbmax = 50 000)

c ... standard mesh arrays (see doc/user_guide.pdf for more detail)
c ... number of points, triangles and boundary edges
      Integer  nv, nt, nb

c ... coordinates of mesh points 
      Double precision   vrt(2,nvmax)

c ... connectivity table for triangles and triangle labels
      Integer  tri(3,ntmax), labelT(ntmax)

c ... connectivity table for boundary edges, and edge labels
      Integer  bnd(2,nbmax), labelB(nbmax)

c ... maxlevel - maximum number of levels for refinement and coarsening
      Integer maxlevel
      Parameter(maxlevel = 150)

c ... work memory (at least 11*ntmax+7)
      Integer    MaxWi
      Parameter(MaxWi = 5 500 000)
      Integer  iW(MaxWi)

c ... history of bisection
      Logical history(maxlevel*ntmax)

c ... number of levels of refinement/coarsening
      integer  nlevel
      
c ... local variables
      Integer  ilevel, iERR

c ======================================================================
c Step 1: generate initial mesh

c ... Define a simple mesh for [0,1]^2 consisting of 2 triangles
      nv = 4
      vrt(1,1) = 0d0
      vrt(2,1) = 0d0
      vrt(1,2) = 0d0
      vrt(2,2) = 2d-1
      vrt(1,3) = 1d0
      vrt(2,3) = 1d0
      vrt(1,4) = 1d0
      vrt(2,4) = 0d0

      nt = 2
      tri(1,1) = 1
      tri(2,1) = 2
      tri(3,1) = 3
      labelT(1) = 1 ! x-y<0 - material 1

      tri(1,2) = 3
      tri(2,2) = 4
      tri(3,2) = 1
      labelT(2) = 2 ! x-y>0 - material 2

      nb = 4
      bnd(1,1) = 1
      bnd(2,1) = 2
      labelB(1) = 1    ! Dirichlet face - label 1

      bnd(1,2) = 2
      bnd(2,2) = 3
      labelB(2) = 1    ! Dirichlet face - label 1

      bnd(1,3) = 3
      bnd(2,3) = 4
      labelB(3) = 2    ! Neumann face - label 2

      bnd(1,4) = 4
      bnd(2,4) = 1
      labelB(4) = 1    ! Dirichlet face - label 1

      Write(*,'(A,2I7)')
     &        'Initial mesh:   numbers of nodes and triangles:',nv, nt

c ... draw mesh; demo has been activated
c     call graph(nv,vrt,nt,tri,'mesh_initial.ps')
      call graph_demo(nv,vrt, nt,tri, 'mesh_initial.ps','Initial mesh')
      
c Step 2: initialize data structure (work memory is at least 11*ntmax+7)

      iERR = 0
      call InitializeRCB(nt, ntmax, vrt, tri, MaxWi, iW, iERR)

      If(iERR.GT.0) Call errMesRCB(iERR, 'main', 
     &                             'error in InitializeRCB')


c Step 3: refine initial mesh nlevel times by local bisection. 
c         The rule is defined in RefineRule

      nlevel = 5

      Do ilevel = 1, nlevel
         Call LocalRefine(
     &        nv, nvmax, nb, nbmax, nt, ntmax,
     &        vrt, tri, bnd, labelB, labelT,
     &        RefineRule, ilevel, maxlevel, history,
     &        MaxWi, iW, iERR)

         If(iERR.GT.0) Call errMesRCB(iERR, 'main',
     &                      'Error in LocalRefine')
      End do
c
      Write(*,'(A,2I7)')
     &        'Refined mesh:   numbers of nodes and triangles:',nv, nt

c ... draw mesh; demo graphics has been activated
c     call graph(nv,vrt,nt,tri,'mesh_final.ps')
      call graph_demo(nv,vrt, nt,tri, 'mesh_final.ps', 
     &               'Locally refined mesh')
      
c Step 4: coarse the refined mesh nlevel times by local coarsening. 
c         The rule is defined in CoarseRule

      Do ilevel = nlevel, 1, -1
         Call LocalCoarse(
     &        nv, nvmax, nb, nbmax, nt, ntmax,
     &        vrt, tri, bnd, labelB,
     &        CoarseRule, ilevel, maxlevel, history,
     &        MaxWi, iW, iERR)

         If(iERR.GT.0) Call errMesRCB(iERR, 'main',
     &                     'Error in LocalCoarse')
      End do
      
      Write(*,'(A,2I7)')
     &        'Coarsened mesh: numbers of nodes and triangles:',nv, nt

c ... draw mesh; demo graphics has been activated
c     call graph(nv,vrt,nt,tri,'mesh_coarse.ps')
      call graph_demo(nv,vrt, nt,tri, 'mesh_coarse.ps','Coarsened mesh')

      Stop 
      End



C ======================================================================
c User routine  defines the rule for local refinement depending on 
c current level.
c On input: nE  current number of elements
c           IPE current connectivity table
c           XYP current coordinates
c           ilevel current level of refinement
c On output: verf marker for refinement of each element
c            (0 - no need to refine, 1 - refine by single bisection,
c             2 - refine by two levels of bisection preserving the shape)
C ======================================================================
      Subroutine RefineRule (nE, IPE, XYP, verf, ilevel)
C ======================================================================
      implicit none
      
      Integer             nE
      Integer             IPE(3,*)
      Double precision    XYP(2,*)
      Integer             verf(*)
      Integer             ilevel

      Integer             i
      Double precision    xy, xy1, xy2, xy3

C ================================================================
      If (ilevel .le. 0) then
        Do i = 1, nE
          verf(i) =  2 ! two levels of bisection (keep the shape)
        End do
      Else  ! refine towards the diagonal y=x
        Do i = 1, nE
          xy1 = XYP(2,IPE(1,i)) - XYP(1,IPE(1,i))
          xy2 = XYP(2,IPE(2,i)) - XYP(1,IPE(2,i))
          xy3 = XYP(2,IPE(3,i)) - XYP(1,IPE(3,i))
          xy  = (xy1 **2 + xy2 **2) *
     &          (xy1 **2 + xy3 **2) *
     &          (xy2 **2 + xy3 **2)
       
          If (xy .eq. 0) then 
            verf(i) =  2 ! two levels of bisection (keep the shape)
           else
            verf(i) =  0 ! no need to refine
          End if
        End do
       End if 

       Return
       End



C ======================================================================
c User routine  defines the rule for local coarsening depending on 
c current level.
c On input: nE  current number of elements
c           IPE current connectivity table
c           XYP current coordinates
c           ilevel current level of refinement
c On output: verf marker for coarsening of each element
c            (0 - no need to coarse, 1 - coarse by single merging,  
c             2 - coarse by two levels of merging   preserving the shape)
C ======================================================================
      Subroutine CoarseRule (nE, IPE, XYP, verf, ilevel)
C ======================================================================
      implicit none
      
      Integer                nE
      Integer                IPE(3,*)
      Double precision       XYP(2,*)
      Integer                verf(*)
      Integer                ilevel

      Integer                i
      Double precision       xy, xy1, xy2, xy3

C ======================================================================
      If (ilevel .le. 0) then
        Do i = 1, nE
          verf(i) =  2 ! two levels of merging (keep the shape)
        End do
      Else ! coarse towards the diagonal y=x
        Do i = 1, nE
          xy1 = XYP(2,IPE(1,i)) - XYP(1,IPE(1,i))
          xy2 = XYP(2,IPE(2,i)) - XYP(1,IPE(2,i))
          xy3 = XYP(2,IPE(3,i)) - XYP(1,IPE(3,i))
          xy  = (xy1 **2 + xy2 **2) *
     &          (xy1 **2 + xy3 **2) *
     &          (xy2 **2 + xy3 **2)
       
          If (xy .eq. 0) then 
            verf(i) =  2 ! two levels of merging (keep the shape)
           else
            verf(i) =  0 ! no need to coarse
          End if
        End do
       End if 

       Return
       End



