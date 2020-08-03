C =====================================================================
      Program Main
C =====================================================================
      implicit none
 
c nvmax - maximum number of mesh nodes
c ntmax - maximum number of mesh triangles
c nbmax - maximum number of boundary edges
c namax - maximum number of non-zero matrix entries
      Integer   nvmax,ntmax,nbmax,namax, nvMetaMax,ntMetaMax
      parameter(nvmax =  10 000, ntmax = 2*nvmax, nbmax = 2 000)
      parameter(namax = 900 000)
      parameter(nvMetaMax = 5 * nvmax, ntMetaMax = 5*ntmax)

c work memory
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 3 000 000, MaxWi = 5 000 000)

      Integer  iW(MaxWi)
      Real*8   rW(MaxWr)


C ======================================================================
C Mesh definition (first mesh)
C ======================================================================
      Integer  nv, nvfix, labelV(nvmax), fixedV(nvmax)
      Real*8   vrt(2,nvmax)

      Integer  nb, nbfix, bnd(2,nbmax), labelB(nbmax), fixedB(nbmax)

      Integer  nc, labelC(nbmax)
      Real*8   crv(2,nbmax)

      Integer  nt, ntfix, tri(3,ntmax), labelT(ntmax), fixedT(ntmax)


C ======================================================================
C Mesh definition (second mesh)
C ======================================================================
      Integer  nv2,nvfix2, labelV2(nvmax), fixedV2(nvmax)
      Real*8   vrt2(2,nvmax)

      Integer  nb2,nbfix2, bnd2(2,nbmax), labelB2(nbmax), fixedB2(nbmax)

      Integer  nc2, labelC2(nbmax)
      Real*8   crv2(2,nbmax)

      Integer  nt2,ntfix2, tri2(3,ntmax), labelT2(ntmax), fixedT2(ntmax)


C ======================================================================
C Meta mesh (intermidiate data structure)
C ======================================================================
      Integer  nv12, nt12, tri12(3,ntMetaMax), parents(2,ntMetaMax)
      Real*8   vrt12(2,nvMetaMax)


C ======================================================================
C For library aniFEM
C ======================================================================
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  IA(nvmax), JA(namax), controlFEM(3), nRow, nCol
      Real*8    A(namax), dDATAFEM(1)
      Integer  iDATAFEM(1)

      EXTERNAL FEM2Dext


C ======================================================================
C For library aniILU
C ======================================================================
      EXTERNAL matvec, prevec2
      Integer  imatvec(1), iprevec(1)

      Integer  iter, info, nunit, verb, UsedWr, UsedWi
      Real*8   resid, tau1,tau2, partlur,partlurout 


C ======================================================================
C For library aniPRJ
C ======================================================================
      Real*8  U1(ntmax), U2(ntmax)


C LOCAL VARIABLEs
      Integer  i, iERR, ipBCG
      Real*8   x, y, errL2, s

C ======================================================================
c ... load the first mesh. The extension must be .ani
      Call loadMani(nv, nvfix, nvmax, vrt, labelV, fixedV,
     &              nb, nbfix, nbmax, bnd, labelB, fixedB,
     &              nc,               crv, labelC, 
     &              nt, ntfix, ntmax, tri, labelT, fixedT,
     &              "../data/simple2.ani")

c ... load the second mesh. The extension must be .ani
      Call loadMani(nv2, nvfix2, nvmax, vrt2, labelV2, fixedV2,
     &              nb2, nbfix2, nbmax, bnd2, labelB2, fixedB2,
     &              nc2,                crv2, labelC2, 
     &              nt2, ntfix2, ntmax, tri2, labelT2, fixedT2,
     &              "../data/simple3.ani")

c === draw the initial mesh (no labels of mesh points is available)
c The name must have extension with .ps
      Call graph_demo(nv,vrt, nt,tri, 'mesh_first.ps', '1st mesh')

      Call graph_demo(nv2,vrt2, nt2,tri2, 'mesh_second.ps', '2nd mesh')


c === no extra data is provided for the user subroutine Ddiff
      dDATAFEM(1) = 0D0
      iDATAFEM(1) = 0

c general sparse matrix in the AMG format 
      controlFEM(1) = IOR(MATRIX_GENERAL, FORMAT_AMG)
      controlFEM(2) = 1

      Call BilinearFormTemplate(
     &     nv, nb, nt, vrt, labelV, bnd, labelB, tri, labelT,
     &     FEM2Dext, dDATAFEM, iDATAFEM, controlFEM,
     &     nvmax, namax, IA, JA, A, U2, nRow, nCol,
     &     MaxWi, MaxWr, iW, rW)


c === library aniPRJ: create a meta-mesh
      Call MetaMesh(nv, vrt, nt, tri, nv2, vrt2, nt2, tri2,
     &              nv12, nvMetaMax, vrt12, 
     &              nt12, ntMetaMax, tri12, parents,
     &              MaxWi, MaxWr, iW, rW, iERR)

      Call graph_demo(nv12,vrt12, nt12,tri12, 
     &                'mesh_final.ps', 'Intersection of two meshes')

c assemble the right-hand side
      Do i = 1, nt2
         U2(i) = 1d0
      End do

      Call assemble_rhs(nv, vrt, nt, tri, nv2, vrt2, nt2, tri2,
     &                  nv12, vrt12, nt12, tri12, parents, 
     &                  IDEN, FEM_P0, IDEN, FEM_P0,
     &                  U1, U2, MaxWi, iW)

      write(*,*)'matrix and rhs: strange rhs?'
      do i = 1, nRow
         write(*,*) A(i), U1(i)
      end do

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

      If(iERR.NE.0) Then
         Write(*,*) 'Initialization(1) of iluoo failed, iERR=', iERR
         Stop
      End if

c iterative solution
      If(UsedWr + 8*nRow.GT.MaxWr) Then
         Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr + 8*nRow
         Stop
      End if

      ipBCG = UsedWr + 1

      iter  = 10000  ! max number of iterations
      resid = 1d-12  ! final residual
      info  = 0      ! no troubles on input
      nunit = 6      ! output to display 
      Do i = 1, nRow ! initial guess
         U2(i) = i !0D0
      End do

      iprevec(1) = nRow
      imatvec(1) = nRow
      Call slpbcgs(prevec2, iprevec, iW, rW,
     &             matvec,  imatvec, IA, JA, A,
     &             rW(ipBCG), nRow, 8,
     &             nRow, U1, U2,
     &             iter, resid, info, nunit)

      If(info.NE.0) Stop 'BiCGStab had failed'
 

c === draw isolines
      errL2 = 0D0
      Do i = 1, nt
         errL2 = errL2 + (1d0 - U2(i)) ** 2
      End do
      Write(*,*) 'Error =', dsqrt(errL2) 

c     Call isolines(U2, nv,vrt, nt,tri, nb,bnd,'isolines.ps',20)
c     Call isolines_demo(U2, nv,vrt, nt,tri, nb,bnd,
c    &          'isolines.ps', 10, 'Solution isolines')


c === scale the meshes back
      Call scale2Meshes(nv, vrt, nv2, vrt2, .FALSE.)

      Stop
      End



C ======================================================================
C  The user defined routines required above
C ======================================================================
c Templated routine for the elemental matrix. It calls standard bilinear
c forms and imposes the boundary conditions using the provided labels.
C ======================================================================
      Subroutine FEM2Dext(XY1, XY2, XY3, 
     &                    lbE, lbF, lbP, dDATA, iDATA, iSYS,
     &                    LDA, A, F, nRow, nCol,
     &                    templateR, templateC)
C ======================================================================
      Include 'fem2Dtri.fd'
      Include 'assemble.fd'

      Real*8  XY1(*), XY2(*), XY3(*)
      
      Integer lbE, lbF(3), lbP(3)
      Real*8  dDATA(*)
      Integer iDATA(*), iSYS(*)

      Real*8  A(LDA, *), F(*)
      Integer templateR(*), templateC(*)

C Local variables
      Integer  ANI_Dnull
      External ANI_Dnull

C ======================================================================
      nRow = 1
      nCol = 1

c ... set up templated degrees of freedom for rows and columns. 
      Do i = 1, 1
         templateR(i) = Edof
      End do

      Do i = 1, 1
         templateC(i) = templateR(i)
      End do

c ... compute the stiffness matrix (M)
      label = lbE

c     A(1:1,1:1) is elemental identity operator;
c     in other words, for the bilinear form <iden(P0), iden(P0)>
      Call fem2Dtri(XY1, XY2, XY3,
     &              IDEN, FEM_P0, IDEN, FEM_P0,
     &              label, ANI_Dnull, dDATA, iDATA, iSYS, 2,
     &              LDA, A, ir, ic)

      Return
      End



