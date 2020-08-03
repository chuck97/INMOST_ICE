C ======================================================================
      Subroutine DGjump_AJ(XYZ, iE, IEE, ks, 
     &                     FEMtypeA, FEMtypeB, 
     &                     Dedge, dDATA, iDATA, iSYS, order,
     &                     LDA, A, nRow, nCol)
C ======================================================================
      Implicit none
      Include 'fem2Dtri.fd'
C ======================================================================
      Real*8  XYZ(2, 3, 2)
      Integer iE, IEE(3), ks(3), FEMtypeA, FEMtypeB(3)

      Integer iDATA(*), iSYS(*), LDA, nRow, nCol, order
      Real*8  dDATA(*), A(LDA, *)

      Integer  Dedge
      EXTERNAL Dedge

C Local variables
      Integer  i,j,n, i1,i2,i3, ir,ic, label, iEt, idx(10), idy(10)
      Real*8   C(10, 10), edge_length, d, s

      Integer  iref(4)
      DATA     iref/1,2,3,1/

C ======================================================================
c ... initialize array by zeros
      Do i = 1, nRow
         Do j = 1, nCol
            A(i, j) = 0D0
         End do
      End do

c ... compute DG jump terms   s [[ u ]]  [[ v ]]  on edge (i1,i2) 
      Do i1 = 1, 3
         iEt = IEE(i1)

         If(iE.LT.iEt) Then
            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            Call enumDOF(FEMtypeA, n, idx, i1)
            Call enumDOF(FEMtypeB(i1), n, idy, 1)
            Do i = 1, n
               idy(i) = idy(i) + ks(i1)
            End do

            d = edge_length(XYZ(1,i1,1), XYZ(1,i2,1))
            s = 1D0 / d

c           + u1 v1
            Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,1),
     &                     IDEN, FEMtypeA, IDEN, FEMtypeA,
     &                     label, Dedge, dDATA, iDATA, iSYS, order,
     &                     10, C, ir, ic)

            Call addBlock(idx, idx, ir, ic, LDA, A, 10, C, s)

c           - u1 v2 - u2 v1
            If(FEMtypeA.NE.FEMtypeB(i1)) Then
               Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,1),
     &                        IDEN, FEMtypeA, IDEN, FEMtypeB(i1),
     &                        label, Dedge, dDATA, iDATA, iSYS, order,
     &                        10, C, ir, ic)
            End if

            Call addBlock(     idx, idy, ir, ic, LDA, A, 10, C, -s)
            Call addBlockTrans(idx, idy, ir, ic, LDA, A, 10, C, -s)

c           + u2 v2
            If(FEMtypeA.NE.FEMtypeB(i1)) Then
               Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,1),
     &                        IDEN, FEMtypeB(i1), IDEN, FEMtypeB(i1),
     &                        label, Dedge, dDATA, iDATA, iSYS, order,
     &                        10, C, ir, ic)
            End if

            Call addBlock(idy, idy, ir, ic, LDA, A, 10, C, s)
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine DGjump_SIP(XYZ, iE, IEE, ks, 
     &                      FEMtypeA, FEMtypeB, 
     &                      Ddiff, dDATA, iDATA, iSYS, order,
     &                      LDA, A, nRow, nCol)
C ======================================================================
      Implicit none
      Include 'fem2Dtri.fd'
C ======================================================================
      Real*8  XYZ(2, 3, 2)
      Integer iE, IEE(3), ks(3), FEMtypeA, FEMtypeB(3)

      Integer iDATA(*), iSYS(*), LDA, nRow, nCol, order
      Real*8  dDATA(*), A(LDA, *)

      Integer  Ddiff
      EXTERNAL Ddiff

C Local variables
      Integer  i,j,n, i1,i2,i3, ir,ic, label, iEt, idx(10), idy(10)
      Real*8   B(10, 10), s

      Integer  iref(4)
      DATA     iref/1,2,3,1/

C ======================================================================
c ... initialize array by zeros
      Do i = 1, nRow
         Do j = 1, nCol
            A(i, j) = 0D0
         End do
      End do


c ... jump terms  -{{ DUDN(v) }} [[ u ]] - [[ v ]] {{ DUDN(u) }}  on edge (i1,i2) 
      Do i1 = 1, 3
         iEt = IEE(i1)

         If(iE.LT.iEt) Then
            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            Call enumDOF(FEMtypeA, n, idx, i1)
            Call enumDOF(FEMtypeB(i1), n, idy, 1)
            Do i = 1, n
               idy(i) = idy(i) + ks(i1)
            End do

            s = 0.5D0

c           - DUDN(u1) v1 - DUDN(v1) u1
            Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,1),
     &                     DUDN, FEMtypeA, IDEN, FEMtypeA,
     &                     label, Ddiff, dDATA, iDATA, iSYS, order,
     &                     10, B, ir, ic)

            Call addBlock(     idx, idx, ir, ic, LDA, A, 10, B, -s)
            Call addBlockTrans(idx, idx, ir, ic, LDA, A, 10, B, -s)

c           DUDN(u1) v2 + DUDN(v1) u2
            If(FEMtypeA.NE.FEMtypeB(i1)) Then
               Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,1),
     &                        DUDN, FEMtypeA, IDEN, FEMtypeB(i1),
     &                        label, Ddiff, dDATA, iDATA, iSYS, order,
     &                        10, B, ir, ic)
            End if

            Call addBlock(     idx, idy, ir, ic, LDA, A, 10, B, s)
            Call addBlockTrans(idx, idy, ir, ic, LDA, A, 10, B, s)

c           DUDN(u2) v1 + DUDN(v2) u1
            Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,2),
     &                     DUDN, FEMtypeB(i1), IDEN, FEMtypeA,
     &                     label, Ddiff, dDATA, iDATA, iSYS, order,
     &                     10, B, ir, ic)

            Call addBlock(     idy, idx, ir, ic, LDA, A, 10, B, s)
            Call addBlockTrans(idy, idx, ir, ic, LDA, A, 10, B, s)

c           - DUDN(u2) v2 - DUDN(v2) u2
            If(FEMtypeA.NE.FEMtypeB(i1)) Then
               Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,2),
     &                        DUDN, FEMtypeB(i1), IDEN, FEMtypeB(i1),
     &                        label, Ddiff, dDATA, iDATA, iSYS, order,
     &                        10, B, ir, ic)
            End if

            Call addBlock(     idy, idy, ir, ic, LDA, A, 10, B, -s)
            Call addBlockTrans(idy, idy, ir, ic, LDA, A, 10, B, -s)
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine DGjump_NIP(XYZ, iE, IEE, ks, 
     &                      FEMtypeA, FEMtypeB, 
     &                      Ddiff, dDATA, iDATA, iSYS, order,
     &                      LDA, A, nRow, nCol)
C ======================================================================
      Implicit none
      Include 'fem2Dtri.fd'
C ======================================================================
      Real*8  XYZ(2, 3, 2)
      Integer iE, IEE(3), ks(3), FEMtypeA, FEMtypeB(3)

      Integer iDATA(*), iSYS(*), LDA, nRow, nCol, order
      Real*8  dDATA(*), A(LDA, *)

      Integer  Ddiff
      EXTERNAL Ddiff

C Local variables
      Integer  i,j,n, i1,i2,i3, ir,ic, label, iEt, idx(10), idy(10)
      Real*8   B(10, 10), s

      Integer  iref(4)
      DATA     iref/1,2,3,1/

C ======================================================================
c ... initialize array by zeros
      Do i = 1, nRow
         Do j = 1, nCol
            A(i, j) = 0D0
         End do
      End do


c ... jump terms  -{{ DUDN(v) }} [[ u ]] - [[ v ]] {{ DUDN(u) }}  on edge (i1,i2) 
      Do i1 = 1, 3
         iEt = IEE(i1)

         If(iE.LT.iEt) Then
            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            Call enumDOF(FEMtypeA, n, idx, i1)
            Call enumDOF(FEMtypeB(i1), n, idy, 1)
            Do i = 1, n
               idy(i) = idy(i) + ks(i1)
            End do

            s = 0.5D0

c           - DUDN(u1) v1 + DUDN(v1) u1
            Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,1),
     &                     DUDN, FEMtypeA, IDEN, FEMtypeA,
     &                     label, Ddiff, dDATA, iDATA, iSYS, order,
     &                     10, B, ir, ic)

            Call addBlock(     idx, idx, ir, ic, LDA, A, 10, B,  s)
            Call addBlockTrans(idx, idx, ir, ic, LDA, A, 10, B, -s)

c           DUDN(u1) v2 - DUDN(v1) u2
            If(FEMtypeA.NE.FEMtypeB(i1)) Then
               Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,1),
     &                        DUDN, FEMtypeA, IDEN, FEMtypeB(i1),
     &                        label, Ddiff, dDATA, iDATA, iSYS, order,
     &                        10, B, ir, ic)
            End if

            Call addBlock(     idx, idy, ir, ic, LDA, A, 10, B, -s)
            Call addBlockTrans(idx, idy, ir, ic, LDA, A, 10, B,  s)

c           DUDN(u2) v1 - DUDN(v2) u1
            Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,2),
     &                     DUDN, FEMtypeB(i1), IDEN, FEMtypeA,
     &                     label, Ddiff, dDATA, iDATA, iSYS, order,
     &                     10, B, ir, ic)

            Call addBlock(     idy, idx, ir, ic, LDA, A, 10, B, -s)
            Call addBlockTrans(idy, idx, ir, ic, LDA, A, 10, B,  s)

c           - DUDN(u2) v2 + DUDN(v2) u2
            If(FEMtypeA.NE.FEMtypeB(i1)) Then
               Call fem2Dedge(XYZ(1,i1,1), XYZ(1,i2,1), XYZ(1,i3,2),
     &                        DUDN, FEMtypeB(i1), IDEN, FEMtypeB(i1),
     &                        label, Ddiff, dDATA, iDATA, iSYS, order,
     &                        10, B, ir, ic)
            End if

            Call addBlock(     idy, idy, ir, ic, LDA, A, 10, B,  s)
            Call addBlockTrans(idy, idy, ir, ic, LDA, A, 10, B, -s)
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine DG_init(XY1, XY2, XY3, iSYS, XYZ, iE, IEE) 
C ======================================================================
      implicit none
      Real*8  XY1(2,2), XY2(2,2), XY3(2,2), XYZ(2,3,2)
      Integer iSYS(*), iE, IEE(3), i, n
C ======================================================================
      Do n = 1, 2
         Do i = 1, 2
            XYZ(i, 1, n) = XY1(i, n)
            XYZ(i, 2, n) = XY2(i, n)
            XYZ(i, 3, n) = XY3(i, n)
         End do
      End do

      iE = iSYS(3)

      IEE(1) = iSYS(13)
      IEE(2) = iSYS(14)
      IEE(3) = iSYS(15)
        
      Return
      End



C ======================================================================
      Subroutine DG_template(iE, IEE, FEMtypeA, FEMtypeB, 
     &                       nRow, templateR, ks) 
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer iE, IEE(3), FEMtypeA, FEMtypeB(3)
      Integer nRow, templateR(*), ks(3), iEt, i, k, n

      Integer dof(3)
      DATA    dof/DGdof12, DGdof23, DGdof31/
C ======================================================================
c     add degrees of freedom for central triangle
      Call listDOF(FEMtypeA, n, templateR)

      Do i = 1, n
         Call addXnode(templateR(i), DGdof)
      End do
      nRow = n

c     add degrees of freedom for neighboring triangles
      Do k = 1, 3
         iEt = IEE(k)

         If(iEt.GT.0) Then
            ks(k) = nRow
            Call listDOF(FEMtypeB(k), n, templateR(nRow+1))

            Do i = nRow + 1, nRow + n
               Call addXnode(templateR(i), dof(k))
            End do
            nRow = nRow + n
         End if
      End do
        
      Return
      End



C ======================================================================
      Subroutine addBlock(idx, idy, ir, ic, LDA, A, LDC, C, factor)
C ======================================================================
      implicit none
      Integer  idx(*), idy(*), ir, ic, LDA, LDC
      Real*8   A(LDA, *), C(LDC, *), factor

      Integer  i, j, i1, j1

      Do i = 1, ir
         i1 = idx(i)
         Do j = 1, ic
            j1 = idy(j)
            A(i1, j1) = A(i1, j1) + C(i, j) * factor
         End do
      End do
    
      Return
      End


 
C ======================================================================
      Subroutine addBlockTrans(idx, idy, ir, ic, LDA, A, LDC, C, factor)
C ======================================================================
      implicit none
      Integer  idx(*), idy(*), ir, ic, LDA, LDC
      Real*8   A(LDA, *), C(LDC, *), factor

      Integer  i, j, i1, j1

      Do i = 1, ir
         i1 = idx(i)
         Do j = 1, ic
            j1 = idy(j)
            A(j1, i1) = A(j1, i1) + C(i, j) * factor
         End do
      End do
    
      Return
      End



C ======================================================================
      Subroutine DG_nodal(nP, nE, IPE, iDATA, Udg, Up)
C ======================================================================
C Routines converts unstructured DG data to vertex data. A simple
C algorithm has been implemented.
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer nP, nE, IPE(3,*), iDATA(*)
      Real*8  Udg(*), Up(*)

      Integer i, k, n, m, iP, dof(MaxDOFs) 
C ======================================================================
      m = 0
      Do n = 1, nE
         Call listDOF(iDATA(n), k, dof)

         Do i = 1, 3
            iP = IPE(i, n)
            Up(iP) = Udg(m + i)
         End do

         m = m + k
      End do

      Return 
      End

