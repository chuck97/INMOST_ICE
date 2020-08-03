C ======================================================================
      Subroutine markDIR(nv, vrt, lbV, nb, bnd, lbB, 
     &                   Dbc, dDATA, iDATA, iSYS)
C ======================================================================
C Mark points as Dirichlet points for edges labled as Dirichlet edges 
C Routine call the user-given routine Dbc
C ======================================================================
      implicit none
      Include 'fem2Dtri.fd'

C ======================================================================
      Integer  nv, lbV(*), nb, bnd(2, *), lbB(*)
      Real*8   vrt(2, *)

      Real*8   dDATA(*)
      Integer  iDATA(*), Dbc, iSYS(*)
      EXTERNAL Dbc

c Local variables
      Integer  i, iv1, iv2, ibc, labelB
      Real*8   x, y, eBC(2)
      Logical  ifXbc

C ======================================================================
      Do i = 1, nv
         lbV(i) = 0
      End do

      Do i = 1, nb
         iv1 = bnd(1, i)
         iv2 = bnd(2, i)

         x = (vrt(1,iv1) + vrt(1,iv2)) / 2
         y = (vrt(2,iv1) + vrt(2,iv2)) / 2

         labelB = lbB(i)
         ibc = Dbc(x, y, labelB, dDATA, iDATA, iSYS, eBC)

         If(ifXbc(ibc, BC_DIRICHLET)) Then
            lbV(iv1) = max(lbV(iv1), lbB(i))
            lbV(iv2) = max(lbV(iv2), lbB(i))
         End if
      End do

      Return
      End


C ======================================================================
      Subroutine applyDIR(LDA, nRow, A, F, k, bc)
C ======================================================================
C Apply the Dirichlet b.c. to the local matrix
C ======================================================================
      implicit none

      Integer LDA, nRow, k, i
      Real*8  A(LDA, *), F(*), bc

C ======================================================================
      Do i = 1, nRow
         F(i) = F(i) - A(i, k) * bc

         A(i, k) = 0D0
         A(k, i) = 0D0
      End do

      F(k) = bc
      A(k, k) = 1D0

      Return
      End


C ======================================================================
      Subroutine applyNEU(Floc, Fglb, k, FEM_X)
C ======================================================================
C Apply the Neumann b.c. on edge k to the local right-hand side vector
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Real*8  Floc(*), Fglb(*)
      Integer k, FEM_X

      Integer iref(4), l, m, n
      DATA    iref/1,2,3,1/

C ======================================================================
      If(FEM_X.EQ.FEM_P1) Then
         l = iref(k + 1)
         Fglb(k) = Fglb(k) + Floc(1)
         Fglb(l) = Fglb(l) + Floc(2)

      Else If(FEM_X.EQ.FEM_P2) Then
         l = iref(k + 1)
         m = k + 3
         Fglb(k) = Fglb(k) + Floc(1)
         Fglb(l) = Fglb(l) + Floc(2)
         Fglb(m) = Fglb(m) + Floc(3)

      Else If(FEM_X.EQ.FEM_P3) Then
         l = iref(k + 1)
         m = k + 3
         n = k + 6
         Fglb(k) = Fglb(k) + Floc(1)
         Fglb(l) = Fglb(l) + Floc(2)
         Fglb(m) = Fglb(m) + Floc(3)
         Fglb(n) = Fglb(n) + Floc(3)
      End if

      Return
      End


C ======================================================================
      Subroutine applyDIRedge(LDA, nRow, A, F, i1, i2, XY1, XY2, FEM_PX,
     &                        label, Dbc, dDATA, iDATA, iSYS)
C ======================================================================
C Apply the Dirichlet b.c. to the local matrix
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer  LDA, nRow, i1, i2, FEM_PX
      Real*8   dDATA(*), A(LDA, *), F(*), XY1(*), XY2(*)
      Integer  iDATA(*), label, Dbc, iSYS(*)
      EXTERNAL Dbc

      Integer  iref(4), i, k, m, n
      Real*8   d, x, y, eBC(1)
      DATA     iref/1,2,3,1/

C ======================================================================
      i1 = k
      i2 = iref(i1 + 1)

      n = 0  ! number of interior points
      If(FEM_PX.EQ.FEM_P2) n = 1
      If(FEM_PX.EQ.FEM_P3) n = 2

      Do i = 1, n + 1
         d = dble(i - 1) / (n + 1)
         m = i1 + 3*(i - 1)

         x = XY1(1) * (1D0-d) + XY2(1) * d
         y = XY1(2) * (1D0-d) + XY2(2) * d

         k = Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
         Call applyDIR(LDA, nRow, A, F, m, eBC(1))
      End do

      x = XY2(1)
      y = XY2(2)

      k = Dbc(x, y, label, dDATA, iDATA, iSYS, eBC)
      Call applyDIR(LDA, nRow, A, F, i2, eBC(1))

      Return
      End


C ======================================================================
      Logical Function ifXbc(clr, iXbc)
C ======================================================================
      Integer clr, iXbc

      ifXbc = IAND(clr, iXbc) .EQ. iXbc

      Return
      End


C ======================================================================
      Logical Function ifXtensor(clr, iXtensor)
C ======================================================================
      Integer clr, iXtensor

      ifXtensor = IAND(clr, iXtensor) .EQ. iXtensor

      Return
      End



