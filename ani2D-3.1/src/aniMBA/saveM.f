C ======================================================================
      Subroutine saveMani(nP, nPfix, XYP, lbP, fixP,
     &                    nF, nFfix, IPF, lbF, fixF,
     &                    nC,        Crv, lbC, 
     &                    nE, nEfix, IPE, lbE, fixE,
     &      fName)
C ======================================================================
C  Routine saves the mesh in the file "fName".
C  The file name must have extension .ani.
C ======================================================================
      implicit none

      Integer  nP, nPfix, lbP(*), fixP(*)
      Real*8  XYP(2, *)

      Integer  nF, nFfix, IPF(2, *), lbF(*), fixF(*)

      Integer  nC, lbC(*)
      Real*8   Crv(2, *)

      Integer  nE, nEfix, IPE(3, *), lbE(*), fixE(*)

      Character*(*) fName


C LOCAL VARIABLES
      Character*30  fNameExt

      Integer  i, j, k, n, kps, kis, kes, kfs, kfe, kcs, kfp, kff

C ======================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.ani')
         i = i + 1
      End do

      fNameExt = fName(1:i) // 'ani'

      Write(*,'(/,A,A)') 'Saving mesh ', fNameExt

c ... save the mesh header
      kps = 10
      kfs = kps + nP + 2
      kes = kfs + nF + 2
      kcs = kes + nE + 2

      kfp = kcs + nC + 2
      kff = kfp + nPfix + 2
      kfe = kff + nFfix + 2

      kis = kfe + nEfix + 2

      Open(10, file=fNameExt, status='UNKNOWN')
      Write(10, '(3(A,I7),A)') 'T points:        ', 
     &          nP, ' (lines ', kps, ' - ', kps + nP - 1, ')'
      Write(10, '(3(A,I7),A)') 'T edges:         ',
     &          nF, ' (lines ', kfs, ' - ', kfs + nF - 1, ')'
      Write(10, '(3(A,I7),A)') 'T elements:      ',
     &          nE, ' (lines ', kes, ' - ', kes + nE - 1, ')'

      If(nC.NE.0) Then
         Write(10, '(3(A,I7),A)') 'T curved edges:  ',
     &             nC, ' (lines ', kcs, ' - ', kcs + nC - 1, ')'
      Else
         Write(10, '(A,I7)') 'T curved faces:  ', nC
      End if

      If(nPfix.NE.0) Then
         Write(10, '(3(A,I7),A)') 'T fixed points:  ',
     &             nPfix, ' (lines ', kfp, ' - ', kfp + nPfix - 1, ')'
      Else
         Write(10, '(A,I7)') 'T fixed points:  ', nPfix
      End if

      If(nFfix.NE.0) Then
         Write(10, '(3(A,I7),A)') 'T fixed edges:   ', 
     &             nFfix, ' (lines ', kff, ' - ', kff + nFfix - 1, ')'
      Else
         Write(10, '(A,I7)') 'T fixed edges:   ', nFfix
      End if

      If(nEfix.NE.0) Then
         Write(10, '(3(A,I7),A)') 'T fixed elements:', 
     &             nEfix, ' (lines ', kfe, ' - ', kfe + nEfix - 1, ')'
      Else
         Write(10, '(A,I7)') 'T fixed elements:', nEfix 
      End if


c ... save the mesh
      Write(10,*)
      Write(10,*) nP, ' # of nodes'
      Do n = 1, nP
         Write(10,*) (XYP(j, n), j = 1, 2), lbP(n)
      End do

      Write(10,*)
      Write(10,*) nF, ' # of edges'
      k = 0
      Do n = 1, nF
         Write(10,*) (IPF(j, n), j=1,2), lbF(n)
      End do

      Write(10,*)
      Write(10,*) nE, ' # of elements'
      Do n = 1, nE
         Write(10,*) (IPE(j, n), j=1,3), lbE(n)
      End do


c ... saving curvilinear edges
      nC = 0
      Do n = 1, nF
         If(lbC(n).GT.0) nC = nC + 1
      End do

      Write(10,*)
      Write(10,*) nC, ' # of curvilinear faces'
      Do n = 1, nF
         If(lbC(n).GT.0) 
     &      Write(10,*) n, (Crv(j, n), j=1,2), lbC(n)
      End do


c ... saving the fixed mesh points, faces and elements
      Write(10,*)
      Write(10,*) nPfix, ' # number of fixed points'
      Do n = 1, nPfix
         Write(10,*) fixP(n)
      End do


      Write(10,*)
      Write(10,*) nFfix, ' # number of fixed faces'
      Do n = 1, nFfix
         Write(10,*) fixF(n)
      End do


      Write(10,*)
      Write(10,*) nEfix, ' # number of fixed elements'
      Do n = 1, nEfix
         Write(10,*) fixE(n)
      End do

      Close(10)

      Return
      End





