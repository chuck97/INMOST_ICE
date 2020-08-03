C ======================================================================
      Subroutine loadMani(nP, nPfix, MaxP, XYP, lbP, fixP,
     &                    nF, nFfix, MaxF, IPF, lbF, fixF,
     &                    nC,              Crv, lbC, 
     &                    nE, nEfix, MaxE, IPE, lbE, fixE,
     &                    fName)
C ======================================================================
C  Routine reads the input mesh from file "fName"
C  The file name MUST have extension .ani.
C ======================================================================
      implicit none

C group (M)
      Integer  nP, nPfix, MaxP, lbP(*), fixP(*)
      Real*8   XYP(2, *)

      Integer  nF, nFfix, MaxF, IPF(2, *), lbF(*), fixF(*)

      Integer  nC, lbC(*)
      Real*8   Crv(2, *)

      Integer  nE, nEfix, MaxE, IPE(3, *), lbE(*), fixE(*)

      Character*(*) fName


C LOCAL VARIABLES
      Logical flagP, flagF, flagE, flagC, flagPV, flagFV, flagEV
      Character*30  fNameExt

      Integer  i, j, k, n
      Real*8   a, b

C ======================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.ani')
         i = i + 1
      End do
      fNameExt = fName(1:i) // 'ani'

      Write(*,'(A,A)') 'Loading mesh ', fNameExt

      Open(10, file=fNameExt, status='OLD', ERR=1000)
      Read(10,*) flagP 
      Read(10,*) flagF
      Read(10,*) flagE
      Read(10,*) flagC
      Read(10,*) flagPV
      Read(10,*) flagFV
      Read(10,*) flagEV

      Read(10,*)
      Read(10,*) nP
      If(nP.GT.MaxP) Then
        Call errMesMBA(1003, 'loadMani', 
     &                       'local parameter MaxP is small')
      End if
      Do n = 1, nP
         Read(10,*) (XYP(j, n), j=1,2), lbP(n)
      End do


c ... reading the boundary edges
      Read(10,*)
      Read(10,*) nF
      If(nF.GT.MaxF) Call errMesMBA(1004, 'loadMani',
     &                          'local parameter MaxF is small')
      Do n = 1, nF
         Read(10,*) (IPF(j, n), j=1,2), lbF(n)
      End do


      Read(10,*)
      Read(10,*) nE
      If(nE.GT.MaxE) Call errMesMBA(1006, 'loadMani',
     &                          'local parameter MaxE is small')
      Do n = 1, nE
         Read(10,*) (IPE(j, n), j=1,3), lbE(n)
      End do


c ... reading parametrization of the curved boundary edges:
c     edge number
c     parametrization of edge end-points
c     parametrization function
      Read(10,*)
      Read(10,*) nC
      If(nC.GT.MaxF) Call errMesMBA(1004, 'loadMani',
     &                            'local parameter MaxF is small')
      Do n = 1, nF
         lbC(n) = 0
      End do

      Do n = 1, nC
         Read(10,*) i, a, b, k

         Crv(1, i) = a 
         Crv(2, i) = b
         lbC(i) = k
      End do


c ... reading the fixed mesh points, faces and elements
      Read(10,*)
      Read(10,*) nPfix
      If(nPfix.GT.MaxP) Call errMesMBA(1003, 'loadMani', 
     &                            'local parameter MaxP is small')
      Do n = 1, nPfix
         Read(10,*) fixP(n)
      End do


      Read(10,*)
      Read(10,*) nFfix
      If(nFfix.GT.MaxF) Call errMesMBA(1004, 'loadMani', 
     &                            'local parameter MaxF is small')
      Do n = 1, nFfix
         Read(10,*) fixF(n)
      End do


      Read(10,*)
      Read(10,*) nEfix
      If(nEfix.GT.MaxE) Call errMesMBA(1006, 'loadMani', 
     &                            'local parameter MaxE is small')
      Do n = 1, nEfix
         Read(10,*) fixE(n)
      End do

      Close(10)

      Return

 1000 Continue
      Call errMesMBA(4001, 'loadMani', 'Input file name is wrong ')
      Return
      End



C ======================================================================
      Subroutine loadMani_24(nP, nPfix, MaxP, XYP, lbP, fixP,
     &                       nF, nFfix, MaxF, IPF, lbF, fixF,
     &                       nC,              Crv, lbC, 
     &                       nE, nEfix, MaxE, IPE, lbE, fixE,
     &                       fName)
C ======================================================================
C  Routine reads the input mesh in format 2.4 from file "fName"
C  The file name MUST have extension .ani.
C ======================================================================
      implicit none

C group (M)
      Integer  nP, nPfix, MaxP, lbP(*), fixP(*)
      Real*8   XYP(2, *)

      Integer  nF, nFfix, MaxF, IPF(2, *), lbF(*), fixF(*)

      Integer  nC, lbC(*)
      Real*8   Crv(2, *)

      Integer  nE, nEfix, MaxE, IPE(3, *), lbE(*), fixE(*)

      Character*(*) fName


C LOCAL VARIABLES
      Logical  flagP, flagF, flagE, flagC, flagPV, flagFV, flagEV
      Character*30  fNameExt

      Integer  i, j, k, n, iCrv, iW(5000)

C ======================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.ani')
         i = i + 1
      End do
      fNameExt = fName(1:i) // 'ani'

      Write(*,'(A,A)') 'Loading mesh ', fNameExt

      Open(10, file=fNameExt, status='OLD', ERR=1000)
      Read(10,*) flagP 
      Read(10,*) flagF
      Read(10,*) flagE
      Read(10,*) flagC
      Read(10,*) flagPV
      Read(10,*) flagFV
      Read(10,*) flagEV

      Read(10,*)
      Read(10,*) nP
      If(nP.GT.MaxP) Call errMesMBA(1003, 'loadMani',
     &                          'local parameter MaxP is small')
      Do n = 1, nP
         Read(10,*) (XYP(j, n), j=1,2)
         lbP(n) = 0
      End do


c ... reading the boundary edges: 
c        edge end-points
c        the corresponding curved edge
c        dummy integer
c        label of the edge 
      Read(10,*)
      Read(10,*) nF
      If(nF.GT.MaxF) Call errMesMBA(1004, 'loadMani',
     &                          'local parameter MaxF is small')
      Do n = 1, nF
         Read(10,*) (IPF(j, n), j=1,2), iW(n), k, lbF(n)
      End do


      Read(10,*)
      Read(10,*) nE
      If(nE.GT.MaxE) Call errMesMBA(1006, 'loadMani',
     &                          'local parameter MaxE is small')
      Do n = 1, nE
         Read(10,*) (IPE(j, n), j=1,3), lbE(n)
      End do


c ... reading parametrization of the curved boundary edges
c     parametrization of edge end-points
c     parametrization function
      Read(10,*)
      Read(10,*) nC
      If(nC.GT.MaxF) Call errMesMBA(1004, 'loadMani',
     &                            'local parameter MaxF is small')
      Do n = 1, nC
         Read(10,*) (Crv(j, n), j=1,2), lbC(n)
      End do

C ... update Crv & lbC
      Do n = nF, 1, -1
         iCrv = iW(n)
         If(iCrv.NE.0) Then
            lbC(n) = lbC(iCrv)
            If(n.NE.iCrv) lbC(iCrv) = 0

            Do i = 1, 2
               Crv(i, n) = Crv(i, iCrv)
            End do
         Else
            lbC(n) = 0
            Do i = 1, 2
               Crv(i, n) = 0D0
            End do
         End if
      End do


c ... reading the fixed mesh points, faces and elements
      Read(10,*) 
      Read(10,*) nPfix
      If(nPfix.GT.MaxP) Call errMesMBA(1003, 'loadMani', 
     &                          'local parameter MaxP is small')
      Do n = 1, nPfix
         Read(10,*) fixP(n)
      End do


      Read(10,*)
      Read(10,*) nFfix
      If(nFfix.GT.MaxF) Call errMesMBA(1004, 'loadMani', 
     &                          'local parameter MaxF is small')
      Do n = 1, nFfix
         Read(10,*) fixF(n)
      End do


      Read(10,*)
      Read(10,*) nEfix
      If(nEfix.GT.MaxE) Call errMesMBA(1006, 'loadMani', 
     &                          'local parameter MaxE is small')
      Do n = 1, nEfix
         Read(10,*) fixE(n)
      End do

      Close(10)

      Return

 1000 Continue
      Call errMesMBA(4001, 'loadMani', 'Input file name is wrong ')
      Return
      End



C ======================================================================
      Subroutine loadManiHeader(nP, nF, nC, nE,
     &                          nPfix, nFfix, nEfix, fName)
C ======================================================================
C  Routines read the header of te input mesh from file fName. 
C  The must have extension .ani.
C ======================================================================
      Integer       nP, nF, nC, nE, nPfix, nFfix, nEfix
      Character*(*) fName

      Logical       flag
      Character*30  fNameExt, text, text2

C ==========================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.ani')
         i = i + 1
      End do
      fNameExt = fName(1:i) // 'ani'

      Write(*,'(A,A)') 'Reading header of mesh ', fNameExt

      Open(10, file=fNameExt, status='OLD', ERR=1000)
      Read(10,*) flag, text, nP
      Read(10,*) flag, text, nF
      Read(10,*) flag, text, nE
      Read(10,*) flag, text, text2, nC
      Read(10,*) flag, text, text2, nPfix
      Read(10,*) flag, text, text2, nFfix
      Read(10,*) flag, text, text2, nEfix

      Close(10)
      Return

 1000 Continue
      Call errMesMBA(4001, 'loadManiHeader', 'File name is wrong')
      Return
      End


