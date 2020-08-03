C ======================================================================
      Subroutine saveMgmv(nP, nE, XYP, IPE, fName)
C ======================================================================
C Routine saves mesh in the GMV file. 
C
C *** Remarks:
C        1. The size of the working memory is nE
C ======================================================================
      implicit none

      Integer nP, nE
      Real*8  XYP(2, *)
      Integer IPE(3, *)
      Character*(*) fName

C group (Local variables)
      Integer      i, n
      Real*8       z
      Character*30 fNameExt

C ======================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.gmv')
         i = i + 1
      End do

      fNameExt = fName(1:i) // 'gmv'
      Write(*,'(A,A)') 'Saving GMV file ', fNameExt

      Open(10, file=fNameExt, status='UNKNOWN')

      Write(10, '(A)') 'gmvinput ascii'
      Write(10, *)
      Write(10, *) 'nodev ', nP

      z = 0D0
      Do n = 1, nP
         Write(10, *) (XYP(i, n), i = 1, 2), z
      End do

c ... save faces
      Write(10, '(A,I10)') 'cells ', nE
      Do n = 1, nE
         Write(10, *) 'general 1'
         Write(10, *) '3 ', (IPE(i, n), i = 1, 3)
      End do

      Write(10, '(A)') 'endgmv'
      Close(10)

      Return
      End



C ======================================================================
      Subroutine saveDgmv(nP, nE, XYP, IPE, U, fName)
C ======================================================================
C Routine saves mesh and data in the GMV file. 
C
C *** Remarks:
C        1. The size of the working memory is nE
C ======================================================================
      implicit none

      Integer nP, nE
      Real*8  XYP(2, *), U(*)
      Integer IPE(3, *)
      Character*(*) fName

C group (Local variables)
      Integer      i, n
      Real*8       z
      Character*30 fNameExt

C ======================================================================
      i = 1
      Do while( fName(i:i+3) .NE. '.gmv')
         i = i + 1
      End do

      fNameExt = fName(1:i) // 'gmv'
      Write(*,'(A,A)') 'Saving GMV file ', fNameExt

      Open(10, file=fNameExt, status='UNKNOWN')

      Write(10, '(A)') 'gmvinput ascii'
      Write(10, *)
      Write(10, *) 'nodev ', nP

      z = 0D0
      Do n = 1, nP
         Write(10, *) (XYP(i, n), i = 1, 2), z
      End do

c ... save faces
      Write(10, '(A,I10)') 'cells ', nE
      Do n = 1, nE
         Write(10, *) 'general 1'
         Write(10, *) '3 ', (IPE(i, n), i = 1, 3)
      End do

c ... save data
      Write(10,*)
      Write(10,'(A)') 'variable'
      Write(10,'(A)') 'solution 1'
      Do i = 1, nP
         Write(10,'(F14.5)') U(i)
      End do
      Write(10,'(A)') 'endvars'

      Write(10, '(A)') 'endgmv'
      Close(10)

      Return
      End



