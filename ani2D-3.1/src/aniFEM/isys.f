C ======================================================================
      Subroutine encodeISYS(iE, iP1, iP2, iP3,  
     &                          iR1, iR2, iR3, 
     &                          iE1, iE2, iE3, nP, nR, nE, iSYS)
C ======================================================================
      implicit none
      include 'assemble.fd'

      Integer iE, iP1,iP2,iP3, iR1,iR2,iR3, iE1,iE2,iE3, nP, nR, nE
      Integer iSYS(MAXiSYS)
C ======================================================================
      iSYS(1) = 0

c ... triangle number
      iSYS(3) = iE

c ... vertices of this triangle 
      iSYS(4) = iP1
      iSYS(5) = iP2
      iSYS(6) = iP3

c ... edges of this triangle 
      iSYS(7) = iR1
      iSYS(8) = iR2
      iSYS(9) = iR3

c ... total number of mesh objects
      iSYS(10) = nP
      iSYS(11) = nR
      iSYS(12) = nE

c ... neighboring triangles
      iSYS(13) = iE1
      iSYS(14) = iE2
      iSYS(15) = iE3

      Return
      End



C ======================================================================
      Subroutine decodeISYS(iE, iP1, iP2, iP3,  
     &                          iR1, iR2, iR3, 
     &                          iE1, iE2, iE3, nP, nR, nE, iSYS)
C ======================================================================
      implicit none
      include 'assemble.fd'

      Integer iE, iP1,iP2,iP3, iR1,iR2,iR3, iE1,iE2,iE3, nP, nR, nE
      Integer iSYS(MAXiSYS)
C ======================================================================
c ... triangle number
      iE = iSYS(3)

c ... vertices of this triangle 
      iP1 = iSYS(4)
      iP2 = iSYS(5)
      iP3 = iSYS(6)

c ... edges of this triangle 
      iR1 = iSYS(7)
      iR2 = iSYS(8)
      iR3 = iSYS(9)

c ... total number of mesh objects
      nP = iSYS(10)
      nR = iSYS(11)
      nE = iSYS(12)

c ... neigboring triangles
      iE1 = iSYS(13)
      iE2 = iSYS(14)
      iE3 = iSYS(15)

      Return
      End


