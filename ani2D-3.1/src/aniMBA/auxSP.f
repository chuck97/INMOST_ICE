C ======================================================================
C  All operations with colors are binary operations. 
C  The same operations go with the input parameter status.
C ======================================================================


C ======================================================================
C     Logical Function ifXstatus(status, iXstatus)
      Logical Function ifXnode(clr, iXnode)
C ======================================================================
      Integer clr, iXnode

      ifXnode = IAND(clr, iXnode) .EQ. iXnode

      Return
      End



C ======================================================================
C     Subroutine addXstatus(status, iXstatus)
      Subroutine addXnode(clr, iXnode)
C ======================================================================
      Integer clr, iXnode

      clr = IOR(clr, iXnode)

      Return
      End



C ======================================================================
C     Subroutine delXstatus(status, iXstatus)
      Subroutine delXnode(clr, iXnode)
C ======================================================================
      Integer clr, iXnode

      clr = clr - IAND(clr, iXnode)

      Return
      End



C ======================================================================
      Integer Function minClr(clr1, clr2)
C ======================================================================
C  The function returns common color for both clr1 and clr2.
C ======================================================================
      Integer clr1, clr2

      minClr = IAND(clr1, clr2)

      Return
      End



C ======================================================================
      Integer Function maxClr(clr1, clr2)
C ======================================================================
C Routine returns minimal color containing both clr1 and clr2.
C ======================================================================
      Integer clr1, clr2

      maxClr = IOR(clr1, clr2)

      Return
      End




