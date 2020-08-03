C ======================================================================
      Subroutine quality_log(nE, L1E, L2E, qE)
C ======================================================================
      Integer L1E(2, *), L2E(*)
      Real*8  qE(*)

C LOCAL VARIABLES
      Integer  nBins, cBins(10)
      Real*8   vBins(10)

C ======================================================================
      nBins = 9
      
      Do n = 1, nBins
         vBins(n) = 1D-1**(nBins - n)
      End do

      Call statistics_bins(nE, L1E, L2E, qE, nBins, vBins, cBins)
      
      Return
      End



C ======================================================================
      Subroutine quality_plain(nE, L1E, L2E, qE)
C ======================================================================
      Integer L1E(2, *), L2E(*)
      Real*8  qE(*)

C LOCAL VARIABLES
      Integer  nBins, cBins(10)
      Real*8   vBins(10)

C ======================================================================
      nBins = 9
      
      Do n = 1, nBins
         vBins(n) = 1D-1 * (n + 1)
      End do

      Call statistics_bins(nE, L1E, L2E, qE, nBins, vBins, cBins)
      
      Return
      End



C ======================================================================
      Subroutine statistics_bins(nE, L1E, L2E, qE, nBins, vBins, cBins)
C ======================================================================
      implicit none

      Integer  nE, L1E(2, *), L2E(*), nBins, cBins(*)
      Real*8   qE(*), vBins(*)

      Integer  i,k,n, iE
C ======================================================================
      Do n = 1, nBins
         cBins(n) = 0
      End do

      iE = L2E(1)
      Do n = 1, nE
         If(qE(iE).LT.vBins(1)) Then
            cBins(1) = cBins(1) + 1
         Else
            Do k = 2, nBins 
               If(vBins(k-1).LE.qE(iE) .AND. qE(iE).LT.vBins(k)) Then
                  cBins(k) = cBins(k) + 1
                  Goto 500
               End if
            End do
         End if

 500     iE = L1E(2, iE)
      End do

      Write(*,5001) (vBins(i), i = 1, nBins)
      Write(*,5002) (cBins(i), i = 1, nBins)

 5001 Format(/, 'Distribution of triangles by quality',/, E7.1,100E9.1)
 5002 Format(I7,100I9)
      Return
      End

