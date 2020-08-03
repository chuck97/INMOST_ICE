C ======================================================================
      Subroutine stream_function(nP, nR, nF, nE, XYP, IPF, IPE,
     &                           U, PSI, MaxWi, iW)
C ======================================================================
C Routine calculates the P1 finite element stream function of P2 finite
C element velocity field U.
C
C The size of working arrya iW should be at least nP.
C ======================================================================
      implicit none

      Integer nP, nR, nF, nE, IPF(2, *), IPE(3, *)
      Real*8  XYP(2, *), U(*), PSI(*)
    
      Integer MaxWi, iW(*)

c Local variables
      Integer i,k,n, iP1,iP2, iux,iuy
      Real*8  tau(2), uavg(2), s
      Logical flag

C ======================================================================
c ... mark all points, initialize the stream function
      Do n = 1, nP
         iW(n)  = 1
         PSI(n) = 0D0
      End do

      iux = 0
      iuy = iux + nP


c ... process one boundary points
      iP1 = IPF(1, 1)
      If(nF.EQ.0) iP1 = IPE(1, 1)
      iW(iP1) = 0


c ... loop until all points are re-calculated
      flag = .TRUE.
      Do while( flag )
         flag = .FALSE.

         Do 200 n = 1, nE
            Do i = 1, 3
               iP1 = IPE(i, n)
               If(iW(iP1).EQ.0) Goto 100
            End do
            Goto 200

 100        Continue
            Do i = 1, 3
               iP2 = IPE(i, n)
               If(iW(iP2).EQ.1) Then
                  iW(iP2) = 0
                  flag = .TRUE. 

                  Do k = 1, 2
                     tau(k) = XYP(k, iP2) - XYP(k, iP1) 
                  End do

                  uavg(1) = (U(iux + iP1) + U(iux + iP2)) / 2
                  uavg(2) = (U(iuy + iP1) + U(iuy + iP2)) / 2

                  s = tau(2) * uavg(1) - tau(1) * uavg(2)
                  PSI(iP2) = PSI(iP1) + s
               End if
            End do
 200     Continue
      End do

      Return
      End


