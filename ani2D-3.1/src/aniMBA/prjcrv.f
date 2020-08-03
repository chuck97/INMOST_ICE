C ======================================================================
      Subroutine prjCrv(XY, prjXY, iFNC, t, CrvFunction,
     &                  L1Et, L2Et, nL2t, tE)
C ======================================================================
C Routine computes a point which lies on the boundary of the domain and 
C is the closest to XY. The later is the point lying on new boundary. 
C ======================================================================
      Real*8  XY(2), prjXY(2), t

      EXTERNAL CrvFunction

      Real*8  tE(*)
      Integer L1Et(2, *), L2Et(*)

      Real*8  XYprev(2), XYnext(2), tnext, tprev
      Integer iPrevL2t, PrevL2, iPrevL1t, PrevL1

C ======================================================================
      If(iFNC.EQ.0) Then
         prjXY(1) = XY(1)
         prjXY(2) = XY(2)
      Else
         iPrevL2t = PrevL2(nL2t, L2Et, tE, t)
         iPrevL1t = PrevL1(L1Et, L2Et, iPrevL2t, tE, t)

         tprev = tE(iPrevL1t)
         tnext = tE(L1Et(2, iPrevL1t))

         Call  aniCrv(tprev, XYprev, iFNC, CrvFunction)
         Call  aniCrv(tnext, XYnext, iFNC, CrvFunction)

         prjXY(1) = (XYprev(1) * (tnext - t) / (tnext - tprev) +
     &               XYnext(1) * (t - tprev) / (tnext - tprev))
         prjXY(2) = (XYprev(2) * (tnext - t) / (tnext - tprev) +
     &               XYnext(2) * (t - tprev) / (tnext - tprev))
      End if

      Return
      End



C ======================================================================
      Subroutine calCrvFnc(nF, lbC, LFnc, nCrvFnc)
C ======================================================================
C Routine computes the number of different functions describing 
C curvilinear boundaries and internal interfaces.
C ======================================================================
      implicit none
      Integer  nF, lbC(*), LFnc(*), nCrvFnc

      Integer  k, n, iC

C ======================================================================
      nCrvFnc = 0

      Do n = 1, nF
         iC = lbC(n)

         If(iC.NE.0) Then
            Call findSE(nCrvFnc, LFnc, iC, k)

            If(k.EQ.0) Then
               nCrvFnc = nCrvFnc + 1
               LFnc(nCrvFnc) = iC
            End if
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine tEMak(tE, nEt, nF, Crv, lbC, iCrvFnc)
C ======================================================================
C Routine creates non-ordered list of parameters for boundary part 
C described by iCrvFnc-th function.
C ======================================================================
      implicit none

      Integer  nEt, nF, lbC(*), iCrvFnc
      Real*8   tE(*), Crv(2, *)


C LOCAL VARIABLES
      Integer  k, n, iC
      Real*8   tmin
      Logical  flag

C ======================================================================
c ... compute minimal parameter for boundary iFnc
      flag = .TRUE.

      Do n = 1, nF
         iC = lbC(n)
 
         If(iC.EQ.iCrvFnc) Then
               Do k = 1, 2
                  If(flag) Then
                  tmin = Crv(k, n)
                     flag = .FALSE.
                  Else
                  tmin = min(tmin, Crv(k, n))
                  End if
               End do
            End if
      End do


c ... collect maxima of two parameters in tE 
      nEt = 1
      tE(nEt) = tmin  

      Do n = 1, nF
         iC = lbC(n)

         If(iC.EQ.iCrvFnc) Then
            nEt = nEt + 1
            tE(nEt) = max(Crv(1, n), Crv(2, n))
         End if
      End do

      Return
      End

