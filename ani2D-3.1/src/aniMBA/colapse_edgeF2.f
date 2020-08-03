C ======================================================================
      Subroutine clpsF2(
C ======================================================================
     &           iwP, iwE,
c group (M)
     &           nP, XYP,
     &           nF, IPF, lbF, Crv, lbC, CrvFunction,
     &           nE, IPE,
c group (M-EXT)
     &           ICP, IEP, IFE, IEE,
     &           IHolP, IHolF, IHolE,
     &           nPw, nEw, XYPw, IPEw, HesPw, iSE, rSE,
C group (Q)
     &           hStar, 
     &           rQuality, HesP, detG,
     &           L1E, L2E, qE, nL2, nStep,
c group (CONTROL)
     &           status, MetricFunction, flagAnalytic,
c group (CRV)
     &           L1Et, L2Et, tE, nL2t, 
     &           LFnc, ILt, nCrvFnc,
C group (W)
     &           lFu, lEu, iFu, iEu, 
     &           flag)
C ======================================================================
      implicit none
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
C ======================================================================
C Routine realizes one of the mesh operations: collapses an edge ending 
C at point iwP of element iwE. It uses the main routine clpsF1.
C ======================================================================
      Integer  iwP, iwE

C group (M)
      Integer  nP
      Real*8   XYP(2, *)

      Integer  nF, IPF(2, *), lbF(*), lbC(*)
      Real*8   Crv(2, *)
      EXTERNAL CrvFunction

      Integer  nE, IPE(3, *)

c group (M-EXT)
      Integer  ICP(*), IEP(*), IFE(3, *), IEE(3, *)
      Integer  IHolP(*), IHolF(*), IHolE(*)

      Integer  nPw, nEw, IPEw(3, *), iSE(*)
      Real*8   XYPw(2, *), HesPw(3, *), rSE(*)

C group (Q)
      Real*8   hStar, rQuality, HesP(3, *), detG(*), qE(*)
      Integer  L1E(2, *), L2E(*), nL2, nStep(4)

c group (CONTROL)
      Integer  MetricFunction, status
      EXTERNAL MetricFunction

      Logical  flagAnalytic

c group (CRV)
      Integer  L1Et(2, *), L2Et(*), nL2t(*), LFnc(*), ILt(*), nCrvFnc
      Real*8   tE(*)

C group (W)
      Integer  lFu, lEu, iFu(*), iEu(*) 
      Logical  flag


C LOCAL VARIABLES
      Integer  lFs, lEs
      Integer  iFs(MaxS), iEs(MaxS), IPFs(2, MaxS), IPEs(3, MaxS)
      Integer  iFx(MaxS), iEx(MaxS), IPFx(2, MaxS), IPEx(3, MaxS)
      Real*8   qEs(MaxS), qEx(MaxS)

      Integer  ip(4)
      Real*8   d, w1, w2
      Logical  ifXnode

      Integer  i, k, i1,i2,i3, iP1,iP2,iP3, iPa,iPb, iwF, itF, kE

C ======================================================================
      flag = .FALSE.

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

C ... checking points which can not be involved
      iP1 = IPE(iwP, iwE)
      If(ifXnode(ICP(iP1), jSnode)) Goto 1000
      If(ifXnode(ICP(iP1), jTnode)) Goto 1000

      i2 = ip(iwP + 1)
      i3 = ip(i2 + 1)
      iPa = IPE(i2, iwE)
      iPb = IPE(i3, iwE)

      Do 10 k = 1, lEu
         kE = iEu(k)
         If(kE.EQ.iwE) Goto 10

         Do i = 1, 3
            If(IPE(i, kE).EQ.iP1) Then
               i1 = i
               Goto 5
            End if
         End do
         Goto 10

 5       i2 = ip(i1 + 1)
         i3 = ip(i2 + 1)

         iP2 = IPE(i2, kE)
         iP3 = IPE(i3, kE)

         d = (XYP(1, iP2) - XYP(1, iP1)) *
     &       (XYP(2, iP3) - XYP(2, iP1)) -
     &       (XYP(1, iP3) - XYP(1, iP1)) *
     &       (XYP(2, iP2) - XYP(2, iP1))

         If(d.LT.0D0) Then
            iwF = min(i1, i3)
            itF = max(i1, i3)
            If(iP3.EQ.iPa .OR. iP3.EQ.iPb) Goto 10
         Else
            iwF = min(i1, i2)
            itF = max(i1, i2)
            If(iP2.EQ.iPa .OR. iP2.EQ.iPb) Goto 10
         End if
         If(iwF.EQ.1 .AND. itF.EQ.3) iwF = 3

         If(iwF.EQ.i1) Then
            w1 = 0D0
            w2 = 1D0
         Else
            w1 = 1D0
            w2 = 0D0
         End if

         Call makSE(kE,  IEP, IPF, IPE, IFE,  IEE,  qE, MaxS,
     &              lFs, lEs, iFs, iEs, IPFs, IPEs, qEs,
     &              status)

         Call clpsF1(
     &           iwF, kE, w1, w2,
c group (M)
     &           nP, XYP, 
     &           nF, IPF, lbF, Crv, lbC, CrvFunction,
     &           nE, IPE,  
c group (M-EXT)
     &           ICP, IEP, IFE, IEE,
     &           IHolP, IHolF, IHolE,
     &           nPw, nEw, XYPw, IPEw, HesPw, iSE, rSE,
C group (Q)
     &           hStar, 
     &           rQuality, HesP, detG,
     &           L1E, L2E, qE, nL2, nStep,
c group (CONTROL)
     &           status, MetricFunction, flagAnalytic,
c group (CRV)
     &           L1Et, L2Et, tE, nL2t, 
     &           LFnc, ILt, nCrvFnc,
C group (W)
     &           lFs, lEs, iFs, iEs, IPFs, IPEs, qEs,
     &                     iFx, iEx, IPFx, IPEx, qEx,
     &           flag)

        If(flag) Goto 1000
 10   Continue

 1000 Return
      End


