C ======================================================================
      Subroutine movePoint(
C ======================================================================
     &           iwP, iwE,
c group (M)
     &           XYP, 
     &           IPF, lbF, 
     &           Crv, lbC, CrvFunction,
     &           nE, IPE,  
c group (M-EXT)
     &           ICP, 
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
     &           lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &                     iFs, iEs, IPFs, IPEs, qEs,
     &           icnt, rMove, flag)
C ======================================================================
      implicit none
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
      include 'magic.fd'
C ======================================================================
C Routine moves a vertex of triangle iwE to increase the quality of 
C mesh elements in the superelement associated with the vertex.
C ======================================================================
      Integer  iwP, iwE

C group (M)
      Real*8   XYP(2, *)

      Integer  IPF(2, *), lbF(*), lbC(*)
      Real*8   Crv(2, *)
      EXTERNAL CrvFunction

      Integer  nE, IPE(3, *)

c group (M-EXT)
      Integer  ICP(*)

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
      Integer  lFu, lEu, iFu(*), iEu(*), IPFu(2, *), IPEu(3, *), icnt
      Integer            iFs(*), iEs(*), IPFs(2, *), IPEs(3, *)
      Real*8   qEu(*), qEs(*), rMove
      Logical  flag


C LOCAL VARIABLES
      Real*8   XYPs(2, 2), HesPs(3), detGs, prjXYPs(2, 2)

C ... for nonlinear minimization procedure
      Real*8   NLnFnc, ZZ(2), U(2), U1(2)

      Integer  ip(4), iCRVs(2), iFNCs(2)
      Real*8   distSP, tri_area, edge_length, heit, v, d
      Real*8   par(4), q(3), qMin, hMin, hMax, xStep, yStep, tStep
      Real*8   tc(2), t1, t2, t3

      Logical  ifXnode, flagTM

      Integer  iP1,iP2,iP3, iPa,iPb, iF1,iF2, nF1,nF2, iE
      Integer  ICPs, iBNDs, iFmove
      Integer  i,j,k,n, i1,i2,i3, lF, lE, nBad, nU, iPend, iref, ir

C ======================================================================
      flag = .FALSE.

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      iP1 = IPE(iwP, iwE)
      ICPs = ICP(iP1)


C ... checking the case when moving is impossible
      If(ifXnode(ICPs, jVnode)) Goto 1000
      If(ifXnode(ICPs, jTnode)) Goto 1000

      Do i = 1, 3
         HesPs(i) = HesP(i, iP1)
      End do


c ... creating a virtual superelement
      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)


c ... calculating maximal relocable distance
      heit = 1D12
      Do 10 n = 1, lE
         iE = iEs(n)

         Do i1 = 1, 3
            If(IPE(i1, iE).EQ.iP1) Then
               i2 = ip(i1 + 1) 
               i3 = ip(i2 + 1)

               iP2 = IPE(i2, iE)
               iP3 = IPE(i3, iE)

               v = tri_area( XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
               d = edge_length(XYP(1, iP2), XYP(1, iP3))
               heit = min(heit, 2D0 * dabs(v) / d)

               Goto 10
            End if
         End do
 10   Continue

      
C ... checking for inverted elements
      flagTM = ifXnode(status, ANIUntangleMesh) 
      if(flagTM) Then
         nBad = 0
         Do n = 1, lE
            If(qEs(n).LE.0D0) nBad = nBad + 1
         End do

         If(heit.LE.0D0) Goto 1000
      End if

      flagTM = flagTM .AND. nBad.GT.0


C ... computing the gradient
      xStep = distSP(iP1,  1D0, 0D0, lE, IPEs, XYP)
      tStep = distSP(iP1, -1D0, 0D0, lE, IPEs, XYP)
      If(tStep.GT.xStep) xStep = -tStep
      xStep = ANIDiscreteGrad * min(dabs(xStep), heit) 
     &                        * dsign(1D0, xStep)

      yStep = distSP(iP1, 0D0,  1D0, lE, IPEs, XYP)
      tStep = distSP(iP1, 0D0, -1D0, lE, IPEs, XYP)
      If(tStep.GT.yStep) yStep = -tStep
      yStep = ANIDiscreteGrad * min(dabs(yStep), heit) 
     &                        * dsign(1D0, yStep)


      If(ICPs.EQ.jInode) Then
         XYPs(1, 1) = XYP(1, iP1) + xStep
         XYPs(2, 1) = XYP(2, iP1)
         XYPs(1, 2) = XYP(1, iP1)
         XYPs(2, 2) = XYP(2, iP1) + yStep

         Do i = 1, 2
            iCRVs(i) = 0
            iFNCs(i) = 0
         End do
      Else if(ifXnode(ICPs, jSnode)) Then
         Call infoP(iP1, iPa, iPb, iF1, iF2, par, IPF, Crv, lF, iFs)
         Do i = 1, 2
            XYPs(i, 1) = XYP(i, iP1) +
     &                   ANIDiscreteGrad * (XYP(i, iPa) - XYP(i, iP1))
            XYPs(i, 2) = XYP(i, iP1) +
     &                   ANIDiscreteGrad * (XYP(i, iPb) - XYP(i, iP1))
         End do

         iCRVs(1) = lbC(iF1)
         iCRVs(2) = lbC(iF2)

         iFNCs(1) = lbC(iF1)
         iFNCs(2) = lbC(iF2)

         tc(1) = par(2) + ANIDiscreteGrad * (par(1) - par(2))
         tc(2) = par(3) + ANIDiscreteGrad * (par(4) - par(3))

         Do i = 1, 2
            Call findSE(nCrvFnc, LFnc, iFNCs(i), k)
            If(k.GT.0) Then
               ir = ILt(k)
               Call prjCrv(XYPs(1, i), prjXYPs(1, i), iFNCs(i), tc(i), 
     &                     CrvFunction,
     &                     L1Et(1, ir), L2Et(ir), nL2t(k), tE(ir))
               Do j = 1, 2
                  XYPs(j, i) = prjXYPs(j, i)
               End do
            End if
         End do
      End if


c ... calculating the search direction
      q(1) = 1D0 - rQuality
      Do i = 1, 2
         q(i + 1) = NLnFnc(
C group (F)
     &       2, XYPs(1, i),
C group (ANI)
     &       XYP, HesP, hStar, 
     &       lE, iEs, prjXYPs, IPEs, HesPs, detGs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       iSE, rSE, iP1, 0, CrvFunction,
     &       L1Et, L2Et, tE,
     &       nL2t, nCrvFnc, LFnc, ILt)


C  ...   updating the spoiled values of qEs
         Do n = 1, lE
            qEs(n) = qEu(n)
         End do
      End do


      nU = 2
      hMin = 0D0
      iref = 1

      Do i = 1, 2
         U(i) = XYP(i, iP1)
      End do
      If(ICPs.EQ.jInode) Then
         ZZ(1) = (q(1) - q(2)) / xStep
         ZZ(2) = (q(1) - q(3)) / yStep
         hMax = distSP(iP1, ZZ(1), ZZ(2), lE, IPEs, XYP)

c  ...   check that the direction is not zero
         If(ZZ(1) ** 2 + ZZ(2) ** 2.LE.0D0) Goto 1000
         If(hMax.LE.0D0)                    Goto 1000

      Else if(ifXnode(ICPs, jSnode)) Then
         If(q(2).LE.q(1)) Then
            If(iCRVs(1).NE.0) Then
               nU = 1

               iFmove = iF1
               iPend = iPa
               hMin = min(par(1), par(2))
               hMax = max(par(1), par(2))
               iref = 1

               t1 = par(2)
               t2 = par(1)
               ZZ(1) = t2 - t1

               U(1) = t1
            Else
               Do i = 1, 2
                  ZZ(i) = XYP(i, iPa) - XYP(i, iP1)
               End do
               hMax = dsqrt(ZZ(1) ** 2 + ZZ(2) ** 2)
            End if
         Else If(q(3).LE.q(1)) Then
            If(iCRVs(2).NE.0) Then
               nU = 1

               iFmove = iF2
               iPend = iPb
               hMin = min(par(3), par(4))
               hMax = max(par(3), par(4))
               iref = 2

               t1 = par(3)
               t2 = par(4)
               ZZ(1) = t2 - t1

               U(1) = t1
            Else
               Do i = 1, 2
                  ZZ(i) = XYP(i, iPb) - XYP(i, iP1)
               End do
               hMax = dsqrt(ZZ(1) ** 2 + ZZ(2) ** 2)
            End if
         Else
            Goto 1000
         End if
      End if

      hMax = hMax - 0.2 * (hMax - hMin)
      If(flagTM) hMax = 2 * hMax


C ... minimizing the functional
      qMin = q(1)

      Call minim(
C group(F)
     &     nU, U, ZZ, hMin, hMax, qMin, U1,
     &     icnt, rMove, flag,
C group (ANI)
     &     XYP, IPE, HesP, hStar,
     &     lE, iEs, XYPs, IPEs, detGs, HesPs, qEs,
     &     nPw, nEw, XYPw, HesPw, IPEw, 
     &     MetricFunction, flagAnalytic,
     &     iSE, rSE, iP1, iFNCs(iref), CrvFunction,
     &     L1Et, L2Et, tE,
     &     nL2t, nCrvFnc, LFnc, ILt,
C group (Tangle)
     &     flagTM, nBad)

      If(rMove.EQ.0D0) Goto 1000


C ... analysing information from the previous routine
      iref = 1
      If(.NOT.flag) Then
         If(nU.EQ.1) Goto 1000

         If(q(2).LT.q(1)) Then
            iref = 1
         Else If(q(3).LT.q(1)) Then
            iref = 2
         Else
            Goto 1000
         End if

         qMin = NLnFnc(
C group (F)
     &       2, XYPs(1, iref),
C group (ANI)
     &       XYP, HesP, hStar,
     &       lE, iEs, prjXYPs, IPEs, HesPs, detGs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       iSE, rSE, iP1, 0, CrvFunction,
     &       L1Et, L2Et, tE,
     &       nL2t, nCrvFnc, LFnc, ILt)
      End if


C ... updating the grid
      flag = .TRUE.
      Call pntUpd(iP1, ICP,  XYP,  HesP,  detG,
     &                 ICPs, XYPs(1, iref), HesPs, detGs)

C ... updating for inverted elements
      If(flagTM) Then
         Do n = 1, lE
            If(iEs(n).GT.0) Then
               Call updQb(n, lE, iEs, XYP, IPEs, qEs)
            End if
         End do
      End if

      If(ifXnode(ICPs, jSnode)) Then
         Call findSE(lF, iFs, iF1, nF1)
         Call findSE(lF, iFs, iF2, nF2)

         IPFs(1, nF1) = iPa
         IPFs(2, nF1) = iP1

         IPFs(1, nF2) = iP1
         IPFs(2, nF2) = iPb

         iBNDs = lbF(iF1)

         t1 = par(1)
         t3 = U(1)
         t2 = par(4)

         Call facUpd(nF1, IPF, lbF, Crv, lbC,
     &        iFs, IPFs, iCRVs(1), iFNCs(1), iBNDs, t1, t3)

         Call facUpd(nF2, IPF, lbF, Crv, lbC,
     &        iFs, IPFs, iCRVs(2), iFNCs(2), iBNDs, t3, t2)
      End if

      Do 20 n = 1, lE
         If(iEs(n).LE.0) Goto 20
         Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(n), qEs(n))
 20   Continue

 1000 Return
      End



C ======================================================
      Real*8 Function distSP(iPo, nx, ny, lE, IPEs, XYP)
C ======================================================
      Integer IPEs(3, *)
      Real*8  XYP(2, *), nx, ny

C group (Local variables)
      Integer iref(4)
      Real*8  ox, oy, ax, ay, bx, by, dt, dta, dtb, dtc, dx, dy
      Real*8  cx, cy, aNorm, bNorm, bet, gam
      Logical flag

C ======================================================
      flag = .FALSE.

      iref(1) = 1
      iref(2) = 2
      iref(3) = 3
      iref(4) = 1

      distSP = 1D24
      Do 20 n = 1, lE
         Do i1 = 1, 3
            If(IPEs(i1, n).EQ.iPo) Then
               i2 = iref(i1 + 1)
               i3 = iref(i2 + 1)

               iPa = IPEs(i2, n)
               iPb = IPEs(i3, n)

               ox = XYP(1, iPo)
               oy = XYP(2, iPo)

               ax = XYP(1, iPa) - ox
               ay = XYP(2, iPa) - oy

               bx = XYP(1, iPb) - ox
               by = XYP(2, iPb) - oy
               Goto 10
            End if
         End do
         Goto 20

c  ...  calculating the bisectris between (ax, ay) and (bx, by)
 10      aNorm = dsqrt(ax ** 2 + ay ** 2)
         bNorm = dsqrt(bx ** 2 + by ** 2)
         cx = ax / aNorm + bx / bNorm
         cy = ay / aNorm + by / bNorm

         dt  = ax * cy - cx * ay
         dtc = nx * cy - cx * ny
         If(dtc / dt.LT.0D0) Goto 15

         dta = ax * ny - nx * ay
         If(dta / dt.GE.0D0) Goto 17

 15      dt = -dt
         If(dtc / dt.LT.0D0) Goto 20

         dtb = bx * ny - nx * by
         If(dtb / dt.LT.0D0) Goto 20

 17      Continue
         If(dabs(ny).LE.1D-12) Then
            dy = oy

            If(dabs(ay - by).LE.1D-12) Goto 20
            dx = ((bx + ox) * ay - (ax + ox) * by) / (ay - by)
         Else
            dt  = nx / ny
            bet = bx - ax + dt * (ay - by)
            gam = (by + oy) * ax - (ay  + oy) * bx +
     &            dt * oy * (by - ay)

            dy = -gam / bet
            dx = ox + dt * (dy - oy)
         End if

         distSP = min(distSP, (dx - ox) ** 2 + (dy - oy) ** 2)
         flag = .TRUE.
 20   Continue
      distSP = dsqrt(distSP)

      If(.NOT.flag) distSP = -distSP
      Return
      End
