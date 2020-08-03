C ======================================================================
      Subroutine insrtPoint(
C ======================================================================
     &           iwF, iwE,
c group (M)
     &           nP, MaxP, XYP,  
     &           nF, MaxF, IPF, lbF, 
     &           Crv, lbC, CrvFunction,
     &           nE, MaxE, IPE, lbE, 
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
     &           lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &                     iFs, iEs, IPFs, IPEs, qEs,
     &           flag)
C ======================================================================
      implicit none
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
      include 'operat.fd'
C ======================================================================
C Routine realizes one of the mesh operations: inserts
C a point at the middle of edge iwF of element iwE.
C
C *** DATA FLOW CHART
C
C  -> check if the edge can be split
C  -> check that 2-arm rule will be preserved
C  -> insert a virtual point 
C  -> project the point onto curved boundary
C
C  -> virtual evaluation of the element quality
C
C  -> check that no traingles were tangled (for curved edge)
C  -> check that "smoothed" element quality has incresed
C
C  -> update the quality of new mesh elements
C  -> update the list of curvilinear edges
C  -> update mesh cross-references for the new elements
C
C ======================================================================
      Integer  iwF, iwE

C group (M)
      Integer  nP, MaxP
      Real*8  XYP(2, *)

      Integer  nF, MaxF, IPF(2, *), lbF(*), lbC(*)
      Real*8   Crv(2, *)
      EXTERNAL CrvFunction

      Integer  nE, MaxE, IPE(3, *), lbE(*)

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
      Integer  lFu, lEu, iFu(*), iEu(*), IPFu(2, *), IPEu(3, *)
      Integer            iFs(*), iEs(*), IPFs(2, *), IPEs(3, *)
      Real*8   qEu(*), qEs(*)
      Logical flag


C LOCAL VARIABLES
      Real*8   prjXYPs(2), XYPs(2), HesPs(3), detGs, XYPt(2)

      Integer  ip(4), iPs(MaxS), iEx(MaxS), info(3)
      Real*8   tri_area2, vol1, vol2, par(6), t1, t2, tc
      Logical  flagBNDs, ifXnode

      Integer  minClr, tri_orient
      EXTERNAL minClr, tri_orient

      Integer  iP1,iP2, iPa, iPb,iPc,iPd, iF,iF1,iF2, iE1,iE2
      Integer  ICPs, iFNCs, iBNDs, iCRVs
      Integer  i,k,n, i1,i2,i3,i4, nF1, nE1,nE2, lP, lF, lE, lEadd, ir

C ======================================================================
      flag = .FALSE.

      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      i1 = iwF
      i2 = ip(i1 + 1)

      iP1 = IPE(i1, iwE)
      iP2 = IPE(i2, iwE)


c ... checking for edges which can not be split
      ICPs = minClr(ICP(iP1), ICP(iP2))
      If(ifXnode(ICPs, jTnode)) Goto 1000

      iE1 = iwE
      iE2 = IEE(iwF, iE1)

      i3 = ip(i2 + 1)
      iPa = IPE(i3, iwE)

      If(iE2.NE.0) Then
         Do i = 1, 3
            If(IPE(i, iE2).NE.iP1 .AND. IPE(i, iE2).NE.iP2) Then
               i4 = i
               iPb = IPE(i, iE2)
            End if
         End do
      End if


C ... checking for surrounding points
      iF = IFE(iwF, iwE)
      If(iF.EQ.0) Then
         ICPs = jInode
      Else
         lP = 2
         iPs(1) = iP1
         iPs(2) = iP2

         If(iE2.NE.0) Then
            ICPs = jInode + jSnode
         Else
            ICPs = jBnode + jSnode
         End if
      End if


C ... skipping the inverted elements
      If(ifXnode(status, ANIUntangleMesh)) Then
         If(qE(iE1).LE.0D0) Goto 1000
         If(iE2.GT.0) Then
            If(qE(iE2).LE.0D0) Goto 1000
         End if
      End if


c ... control parameyters for the interpolation routine
      info(1) = 0
      info(2) = 1


c ... making the point for inserting
      iFNCs = 0
      Do i = 1, 3
         HesPs(i) = (HesP(i, iP1) + HesP(i, iP2)) / 2
      End do

      If(ifXnode(ICPs, jSnode)) Then
         Call infoF(iF, iP1, iP2, iF1, iF2, iPc, iPd,
     &              par, IPF, lbF, Crv, lF, iFs)

         iCRVs = lbC(iF)
         If(iCRVs.NE.0) Then
            t1 = par(3)
            t2 = par(4)
            tc = (t1 + t2) / 2

            iFNCs = iCRVs
            Call aniCrv(tc, XYPs, iFNCs, CrvFunction)

            Call findSE(nCrvFnc, LFnc, iFNCs, k)
            ir = ILt(k)
            Call prjCrv(XYPs, prjXYPs, iFNCs, tc, CrvFunction,
     &                  L1Et(1, ir), L2Et(ir), nL2t(k), tE(ir))

            If( .NOT.flagAnalytic ) Then
               Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1,
     &                     prjXYPs, HesPs, iSE, rSE, info)
            Else
               Call iniQanalytic(1, prjXYPs, MetricFunction, HesPs)
            End if

            Do i = 1, 2
               XYPt(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do
         Else
            Do i = 1, 2
               XYPs(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do

            If( .NOT.flagAnalytic ) Then
               Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1,
     &                     XYPs, HesPs, iSE, rSE, info)
            Else
               Call iniQanalytic(1, XYPs, MetricFunction, HesPs)
            End if
         End if
      Else
         iCRVs = 0
         Do i = 1, 2
            XYPs(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
         End do

         If( .NOT.flagAnalytic ) Then
            Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1,
     &                  XYPs, HesPs, iSE, rSE, info)
         Else
            Call iniQanalytic(1, XYPs, MetricFunction, HesPs)
         End if
      End if


      Call calDet(HesPs, detGs)


C ... making a virtual evaluation of the quality
      lEadd = 0
      Call findSE(lE, iEs, iE1, nE1)
      Call calQE(
     &     HesP(1, iP1), XYP(1, iP1),
     &     HesP(1, iPa), XYP(1, iPa),
     &     HesPs,        XYPs,
     &     hStar, qEs(nE1))

      If(qEs(nE1).LE.rQuality) Goto 1000

      lE = lE + 1
      Call calQE(
     &     HesP(1, iP2), XYP(1, iP2),
     &     HesP(1, iPa), XYP(1, iPa),
     &     HesPs,        XYPs,
     &     hStar, qEs(lE))

      If(qEs(lE).LE.rQuality) Goto 1000


C  ...  checking for the orientation
      If(iCRVs.NE.0) Then
c        vol1 = tri_area2(XYP(1, iPa), XYP(1, iP1), XYPs)
c        vol2 = tri_area2(XYP(1, iPa), XYP(1, iP1), XYPt)
c        If(vol1 * vol2.LE.0D0) Goto 1000

c        vol1 = tri_area2(XYP(1, iPa), XYP(1, iP2), XYPs)
c        vol2 = tri_area2(XYP(1, iPa), XYP(1, iP2), XYPt)
c        If(vol1 * vol2.LE.0D0) Goto 1000

         i1 = tri_orient(XYP(1, iPa), XYP(1, iP1), XYPs)
         i2 = tri_orient(XYP(1, iPa), XYP(1, iP1), XYPt)
         If(i1 * i2.LE.0) Goto 1000

         i1 = tri_orient(XYP(1, iPa), XYP(1, iP2), XYPs)
         i2 = tri_orient(XYP(1, iPa), XYP(1, iP2), XYPt)
         If(i1 * i2.LE.0) Goto 1000
      End if


      If(iE2.NE.0) Then
         lEadd = 1
         Call findSE(lE, iEs, iE2, nE2)
         Call calQE(
     &        HesP(1, iP1), XYP(1, iP1),
     &        HesP(1, iPb), XYP(1, iPb),
     &        HesPs,        XYPs,
     &        hStar, qEs(nE2))

         If(qEs(nE2).LE.rQuality) Goto 1000

         lE = lE + 1
         Call calQE(
     &        HesP(1, iP2), XYP(1, iP2),
     &        HesP(1, iPb), XYP(1, iPb),
     &        HesPs,        XYPs,
     &        hStar, qEs(lE))

         If(qEs(lE).LE.rQuality) Goto 1000


C  ...  checking for the orientation
         If(iCRVs.NE.0) Then
c           vol1 = tri_area2(XYP(1, iPb), XYP(1, iP1), XYPs)
c           vol2 = tri_area2(XYP(1, iPb), XYP(1, iP1), XYPt)
c           If(vol1 * vol2.LE.0D0) Goto 1000
c
c           vol1 = tri_area2(XYP(1, iPb), XYP(1, iP2), XYPs)
c           vol2 = tri_area2(XYP(1, iPb), XYP(1, iP2), XYPt)
c           If(vol1 * vol2.LE.0D0) Goto 1000

            i1 = tri_orient(XYP(1, iPb), XYP(1, iP1), XYPs)
            i2 = tri_orient(XYP(1, iPb), XYP(1, iP1), XYPt)
            If(i1 * i2.LE.0) Goto 1000

            i1 = tri_orient(XYP(1, iPb), XYP(1, iP2), XYPs)
            i2 = tri_orient(XYP(1, iPb), XYP(1, iP2), XYPt)
            If(i1 * i2.LE.0) Goto 1000
         End if
      End if


C ... checking for the triangles orientation
      If(iF.NE.0 .AND. iCRVs.NE.0) Then
c        vol1 = tri_area2(XYP(1, iPa), XYPs, XYP(1, iP2))
c        vol2 = tri_area2(XYP(1, iPa), XYPs, XYP(1, iP1))
c        If(vol1 * vol2.GE.0D0) Goto 1000

         i1 = tri_orient(XYP(1, iPa), XYPs, XYP(1, iP2))
         i2 = tri_orient(XYP(1, iPa), XYPs, XYP(1, iP1))
         If(i1 * i2.GE.0) Goto 1000

         If(iE2.NE.0) Then
c           vol1 = tri_area2(XYP(1, iPb), XYPs, XYP(1, iP2))
c           vol2 = tri_area2(XYP(1, iPb), XYPs, XYP(1, iP1))
c           If(vol1 * vol2.GE.0D0) Goto 1000

            i1 = tri_orient(XYP(1, iPb), XYPs, XYP(1, iP2))
            i2 = tri_orient(XYP(1, iPb), XYPs, XYP(1, iP1))
            If(i1 * i2.GE.0) Goto 1000
         End if
      End if


c ... final check (will be moved up to the driver)
      If(nP.EQ.MaxP)               Goto 9100
      If(nF.EQ.MaxF .AND. iF.NE.0) Goto 9200
      If(nE.GE.MaxE-lEadd)         Goto 9300


C ... updating the grid
      flag = .TRUE.

      Call pntAdd(iPc, nP, ICP,  XYP,  HesP,  detG, IHolP,
     &                           ICPs, XYPs, HesPs, detGs)

      IPEs(1, nE1) = iP1
      IPEs(2, nE1) = iPc
      IPEs(3, nE1) = iPa

      IPEs(1, lE - lEadd) = iP2
      IPEs(2, lE - lEadd) = iPc
      IPEs(3, lE - lEadd) = iPa


      If(iF.NE.0) Then
         Call facAdd(iF2, nF, MaxF, IHolF)

         Call findSE(lF, iFs, iF, nF1)
         IPFs(1, nF1) = iP1
         IPFs(2, nF1) = iPc

         lF = lF + 1
         iFs(lF) = iF2
         IPFs(1, lF) = iPc
         IPFs(2, lF) = iP2

         iCRVs = lbC(iF)
         iFNCs = iCRVs
         iBNDs = lbF(iF)
         Call facUpd(nF1, IPF, lbF, Crv, lbC,
     &               iFs, IPFs, iCRVs, iFNCs, iBNDs, par(3), tc)

         If(iCRVs.NE.0) iCRVs = iF2
         Call facUpd(lF,  IPF, lbF, Crv, lbC,
     &               iFs, IPFs, iCRVs, iFNCs, iBNDs, tc, par(4))
      End if


      Call eleDel(iEs(nE1), IPE, IEE)

      If(iE2.NE.0) Then
         IPEs(1, nE2) = iP1
         IPEs(2, nE2) = iPc
         IPEs(3, nE2) = iPb

         IPEs(1, lE) = iP2
         IPEs(2, lE) = iPc
         IPEs(3, lE) = iPb

         Call eleDel(iEs(nE2), IPE, IEE)
      End if

      Do n = lE - lEadd, lE
         Call lstAdd(nE, L1E, nL2, L2E, nStep, IHolE,
     &               qE, qEs(n), iEs(n))
         Call eleDel(iEs(n), IPE, IEE)
      End do

      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE1), qEs(nE1))
      If(iE2.NE.0) Then
         Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE2), qEs(nE2))
      End if

      Call eleUpd(nE1, IEP, IPE, IFE, IEE,
     &            lF, lE, iFs, iEs, IPFs, IPEs, iEx)
      Call eleUpd(lE, IEP, IPE, IFE, IEE,
     &            lF, lE, iFs, iEs, IPFs, IPEs, iEx)

      lbE(iEs(lE - lEadd)) = lbE(iE1)

      If(iE2.NE.0) Then
         Call eleUpd(nE2, IEP, IPE, IFE, IEE,
     &               lF, lE, iFs, iEs, IPFs, IPEs, iEx)
         Call eleUpd(lE - 1, IEP, IPE, IFE, IEE,
     &               lF, lE, iFs, iEs, IPFs, IPEs, iEx)

         lbE(iEs(lE)) = lbE(iE2)
      End if

 1000 Continue
      Return

 9100 Continue
      Call errMesMBA(1003, 'insrtPoint', 'variable MaxP is too small')

 9200 Continue
      Call errMesMBA(1004, 'insrtPoint', 'variable MaxF is too small')

 9300 Continue
      Call errMesMBA(1006, 'insrtPoint', 'variable MaxE is too small')
      Return
      End


