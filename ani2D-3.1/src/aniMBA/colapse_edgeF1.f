C ======================================================================
      Subroutine clpsF1(
C ======================================================================
     &           iwF, iwE, w1, w2,
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
C Routine realizes one of the mesh operations: collapses
C edge iwF of element iwE to its middle point.
C
C *** DATA FLOW CHART
C
C  -> simple check if the edge can be collapsed
C  -> collect information about the edge
c
c  -> define a point to collapse the edge
c  ----> mid-point of the edge 
c  ----> projection onto the curved boundary
c  ----> one of the terminal points of the edge 
c
c  -> virtual evaluation of the element quality
c
c  -> check that no boundary triangles were created
c  -> check that 2-arm rule helds
c  -> check that no traingles were tangled
c
c  -> update the quality of new mesh elements
c  -> update the list of curvilinear edges
c  -> update mesh cross-references for the new elements  
C
C ======================================================================
      Integer  iwF, iwE
      Real*8   w1, w2

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
      Integer  lFu, lEu, iFu(*), iEu(*), IPFu(2, *), IPEu(3, *)
      Integer            iFs(*), iEs(*), IPFs(2, *), IPEs(3, *)
      Real*8   qEu(*), qEs(*)
      Logical  flag


C LOCAL VARIABLES
      Real*8   prjXYPs(2), XYPs(2), HesPs(3), detGs

      Integer  ip(4), iPs(MaxS), iOs(MaxS), iEx(MaxS), info(3)
      Real*8   par(6), t1, t2, tc
      Logical  flagBNDs, flagOrient, flagTM
      
      Logical  chkTangled, ifXnode
      Integer  minClr, maxClr
      EXTERNAL minClr, maxClr

      Integer  iPa,iPb,iPc,iPd, iPt,iP1,iP2,iP3, iF1,iF2, iE,iEt,iE1,iE2
      Integer  lP,lF, lE,kE, iF, ICP1,ICP2,ICP3, ICPt, ICPs
      Integer  iBNDs, iCRVs, iFNCs
      Integer  i,k,n, ir, i1,i2,i3, nFd, nF1,nF2, nBad,mBad

C ======================================================================
      flag = .FALSE.

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      i1 = iwF
      i2 = ip(i1 + 1)
      i3 = ip(i2 + 1)

      iP1 = IPE(i1, iwE)
      iP2 = IPE(i2, iwE)
      iP3 = IPE(i3, iwE)

      iF = IFE(iwF, iwE)

      iE1 = iwE
      iE2 = IEE(iwF, iE1)

C ... checking for faces which can not be collapsed
      ICP1 = ICP(iP1)
      ICP2 = ICP(iP2)
      ICP3 = ICP(iP3)

      ICPs = minClr(ICP1, ICP2)
      If(ifXnode(ICPs, jVnode)) Goto 1000
      If(ifXnode(ICPs, jTnode)) Goto 1000

      If(ifXnode(ICPs, jSnode) .AND. iF.EQ.0) Goto 1000

      ICPt = maxClr(ICP1, ICP2)
      If(ifXnode(ICPt, jVnode) .AND. ifXnode(ICPt, jTnode)) Goto 1000


C ... gathering information about the face
      Do i = 1, 6
         par(i) = 0D0
      End do

      If(iF.NE.0) Then
         Call infoF(iF, iP1, iP2, iF1, iF2, iPc, iPd,
     &              par, IPF, lbF, Crv, lFu, iFu)

         If(iPc.LE.0 .OR. iPd.LE.0) Goto 1000
      End if


C ... check for inverted elements
      flagTM = ifXnode(status, ANIUntangleMesh)

      if(flagTM) Then
         nBad = 0
         Do n = 1, lEu
            If(qEu(n).LE.0D0) nBad = nBad + 1
         End do
      End if

      flagTM = flagTM .AND. nBad.GT.0
      

c ... creating a virtual superelement
      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)


c ... control parameyters for the interpolation routine
      info(1) = 0
      info(2) = 1


C ... find a point to which we collapse the edge
      iFNCs = 0
      If(ICP1.EQ.jInode .AND. ICP2.EQ.jInode) Then
         ICPs = jInode
         Do i = 1, 2
            XYPs(i) = XYP(i, iP1) * w1 + XYP(i, iP2) * w2
         End do

         Do i = 1, 3
            HesPs(i) = HesP(i, iP1) * w1 + HesP(i, iP2) * w2
         End do

         If(.NOT.flagAnalytic) Then
            Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &                  XYPs, HesPs, iSE, rSE, info)
         Else
            Call iniQanalytic(1, XYPs, MetricFunction, HesPs)
         End if

         Call calDet(HesPs, detGs)

      Else If(ifXnode(ICP1, jSnode) .AND. ICP2.EQ.jInode .OR.
     &        ifXnode(ICP1, jTnode) .AND. ICP2.EQ.jInode .OR.
     &        ifXnode(ICP1, jTnode) .AND. ifXnode(ICPs, jBnode) .OR.
     &        ifXnode(ICP1, jVnode)) Then
         ICPs = ICP1
         Do i = 1, 2
            XYPs(i) = XYP(i, iP1)
         End do

         Do i = 1, 3
            HesPs(i) = HesP(i, iP1)
         End do
         detGs = detG(iP1)

         t1 = par(2)
         t2 = par(3)

      Else If((ifXnode(ICP2, jSnode) .AND. ICP1.EQ.jInode) .OR.
     &        (ifXnode(ICP2, jTnode) .AND. ICP1.EQ.jInode) .OR.
     &        (ifXnode(ICP2, jTnode) .AND. ifXnode(ICPs, jBnode)) .OR.
     &        (ifXnode(ICP2, jVnode))) Then
         ICPs = ICP2
         Do i = 1, 2
            XYPs(i) = XYP(i, iP2)
         End do

         Do i = 1, 3
            HesPs(i) = HesP(i, iP2)
         End do
         detGs = detG(iP2)

         t1 = par(4)
         t2 = par(5)

C ...    changing order of points to be consistent with the previous case
         i = iP2
         iP2 = iP1
         iP1 = i

      Else If(ifXnode(ICPs, jSnode)) Then
c!       ICPs = minClr(ICP1, ICP2)
         iCRVs = lbC(iF)
         If(iCRVs.NE.0) Then
            t1 = par(3)
            t2 = par(4)
            tc = (t1 + t2) / 2

            iFNCs = iCRVs
            Call aniCrv(tc, XYPs, iFNCs, CrvFunction)

            t1 = tc
            t2 = tc

            Call findSE(nCrvFnc, LFnc, iFNCs, k)
            ir = ILt(k)
            Call prjCrv(XYPs, prjXYPs, iFNCs, tc, CrvFunction,
     &                  L1Et(1, ir), L2Et(ir), nL2t(k), tE(ir))
  
            If(.NOT.flagAnalytic) Then
               Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &                     prjXYPs, HesPs, iSE, rSE, info)
            Else
               Call iniQanalytic(1, prjXYPs, MetricFunction, HesPs)
            End if
         Else
            Do i = 1, 2
               XYPs(i) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do

            If(.NOT.flagAnalytic) Then
               Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &                     XYPs, HesPs, iSE, rSE, info)
            Else
               Call iniQanalytic(1, XYPs, MetricFunction, HesPs)
            End if
         End if

         Call calDet(HesPs, detGs)
      Else
         Goto 1000
      End if


C ... virtual evaluation of the superelement quality
      kE = 0
      Do 10 n = 1, lE
         iE = iEs(n)
         If(iE.EQ.iE1 .OR. iE.EQ.iE2) Then
            iEs(n) = -iEs(n)
            Goto 10
         End if

         Do i1 = 1, 3
            If(IPEs(i1, n).EQ.iP1 .OR. IPEs(i1, n).EQ.iP2) Then
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)

               IPEs(i1, n) = iP1
               iPa = IPEs(i2, n)
               iPb = IPEs(i3, n)

               Call calQE(
     &              HesP(1, iPa), XYP(1, iPa),
     &              HesP(1, iPb), XYP(1, iPb),
     &              HesPs,        XYPs,
     &              hStar, qEs(n))

               If(qEs(n).LE.rQuality) Goto 1000

               kE = kE + 1
               Goto 10
            End if
         End do

         iEs(n) = 0
 10   Continue


C ... check for number of elements in a superelement
c     If(kE.GE.9) Goto 1000


C ... check for boundary triangles
      If(ifXnode(status, ANIForbidBoundaryElements)) Then
         Do 20 n = 1, lE
            iE = iEs(n)
            If(iE.LE.0) Goto 20

            Do i = 1, 3
               iPt = IPEs(i, n)
               If(ifXnode(ICP(iPt), jInode)) Goto 20
            End do

            Goto 1000
 20      Continue
      End if


C ... checking for orientation of triangles
      Call calSO(XYP, IPE, lE, iEs, iOs)
      Call chkSO(iP1, iP2, XYPs, XYP, IPE, lE, iEs, iOs, flagOrient)
      If(.NOT.flagOrient) Goto 1000


C ... checking for inverted elements
      If(flagTM) Then
         Do n = 1, lE
            If(iEs(n).GE.0 .AND. qEs(n).GT.0D0) Then
               Call updQb(n, lE, iEs, XYP, IPEs, qEs)
            End if
         End do

         mBad = 0
         Do n = 1, lE
            If(iEs(n).GE.0 .AND. qEs(n).LE.0D0) mBad = mBad + 1
         End do

         If(mBad.GE.nBad) Goto 1000

C  ...  colapsing may result in topologically wrong mesh
        flagTM = chkTangled(lE, iEs, IPEs)
        If(flagTM) Goto 1000
      End if


C ... updating the grid
      flag = .TRUE.

      Call pntUpd(iP1, ICP,  XYP,  HesP,  detG,
     &                 ICPs, XYPs, HesPs, detGs)
      Call pntDel(iP2, nP, ICP, IHolP)


      If(ifXnode(ICP1, jSnode) .AND. ifXnode(ICP2, jSnode)) Then
         Call findSE(lF, iFs, iF, nFd)
         Call facDel(iF, nF, IPF, lbC, IHolF)
         iFs(nFd) = -iFs(nFd)

         Call findSE(lF, iFs, iF1, nF1)
         Call findSE(lF, iFs, iF2, nF2)

         IPFs(1, nF1) = iPc
         IPFs(2, nF1) = iP1

         IPFs(1, nF2) = iP1
         IPFs(2, nF2) = iPd

         iCRVs = lbC(iF1)
         iFNCs = iCRVs
         iBNDs = lbF(iF1)

         Call facUpd(nF1, IPF, lbF, Crv, lbC,
     &               iFs, IPFs, iCRVs, iFNCs, iBNDs, par(1), t1)

         iCRVs = lbC(iF2)
         iFNCs = iCRVs
         iBNDs = lbF(iF2)
         Call facUpd(nF2, IPF, lbF, Crv, lbC,
     &               iFs, IPFs, iCRVs, iFNCs, iBNDs, t2, par(6))
      End if


      Do n = 1, lE
         iEt = iEs(n)
         If(iEt.LT.0) Then
            iEt = -iEt
            Call lstDel(nE, L1E, nL2, L2E, nStep, IHolE, qE, iEt)
            Call eleDel(iEt, IPE, IEE)
         Else If(iEt.GT.0) Then
            Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEt, qEs(n))

            Call eleUpd(n, IEP, IPE, IFE, IEE,
     &                  lF, lE, iFs, iEs, IPFs, IPEs, iEx)
         End if
      End do

1000  Return
      End





