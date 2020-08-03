C ======================================================================
      Subroutine swapEdge(
C ======================================================================
     &           iwF, iwE,
c group (M)
     &           XYP, nE, IPE,  
c group (M-EXT)
     &           ICP, IEP, IFE, IEE,
C group (Q)
     &           hStar, rQuality, HesP, 
     &           L1E, L2E, qE, nL2, nStep,
c group (CONTROL)
     &           status, 
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
      Integer  iwF, iwE

C group (M)
      Integer  nE, IPE(3, *)
      Real*8  XYP(2, *)

C group (M-EXT)
      Integer  ICP(*), IEP(*), IFE(3, *), IEE(3, *)

C group (Q)
      Real*8   hStar, rQuality, HesP(3, *), qE(*) 
      Integer  L1E(2, *), L2E(*), nL2, nStep(4)

c group (CONTROL)
      Integer  status

C group (W)
      Integer  lFu, lEu, iFu(*), iEu(*), IPFu(2, *), IPEu(3, *)
      Integer            iFs(*), iEs(*), IPFs(2, *), IPEs(3, *)
      Real*8   qEu(*), qEs(*)
      Logical  flag

C LOCAL VARIABLES
      Integer  ip(4), iPs(MaxS), iEx(MaxS)
      Real*8   t1, t2, tri_area2
      Logical  flagBNDs, ifXnode

      Integer  minClr, tri_orient
      EXTERNAL minClr, tri_orient

      Integer  lP, lF, lE, iP1,iP2, iPa,iPb, iF, iE1,iE2, nE1, nE2
      Integer  i, i1,i2,i3,i4, ICPs

C ======================================================================
      flag = .FALSE.

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      iF = IFE(iwF, iwE)

      i1 = iwF
      i2 = ip(i1 + 1)

      iP1 = IPE(i1, iwE)
      iP2 = IPE(i2, iwE)

      iE1 = iwE
      iE2 = IEE(iwF, iE1)


c ... checking the case when swapping is impossible
      If(iE2.EQ.0) Goto 1000
      If(iF.NE.0)  Goto 1000

      ICPs = minClr(ICP(iP1), ICP(iP2))
      If(ifXnode(ICPs, jTnode)) Goto 1000


      i3 = ip(i2 + 1)
      iPa = IPE(i3, iwE)

      Do i = 1, 3
         If(IPE(i, iE2).NE.iP1 .AND. IPE(i, iE2).NE.iP2) Then
            i4 = i
            iPb = IPE(i, iE2)
         End if
      End do


C ... checking for boundary triangles
      If(ifXnode(status, ANIForbidBoundaryElements)) Then
         If(ifXnode(ICP(iPa), jInode)) Goto 100
         If(ifXnode(ICP(iPb), jInode)) Goto 100

         If(ifXnode(ICP(iP1), jBnode)) Goto 1000
         If(ifXnode(ICP(iP2), jBnode)) Goto 1000

 100     Continue
      End if


c ... checking the case when swapping is impossible
c     t1 = tri_area2(XYP(1, iPa), XYP(1, iPb), XYP(1, iP1))
c     t2 = tri_area2(XYP(1, iPa), XYP(1, iPb), XYP(1, iP2))
c     If(t1 * t2.GE.0D0) Goto 1000

      i1 = tri_orient(XYP(1, iPa), XYP(1, iPb), XYP(1, iP1))
      i2 = tri_orient(XYP(1, iPa), XYP(1, iPb), XYP(1, iP2))
      If(i1 * i2.GE.0) Goto 1000


C ... skipping the inverted elements
      if(ifXnode(status, ANIUntangleMesh)) Then
         If(qE(iE1).LE.0D0) Goto 1000
         If(qE(iE2).LE.0D0) Goto 1000
      End if


c ... creating a virtual superelement
      Call copySE(lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            lF,  lE,  iFs, iEs, IPFs, IPEs, qEs)


C ... making a virtual evaluation of the quality
      Call findSE(lE, iEs, iE1, nE1)
      Call calQE(
     &     HesP(1, iP1), XYP(1, iP1),
     &     HesP(1, iPa), XYP(1, iPa),
     &     HesP(1, iPb), XYP(1, iPb),
     &     hStar, qEs(nE1))

      If(qEs(nE1).LE.rQuality) Goto 1000

      Call findSE(lE, iEs, iE2, nE2)
      Call calQE(
     &     HesP(1, iP2), XYP(1, iP2),
     &     HesP(1, iPa), XYP(1, iPa),
     &     HesP(1, iPb), XYP(1, iPb),
     &     hStar, qEs(nE2))

      If(qEs(nE2).LE.rQuality) Goto 1000


C ... updating the grid
      flag = .TRUE.

      IPEs(1, nE1) = iP1
      IPEs(2, nE1) = iPa
      IPEs(3, nE1) = iPb

      IPEs(1, nE2) = iP2
      IPEs(2, nE2) = iPa
      IPEs(3, nE2) = iPb


      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE1), qEs(nE1))
      Call lstUpd(nE, L1E, nL2, L2E, nStep, qE, iEs(nE2), qEs(nE2))

      Call eleDel(iEs(nE1), IPE, IEE)
      Call eleDel(iEs(nE2), IPE, IEE)

      Call eleUpd(nE1, IEP, IPE, IFE, IEE,
     &            lF, lE, iFs, iEs, IPFs, IPEs, iEx)
      Call eleUpd(nE2, IEP, IPE, IFE, IEE,
     &            lF, lE, iFs, iEs, IPFs, IPEs, iEx)

 1000 Return
      End


