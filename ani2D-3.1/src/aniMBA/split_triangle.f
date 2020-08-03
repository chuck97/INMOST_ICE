C ======================================================================
      Subroutine splitTriangle(iwE,
C ======================================================================
c group (M)
     &                  nP, MaxP, XYP, 
     &                  nE, MaxE, IPE, lbE, 
c group (M-EXT)
     &            ICP, IEP, IFE, IEE,
     &            IHolP, IHolE,
C group (Q)
     &            HesP, detG, qE,
C group (W)
     &                  lFu, lEu, iFu, iEu, IPFu, IPEu, qEu,
     &            flag)
C ======================================================================
      implicit none
      include 'makS.fd'
      include 'colors.fd'
      include 'operat.fd'
C ======================================================================
C Routine splits triangle into 3 triangles by inserting
C one interior point at the center of mass.
C
C *** Remarks:
C        1. The mesh quality may drop down after the splitting.
C ======================================================================
      Integer  iwE

C group (M)
      Integer  nP, MaxP
      Real*8  XYP(2, *)

      Integer  nE, MaxE, IPE(3, *), lbE(*)

C group (M-EXT)
      Integer  ICP(*), IEP(*), IFE(3, *), IEE(3, *)
      Integer  IHolP(*), IHolE(*)

C group (Q)
      Real*8  HesP(3, *), detG(*), qE(*)

C group (W)
      Integer  lFu, lEu, iFu(*), iEu(*), IPFu(2, *), IPEu(3, *)
      Real*8   qEu(*)
      Logical flag


C LOCAL VARIABLES
      Integer  ip(4), minClr, iEx(MaxS)
      Real*8   XYPs(2), HesPs(3), detGs
      Logical ifXnode

      Integer  i,n, i1,i2, iPs, iP1,iP2,iP3, ICPs, lF,lE,kE, nE1
      Integer  lEold, lEadd

C ======================================================================
      flag = .FALSE.

      lF = lFu
      lE = lEu

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      Call findSE(lE, iEu, iwE, nE1)

      iP1 = IPE(1, iwE)
      iP2 = IPE(2, iwE)
      iP3 = IPE(3, iwE)


c ... check for faces which can not be split
      ICPs = minClr(ICP(iP1), ICP(iP2))
      ICPs = minClr(ICPs, ICP(iP3))
      If(ifXnode(ICPs, jTnode)) Goto 1000


c ... create the interior point
      ICPs = jInode
      Do i = 1, 2
         XYPs(i) = (XYP(i, iP1) + XYP(i, iP2) + XYP(i, iP3)) / 3
      End do

      Do i = 1, 3
         HesPs(i) = (HesP(i, iP1) + HesP(i, iP2) + HesP(i, iP3)) / 3
      End do

      Call calDet(HesPs, detGs)


c ... final check (will be moved in the driver)
      If(nP.EQ.MaxP)   Goto 9100
      If(nE+3.GE.MaxE) Goto 9300

      Call pntAdd(iPs, nP, ICP,  XYP,  HesP,  detG, IHolP,
     &                           ICPs, XYPs, HesPs, detGs)

c ... create 3 elements
      lEold = lE
      Do i1 = 1, 3
         i2 = ip(i1 + 1)

         iP1 = IPE(i1, iwE)
         iP2 = IPE(i2, iwE)

         If(i1.EQ.1) Then
            kE = nE1
         Else
            lE = lE + 1
            If(lE.GT.MaxS) Goto 9000
            kE = lE
         End if

c        Call calQE(
c    &        HesP(1, iP1), detG(iP1), XYP(1, iP1),
c    &        HesP(1, iPa), detG(iPa), XYP(1, iPa),
c    &        HesPs, detGs, XYPs,
c    &        hStar, qEu(kE))

         IPEu(1, kE) = iP1 
         IPEu(2, kE) = iP2
         IPEu(3, kE) = iPs 
      End do


C ... update the grid
      flag = .TRUE.

C!!!  next line simulates lstAdd
      qE(iEu(nE1)) = qEu(nE1)
      Call eleDel(iEu(nE1), IPE, IEE)

      Do n = lEold + 1, lE
C!!!     3 next lines simulate lstAdd
         nE = nE + 1
         iEu(n) = nE 
         qE(iEu(n)) = qEu(n)
         Call eleDel(iEu(n), IPE, IEE)
      End do

      Call eleUpd(nE1, IEP, IPE, IFE, IEE,
     &            lF, lE, iFu, iEu, IPFu, IPEu, iEx)

      Do n = lEold + 1, lE
         Call eleUpd(n, IEP, IPE, IFE, IEE,
     &               lF, lE, iFu, iEu, IPFu, IPEu, iEx)

         lbE(iEu(n)) = lbE(iwE)
      End do

 1000 Continue
      Return

 9000 Continue
      Call errMesMBA(1007, 'splitTriangle', 'parameter MaxS is small')

 9100 Continue
      Call errMesMBA(1003, 'splitTriangle', 'variable MaxP is small')

 9300 Continue
      Call errMesMBA(1006, 'splitTriangle', 'variable MaxE is small')
      End


