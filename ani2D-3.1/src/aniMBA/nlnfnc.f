C ======================================================================
      Double Precision Function NLnFnc(
C ======================================================================
C group (F)
     &       nU, U,
C group (ANI)
     &       XYP, HesP, hStar, 
     &       lE, iEs, XYPs, IPEs, HesPs, detGs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw,
     &       MetricFunction, flagAnalytic,
     &       iSE, rSE, iP1, iFNCs, CrvFunction,
     &       L1Et, L2Et, tE,
     &       nL2t, nCrvFnc, LFnc, ILt)
C ======================================================================
C group (F)
      Real*8 U(*)

C group (ANI)
      Integer iEs(*), IPEs(3, *), IPEw(3, *), iSE(*)
      Real*8  XYP(2, *), HesP(3, *), hStar

      Real*8  XYPs(2),   HesPs(*),  detGs,  qEs(*)
      Real*8  XYPw(2, *), HesPw(3, *), rSE(*)
      Real*8  tE(*)

      EXTERNAL CrvFunction

      Integer L1Et(2, *), L2Et(*)
      Integer nL2t(*)
      Integer LFnc(*), ILt(*)

      Logical flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (Local variables)
      Integer  ip(4), info(3)
      Real*8   prjXYPs(2), XYPo(2), tc

C ======================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

c ... control parameters for the interpolation routine
      info(1) = 0
      info(2) = 1


      If(nU.EQ.2) Then
         Do i = 1, 2
            XYPs(i) = U(i)
         End do
      Else if(nU.EQ.1) Then
         tc = U(1)

         Call aniCrv(tc, XYPs, iFNCs, CrvFunction)
      End if

      Call findSE(nCrvFnc, LFnc, iFNCs, k)
      If(k.GT.0) Then
         ir = ILt(k)
         Call prjCrv(XYPs, prjXYPs, iFNCs, tc, CrvFunction,
     &               L1Et(1, ir), L2Et(ir), nL2t(k), tE(ir))

         If(.NOT.flagAnalytic) Then
            Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &                  prjXYPs, HesPs, iSE, rSE, info)
         Else
c           Call scaleBack(prjXYPs, XYPo)
            Call iniQanalytic(1, prjXYPs, MetricFunction, HesPs)
         End if
      Else 
         If(.NOT.flagAnalytic) Then
            Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 1, 
     &                  XYPs, HesPs, iSE, rSE, info)
         Else
c           Call scaleBack(XYPs, XYPo)
            Call iniQanalytic(1, XYPs, MetricFunction, HesPs)
         End if
      End if

      Call calDet(HesPs, detGs)


      Do 50 n = 1, lE
         iE = iEs(n)
         If(iE.LE.0) Goto 50

         Do i1 = 1, 3
            If(IPEs(i1, n).EQ.iP1) Then
               i2 = ip(i1 + 1)
               i3 = ip(i2 + 1)

               iPa = IPEs(i2, n)
               iPb = IPEs(i3, n)

               Call calQE(
     &              HesPs,        XYPs,
     &              HesP(1, iPa), XYP(1, iPa),
     &              HesP(1, iPb), XYP(1, iPb),
     &              hStar, qEs(n))

               Goto 50
            End if
         End do
 50   Continue

      NLnFnc = 1D0
      Do n = 1, lE
         NLnFnc = min(NLnFnc, qEs(n))
      End do
      NLnFnc = 1D0 - NLnFnc
      Return
      End
