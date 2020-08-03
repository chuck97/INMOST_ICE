C*************************************************************
C     MINIMIZATION OF THE NONLINEAR FUNCTIONAL FUNC(U)
C     IN GIVEN DIRECTION  ZZ
C*************************************************************
      Subroutine minim(
C group(F)
     &       nU, U, ZZ, Lamda1, Lamda2, fMin, U1,
     &       icnt, rMove, flagResult,
C group (ANI)
     &       XYP, IPE, HesP, hStar, 
     &       lE, iEs, XYPs, IPEs, detGs, HesPs, qEs,
     &       nPw, nEw, XYPw, HesPw, IPEw, 
     &       MetricFunction, flagAnalytic,
     &       iSE, rSE, iP1, iFNCs, CrvFunction,
     &       L1Et, L2Et, tE,
     &       nL2t, nCrvFnc, LFnc, ILt,
C group (Tangle)
     &       flagTM, nBad)
C ===========================================================
      include 'makS.fd'
C ===========================================================
C group (F)
      Real*8  U(*), ZZ(*), U1(*)
      Real*8  Lamda1, Lamda2, fMin, rMove
      Logical flagResult

C group (ANI)
      Integer IPE(3, *)
      Integer iEs(*), IPEs(3, *), IPEw(3, *), iSE(*)
      Real*8  XYP(2, *),  HesP(3, *), hStar
      
      Real*8  XYPs(2, *), HesPs(*),   detGs, qEs(*)
      Real*8  XYPw(2, *), HesPw(3, *), rSE(*)
      Real*8  tE(*)

      EXTERNAL CrvFunction

      Integer L1Et(2, *), L2Et(*)
      Integer nL2t(*)
      Integer LFnc(*), ILt(*)

      Logical flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (Tangle)
      Logical flagTM
      Integer nBad

C group (Local variables)
      Real*8  NLnFnc
      Real*8  Lamda, dLamda, step, x, f0, f1
      Real*8  aStep

      Real*8  qEt(MaxS), XYPt(2), HesPt(3), detGt
      Real*8  XYPu(2)

      Integer iOs(MaxS)
      Logical flagOrient

C ===========================================================
      flagResult = .FALSE.

      x = 0D0
      Do i = 1, nU
         x = x + ZZ(i) ** 2
      End do
      x = dsqrt(x)

      Do i = 1, nU
         ZZ(i) = ZZ(i) / x
      End do

      x=Lamda1
      Lamda1 = min(x, Lamda2)
      Lamda2 = max(x, Lamda2)

      icnt = 0
      step = 5D-3
      aStep = 0D0
      dLamda = Lamda2 - Lamda1

      f0 = fMin

      Call copySQ(lE, qEs, XYPs, HesPs, detGs,
     &                qEt, XYPt, HesPt, detGt)


 10   Lamda = step  * dLamda
      If(Lamda.GT.(1D0 - aStep - step) * dLamda) Then
         f1 = f0 + 1D0
      Else
         Do i = 1, nU
            U1(i) = U(i) + Lamda * ZZ(i)
         End do


         icnt = icnt + 1
         f1 = NLnFnc(nU, U1,
C group (ANI)
     &        XYP, HesP, hStar,
     &        lE, iEs, XYPs, IPEs, HesPs, detGs, qEs,
     &        nPw, nEw, XYPw, HesPw, IPEw,
     &        MetricFunction, flagAnalytic,
     &        iSE, rSE, iP1, iFNCs, CrvFunction,
     &        L1Et, L2Et, tE,
     &        nL2t, nCrvFnc, LFnc, ILt)


c  ...  checking for inverted elements
         If(flagTM) Then
            Do i = 1, 2
               XYPu(i) = XYP(i, iP1)
               XYP(i, iP1) = XYPs(i, 1)
            End do

            Do n = 1, lE
               If(qEs(n).GT.0D0) Then
                  Call updQb(n, lE, iEs, XYP, IPEs, qEs)
               End if
            End do

            Do i = 1, 2
               XYP(i, iP1) = XYPu(i)
            End do

            mBad = 0
            Do n = 1, lE
               If(qEs(n).LE.0D0) Then
                  mBad = mBad + 1
                  f1 = max(f1, 1D0 - qEs(n))
               End if
            End do
            If(mBad.GT.nBad) f1 = f0 + 1D0

c  ...  checking orientation of triangles
         Else
            Call calSO(XYP, IPE, lE, iEs, iOs)
            Call chkSO(iP1, iP1, XYPs, XYP, IPE, lE, iEs, iOs, 
     &                 flagOrient)
            If(.NOT.flagOrient) f1 = f0 + 1D0
         End if
      End if

      If(f1.LT.f0) Then
         f0 = f1
         Do i = 1, nU
            U(i) = U1(i)
         End do

         flagResult = .TRUE.

         Call copySQ(lE, qEs, XYPs, HesPs, detGs,
     &                   qEt, XYPt, HesPt, detGt)

         aStep = aStep + step
         step = 4D0 * step
      Else
         step = step / 4D0
         If(step.LT.2D-5) Then
            fMin = f0
            rMove = aStep * dLamda
            Goto 1000
         End if
      End if
      Goto 10

 1000 Call copySQ(lE, qEt, XYPt, HesPt, detGt,
     &                qEs, XYPs, HesPs, detGs)

      Return
      End




