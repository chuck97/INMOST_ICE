C ======================================================================
      Subroutine mbaFixShape(
C ======================================================================
c group (M)
     &      nP, nPfix, MaxP, XYP, lbP, fixP,
     &      nF, nFfix, MaxF, IPF, lbF, fixF,
     &      nC,              Crv, lbC, CrvFunction,
     &      nE, nEfix, MaxE, IPE, lbE, fixE,
c group (CONTROL)
     &      nEStar, Quality, control, MetricFunction,
c group (W)
     &      MaxWr, MaxWi, rW, iW, iERR)
C ======================================================================
      implicit none
      include 'lintrp.fd'
      include 'status.fd'
C ======================================================================
C  Routine is similar to mbaAnalytic but the element quality is 
C  defined only by shape regularity (not size!) of the element.
C ======================================================================
C  VARIABLES & PARAMETER are described in mba_nodal.f
C
C  MetricFunction - integer function created by the user (see 
C                   example in file forlibmba.f)
C
C    Integer Function MetricFunction(x, y, z, Metric)
C
C  This routine creates a metric at the given point (x,y, z). The
C  metric is a 2x2 positive definite symmetric tensor:
C
C                M11   M12
C      Metric =     
C                M12   M22
C
C  Only the upper triangular part of array Metric must be defined.
C
C
C *** Authors: K. Lipnikov (lipnikov@gmail.com)
C              Y. Vassilevski (yuri.vasilevski@gmail.com)
C ======================================================================
C group (M)
      Integer  nP, nPfix, MaxP, lbP(*), fixP(*)
      Real*8  XYP(2, *)

      Integer  nF, nFfix, MaxF, IPF(2, *), lbF(*), fixF(*)
      
      Integer  nC, lbC(*)
      Real*8   Crv(2, *)
      EXTERNAL CrvFunction

      Integer  nE, nEfix, MaxE, IPE(3, *), lbE(*), fixE(*)

C group (CONTROL)
      Integer  nEStar, control(*), MetricFunction
      Real*8   Quality
      EXTERNAL MetricFunction

C group (W)
      Integer  MaxWr, MaxWi, iW(*), iERR
      Real*8  rW(*)


C LOCAL VARIABLES
      Integer  MaxSkipE, MaxQItr, status, nQItr, iPrint, iErrMesg
      Integer flagFixShape
      Real*8   hStar, rQuality
      Logical  flagAnalytic, flagAuto

      Integer  iICP, iIEP, iIFE, iIEE, iXYPw, iIPEw, iIEPw, inEPw
      Integer  iIHolP, iIHolF, iIHolE
      Integer  iL1E, iL2E, iL1Et, iL2Et, inL2t, inStept
      Integer  iqE, iLFnc, iILt, inEt, itE, irSE, iiSE
      Integer  iHesP, iHesPw, idG

      Integer  i, n, nWr, nWi, mrLINTRP, miLINTRP


C ======================================================================
C group (Common blocks)
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
 
      Real*8  refXYP(2), scaXYP(2)
      Integer isON
      Common /rescale/refXYP, scaXYP, isON

C ======================================================================
      iERR = 0
      nEStar = nE

c ... unpack control
      MaxSkipE = control(1)
      MaxQItr  = control(2)
      status   = control(3)
      flagAuto = control(4).GT.0
      iPrint   = control(5)
      iErrMesg = control(6)

      If(MaxSkipE.LE.0) MaxSkipE = max(100,  nE / 20)
      If(MaxQItr .LE.0) MaxQItr  = max(1000, nE * 3)


c ... default scaling (no scaling)
      isON = 0
      Do i = 1, 2
         refXYP(i) = 0D0
         scaXYP(i) = 1D0
      End do


      miLINTRP = 10 * nP + 3 * nE + 6
      mrLINTRP =  4 * nP + MaxH + 4

      inEt = 1
      inStept = inEt + MaxF
      inL2t = inStept + 4 * MaxF
      iLFnc = inL2t + MaxF
      iILt  = iLFnc + MaxF
      iL1Et = iILt + MaxF
      iL2Et = iL1Et + 2 * MaxF
      iIHolP = iL2Et + 2 * MaxF
      iIHolF = iIHolP + MaxP
      iIHolE = iIHolF + MaxF
      iICP = iIHolE + MaxE
      iIEP = iICP + MaxP
      iIFE = iIEP + MaxP
      iIEE = iIFE + 3 * MaxE
      iL1E = iIEE + 3 * MaxE
      iL2E = iL1E + 2 * MaxE
      iIPEw = iL2E + 2 * MaxE
      iiSE  = iIPEw + 3 * nE
      iIEPw = iiSE + miLINTRP
c ... we need twice less memory for backReferences
      inEPw = iIEPw + max(6 * nE, 4 * MaxF)
      nWi   = inEPw + max(3 * MaxP, 2 * MaxF)


      iHesP = 1
      itE = iHesP + 3 * MaxP
      idG = itE + MaxF
      iqE = idG + MaxP
      iHesPw = iqE + MaxE
      iXYPw = iHesPw + 3 * nP
      irSE = iXYPw + 2 * nP
      nWr  = irSE + max(mrLINTRP, max(nE, MaxF))


      iW(1) = nWi
      iW(2) = nWr
      If(nWi.GT.MaxWi) Then
         iERR = 1001
         Goto 1000
      End if

      If(nWr.GT.MaxWr) Then
         iERR = 1002
         Goto 1000
      End if


      Do n = 1, nWr
         rW(n) = 0D0
      End do

      Do n = 1, nWi
         iW(n) = 0
      End do


c ... compute the analytic metric
      Call iniQanalytic(nP, XYP, MetricFunction, rW(iHesP))


c ... scale geometry to unit cube
      Call scale2Square(nP, nF, XYP, lbC, .TRUE.)


c ... print Ani2D header
      If(iPrint.GE.1) Write(*, 5004) Quality, nEStar, MaxSkipE, MaxQItr


c ... set up default status
      Call setStatus(flagAuto, status, iPrint)


c ... call the main module 
      flagAnalytic = .TRUE.
      flagFixShape = 1       

      Call ani2(
c group (M)
     &      nP, nPfix, MaxP, XYP, lbP, fixP,
     &      nF, nFfix, MaxF, IPF, lbF, fixF,
     &      nC,              Crv, lbC, CrvFunction,
     &      nE, nEfix, MaxE, IPE, lbE, fixE,
c group (M-EXT)
     &      iW(iICP), iW(iIEP), iW(iIFE), iW(iIEE),
     &      iW(iIHolP), iW(iIHolF), iW(iIHolE),
     &      rW(iXYPw), iW(iIPEw), rW(iHesPw), 
     &      iW(iIEPw), iW(inEPw), iW(iiSE), rW(irSE),
c group (CRV)
     &      iW(iL1Et), iW(iL2Et), rW(itE),
     &      iW(inL2t), iW(inStept), iW(inEt),
     &      iW(iLFnc), iW(iILt),
c group (Q)
     &      nEStar, hStar, 
     &      Quality, rQuality, rW(iHesP), rW(idG),
     &      iW(iL1E), iW(iL2E), rW(iqE),
c group (CONTROL)
     &      flagFixShape, flagAuto, status,
     &      MetricFunction, flagAnalytic,
     &      MaxSkipE, MaxQItr, nQItr,
     &      iPrint, iERR)


c ... rescale geometry back
      Call scale2Square(nP, nF, XYP, lbC, .FALSE.)


c ... returning sadditional information
      Do n = 1, nP
         iW(n) = iW(iICP + n - 1)
      End do
      iW(nP + 1) = nQItr

      rW(1) = 0D0
      rW(2) = rQuality
      rW(3) = hStar


 1000 If(iERR.EQ.0 .OR. iERR.EQ.1000) Return
      Call errMesMBA(iERR, 'mbaFixShape', 'See error.f for more detail')

      Return

 5004 Format('MBA: STONE FLOWER! (1997-2014), version 3.1', /,
     &       5X,'Target Quality=', F4.2, ' (nEStar:', I7, 
     &       ', SkipE:', I6, ', maxITR:', I8,')') 
      End




