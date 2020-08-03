C ======================================================================
      Subroutine ani2(
C ======================================================================
c group (M)
     &           nP, nPfix, MaxP, XYP, lbP, fixP,
     &           nF, nFfix, MaxF, IPF, lbF, fixF,
     &           nC,              Crv, lbC, CrvFunction,
     &           nE, nEfix, MaxE, IPE, lbE, fixE,
c group (M-EXT)
     &           ICP, IEP, IFE, IEE,
     &           IHolP, IHolF, IHolE,
     &           XYPw, IPEw, HesPw, 
     &           IEPw, nEPw, iSE, rSE,
c group (CRV)
     &           L1Et, L2Et, tE,
     &           nL2t, nStept, nEt,
     &           LFnc, ILt,
c group (Q)
     &           nEStar, hStar, 
     &           Quality, rQuality, HesP, detG,
     &           L1E, L2E, qE,
c group (CONTROL)
     &           flagFixShape, flagAuto, status,
     &           MetricFunction, flagAnalytic,
     &           MaxSkipE, MaxQItr, nQItr,
     &           iPrint, iERR)
C ======================================================================
      implicit none 

      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
      include 'operat.fd'
      include 'magic.fd'
C ======================================================================
C The main driver routine. 
C
C *** DATA FLOW CHART
C
C  -> put points in increasing order
C  -> dublicate initial metric
C  -> build various mesh structures (cross maps)
C  -> check initial mesh and self-check other mesh structures
C  -> add mesh points to satify 2-arm rule  
C
C  -> compute qualities of mesh elements
C  -> create a list of triangle qualities
C  -> loop of curvilinear boundaries
C  ->  -> created lists of curvilinear edges
C
C  -> create quad-tree
C  -> while( #loop <= #Baskets )
C  -----> begin infinite basket-loop
C  ---------> take an element (E) with the worst quality
C  ---------> order edges of element E
C  ---------> collect a superelement around E
C
C  ---------> loop over edges of E
C  -------------> collapse the edge
C  -------------> insert point in the middle of the edge
C  -------------> swap the edge
C
C  ---------> loop over vertices of E
C  -------------> swap an edge ending at the vertex of E
C  -------------> move the vertex of E
C  
C  ---------> add E to the basket is the operations have failed
C  ---------> terminate basket-loop if the basket is full
C  -----> clean the basket
C
C  -> calculate the distribution of bad elements
C  -> compress the mesh data structure (clean holes)
C  -> check the final mesh
C  -> check the area and the perimetr of the meshed area
C 
C ======================================================================
C group (M)
      Integer  nP, nPfix, MaxP, lbP(*), fixP(*)
      Real*8  XYP(2, *)

      Integer  nF, nFfix, MaxF, IPF(2, *), lbF(*), fixF(*)

      Integer  nC, lbC(*)
      Real*8   Crv(2, *)
      EXTERNAL CrvFunction

      Integer  nE, nEfix, MaxE, IPE(3, *), lbE(*), fixE(*)

C group (M-EXT)
      Integer ICP(*), IEP(*), IFE(3, *), IEE(3, *)
      Integer IHolP(*), IHolF(*), IHolE(*)

      Integer IPEw(3, *), IEPw(*), nEPw(*), iSE(*)
      Real*8  XYPw(2, *), HesPw(3, *), rSE(*)

C group (CRV)
      Integer L1Et(2, *), L2Et(*)
      Real*8  tE(*)

      Integer nL2t(*), nStept(4, *), nEt(*)
      Integer LFnc(*), ILt(*)

C group (Q)
      Integer nEStar
      Real*8  hStar, Quality, rQuality, HesP(3, *), detG(*), qE(*)
      Integer L1E(2, *), L2E(*), nStep(4)

c group (CONTROL)
      Integer  flagFixShape, status
      Logical  flagAuto

      Integer  MetricFunction
      EXTERNAL MetricFunction
      Logical  flagAnalytic

      Integer  MaxSkipE, MaxQItr, nQItr
      Integer  iPrint, iERR


C LOCAL VARIABLES
      Integer iw(3)
      Integer iFu(MaxS), iEu(MaxS), IPFu(2, MaxS), IPEu(3, MaxS)
      Integer iFs(MaxS), iEs(MaxS), IPFs(2, MaxS), IPEs(3, MaxS)
      Real*8  qEu(MaxS), qEs(MaxS)

      Real*8  XYPs(2), HesPs(3), rMove
      Logical flagTM, flagTM0, flagTest
      Logical flagMOVE, flagSWAP, flagCLPS1, flagCLPS2, flagINSRT

      Integer nCLPS1, nCLPS2, nINSRT, nMOVE, nSWAP, nNOTHING
      Integer mCLPS1, mCLPS2, mINSRT, mSWAP, mMOVE
      Real*8  tm1, tm2

      Real*8  domainArea, domainPerimetr, dao, dpo, rRmax, rRavg
      Real*8  dao_new, dpo_new
      Real*8  avgQ, aQuality

      Logical ifXnode

      Integer iDummy(1), info(3)
      Real    ANItime, tmdata(2)

      Integer nL2, nCrvFnc, nQItrAdd, nQItrBig
      Integer nSkipE

      Integer lF, lE, kE, nPo, nPVo, nFo, nPw, nEw
      Integer i, n, ir, iP1,iP2,iP3, iPt, iwE, icnt

C ======================================================================
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
C ======================================================================
      tm1 = ANItime(tmdata)

      nCLPS1 = 0
      nCLPS2 = 0
      nINSRT = 0
      nSWAP  = 0
      nMOVE  = 0
      nNOTHING = 0

      mCLPS1 = 0
      mCLPS2 = 0
      mINSRT = 0
      mSWAP  = 0
      mMOVE  = 0

      flagTM = .TRUE.

      nQItr = 0

c ... put initial data in increasing order
c     Do n = 1, nE
c        Call orderijk(IPE(1, n), IPE(2, n), IPE(3, n))
c     End do

c ... duplicating initial data for accurate data interpolation
      nPw = nP
      Do n = 1, nPw
         Do i = 1, 3
            HesPw(i, n) = HesP(i, n)
         End do

         Do i = 1, 2
            XYPw(i, n) = XYP(i, n)
         End do
      End do

      nEw = nE
      Do n = 1, nEw
         Do i = 1, 3
            IPEw(i, n) = IPE(i, n)
         End do
      End do


c ... building mesh data structure
      nPVo = nPfix
      nPo = nP
      nFo = nF
      Call makeMeshData(
c group (M)
     &     nP, nPfix,       XYP,      fixP,
     &     nF, nFfix, MaxF, IPF, lbF, fixF, Crv, lbC, 
     &     nE, nEfix, MaxE, IPE, lbE, fixE,
c group (M-EXT)
     &     ICP, IEP, IFE, IEE,
     &     IHolP, IHolF, IHolE,
     &     IEPw, nEPw,
c group (CONTROL)
     &     status, iERR)
      If(iERR.NE.0) Goto 9000


      If(nPVo.NE.nPfix) Then
         If(.NOT.flagAuto) Call errMesMBA(4001, 'ani2', 
     &                                'inconsistent input data')
         If(iPrint.GE.1) Write(*, 5007) nPfix - nPVo
      End if

      If(nFo.NE.nF) Then
         If(.NOT.flagAuto) Call errMesMBA(4001, 'ani2', 
     &                                'inconsistent input data')
         If(iPrint.GE.1) Write(*, 5008) nF - nFo
      End if


c ... checking (and self-checking) of the initial mesh
      Call meshAudit(nP, XYP,
     &               nF, IPF, Crv, lbC, CrvFunction,
     &               nE, IPE, lbE, 
     &               ICP, IFE, IEE, rRmax, rRavg, status)


c ... remove boundary triangles
      If(ifXnode(status, ANIForbidBoundaryElements)) Then
         kE = nE
         Do 20 n = 1, kE
            Do i = 1, 3
               iPt = IPE(i, n)
               If(ifXnode(ICP(iPt), jInode)) Goto 20
            End do

            nQItr = nQItr + 1

            Call makSE(n, IEP, IPF, IPE, IFE, IEE, qE, MaxS,
     &                 lF, lE, iFu, iEu, IPFu, IPEu, qEu,
     &                 status)

            Call splitTriangle(n,
     &                         nP, MaxP, XYP,  ! group M
     &                         nE, MaxE, IPE, lbE, 
     &                         ICP, IEP, IFE, IEE,  ! group M-EXT
     &                         IHolP, IHolE,
     &                         HesP, detG, qE,  ! group Q
     &                         lF, lE, iFu, iEu, IPFu, IPEu, qEu,  ! group W
     &                         flagTest)

            If(.NOT.flagTest) Call errMesMBA(4103, 'ani2.f', 
     &                       'The input data contradicts status')
 20      Continue
      End if


c ... compute quality of elements
      Call makQ(nP, XYP, 
     &          nE, nEfix, IPE, fixE,
     &          nEStar, hStar, HesP, detG, qE,
     &          flagFixShape)


c ... update qualities of tangled elements
      If(ifXnode(status, ANIUntangleMesh)) Then
         Do n = 1, nE
            Call updQa(n, XYP, IPE, IEE, qE)
         End do
      End if


c ... initilize the ordered list of element qualities
      Do n = 1, nE
         rSE(n) = qE(n)
         L2E(n) = n
      End do
      Call DSORT(rSE, L2E, nE, 2, iERR)
      If(iERR.NE.0) Goto 9000

c nStep(1) - typical interval length 
c nStep(2) - rank of interval length is [nStep(1)-nStep(2),nStep(1)+nStep(2)]
c nStep(3) - ipos2 of the second part of L2E
c nStep(4) - output channel in case of debugging (=0 for no debugging)
c     nStep(1) = 7 * log(real(nEStar)) 
      nStep(1) = 2 * log(real(nEStar)) 
c     nStep(1) = sqrt(real(nEStar)) 
      nStep(2) = nStep(1) / 4 
      nStep(3) = MaxE
      nStep(4) = 1

      Call lstMak(nE, L1E, L2E, nL2, nStep, IHolE)


c ... output of statistics
      If(iPrint.GE.3) Then
         Call quality_log(nE, L1E, L2E, qE)
      End if


c ... count the number of different curvilinear functions
      Call calCrvFnc(nF, lbC, LFnc, nCrvFnc)


c ... creating ordered (and fixed) lists of initial parametrizations
      ir = 1
      Do n = 1, nCrvFnc
         ILt(n) = ir
         Call tEMak(tE(ir), nEt(n), nF, Crv, lbC, LFnc(n))

         If(nEt(n).GT.0) Then
            Do i = 1, nEt(n)
               rSE(i) = tE(ir + i - 1)
               L2Et(ir + i - 1) = i
            End do
            Call DSORT(rSE, L2Et(ir), nEt(n), 2, iERR)
            If(iERR.NE.0) Goto 9000

            nStept(1, n) = sqrt(real(nEt(n)))
            nStept(2, n) = 0
            nStept(3, n) = MaxF
            nStept(4, n) = 0 
            Call lstMak(nEt(n), L1Et(1, ir), L2Et(ir),
     &                  nL2t(n), nStept(1, n), iDummy)
         End if
         ir = ir + nEt(n)
      End do


c ... initialize and populate the 4-tree
      If( .NOT.flagAnalytic ) Then
         info(1) = 1
         info(2) = 0
         Call LINTRP(nEw, IPEw, nPw, XYPw, 3, HesPw, 0, XYPs,
     &               HesPs,  iSE, rSE, info)
      End if


c ... output of the initial mesh quality
      tm2 = ANItime(tmdata)
      If(iPrint.GE.1) Then
         aQuality = avgQ(nE, qE, L1E, L2E)

         Write(*, 5006) aQuality, rRmax, rRavg, status
         Write(*, 5000) nQItr, qE(L2E(1)), nP, nF, nE, tm2 - tm1
      End if

      flagTM0 = .FALSE.
      If(qE(L2E(1)).LT.0D0) flagTM0 = .TRUE.


C ... START OF THE MAIN LOOP
      nQItrAdd = 0
      nQItrBig = 0
      iERR = 0

 100  nQItrBig = nQItrBig + 1
      If(nQItrBig.GT.ANIMaxBaskets) Then
         iERR = 1000
         Goto 1000
      End if


c ... output of the intermediate mesh quality
      tm2 = ANItime(tmdata)
      If(iPrint.GE.4 .AND. nSkipE.GT.0)
     &   Write(*, 5000) nQItr, qE(L2E(1)), nP, nF, nE, tm2 - tm1

      If(iPrint.GE.4 .AND. nQItrBig.GT.1 .AND. nSkipE.GT.0) 
     &   Write(*, 5001) nQItrBig - 1, nSkipE


      nSkipE = 0
 300  nQItr = nQItr + 1

      If(nQItr.GT.MaxQItr + nQItrAdd) Then
         If(nQItrBig.GE.1) Then
            nQItrAdd = nSkipE
            Goto 100
         End if

         iERR = 1000
         Goto 1000
      End if

      iwE = L2E(1)
      Do i = 1, nSkipE
         iwE = L1E(2, iwE)
      End do


c ... check if we need to removed the flag of tangled mesh
      If(qE(L2E(1)).GT.0D0 .AND. flagTM) Then
         flagTM = .FALSE.

         Call delXnode(status, ANIUntangleMesh)

         dao = domainArea(nE + IHolE(1), XYP, IPE)
         dpo = domainPerimetr(nF + IHolF(1), XYP, IPF)
      End if


c ... output of the intermediate mesh quality
      If((nQItr / 5000) * 5000.EQ.nQItr) Then
         tm2 = ANItime(tmdata)
         If(iPrint.GE.2)
     &      Write(*, 5000) nQItr, qE(L2E(1)), nP, nF, nE, tm2 - tm1
      End if

      rQuality = qE(iwE)
      If(rQuality.GT.Quality) Then
         If(nSkipE.EQ.0) Goto 1000
         Goto 100
      End if

c ... order edges of the triangle whose quality we want to increase
      iP1 = IPE(1, iwE)
      iP2 = IPE(2, iwE)
      iP3 = IPE(3, iwE)

      Call calQF(HesP(1, iP1), detG(iP1), XYP(1, iP1),
     &           HesP(1, iP2), detG(iP2), XYP(1, iP2),
     &           HesP(1, iP3), detG(iP3), XYP(1, iP3),
     &           hStar, iw)

      Call makSE(iwE, IEP, IPF, IPE, IFE, IEE, qE, MaxS,
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
     &           status)


c ... setup flags for the local mesh modifications
      flagMOVE  = .FALSE.
      flagSWAP  = .FALSE.
      flagCLPS1 = .FALSE. 
      flagCLPS2 = .FALSE. 
      flagINSRT = .FALSE.

      Do i = 1, 3
         If(iCLPSop.EQ.1) Then
            Call clpsF1(
     &           iw(i), iwE, 5D-1, 5D-1,
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
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
     &                   iFs, iEs, IPFs, IPEs, qEs,
     &           flagTest)

            mCLPS1 = mCLPS1 + 1

            If(flagTest) Then
               nCLPS1 = nCLPS1 + 1
               flagCLPS1 = .TRUE. 
               Goto 400
            End if
         End if

         If(iINSRTop.EQ.1) Then
            Call insrtPoint(
     &           iw(i), iwE,
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
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
     &                   iFs, iEs, IPFs, IPEs, qEs,
     &           flagTest)

            mINSRT = mINSRT + 1
            If(flagTest) Then
               nINSRT = nINSRT + 1
               flagINSRT = .TRUE.
               Goto 400
            End if
         End if

         If(iSWAPop.EQ.1) Then
            Call swapEdge(
     &           iw(i), iwE,
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
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
     &                   iFs, iEs, IPFs, IPEs, qEs,
     &           flagTest)

            mSWAP = mSWAP + 1

            If(flagTest) Then
               nSWAP = nSWAP + 1
               flagSWAP  = .TRUE.
               Goto 400
            End if
         End if
      End do


      Do i = 1, 3
         If(iCLPSop.EQ.1) Then
            Call clpsF2(
     &           i, iwE,
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
     &           lF, lE, iFu, iEu,
     &           flagTest)

            mCLPS2 = mCLPS2 + 1

            If(flagTest) Then
               nCLPS2 = nCLPS2 + 1
               flagCLPS2 = .TRUE. 
               Goto 400
            End if
         End if

         If(iMOVEop.EQ.1) Then
            Call movePoint(
     &           i, iwE,
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
     &           lF, lE, iFu, iEu, IPFu, IPEu, qEu,
     &                   iFs, iEs, IPFs, IPEs, qEs,
     &           icnt, rMove, flagTest)

            mMOVE = mMOVE + 1

            If(flagTest) Then
               nMOVE = nMOVE + 1
               flagMOVE = .TRUE.
               Goto 400
            End if
         End if
      End do

      nNOTHING = nNOTHING + 1

 400  Continue 
c ... debug section
c     Call debug(nP, nF, nE,
c    &           XYP, IPF, lbC, IPE, IFE, IEE, lbE,
c    &           CrvFunction, Crv, lbC,
c    &           IHolP, IHolF, IHolE,
c    &           ICP, status, nQItr,
c    &           flagMOVE, flagSWAP, flagCLPS1, flagCLPS2, flagINSRT)
c ... end debug section
      If(flagTest) Goto 300


      nSkipE = nSkipE + 1
      If(nSkipE.GT.MaxSkipE .OR. nSkipE.GE.nE) Goto 100
      Goto 300

 1000 rQuality = qE(L2E(1))

      tm2 = ANItime(tmdata)
      If(iPrint.GE.1) 
     &   Write(*, 5000) nQItr - 1, qE(L2E(1)), nP, nF, nE, tm2 - tm1


C ... calculating the number of bad triangles
      If(iPrint.GE.3) Then
         Write(*,5009) tmdata(1), tmdata(2)

         Write(*,5003) nCLPS1, nINSRT, nSWAP, nCLPS2, nMOVE, nNOTHING,
     &                 mCLPS1, mINSRT, mSWAP, mCLPS2, mMOVE

         Call quality_plain(nE, L1E, L2E, qE)
      End if


C ... remove 'holes' from the final grid
      nFo = nF 
      Call updM(
c group (M)
     &     nP, nPfix, XYP, fixP,
     &     nF,      IPF, lbF, Crv, lbC, 
     &     nE,      IPE, lbE,
c group (M-EXT)
     &     ICP, IFE, IEE, 
     &            IHolP, IHolF, IHolE,
c group (MISC)
     &     status, qE, nEPw)

c     Call draw(nP, nF, nE, XYP, ICP, IPF, IPE, 'fin.ps')
c     Call draw_Q(nP, nE, XYP, IPE, qE, Quality, 'qE.ps')

      Call meshAudit(nP, XYP,
     &               nF, IPF, Crv, lbC, CrvFunction,
     &               nE, IPE, lbE, 
     &               ICP, IFE, IEE, rRmax, rRavg, status)

      If(iPrint.GE.1) Then
         aQuality = avgQ(nE, qE, L1E, L2E)
         Write(*, 5006) aQuality, rRmax, rRavg, status
      End if

      If(nCrvFnc.EQ.0 .AND. nF.GT.0 .AND. .NOT.flagTM0) Then 
         dao_new = domainArea(nE, XYP, IPE)
         dpo_new = domainPerimetr(nF, XYP, IPF)

         If(iPrint.GE.2) Write(*, 5005) dao, dpo, dao_new, dpo_new

         dao = dabs(dao - dao_new) / dao
         dpo = dabs(dpo - dpo_new) / dpo

         If(dao.GT.1D-8) Call errMesMBA(6005, 'ani2.f', 
     &                       'Lose of the domain area')

         If(dpo.GT.1D-8 .AND. nF.EQ.nFo) Call errMesMBA(6005, 'ani2.f', 
     &                       'Lose of the total boundary lenght')
      End if


 5000 Format('ITRs:', I7, '  Q=', E10.4, '   #V#B#T:', 3I8,
     &     '   tm=', F7.2, 's')

 5001 Format('Cleaning basket #', I3, ' with', I4, ' elements')

 5003 Format(/,
     & 'Collapse I-B   Insert   Swapping   Collapse I-E   Moving  Nothin
     &g',/, I12,I9,I11,I15,I9,I9,/,I12,I9,I11,I15,I9)

 5005 Format('Domain area and total      boundary: ', 2E14.6,/,
     &       'Domain area and user-given boundary: ', 2E14.6)

 5006 Format('Avg Quality =', E11.4, ',  R/r(max,avg):', 2E10.3,
     &       ', status =', I4)

 5007 Format('Warning:', I6, ' new fix vertices have been added') 
 5008 Format('Warning:', I6, ' new edges have been added') 

 5009 Format('User time:', F7.2, ',  system time:', F7.2)

 9000 Return
      End


  

