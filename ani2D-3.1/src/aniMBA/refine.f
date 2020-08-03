C ======================================================================
      Subroutine uniformRefinement(nP, MaxP, nF, MaxF, nE, MaxE,  
     &                             XYP, IPF, lbF, IPE, lbE,
     &                             CrvFunction, ParCrv, lbC, IRE,
     &                             F, LDF, iW, MaxWi)
C ======================================================================
C Routine refines the input mesh and interpolates (linearly) the
C nodal function F(LDF, *).
C 
C *** Remarks:
C        1. The size of working memory is 3 * nE + nP
C ======================================================================
      Real*8  XYP(2, *)
      Integer IPE(3, *), IPF(2, *), lbF(*), lbE(*)

      EXTERNAL CrvFunction
      Real*8   ParCrv(2, *)
      Integer  lbC(*)

      Integer IRE(3, *), iW(*)
      Integer LDF
      Real*8  F(LDF, *)

C (Local variables)
      Integer iref(4), iEt(8)
      Real*8  t1, t2, t3, s
      Logical cmpE, tangled, check22, flagE

      Real*8  refXYP(2), scaXYP(2)
      Integer isON
      Common /rescale/refXYP, scaXYP, isON

      DATA    iref /1,2,3,1/
C ======================================================================
      If(scaXYP(1).LE.0D0) Then
         isON = 0
         Do i = 1, 2
            refXYP(i) = 0D0
            scaXYP(i) = 1D0
         End do
      End if

      inEP = 1
      iIEP = inEP + nP
      iEnd = iIEP + 3 * nE

      If(iEnd.GT.MaxWi) Call errMesMBA(1001, 'uniformRefinement',
     &                             'not enough working memory')


c ... compute map E -> F
      Call listE2R(nP, nFtot, nE, IPE, IRE, iW(inEP), iW(iIEP))

      nPo = nP
      nFo = nF
      nEo = nE

      nP = nPo + nFtot
      nF = 2 * nFo
      nE = 4 * nEo

      If(nP.GT.MaxP) Call errMesMBA(1003, 'uniformRefinement',
     &                          'local parameter MaxP is small')
      If(nF.GT.MaxF) Call errMesMBA(1004, 'uniformRefinement',
     &                          'local parameter MaxF is small')
      If(nE.GT.MaxE) Call errMesMBA(1006, 'uniformRefinement',
     &                          'local parameter MaxE is small')


c ... split edges
      Do 10 n = 1, nFo
         iP1 = IPF(1, n)
         iP2 = IPF(2, n)

         flagE = cmpE(iP1, iP2, iW(iIEP), iW(inEP), 0, iE)

         Do i1 = 1, 3
            i2 = iref(i1 + 1)

            jP1 = IPE(i1, iE)
            jP2 = IPE(i2, iE)

            If(check22(iP1, iP2, jP1, jP2)) Then
               kP1 = nPo + IRE(i1, iE)

               IPF(1, nFo + n) = iP1
               IPF(2, nFo + n) = kP1

               lbC(nFo + n) = 0
               lbF(nFo + n) = lbF(n)

               IPF(1, n) = kP1
               Goto 10
            End if
         End do
  10  Continue


c ... split elements
      kE = nEo 
      Do n = 1, nEo
         Do i1 = 1, 3
            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            jP1 = nPo + IRE(i1, n)
            jP3 = nPo + IRE(i3, n)

            kE = kE + 1
            IPE(1, kE) = iP1
            IPE(2, kE) = jP1
            IPE(3, kE) = jP3

            lbE(kE) = lbE(n)

            Do i = 1, 2
               XYP(i, jP1) = (XYP(i, iP1) + XYP(i, iP2)) / 2
            End do

            Do i = 1, LDF
               F(i, jP1) = (F(i, iP1) + F(i, iP2)) / 2
            End do
         End do

         Do i = 1, 3
            IPE(i, n) = nPo + IRE(i, n)
         End do
      End do


c ... split curved edges
      Do n = 1, nFo
         nC = lbC(n) 

         If(nC.GT.0) Then
            s = 1.0D0

            iloop = 0
 20         iloop = iloop + 1
            If(iloop.GT.3) Call errMesMBA(6004, 'uniformRefinement', 
     &                          'Mesh is tangled after refinement')

            t1 = ParCrv(1, n)
            t2 = ParCrv(2, n)

            s = s / 2
            t3 = s * t1 + (1 - s) * t2
            ParCrv(1, n) = t3

            kP1 = IPF(1, n)
            Call aniCrv(t3, XYP(1, kP1), nC, CrvFunction)

c  ...  check for inverted elements    
            iP1 = IPF(1, nFo + n) 
            iP2 = IPF(2, n) 

            nEt = 0
            iE1 = 0
            Do k = 1, 2
               If(cmpE(iP1, iP2, iW(iIEP), iW(inEP), iE1, iE2)) Then
                  Do i = 1, 3 
                     iEt(nEt + i) = nEo + 3 * (iE2 - 1) + i 
                  End do
                  iEt(nEt + 4) = iE2
                  nEt = nEt + 4
               End if

               iE1 = iE2
            End do

            Do i = 1, nEt
               Do j = max(i + 1, 4), nEt
                  If(tangled(iEt(i), iEt(j), XYP, IPE)) Goto 20
               End do
            End do

c  ...  add new curved edge
            lbC(nFo + n) = nC
            ParCrv(1, nFo + n) = t1
            ParCrv(2, nFo + n) = t3
         End if
      End do

      Return
      End

