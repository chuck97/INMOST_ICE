C ====================================================================
      Subroutine makeMeshData(
C ====================================================================
c group (M)
     &           nP, nPfix,       XYP,      fixP,
     &           nF, nFfix, MaxF, IPF, lbF, fixF, Crv, lbC, 
     &           nE, nEfix, MaxE, IPE, lbE, fixE,
c group (M-EXT)
     &           ICP, IEP, IFE, IEE,
     &           IHolP, IHolF, IHolE,
     &           IEPw, nEPw,
c group (CONTROL)
     &           status, iERR)
C ====================================================================
      implicit none
      include 'makS.fd'
      include 'colors.fd'
      include 'status.fd'
C ====================================================================
C group (M)
      Integer  nP, nPfix, fixP(*)
      Real*8  XYP(2, *)

      Integer  nF, nFfix, MaxF, IPF(2, *), lbF(*), fixF(*), lbC(*)
      Real*8   Crv(2, *)

      Integer  nE, nEfix, MaxE, IPE(3, *), lbE(*), fixE(*)

C group (M-EXT)
      Integer ICP(*), IEP(*)
      Integer IFE(3, *), IEE(3, *)
      Integer IHolP(*), IHolF(*), IHolE(*)
      Integer IEPw(*), nEPw(*)

c group (CONTROL)
      Integer status
      Integer iERR


C LOCAL VARIABLES
      Integer Mlist(2, MaxS), Clist(MaxS)
      Integer IPFs(2)

      Integer ip(4) 
      Logical ifXnode, cmpP, cmpE, sharpAngle

      Integer i, k, n, i1, i2, j1
      Integer ip1, ip2, iPt, iF, iFt, iE, iE2, iEt, ilist
      Integer  iBNDs, ic1, ic2, icmat, icbnd, icmax

C ====================================================================
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
C ====================================================================
      iERR = 0

      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      IHolP(1) = 0
      IHolF(1) = 0
      IHolE(1) = 0


C ... create an auxiliary structure
      Call backReferences(nP, nE, 3, 3, IPE, nEPw, IEPw)

      Do n = 1, MaxE
         Do i = 1, 3
            IEE(i, n) = 0
            IFE(i, n) = 0
         End do
      End do

C ... create IEE & IEP
      Do n = 1, nE
         Do i1 = 1, 3
            i2 = ip(i1 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)

            IEP(ip1) = n

            If(cmpE(ip1, ip2, IEPw, nEPw, n, iE2)) Then
               IEE(i1, n) = iE2
            End if
         End do
      End do


c ... creating an auxiliary structure
      Call backReferences(nP, nF, 2, 2, IPF, nEPw, IEPw)

      Do n = 1, nP
         ICP(n) = 0
      End do

C ... create IFE  
      Do n = 1, nE
         Do i1 = 1, 3
            i2 = ip(i1 + 1)

            ip1 = IPE(i1, n)
            ip2 = IPE(i2, n)

            If(cmpE(ip1, ip2, IEPw, nEPw, 0, iF)) Then
               IFE(i1, n) = iF
            End if
         End do
      End do


C ... count existing surfaces
      ilist = 0
      icmax = 0
      Do n = 1, nF
         ic1 = lbF(n)
         Call findSE(ilist, Clist, ic1, k)

         If(k.EQ.0) Then
            ilist = ilist + 1
            If(ilist.GT.MaxS) Call errMesMBA(1007, 'makeMeshData',
     &                             'local parameter MaxS is small')

c  ...  user-given interface
            Do i = 1, 2
               Mlist(i, ilist) = 0
            End do

            Clist(ilist) = ic1
            icmax = max(ic1, icmax) 
         End if
      End do


C ... create temporary surfaces (material and boundary)
      iDomBnd = icmax + 1
      iMatBnd = iDomBnd + 3 * nE

      icbnd = 0
      icmat = 0

      Do n = 1, nE
         Do 100 i1 = 1, 3
            iFt = IFE(i1, n)
            iEt = IEE(i1, n)

            If(iEt.GT.n .OR. iFt.GT.0) Goto 100

            ic1 = lbE(n)
            If(iEt.GT.0) Then
               ic2 = lbE(iEt)
            Else If(iFt.EQ.0) Then
               ic2 = 0
            Else
               ic2 = ic1
            End if

c  ...  order materials such that ic1 > ic2
            If(ic1.LT.ic2) Call swapii(ic1, ic2)

            If(ic1.NE.ic2) Then
               Do i = 1, ilist
                  If(Mlist(1, i).EQ.ic1 .AND.
     &               Mlist(2, i).EQ.ic2) Then
                     iBNDs = Clist(i)
                     Goto 10
                  End if
               End do

               If(ic2.EQ.0) Then
                  icbnd = icbnd + 1
                  iBNDs = iDomBnd + icbnd
               Else
                  icmat = icmat + 1
                  iBNDs = iMatBnd + icmat
               End if

               ilist = ilist + 1
               If(ilist.GT.MaxS) Call errMesMBA(1007, 'makeMeshData',
     &                                'local parameter MaxS is small')
               Mlist(1, ilist) = ic1
               Mlist(2, ilist) = ic2
               Clist(ilist) = iBNDs

 10            Continue
               If(nF.EQ.MaxF) Goto 9200
               Call facAdd(iF, nF, MaxF, IHolF)

               i2 = ip(i1 + 1)

               IPFs(1) = IPE(i1, n)  
               IPFs(2) = IPE(i2, n)

               Call facUpd(1, IPF, lbF, Crv, lbC,
     &                     iF, IPFs, 0, 0, iBNDs, 0D0, 0D0)

c  ...  update IFE
               IFE(i1, n) = iF
               If(iEt.GT.0) Then
                  Do j1 = 1, 3
                     If(IEE(j1, iEt).EQ.n) IFE(j1, iEt) = iF
                  End do
               End if
            End if
 100     Continue
      End do
      

C ... create ICP
      Do n = 1, nF
         Do i = 1, 2
            iP1 = IPF(i, n)
            Call addXnode(ICP(iP1), jSnode)
         End do
      End do


      Do n = 1, nE
         Do i1 = 1, 3
            iE = IEE(i1, n)
            If(iE.EQ.0) Then
               i2 = ip(i1 + 1)

               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)

               Call addXnode(ICP(iP1), jBnode)
               Call addXnode(ICP(iP2), jBnode)
            End if
         End do
      End do

      Do n = 1, nP
         If(.NOT.ifXnode(ICP(n), jBnode)) Call addXnode(ICP(n), jInode)
      End do

      Do n = 1, nPfix
         iP1 = fixP(n)
         Call addXnode(ICP(iP1), jVnode)
      End do


c ... first realization of the fixed edges and triangles
      Do n = 1, nFfix
         iFt = fixF(n)

         Do i = 1, 2
            iPt = IPF(i, iFt)
            Call addXnode(ICP(iPt), jTnode)
         End do
      End do

      Do n = 1, nEfix
         iEt = fixE(n)

         Do i = 1, 3
            iPt = IPE(i, iEt)
            Call addXnode(ICP(iPt), jTnode)
         End do
      End do


c ... color the points (vertices)
      Call backReferences(nP, nF, 2, 2, IPF, nEPw, IEPw)

      Do n = 1, nP
         If(cmpP(n, IEPw, nEPw, lbF)) Call addXnode(ICP(n), jVnode)
      End do


c ... color sharp angles
      Do n = 1, nP
         If(.NOT.ifXnode(ICP(n), jVnode)) Then
            If(sharpAngle(n, IEPw, nEPw, XYP, IPF)) 
     &         Call addXnode(ICP(n), jVnode)
         End if
      End do


c ... color the points (fixed boundary points)
      If(ifXnode(status, ANIFixBoundaryPoints)) Then
         Do n = 1, nE
            Do i1 = 1, 3
               If(IEE(i1, n).EQ.0) Then
                  i2 = ip(i1 + 1)

                  iP1 = IPE(i1, n)
                  iP2 = IPE(i2, n)
                  Call addXnode(ICP(iP1), jVnode)
                  Call addXnode(ICP(iP2), jVnode)
               End if
            End do
         End do
      End if


c ... color the points (fixed boundary edges)
      If(ifXnode(status, ANIFixBoundaryEdges)) Then
         Do n = 1, nE
            Do i1 = 1, 3
               If(IEE(i1, n).EQ.0) Then
                  i2 = ip(i1 + 1)

                  iP1 = IPE(i1, n)
                  iP2 = IPE(i2, n)
                  Call addXnode(ICP(iP1), jTnode)
                  Call addXnode(ICP(iP2), jTnode)
               End if
            End do
         End do
      End if

      Return

 9200 Continue
      Call errMesMBA(1004, 'makeMeshData', 
     &               'local variable MaxF is too small')
      End



C ======================================================================
      Subroutine updM(
C ======================================================================
c group (M)
     &           nP, nPfix, XYP, fixP,
     &           nF,        IPF, lbF, Crv, lbC, 
     &           nE,        IPE, lbE,
c group (M-EXT)
     &           ICP, IFE, IEE, 
     &            IHolP, IHolF, IHolE,
c group (MISC)
     &           status, qE, IPw)
C ======================================================================
      implicit none
      include 'status.fd'
C ======================================================================
C *** Remarks:
C     1. size of working memory is equal to MaxP, IPw(MaxP)
C ======================================================================
C group (M)
      Integer  nP, nPfix, fixP(*)
      Real*8  XYP(2, *)

      Integer  nF, IPF(2, *), lbF(*), lbC(*)
      Real*8   Crv(2, *)

      Integer  nE, IPE(3, *), lbE(*)

C group (M-EXT)
      Integer  ICP(*), IFE(3, *), IEE(3, *)
      Integer IHolP(*), IHolF(*), IHolE(*)

c group (MISC)
      Integer  status, IPw(*)
      Real*8  qE(*)


C LOCAL VARIABLES
      Integer  iP, iF, iE, iFt, iEt, lP,lF,kF, lE,kE, mP, mF, mE
      Integer  i, j, k, m, n
      Logical ifXnode
 
C ======================================================================
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
C ======================================================================
c ... delete references to material or fictitious faces (Id >= iDomBnd)
      If(ifXnode(status, ANIDeleteTemporaryEdges)) Then
         mF = nF + IHolF(1)
         mE = nE + IHolE(1)

         Do n = 1, mE
            If(IPE(1, n).GT.0) Then
               Do i = 1, 3
                  iFt = IFE(i, n)
                  If(iFt.GT.0) Then
                     If(lbF(iFt).GE.iDomBnd) IFE(i, n) = 0
                  End if
               End do
            End if
        End do

         Do n = 1, mF
            If(IPF(1, n).GT.0 .AND. lbF(n).GE.iDomBnd) Then
               Call facDel(n, nF, IPF, lbC, IHolF)
            End if
         End do
      End if


      nP = nP + IHolP(1)
      nF = nF + IHolF(1)
      nE = nE + IHolE(1)


c ... fill in holes in points
      Do n = 1, nP
         IPw(n) = 0
      End do

      lP = IHolP(1)
      Do n = 1, lP
         iP = IHolP(n + 1)
         IPw(iP) = -1
      End do

      mP = 0
      Do n = 1, nP
         If(IPw(n).EQ.0) Then
            mP = mP + 1
            IPw(n) = mP

            Do i = 1, 2
               XYP(i, mP) = XYP(i, n)
            End do

            ICP(mP) = ICP(n)
         End if
      End do


      Do n = mP + 1, nP
         ICP(n) = 0
      End do

      Do n = 1, nPfix
         fixP(n) = IPw(fixP(n))
      End do

      Do n = 1, nF
         If(IPF(1, n).GT.0) Then
            Do i = 1, 2
               IPF(i, n) = IPw(IPF(i, n))
            End do
         End if
      End do

      Do n = 1, nE
         If(IPE(1, n).GT.0) Then
            Do i = 1, 3
               IPE(i, n) = IPw(IPE(i, n))
            End do
         End if
      End do

      nP = mP


c ... fill in holes in edges
      lF = IHolF(1)
      Do 200 n = 1, lF
         iF = IHolF(n + 1)

         Do m = nF, iF + 1, -1
            If(IPF(1, m).NE.0) Then
               kF = m
               Goto 20
            End if
         End do
         Goto 200

 20      Do i = 1, 2
            IPF(i, iF) = IPF(i, kF)
            Crv(i, iF) = Crv(i, kF)
         End do

         lbF(iF) = lbF(kF)
         lbC(iF) = lbC(kF)

         Do k = 1, nE
            Do i = 1, 3
               If(IFE(i, k).EQ.kF) IFE(i, k) = iF
            End do
         End do

         IPF(1, kF) = 0
 200     nF = nF - 1


c ... fill in holes in elements
      lE = IHolE(1)
      Do 300 n = 1, lE
         iE = IHolE(n + 1)

         Do m = nE, iE + 1, -1
            If(IPE(1, m).NE.0) Then
               kE = m
               Goto 30
            End if
         End do
         Goto 300

 30      Do i = 1, 3
            IPE(i, iE) = IPE(i, kE)
         End do

         qE(iE) = qE(kE)

         lbE(iE) = lbE(kE)

c  ...   update auxiliary structures
         Do i = 1, 3
            IFE(i, iE) = IFE(i, kE)
            IEE(i, iE) = IEE(i, kE)
            iEt = IEE(i, iE)
            If(iEt.GT.0) Then
              Do j = 1, 3
                 If(IEE(j, iEt).EQ.kE) IEE(j, iEt) = iE
              End do
            End if
         End do


         IPE(1, kE) = 0
 300     nE = nE - 1

      Return
      End



C ======================================================================
      Logical Function cmpE(i1, i2, IEP, nEP, iE1, iE2)
C ======================================================================
C Search for element iE2, other than iE1, in the intersection of
C sublists associated with i1 and i2. The routine returns value
C .FALSe. and iE2 = 0 when thre is no such element. 
C ======================================================================
      Integer IEP(*), nEP(*)

C group (Local variables)
      Integer ib(2), ie(2), ip(2)

      ip(1) = i1
      ip(2) = i2
      Do i = 1, 2
         If(ip(i).EQ.1) Then
            ib(i) = 1
         Else
            ib(i) = nEP(ip(i) - 1) + 1
         End if
         ie(i) = nEP(ip(i))
      End do

      Do 10 i = ib(1), ie(1)
         iE2 = IEP(i)
         If(iE2.EQ.iE1) Goto 10
         Do j = ib(2), ie(2)
            If(iE2.EQ.IEP(j)) Then
               cmpE = .TRUE.
               Goto 1000
            End if
         End do
 10   Continue

      iE2 = 0
      cmpE = .FALSE.
 1000 Return
      End



C ======================================================================
      Logical Function cmpF(i1, IFP, nFP, iF1, iF2)
C ======================================================================
C Search of one coincident (different from iF1) in two lists
C ======================================================================
      Integer IFP(*), nFP(*)

      ip = i1
      If(ip.EQ.1) Then
         ib = 1
      Else
         ib = nFP(ip - 1) + 1
      End if
      ie = nFP(ip)

      Do 10 i = ib, ie
         iF2 = IFP(i)
         If(iF2.EQ.iF1) Goto 10

         cmpF = .TRUE.
         Goto 1000
 10   Continue

      cmpF = .FALSE.
 1000 Return
      End



C ======================================================================
      Logical Function cmpP(iP, IFP, nFP, lbF)
C ======================================================================
C cmpP = TRUE if point iP belongs to two faces with different colors. 
C Otherwise cmpP = FALSE.
C ======================================================================
      Integer IFP(*), nFP(*), lbF(*)

      If(iP.EQ.1) Then
         ib = 1
      Else
         ib = nFP(iP - 1) + 1
      End if
      ie = nFP(iP)

      Do i = ib, ie
         ICF1 = lbF(IFP(i))
         Do j = i + 1, ie
            ICF2 = lbF(IFP(j))
            If(ICF1.NE.ICF2) Then
               cmpP = .TRUE.
               Goto 1000
            End if
         End do
      End do

      cmpP = .FALSE.
 1000 Return
      End



C ======================================================================
      Logical Function sharpAngle(n, IFP, nFP, XYP, IPF)
C ======================================================================
C Routine checks that the angle between two surface edges is
C less than ANIMaxSharpAngle.
C ======================================================================
      implicit none
      include 'magic.fd'
C ======================================================================
      Integer IFP(*), nFP(*), IPF(2, *), n
      Real*8  XYP(2, *)

C LOCAL VARIABLES
      Real*8  dotMul, calNorm
      Real*8  v1(2), v2(2), angle, norms

      Integer i, ib, ie, iP1, iP2, iPt, iF1, iF2, nFt
C ======================================================================
      sharpAngle = .FALSE.

      If(n.EQ.1) Then
         ib = 1
      Else
         ib = nFP(n - 1) + 1
      End if
      ie = nFP(n)

      nFt = ie - ib + 1
      If(nFt.LT.2) Goto 1000
      If(nFt.GT.2) Goto 500


c ... find two neighboring points for point n
      iF1 = IFP(ib)
      iF2 = IFP(ie)

      Do i = 1, 2
         iPt = IPF(i, iF1)
         If(iPt.NE.n) iP1 = iPt

         iPt = IPF(i, iF2)
         If(iPt.NE.n) iP2 = iPt
      End do

      Do i = 1, 2
         v1(i) = XYP(i, iP1) - XYP(i, n)
         v2(i) = XYP(i, iP2) - XYP(i, n)
      End do  

      angle = dotMul(v1, v2)
      norms = calNorm(v1) * calNorm(v2)

      If(angle.LE.ANIMaxSharpAngle*norms) Goto 1000
   
  500 sharpAngle = .TRUE.

 1000 Return
      End



C ======================================================================
      Subroutine RandR(XY1, XY2, XY3, rOut, rIn)
C ======================================================================
C Routine computes curcumscribed and inscribed radii for the 
c triangle defined by three vertices.
C ======================================================================
      Real*8  XY1(3), XY2(3), XY3(3)
      Real*8  rOut, rIn, edge_length

      Real*8  a, b, c, s
C ======================================================================
      a = edge_length(XY1, XY2)
      b = edge_length(XY1, XY3)
      c = edge_length(XY2, XY3)

      s = (a + b + c) / 2
 
      rOut = a * b * c / (4 * dsqrt(s * (a + b - s) * 
     &                                  (a + c - s) * (b + c - s)))

      rIn = a * b * c / (4 * ROut * s)

      Return
      End



C ======================================================================
      Real*8  Function domainArea(nE, XYP, IPE)
C ======================================================================
      Real*8  XYP(2, *)
      Integer IPE(3, *)

c (Local variables)
      Real*8 s

      s = 0D0
      Do 10 n = 1, nE
         If(IPE(1, n).EQ.0) Goto 10

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         s = s + dabs((XYP(1, iP2) - XYP(1, iP1)) *
     &                (XYP(2, iP3) - XYP(2, iP1)) -
     &                (XYP(1, iP3) - XYP(1, iP1)) *
     &                (XYP(2, iP2) - XYP(2, iP1)))
 10   Continue

      domainArea = 5D-1 * s 
      Return
      End



C ======================================================================
      Real*8  Function domainPerimetr(nF, XYP, IPF)
C ======================================================================
      Real*8  XYP(2, *)
      Integer IPF(2, *)

c (Local variables)
      Real*8 s
C ======================================================================
      s = 0D0
      Do 10 n = 1, nF
         If(IPF(1, n).EQ.0) Goto 10

         iP1 = IPF(1, n)
         iP2 = IPF(2, n)

         s = s + dsqrt((XYP(1, iP2) - XYP(1, iP1)) ** 2 +
     &                 (XYP(2, iP2) - XYP(2, iP1)) ** 2)
 10   Continue

      domainPerimetr = s 
      Return
      End



