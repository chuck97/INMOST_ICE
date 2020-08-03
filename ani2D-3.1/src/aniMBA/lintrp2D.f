C ======================================================================
      Subroutine LINTRP(nt,tri, nv,vrt, LDF,F, NXY,XY, G, iW,rW, info)
C ======================================================================
      implicit none
      Include 'lintrp.fd'
C ======================================================================
C Routine interpolates values of a piecewise linear finite element at
C a given set of points inside the triangulation.
C
C   nt         - [input] number of mesh triangles
C   tri(3, nt) - [input] list of triangles 
C   nv         - [input] number of mesh vertices
C   vrt(2, nv) - [input] coordinates of vertices
C   LDF        - [input] leading dimension of vector function F and G
C   F(LDF, nv) - [input] values of function at mesh vertices
C
C   NXY        - [input] number of points where the value of function 
C                        to be determined
C   xy(2, NXY) - [input] coordnates of these points
C   G(LDF, NXY)- [ouput] the interpolated values
C 
C   iW(4*NQ + 3*nt + nv + 1) - auxiliary array of integers
C   rW(2*NQ)                 - auxiliary array of double precisions
C
C   info - [input/output] control array
C              flagInit    = info(1)>0;  If .TRUE. will build the tree
C              interpolate = info(2)>0;  If .TRUE. will interpolate
C              info(3) = iT returns the triangle number for the 1st point
C
C   Remarks:
C       The constant NQ lies in range   1/3 * nv <  NQ  < 2/3 * nv  and
C       dependes on spacing of triangulaion. To guaranty the robustness
C       of subroutine set NQ to  2/3 * nv. The exact sizes of auxiliary
C       storage been required  by  subroutine are returned in the first
C       elements of correspondent working arrays.
C
C ======================================================================
      Integer nt, nv, LDF, NXY, tri(3, nt)
      Real*8  vrt(2, nv), F(LDF, nv), xy(2, NXY), G(LDF, NXY)

      Integer iW(*), info(3)
      Real*8  rW(*)


C LOCAL VARIABLEs
      Integer i, nQT, iQT, inEP, iIEP
      Real*8  h(MaxH)
      Logical flagInit, flagCvxCrvBnd

C ======================================================================
      flagInit = info(1).GT.0

      If(flagInit) Then
         iQT = 5
         Call initQT(nQT, iW(iQT), h, rW(MaxH + 1))

         Do i = 1, nv
            Call addQT(nQT, iW(iQT), h, rW(MaxH + 1), vrt, i)
         Enddo

         inEP = iQT + 4 * nQT
         iIEP = inEP + nv + 1
      Else
         iQT  = iW(1)
         nQT  = iW(2)
         inEP = iW(3)
         iIEP = iW(4)

         Do i = 1, MaxH
            h(i) = rW(i)
         End do
      End if

        flagCvxCrvBnd = .true.
C This flag indicates that the domain can be non-convex. 
C In this case, an error will be issued if no triangle is found for
C a given point and tolerance PREC=1e-6.
C
C If the domain is convex and boundary is curved, it is possible that 
C a given point is inside the domain but outside of the mesh. In this 
C case, the tolerance PREC will be relaxed until an element is found. 

      Call interpolate(nt, tri, nv, vrt, LDF, F, NXY, xy, G,
     &                 iW(iQT), rW(MaxH + 1), iW(inEP), iW(iIEP),
     &                 h, info, flagCvxCrvBnd)

c     iW(1) = nQT * 4 + 3 * NT + NV + 1
c     rW(1) = nQT * 2

      If(flagInit) Then
         iW(1) = iQT
         iW(2) = nQT
         iW(3) = inEP
         iW(4) = iIEP

         Do i = 1, MaxH
            rW(i) = h(i)
         End do
      End if

      Return
      End



C ======================================================================
      subroutine initQT(nQT, QT, h, XYc)
C ======================================================================
C Routine initializes a quadtree sructure
C ======================================================================
      implicit none
      Include 'lintrp.fd'

      Integer  nQT, QT(4, *)
      Real*8   XYc(2, *), h(*)

      Integer  i, k

C ======================================================================
        h(1) = 0.5

      Do k = 2, MaxH
         h(k) = h(k-1) / 2
      End do

      nQT = 1

      Do i = 1, 4
         QT(i, nQT) = 0
      End do

      XYc(1, nQT) = 0.5
      XYc(2, nQT) = 0.5

      Return
      End



C ======================================================================
      Subroutine setIJ(XYc, XY, I, J)
C ======================================================================
C Finds the quadrant of point XY with respect to the center XYc
C ======================================================================
      implicit none
      Integer I, J
      Real*8  XYc(2), XY(2)

C ======================================================================
      If(XY(1).GT.XYc(1)) Then
          I = 2
      Else
          I = 1
      End if

      If(XY(2).GT.XYc(2)) Then
          J = 2
      Else
          J = 1
      End if

      Return
      End



C ======================================================================
      subroutine addQT(nQT, QT, h, XYc, XY, idx)
C ======================================================================
C  Arrangement of a new point XY, in the quadtree QT
C ======================================================================
      implicit none
      include 'lintrp.fd'

      Integer nQT, idx, QT(2, 2, *)
      Real*8  XYc(2, *), XY(2, *), h(*)

      Integer i,j, i1,i2, L, ip, ptr, dir(2)
C ======================================================================
      dir(1) =-1
      dir(2) = 1

      ip = 1
      L = 1
      Do while (.TRUE.)
         L = L + 1
         If(L.GT.MaxH) 
     &      Call errMesMBA(6201, 'addQT', 'No memory for quadtree')

         Call setIJ(XYc(1, ip), XY(1, idx), i, j)

         ptr = QT(i, j, ip)
         If(ptr.EQ.0) Then
            QT(i, j, ip) = idx
            Return
         Else If(ptr.LT.0) Then
            ip = -ptr
         Else If(ptr.GT.0) Then
            nQT = nQT + 1
            Do i1 = 1, 2
               Do i2 = 1, 2
                  QT(i1, i2, nQT) = 0
               End do
            End do

            QT(i, j, ip) = -nQT

            XYc(1, nQT) = XYc(1, ip) + dir(i) * h(L)
            XYc(2, nQT) = XYc(2, ip) + dir(j) * h(L)

            Call setIJ(XYc(1, nQT), XY(1, ptr), i, j)
  
            QT(i, j, nQT) = ptr
            ip = nQT
         Endif
      Enddo

      Return
      End



C ======================================================================
      subroutine order2(XYc, h, XY, ord, sqrd)
C ======================================================================
C  Ordering of subquads of quadtree QT by their distance from point XY
C ======================================================================
      implicit none
      Integer  ord(4)
      Real*8   XYc(2), XY(2), h, sqrd(4)

      Integer  i, j, ofs
      Real*8   dist(2), sqr(2, 2)

C ======================================================================
        ofs(i, j) = i + 2 * (j - 1)

      Do i = 1, 2
         dist(i)   = dabs(XY(i) - XYc(i))
          sqr(i, 1) = dist(i)**2
         sqr(i, 2) = ddim(dist(i), h)**2
      End do

      Call setIJ(XYc, XY, i, j)

        ord(1) = ofs(i, j)
        sqrd(1) = sqr(1, 1) + sqr(2, 1)

      If(dist(1).LT.dist(2)) Then
          ord(2) = ofs(3-i, j)
          ord(3) = ofs(i, 3-j)

          sqrd(2) = sqr(1, 1) + sqr(2, 2)
          sqrd(3) = sqr(1, 2) + sqr(2, 1)
      Else
          ord(2) = ofs(i, 3-j)
          ord(3) = ofs(3-i, j)

          sqrd(2) = sqr(1, 2) + sqr(2, 1)
          sqrd(3) = sqr(1, 1) + sqr(2, 2)
      Endif

        ord(4) = ofs(3-i, 3-j)
        sqrd(4) = sqr(1, 2) + sqr(2, 2)

      Return
      End



C ======================================================================
      Integer Function nearestQT(QT, XYc, xy, point, h)
C ======================================================================
C  This function returns the index of the quadtree point which is
C  the nearest to given one
C ======================================================================
      implicit none
      include 'lintrp.fd'

      Integer QT(4, *)
      REal*8  XYc(2, *), xy(2, *), point(2), h(*)

      Integer L, ip(MaxH), k(MaxH), ord(4, MaxH), ptr
      Real*8  sqrd(4, MaxH), sqrEdge, sqdist, min

C ======================================================================
        L   = 1
      ptr =-1
        min = 2.0
      nearestQT = 0

      Do while( .TRUE. )
         If(ptr.LT.0) Then
            L = L + 1
            If(L.GT.MaxH) Call errMesMBA(6202, 'nearestQT', 
     &                                         'No memory for quadtree')

            ip(L) = -ptr
            Call order2(XYc(1,-ptr), h(L), point, ord(1,L), sqrd(1,L))
            k(L) = 1
         Else
            If(ptr.GT.0) Then
               sqdist = sqrEdge(point, xy(1, ptr))
               If(min.GT.sqdist) Then
                min = sqdist
                  nearestQT = ord(k(L),L) + 4*(ip(L) - 1)
               End if
            End if

            Do while( .TRUE. )
              k(L) = k(L) + 1
               If(k(L).LE.4) Then
                 If(sqrd(k(L), L).LT.min) Goto 100
               End if

              L = L - 1
               If(L.EQ.1) Goto 200
            End do
         Endif

 100     Continue
          ptr = QT(ord(k(L), L), ip(L))
      End do

 200  Continue
      Return
      End



C ======================================================================
      Integer Function searchQT(QT, XYc,vrt, nt,tri, nEP,IEP, XY,h, TOL)
C ======================================================================
C  Routine determines a triangle containing point XY with the relative 
C  tolerance TOL.
C ======================================================================
      implicit none
      include 'lintrp.fd'

      Integer QT(*), nt, tri(3, *), nEP(*), IEP(*)
      Real*8  XYc(2, *), vrt(2, *), XY(2), h(*), TOL

      Integer buf(2, MaxBuf), tbuf(MaxTrSSEl), nearestQT, MaxTrials
      Integer i,j,k,m,n, i1,i2, m1,m2, idx,idx2, ip, iT,mT, nbuf, nTs
      Logical encloseQT

C ======================================================================
      searchQT = 0
      MaxTrials = max(int(nt * RelativeTrials), 20)
        k = 0

      Do while( k.lt.MaxTrials )
         nbuf = 0

         ip = nearestQT(QT, XYc, vrt, XY, h)
         If(ip.EQ.0) Goto 1000

          idx = QT(ip)

         If(idx.EQ.1) Then
            i1 = 1
         Else
            i1 = nEP(idx - 1) + 1
         End if
         i2 = nEP(idx)

         Do 100 i = i1, i2
            iT = IEP(i)
            Call findSE(nbuf, tbuf, iT, nTs)
            If(nTs.GT.0) Goto 100

            If(encloseQT(XY, vrt, tri(1,iT), TOL)) Then
               searchQT = iT
               Goto 1000
            Endif
 
            If(nbuf.EQ.MaxTrSSEl) Goto 9000

            nbuf = nbuf + 1
            tbuf(nbuf) = iT
 100     Continue

C ...  checking for neighbooring triangles   
         Do i = i1, i2
            iT = IEP(i)

            Do j = 1, 3
               idx2 = iabs(tri(j, iT))

               If(idx2.EQ.1) Then
                  m1 = 1
               Else
                  m1 = nEP(idx2 - 1) + 1
               End if
               m2 = nEP(idx2)

               Do 200 m = m1, m2
                  mT = IEP(m)
                  Call findSE(nbuf, tbuf, mT, nTs)
                  If(nTs.GT.0) Goto 200

                  If(encloseQT(XY, vrt, tri(1,mT), TOL)) Then
                     searchQT = mT
                     Goto 1000
                  Endif

                  If(nbuf.EQ.MaxTrSSEl) Goto 9000

                  nbuf = nbuf + 1
                  tbuf(nbuf) = mT
 200           Continue
               End do   
            End do

          k = k + 1
         If(k.GT.MaxBuf) Goto 9001

          buf(1, k) = ip
          buf(2, k) = idx
          QT(ip) = 0
      End do


 1000 Continue
      Do n = k, 1, -1
         QT(buf(1, n)) = buf(2, n)
      End do

      Return

c ... error messages
 9000 Continue
      Call errMesMBA(6203,'searchQT', 'MaxTrSSEl in lintrp.fd is short')
 9001 Continue
      Call errMesMBA(6203,'searchQT', 'MaxBuf in lintrp.fd is short')

      Return
      End



C ======================================================================
      Subroutine interpolate(nt,tri, nv,vrt, LDF, F, NXY,xy, G, 
     &                       QT, XYc, nEP, IEP, h, info, flagCvxCrvBnd)
C ======================================================================
C Routine builds the reversed map points -> triangles and interpolates
C the values of function F at the points xy
C ======================================================================
      implicit none

      Integer  nv, nt, LDF, NXY
      Integer  QT(2,2, *), tri(3, nt), nEP(*), IEP(*), info(3)
      Real*8   XYc(2, *), vrt(2, nv), xy(2, NXY)
      Real*8   F(LDF, nv), G(LDF, NXY), h(*)
      Logical  flagCvxCrvBnd

c Tolerance for checking whether the point belongs to an element
      Real*8    PREC
      Parameter(PREC = 1D-6)
c Number of possible relaxations of the above tolerance (6 orders of 10, here)
      Integer   MaxRelax
      Parameter(MaxRelax = 7)

      Integer  i,k, idx, iv1,iv2,iv3, searchQT, irel
      Logical  flagInit, flagInterpolate
      Real*8   tri_area0, a,b,c,d, TOL, xy1(2), xy2(2), xy3(2)
      EXTERNAL tri_area0

C ======================================================================
      flagInit        = info(1).GT.0
      flagInterpolate = info(2).GT.0

      If(flagInit) Then
         Call backReferences(nv, nt, 3, 3, tri, nEP, IEP)
      End if

      Do i = 1, NXY
         TOL = PREC
          irel = 0

 100     Continue
         idx = searchQT(QT,XYc, vrt, nt,tri, nEP, IEP, xy(1,i), h, TOL)

c If curvelinear boundaries do exist and/or the domain is convex, 
c then the node could be outside of the computational mesh
c In such a case, the relaxation of PREC is used
         If(idx.LE.0) Then
            If(flagCvxCrvBnd) Then
               If(TOL.LT.5D-2) Then
                  TOL = TOL * 10
               Else
                  TOL = TOL * 3.16
               End if

               If(irel.LT.MaxRelax) Then
                  irel = irel + 1
                  Goto 100
               Else
                  Call errMesMBA(6204, 'lintrp2D',
     .                 'Max number of relaxations for PREC is reached')
               End if
c If there are no curvelinear boundaries, then there is an error or 
c some inconsistency of parameters in lintrp.fd. 
            Else
               Call errMesMBA(6204, 'lintrp2D',
     .             'Failed to find element within PREC tolerance')
            End if
         End if

         info(3) = idx
         If(flagInterpolate) Then
            iv1 = tri(1, idx)
            iv2 = tri(2, idx)
            iv3 = tri(3, idx)

            Do k = 1, 2
               xy1(k) = vrt(k, iv1) - xy(k, i)
               xy2(k) = vrt(k, iv2) - xy(k, i)
               xy3(k) = vrt(k, iv3) - xy(k, i)
            End do

            a = tri_area0(xy2, xy3)
            b = tri_area0(xy3, xy1)
            c = tri_area0(xy1, xy2)
            d = a + b + c

            If(d.EQ.0D0) Then
             a = 1D0
             d = 1D0
            End if

            Do k = 1, LDF
               G(k, i) = (a * F(k,iv1) + b * F(k, iv2) 
     &                                 + c * F(k, iv3)) / d
            End do
         End if
      End do

      Return
      End



C ======================================================================
      Logical Function encloseQT(XY, vrt, TRI, tol)
C ======================================================================
C Routine determines if the point XY belongs to the triangle TRI
C ======================================================================
      implicit none
      Real*8   XY(2), vrt(2, *), tol
      Integer  TRI(3)

      Integer i, j, k
      Real*8  a, b, c, alph, beta, x(3), y(3), F, frac

C ======================================================================
      F(alph, beta) = a * alph + b * beta + c

      Do i = 1, 3
        x(i) = vrt(1, TRI(i))
        y(i) = vrt(2, TRI(i))
      End do

      k = 3
      j = 2
      Do i = 1, 3
         a = y(j) - y(i)
         b = x(i) - x(j)
         c = x(j) * y(i) - x(i) * y(j)

         frac = F(XY(1), XY(2)) / F(x(k), y(k))

         If(frac.LE.-tol) Then
            encloseQT = .FALSE.
            Return
         Endif
         j = k
         k = i
      End do

      encloseQT = .TRUE.

      Return
      End



C ======================================================================
      Double Precision Function sizeQT(point, iW, rW)
C ======================================================================
      include 'lintrp.fd'
C ======================================================================
C  The function returns the size of the quadtree cell which contains
C  the given point with coords point(2).
C
C  PARAMETERS:
C     point(2) - user given point inside the unit square (0,1)^2
C
C     rW(*)  - real*8  working memory of LINTRP2D
C     iW(*)  - integer working memory of LINTRP2D
C
C  REMARK 1. This function is the wrapping for function SizeHost.
C
C ======================================================================
      Real*8   point(2), rW(*)
      Integer  iW(*)

      Real*8   SizeHost

      iQT = iW(1)
      sizeQT = SizeHost(iW(iQT), rW(MaxH + 1), point, rW)

      Return
      End



C ======================================================================
      Real*8  Function SizeHost(QT, XYc, point, h)
C ======================================================================
      implicit none
      include 'lintrp.fd'

      Integer QT(4, *)
      Real*8  XYc(2, *), point(2), h(*)

      Integer L, ip(MaxH), ord(4, MaxH), ptr
      Real*8  sqrd(4, MaxH)

C ======================================================================
      L   = 1
      ptr = -1
      SizeHost = 1D0

      Do while( .TRUE. )
         If(ptr.LT.0) Then
            L = L + 1
            If(L.GT.MaxH) Call errMesMBA(1009, 'SizeHost', 
     &                            'local parameter MaxH is small')

            ip(L) = -ptr
            Call order2(XYc(1,-ptr), h(L), point, ord(1,L), sqrd(1,L))
         Else
            SizeHost = h(L)
            Return
         Endif

         ptr = QT(ord(1, L), ip(L))
      End do

      Call errMesMBA(6204, 'SizeHost', 'Host cell was not found')

      Return
      End


C*************************************************************************
      double precision function detf(i,j,ii,ip,xy,vrt)
C*************************************************************************
      implicit none
      integer i,j,ii,ip(*)
      double precision xy(2,*),vrt(2,*)

        detf = (vrt(1, ip(i)) - xy(1,ii))*(vrt(2, ip(j)) - xy(2,ii))
     &       - (vrt(2, ip(i)) - xy(2,ii))*(vrt(1, ip(j)) - xy(1,ii))
      return
      end


