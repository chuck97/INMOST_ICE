C ======================================================================
      Subroutine Delaunay(nP, nE, XYP, IPE, MaxWi, iW)
C ======================================================================
      implicit none
      include 'makS.fd'
C ======================================================================
C  The routine builds the Delaunay triangulation from the existing
C  triangution by swapping edges in pairs of triangles.
C
C  Data flow:
C      1. calculate map E->E
C      2. loop over triangles
C         2.1 check criterium for sum of opposite angles 
C         2.2 swap triangles
C         2.3 update the map
C
C  *** Remarks 
C ======================================================================
      Integer  nP, nE

      Real*8   XYP(2, *)
      Integer  IPE(3, *)

      Integer  MaxWi, iW(*)

      Integer  iref(4)
      Logical  check22, DelonePair, flagREPEAT, flagDELONE

      Integer  i1,i2,i3, j1,j2,j3, k, n
      Integer  iP1,iP2,iP3, jP1,jP2,jP3, iE1,iE2, jE2,jE3
      Integer  inEp, iIEP, iIEE, iEnd, iEmem, jEmem, kEmem
      Integer  nswap, nswapold, nloop
      Real*8   sa1, sa2

      DATA     iref /1,2,3,1/

C ======================================================================
      iIEE = 1
      inEP = iIEE + 3 * nE
      iIEP = inEP + nP
      iEnd = iIEP + 3 * nE

      If(iEnd.GT.MaxWi) Call errMesMBA(1001,
     &      'delaunay', 'MaxWi is too small')

      Call listE2E(nP, nE, IPE, iW(iIEE), iW(inEP), iw(iIEP))

c ... initialize loops
      nswapold = nE
      nloop = 0

c ... new loop
 100  flagREPEAT = .FALSE.
      nswap = 0

      Do n = 1, nE
         iEmem = iIEE + 3 * (n - 1) - 1

         Do 30 i1 = 1, 3
            iE1 = iW(iEmem + i1)
            If(iE1.LE.n) Goto 30

            i2 = iref(i1 + 1)
            i3 = iref(i2 + 1)

            iP1 = IPE(i1, n) 
            iP2 = IPE(i2, n) 
            iP3 = IPE(i3, n) 

            Do j1 = 1, 3
               j2 = iref(j1 + 1) 

               jP1 = IPE(j1, iE1)
               jP2 = IPE(j2, iE1)
               If(check22(iP1, iP2, jP1, jP2)) Then
                  j3 = iref(j2 + 1) 
                  jP3 = IPE(j3, iE1)
                  Goto 10 
               End if
            End do

 10         Continue
            flagDELONE = DelonePair(XYP(1, iP1), XYP(1, iP2), 
     &                              XYP(1, iP3), XYP(1, jP3), sa1)

            If(.NOT.flagDELONE) Then
               flagREPEAT = .TRUE.

               flagDELONE = DelonePair(XYP(1, iP3), XYP(1, jP3), 
     &                                 XYP(1, iP1), XYP(1, jP2), sa2)
               If(flagDELONE .OR. sa2.GT.sa1) Then
                  nswap = nswap + 1
                  iE2 = iW(iEmem + i2) 

                  jEmem = iIEE + 3 * (iE1 - 1) - 1
                  jE2 = iW(jEmem + j2)
                  jE3 = iW(jEmem + j3)
                  If(iP1.NE.jP1) Call swapii(jE3, jE2)

                  IPE(i1, n) = iP3
                  IPE(i2, n) = jP3
                  IPE(i3, n) = iP1

                  IPE(j1, iE1) = iP3
                  IPE(j2, iE1) = jP3
                  IPE(j3, iE1) = iP2

                  iW(iEmem + i2) = jE3
                  If(jE3.GT.0) Then
                     kEmem = iIEE + 3 * (jE3 - 1) - 1 
                     Do k = 1, 3
                        If(iW(kEmem + k).EQ.iE1) Then
                           iW(kEmem + k) = n
                           Goto 20
                        End if
                     End do

c                    Call draw_T(nP, 0, nE, XYP, iW, iW, IPE, 'fin.ps')
                     Call errMesMBA(6006, 
     &                   'delaunay.f', 'Cannot build Delaunay mesh')
                  End if

 20               iW(jEmem + j2) = jE2
                  iW(jEmem + j3) = iE2
                  If(iE2.GT.0) Then
                     kEmem = iIEE + 3 * (iE2 - 1) - 1
                     Do k = 1, 3
                        If(iW(kEmem + k).EQ.n) Then
                           iW(kEmem + k) = iE1 
                           Goto 30
                        End if
                     End do

c                    Call draw_T(nP, 0, nE, XYP, iW, iW, IPE, 'fin.ps')
                     Call errMesMBA(6006, 
     &                   'delaunay.f', 'Cannot build Delaunay mesh')
                  End if
               End if
            End if
 30      Continue
      End do

      If(flagREPEAT) Then
         nloop = nloop + 1
 
         If(nloop.GT.nE) Then
c          Call draw_T(nP, 0, nE, XYP, iW, iW, IPE, 'fin.ps')
           Call wrnMesMBA(6006, 
     &                   'delaunay.f', 'Cannot build Delaunay mesh')
           Return
         End if

         nswapold = nswap
         Goto 100
      End if
 
      Return
      End



C ======================================================================
      Logical Function DelonePair(xy1, xy2, xy3, xy4, sa)
C ======================================================================
C  This roune checks that the pair {1,2,3} and {1,2,4} is the
C  Delaunay pair. It returns .FALSE.otherwise.
C ======================================================================
      Real*8  xy1(2), xy2(2), xy3(2), xy4(2)

      Real*8  sa, sb, ca, cb

C ======================================================================
      DelonePair = .FALSE.

      ca = 0
      cb = 0

      Do i = 1, 2 
         ca = ca + (xy4(i) - xy1(i)) * (xy4(i) - xy2(i)) 
         cb = cb + (xy3(i) - xy1(i)) * (xy3(i) - xy2(i)) 
      End do
            
c ... Delaunay condition fails
      If(ca.LT.0 .AND. cb.LT.0) Goto 9000

c ... Delaunay condition true
      If(ca.GE.0 .AND. cb.GE.0) Then
         DelonePair = .TRUE.
         Goto 9000
      End if

c ... the rest of the formula
      sa = (xy4(1) - xy1(1)) * (xy4(2) - xy2(2)) 
     &   - (xy4(1) - xy2(1)) * (xy4(2) - xy1(2))

      sb = (xy3(1) - xy1(1)) * (xy3(2) - xy2(2)) 
     &   - (xy3(1) - xy2(1)) * (xy3(2) - xy1(2)) 

c ... checking when swapping is not possible
c     If(sa * sb.GE.0) Goto 9000

      sa = dabs(sa) * cb + ca * dabs(sb)
      If(sa.GE.0) DelonePair = .TRUE.
        
 9000 Return
      End



C ======================================================================
      Subroutine RandXY(xy1, xy2, xy3, xyc, r)
C ======================================================================
C Routine computes the center and radius of the circle 
c curcumscribed around triangle defined by three vertices.
C ======================================================================
      Real*8  xy1(2), xy2(2), xy3(2)
      Real*8  xyc(2), r

      Real*8  a, b, c, x1, y1, x2, y2, r1, r2
C ======================================================================
      x1 = xy1(1) - xy3(1)
      y1 = xy1(2) - xy3(2)

      x2 = xy2(1) - xy3(1)
      y2 = xy2(2) - xy3(2)

      r1 = x1 ** 2 + y1 ** 2
      r2 = x2 ** 2 + y2 ** 2

      a = x1 * y2 - y1 * x2
      b = r1 * y2 - r2 * y1
      c = r1 * x2 - r2 * x1  

      xyc(1) = b / (2 * a)
      xyc(2) =-c / (2 * a)

      r = dsqrt(xyc(1) ** 2 + xyc(2) ** 2)

      xyc(1) = xyc(1) + xy3(1)
      xyc(2) = xyc(2) + xy3(2)

      Return
      End


