C ======================================================================
      Subroutine makQ(nP, XYP, 
     &                nE, nEv, IPE, fixE,
     &                nEStar, hStar, HesP, detG, qE,
     &                flagFixShape)
C ======================================================================
C Routine computes hStar and Quality of elements in the initial mesh. 
C ======================================================================
      implicit none

      Integer  nP, nE, nEv, IPE(3, *), fixE(*)
      Real*8  XYP(2, *)

      Integer  nEStar
      Real*8   hStar, HesP(3, *), detG(*), qE(*)
      Integer flagFixShape


C LOCAL VARIABLES
      Integer  ip(4)
      Real*8   HesAvg(3), detAvg, dsum, tri_area2, VStar, Vn

      Integer  i,n, i1,i2,i3, iP1,iP2,iP3

C ======================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

      Do n = 1, nP
         Call calDet(HesP(1, n), detG(n))
      End do

      If(flagFixShape.LE.0) Then
         VStar = 0D0
         Do n = 1, nE
            dsum = 0D0

            Do i1 = 1, 3 
               i2 = ip(i1 + 1)
               iP1 = IPE(i1, n)
               iP2 = IPE(i2, n)

               Do i = 1, 3
                  HesAvg(i) = (HesP(i, iP1) + HesP(i, iP2)) / 2
               End do

               Call calDet(HesAvg, detAvg)

               dsum = dsum + dsqrt(detAvg)
            End do
  
            i3  = ip(i2 + 1)
            iP3 = IPE(i3, n)

            Vn = tri_area2(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
            VStar = VStar + dabs(Vn) * dsum / 6
         End do

         hStar = dsqrt(VStar / nEStar * 4D0 / dsqrt(3D0))
      Else
         hStar = -1D0
      End if

      Do n = 1, nE
         i1 = IPE(1, n)
         i2 = IPE(2, n)
         i3 = IPE(3, n)

         Call calQE(
     &        HesP(1, i1), XYP(1, i1),
     &        HesP(1, i2), XYP(1, i2),
     &        HesP(1, i3), XYP(1, i3),
     &        hStar, qE(n))
      End do


c ... set quality of fixed elements to 1
      Do n = 1, nEv
         qE(fixE(n)) = 1D0
      End do

      Return
      End



C ======================================================================
      Subroutine calDet(HesP, detG)
C ======================================================================
C Routine computes determinat det(H) of matrix H where
C       | HesP(1)  HesP(3) |
C   H = |                  |
C       | HesP(3)  HesP(2) |  
C ======================================================================
      Real*8  HesP(3), detG


C ======================================================================
      detG = HesP(1) * HesP(2) - HesP(3) ** 2

      If(detG.LE.0D0) Then
         Call SpectralModule(HesP, detG)
      End if

      Return
      End



C ======================================================================
      Subroutine calQE(Hes1, xy1, Hes2, xy2, Hes3, xy3, hStar, qE)
C ======================================================================
C Routines computes quality of triangle {xy1, xy2, xy3} assuming
C that the metric is linear.
C
C *** Remarks:
C        1. round-off errors require to use |Lk|
C        2. moved division by 2 in HesAvg in calculation of Lk 
C ======================================================================
      Real*8 Hes1(3), xy1(2)
      Real*8 Hes2(3), xy2(2)
      Real*8 Hes3(3), xy3(2)
      Real*8 hStar, qE

C LOCAL VARIABLES
      Real*8 HesAvg(3), det12, det13, det23, dsum
      Real*8 Pk, Vk, Lk
      Real*8 F, x1, y1, x2, y2

C group (FUnction)
      F(x1) = (x1 * (2D0 - x1)) ** 3

C ======================================================================
      Do i = 1, 3
         HesAvg(i) = Hes2(i) + Hes3(i)
      End do
      Call calDet(HesAvg, det23)

      x1 = xy3(1) - xy2(1)
      y1 = xy3(2) - xy2(2)
      Lk = dabs(HesAvg(1) * x1 * x1 + HesAvg(2) * y1 * y1 +
     &                            2 * HesAvg(3) * x1 * y1) / 2
      Pk = dsqrt(Lk)

      Do i = 1, 3
         HesAvg(i) = Hes1(i) + Hes2(i)
      End do
      Call calDet(HesAvg, det12)

      x1 = xy1(1) - xy2(1)
      y1 = xy1(2) - xy2(2)
      Lk = dabs(HesAvg(1) * x1 * x1 + HesAvg(2) * y1 * y1 +
     &                            2 * HesAvg(3) * x1 * y1) / 2
      Pk = Pk + dsqrt(Lk)

      Do i = 1, 3
         HesAvg(i) = Hes1(i) + Hes3(i)
      End do
      Call calDet(HesAvg, det13)

      x2 = xy1(1) - xy3(1)
      y2 = xy1(2) - xy3(2)
      Lk = dabs(HesAvg(1) * x2 * x2 + HesAvg(2) * y2 * y2 +
     &                            2 * HesAvg(3) * x2 * y2) / 2
      Pk = Pk + dsqrt(Lk)

      dsum = (dsqrt(det12) + dsqrt(det23) + dsqrt(det13)) / 2

c     Vk = dabs(x1 * y2 - y1 * x2) * 5D-1 * dsum / 3
c     qE = 20.784619D0 * Vk / (Pk ** 2) * F(x1)

      Vk = dabs(x1 * y2 - y1 * x2) * dsum

      qE = 3.46410316666D0 * Vk / (Pk ** 2)

      If(hStar.GT.0D0) Then
         x1 = Pk / (3 * hStar)
         If(x1.GT.1D0) x1 = 1D0 / x1

         qE = qE * F(x1)
      End if

      Return
      End



C ======================================================================
      Subroutine calQF(
C ======================================================================
     &      Hes1, det1, xy1,
     &      Hes2, det2, xy2,
     &      Hes3, det3, xy3,
     &      hStar, iw)
C ======================================================================
C Routine orders mesh edges in assending of their quality 
C ======================================================================
      Real*8  Hes1(3), det1, xy1(2)
      Real*8  Hes2(3), det2, xy2(2)
      Real*8  Hes3(3), det3, xy3(2)
      Real*8  hStar

      Integer iw(3)

C group (Local variables)
      Real*8  HesMax(3), detMax
      Real*8  qF(3)
      Real*8  F, x1, y1

C group (Function)
      F(x1) = (x1 * (2D0 - x1)) ** 3

C ======================================================================
      detMax = det1
      Do i = 1, 3
         HesMax(i) = Hes1(i)
      End do

      If(det2.GT.detMax) Then
         detMax = det2
         Do i = 1, 3
            HesMax(i) = Hes2(i)
         End do
      End if

      If(det3.GT.detMax) Then
         detMax = det3
         Do i = 1, 3
            HesMax(i) = Hes3(i)
         End do
      End if

      x1 = xy1(1) - xy2(1)
      y1 = xy1(2) - xy2(2)
      x1 = dsqrt(HesMax(1) * x1 * x1 + HesMax(2) * y1 * y1 +
     &                             2 * HesMax(3) * x1 * y1) / hStar
      x1 = min(x1, 1D0 / x1)
      qF(1) = F(x1)


      x1 = xy3(1) - xy2(1)
      y1 = xy3(2) - xy2(2)
      x1 = dsqrt(HesMax(1) * x1 * x1 + HesMax(2) * y1 * y1 +
     &                             2 * HesMax(3) * x1 * y1) / hStar
      x1 = min(x1, 1D0 / x1)
      qF(2) = F(x1)

      x1 = xy1(1) - xy3(1)
      y1 = xy1(2) - xy3(2)
      x1 = dsqrt(HesMax(1) * x1 * x1 + HesMax(2) * y1 * y1 +
     &                             2 * HesMax(3) * x1 * y1) / hStar
      x1 = min(x1, 1D0 / x1)
      qF(3) = F(x1)

      iw(1) = 1
      iw(2) = 2
      iw(3) = 3

      Do i = 1, 2
         iMin = i
         Do j = i + 1, 3
            If(qF(j).LT.qF(iMin)) iMin = j
         End do

         x1 = qF(i)
         qF(i) = qF(iMin)
         qF(iMin) = x1

         k = iw(i)
         iw(i) = iw(iMin)
         iw(iMin) = k
      End do

      Return
      End



C ======================================================================
      Real*8 Function avgQ(nE, qE, L1E, L2E)
C ======================================================================
      Real*8  qE(*)
      Integer L2E(*), L1E(2, *) 

      avgQ = 0D0

      iE = L2E(1)
      Do n = 1, nE
         avgQ = avgQ + qE(iE)
         iE = L1E(2, iE)
      End do

      avgQ = avgQ / nE

      Return
      End

