C ======================================================================
      Subroutine MetricBV(nP, nE, XYP, IPE, Metric, rW, iW)
C ======================================================================
C Routine calculates metric jumps across mesh edges. 
C ======================================================================
      Integer IPE(3, *), iW(*)
      Real*8  XYP(2, *), Metric(3, *), rW(*)

      Integer iref(4)
      Real*8  x, y, vb, vn, vol, tri_area, det, lenght

      DATA    iref/1,2,3,1/
C ======================================================================
      iIRE = 1
      inEP = iIRE + 3 * nE
      iIEP = inEP + nP
      iIRE = iIEP + 3 * nE
      iMrk = iIRE + 3 * nE

      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

c ... mark interior edges
      Do n = 1, nR
         iW(iMrk + n - 1) = -1
      End do

      Do n = 1, nE
         Do i = 1, 3
            iRt = iW(iIRE + 3*(n-1) + (i-1))
            k = iMrk + iRt - 1
            iW(k) = iW(k) + 1
         End do
      End do


c ... calculate jump for interior edges (Mark=1)
      nRi = 0
      Do n = 1, nR
         rW(n) = 0D0
      End do

      lenght = 1D+24
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         vol = tri_area(XYP(1,iP1), XYP(1,iP2), XYP(1,iP3))
         Call calDet(Metric(1,n), det)
 
         vol = dabs(vol) * dsqrt(det)

         Do i1 = 1, 3
            i2 = iref(i1 + 1)

            iRt = iW(iIRE + 3*(n-1) + (i1-1))
            iVt = iRt + nR
            k = iMrk + iRt - 1

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            x = XYP(1, iP1) - XYP(1, iP2)
            y = XYP(2, iP1) - XYP(2, iP2)

            vb = Metric(1,n) * x*x + Metric(2,n) * y*y 
     &                             + Metric(3,n) * x*y
            vb = dsqrt(vb) 
            lenght = min(lenght, dsqrt(x*x + y*y))

            If(iW(k).EQ.1) Then 
               rW(iRt) = vb
               rW(iVt) = vol
               iW(k) = iW(k) + 1
            Else If(iW(k).EQ.2) Then 
               vol = vol + rW(iVt)

               vn = rW(iRt)
               rW(iRt) = ((vn - vb) ** 2) * vol
               rW(iVt) = ((vn + vb) ** 2) * vol

               nRi = nRi + 1
            End if
         End do
      End do

      vb = 0D0
      vn = 0D0
      Do n = 1, nR
         vb = vb + rW(n)
         vn = vn + rW(n + nR)
      End do

      Write(*,'(2(A,E12.6))') '   Metric norm =', dsqrt(vb / vn),
     &                        '   minimal edge lenght =', lenght

      Return
      End


