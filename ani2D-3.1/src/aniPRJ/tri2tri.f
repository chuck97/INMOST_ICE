C =====================================================================
      Subroutine tri2tri(p1, p2, p3, pp1, pp2, pp3, N, XYP)
C =====================================================================
      implicit none
      Integer  N
      Real*8   p1(2),p2(2),p3(2), pp1(2),pp2(2),pp3(2), XYP(2,6)

      Integer  i,j, i1,i2, k1,k2, N1,M, iref(4), next(9)
      Real*8   XY1(2, 3), XY2(2, 9), normal(2), dist(9), a
      DATA     iref/1,2,3,1/
C =====================================================================
      Do i = 1, 2
         XY1(i, 1) = p1(i)
         XY1(i, 2) = p2(i)
         XY1(i, 3) = p3(i)
      End do

      Do i = 1, 2
         XY2(i, 1) = pp1(i)
         XY2(i, 2) = pp2(i)
         XY2(i, 3) = pp3(i)
      End do

      next(1) = 2
      next(2) = 3
      next(3) = 1

      N  = 3
      M  = 3
      k1 = 1

      Do i1 = 1, 3
         i2 = iref(i1 + 1)

         Call edge_normal(XY1(1,i1), XY1(1,i2), normal)

         Do i = 1, N
            dist(k1) = (XY1(1,i1) - XY2(1,k1)) * normal(1) +
     &                 (XY1(2,i1) - XY2(2,k1)) * normal(2)
            If(dabs(dist(k1)).LT.1D-10) dist(k1) = 0
            k1 = next(k1)
         End do

         N1 = N
         Do i = 1, N1
            k2 = next(k1)

            If(dist(k1)*dist(k2).LT.0) Then
               N = N + 1
               M = M + 1 

               a = dist(k1) / dist(k2)
               Do j = 1, 2
                  XY2(j, M) = (XY2(j, k1) - a * XY2(j, k2)) / (1 - a)
               End do
               next(k1) = M
               next(M)  = k2
               dist(M)  = 0
            End if

            k1 = k2
         End do


         N1 = N
         Do i = 1, N1
            k2 = next(k1)
            If(dist(k2).GT.0) Then
               N = N - 1
               If(N.EQ.2) Return

               next(k1) = next(k2)
            Else
               k1 = next(k1)
            End if
         End do
      End do


c ... contigious arrays of vertices 
      Do i = 1, N
         Do j = 1, 2
            XYP(j, i) = XY2(j, k1)
         End do

         k1 = next(k1)
      End do

      Return
      End

