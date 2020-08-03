C ======================================================================
      Real*8 Function tri_area(xy1, xy2, xy3)
C ======================================================================
C Routine computes the oriented triangle area
C ======================================================================
      Real*8 xy1(2), xy2(2), xy3(2)

      tri_area = ((xy1(1) - xy3(1)) * (xy2(2) - xy3(2)) -
     &            (xy1(2) - xy3(2)) * (xy2(1) - xy3(1))) / 2
      Return
      End



C ======================================================================
      Real*8 Function tri_area2(xy1, xy2, xy3)
C ======================================================================
C Routine computes twice the triangle area
C ======================================================================
      Real*8 xy1(2), xy2(2), xy3(2)

      tri_area2 = (xy1(1) - xy3(1)) * (xy2(2) - xy3(2)) -
     &            (xy1(2) - xy3(2)) * (xy2(1) - xy3(1))
      Return
      End



C ======================================================================
      Integer Function tri_orient(xy1, xy2, xy3)
C ======================================================================
C Routine computes sign of the algebraic area
C ======================================================================
      Real*8 xy1(2), xy2(2), xy3(2)
      Real*8 a, b, c, d, ab, cd, area, eps

      a = xy1(1) - xy3(1)
      b = xy2(2) - xy3(2)
      c = xy1(2) - xy3(2)
      d = xy2(1) - xy3(1)

      ab = a * b
      cd = c * d
      area = ab - cd
      eps = 1D-16 * (dabs(a) + dabs(b) + dabs(c) + dabs(d) 
     &                       + dabs(ab) + dabs(cd))

      If(area.GT.eps) Then
         tri_orient = 1
      Else If(area.LT.-eps) Then
         tri_orient = -1
      Else
         tri_orient = 0
      End if

      Return
      End


C ======================================================================
      Real*8 Function tri_area0(xy1, xy2)
C ======================================================================
C Routine computes the oriented triangle area when xy3 = (0,0)
C ======================================================================
      Real*8 xy1(2), xy2(2)

      tri_area0 = (xy1(1) * xy2(2) - xy1(2) * xy2(1)) / 2
      Return
      End



c ======================================================================
      Real*8 Function tri_diameter(xy1, xy2, xy3)
c ======================================================================
c  Routines computes diameter of the triangle 
c ======================================================================
      implicit none
      Real*8   xy1(2), xy2(2), xy3(2)
      Real*8   sqrEdge, d1, d2, d3
c ======================================================================
      d1 = sqrEdge(xy1, xy2)
      d2 = sqrEdge(xy1, xy3)
      d3 = sqrEdge(xy2, xy3)

      d1 = max(d1, d2)
      d1 = max(d1, d3)

      tri_diameter = dsqrt(d1)

      Return
      End

 

c ======================================================================
      Subroutine tri_center(xy1, xy2, xy3, x, y)
c ======================================================================
c  Routines computes center of the triangle 
c ======================================================================
      implicit none
      Real*8   xy1(2), xy2(2), xy3(2), x, y
c ======================================================================
      x = (xy1(1) + xy2(1) + xy3(1)) / 3
      y = (xy1(2) + xy2(2) + xy3(2)) / 3

      Return
      End



c ======================================================================
      Subroutine tri_heights(xy1, xy2, xy3, H)
c ======================================================================
c  Routines computes three heights of the ttriangle
c ======================================================================
      implicit none
      Real*8   xy1(2), xy2(2), xy3(2), H(3)
      Real*8   tri_area, edge_length, v, e

c ======================================================================
      v = tri_area(xy1, xy2, xy3)
      v = dabs(v)

      e = edge_length(xy1, xy2)
      H(1) = 2 * v / e

      e = edge_length(xy2, xy3)
      H(2) = 2 * v / e

      e = edge_length(xy3, xy1)
      H(3) = 2 * v / e

      Return
      End



c ======================================================================
      Subroutine tri_angles(xy1, xy2, xy3, angles)
c ======================================================================
c  Routines computes three angles of the ttriangle
c ======================================================================
      implicit none
      Real*8   xy1(2), xy2(2), xy3(2), angles(3)
      Real*8   tri_area, edge_length, e12,e13,e23, c, pi

c ======================================================================
      pi = 3.14159265358979323846D0

      e12 = edge_length(xy1, xy2)
      e13 = edge_length(xy1, xy3)
      e23 = edge_length(xy2, xy3)

      c = (xy2(1) - xy1(1)) * (xy3(1) - xy1(1)) +
     &    (xy2(2) - xy1(2)) * (xy3(2) - xy1(2)) 

      angles(1) = c / (e12 * e13)
      angles(1) = dacos(angles(1)) / pi * 180D0

      c = (xy1(1) - xy2(1)) * (xy3(1) - xy2(1)) +
     &    (xy1(2) - xy2(2)) * (xy3(2) - xy2(2)) 

      angles(2) = c / (e12 * e23)
      angles(2) = dacos(angles(2)) / pi * 180D0

      angles(3) = 180 - angles(1) - angles(2)

      Return
      End

