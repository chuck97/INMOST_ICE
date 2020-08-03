C ======================================================================
      Real*8 Function sqrEdge(xy1, xy2)
C ======================================================================
C Routine computes sqruare distance between two points
C ======================================================================
      Real*8 xy1(2), xy2(2)

      sqrEdge = (xy1(1) - xy2(1)) ** 2 + (xy1(2) - xy2(2)) ** 2
      Return
      End



C ======================================================================
      Real*8 Function edge_length(xy1, xy2)
C ======================================================================
C Routine computes distance between two points
C ======================================================================
      Real*8 xy1(2), xy2(2)

      edge_length = dsqrt((xy1(1) - xy2(1)) ** 2 +
     &                    (xy1(2) - xy2(2)) ** 2)
      Return
      End



C ======================================================================
      Subroutine edge_normal(xy1, xy2, xyn)
C ======================================================================
C Routines compute a normal vector to the edge {xy1, xy2}
C ======================================================================
      Real*8 xy1(2), xy2(2), xyn(2)
      Real*8 x, y, d

      x = xy2(1) - xy1(1)
      y = xy2(2) - xy1(2)

      d = dsqrt(x * x + y * y)
 
      xyn(1) = -y / d
      xyn(2) =  x / d

      Return
      End


C ======================================================================
      Subroutine extNormal(xy1, xy2, xy3, xyn)
C ======================================================================
C Routines compute external normal vector to the edge {xy1, xy2}
C of triangle {xy1, xy2, xy3}
C ======================================================================
      Real*8 xy1(2), xy2(2), xy3(2), xyn(2)
      Real*8 x, y, d

      x = xy2(1) - xy1(1)
      y = xy2(2) - xy1(2)

      d = dsqrt(x * x + y * y)
 
      xyn(1) = -y / d
      xyn(2) =  x / d

c ... orientation
      x = xy3(1) - xy1(1)
      y = xy3(2) - xy1(2)

      If( x*xyn(1) + y*xyn(2).GT.0D0 ) Then
         xyn(1) = -xyn(1)
         xyn(2) = -xyn(2)
      End if 

      Return
      End









