C ======================================================================
c Example of user defined mesh size function
C ======================================================================
      double precision Function meshSize(XY)
C ======================================================================
      double precision XY(2)  ! point where to evaluate mesh size
      double precision R, m0, c0

c distance to the boundary of the unit square
      R = min(XY(1), XY(2), 1D0-XY(1), 1D0-XY(2))

c minimal mesh size
      m0 = 0.005  ! if  Coef = 3D-4  in function Ddiff (forlibfem.f)
c     m0 = 0.003  ! if  Coef = 2D-4  in function Ddiff (forlibfem.f)

c coarsening speed away from the boundary
      c0 = 0.13

      meshSize = dsqrt(m0 ** 2 + R**1.5d0 * c0**2)

      Return
      End

