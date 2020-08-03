c ======================================================================
      Subroutine GradU(xy1, xy2, xy3, u1, u2, u3, GRAD, area)
c ======================================================================
c Routine computes gradient of a linear function given
c by value u_i at points xy_i. 
c ======================================================================
      Real*8  xy1(*), xy2(*), xy3(*), GRAD(*)
      Real*8  u1, u2, u3, area

      Real*8  v12(2), v13(2), VecMul, det
c ======================================================================
      Do i = 1, 2
         v12(i) = xy1(i) - xy2(i)
         v13(i) = xy1(i) - xy3(i)
      End do

      det = VecMul(v12, v13)
      area = dabs(det) / 2

      GRAD(1) = ((u1 - u2) * v13(2) - (u1 - u3) * v12(2)) / det 
      GRAD(2) =-((u1 - u2) * v13(1) - (u1 - u3) * v12(1)) / det 

      Return
      End

