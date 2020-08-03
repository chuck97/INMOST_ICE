c ======================================================================
      Subroutine GradUtri(xy1, xy2, xy3, u1, u2, u3, GRAD)
c ======================================================================
c Routine computes gradient of a linear function given by values u1, u2,
c and u3 at points xy1, xy2, and xy3
c ======================================================================
      implicit none

      Real*8  xy1(*), xy2(*), xy3(*), GRAD(*)
      Real*8  u1, u2, u3

      Real*8  v12(2), v13(2), VecMul, det
      Integer i 
c ======================================================================
      Do i = 1, 2
         v12(i) = xy1(i) - xy2(i)
         v13(i) = xy1(i) - xy3(i)
      End do

      det = VecMul(v12, v13)

      GRAD(1) = ((u1 - u2) * v13(2) - (u1 - u3) * v12(2)) / det 
      GRAD(2) =-((u1 - u2) * v13(1) - (u1 - u3) * v12(1)) / det 

      Return
      End



C ======================================================================
C  The routine takes an input matrix D and rotates it counter clockwise
c  by angle theta. The angle theta is measured in degrees.
C ======================================================================
      Subroutine rotateTensor(LDD, D, theta)
C ======================================================================
      implicit none

      Real*8    PI
      PARAMETER(PI = 3.14159265358979D0)

      Integer  LDD, k
      Real*8   D(LDD, 2), theta, c, s, a, b
C ======================================================================
      a = theta * PI / 180 
      c = dcos(a)
      s = dsin(a)

      Do k = 1, 2
         a = D(k, 1)
         b = D(k, 2)

         D(k, 1) = a * c - b * s
         D(k, 2) = a * s + b * c
      End do

      Do k = 1, 2
         a = D(1, k)
         b = D(2, k)

         D(1, k) = a * c - b * s
         D(2, k) = a * s + b * c
      End do

      Return
      End



C ======================================================================
      Subroutine order2D(L, nE, IPE)
C ======================================================================
C  The routine orders each column of 2D array
C ======================================================================
      implicit none

      Integer  L, nE, IPE(L, *)
      Integer  i, j, n, im, mv

C ======================================================================
      Do n = 1, nE
         Do i = 1, L
            im = i
            mv = IPE(i, n)
            Do j = i + 1, L
               If(mv.GT.IPE(j, n)) Then
                  im = j
                  mv = IPE(j, n)
               End if
            End do
            IPE(im, n) = IPE(i,  n)
            IPE(i,  n) = mv
         End do
      End do

      Return
      End


