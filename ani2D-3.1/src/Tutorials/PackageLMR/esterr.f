C ======================================================================
      Subroutine triIntL8(U, nv, vrt, nt, tri, nx, error)
C ======================================================================
C Routine calculates the maximum norm of the interpolation
C error in triangles by spliting them in smaller subtriangles.
C ======================================================================
      Integer   Nstpan,Nstpsz
      Parameter(Nstpan=20, Nstpsz=20)
C ======================================================================
      Real*8   U(*), vrt(2, *), error(*)
      Integer  tri(3,*)

      Real*8   dx21,dy21,dx31,dy31, da,an,dl, lx,ly
      Real*8   xy(2), a,b,c,d, u1,u2,u3,uint, diff
      Integer  ip(3)
      
      Real*8   Func, detff
      EXTERNAL Func, detff

C ======================================================================
      nx = nt
      da = 1D0 / Nstpan
      dl = 1D0 / Nstpsz

      Do n = 1, nt
         error(n) = 0D0

         iP1 = tri(1, n)
         iP2 = tri(2, n)
         iP3 = tri(3, n)

         dx21 = vrt(1, iP2) - vrt(1, iP1)
         dy21 = vrt(2, iP2) - vrt(2, iP1)
         dx31 = vrt(1, iP3) - vrt(1, iP1)
         dy31 = vrt(2, iP3) - vrt(2, iP1)

         ip(1) = iP1
         ip(2) = iP2
         ip(3) = iP3

         u1 = U(iP1)
         u2 = U(iP2)
         u3 = U(iP3)

         Do j = 0, Nstpan
            an = j * da

            Do k = 0, Nstpsz
               lx = dx21*(1 - an) + dx31 * an
               ly = dy21*(1 - an) + dy31 * an

               xy(1) = vrt(1, iP1) + (k*dl) * lx 
               xy(2) = vrt(2, iP1) + (k*dl) * ly 

               a = detff(2, 3, ip,xy,vrt)
               b = detff(3, 1, ip,xy,vrt)
               c = detff(1, 2, ip,xy,vrt)
               d = a + b + c

               uint = (a * u1 + b * u2 + c * u3) / d
               diff = dabs(Func(xy) - uint)
               error(n) = max(error(n), diff)
            End do
         End do
      End do

      Return
      End


         
C ======================================================================
      Subroutine edgeIntL8(U, nv, vrt, nt, tri, nx, error)
C ======================================================================
C Routine calculates the maximum norm of the interpolation
C error on mesh edges by spltting edges into small subedges.
C ======================================================================
      Integer   Nsteps
      Parameter(Nsteps = 20)
C ======================================================================
      Real*8   U(*), vrt(2, *), error(3, *)
      Integer  tri(3, *)

      Real*8   xy(2), da, dx,dy, u1,u2,uint, diff
      Integer  iref(4)

      Real*8   Func
      EXTERNAL Func

      DATA iref/1,2,3,1/

C ======================================================================
      nx = 3 * nt
      da = 1D0 / Nsteps

      Do n = 1, nt
         Do i1 = 1, 3
            error(i1, n) = 0D0

            i2 = iref(i1 + 1)

            iP1 = tri(i1, n)
            iP2 = tri(i2, n)

            dx = (vrt(1, iP2) - vrt(1, iP1)) / Nsteps
            dy = (vrt(2, iP2) - vrt(2, iP1)) / Nsteps

            u1 = U(iP1)
            u2 = U(iP2)

            Do i = 1, Nsteps - 1
               xy(1) = vrt(1, iP1) + dx * i
               xy(2) = vrt(2, iP1) + dy * i

               a = da * i

               uint = a * u2 + (1 - a) * u1
               diff = dabs(Func(xy) - uint)
               error(i1, n) = max(error(i1, n), diff)
            End do
         End do
      End do

      Return
      End



C ======================================================================
      Real*8 Function detff(i,j,ip,xy,vrt)
C ======================================================================
C Routines will be replaces by tri_area(...)
C ======================================================================
      implicit none

      Integer i,j, ip(*)
      Real*8  xy(2), vrt(2,*)

C ======================================================================
      detff = (vrt(1, ip(i)) - xy(1))*(vrt(2, ip(j)) - xy(2))
     &      - (vrt(2, ip(i)) - xy(2))*(vrt(1, ip(j)) - xy(1))

      Return
      End


