C ======================================================================
      Subroutine isolines(U, nv, vrt, nt, tri, nb, bnd, fName, ni)
C ======================================================================
C  Routine draws solution isolines. The number of isolines is ni.
C  The file fName must have extension '.ps'.
C ======================================================================
      implicit none

      Character*(*) fName
      Integer       nv, nt, nb, ni
      Real*8        U(*), vrt(2,*)
      Integer       tri(3,*), bnd(2,*)

      Real*8        xymax(2), xymin(2), umin, umax, ucur, x(2), y(2)
      Real*8        vi, kx, tet, scale
      Integer       i, j, k, is, icount, jp, kp
      Character*30  fNameExt

C ======================================================================
      i = 1
      Do while( fName(i:i+2) .NE. '.ps')
         i = i + 1
      End do
      fNameExt = fName(1:i+2)

      xymin(1) = vrt(1,1)
      xymin(2) = vrt(2,1)
      xymax(1) = vrt(1,1)
      xymax(2) = vrt(2,1)

      umin = u(1)
      umax = u(1)
      Do i = 2, nv
         xymin(1) = min(xymin(1), vrt(1,i))
         xymin(2) = min(xymin(2), vrt(2,i))

         xymax(1) = max(xymax(1), vrt(1,i))
         xymax(2) = max(xymax(2), vrt(2,i))

         umax = max(umax, u(i))
         umin = min(umin, u(i))
      End do

      scale = max(xymax(1)-xymin(1), xymax(2)-xymin(2))
      kx = 500.0 / scale

      Open(UNIT=1, FILE=fNameExt )
      Write(1,'(A)') '%!PS-Adobe-2.0 EPSF-2.0'
      Write(1,'(A)') '%%BoundingBox: 0 0  520 520'
      Write(1,'(A)') '%%EndComments'
      Write(1,'(A)') ' 10 10 translate 0 setlinewidth'
      Write(1,'(A)') ' /v{moveto lineto stroke}def'

      Write(1,'(A)') '/Times-Roman findfont 12 scalefont setfont'
      Write(1,'(A)') ' 0 -10 moveto (solution isolines) show'

      Do i = 1, nt
         Do is = 1, ni
            icount = 0
            ucur = umin + (is-0.5)*(umax-umin)/(ni)

            Do j = 1, 3
               k = j - 1
               If(k.EQ.0) k = 3

               jp = tri(j, i)
               kp = tri(k, i)

               If(U(jp).EQ.U(kp) .AND. U(jp).EQ.ucur 
     &                           .AND. icount.LT.2) Then
                  x(1) = vrt(1, jp)
                  y(1) = vrt(2, jp)

                  x(2) = vrt(1, kp)
                  y(2) = vrt(2, kp)

                  vi = ucur
                  icount = 2
               End if

               If(((U(jp).LT.ucur .AND. ucur.LE.U(kp)) .OR.
     &             (U(jp).GE.ucur .AND. ucur.GT.U(kp))) .AND. 
     &            icount.LT.2) Then

                  icount = icount + 1

                  tet = (ucur - U(jp)) / (U(kp) - U(jp))
                  x(icount) = vrt(1,jp) + tet * (vrt(1,kp) - vrt(1,jp))
                  y(icount) = vrt(2,jp) + tet * (vrt(2,kp) - vrt(2,jp))

                  vi = ucur
               End if
            End do

            If(icount.NE.0) Then
              Write(1,'(4F7.1,A)')
     &           (x(1)-xymin(1)) * kx, (y(1)-xymin(2)) * kx,
     &           (x(2)-xymin(1)) * kx, (y(2)-xymin(2)) * kx, ' v'
            End if
         End do
      End do

c    draw boundary
      Do  i = 1, nb
          jp = bnd(1, i)
          kp = bnd(2, i)

          Write(1,'(4F7.1,A)') 
     &      (vrt(1,jp)-xymin(1)) * kx, (vrt(2,jp)-xymin(2)) * kx,
     &      (vrt(1,kp)-xymin(1)) * kx, (vrt(2,kp)-xymin(2)) * kx, ' v'
      End do

      Write(1,*) ' showpage'
      Close(1)

      Return
      End

