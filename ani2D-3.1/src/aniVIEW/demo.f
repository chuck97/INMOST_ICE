C ======================================================================
      Subroutine graph_demo(nv, vrt, nt, tri, fName, demo_message)
C ======================================================================
c  Make a simple ps-figure of the triangulation. The file name
C  MUST have extension '.ps'!
C ======================================================================
      Integer     nv, nt
      Real*8      vrt(2,*)
      Integer     tri(3,*)
      Character*(100) fName, demo_message

      Real*8       xymax(2), xymin(2), kx, scale
      Integer      i, k, l
      Character*30 fNameExt

C ======================================================================
      i = 1
      Do while(fName(i:i+2) .NE. '.ps')
         i = i + 1
      End do
      fNameExt = fName(1:i+2)

      xymin(1) = vrt(1,1)
      xymin(2) = vrt(2,1)
      xymax(1) = vrt(1,1)
      xymax(2) = vrt(2,1)
      Do i = 2, nv
         xymin(1) = min(xymin(1),vrt(1,i))
         xymin(2) = min(xymin(2),vrt(2,i))
         xymax(1) = max(xymax(1),vrt(1,i))
         xymax(2) = max(xymax(2),vrt(2,i))
      End do

      scale = max(xymax(1)-xymin(1),xymax(2)-xymin(2))
      kx = 500.0 / scale

      Open( UNIT=1, FILE=fName )
      Write(1,'(A)') '%!PS-Adobe-2.0 EPSF-2.0'
      Write(1,'(A)') '%%BoundingBox: 0 0  520 520'
      Write(1,'(A)') '%%EndComments'
      Write(1,'(A)') ' 10 10 translate 0 setlinewidth'
      Write(1,'(A)') 
     &   ' /t{newpath moveto lineto lineto closepath stroke}def'

      Do k = 1,nt
         Write(1,*) (((vrt(i,tri(l,k))-xymin(i))*kx,i=1,2),l=1,3),' t'
      End do

c ... demo part
      Call headerPS_demo(1)

      Write(1,'(3A)') '250 450 (', demo_message, ') clearANDctext'
      Write(1,*) ' showpage'
      Close(1)

      Return
      End



C ======================================================================
      Subroutine isolines_demo(U, nv, vrt, nt, tri, nb, bnd, fName, ni,
     &                         demo_message)
C ======================================================================
C  Routine draws solution isolines. The number of isolines is ni.
C  The file fName must have extension '.ps'.
C ================================================================
      implicit none

      Character*(*) fName, demo_message
      Integer       nv, nt, nb, ni
      Real*8        U(*), vrt(2,*)
      Integer       tri(3,*), bnd(2,*)

      Real*8        xymax(2), xymin(2), umin, umax, ucur, x(2), y(2)
      Real*8        vi, kx, tet, scale
      Integer       i, j, k, is, icount, jp, kp
      Character*30  fNameExt

C ================================================================
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

      Do i = 1, nt
         Do is = 1, ni
            icount = 0
            ucur = umin + (is-1)*(umax-umin)/(ni-1)

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

c ... demo part
      Call headerPS_demo(1)

      Write(1,'(3A)') '250 450 (', demo_message, ') clearANDctext'
      Write(1,*) ' showpage'
      Close(1)

      Return
      End



C ======================================================================
      Subroutine headerPS_demo(io)
C ======================================================================
C The routines writes a header for a PS demo file 
C ======================================================================
      Write(io,'(A)') '/l{lineto}def /m{moveto}def /s{l stroke}def'
      Write(io,'(A)') '/x{exch}def /rm{rmoveto}def'
      Write(io,'(A,/)') '/np{newpath}def /cp{closepath}def'

      Write(io,'(A)') '/fsize 27 def'
      Write(io,'(A,/)') '/Times-Roman findfont fsize scalefont setfont'

      Write(io,'(A)') '/ctext{3 1 roll m dup stringwidth fsize 0.33' 
      Write(io,'(A,/)') 'mul sub x 2 div 0 x sub x rm show}def'

      Write(io,'(A)') '/clearField{/yUp x def /xRight x def'
      Write(io,'(A)') '/yDown x def /xLeft x def'
      Write(io,'(A)') 'gsave np xLeft yDown m xRight yDown'
      Write(io,'(A)') 'l xRight yUp l xLeft yUp l cp 1 setgray fill'
      Write(io,'(A,/)') 'grestore}def'

      Write(io,'(A)') '/clearANDctext{3 copy stringwidth /dY fsize 2'
      Write(io,'(A)') 'div def pop /dX x 2 div def /yDown x def /xLeft'
      Write(io,'(A)') 'x def xLeft dX sub yDown dY sub xLeft dX add'
      Write(io,'(A,/)') 'yDown dY add clearField ctext}def'

      Write(io,'(A)') '0. 0. 1. setrgbcolor'

      Return
      End

