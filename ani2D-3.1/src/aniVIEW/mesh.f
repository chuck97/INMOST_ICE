C ================================================================
      Subroutine graph(nv, vrt, nt, tri, fName)
C ================================================================
c  Make a simple ps-figure of the triangulation. The file name
C  MUST have extension '.ps'!
C ================================================================
      Integer     nv, nt
      Real*8      vrt(2,*)
      Integer     tri(3,*)
      Character*(100) fName

      Real*8       xymax(2), xymin(2), kx, scale
      Integer      i, k, l
      Character*30 fNameExt

C ================================================================
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

      Write(1,'(A)') '/Times-Roman findfont 12 scalefont setfont'
      Write(1,'(A)') ' 0 -10 moveto (mesh) show'

      Do k = 1,nt
         Write(1,*) (((vrt(i,tri(l,k))-xymin(i))*kx,i=1,2),l=1,3),' t'
      End do

      Write(1,*) ' showpage'
      Close(1)

      Return
      End



C ================================================================
      Subroutine graph_front(nv, vrt, nb, bnd, fName)
C ================================================================
c  Make a simple ps-figure of the triangulation. The file name
C  MUST have extension '.ps'!
C ================================================================
      Integer     nv, nb
      Real*8      vrt(2,*)
      Integer     bnd(2,*)
      Character*(*) fName

      Real*8       xymax(2), xymin(2), kx, scale
      Integer      i, k, l
      Character*30 fNameExt

C ================================================================
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
     &   ' /t{moveto lineto stroke}def'

      Write(1,'(A)') '/Times-Roman findfont 12 scalefont setfont'
      Write(1,'(A)') ' 0 -10 moveto (mesh) show'

      Do k = 1,nb
         Write(1,*) (((vrt(i,bnd(l,k))-xymin(i))*kx,i=1,2),l=1,2),' t'
      End do

      Write(1,*) ' showpage'
      Close(1)

      Return
      End



C ================================================================
      Subroutine graph_py(nv, vrt, nt, tri, fName)
C ================================================================
c  Make a simple ps-figure of the triangulation. The file name
C  MUST have extension '.ps'!
C ================================================================
      Integer     nv, nt
      Real*8      vrt(2,*)
      Integer     tri(3,*)
      Character*(*) fName

      Real*8       xymax(2), xymin(2), kx, scale
      Integer      i, k, l, is
      Character*30 fNameExt

C ================================================================
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
      Write(1,*) nt

      is = 10
      Do k = 1,nt
         Write(1,*) (((vrt(i,tri(l,k))-xymin(i))*kx+is,i=1,2),l=1,3)
      End do

      Write(1,*)
      Close(1)

      Return
      End



C ================================================================
