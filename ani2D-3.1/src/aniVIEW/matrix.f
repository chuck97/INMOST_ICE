C ================================================================
      Subroutine draw_matrix(N, IA, JA, fName)
C ================================================================
c  Make a simple ps-figure of the matrix sparcity.
C  MUST have extension '.ps'!
C ================================================================
      Integer     N, IA(*), JA(*)
      Character*(*) fName

      Real*8       xymax(2), xymin(2), size, kx, dot
      Integer      i, k, l
      Character*30 fNameExt

C ================================================================
      i = 1
      Do while(fName(i:i+2) .NE. '.ps')
         i = i + 1
      End do
      fNameExt = fName(1:i+2)

      size = 500.0
      dot = max(0.1, 2D1 / N)
      kx = size / N

      Open( UNIT=1, FILE=fName )
      Write(1,'(A)') '%!PS-Adobe-2.0 EPSF-2.0'
      Write(1,'(A)') '%%BoundingBox: 0 0  520 520'
      Write(1,'(A)') '%%EndComments'
      Write(1,'(A)') ' 10 10 translate 0 setlinewidth'
      Write(1,'(A)') 
     &   ' /t{newpath moveto lineto lineto closepath stroke}def'
      Write(1,'(A,F4.2,A)') '/fillC{gsave translate newpath 0 0 ', dot, 
     &                ' 0 360 arc closepath 0 setgray fill grestore}def'

      Write(1,'(A)') '/Times-Roman findfont 12 scalefont setfont'

      Write(1,'(A)') ' 0 -10 moveto (matrix) show'

      Do i = 1, N 
         y = size - i * kx
         Do j = IA(i), IA(i+1)-1
            x = JA(j) * kx
            Write(1,*) x, y,' fillC'
         End do
      End do

      Write(1,*) ' showpage'
      Close(1)

      Return
      End



C ================================================================
      Subroutine draw_matrix_tol(N, IA, JA, A, tol, fName)
C ================================================================
c  Make a simple ps-figure of the matrix sparcity.
C  MUST have extension '.ps'!
C ================================================================
      Integer     N, IA(*), JA(*)
      Real*8      A(*), tol
      Character*(*) fName

      Real*8       xymax(2), xymin(2), size, kx, dot
      Integer      i, k, l
      Character*30 fNameExt

C ================================================================
      i = 1
      Do while(fName(i:i+2) .NE. '.ps')
         i = i + 1
      End do
      fNameExt = fName(1:i+2)

      size = 500.0
      dot = max(0.1, 2D1 / N)
      kx = size / N

      Open( UNIT=1, FILE=fName )
      Write(1,'(A)') '%!PS-Adobe-2.0 EPSF-2.0'
      Write(1,'(A)') '%%BoundingBox: 0 0  520 520'
      Write(1,'(A)') '%%EndComments'
      Write(1,'(A)') ' 10 10 translate 0 setlinewidth'
      Write(1,'(A)') 
     &   ' /t{newpath moveto lineto lineto closepath stroke}def'
      Write(1,'(A,F4.2,A)') '/fillC{gsave translate newpath 0 0 ', dot, 
     &                ' 0 360 arc closepath 0 setgray fill grestore}def'

      Write(1,'(A)') '/Times-Roman findfont 12 scalefont setfont'

      Write(1,'(A)') ' 0 -10 moveto (matrix) show'

      Do i = 1, N 
         y = size - i * kx
         Do j = IA(i), IA(i+1)-1
            x = JA(j) * kx
            If(dabs(A(j)).GT.tol) Write(1,*) x, y,' fillC'
         End do
      End do

      Write(1,*) ' showpage'
      Close(1)

      Return
      End
