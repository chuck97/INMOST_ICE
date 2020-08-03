C ======================================================================
      Subroutine draw(nP, nF, nE, XYP, ICP, IPF, lbF, IPE, fName)
C ======================================================================
      include 'colors.fd'
C ======================================================================
C Routine converts the mesh into a PostScript file fName.
C The file name MUST have extension '.ps'. 
C ======================================================================
C group (M)
      Real*8  XYP(2, *)
      Integer ICP(*), IPF(2, *), lbF(*), IPE(3, *)

C group (File)
      Character*(*) fName

C group (Local variables)
      Real  mmTOpt, kx, ky
      Real  x1, y1, x2, y2, x3, y3, d
      Real  fxmax, fymax

      Real*8  edge_length
      Logical ifXnode

      Character*30 fNameExt

C ======================================================================
      i = 1
      Do while(fName(i:i+2) .NE. '.ps')
         i = i + 1
      End do
      fNameExt = fName(1:i+2)


      mmTOpt = 72.0 / 25.4
      Pi = 4.0*atan(1.0)
      rdTOdg = 180.0 / Pi
      xBig = 180.0
      yBig = 180.0

      fxmax = XYP(1, 1)
      fxmin = fxmax
      fymax = XYP(2, 1)
      fymin = fymax

      Do n = 2, nP
         fxmin = min(fxmin, real(XYP(1, n)))
         fxmax = max(fxmax, real(XYP(1, n)))
         fymin = min(fymin, real(XYP(2, n)))
         fymax = max(fymax, real(XYP(2, n)))
      End do

      kx = xBig / (fxmax-fxmin) * mmTOpt
      ky = yBig / (fymax-fymin) * mmTOpt

      If(kx.GT.ky) kx = ky

      ibx = (fxmax-fxmin) * kx + 10
      iby = (fymax-fymin) * kx + 10

      Open(30, file = fNameExt, status='UNKNOWN')

      Call headerPS(30, ibx, iby)

      Do n = 1, nP
         x1 = (XYP(1, n) - fxmin) * kx
         y1 = (XYP(2, n) - fymin) * kx
c        Write(30,*) x1,y1, ' m', x1, y1, n, ' intTOtext ctext'
      End do 

      Do 10 n = 1, nE
         If(IPE(1, n).EQ.0) Goto 10
         x1 = (XYP(1, IPE(1, n)) - fxmin) * kx
         y1 = (XYP(2, IPE(1, n)) - fymin) * kx

         x2 = (XYP(1, IPE(2, n)) - fxmin) * kx
         y2 = (XYP(2, IPE(2, n)) - fymin) * kx

         x3 = (XYP(1, IPE(3, n)) - fxmin) * kx
         y3 = (XYP(2, IPE(3, n)) - fymin) * kx

         Write(30,*) 'newpath'
         Write(30,*) x1,y1, ' m', x2,y2, ' l', x3,y3, ' l'
         Write(30,'(A)') ' closepath cBlack stroke'

c  ...   mark fix points 
         rad = fxmax - fxmin
         Do i = 1, 3
            iP1 = IPE(i, n)
            Do j = i + 1, 3
               iP2 = IPE(j, n)

               d = edge_length(XYP(1, iP1), XYP(1, iP2)) 
               rad = min(rad, d)
            End do
         End do
         rad = rad * kx / 4
         rad = min(rad, 1.0)

         Do i = 1, 3
            iP1 = IPE(i, n)
            If(ifXnode(ICP(iP1), jVnode)) Then
               x1 = (XYP(1, iP1) - fxmin) * kx
               y1 = (XYP(2, iP1) - fymin) * kx
               Write(30,*) x1,y1, rad, ' cRed fillCircle'
            End if
         End do
 10   Continue


c ... draw boundaries using thick lines
      Write(30,*) '0.4 slw'
      Do 20 n = 1, nF
         If(IPF(1, n).LE.0) Goto 20
         iP1 = IPF(1, n)
         iP2 = IPF(2, n)

         x1 = (XYP(1, iP1) - fxmin) * kx
         y1 = (XYP(2, iP1) - fymin) * kx

         x2 = (XYP(1, iP2) - fxmin) * kx
         y2 = (XYP(2, iP2) - fymin) * kx

         Write(30,*) x1,y1, ' m', x2,y2, ' l '
c        Write(30,*) x1,y1, x2,y2, ' emptyArrow '

         ic = lbF(n) / 5
         ic = max(0, lbF(n) - 5 * ic) 
         Write(30,'(A,I1,A)') 'c', ic, ' stroke'

         x1 = (x1 + x2) / 2
         y1 = (y1 + y2) / 2
c        Write(30,*) x1,y1, ' m', x1, y1, n, ' intTOtext ctext'
 20   Continue

      Write(30,'(A)') 'showpage '
      Close(30)
      Return
      End



C ======================================================================
      Subroutine draw_Q(nP, nE, XYP, IPE, qE, rQuality, fName)
C ======================================================================
      include 'colors.fd'
C ======================================================================
C Routine converts the mesh into a postcript file fName.
C ======================================================================
C group (M)
      Real*8  XYP(2, *), qE(*), rQuality
      Integer IPE(3, *)

C group (File)
      Character*(*) fName

C group (Local variables)
      Real  mmTOpt, kx, ky
      Real  x1, y1, x2, y2, x3, y3
      Real  fxmax, fymax
      Character*30 fNameExt

C ==========================================================
      i = 1
      Do while( fName(i:i+2) .NE. '.ps')
         i = i + 1
      End do
      fNameExt = fName(1:i+2)

      mmTOpt = 72.0 / 25.4
      Pi = 4.0*atan(1.0)
      rdTOdg = 180.0 / Pi
      xBig = 180.0
      yBig = 180.0

      fxmax = XYP(1, 1)
      fxmin = fxmax
      fymax = XYP(2, 1)
      fymin = fymax

      Do n = 2, nP
         fxmin = min(fxmin, real(XYP(1, n)))
         fxmax = max(fxmax, real(XYP(1, n)))
         fymin = min(fymin, real(XYP(2, n)))
         fymax = max(fymax, real(XYP(2, n)))
      End do

      kx = xBig / (fxmax-fxmin) * mmTOpt
      ky = yBig / (fymax-fymin) * mmTOpt

      If(kx.GT.ky) kx = ky

      ibx = (fxmax-fxmin) * kx + 10
      iby = (fymax-fymin) * kx + 10

      Open(30, file = fNameExt, status='UNKNOWN')

      Call headerPS(30, ibx, iby)

      Do 10 n = 1, nE
         If(IPE(1, n).EQ.0) Goto 10
         x1 = (XYP(1, IPE(1, n)) - fxmin) * kx
         y1 = (XYP(2, IPE(1, n)) - fymin) * kx

         x2 = (XYP(1, IPE(2, n)) - fxmin) * kx
         y2 = (XYP(2, IPE(2, n)) - fymin) * kx

         x3 = (XYP(1, IPE(3, n)) - fxmin) * kx
         y3 = (XYP(2, IPE(3, n)) - fymin) * kx

         Write(30,*) 'newpath'
         Write(30,*) x1,y1, ' m', x2,y2, ' l', x3,y3, ' l'
         Write(30,'(A)') ' closepath'

c  ...   color all elements
c        Write(30,*) qE(n), ' setgray fill stroke'

c  ...   color only low-quality elements
         If(qE(n).LE.rQuality) Then
            Write(30,*) ' 0 fill stroke'
         Else
            Write(30,*) ' cBlack stroke'
         End if
 10   Continue

      Write(30,'(A)') 'showpage '
      Close(30)
      Return
      End



C ======================================================================
      Subroutine draw_D(nP, nF, nE, XYP, ICP, IPF, lbF, IPE, fName)
C ======================================================================
      include 'colors.fd'
C ======================================================================
C Routine converts the Delone mesh into a postcript file fName.
C ======================================================================
C group (M)
      Real*8  XYP(2, *)
      Integer ICP(*), IPF(2, *), lbF(*), IPE(3, *)

C group (File)
      Character*(*) fName

C group (Local variables)
      Real  mmTOpt, kx, ky
      Real  x1, y1, x2, y2, x3, y3, rad
      Real  fxmax, fymax

      Real*8  edge_length
      Logical ifXnode

      Real*8  xyc(2), rOut

C ==========================================================
      mmTOpt = 72.0 / 25.4
      Pi = 4.0*atan(1.0)
      rdTOdg = 180.0 / Pi
      xBig = 150.0
      yBig = 150.0

      fxmax = XYP(1, 1)
      fxmin = fxmax
      fymax = XYP(2, 1)
      fymin = fymax

      Do n = 2, nP
         fxmin = min(fxmin, real(XYP(1, n)))
         fxmax = max(fxmax, real(XYP(1, n)))
         fymin = min(fymin, real(XYP(2, n)))
         fymax = max(fymax, real(XYP(2, n)))
      End do

      kx = xBig / (fxmax-fxmin) * mmTOpt
      ky = yBig / (fymax-fymin) * mmTOpt

      If(kx.GT.ky) kx = ky

      ibx = (fxmax-fxmin) * kx + 10
      iby = (fymax-fymin) * kx + 10

      Open(30, file = fName, status='UNKNOWN')

      Call headerPS(30, ibx, iby)

      Do n = 1, nP
         x1 = (XYP(1, n) - fxmin) * kx
         y1 = (XYP(2, n) - fymin) * kx
c        Write(30,*) x1,y1, ' m', x1, y1, n, ' intTOtext ctext'
      End do 

      Do 10 n = 1, nE
         If(IPE(1, n).EQ.0) Goto 10
         x1 = (XYP(1, IPE(1, n)) - fxmin) * kx
         y1 = (XYP(2, IPE(1, n)) - fymin) * kx

         x2 = (XYP(1, IPE(2, n)) - fxmin) * kx
         y2 = (XYP(2, IPE(2, n)) - fymin) * kx

         x3 = (XYP(1, IPE(3, n)) - fxmin) * kx
         y3 = (XYP(2, IPE(3, n)) - fymin) * kx

         Write(30,*) 'newpath'
         Write(30,*) x1,y1, ' m', x2,y2, ' l', x3,y3, ' l'
         Write(30,'(A)') ' closepath cBlack stroke'

c  ...   mark fix points 
         rad = fxmax - fxmin
         Do i = 1, 3
            iP1 = IPE(i, n)
            Do j = i + 1, 3
               iP2 = IPE(j, n)

               rad = min(rad, edge_length(XYP(1, iP1), XYP(1, iP2)))
            End do
         End do
         rad = rad * kx / 4
         rad = min(rad, 1.0)

         Do i = 1, 3
            iP1 = IPE(i, n)
            If(ifXnode(ICP(iP1), jVnode)) Then
               x1 = (XYP(1, iP1) - fxmin) * kx
               y1 = (XYP(2, iP1) - fymin) * kx
               Write(30,*) x1,y1, rad, ' cRed fillCircle'
            End if
         End do
 10   Continue


c ... draw boundaries using thick lines
      Write(30,*) 'blw'
      Do 20 n = 1, nF
         If(IPF(1, n).LE.0) Goto 20
         iP1 = IPF(1, n)
         iP2 = IPF(2, n)

         x1 = (XYP(1, iP1) - fxmin) * kx
         y1 = (XYP(2, iP1) - fymin) * kx

         x2 = (XYP(1, iP2) - fxmin) * kx
         y2 = (XYP(2, iP2) - fymin) * kx

         Write(30,*) x1,y1, ' m', x2,y2, ' l '
c        Write(30,*) x1,y1, x2,y2, ' emptyArrow '

         ic = lbF(n) / 5
         ic = max(0, lbF(n) - 5 * ic) 
         Write(30,'(A,I1,A)') 'c', ic, ' stroke'

         x1 = (x1 + x2) / 2
         y1 = (y1 + y2) / 2
c        Write(30,*) x1,y1, ' m', x1, y1, n, ' intTOtext ctext'
 20   Continue


c ... draw the circumcribed circle
      Do 30 n = 1, nE
         If(IPE(1, n).EQ.0) Goto 30
            
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call RandXY(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), xyc, rOut)

         x1  = (xyc(1) - fxmin) * kx
         y1  = (xyc(2) - fymin) * kx
         rad = rOut * kx

         Write(30,*) ' gsave', x1, y1, 
     &               ' translate np 0 0 ', rad, 
     &               ' 0 360 arc cp cGray stroke grestore'
 30   Continue

      Write(30,'(A)') 'showpage '
      Close(30)
      Return
      End



C ======================================================================
      Subroutine draw_T(nP, nF, nE, XYP, IEE, IPF, lbF, IPE, fName)
C ======================================================================
      include 'colors.fd'
C ======================================================================
C Routine converts the mesh into a postcript file fName.
C ======================================================================
C group (M)
      Real*8  XYP(2, *)
      Integer IEE(3, *), IPF(2, *), lbF(*), IPE(3, *)

C group (File)
      Character*(*) fName

C group (Local variables)
      Real  mmTOpt, kx, ky
      Real  x1, y1, x2, y2, x3, y3
      Real  fxmax, fymax

      Character*30 fNameExt

C ======================================================================
      i = 1
      Do while( fName(i:i+2) .NE. '.ps')
         i = i + 1
      End do
      fNameExt = fName(1:i+2)

      mmTOpt = 72.0 / 25.4
      Pi = 4.0*atan(1.0)
      rdTOdg = 180.0 / Pi
      xBig = 180.0
      yBig = 180.0

      fxmax = XYP(1, 1)
      fxmin = fxmax
      fymax = XYP(2, 1)
      fymin = fymax

      Do n = 2, nP
         fxmin = min(fxmin, real(XYP(1, n)))
         fxmax = max(fxmax, real(XYP(1, n)))
         fymin = min(fymin, real(XYP(2, n)))
         fymax = max(fymax, real(XYP(2, n)))
      End do

      kx = xBig / (fxmax-fxmin) * mmTOpt
      ky = yBig / (fymax-fymin) * mmTOpt

      If(kx.GT.ky) kx = ky

      ibx = (fxmax-fxmin) * kx + 10
      iby = (fymax-fymin) * kx + 10

      Open(30, file = fNameExt, status='UNKNOWN')
      Call headerPS(30, ibx, iby)

      Do n = 1, nP
         x1 = (XYP(1, n) - fxmin) * kx
         y1 = (XYP(2, n) - fymin) * kx
c        Write(30,*) x1,y1, ' m', x1, y1, n, ' intTOtext ctext'
      End do 

      Do 10 n = 1, nE
         If(IPE(1, n).EQ.0) Goto 10
         x1 = (XYP(1, IPE(1, n)) - fxmin) * kx
         y1 = (XYP(2, IPE(1, n)) - fymin) * kx

         x2 = (XYP(1, IPE(2, n)) - fxmin) * kx
         y2 = (XYP(2, IPE(2, n)) - fymin) * kx

         x3 = (XYP(1, IPE(3, n)) - fxmin) * kx
         y3 = (XYP(2, IPE(3, n)) - fymin) * kx

         Write(30,*) 'newpath'
         Write(30,*) x1,y1, ' m', x2,y2, ' l', x3,y3, ' l'
         Write(30,'(A)') ' closepath cBlack stroke'

         Do i = 1, 3
            j = i + 1
            If(j.GT.3) j = 1

            k = j + 1
            If(k.GT.3) k = 1

            iP1 = IPE(i, n)
            iP2 = IPE(j, n)
            iP3 = IPE(k, n)

            x1 = (XYP(1,iP1) + XYP(1,iP2)) * 0.45 + XYP(1,iP3) * 0.1
            y1 = (XYP(2,iP1) + XYP(2,iP2)) * 0.45 + XYP(2,iP3) * 0.1

            x1 = (x1 - fxmin) * kx
            y1 = (y1 - fxmin) * kx

            iE = IEE(i, n)
            Write(30,*) x1,y1, ' m', x1, y1, iE, ' intTOtext ctext'
         End do
 10   Continue


c ... draw boundaries using thick lines
      Write(30,*) '0.0 slw'
      Do 20 n = 1, nF
         If(IPF(1, n).LE.0) Goto 20
         iP1 = IPF(1, n)
         iP2 = IPF(2, n)

         x1 = (XYP(1, iP1) - fxmin) * kx
         y1 = (XYP(2, iP1) - fymin) * kx

         x2 = (XYP(1, iP2) - fxmin) * kx
         y2 = (XYP(2, iP2) - fymin) * kx

         Write(30,*) x1,y1, ' m', x2,y2, ' l '
c        Write(30,*) x1,y1, x2,y2, ' emptyArrow '

         ic = lbF(n) / 5
         ic = max(0, lbF(n) - 5 * ic) 
         Write(30,'(A,I1,A)') 'c', ic, ' stroke'

         x1 = (x1 + x2) / 2
         y1 = (y1 + y2) / 2
c        Write(30,*) x1,y1, ' m', x1, y1, n, ' intTOtext ctext'
 20   Continue

      Write(30,'(A)') 'showpage '
      Close(30)
      Return
      End



C ======================================================================
      Subroutine headerPS(io, ibx, iby)
C ======================================================================
C The routines writes a header for a PS file 
C ======================================================================
      Write(io,'(A)') '%!PS-Adobe-2.0 EPSF-2.0'
      Write(io,'(A,2I5)') '%%BoundingBox: 0 0 ', ibx, iby 
      Write(io,'(A)') '%%EndComments'

      Write(io, '(A)') '/m{moveto} def'
      Write(io, '(A)') '/l{lineto} def'

      Write(io, '(A)') '/cBlack{ 0.0  0.0  0.0 setrgbcolor}def'
      Write(io, '(A)' ) '/cGreen{ 0.0  1.0  0.0 setrgbcolor}def'
      Write(io, '(A)' ) '/cBlue{  0.0  0.0  1.0 setrgbcolor}def'
      Write(io, '(A)' ) '/cGelb{  1.0  1.0  0.0 setrgbcolor}def'
      Write(io, '(A)' ) '/cRed{   1.0  0.0  0.0 setrgbcolor}def'
      Write(io, '(A)' ) '/cViole{ 1.0  0.0  1.0 setrgbcolor}def'
      Write(io, '(A)' ) '/cGray{  0.5  0.5  0.5 setrgbcolor}def'


      Write(io, '(A)') '/c0{ cBlue } def'
      Write(io, '(A)') '/c1{ cGreen }def'
      Write(io, '(A)') '/c2{ cGelb } def'
      Write(io, '(A)') '/c3{ cRed }  def'
      Write(io, '(A)') '/c4{ cViole }def'

      Write(io, '(A)') '/slw{setlinewidth}def'

      Write(io, '(5(A,/))') 
     &   '%fillCircle : xCenter yCenter Radius -> NULL',
     &   '/fillCircle{',
     &   '  /rad exch def',
     &   '  gsave translate newpath 0 0 rad 0 360 arc closepath',
     &   '  fill grestore}def'


      Write(io,*) '0 slw'
      Write(io,*) '5 5 translate'
      
      Return
      End


