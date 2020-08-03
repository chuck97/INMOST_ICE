C ======================================================================
      Subroutine DG2P1(nP, nE, XYP, IPE, IRE, fDG, fP1, MaxWi, iW, iERR)
C ======================================================================
      implicit integer (i)
      include 'makS.fd'
C ======================================================================
C  The routine maps a discontinuous piece-wise linear function with
C  d.o.f on edges onto a continuous piece-wise linear function with 
C  d.o.f. at vertices. We use the ZZ method for that.
C ======================================================================
C  *** Input
C         nP  - the number of mesh nodes
C         nE  - the number of triangles
C
C         XYP - coordinates of the mesh nodes
C         IPE - map: triangle -> vertices
C         IRE - map: triangle -> edges
C
C         fDG - discontinuous piece-wise linear function
C
C  *** Output: 
C         fP1 - continuous piece-wice linear function
C                    
C  *** Work memory: 
C          iW - integer array of size MaxWi
C
C  *** Remarks:
C          iW(1) returns the required memory
C ======================================================================
      Real*8   XYP(2, *)
      Integer  IPE(3, *), IRE(3, *), nP, nE

      Real*8   fDG(*), fP1(*)

      Integer  MaxWi, iW(*)

C ... external functions
      Real*8   LSvalue

C Local variables
      Real*8   xy(2, MaxS), values(MaxS)
      Integer  kbuf(MaxS), iref(4), iERR

      DATA  iref/1,2,3,1/

C ======================================================================
      iERR = 0

C ... distribute memory
      inEP = 0
      iIEP = inEP + nP
      iiW  = iIEP + 3 * nE

      If(MaxWi.LT.iiW) Then
         iERR = 1001
         Goto 9000
      End if


C ... create maps: vertix -> triangles
      Call backReferences(nP, nE, 3, 3, IPE, iW(inEP+1), iW(iIEP+1))


C ... interpolate at mesh points with least squares
      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = iW(inEP + n)

         If(i2 - i1 + 1.GT.MaxS) Call errMesMBA(1007,
     &      'DG2P1', 'Local parameter MaxS is small')

         kR = 0
         Do i = i1, i2
            iE = iW(iIEP + i)
   
            Do 100 j1 = 1, 3
               j2 = iref(j1 + 1) 

               iP1 = IPE(j1, iE)
               iP2 = IPE(j2, iE)

               If(iP1.NE.n .AND. iP2.NE.n) Goto 100

               iR1 = IRE(j1, iE) 

               Call findSE(kR, kbuf, iR1, m)
               If(m.GT.0) Goto 100 

               kR = kR + 1
               kbuf(kR) = iR1

               Do j = 1, 2
                  xy(j, kR) = (XYP(j, iP1) + XYP(j, iP2)) / 2
               End do
               values(kR) = fDG(iR1)
  100      Continue
         End do

         If(kR.EQ.2) Then
            fP1(n) = (values(1) + values(2)) / 2
         Else   
            fP1(n) = LSvalue(kR, xy, values, XYP(1, n))
         End if 
      End do

 9000 Continue
      iW(1) = iiW

      Return
      End



C ======================================================================
      Subroutine P02P1(nP, nE, XYP, IPE, fP0, fP1, MaxWi, iW, iERR)
C ======================================================================
      implicit none
      include 'makS.fd'
C ======================================================================
C  The routine maps a discontinuous piece-wise constant function 
C  with d.o.f in elements onto a continuous piece-wise linear 
C  function with d.o.f. at vertices. We use the ZZ method for that.
C
C  *** Remarks 
C      1. Extrapolation at boundary nodes adds some diffusion for
C         functions with sharp gradients near boundary.
C
C      2. The required memory is returned in iW(1) 
C ======================================================================
C  *** Input
C         nP  - the number of mesh nodes
C         nE  - the number of triangles
C
C         XYP - coordinates of the mesh nodes
C         IPE - map: triangle -> vertices
C
C         fP0 - discontinuous piece-wise constant function
C
C  *** Output: 
C         fP1 - continuous piece-wice linear function
C                    
C  *** Work memory: 
C          iW - integer array of size MaxWi, minimum is 3*(nP+nE)
C ======================================================================
      Real*8   XYP(2, *)
      Integer  IPE(3, *), nP, nE, iERR

      Real*8   fP0(*), fP1(*)

      Integer  MaxWi, iW(*)

C ... external functions
      Real*8   LSvalue

C ... local variables
      Integer  i,j,k, n, i1,i2, iP1,iP2,iP3, iE, inEP, iIEP, iiW
      Real*8   xy(2, MaxS), values(MaxS)

C ======================================================================
      iERR = 0

C ... distribute memory
      inEP = 0
      iIEP = inEP + nP
      iiW  = iIEP + 3 * nE

      If(MaxWi.LT.iiW) Then
         iERR  = 1001
         Goto 9000
      End if


C ... create maps: vertix -> triangles
      Call backReferences(nP, nE, 3, 3, IPE, iW(inEP+1), iW(iIEP+1))


C ... interpolate at mesh points with least squares
      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = iW(inEP + n)

         If(i2 - i1 + 1.GT.MaxS) Call errMesMBA(1007,
     &      'DG2P1', 'Local parameter MaxS is small')

         Do i = i1, i2
            iE = iW(iIEP + i)
   
            iP1 = IPE(1, iE)
            iP2 = IPE(2, iE)
            iP3 = IPE(3, iE)

            k = i - i1 + 1
            values(k) = fP0(iE)
            Do j = 1, 2
               xy(j, k) = (XYP(j,iP1) + XYP(j,iP2) + XYP(j,iP3)) / 3
            End do
         End do

         If(k.EQ.1) Then
            fP1(n) = values(1)
         Else if(k.EQ.2) Then
            fP1(n) = (values(1) + values(2)) / 2
         Else
            k = i2 - i1 + 1
            fP1(n) = LSvalue(k, xy, values, XYP(1, n))
         End if
      End do

 9000 Continue 
      iW(1) = iiW

      Return
      End



C ======================================================================
      Real*8 function LSvalue(k, xy, values, xy0)
C ======================================================================
C  This routine uses least square linear fit to points xy(2,*) and
C  evaluates the value of the linear function at point xy0.
C ======================================================================
      implicit none

      Integer  k
      Real*8   xy(2,*), values(*), xy0(2)

c (local variables)
      Real*8   A(3, 3), S(3), work(30)
      Integer  i, j, ipiv(3), info     
      
C ======================================================================
      Do i = 1, 3 
         Do j = 1, 3 
            A(i, j) = 0D0
         End do
         S(i) = 0D0
      End do

c ... generate the least squares matrix
      Do i = 1, k
         A(1,1) = A(1,1) + xy(1, i) * xy(1, i)
         A(1,2) = A(1,2) + xy(1, i) * xy(2, i)
         A(1,3) = A(1,3) + xy(1, i) 

         A(2,2) = A(2,2) + xy(2, i) * xy(2, i)
         A(2,3) = A(2,3) + xy(2, i) 
      End do
      A(3,3) = k

c ... generate the RHS
      do i = 1, k
         S(1) = S(1) + xy(1, i) * values(i)
         S(2) = S(2) + xy(2, i) * values(i)
         S(3) = S(3) + values(i)
      end do

      Call dsysv('U', 3, 1, A, 3, ipiv, S, 3, work, 30, info)

      If(info.NE.0) Call errMesMBA(3011, 'LSvalue',
     &                   'Error in Lapack routine dsysv')

c ... evaluate the linear function
      LSvalue = S(1) * xy0(1) + S(2) * xy0(2) + S(3)

      Return
      End

