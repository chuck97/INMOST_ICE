C ======================================================================
      Subroutine smoothingMesh(nP, nE, XYP, IPE, MaxWi, iW)
C ======================================================================
      include 'makS.fd'
C ======================================================================
C Routines smoothes the mesh with Laplacian smoother.
C ======================================================================
      Real*8  XYP(2, *)
      Integer IPE(3, *)

      Integer MaxWi
      Integer iW(*)

C group (Local variables)
      Real*8  XYPs(2)
      Logical flagOrient
      Character*80 message
C ======================================================================
      inEP = 0
      iIEP = inEP + nP
      iMrk = iIEP + 3 * nE
      iiOs = iMrk + nP
      iEnd = iiOs + MaxS

      If(iEnd.GT.MaxWi) Then
         iERR = 1001
         Write(message,'(A,I10)')
     &        'The approximate size of iW is ', iEnd
         Call errMesMBA(iERR, 'smoothingMesh', message)
      End if


C ... unmarked points
      Do n = 1, nP
         iW(iMrk + n) = 0
      End do

C ... create a map P->E
      Call backReferences(nP, nE, 3, 3, IPE, iW(inEP + 1), iW(iIEP + 1))

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = iW(inEP + n)
  
c  ...  calculate mid-point of the superelement
         Do j = 1, 2
            XYPs(j) = 0D0
         End do

         m = 0
         Do i = i1, i2
            iE = iW(iIEP + i)

            Do k = 1, 3
               iP = IPE(k, iE)
               If(iP.NE.n .AND. iW(iMrk + iP).EQ.0) Then
                  iW(iMrk + iP) = 1

                  m = m + 1
                  Do j = 1, 2
                     XYPs(j) = XYPs(j) + XYP(j, iP)
                  End do
               End if
            End do 
         End do

         Do j = 1, 2
            XYPs(j) = XYPs(j) / m
         End do

c  ...  check for interior point
         lE = i2 - i1 + 1
         If(m.NE.lE) Goto 100

c  ...  check orientation
         Call calSO(XYP, IPE, lE, iW(iIEP + i1), iW(iiOs + 1))
         Call chkSO(n, n, XYPs, XYP, IPE, 
     &                    lE, iW(iIEP + i1), iW(iiOs + 1), flagOrient)
         If(.NOT.flagOrient) Goto 100

c  ...  update the mesh
         Do j = 1, 2
            XYP(j, n) = XYPs(j)
         End do

c  ...  upmark points
  100    Continue
         Do i = i1, i2
            iE = iW(iIEP + i)

            Do k = 1, 3
               iP = IPE(k, iE)
               iW(iMrk + iP) = 0
            End do 
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine smoothingP1(
C ======================================================================
C Routine smoothes the piecewise linear function Sol defined at mesh
C points. The value of the smoothed function at a mesh point P is equal
C to the integral average over the superelement associated with point P. 
C This method is known as the Steklov smoothing.
C ======================================================================
     &           nP, nE, XYP, IPE, Sol,
     &           MaxWr, MaxWi, rW, iW)
C ======================================================================
      Real*8  XYP(2, *), Sol(*)
      Integer IPE(3, *)

      Integer MaxWi, MaxWr
      Integer iW(*)
      Real*8  rW(*)

c group (Local variables)
      Real*8  v, tri_area
      Character*80 message

C ======================================================================
      inEP = 0
      iIEP = inEP + nP
      iEnd = iIEP + 3 * nE

      If(iEnd.GT.MaxWi) Then
         iERR = 1001
         Write(message,'(A,I10)')
     &        'The approximate size of iW is ', iEnd
         Call errMesMBA(iERR, 'smoothingP1', message)
      End if

      iVol  = 0
      iSup = iVol + nE
      iSol = iSup + nP
      iEnd  = iSol + nP

      If(iEnd.GT.MaxWr) Then
         iERR = 1002
         Write(message,'(A,I10)')
     &        'The approximate size of rW is ', iEnd
         Call errMesMBA(iERR, 'smoothingP1', message)
      End if


C ... creating an auxiliary structure
      Call backReferences(nP, nE, 3, 3, IPE, iW(inEP + 1), iW(iIEP + 1))


C ... computing volumes of elements and superelements
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         v = tri_area(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))

         rW(iVol + n) = dabs(v)
      End do

      i2 = 0
      Do n = 1, nP
         i1 = i2 + 1
         i2 = iW(inEP + n)

         v = 0D0
         Do i = i1, i2
            iE = iW(iIEP + i)
            v = v + rW(iVol + iE)
         End do
         rW(iSup + n) = v
      End do
      

C ... itegrating the piecewise linear function SOL
      Do n = 1, nP
         rW(iSol + n) = 0D0
      End do

      Do n = 1, nE
         v = 0D0
         Do i = 1, 3
            iP1 = IPE(i, n)
            v = v + Sol(iP1) 
         End do
         v = v * rW(iVol + n) / 3 

         Do i = 1, 3
            iP1 = IPE(i, n)
            rW(iSol + iP1) = rW(iSol + iP1) + v  / rW(iSup + iP1)
         End do
      End do
      
      Do n = 1, nP
         Sol(n) = rW(iSol + n)
      End do

      Return
      End


