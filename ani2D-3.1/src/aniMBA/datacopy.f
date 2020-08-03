C ======================================================================
      Subroutine copyEdge2Triangle(nP, nE, IPE, Rdata, Edata, iW)
C ======================================================================
C Routines copies data from edges to triangles 
C The working memory is nP + 3*nE + nR < nP + 6*nE
C ======================================================================
      implicit none

      Integer  nP, nE, IPE(3, *), iW(*)
      Real*8  Rdata(*), Edata(3, *)

C LOCAL VARIABLES
      Integer  i, n, inEP, iIEP, iIRE, iR, nR

C ======================================================================
      inEP = 1
      iIEP = inEP + nP
      iIRE = iIEP + 3 * nE

      Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

      Do n = 1, nE
         Do i = 1, 3
            iR = iW(iIRE + 3 * (n-1) + i-1)
            Edata(i, n) = Rdata(iR)
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine packMesh(nP,XYP,lbP, nF,IPF,lbF, nE,IPE,lbE, iW, rW)
C ======================================================================
C Routine packs the mesh data into two arrays iW(11:Ni) and rW(1:Nr).
C The numbers Ni,Nr are saved in iW(1:2), respectively.
C The numbers nP,nF,nE are saved in iW(3:5), respectively.
C Positions iW(6:10) are reserved for the user.
C ======================================================================
      implicit none 

      Integer  nP, lbP(*), nF, IPF(2, *), lbF(*), nE, IPE(3, *), lbE(*)
      Real*8   XYP(2, *)

      Integer  iW(*), i, k, n
      Real*8   rW(*)

c pack mesh numbers
      iW(1) = 10 + nP + 3 * nF + 4 * nE
      iW(2) = 2 * nP
      iW(3) = nP
      iW(4) = nF
      iW(5) = nE

c pack coordinates of points
      k = 0
      Do n = 1, nP
         Do i = 1, 2
            k = k + 1
            rW(k) = XYP(i, n)
         End do
      End do

c pack labels of points        
      k = 10
      Do n = 1, nP
         k = k + 1
         iW(k) = lbP(n)
      End do

c pack boundary edges
      Do n = 1, nF
         Do i = 1, 2
            k = k + 1
            iW(k) = IPF(i, n)
         End do
      End do

c pack boundary labels
      Do n = 1, nF
         k = k + 1
         iW(k) = lbF(n)
      End do

c pack connectivity table for triangles
      Do n = 1, nE
         Do i = 1, 3
            k = k + 1
            iW(k) = IPE(i, n)
         End do
      End do

c pack labels of triangles 
      Do n = 1, nE
         k = k + 1
         iW(k) = lbE(n)
      End do

      Return
      End



C ======================================================================
      Subroutine unpackMesh(nP,ipXYP,iplbP, nF,ipIPF,iplbF, 
     &                                      nE,ipIPE,iplbE, iW, rW)
C ======================================================================
C Routine returns pointers to mesh packed in arrays iW(3:,Ni), rW(1:Nr).
C ======================================================================
      implicit none 

      Integer  nP,ipXYP,iplbP, nF,ipIPF,iplbF, nE,ipIPE,iplbE, iW(*)
      Real*8   rW(*)

c unpack mesh numbers
      nP = iW(3)
      nF = iW(4)
      nE = iW(5)

c extract mesh pointers
      ipXYP = 1
      iplbP = 11
      ipIPF = iplbP + nP
      iplbF = ipIPF + 2*nF
      ipIPE = iplbF + nF
      iplbE = ipIPE + 3*nE

      Return
      End


