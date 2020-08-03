C ======================================================================
      Subroutine orientBoundary(nP, nF, nE, XYP, IPF, IPE, iW, MaxWi)
C ======================================================================
C Routine orients the external boundary of a given mesh in such 
C a way that the domain is on the left of an edge. In orwer words, 
C IPF(1, *) and IPF(2, *) are flipped if neccessary.
C
C *** Remarks:
C        1. The size of working memory is 3 * nE + 2 * nF + nP
C ======================================================================
      Real*8  XYP(2, *)
      Integer IPE(3, *), IPF(2, *), iW(*)

C (Local variables)
      Real*8  v, tri_area2
      Integer iref(4)
      DATA    iref /1,2,3,1/

C ======================================================================
      inEP = 1
      iIEP = inEP + nP
      iIFE = iIEP + 2 * nF
      iEnd = iIFE + 3 * nE

      If(iEnd.GE.MaxWi) Call errMesMBA(1001, 'orientBoundary',
     &                             'not enough working memory')

c ... compute map E -> F
      Call listE2F(nP, nF, nE, IPF, IPE, iW(iIFE), iW(inEP), iW(iIEP))

      Do n = 1, nE
         Do i1 = 1, 3
            ife = iW(iIFE + 3 * (n - 1) + i1 - 1)
            If(ife.GT.0) Then 
               i2 = iref(i1 + 1)
               i3 = iref(i2 + 1)
 
               iP1 = IPF(1, ife)  
               iP2 = IPF(2, ife)  
               iP3 = IPE(i3, n)  

               v = tri_area2(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
               If(v.GT.0D0) Call swapii(IPF(1, ife), IPF(2, ife))
            End if
         End do
      End do      

      Return
      End

