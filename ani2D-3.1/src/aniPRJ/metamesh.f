C =====================================================================
      Subroutine MetaMesh(nv, vrt, nt, tri, nv2, vrt2, nt2, tri2,
     &                    nv12, nv12max, vrt12, 
     &                    nt12, nt12max, tri12, parents,
     &                    MaxWi, MaxWr, iW, rW, iERR)
C =====================================================================
      implicit none

      Integer  nv, nt, tri(3, *), nv2, nt2, tri2(3, *)
      Real*8   vrt(2, *), vrt2(2, *)

      Integer  nv12, nv12max, nt12, nt12max, tri12(3, *), parents(2, *)
      Real*8   vrt12(2, *)

      Integer  iW(*), MaxWi, MaxWr, iERR
      Real*8   rW(*)


c LOCAL VARIABLEs
      Real*8   F(1), G(1), xy(2), xc,yc, XYP(2, 6), tri_area2, v,v1,dV
      Integer  i,l, n,n2,n2t, m, l1,l2,l2t, iv1,iv2,iv3, jv1,jv2,jv3
      Integer  iLST, iMRK, iIEE, inEP, iIEP, iLIN, iEnd, info(3)
      Integer  rVOL, rLIN, rEnd, nOverlaps
      Logical  flagOVERLAP

C =====================================================================
      iERR = 0  ! reserved for future

c ... memory management
      iLST = 1
      iMRK = iLST + nt2
      iIEE = iMRK + nt2
      inEP = iIEE + 3 * nt2
      iIEP = inEP + nv2
      iEnd = iIEP + 3 * nt2

      rVOL = 1
      rEnd = rVOL + nt
      
      iLIN = inEP
      rLIN = rEnd

      If(iEnd + 3*nv2.GT.MaxWi) Then
         iERR = 1001 
         Call errMesg(iERR, 'MetaMesh', 'MaxWi is too small')
      End if
      If(rEnd + 2*nv2.GT.MaxWr) Then
         iERR = 1002
         Call errMesg(iERR, 'MetaMesh', 'MaxWr is too small')
      End if


c ... create a map E2 -> E2
      Call listE2E(nv2, nt2, tri2, iW(iIEE), iW(inEP), iW(iIEP))


c ... populate the quad-tree structure
      Call scale2Meshes(nv, vrt, nv2, vrt2, .TRUE.)

      info(1) = 1
      info(2) = 0
      Call LINTRP(nt2,tri2, nv2,vrt2, 1,F, 0,xy, G, 
     &            iW(iLIN), rW(rLIN), info)


c ... clean the list
      Do n2 = 1, nt2
         iW(iMRK + n2-1) = 0
      End do


c ... order the triangles counterclock wise
      Do n = 1, nt
         iv1 = tri(1, n)
         iv2 = tri(2, n)
         iv3 = tri(3, n)

         v = tri_area2(vrt(1,iv1), vrt(1,iv2), vrt(1,iv3))

         rW(rVOL + n-1) = dabs(v)

         If(v.LT.0) Then
            tri(1, n) = iv2 
            tri(2, n) = iv1 
         End if
      End do

      Do n = 1, nt2
         iv1 = tri2(1, n)
         iv2 = tri2(2, n)
         iv3 = tri2(3, n)

         v = tri_area2(vrt2(1,iv1), vrt2(1,iv2), vrt2(1,iv3))

         If(v.LT.0) Then
            tri2(1, n) = iv2 
            tri2(2, n) = iv1 
         End if
      End do


c ... loop over elements of the first mesh and create a meta-mesh
      info(1) = 0
      info(2) = 0

      nv12 = 0
      nt12 = 0
      nOverlaps = 0
      dV = 0

      Do n = 1, nt
         v1 = rW(rVOL + n-1)

c find triangle #n2 in mesh 2 that intersects triangle #n in mesh 1
         iv1 = tri(1, n)
         iv2 = tri(2, n)
         iv3 = tri(3, n)

         Call tri_center(vrt(1,iv1), vrt(1,iv2), vrt(1,iv3), xc, yc)

         xy(1) = xc
         xy(2) = yc
         Call LINTRP(nt2,tri2, nv2,vrt2, 1,F, 1,xy, G, 
     &               iW(iLIN), rW(rLIN), info)

         n2 = info(3)
         l1 = 1
         l2 = 1
         iW(iLST) = n2 
         iW(iMRK + n2-1) = 1

 100     Continue
c loop over neighboors of the previous triangles
         flagOVERLAP = .FALSE.
         Do l = l1, l2
            n2 = iW(iLST + l-1)

            jv1 = tri2(1, n2)
            jv2 = tri2(2, n2)
            jv3 = tri2(3, n2)
            Call tri2tri(vrt( 1,iv1), vrt( 1,iv2), vrt( 1,iv3), 
     &                   vrt2(1,jv1), vrt2(1,jv2), vrt2(1,jv3), m, XYP)
            nOverlaps = nOverlaps + 1

            If(nt12 + m.GT.nt12max) Then  
               Call errMesFEM(1001, 'metamesh', 
     &                        'Not enough memory for meta-mesh')
            End if

            Do i = 3, m
               nt12 = nt12 + 1
               tri12(1, nt12) = nv12 + 1
               tri12(2, nt12) = nv12 + i-1
               tri12(3, nt12) = nv12 + i

               parents(1, nt12) = n
               parents(2, nt12) = n2

               v  = tri_area2(XYP(1,1), XYP(1,i-1), XYP(1,i))
               v1 = v1 - v
            End do

            If(m.GE.3) Then
               flagOVERLAP = .TRUE. 
               If(nv12 + m.GT.nv12max) Then  
                  Call errMesFEM(1001, 'metamesh', 
     &                           'Not enough memory for meta-mesh')
               End if

               Do i = 1, m 
                  nv12 = nv12 + 1
                  vrt12(1, nv12) = XYP(1, i) 
                  vrt12(2, nv12) = XYP(2, i) 
               End do 
            End if
         End do

c add more neighboors if the intersection is not complete
c otherwise, clean the marked triangles
         If(flagOVERLAP) Then
            l2t = l2
            Do l = l1, l2t
               n2 = iW(iLST + l-1)

               Do i = 1, 3
                  n2t = iW(iIEE + 3*(n2-1) + i-1)
                  If(n2t.GT.0 .AND. iW(iMRK + n2t-1).EQ.0) Then
                     l2 = l2 + 1
                     iW(iLST + l2-1)  = n2t
                     iW(iMRK + n2t-1) = 1
                  End if
               End do
            End do

            l1 = l2t + 1
            Goto 100
         Else
            dV = dV + v1
            Do l = 1, l2
               n2 = iW(iLST + l-1)
               iW(iMRK + n2-1) = 0
            End do
         End if
      End do

c scale the meshes back (the flag should be false)
      Call scale2Meshes(nv, vrt, nv2, vrt2, .FALSE.)
      Call scale2Meshes(nv12, vrt12, 0, vrt2, .FALSE.)

      m = 100 * nOverlaps / (nt * nt2)

      Write(*,'(A,2I7)') 'PRJ: 1st mesh: nv/nt:', nv,  nt
      Write(*,'(A,2I7)') '     2nd mesh: nv/nt:', nv2, nt2
      Write(*,'(A,2I7,A,I8,A,I2,A)') '     metamesh: nv/nt:', nv12,nt12,
     &         '   after', nOverlaps, ' overlaps (', m, '%)'
      Write(*,'(A,E12.4)') '     area mismatch is', dV

      Return
      End




