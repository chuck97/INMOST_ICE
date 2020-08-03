C ======================================================================
      Subroutine AssembleFromTemplate(
C ======================================================================
     &           nP, nF, nE, XYP, lbP, IPF, lbF, IPE, lbE, 
     &           FEM2Dext, DDATA, IDATA, control, 
     &           MaxF, MaxA, IA, JA, A, F, nRow, nCol,
     &           MaxWi, MaxWr, iW, rW)
C ======================================================================
      implicit none

      Include 'fem2Dtri.fd'
      Include 'assemble.fd'
C ======================================================================
C  The routine assembles a stiffness matrix for a bilinear form B.
C  Here we describe a new way to assemble elemental matrices. 
C
C  Parameter are described in template.f
C ======================================================================
C
C *** Authors: K. Lipnikov (lipnikov@gmail.com)
C              Yu. Vassilevski (yuri.vasilevski@gmail.com)
C *** Date:   2005 - 2011
C *** Complains & Comments: lipnikov@gmail.com
C *** External routines: FEM2Dext
C
C ======================================================================
      Real*8   XYP(2, *)
      Integer  lbP(*), IPF(2, *), lbF(*), IPE(3, *), lbE(*)
      Integer  nP, nF, nE

      Integer  IDATA(*), iSYS(MAXiSYS), control(3)
      Real*8   DDATA(*)
      EXTERNAL FEM2Dext

      Integer  IA(*), JA(*), MaxF, MaxA, nRow, nCol
      Real*8   A(*), F(*)

      Integer  MaxWi, MaxWr, iW(*)
      Real*8   rW(*)

C ======================================================================
C Local variables
      Integer templateR(MaxSizeDG), templateC(MaxSizeDG)
      Real*8  Aloc(MaxSizeDG, MaxSizeDG), Floc(MaxSizeDG)
      Integer indR(MaxSizeDG), indC(MaxSizeDG)
 
      Integer assembleStatus, iPrint

      Integer lbFloc(3,3), lbPloc(3,2)
      Logical flagMS, flagMG, flagFA, flagFR, flagFC

      Integer nEdofR, nEdofC, nR
      Integer irp,irr,ire,ire_dg, icp,icr,ice,ice_dg, irr_half,icr_half

      Integer imER, imEC, inEP, iIEP, iIFE, iIRE, inER, iIER, iIEE
      Integer iEnd, jEnd, iIA, iJA
      Integer ip1,ip2,ip3, ir1,ir2,ir3, ie1,ie2,ie3
      Integer ir,ic, iP,iQ, iPt,iFt,iRt,iEt
      Integer i,j,k,m,n,p, i1,i2, vecC, vecR, non_zero, iusage, rusage

      Integer map1(6),map2(6),map3(6), smap(MaxSizeDG),tmap(MaxSizeDG)
      Real*8  XY1(2,2), XY2(2,2), XY3(2,2), s

      Integer findY
      Character*80 message
      Logical ifXnode

C ======================================================================
c ... check the input parameters
      If(nP.LE.0   .OR. nF.LT.0 .OR. nE.LE.0 .OR.
     &   MaxF.LE.0 .OR. MaxA.LE.0 . OR. MaxWi.LE.0 .OR. MaxWr.LE.0) 
     &   Call errMesFEM(4101, 'AssembleFromTemplate', 
     &                        'mesh numbers are negative')

c ... unpack control parameters
      assembleStatus = control(1)
      iPrint = control(2)

c ... default values for flags
      flagMG = .TRUE.

      flagMS = IAND(assembleStatus, MATRIX_SYMMETRIC).NE.0
      flagMG = IAND(assembleStatus, MATRIX_GENERAL).NE.0
      If(flagMG) Then
         flagMS = .FALSE.
      End if

      flagFA  = IAND(assembleStatus, FORMAT_AMG).NE.0
      flagFR  = IAND(assembleStatus, FORMAT_CSR).NE.0
      flagFC  = IAND(assembleStatus, FORMAT_CSC).NE.0
      If(flagFC) Then
         flagFA = .FALSE.
         flagFR = .FALSE.
         If(iPrint.GT.0) 
     &      Write(*,'(A)') 'FEM: sparse compressed column format'
      Else If(flagFR) Then
         flagFA = .FALSE.
         If(iPrint.GT.0) 
     &      Write(*,'(A)') 'FEM: sparse compressed row format'
      Else 
         If(iPrint.GT.0) 
     &      Write(*,'(A)') 'FEM: sparse compressed AMG row format'
      End if
         

c ... memory allocation
c     pointers to arrays collecting number of dofs of type Edof 
      imER = 1
      imEC = imER + nE + 1
      inEP = imEC + nE + 1
      iIEP = inEP + nP
      iIEE = iIEP + max(3 * nE, 2 * nF)
      iIFE = iIEE + 3 * nE 
      iEnd = iIFE + 3 * nE
      If(iEnd.GT.MaxWi) Goto 500

      Call listE2E(nP, nE, IPE, iW(iIEE), iW(inEP), iW(iIEP))

      Call listE2F(nP, nF, nE, IPF, IPE, iW(iIFE), iW(inEP), iW(iIEP))


c ... get number of dofs that are strictly inside the element
      iW(imER) = 1
      iW(imEC) = 1

      nEdofR = 0
      nEdofC = 0

      Do i = 1, 3
         lbFloc(i,1) = 0
         lbPloc(i,1) = 0
      End do

      Do n = 1, nE
         ip1 = IPE(1, n)
         ip2 = IPE(2, n)
         ip3 = IPE(3, n)

         ie1 = iW(iIEE + 3*(n-1))
         ie2 = iW(iIEE + 3*(n-1) + 1)
         ie3 = iW(iIEE + 3*(n-1) + 2)

c        not all information is available at this moment to populate iSYS
         nR = 0
         Call encodeISYS(n, ip1,ip2,ip3, 0,0,0, 
     &                      ie1,ie2,ie3, nP,nR,nE, iSYS)

c        signaling a special request
         iSYS(1) = -1
         Call FEM2Dext(XYP(1, ip1), XYP(1, ip2), XYP(1, ip3),
     &                 lbE(n), lbFloc, lbPloc, DDATA, IDATA, iSYS,
     &                 MaxSizeDG, Aloc, Floc, ir, ic, 
     &                 templateR, templateC)
 
         If(ir.GT.MaxSizeDG .OR. ic.GT.MaxSizeDG) 
     &      Call errMesFEM(2001, 'AssembleFromTemplate', 
     &                     'MaxSizeDG is small')

c        extracting vector information
         Call vectorDOF(ir, templateR, vecR)
         Call vectorDOF(ic, templateC, vecC)

c        analyze rows
         irp = 0
         irr = 0
         ire = 0
         ire_dg = 0
         Do i = 1, ir
            If(templateR(i).EQ.Vdof) Then
               irp = irp + 1
            Else If(templateR(i).EQ.Rdof .OR. 
     &              templateR(i).EQ.RdofOrient) Then
               irr = irr + 1
            Else If(templateR(i).EQ.Edof .OR.
     &              ifXnode(templateR(i), DGdof)) Then
               ire = ire + 1
            Else
               ire_dg = ire_dg + 1
            End if
         End do
         irp = irp / 3 
         irr = irr / 3 
         irr_half = irr / (vecR + 1) 

         nEdofR = nEdofR + ire
         iW(imER + n) = nEdofR + 1

c        analyze columns
         icp = 0
         icr = 0
         ice = 0
         ice_dg = 0
         Do i = 1, ic
            If(templateC(i).EQ.Vdof) Then
               icp = icp + 1
            Else If(templateC(i).EQ.Rdof .OR.
     &              templateC(i).EQ.RdofOrient) Then
               icr = icr + 1
            Else If(templateC(i).EQ.Edof .OR.
     &              ifXnode(templateC(i), DGdof)) Then
               ice = ice + 1
            Else
               ice_dg = ice_dg + 1
            End if
         End do
         icp = icp / 3
         icr = icr / 3
         icr_half = icr / (vecC + 1) 

         nEdofC = nEdofC + ice
         iW(imEC + n) = nEdofC + 1
      End do


c ... optional memory allocation for map E->R and its reverse R->E
      If(irr.GT.0 .OR. icr.GT.0) Then      
         iIRE = iEnd
         iEnd = iIRE + 3 * nE
         If(iEnd.GT.MaxWi) Goto 500

         Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

         inER = iEnd
         iIER = inER + nR
         iEnd = iIER + 3 * nE 
         If(iEnd.GT.MaxWi) Goto 500

         Call backReferences(nR, nE, 3,3, iW(iIRE), iW(inER), iW(iIER))
      Else
         iIRE = 1
         iIER = 1
      End if

 500  Continue
      If(iEnd.GT.MaxWi) Then
         Write(message,'(A,I10)')
     &      'The approximate size of iW is ', iEnd
         Call errMesFEM(1001, 'AssembleFromTemplate', message)
      End if

      nRow = irp * nP + irr * nR + nEdofR
      nCol = icp * nP + icr * nR + nEdofC

      If(nRow.GT.MaxF) Call errMesFEM(1013, 'AssembleFormTemplate',
     &                        'local parameter MaxF is small')

      iusage = (iEnd * 100) / MaxWi
      rusage = 0


c ... logical assembling of the matrix
c     rough estimate of the number of non-zero elements
      Do i = 1, nRow + 1
         IA(i) = 1
      End do

      Do n = 1, nE
         ip1 = IPE(1, n)
         ip2 = IPE(2, n)
         ip3 = IPE(3, n)

         If(irr.GT.0 .OR. icr.GT.0) Then
            ir1 = iW(iIRE + 3*(n-1))
            ir2 = iW(iIRE + 3*(n-1) + 1)
            ir3 = iW(iIRE + 3*(n-1) + 2)
         Else
            ir1 = 0
            ir2 = 0
            ir3 = 0
         End if

         ie1 = iW(iIEE + 3*(n-1))
         ie2 = iW(iIEE + 3*(n-1) + 1)
         ie3 = iW(iIEE + 3*(n-1) + 2)

         Call encodeISYS(n, ip1,ip2,ip3, ir1,ir2,ir3, 
     &                      ie1,ie2,ie3, nP,nR,nE, iSYS)

         iSYS(1) = -1
         Call FEM2Dext(XYP(1, ip1), XYP(1, ip2), XYP(1, ip3),
     &                 lbE(n), lbFloc, lbPloc, DDATA, IDATA, iSYS,
     &                 MaxSizeDG, Aloc, Floc, ir, ic, 
     &                 templateR, templateC)

c        extracting vector information
         Call vectorDOF(ir, templateR, vecR)
  
         Call globalDOF(n, nP, nR, IPE, iW(iIRE), iW(iIEE), templateR,
     &                  irp, irr, irr_half, iW(imER), ire_dg, indR)
 
         ic = 3 * icp + 3 * icr + iW(imEC + n) - iW(imEC + n-1) 
c        add correction for DG methods
         If(ice_dg.GT.0) Then
            Do k = 1, 3
               iEt = iW(iIEE + 3*(n-1) + k-1)
               If(iEt.GT.0) ic = ic + iW(imEC + iEt) - iW(imEC + iEt-1)
            End do
         End if

         Do i = 1, ir
            iP = indR(i) + 1
            IA(iP) = IA(iP) + ic - 1
         End do
      End do

      Do i = 1, nRow
         IA(i + 1) = IA(i + 1) + IA(i)
      End do
      IA(1) = 1

      iEnd = IA(nRow + 1)
      If(iEnd.GT.MaxA) Call errMesFEM(1014,
     &     'AssembleFromTemplate', 'local parameter MaxA is small')


C ... fill in the sparsity structure
      Do n = 1, nRow
         F(n) = 0D0
      End do

      Do n = 1, IA(nRow + 1) - 1
         JA(n) = 0
          A(n) = 0D0
      End do

      Do n = 1, nE
c        provide the user global numbering of local objects
         ip1 = IPE(1, n)
         ip2 = IPE(2, n)
         ip3 = IPE(3, n)

         If(irr.GT.0 .OR. icr.GT.0) Then
            ir1 = iW(iIRE + 3*(n-1))
            ir2 = iW(iIRE + 3*(n-1) + 1)
            ir3 = iW(iIRE + 3*(n-1) + 2)
         Else
            ir1 = 0
            ir2 = 0
            ir3 = 0
         End if

         ie1 = iW(iIEE + 3*(n-1))
         ie2 = iW(iIEE + 3*(n-1) + 1)
         ie3 = iW(iIEE + 3*(n-1) + 2)
 
         Call encodeISYS(n, ip1,ip2,ip3, ir1,ir2,ir3, 
     &                      ie1,ie2,ie3, nP,nR,nE, iSYS)

c        save point and edge labels
         Do i = 1, 3
            iPt = IPE(i, n)
            lbPloc(i,1) = lbP(iPt)
         End do

         Do i = 1, 3
            iFt = iW(iIFE + 3 * (n-1) + i-1)
            lbFloc(i,1) = 0
            If(iFt.GT.0) lbFloc(i,1) = lbF(iFt)
         End do

c        save vertices
         Do i = 1, 2
            XY1(i, 1) = XYP(i, ip1)
            XY2(i, 1) = XYP(i, ip2)
            XY3(i, 1) = XYP(i, ip3)
         End do

c        DG extension: opposite points and faces in neighboring triangles
         Call oppositeV(1, XYP, IPE, iW(iIFE), lbP, lbF, n,ie1, 
     &                  XY3, lbPloc, lbFloc, map1)
         Call oppositeV(2, XYP, IPE, iW(iIFE), lbP, lbF, n,ie2, 
     &                  XY1, lbPloc, lbFloc, map2)
         Call oppositeV(3, XYP, IPE, iW(iIFE), lbP, lbF, n,ie3, 
     &                  XY2, lbPloc, lbFloc, map3)

         Call FEM2Dext(XY1, XY2, XY3,
     &                 lbE(n), lbFloc, lbPloc, DDATA, IDATA, iSYS,
     &                 MaxSizeDG, Aloc, Floc, ir, ic, 
     &                 templateR, templateC)

         Call vectorDOF(ir, templateR, vecR)
         Call vectorDOF(ic, templateC, vecC)

c        correct orientation of the basis functions in the local matrix
c        R-unknowns: the global orientation is from bigger E to smaller E
         i = 1
         Do While(i.LE.ir) 
            If(templateR(i).EQ.RdofOrient) Then
               Do k = 0, 2
                  iRt = iW(iIRE + 3*(n-1) + k)
                  iEt = findY(iRt, iW(iIER), iW(inER), n)

                  If(n.LT.iEt) Then
                     Do j = 1, ic
                        Aloc(i, j) = -Aloc(i, j)
                     End do
                     Floc(i) = -Floc(i)
                  End if
                  i = i + 1
               End do
            Else
               i = i + 1
            End if
         End do

         i = 1
         Do While(i.LE.ic) 
            If(templateC(i).EQ.RdofOrient) Then
               Do k = 0, 2
                  iRt = iW(iIRE + 3*(n-1) + k)
                  iEt = findY(iRt, iW(iIER), iW(inER), n)

                  If(n.LT.iEt) Then
                     Do j = 1, ic
                        Aloc(j, i) = -Aloc(j, i)
                     End do
                  End if
                  i = i + 1
               End do
            Else
               i = i + 1
            End if
         End do

c        create arrays of row and column indices 
         Call globalDOF(n, nP, nR, IPE, iW(iIRE), iW(iIEE), templateR,
     &                  irp, irr, irr_half, iW(imER), ire_dg, indR)

         Call globalDOF(n, nP, nR, IPE, iW(iIRE), iW(iIEE), templateC,
     &                  icp, icr, icr_half, iW(imEC), ice_dg, indC)

c        create supermaps
         k = 3 * irp + 3 * irr + iW(imER + n) - iW(imER + n-1)
         Call supermap(k, ir, map1, map2, map3, templateR, smap)

         k = 3 * icp + 3 * icr + iW(imEC + n) - iW(imEC + n-1)
         Call supermap(k, ic, map1, map2, map3, templateC, tmap)

c        assemble the right hand side
         Do i = 1, ir
            iP = indR(smap(i))
            F(iP) = F(iP) + Floc(i)
         End do

c        assemble the elemental matrix
         Do i = 1, ir
            iP = indR(smap(i))

            Do 600 j = 1, ic
               iQ = indC(tmap(j))

               Do k = IA(iP), IA(iP + 1) - 1
                  If(JA(k).EQ.0) Then
                     JA(k) = iQ
                      A(k) = Aloc(i, j)
                     Goto 600
                  Else If(JA(k).EQ.iQ) Then
                     A(k) = A(k) + Aloc(i, j)
                     Goto 600
                  End if
               End do

               Call errMesFEM(6001, 'AssembleFromTemplate', 
     &                       'System error: logical assembling failed')
 600        Continue
         End do
      End do


C ... cleaning sparcity defects
      n  = 1
      i2 = 1 
      Do iP = 1, nRow
         i1 = i2
         i2 = IA(iP + 1)

         Do k = i1, i2 - 1
            If(JA(k).GT.0) Then
               If(A(k).NE.0D0) Then
                  JA(n) = JA(k)
                   A(n) =  A(k)
                  n = n + 1
               End if
            Else
              Goto 700
            End if
         End do

 700     Continue
         IA(iP + 1) = n
      End do


C ... comply with the AMG format
      If(flagFA) Then
         Do 800 n = 1, nRow
            Do j = IA(n), IA(n + 1) - 1
               If(JA(j).EQ.n) Then
                  Call swapii(JA(j), JA(IA(n)))
                  Call swapdd( A(j),  A(IA(n)))
                  Goto 800
               End if
            End do
 800     Continue
      End if


c ... convert row format to column format using ugly method
      non_zero = IA(nRow + 1) - 1 

      If(flagFC) Then
         iIA  = 1
         iJA  = iIA + nRow + 1
         iEnd = iJA + non_zero

         jEnd = non_zero

         If(iEnd.GT.MaxWi) Goto 500 
         If(jEnd.GT.MaxWr) Call errMesFEM(1001, 
     &     'AssembleFormTemplate', 'Not enough Real*8 memory')

         Do i = 1, nRow + 1
            iW(i) = IA(i)
         End do

         Do i = 1, non_zero
            iW(iJA + i-1) = JA(i)
            rW(i)         =  A(i)
         End do

         Call AMG2CSC(nRow, iW(iIA), iW(iJA), rW, nCol, IA, JA, A)

         ir   = nRow   
         nRow = nCol
         nCol = ir

         iusage = max(iusage, iEnd * 100 / MaxWi)
         rusage = jEnd * 100 / MaxWr
      End if
         

c ... return some useful information
      iW(1) = nP
      iW(2) = nR
      iW(3) = nE

      If(iPrint.GT.0) Write(*,'(5X,A,I6,A,I6,A,I9,A)') 
     &        'matrix size:', nRow,' x', nCol, 
     &        '  (', non_zero, ' non-zero entries)'

      If(iPrint.GE.9) 
     &   Write(*,'(A,3(I2,A))') '     work memory usage: ', 
     &        iusage, '% (Integer) and ', rusage, '% (Real*8)'

      Return
      End



C ======================================================================
      Subroutine globalDOF(n, nP, nR, IPE, IRE, IEE, templateC, 
     &                     icp, icr, icr_half, nCE, ice_dg, indC)
C ======================================================================
C  Technical subroutine for calculating global indexing 
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer  n, nP, nR, IPE(3, *), IRE(3, *), IEE(3, *), templateC(*)
      Integer  icp, icr, icr_half, nCE(*), ice_dg, indC(*)

      Integer i,j,k, ip1,ip2,ip3, ir1,ir2,ir3, jcp,jcr,jce, ie1,ie2,ie3
      Integer iscR, iscE, jce12, jce23, jce31, icr_copy, icr_vecx, ic
      Integer orient(3), ixc(3)

      Logical ifXnode

C ======================================================================
      iscR = icp * nP
      iscE = iscR + icr * nR

      ip1 = IPE(1, n)
      ip2 = IPE(2, n)
      ip3 = IPE(3, n)

      If(icr.GT.0) Then
         ir1 = IRE(1, n)
         ir2 = IRE(2, n)
         ir3 = IRE(3, n)
      End if

      ie1 = IEE(1, n)
      ie2 = IEE(2, n)
      ie3 = IEE(3, n)

      ic = 3 * icp + 3 * icr + (nCE(n + 1) - nCE(n))
      If(ice_dg.GT.0) Then
         If(ie1.GT.0) ic = ic + (nCE(ie1 + 1) - nCE(ie1))
         If(ie2.GT.0) ic = ic + (nCE(ie2 + 1) - nCE(ie2))
         If(ie3.GT.0) ic = ic + (nCE(ie3 + 1) - nCE(ie3))
      End if

c     orientation of edges from smaller to higher
      Do i = 1, 3
         orient(i) = 1
      End do
      If(ip1.GT.ip2) orient(1) = 0
      If(ip2.GT.ip3) orient(2) = 0
      If(ip3.GT.ip1) orient(3) = 0

c     populate all point-based dofs 
      icr_copy = icr_half - 1
      icr_vecx = 0

      jcp = 0
      jcr = 0
      jce = 0
      jce12 = 0
      jce23 = 0
      jce31 = 0

      i = 1
      Do While(i.LE.ic) 
c        populate vertex-based dofs
         If(templateC(i).EQ.Vdof) Then
            indC(i)     = jcp * nP + ip1 
            indC(i + 1) = jcp * nP + ip2 
            indC(i + 2) = jcp * nP + ip3 

            i = i + 3
            jcp = jcp + 1

c        populate edge-based dofs for vector component X
         Else If(templateC(i).EQ.Rdof .OR.
     &           templateC(i).EQ.RdofOrient) Then
            Do k = 1, 3
               ixc(k) = icr_vecx 
     &                + icr_copy * (1 - orient(k)) + jcr * orient(k)
            End do

            indC(i)     = iscR + ixc(1) * nR + ir1
            indC(i + 1) = iscR + ixc(2) * nR + ir2
            indC(i + 2) = iscR + ixc(3) * nR + ir3

            i = i + 3

            jcr = jcr + 1
            icr_copy = icr_copy - 1
            If(icr_copy.LT.0) Then
               icr_vecx = icr_half
               icr_copy = icr_half - 1
               jcr = 0
            End if

c        populate element-based dofs
         Else If(templateC(i).EQ.Edof .OR.
     &           ifXnode(templateC(i), DGdof)) Then
            indC(i) = iscE + nCE(n) + jce

            i = i + 1
            jce = jce + 1

c        populate exterior element-based dofs
         Else If(ifXnode(templateC(i), DGdof12)) Then
            indC(i) = iscE + nCE(ie1) + jce12

            i = i + 1
            jce12 = jce12 + 1
         Else If(ifXnode(templateC(i), DGdof23)) Then
            indC(i) = iscE + nCE(ie2) + jce23

            i = i + 1
            jce23 = jce23 + 1
         Else If(ifXnode(templateC(i), DGdof31)) Then
            indC(i) = iscE + nCE(ie3) + jce31

            i = i + 1
            jce31 = jce31 + 1
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine oppositeV(i, XYP, IPE, IFE, lbP, lbF, n1, n2, 
     &                     XYi, lbPi, lbFi, map) 
C ======================================================================
C  Routine calculates opposite vertices of two triangles n1 and n2,
C  saves them in in array XYi and creates a map for triangle n2.
C ======================================================================
      implicit none
      Integer i, IPE(3, *), IFE(3, *), lbP(*), lbF(*), n1, n2
      Integer lbPi(3,2), lbFi(3,3), map(6)
      Real*8  XYP(2, *), XYi(2, 2)

      Integer iP1,iP2,iP3, iQ1,iQ2, iPt,iFt, i1,i2,i3, k
      Logical check22

      Integer iref(4)
      DATA    iref/1,2,3,1/
C ======================================================================
      i2 = iref(i  + 1)
      i3 = iref(i2 + 1)

      iP1 = IPE(i,  n1)
      iP2 = IPE(i2, n1)

      If(n2.EQ.0) Then
         map(1) = 0
         lbPi(i3,2) = 0
         return
      End if

c     find opposite vertex
      Do k = 1, 3
         iPt = IPE(k, n2)

         If(iPt.EQ.iP1) Then
            map(1) = k
         Else If(iPt.EQ.iP2) Then
            map(2) = k
         Else
            map(3) = k 
            XYi(1, 2) = XYP(1, iPt)
            XYi(2, 2) = XYP(2, iPt)

            iP3 = iPt
            lbPi(i3, 2) = lbP(iPt)
         End if
      End do

c     find faces
      lbFi(i, 2) = 0
      lbFi(i, 3) = 0

      Do i1 = 1, 3
         i2 = iref(i1 + 1)

         iQ1 = IPE(i1, n2)
         iQ2 = IPE(i2, n2)

         iFt = IFE(i1, n2)

         If(check22(iQ1,iQ2, iP2,iP3)) Then
            map(5) = i1
            If(iFt.GT.0) lbFi(i, 2) = lbF(iFt)  
         Else If(check22(iQ1,iQ2, iP1,iP3)) Then
            map(6) = i1
            If(iFt.GT.0) lbFi(i, 3) = lbF(iFt) 
         Else
            map(4) = i1
         End if
      End do

      Return 
      End



C ======================================================================
      Subroutine supermap(k, ir, map1, map2, map3, template, smap)
C ======================================================================
C  Routine reconstructs global map from elemental maps of triangles
C ======================================================================
      implicit none
      include "fem2Dtri.fd"

      Integer k, ir, map1(6), map2(6), map3(6), template(*), smap(*)

      Integer i,n,m, n12,n23,n31, i12,i23,i31, m12,m23,m31
      Logical ifXnode
C ======================================================================
      Do i = 1, k
         smap(i) = i
      End do
      If(k.EQ.ir) Return

      i12 = 0
      i23 = 0
      i31 = 0

c ... populate supermap
      n = k
      m = k + 1
      Do While(m.LE.ir) 
         If(ifXnode(template(m), DGdof12)) Then
            If(ifXnode(template(m), Vdof)) Then
               Do i = 1, 3
                  smap(n + i) = n + map1(i)
               End do
               n = n + 3
               m = m + 3
            Else If(ifXnode(template(m), Rdof)) Then
               n12 = n
               i12 = i12 + 1
               Do i = 1, 3
                  smap(n + i) = n + map1(3+i)
               End do
               n = n + 3
               m = m + 3
            Else If(ifXnode(template(m), Edof)) Then
               n = n + 1
               m = m + 1
               smap(n) = n
            End if
         Else
            m = m + 1
         End if
      End do 

      m = k + 1
      Do While(m.LE.ir) 
         If(ifXnode(template(m), DGdof23)) Then
            If(ifXnode(template(m), Vdof)) Then
               Do i = 1, 3
                  smap(n + i) = n + map2(i)
               End do
               n = n + 3
               m = m + 3
            Else If(ifXnode(template(m), Rdof)) Then
               n23 = n
               i23 = i23 + 1
               Do i = 1, 3
                  smap(n + i) = n + map2(3+i)
               End do
               n = n + 3
               m = m + 3
            Else If(ifXnode(template(m), Edof)) Then
               n = n + 1
               m = m + 1
               smap(n) = n
            End if
         Else
            m = m + 1
         End if
      End do 

      m = k + 1
      Do While(m.LE.ir) 
         If(ifXnode(template(m), DGdof31)) Then
            If(ifXnode(template(m), Vdof)) Then
               Do i = 1, 3
                  smap(n + i) = n + map3(i)
               End do
               n = n + 3
               m = m + 3
            Else If(ifXnode(template(m), Rdof)) Then
               n31 = n
               i31 = i31 + 1
               Do i = 1, 3
                  smap(n + i) = n + map3(3+i)
               End do
               n = n + 3
               m = m + 3
            Else If(ifXnode(template(m), Edof)) Then
               n = n + 1
               m = m + 1
               smap(n) = n
            End if
         Else
            m = m + 1
         End if
      End do 

c ... reorder DG edge-based unknowns (DG + Rdof)
c     we assume that edge unknowns are contiguous
      If(map1(1).NE.map1(4)) Then
         i12 = i12 - 1
         m12 = n12 - 3 * i12

         Do m = 1, i12
            Do i = 1, 3
               Call swapii(smap(m12 + i), smap(n12 + i))
            End do
            m12 = m12 + 3
            n12 = n12 - 3
         End do
      End if

      If(map2(1).NE.map2(4)) Then
         i23 = i23 - 1
         m23 = n23 - 3 * i23

         Do m = 1, i23
            Do i = 1, 3
               Call swapii(smap(m23 + i), smap(n23 + i))
            End do
            m23 = m23 + 3
            n23 = n23 - 3
         End do
      End if

      If(map3(1).NE.map3(4)) Then
         i31 = i31 - 1
         m31 = n31 - 3 * i31

         Do m = 1, i31
            Do i = 1, 3
               Call swapii(smap(m31 + i), smap(n31 + i))
            End do
            m31 = m31 + 3
            n31 = n31 - 3
         End do
      End if

      Return
      End

