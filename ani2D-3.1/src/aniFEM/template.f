C ======================================================================
      Subroutine BilinearFormTemplate(
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
C  Here we describe new way to assemble elemental matrices. 
C ======================================================================
C *** PROBLEM INFORMATION ***
C  The routines assemble elemental matrices for the bilinear form
C
C              B(u, v)
C
C  specified by the user. Each elemental matrix is a combination of 
C  fem2Dtri calls because FEM2Dext may consists of a few simple bilinear 
C  forms. The elemental matrix is described by array 'template' and
C  its contruction may use arrays DDATA and IDATA if neccessary. 
C
C
C *** SUBROUTINE FEM2Dext ***
C  
C     SUBROUTINE FEM2Dext(XY1, XY2, XY3, 
C    &                    lbEloc, lbFloc, lbPloc, DDATA, IDATA, iSYS,
C    &                    LDA, A, F, nRow, nCol,
C    &                    templateR, templateC)
C
C     Parameters:
C              XY1(2) - Real*8 Cartesian coordinates of the first
C                         vertex of the triangle.
C              XY2(2) - second vertex of the triangle
C              XY3(2) - third vertex of the triangle
C
C              lbEloc    - label of the triangle (material label)
C              lbFloc(3) - labels of the triangle edges (boundary label),
C                          lbFloc(i) = 0 for internal edges
C              lbPloc(3) - labels of points (copied from global lbP), 
C
C              DDATA  - Real*8  user given data (a number or an array)
C              IDATA  - Integer user given data (a number or an array)
C
C              iSYS(12) - system buffer to exchange information about mesh
C                       iSYS(3) : triangle number
C                       iSYS(4) : 1st vertex number
C                       iSYS(5) : 2nd vertex number
C                       iSYS(6) : 3rd vertex number
C
C                       iSYS(7) : 1st edge number, connects vertices 1 & 2
C                       iSYS(8) : 2nd edge number, connects vertices 2 & 3
C                       iSYS(9) : 3rd edge number, connects vertices 3 & 1
C
C                       iSYS(10) : total number of points,    nP
C                       iSYS(11) : total number of edges,     nR
C                       iSYS(12) : total number of triangles, nE
C
C                       iSYS(13) : neiboring triangle to edge 12
C                       iSYS(14) : neiboring triangle to edge 23
C                       iSYS(15) : neiboring triangle to edge 31
C
C               iSYS(5:7) and iSYS(9) may be empty if the edge degrees
C               of freedom are not used. Two routines in isys.f will help 
C               to encode and to decode iSYS().
C               iSYS(13:15) will be populated only DG method is detected.
C
C              LDA       - leading dimention of matrix A
C              A(LDA, *) - elemental matrix, the order of degrees of freedom
C                          should be in an agreement with array template
C              F(nRow)   - the vector of the elemental right-hand side     
C 
C              nRow   - the number of rows in A
C              nCol   - the number of columns in A
C
C         templateR(nRow) - array of degrees of freedom (DOF) for rows;
C                the DOF must follows in groups, e.g. three for vertices, 
C                three for edges, etc.
C                The list of availabel DOF is defined in file fem2Dtri.fd
C                         and includes Vdof (for vertices), Rdof(for edges),
C                         Edof (for elements), RdofOrient (for edges).
C                   
C                The difference between Rdof and RdofOrient is that the
C                former corresponds to scalar unknowns (e.g. Lagrange
C                multipliers) that have no orientation while the latter
C                correponds to vector unknowns (e.g. RT finite elements) 
C                that have an orientation.
C
C              templateC(nCol) - array of degrees of freedom for columns. 
C
C
C *** SPARSE MATRIX ***
C   By default, our sparse matrix is assembled in one of the row formats
C   (IA, JA, A) specified by assembleStatus. Other formats are supported
C   thorough the library of format transformations.
C
C     control(3) - integer array with control parameters:
C 
C         assembleStatus = control(1) - information about the matrix A
C         iPrint         = control(2)
C         iErrMesg       = control(3)
C
C         REMARKS on assembleStatus:
C              MATRIX_SYMMETRIC - symmetric matrix
C              MATRIX_GENERAL   - general matrix
C              FORMAT_AMG       - format used in AMG
C              FORMAT_CSR       - compressed row format (diagonal is not saved)
C
C              Any logical AND can be used to defined the variable, 
C              possible contradictions will  be checked by the code.
C
C
C     MaxF  - the maximal size of matrix A and vector F
C     MaxA  - the maximal number of nonzero entries in matrix A
C
C     IA, JA, A - sparcity structure of matrix A:
C
C                IA(k + 1) - IA(k) equals to the number of
C                            nonzero entries in the k-th row
C
C                JA(M) - list of column indexes of non-zero
C                        entries ordered by rows; M = IA(nRow + 1) - 1
C
C                A(M)  - non-zero entries ordered as in JA
C
C     nRow      - the number of rows in A
C     nCol      - the number of columns in A, nCol = nRow for
C                 symmetric bilinear forms.
C
C
C  The following default rules are applied for numbering of unknowns:
C      A) 1st, unknowns associated with vertices are numerated. 
C
C      B) 2nd, unknowns associated with edges are numerated.
C
C      C) The unknowns corresponding to vector functions (i.e. 
C         with 2 degrees of freedom per a mesh object (vertex or 
C         edge) are enumerated first by the corresponding 
C         mesh objects and then by the space coordinates, x and y.
C
C      D) Unknowns associated with elements are enumerated by 
C         elements regadless of the space coordinate.
C
C
C *** WORKING MEMORY ***
C     MaxWi  - the size of the working integer array
C     MaxWr  - the size of the working Real*8  array
C
C     iW(MaxWi)  - the Integer working array.
C     rW(MaxWr)  - the Real*8  working array.
C
C ======================================================================
C *** Note:
C       Input parameters:  nP, nE, XYP, IPE, lbE,
C                          FEM2Dext, DDATA, IDATA,
C                          MaxF, MaxA, MaxWi, MaxWr
C
C       Input / output:    iSYS, control, iW, rW
C
C       Output parameters: IA, JA, A, F, nRow, nCol
C
C ======================================================================
C
C *** Authors: K. Lipnikov (lipnikov@gmail.com)
C              Yu. Vassilevski (yuri.vasilevski@gmail.com)
C *** Date:   2005 - 2007
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
      Integer templateR(MaxSize), templateC(MaxSize)
      Integer dummyR(MaxSize), dummyC(MaxSize)

      Real*8  Aloc(MaxSize, MaxSize), Floc(MaxSize)
      Integer indR(MaxSize), indC(MaxSize)
 
      Integer gpr(MaxNumGroups), gnr(MaxNumGroups), gur(MaxNumGroups)
      Integer gpc(MaxNumGroups), gnc(MaxNumGroups), guc(MaxNumGroups)
      Integer dpr(MaxSize), dpc(MaxSize)

      Integer assembleStatus, iPrint
      Integer lbFloc(3), lbPloc(3), IEE(1)
      Logical flagMS, flagMG, flagFA, flagFR, flagFC

      Integer nR, nX, MaxL, MaxR, MaxX
      Integer iIPP,iIPR, iIFP,iIFE, iIRP,iIRR,iIRE, iIEP,jIEP,iIER
      Integer inPP,inPR, inFP,      inRP,inRR,inEP, inER,jnEP
      Integer ip1,ip2,ip3, ir1,ir2,ir3, ie1,ie2, ir,ic
      Integer irp,irr,ire, icp,icr,ice, ixr(3), ire_dg,ice_dg
      Integer iPt,iRt,iFt,iEt, iP,iQ, iiW
      Integer iIA, iJA, iEnd, jEnd, non_zero, iusage, rusage
      Integer i,i1,i2, j,k,l,m,n,p, nD, label, group, ishift

      Integer findY, orient(3)
      Character*80 message
      Logical ifXnode, flagR, flagC

C ======================================================================
c ... check the input parameters
      If(nP.LE.0   .OR. nF.LT.0 .OR. nE.LE.0 .OR.
     &   MaxF.LE.0 .OR. MaxA.LE.0 . OR. MaxWi.LE.0 .OR. MaxWr.LE.0) 
     &   Call errMesFEM(4101, 'BilinearFormTemplate', 
     &                        'mesh numbers are negative')

      nR = 0
      MaxR = 3 * nE

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
         
c ... default pointers
      iIRE = 1


c ... order vertex indices in the connectivity table
c     Call order2D(3, nE, IPE)


c ... get template array
      iP1 = IPE(1, 1)
      iP2 = IPE(2, 1)
      iP3 = IPE(3, 1)

      Do i = 1, 3
         lbFloc(i) = 0
         lbPloc(i) = 0
      End do

c     not all information is available at this moment to populate iSYS
      nR = 0
      Call encodeISYS(1, iP1,iP2,iP3, 0,0,0, 0,0,0, nP,nR,nE, iSYS)

c     signaling a special request
      iSYS(1) = -1
      Call FEM2Dext(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3),
     &              lbE(1), lbFloc, lbPloc, DDATA, IDATA, iSYS,
     &              MaxSize, Aloc, Floc, ir, ic, 
     &              templateR, templateC)
 
      If(ir.GT.MaxSize .OR. ic.GT.MaxSize) 
     &  Call errMesFEM(2001, 'BilinearFormTemplate', 'MaxSize is small')

c     verifying templates structure
      Call verifyDOF(ir, templateR, flagR)
      Call verifyDOF(ic, templateC, flagC)

      If(.NOT.flagR .OR. .NOT.flagC) 
     &   Call errMesFEM(4101, 'BilinearFormTemplate', 
     &       'template[R,C] is invalid or breaks ordering rules')

c ... analyze rows
      irp = 0
      irr = 0
      ire = 0
      ire_dg = 0
      Do i = 1, ir
         If(ifXnode(templateR(i), Vdof)) Then
            irp = irp + 1
         Else If(ifXnode(templateR(i), Rdof) .OR. 
     &           ifXnode(templateR(i), RdofOrient)) Then
            irr = irr + 1
         Else If(ifXnode(templateR(i), Edof)) Then
            ire = ire + 1
         Else
            ire_dg = ire_dg + 1
         End if
      End do
      irp = irp / 3 
      irr = irr / 3 


c ... analyze columns
      icp = 0
      icr = 0
      ice = 0
      ice_dg = 0
      Do i = 1, ic
         If(ifXnode(templateC(i), Vdof)) Then
            icp = icp + 1
         Else If(ifXnode(templateC(i), Rdof) .OR.
     &           ifXnode(templateC(i), RdofOrient)) Then
            icr = icr + 1
         Else If(ifXnode(templateC(i), Edof)) Then
            ice = ice + 1
         Else
            ice_dg = ice_dg + 1
         End if
      End do
      icp = icp / 3
      icr = icr / 3


c ... memory distribution (lists P -> X)
      inEP  = 1
      iIEP  = inEP + nP
      iEnd  = iIEP + 3 * nE

      If(irp.GT.0 .AND. icp.GT.0) Then
         MaxL = 10 * nP

         inPP = iEnd
         iIPP = inPP + nP
         iiW  = iIPP + MaxL
         iEnd = iiW + 3 * nE + nP
         If(iEnd.GT.MaxWi) Goto 500

         Call listP2P(nP, nE, MaxL, IPE, iW(inPP), iW(iIPP), iW(iiW))

         iEnd = iIPP + iW(inPP + nP - 1)
      Else
         inPP = 1
         iIPP = 1
      End if

      If((irp.GT.0 .AND. icr.GT.0) .OR.
     &   (irr.GT.0 .AND. icp.GT.0)) Then
         iIRE = iEnd
         iEnd = iIRE + 3 * nE
         If(iEnd.GT.MaxWi) Goto 500

         Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

         inRP = iEnd
         iIRP = inRP + nP
         iiW  = iIRP + 4 * nR  
         iEnd = iiW + nR
         If(iEnd.GT.MaxWi) Goto 500

         MaxX = 4 * nR

         Call listConv(nP, nR, nE, iW(inEP), iW(iIEP), iW(iIRE), 
     &                 nX, MaxX, iW(inRP), iW(iIRP), iW(iiW))

         iEnd = iIRP + nX 

         inPR = iEnd
         iIPR = inPR + nR
         iEnd = iIPR + nX
         Call reverseMap(nP, nR, iW(inRP),iW(iIRP), iW(inPR),iW(iIPR))
      Else
         inRP = 1
         iIRP = 1

         inPR = 1
         iIPR = 1
      End if

      If(irp.GT.0 .AND. ice.GT.0) Then
         Call backReferences(nP, nE, 3,3, IPE, iW(inEP), iW(iIEP))
      End if


c ... memory distribution (lists R -> X)
      If(irr.GT.0 .AND. icr.GT.0) Then
         MaxL = 9 * nE

         iIRE = iEnd
         iEnd = iIRE + 3 * nE
         If(iEnd.GT.MaxWi) Goto 500

         Call listE2R(nP, nR, nE, IPE, iW(iIRE), iW(inEP), iW(iIEP))

         inRR = iEnd
         iIRR = inRR + MaxR
         iiW  = iIRR + MaxL
         iEnd = iiW + 9 * nE 
         If(iEnd.GT.MaxWi) Goto 500

         Call listR2R(nP, nR, nE, MaxL, IPE, iW(inRR),iW(iIRR), iW(iiW))

         iEnd = iIRR + iW(inRR + nR - 1)
      End if

      If((irr.GT.0 .AND. ice.GT.0) .OR.
     &   (ire.GT.0 .AND. icr.GT.0)) Then
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
         inER = 1
         iIER = 1
      End if


c ... create map E -> F
      iIFE = iEnd
      jnEP = iIFE + 3 * nE
      jIEP = jnEP + nP
      iEnd = jIEP + 3 * nE 
      If(iEnd.GT.MaxWi) Goto 500
      
      Call listE2F(nP, nF, nE, IPF, IPE, iW(iIFE), iW(jnEP), iW(jIEP))


c ... create map P -> F
      iIFP = iEnd
      inFP = iIFP + 2 * nF
      iEnd = inFP + nP
      If(iEnd.GT.MaxWi) Goto 500

      Call backReferences(nP, nF, 2, 2, IPF, iW(inFP), iW(iIFP))


      nRow = irp * nP + irr * nR + ire * nE
      nCol = icp * nP + icr * nR + ice * nE


 500  Continue
      If(iEnd.GT.MaxWi) Then
         Write(message,'(A,I10)')
     &      'The approximate size of iW is ', iEnd
         Call errMesFEM(1001, 'BilinearFormTemplate', message)
      End if

      If(nRow.GT.MaxF) Call errMesFEM(1013, 'BilinearFormTemplate',
     &                        'local parameter MaxF is small')

      iusage = (iEnd * 100) / MaxWi
      rusage = 0


c     setting global group (g) and degrees (d) pointers for rows
      iEnd = 1
      group = Group1
      Do n = 1, MaxNumGroups
         Call findGroupDOF(ir, templateR, group, i1, i2, label)
         If(i1.GT.0) dpr(i1) = iEnd
         gpr(n) = iEnd
         gnr(n) = 1
         gur(n) = label

         If(label.EQ.Vdof) Then
            iEnd = iEnd + nP
         Else If(label.EQ.Rdof .OR. label.EQ.RdofOrient) Then
            m = (i2 - i1 + 1) / 3
            gnr(n) = m
            iEnd = iEnd + nR * m
         Else If(label.EQ.Edof) Then
            m = i2 - i1 + 1
            iEnd = iEnd + nE * m
         End if
         group = group * 2
      End do

c     setting global group (g) and degrees (d) pointers for columns
      iEnd = 1
      group = Group1
      Do n = 1, MaxNumGroups
         Call findGroupDOF(ir, templateC, group, i1, i2, label)
         If(i1.GT.0) dpc(i1) = iEnd
         gpc(n) = iEnd
         gnc(n) = 1
         guc(n) = label

         If(label.EQ.Vdof) Then
            iEnd = iEnd + nP
         Else If(label.EQ.Rdof .OR. label.EQ.RdofOrient) Then
            m = (i2 - i1 + 1) / 3
            gnc(n) = m
            iEnd = iEnd + nR * m
         Else If(label.EQ.Edof) Then
            m = i2 - i1 + 1
            iEnd = iEnd + nE * m
         End if
         group = group * 2
      End do


c ... the symbolic assemble of the matrix 
      IA(1) = 1
      group = Group1
      Do m = 1, MaxNumGroups
         k = gpr(m)
         label = gur(m)

c  ...  adding connections points->others 
         If(label.EQ.Vdof) Then
            ip1 = 0
            ir1 = 0
            ie1 = 0

            Do n = 1, nP 
               ip2 = iW(inPP + n - 1) 
               ir2 = iW(inRP + n - 1) 
               ie2 = iW(inEP + n - 1) 

               iEnd = IA(k) + icp * (ip2 - ip1) 
     &                      + icr * (ir2 - ir1)
     &                      + ice * (ie2 - ie1)
               If(iEnd.GT.MaxA) Goto 600

               p = IA(k)
               Call fillJA(p, JA, gpc, gnc, guc, nR,
     &                     iW(iIPP), iW(iIRP), iW(iIEP),
     &                     icp, icr, ice, 
     &                     ip1, ir1, ie1, ip2, ir2, ie2)

               k = k + 1
               IA(k) = iEnd

               ip1 = ip2
               ir1 = ir2
               ie1 = ie2
            End do

c ... adding connections edges->others 
         Else If(label.EQ.Rdof .OR. label.EQ.RdofOrient) Then
            Do i = 1, gnr(m)
               k = gpr(m) + (i - 1) * nR
               ip1 = 0
               ir1 = 0
               ie1 = 0

               Do n = 1, nR 
                  ip2 = iW(inPR + n - 1)
                  ir2 = iW(inRR + n - 1) 
                  ie2 = iW(inER + n - 1) 

                  iEnd = IA(k) + icp * (ip2 - ip1)
     &                         + icr * (ir2 - ir1)
     &                         + ice * (ie2 - ie1)
                  If(iEnd.GT.MaxA) Goto 600

                  p = IA(k)
                  Call fillJA(p, JA, gpc, gnc, guc, nR,
     &                        iW(iIPR), iW(iIRR), iW(iIER),
     &                        icp, icr, ice, 
     &                        ip1, ir1, ie1, ip2, ir2, ie2)

                  k = k + 1
                  IA(k) = iEnd

                  ip1 = ip2
                  ir1 = ir2
                  ie1 = ie2
               End do
            End do

c ... adding connections elements->others 
         Else If(label.EQ.Edof) Then
            Do n = 1, nE
               IEE(1) = n

               Do l = 1, ire
                  ip1 = 3 * (n - 1)
                  ir1 = 3 * (n - 1)
                  ip2 = ip1 + 3
                  ir2 = ir1 + 3

                  iEnd = IA(k) + icp * 3 + icr * 3 + ice
                  If(iEnd.GT.MaxA) Goto 600

                  p = IA(k)
                  Call fillJA(p, JA, gpc, gnc, guc, nR,
     &                  IPE, iW(iIRE), IEE,
     &                  icp, icr, ice, 
     &                  ip1, ir1, 0, ip2, ir2, 1)

                  k = k + 1
                  IA(k) = iEnd
               End do
            End do
         End if
      End do

      iEnd = IA(nRow + 1) - 1
 600  Continue 
      If(iEnd.GT.MaxA) Call errMesFEM(1014,
     &     'BilinearFormTemplate', 'local parameter MaxA is small')


C ... comply with the AMG format
      If(flagFA) Then
         Do 700 n = 1, nRow
            Do j = IA(n), IA(n + 1) - 1
               If(JA(j).EQ.n) Then
                  Call swapii(JA(j), JA(IA(n)))
                  Goto 700
               End if
            End do
 700     Continue
      End if


C ... fill in the sparsity structure
      Do n = 1, nRow
         F(n) = 0D0
      End do

      Do n = 1, IA(nRow + 1) - 1
         A(n) = 0D0
      End do

      Do n = 1, nE
         ip1 = IPE(1, n)
         ip2 = IPE(2, n)
         ip3 = IPE(3, n)

         If(iIRE.GT.1) Then
            ir1 = iW(iIRE + 3*(n-1))
            ir2 = iW(iIRE + 3*(n-1) + 1)
            ir3 = iW(iIRE + 3*(n-1) + 2)
         Else
            ir1 = 0
            ir2 = 0
            ir3 = 0
         End if

c  ...   save edge labels
         Do i = 1, 3
            iFt = iW(iIFE + 3 * (n - 1) + i - 1)
            lbFloc(i) = 0
            If(iFt.GT.0) lbFloc(i) = lbF(iFt)
         End do

         Do i = 1, 3
            iPt = IPE(i, n)
            lbPloc(i) = lbP(iPt)
         End do

         Call encodeISYS(n, ip1,ip2,ip3, ir1,ir2,ir3, 
     &                   0,0,0, nP,nR,nE, iSYS)

         Call FEM2Dext(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3),
     &                 lbE(n), lbFloc, lbPloc, DDATA, IDATA, iSYS,
     &                 MaxSize, Aloc, Floc, ir, ic, 
     &                 dummyR, dummyC)

c  ...  calculate orientation of edges
         Do i = 1, 3
            orient(i) = 1
         End do
         If(iP1.GT.iP2) orient(1) = 0
         If(iP2.GT.iP3) orient(2) = 0
         If(iP3.GT.iP1) orient(3) = 0


c  ...  create arrays of row indices 
         group = Group1
         Do m = 1, MaxNumGroups
            Call findGroupDOF(ir, templateR, group, i1, i2, label)
            ishift = gpr(m) - 1
            If(label.EQ.Vdof) Then
               indR(i1)     = ishift + iP1 
               indR(i1 + 1) = ishift + iP2 
               indR(i1 + 2) = ishift + iP3 
            Else If(label.EQ.Rdof .OR. label.EQ.RdofOrient) Then
               l = (i2 - i1 + 1) / 3
               Do i = 1, l 
                  Do k = 1, 3
                     ixr(k) = (i - 1) * orient(k)
     &                      + (l - i) * (1 - orient(k)) 
                  End do

                  k = iIRE + 3 * (n - 1)
                  indR(i1)     = ishift + ixr(1) * nR + iW(k) 
                  indR(i1 + 1) = ishift + ixr(2) * nR + iW(k + 1)
                  indR(i1 + 2) = ishift + ixr(3) * nR + iW(k + 2)
                  i1 = i1 + 3
               End do
            Else If(label.EQ.Edof) Then
               l = i2 - i1 + 1
               Do i = i1, i2
                  indR(i) = ishift + (n - 1) * l + (i - i1 + 1)
               End do
            End if
            group = group * 2
         End do

c  ...  create arrays of column indices 
         group = Group1
         Do m = 1, MaxNumGroups
            Call findGroupDOF(ic, templateC, group, i1, i2, label)
            ishift = gpc(m) - 1
            If(label.EQ.Vdof) Then
               indC(i1)     = ishift + iP1 
               indC(i1 + 1) = ishift + iP2 
               indC(i1 + 2) = ishift + iP3 
            Else If(label.EQ.Rdof .OR. label.EQ.RdofOrient) Then
               l = (i2 - i1 + 1) / 3
               Do i = 1, l 
                  Do k = 1, 3
                     ixr(k) = (i - 1) * orient(k)
     &                      + (l - i) * (1 - orient(k)) 
                  End do

                  k = iIRE + 3 * (n - 1)
                  indC(i1)     = ishift + ixr(1) * nR + iW(k) 
                  indC(i1 + 1) = ishift + ixr(2) * nR + iW(k + 1)
                  indC(i1 + 2) = ishift + ixr(3) * nR + iW(k + 2)
                  i1 = i1 + 3
               End do
            Else If(label.EQ.Edof) Then
               l = i2 - i1 + 1
               Do i = i1, i2
                  indC(i) = ishift + (n - 1) * l + (i- i1 + 1) 
               End do
            End if
            group = group * 2
         End do


c  ...  correct orientation of the basis functions in the local matrix
c  ...  R-unknowns: the global orientation is from bigger E to smaller E
         i = 1
         Do While(i.LE.ir)
            If(ifXnode(templateR(i), RdofOrient)) Then
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
            If(ifXnode(templateC(i), RdofOrient)) Then
               Do k = 0, 2
                  iRt = iW(iIRE + 3*(n-1) + k)
                  iEt = findY(iRt, iW(iIER), iW(inER), n)

                  If(n.LT.iEt) Then
                     Do j = 1, ir
                        Aloc(j, i) = -Aloc(j, i)
                     End do
                  End if

                  i = i + 1
               End do
            Else
               i = i + 1
            End if
         End do


c  ...  assemble the right hand side
         Do i = 1, ir
            iP = indR(i)
            F(iP) = F(iP) + Floc(i)
         End do


c  ...  assemble the elemental matrix
         Do i = 1, ir
            iP = indR(i)

            Do 200 j = 1, ic
               iQ = indC(j)

               Do k = IA(iP), IA(iP + 1) - 1
                  If(JA(k).EQ.iQ) Then
                     A(k) = A(k) + Aloc(i, j)
                     Goto 200
                  End if
               End do
             
               Call errMesFEM(6001, 'BilinearFormTemplate', 
     &                       'System error: symbolic assembling wrong')
 200        Continue
         End do
      End do


c ... convert row format to column format using ugly method
      non_zero = IA(nRow + 1) - 1 

      If(flagFC) Then
         iIA  = 1
         iJA  = iIA + nRow + 1
         iEnd = iJA + non_zero

         jEnd = non_zero

         If(iEnd.GT.MaxWi) Goto 500 
         If(jEnd.GT.MaxWr) Call errMesFEM(1001, 
     &     'BilinearFormTemplate', 'Not enough Real*8 memory')

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
      Subroutine fillJA(p, JA, 
     &                  gpc, gnc, guc, nR, IPX, IRX, IEX,
     &                  icp, icr, ice, 
     &                  ip1, ir1, ie1, ip2, ir2, ie2)
C ======================================================================
C  Technical subroutine for incorporating a part of the local 
C  sparsity structure into the global one.
C ======================================================================
      implicit none
      Include 'fem2Dtri.fd'

      Integer p, JA(*), gpc(*), gnc(*), guc(*), icp, icr, ice
      Integer nR, IPX(*), IRX(*), IEX(*), ip1, ir1, ie1, ip2, ir2, ie2

      Integer j, k, l, n, ishift, label

C ======================================================================
      Do n = 1, MaxNumGroups
         ishift = gpc(n) - 1
         label = guc(n)
   
         If(label.EQ.Vdof) Then
            Do j = ip1 + 1, ip2
               JA(p) = ishift + IPX(j)
               p = p + 1
            End do

         Else If(label.EQ.Rdof .OR. label.EQ.RdofOrient) Then
            Do k = 1, gnc(n) 
               ishift = gpc(n) + (k - 1) * nR - 1

               Do j = ir1 + 1, ir2
                  JA(p) = ishift + IRX(j)
                  p = p + 1
               End do
            End do

         Else If(label.EQ.Edof) Then
            Do j = ie1 + 1, ie2
               Do l = 1, ice
                  JA(p) = ishift + (IEX(j) - 1) * ice + l
                  p = p + 1
               End do
            End do
         End If
      End do

      Return
      End



C ================================================================
      Integer Function findY(iX, IYX, nYX, iY)
C ================================================================
C findY returns pointer to object Y, other than iY, in the sublist 
c of X(iX) -> {iY1, iY2, ..., iYk} of map X -> Y.
C If not found, findY = 0.
C ================================================================
      implicit none
      Integer  iX, iY, IYX(*), nYX(*)
      Integer  i, ib, ie, iYt

      If(iX.EQ.1) Then
         ib = 1
      Else
         ib = nYX(iX - 1) + 1
      End if
      ie = nYX(iX)

      Do i = ib, ie
         iYt = IYX(i)
         If(iYt.NE.iY) Then
            findY = iYt
            Goto 1000
         End if
      End do

      findY = 0
 1000 Return
      End


