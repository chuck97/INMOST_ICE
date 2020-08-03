C ======================================================================
      Subroutine assemble_rhs(nv, vrt, nt, tri, nv2, vrt2, nt2, tri2,
     &                        nv12, vrt12, nt12, tri12, parents, 
     &                        operatorA, FEMtypeA, operatorB, FEMtypeB,
     &                        U1, U2, MaxWi, iW)
C ======================================================================
C Routine calculates the right-hand side vector U1 = M12 U2, where
C Ui is the function on the i-th mesh and M12 is the mass or stiffness
C (for interpolation in general normed spaces) matrix corresponing
C to cross-integration of basis functions v1_i and v2_j on different 
C meshes:
C
C    \Int_\Omega v1_i v2_j dx
C
C The matrix M12 is not assembled and multiplication is done using the
C meta-mesh.
C
C At the moment, the integration weigth is a constant, but it may be 
C generalized to spatially-varying weigths.
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'
      include 'assemble.fd'

      Integer  nv, nt, tri(3, *), nv2, nt2, tri2(3, *)
      Real*8   vrt(2, *), vrt2(2, *)

      Integer  nv12, nt12, tri12(3, *), parents(2, *)
      Real*8   vrt12(2, *), U1(*), U2(*)

      Integer  operatorA, FEMtypeA, operatorB, FEMtypeB
      Integer  MaxWi, iW(*)


C LOCAL VARIABLES
      Integer  iDATA(1), iSYS(MAXiSYS), label, order, nRow, nCol
      Real*8   dDATA(1), M12(MaxSize, MaxSize), s, tvol

      Integer  dof(MaxSize), dof2(MaxSize), ip(MaxSize), iq(MaxSize)
      Integer  ANI_Dnull, orderDOF
      EXTERNAL ANI_Dnull, orderDOF

      Integer  i,j,k1,k2, n,n1,n2, iv1,iv2,iv3, ip1,ip2,ip3, iq1,iq2,iq3
      Integer  nr, nr2, iIRE, iIRE2, inEP, iIEP, iEnd
      
C ======================================================================
      label = 1   !dummy: is not used by the ANI_Dnull routine
      order = orderDOF(FEMtypeA) + orderDOF(FEMtypeB)  !maximal quadrature

c memory allocation for mapping
      iIRE  = 1
      iIRE2 = iIRE  + 3 * nt
      inEP  = iIRE2 + 3 * nt2
      iIEP  = inEP  + max(nv, nv2)
      iEnd  = iIEP  + 3 * max(nt, nt2)

      If(iEnd.GT.MaxWi) Call errMesMBA(1001, 
     &   'assemble_rhs', 'Not enough integer memory')


c maps E->R for both meshes
      Call listE2R(nv,  nr,  nt,  tri,  iW(iIRE),  iW(inEP), iW(iIEP))
      Call listE2R(nv2, nr2, nt2, tri2, iW(iIRE2), iW(inEP), iW(iIEP))


c create local -> global maps
      Call listDOF(FEMtypeA, k1, dof)  
      Call listDOF(FEMtypeB, k2, dof2)  


c assemble the righ-hand side
      Call dof2nRow(k1, dof, nv, nr, nt, nRow)

      tvol = 0.0d0
      Do n = 1, nRow
         U1(n) = 0
      End do

      Do n = 1, nt12
         iv1 = tri12(1, n)
         iv2 = tri12(2, n)
         iv3 = tri12(3, n)

         n1 = parents(1, n)
         n2 = parents(2, n)

c vertices of a triangle on mesh #1
         ip1 = tri(1, n1)
         ip2 = tri(2, n1)
         ip3 = tri(3, n1)

c vertices a triangle on mesh #2
         iq1 = tri2(1, n2)
         iq2 = tri2(2, n2)
         iq3 = tri2(3, n2)

c calling integration routine from library aniFEM
         Call fem2Dsub(vrt12(1,iv1), vrt12(1,iv2), vrt12(1,iv3),  
     &                 vrt(  1,ip1), vrt(  1,ip2), vrt(  1,ip3),  
     &                 vrt2( 1,iq1), vrt2( 1,iq2), vrt2( 1,iq3),  
     &                 operatorA, FEMtypeA, operatorB, FEMtypeB,
     &                 label, ANI_Dnull, dDATA, iDATA, iSYS, order,
     &                 MaxSize, M12, nRow, nCol)

c        Write(*,'(A,2I4)') 'Elemental matrix for:', FEMtypeA, FEMtypeB
c        Do i = 1, nRow
c           Write(*, '(100E14.6)') (M12(i, j), j = 1, nCol)
c        End do

c integrating over triangle iv1-iv2-iv3
         Call dof2global(k1, dof,  nv,  nr,  n1, tri,  iW(iIRE),  ip)
         Call dof2global(k2, dof2, nv2, nr2, n2, tri2, iW(iIRE2), iq)

         Do i = 1, nRow
            s = 0D0
            Do j = 1, nCol
               s = s + M12(i, j) * U2(iq(j))
               tvol = tvol + M12(i, j)
            End do
            U1(ip(i)) = U1(ip(i)) + s
         End do
      End do
c     Write(*,'(A,E14.8)') 'PRJ: DEBUG of tvol =', tvol

      Return
      End



C ======================================================================
      Subroutine dof2global(N, dof, nP, nR, kE, IPE, IRE, ind)
C ======================================================================
C Default rules are used for local->global mapping of DOF for element iE
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer N, dof(*), nP, nR, IPE(3, *), IRE(3, *), ind(*), kE

      Integer i,j, ip,ir,ie, isR,isE

C ======================================================================
c cound gropus of unknowns
      ip = 0
      ir = 0
      ie = 0
      Do i = 1, N
         If(dof(i).EQ.Vdof) Then
            ip = ip + 1
         Else If(dof(i).EQ.Rdof .OR. dof(i).EQ.RdofOrient) Then
            ir = ir + 1
         Else If(dof(i).EQ.Edof) Then
            ie = ie + 1
         End if
      End do
      ip = ip / 3 
      ir = ir / 3 


c define global shifts
      isR = ip * nP
      isE = isR + ir * nR

      ip = 0
      ir = 0

      i = 1
      Do While(i.LE.N)
         If(dof(i).EQ.Vdof) Then
            ind(i)     = ip * nP + IPE(1, kE) 
            ind(i + 1) = ip * nP + IPE(2, kE)
            ind(i + 2) = ip * nP + IPE(3, kE)

            i  = i + 3
            ip = ip + 1

         Else If(dof(i).EQ.Rdof .OR. dof(i).EQ.RdofOrient) Then
            ind(i)     = isR + ir * nR + IRE(1, kE)
            ind(i + 1) = isR + ir * nR + IRE(2, kE)
            ind(i + 2) = isR + ir * nR + IRE(3, kE)

            i  = i + 3
            ir = ir + 1

         Else If(dof(i).EQ.Edof) Then
            Do j = 1, ie
               ind(i + j - 1) = isE + ie * (kE - 1) + j 
            End do

            i = i + ie
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine dof2nRow(N, dof, nP, nR, nE, nRow)
C ======================================================================
C Calculates the size of the array U1.
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer N, dof(*), nP, nR, nE, nRow
      Integer i, ip,ir,ie
C ======================================================================
c cound gropus of unknowns
      ip = 0
      ir = 0
      ie = 0
      Do i = 1, N
         If(dof(i).EQ.Vdof) Then
            ip = ip + 1
         Else If(dof(i).EQ.Rdof .OR. dof(i).EQ.RdofOrient) Then
            ir = ir + 1
         Else If(dof(i).EQ.Edof) Then
            ie = ie + 1
         End if
      End do
      ip = ip / 3 
      ir = ir / 3 

      nRow = ip * nP + ir * nR + ie * nE

      Return
      End



