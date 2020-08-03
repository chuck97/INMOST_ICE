C ======================================================================
      Subroutine vectorDOF(N, dof, vector)
C ======================================================================
C Pre-process input degrees of freedom
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer  N, dof(*), vector, i
      Logical  ifXnode
 
      vector = 0

      Do i = 1, N
         if(ifXnode(dof(i), VectorY)) vector = 1
         Call delXnode(dof(i), VectorY)
      End do

      Return
      End



C ======================================================================
      Subroutine verifyDOF(N, dof, flag)
C ======================================================================
C Pre-process input degrees of freedom
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer  N, dof(*), i,i1,i2, k, label, group, Xdof, GroupX
      Logical  flag, ifXnode
 
      flag = .FALSE.

c     check that dofs have valid labels
      Xdof = Vdof + Rdof + Edof + RdofOrient
      Do i = 1, N
         If(IAND(dof(i), Xdof).EQ.0) Return
      End do

c     check groups are either not defined or defined everywhere
      GroupX = Group1 * (2 ** MaxNumGroups) - Group1
      group = IAND(dof(1), GroupX)
      Do i = 1, N
         i1 = IAND(dof(i), GroupX)
         If(group.EQ.0 .AND. i1.GT.0) Return
         If(group.GT.0 .AND. i1.EQ.0) Return
      End do
      If(group.EQ.0) group = Group1 / 2

c     check length of groups
      i1 = 1
      Do While(i1.LE.N)
         Call groupDOF(N, dof, i1, i2, label, group)
         If(label.NE.Edof .AND. mod(i2 - i1 + 1, 3).NE.0) Return
         i1 = i2 + 1
      End do

c     groups are created or given. Now we check their uniqueness.
      i1 = 1
      Do While(i1.LE.N)
         Call groupDOF(N, dof, i1, i2, label, group)
         If(.NOT.ifXnode(GroupX, group)) Return 
         Call delXnode(GroupX, group)
         i1 = i2 + 1
      End do

      flag = .TRUE.
      Return
      End



C ======================================================================
      Subroutine groupDOF(N, dof, i1, i2, label, group)
C ======================================================================
C Identify a group of common degrees of freedom.
C Parameter group enters with the previous id and returns new id.
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer  N, dof(*), label, group, i, i1, i2, Xdof, GroupX
      Logical  ifXnode
 
      Xdof = Vdof + Rdof + Edof + RdofOrient
      label = IAND(dof(i1), Xdof)

c     mark the next group if none exists
      GroupX = Group1 * (2 ** MaxNumGroups) - Group1
      i = group
      group = IAND(dof(i1), GroupX)
      If(group.EQ.0) Then
         group = i * 2 

         If(label.EQ.Edof) Then
            i2 = i1
         Else
            i2 = i1 + 2
         End if

         Do i = i1, i2
            Call addXnode(dof(i), group)
         End do
         Return
      End if

c     identify the existing group
      i2 = N
      Do i = i1 + 1, N 
         If(dof(i).NE.dof(i1)) Then
            i2 = i - 1
            Return
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine findGroupDOF(N, dof, group, i1, i2, label)
C ======================================================================
C Find a group of DOFs with the given id. 
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Integer  N, dof(*), label, group, i, i1, i2, Xdof
      Logical  ifXnode

      Xdof = Vdof + Rdof + Edof + RdofOrient

      i1 = 0
      Do i = 1, N
         If(ifXnode(dof(i), group)) Then
            i1 = i
            Goto 100
         End if
      End do

      i2 = 0
      label = 0
      Return

 100  Continue
      label = IAND(dof(i1), Xdof)
      i2 = N
      Do i = i1 + 1, N
         If(dof(i).NE.dof(i1)) Then
            i2 = i - 1
            Return
         End if
      End do

      Return
      End



C ======================================================================
      Subroutine listDOF(FEMtype, N, dof)
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'
 
      Integer FEMtype, N, dof(*), i

      If(FEMtype.EQ.FEM_P0) Then
         N = 1
         dof(1) = Edof

      Else If(FEMtype.EQ.FEM_P1) Then
         N = 3
         Do i = 1, 3
            dof(i) = Vdof
         End do

      Else If(FEMtype.EQ.FEM_P2) Then
         N = 6
         Do i = 1, 3
            dof(i)   = Vdof
            dof(i+3) = Rdof
         End do

      Else If(FEMtype.EQ.FEM_P3) Then
         N = 10
         Do i = 1, 3
            dof(i)   = Vdof
            dof(i+3) = Rdof
            dof(i+6) = Rdof
         End do
         dof(10) = Edof
       
      Else If(FEMtype.EQ.FEM_P4) Then
         N = 15
         Do i = 1, 3
            dof(i)    = Vdof
            dof(i+ 3) = Rdof
            dof(i+ 6) = Rdof
            dof(i+ 9) = Rdof
            dof(i+12) = Edof
         End do
       
      Else If(FEMtype.EQ.FEM_P1vector) Then
         N = 6
         Do i = 1, 6
            dof(i) = Vdof
         End do

      Else If(FEMtype.EQ.FEM_P2vector) Then
         N = 12
         Do i = 1, 3
            dof(i)   = Vdof
            dof(i+3) = Rdof
            dof(i+6) = Vdof
            dof(i+9) = Rdof
         End do

      Else If(FEMtype.EQ.FEM_P2reduced) Then
         N = 9
         Do i = 1, 3
            dof(i)   = Vdof
            dof(i+3) = Rdof
            dof(i+6) = Vdof
         End do

      Else If(FEMtype.EQ.FEM_MINI) Then
         N = 8
         Do i = 1, 3
            dof(i)   = Vdof
            dof(i+4) = Vdof
         End do
         dof(4) = Edof
         dof(8) = Edof

      Else If(FEMtype.EQ.FEM_RT0) Then
         N = 3
         Do i = 1, 3
            dof(i) = RdofOrient
         End do

      Else If(FEMtype.EQ.FEM_BDM1) Then
         N = 6
         Do i = 1, 6
            dof(i) = RdofOrient
         End do

      Else If(FEMtype.EQ.FEM_CR1) Then
         N = 3
         Do i = 1, 3
            dof(i) = Rdof
         End do

      Else If(FEMtype.EQ.FEM_CR1vector) Then
         N = 6
         Do i = 1, 6
            dof(i) = Rdof
         End do

      Else
         N = 0
      End if

      Return 
      End



C ======================================================================
      Subroutine enumDOF(FEMtype, N, dof, i1)
C ======================================================================
C Local enumeration of degrees of freedom starting from the vertex i1.
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'
 
      Integer FEMtype, N, dof(*), i1
      Integer iref(4), i2, i
      DATA    iref/1,2,3,1/

      i2 = i1
 
      If(FEMtype.EQ.FEM_P1) Then
         N = 3
         Do i = 1, 3
            dof(i) = i2
            i2 = iref(i2 + 1)
         End do

      Else If(FEMtype.EQ.FEM_P2) Then
         N = 6
         Do i = 1, 3
            dof(i)   = i2
            dof(i+3) = i2 + 3
            i2 = iref(i2 + 1)
         End do

      Else If(FEMtype.EQ.FEM_P3) Then
         N = 10
         Do i = 1, 3
            dof(i)   = i2
            dof(i+3) = i2 + 3
            dof(i+6) = i2 + 6
            i2 = iref(i2 + 1)
         End do
         dof(10) = 10
       
      Else If(FEMtype.EQ.FEM_P4) Then
         N = 15
         Do i = 1, 3
            dof(i)    = i2
            dof(i+ 3) = i2 + 3
            dof(i+ 6) = i2 + 6
            dof(i+ 9) = i2 + 9
            dof(i+12) = i + 12
            i2 = iref(i2 + 1)
         End do

      Else
         N = 0
      End if

      Return 
      End



C ======================================================================
      Integer Function orderDOF(FEMtype)
C ======================================================================
C Calculate maximal order for integration
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'
 
      Integer FEMtype

      If(FEMtype.EQ.FEM_P0) Then
         orderDOF = 1

      Else If(FEMtype.EQ.FEM_P1) Then
         orderDOF = 1

      Else If(FEMtype.EQ.FEM_P2) Then
         orderDOF = 2

      Else If(FEMtype.EQ.FEM_P3) Then
         orderDOF = 3
       
      Else If(FEMtype.EQ.FEM_P4) Then
         orderDOF = 4
       
      Else If(FEMtype.EQ.FEM_P1vector) Then
         orderDOF = 1

      Else If(FEMtype.EQ.FEM_P2vector) Then
         orderDOF = 2

      Else If(FEMtype.EQ.FEM_P2reduced) Then
         orderDOF = 2

      Else If(FEMtype.EQ.FEM_MINI) Then
         orderDOF = 3

      Else If(FEMtype.EQ.FEM_RT0) Then
         orderDOF = 1

      Else If(FEMtype.EQ.FEM_BDM1) Then
         orderDOF = 2

      Else If(FEMtype.EQ.FEM_CR1) Then
         orderDOF = 1

      Else If(FEMtype.EQ.FEM_CR1vector) Then
         orderDOF = 1

      Else
         orderDOF = 1
      End if

      Return 
      End


