C ======================================================================
      Logical function tangled(iE1, iE2, XYP, IPE)
C ======================================================================
C Routines returns .TRUE. is triangles iE1 and iE2 are inverted.
C ======================================================================
      implicit none
 
      Real*8  XYP(2, *)
      Integer iE1, iE2, IPE(3, *)
      
      Real*8  tri_area2, v1, v2
      Logical check22

      Integer  tri_orient
      EXTERNAL tri_orient

C LOCAL VARIABLES
      Integer  i1,i2,i3, j1,j2,j3, k1,k2, iP1,iP2,iP3, jP1,jP2,jP3
      Integer  iref(4)
      DATA     iref /1,2,3,1/

C ======================================================================
      tangled = .FALSE.

      Do 20 i1 = 1, 3
         i2 = iref(i1 + 1)

         iP1 = IPE(i1, iE1)
         iP2 = IPE(i2, iE1)

         Do j1 = 1, 3
            j2 = iref(j1 + 1)

            jP1 = IPE(j1, iE2)
            jP2 = IPE(j2, iE2)

            If(check22(iP1, iP2, jP1, jP2)) Then
               i3  = iref(i2 + 1)
               iP3 = IPE(i3, iE1)

               j3  = iref(j2 + 1)
               jP3 = IPE(j3, iE2)

c              v1 = tri_area2(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
c              v2 = tri_area2(XYP(1, iP1), XYP(1, iP2), XYP(1, jP3))
               k1 = tri_orient(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3))
               k2 = tri_orient(XYP(1, iP1), XYP(1, iP2), XYP(1, jP3))

               If(k1 * k2.GE.0) Then
                  tangled = .TRUE.
                  Goto 9000 
               End if
            End if
         End do
 20   Continue

 9000 Return
      End



C ======================================================================
      Subroutine orderijk(i, j, k)
C ======================================================================
C This routines sets i, j, k in non-decreasing order
C ======================================================================
      If(i.LT.j) Then
         n = j
         j = i
         i = n
      End if

      If(j.LT.k) Then
         n = k
         k = j
         j = n
      End if

      If(i.LT.j) Then
         n = j
         j = i
         i = n
      End if
      Return
      End

