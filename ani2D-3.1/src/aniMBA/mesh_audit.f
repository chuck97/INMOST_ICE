C ======================================================================
      Subroutine meshAudit(nP, XYP, 
     &                     nF, IPF, Crv, lbC, CrvFunction,
     &                     nE, IPE, lbE, 
     &                     ICP, IFE, IEE, rRmax, rRavg, status)
C ======================================================================
      implicit none
      include 'colors.fd'
      include 'status.fd'
C ======================================================================
C Routines checks topology of the input and the oputput meshes.
C ======================================================================
C group (M)
      Integer  nP
      Real*8   XYP(2, *)

      Integer  nF, IPF(2, *), lbC(*)
      Real*8   Crv(2, *)
      EXTERNAL CrvFunction

      Integer  nE, IPE(3, *), lbE(*)

C group (M-EXT)
      Integer  ICP(*), IFE(3, *), IEE(3, *), status
      Real*8   rRmax, rRavg

C LOCAL VARIABLES
      Integer  ip(4)
      Real*8   XYPs(2), err, rOut, rIn, sqrEdge, tri_area2
      Logical  ifXnode, check22

      Integer  i,k,n, i1,i2,i3, j1,j2,j3, iPt, iP1,iP2,iP3, jP1,jP2,jP3
      Integer  k1,k2, iF, iF1,iF2,iF3, iC, iE, iE1,iE2,iE3, nTet, iERR
      Real*8   v1, v2, rR

      Integer  tri_orient
      EXTERNAL tri_orient

C ======================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

c ... check color of mesh points
      Do n = 1, nP
         If(ICP(n).LE.0) 
     &      Call errMesMBA(5001, 'meshAudit', 'wrong point color')

         If(ifXnode(ICP(n), jInode) .AND. 
     &      ifXnode(ICP(n), jBnode)) 
     &      Call errMesMBA(5015, 'meshAudit', 'color contradiction')
      End do


c ... check labels of mesh faces
      Do n = 1, nF
         If(IPF(1, n).LE.0) Call errMesMBA(5002, 'meshAudit', 
     &                                     'wrong map edge -> points')
      End do


c ... check elements
      Do n = 1, nE
         nTet = n

         Do i = 1, 3
            If(IPE(1, n).LE.0) Call errMesMBA(5004, 'meshAudit', 
     &                              'wrong map element -> points')
         End do

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         If(iP1.EQ.iP2 .OR. iP1.EQ.iP3 .OR. iP2.EQ.iP3) Then
            iERR = 5006
            Goto 500
         End if

         If(lbE(n).LE.0) Call errMesMBA(5016, 'meshAudit', 
     &                  'element identificator (lbE) is not positive')

         iF1 = IFE(1, n)
         iF2 = IFE(2, n)
         iF3 = IFE(3, n)

         If(iF1.EQ.iF2 .AND. iF1.GT.0 .OR.
     &      iF1.EQ.iF3 .AND. iF1.GT.0 .OR.
     &      iF2.EQ.iF3 .AND. iF2.GT.0) Then
            iERR = 5007
            Goto 500
         End if


         iE1 = IEE(1, n)
         iE2 = IEE(2, n)
         iE3 = IEE(3, n)

         If(iE1.EQ.iE2 .AND. iE1.NE.0 .OR.
     &      iE1.EQ.iE3 .AND. iE1.NE.0 .OR.
     &      iE2.EQ.iE3 .AND. iE2.NE.0) Then
            iERR = 5008
            Goto 500
         End if


         Do 20 i1 = 1, 3
            iF = IFE(i1, n)
            iE = IEE(i1, n)
            If(iF.EQ.0 .AND. iE.EQ.0 .AND. nF.GT.0) Then
               iERR = 5009
               Goto 500
            End if

            i2 = ip(i1 + 1)

            iP1 = IPE(i1, n)
            iP2 = IPE(i2, n)

            If(iF.NE.0) Then
               jP1 = IPF(1, iF)
               jP2 = IPF(2, iF)

               If(.NOT.check22(iP1, iP2, jP1, jP2)) Then
                  iERR = 5010
                  Goto 500
               End if
            End if

            If(iE.NE.0) Then
               If(iF.NE.0) Then
                  Do j1 = 1, 3
                     If(IFE(j1, iE).EQ.iF) Goto 10
                  End do

                  iERR = 5011
                  Goto 500
               End if

  10           Do j1 = 1, 3
                  j2 = ip(j1 + 1)

                  jP1 = IPE(j1, iE)
                  jP2 = IPE(j2, iE)

                  If(check22(iP1, iP2, jP1, jP2)) Then
                     If(IEE(j1, iE).NE.n) Then
                        iERR = 5012
                        Goto 500
                     End if

                     i3  = ip(i2 + 1)
                     iP3 = IPE(i3, n)

                     j3  = ip(j2 + 1)
                     jP3 = IPE(j3, iE)

c                    v1 = tri_area2(XYP(1,iP1), XYP(1,iP2), XYP(1,iP3))
c                    v2 = tri_area2(XYP(1,iP1), XYP(1,iP2), XYP(1,jP3))
                     k1 = tri_orient(XYP(1,iP1), XYP(1,iP2), XYP(1,iP3))
                     k2 = tri_orient(XYP(1,iP1), XYP(1,iP2), XYP(1,jP3))

c  ...   check for inverted elements
                     If(k1 * k2.GE.0 .AND.
     &                  .NOT.ifXnode(status, ANIUntangleMesh)) Then
                        iERR = 5013
                        Goto 500
                     End if
                     Goto 20
                  End if
               End do

               iERR = 5014
               Goto 500
            End if
 20      Continue
      End do


c ... check the parametrization
      Do n = 1, nF
         iC = lbC(n)

         If(iC.GT.0) Then
            Do i = 1, 2
               iPt = IPF(i, n)

               Call aniCrv(Crv(i, n), XYPs, iC, CrvFunction)
               err = sqrEdge(XYP(1, iPt), XYPs)
c              If(err.GE.1D-16) 
c    &            write(*,*) n,i,iC,XYP(1,iPt),XYP(2,iPt),XYPs

               If(err.GE.1D-16)  Call errMesMBA(5006, 'meshAudit', 
     &                                          'wrong parametrization')
            End do
         End if
      End do   


c ... compute R / r
      rRmax = 0D0
      rRavg = 0D0
      Do n = 1, nE
         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)

         Call RandR(XYP(1, iP1), XYP(1, iP2), XYP(1, iP3), rOut, rIn)

         rR = rOut / rIn
         rRmax = max(rRmax, rR)
         rRavg = rRavg + rR
      End do
      rRavg = rRavg / nE

      Return


 500  Write(*, 5000) nTet, iERR, (IPE(i, nTet), i = 1, 3),
     &                           (IFE(i, nTet), i = 1, 3),
     &                           (IEE(i, nTet), i = 1, 3),
     &                           (ICP(IPE(i, nTet)), i = 1, 3)
      Do k = 1, 3
         iE = IEE(k, nTet)
         If(iE.GT.0) Then
            Write(*, 5000) iE, iERR, (IPE(i, iE), i = 1, 3),
     &                               (IFE(i, iE), i = 1, 3),
     &                               (IEE(i, iE), i = 1, 3),
     &                               (ICP(IPE(i, iE)), i = 1, 3)
         End if
      End do
      Do k = 1, 3
         iF = IFE(k, nTet)
         If(iF.GT.0) Then
            Write(*, 5004) iF, nF, (IPF(i, iF), i = 1, 2),
     &                             (ICP(IPF(i, max(1, iF))), i = 1, 2)
         End if
      End do

      Call errMesMBA(iERR, 'meshAudit', 'triangles are wrong')

 5000 Format('Error in checking triangle =', I6, '   iERR=', I5, /,
     &       'Points    =', 3I6, /,
     &       'Edges     =', 3I6, /,
     &       'Triangles =', 3I6, /,
     &       'colors    =', 3I6, /)

 5004 Format('Edge ', I6, '  (', I5, ')  of the bad triangle', /,
     &       'Points =', 2I6, /,
     &       'colors =', 2I6, /)

      End
