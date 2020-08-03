C ======================================================================
      Subroutine debug(
     &           nP, nF, nE,
     &           XYP, IPF, lbC, IPE, IFE, IEE, lbE,
     &           calCrv, parCrv, iFnc,
     &           IHolP, IHolF, IHolE,
     &           ICP, status, nQItr,
     &           flagMOVE, flagSWAP, flagCLPS1, flagCLPS2, flagINSRT)
C ======================================================================
      include 'colors.fd'
      include 'status.fd'
C ======================================================================
C Routines checks topology of the input and the oputput meshes. The
C size of working array is max(nP+nPholes, nF+nFholes, nE+nEholes).
C ======================================================================
C group (M)
      Integer IPF(2, *), lbC(*), IPE(3, *), IFE(3, *), IEE(3, *), lbE(*)
      Real*8  XYP(2, *)

      EXTERNAL calCrv
      Integer  iFnc(*)
      Real*8   parCrv(2, *)

      Integer  IHolP(*), IHolF(*), IHolE(*)

      Integer ICP(*), status

      Integer nQItr
      Logical flagMOVE, flagSWAP, flagCLPS1, flagCLPS2, flagINSRT

C group (Local variables)
      Integer ip(4), iW(10000)
      Real*8  XYPs(2), err
      Real*8  sqrEdge, tri_area2

      Logical ifXnode, check22
      Character*5 text(3)

C ======================================================================
      ip(1) = 1
      ip(2) = 2
      ip(3) = 3
      ip(4) = 1

c ... check color of mesh points
      nPt = nP + IHolP(1)
  
      Do n = 1, nPt
         iW(n) = 0
      End do 

      Do n = 1, IHolP(1)
         iPt = IHolP(n + 1)
         iW(iPt) = -1
      End do

      Do 100 n = 1, nPt
         If(iW(n).LT.0) Goto 100

         If(ICP(n).LE.0) 
     &      Call errMesMBA(5001, 'DEBUG', 'wrong point color')

         If(ifXnode(ICP(n), jInode) .AND. 
     &      ifXnode(ICP(n), jBnode)) 
     &      Call errMesMBA(5015, 'DEBUG', 'color contradiction')
 100  Continue


c ... check labels of mesh faces
      nFt = nF + IHolF(1)
  
      Do n = 1, nFt
         iW(n) = 0
      End do 

      Do n = 1, IHolF(1)
         iFt = IHolF(n + 1)
         iW(iFt) = -1
      End do

      Do 200 n = 1, nFt
         If(iW(n).LT.0) Goto 200

         If(IPF(1, n).LE.0) 
     &      Call errMesMBA(5002, 'DEBUG', 'wrong map edge -> points')

         If(lbC(n).GT.0) Then
            iCrv = lbC(n)

            Do i = 1, 2
               iPt = IPF(i, n)
               Call aniCrv(parCrv(i, iCrv), XYPs, iFnc(iCRV), calCrv)
               err = sqrEdge(XYP(1, iPt), XYPs)

               If(err.GE.1D-16)
     &            Call errMesMBA(5006, 'DEBUG', 'wrong parametrization')
            End do
         End if
 200  Continue


c ... check elements
      nEt = nE + IHolE(1)
  
      Do n = 1, nEt
         iW(n) = 0
      End do 

      Do n = 1, IHolE(1)
         iEt = IHolE(n + 1)
         iW(iEt) = -1
      End do

      Do 300 n = 1, nEt
         If(iW(n).LT.0) Goto 300

         nTet = n

         Do i = 1, 3
            If(IPE(1, n).LE.0) Call errMesMBA(5004, 'DEBUG', 
     &                              'wrong map element -> points')
         End do

         iP1 = IPE(1, n)
         iP2 = IPE(2, n)
         iP3 = IPE(3, n)
         If(iP1.EQ.iP2 .OR. iP1.EQ.iP3 .OR. iP2.EQ.iP3) Then
            iERR = 5006
            Goto 500
         End if

         If(lbE(n).LE.0) Call errMesMBA(5016, 'DEBUG', 
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


         Do 220 i1 = 1, 3
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
                     If(IFE(j1, iE).EQ.iF) Goto 210
                  End do

                  iERR = 5011
                  Goto 500
               End if

  210          Do j1 = 1, 3
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

                     v1 = tri_area2(XYP(1,iP1), XYP(1,iP2), XYP(1,iP3))
                     v2 = tri_area2(XYP(1,iP1), XYP(1,iP2), XYP(1,jP3))

c  ...   check for inverted elements
                     If(v1 * v2.GE.0D0 .AND.
     &                  .NOT.ifXnode(status, ANIUntangleMesh)) Then
                        iERR = 5013
                        Goto 500
                     End if
                     Goto 220
                  End if
               End do

               iERR = 5014
               Goto 500
            End if
 220     Continue
 300  Continue

      Return


 500  Continue
      Write(*,'(/,A,I5,/,5(A,L1),/)') 'Iteration # ', nQItr, 
     &  'Flags: moveP=', flagMOVE,  ' clpsF1=', flagCLPS1,
     &       ' clpsF2=', flagCLPS2, ' insrtP=', flagINSRT,
     &        ' swapE=', flagSWAP

      Do i = 1, 3
         Call color2text(ICP(IPE(i, nTet)), text(i))
      End do

      Write(*, 5000) nTet, iERR, (IPE(i, nTet), i = 1, 3),
     &                           (IFE(i, nTet), i = 1, 3),
     &                           (IEE(i, nTet), i = 1, 3),
     &                           (text(i),      i = 1, 3)

      Do k = 1, 3
         iE = IEE(k, nTet)
         If(iE.GT.0) Then
            Do i = 1, 3
               Call color2text(ICP(IPE(i, nTet)), text(i))
            End do

            Write(*, 5000) iE, iERR, (IPE(i, iE), i = 1, 3),
     &                               (IFE(i, iE), i = 1, 3),
     &                               (IEE(i, iE), i = 1, 3),
     &                               (text(i),    i = 1, 3)
         End if
      End do
      Do k = 1, 3
         iF = IFE(k, nTet)
         If(iF.GT.0) Then
            Write(*, 5004) iF, nF, (IPF(i, iF), i = 1, 2),
     &                             (ICP(IPF(i, max(1, iF))), i = 1, 2)
         End if
      End do

      Call errMesMBA(iERR, 'DEBUG', 'triangles are wrong')

 5000 Format('Error in checking triangle =', I6, '   iERR=', I5, /,
     &       'Points    =', 3I6, /,
     &       'Edges     =', 3I6, /,
     &       'Triangles =', 3I6, /,
     &       'colors    =', 3A6, /)

 5004 Format('Edge ', I6, '  (', I5, ')  of the bad triangle', /,
     &       'Points =', 2I6, /,
     &       'colors =', 2I6, /)

      End

    

C ======================================================================
      Subroutine color2text(iclr, text)
C ======================================================================
      include 'colors.fd'

      Integer     iclr
      Character*5 text
 
      Logical     ifXnode

C ======================================================================
      n = 5
      text = '     '

      If(ifXnode(iclr, jTnode)) Then  
         text(n:n) = 'T'
         n = n - 1
      End if

      If(ifXnode(iclr, jInode)) Then  
         text(n:n) = 'I'
         n = n - 1
      End if

      If(ifXnode(iclr, jSnode)) Then  
         text(n:n) = 'S'
         n = n - 1
      End if

      If(ifXnode(iclr, jBnode)) Then  
         text(n:n) = 'B'
         n = n - 1
      End if

      If(ifXnode(iclr, jVnode)) Then  
         text(n:n) = 'V'
         n = n - 1
      End if

      Return
      End

 

C ======================================================================
      Subroutine saveMerr(
C ======================================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE, nPv,
     &      XYP, IPF, IPE, IPV,
     &      CrvFunction, ParCrv, iFnc,
     &      nEStar, 
c group (D)
     &      nFv, nEv, IFV, IEV, lbE,
     &      flagAuto, status,
c group (Q)
     &      MaxSkipE, MaxQItr,
     &      Metric, Quality, rQuality, 
c group (W)
     &      MaxWr, MaxWi, iPrint)
C ======================================================================
C group (M)
      Integer MaxP, MaxF, MaxE
      Real*8  XYP(2, *)
      Integer IPE(3, *), IPF(4, *), IPV(*)

      EXTERNAL CrvFunction
      Real*8   ParCrv(2, *)
      Integer  iFnc(*)

C group (D)
      Integer nFv, nEv
      Integer IFV(*), IEV(*), lbE(*)
      
      Logical flagAuto
      Integer status

C group (Q)
      Integer MaxSkipE, MaxQItr
      Real*8  Metric(3, *)
      Real*8  Quality, rQuality
      Logical flagAnalytic

c group (W)
      Integer MaxWr, MaxWi, iPrint

C ======================================================================
      Open(10, file='error.ani')
      Write(10,*) nP, MaxP, nF, MaxF, nE, MaxE, nPv
      Write(10,*) ((XYP(i,n),i=1,2),n=1,nP)
      Write(10,*) ((IPF(i,n),i=1,4),n=1,nF)
      Write(10,*) ((IPE(i,n),i=1,3),n=1,nE)
      Write(10,*) (IPV(n),n=1,nPv)
      Write(10,*) ((ParCrv(i,n),i=1,2),n=1,nF)
      Write(10,*) (iFnc(n),n=1,nF)
      Write(10,*) nEStar
c 
      Write(10,*) nFv, nEv
      Write(10,*) (IFV(i), i=1, nFv)
      Write(10,*) (IEV(i), i=1, nEv)
      Write(10,*) (lbE(i), i=1, nE)
      Write(10,*) flagAuto, status
c 
      Write(10,*) MaxSkipE, MaxQItr
      Write(10,*) ((Metric(i,n),i=1,3),n=1,nE)
      Write(10,*) Quality, rQuality, flagAnalytic
      Write(10,*) MaxWr, MaxWi, iPrint
      Close(10)

      Return 
      End

