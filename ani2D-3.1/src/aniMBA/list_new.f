C ======================================================================
C The routines below work with a weakly ordered list of real numbers.
C The available operations are:
C     (a) intialize the list
C     (b) add new element in the list
C     (c) update an element value
C     (d) exclude an element from the list
C
C The dynamic order is given by array L2E(2, MaxE). 
C The minimal interval length is max(AniMinInterval, iCntl(1) / 4). 
C An interval is removed if it is smaller than that. 
C The maximal interval length is 2 iCntl(1). 
C An interval is split into two is it is biggler than that.
C
C iCntl(1) is typical interval length, we recommend 7*log(nEStar).
C
C ======================================================================
      Subroutine lstMak(nE, L1E, L2E, nL2, iCntl, LHol)
C ======================================================================
      Implicit none
      Include 'list.fd'

      Integer L1E(2,*), L2E(2,*), LHol(*), nE, nL2, iCntl(4)

      Integer i, i1, i2, n, iE1, iE2, nstep

C ======================================================================
      Call weirdCopy(nE, L2E, L2E)

      iE1 = L2E(1,1)
      iE2 = L2E(1,nE)

      L1E(1,iE1) = 0  ! pointers for the 1st and last elements are 0 
      L1E(2,iE2) = 0

      i2 = iE1
      Do n = 2, nE
         i1 = i2
         i2 = L2E(1,n)
         L1E(2,i1) = i2
         L1E(1,i2) = i1
      End do

c ... mark the array with equally distributed pointers
      nstep = min(iCntl(1), AniMinInterval)

      nL2 = 0
      Do n = 1, nE, nstep
         nL2 = nL2 + 1
         L2E(1,nL2) = L2E(1,n)
         L2E(2,nL2) = nstep
      End do
      L2E(2,nL2) = nE - nstep * (nL2-1)

      LHol(1) = 0

      Return
      End



C ======================================================================
      Subroutine lstUpd(nE, L1E, nL2, L2E, iCntl, qE, iE, qiE)
C ======================================================================
C Updates element iE in the list {qE, L2E,L2E}. We ellimate the
C corresponding interval if it becomes small and split the inteval into
C two if it becomes too big.
C ======================================================================
      Implicit none
      Integer L1E(2,*), L2E(2,*), nE, nL2, iE, iCntl(4)
      Real*8  qE(*), qiE

      Integer i1, i2, LHol(2)
      Logical flag

C ======================================================================
      i1 = L1E(1,iE)
      i2 = L1E(2,iE)
      
      flag = .TRUE.
      If(i1.GT.0) flag = qE(i1).LE.qiE

      If(flag .AND. i2.GT.0) flag = qE(i2).GE.qiE

      If(.NOT.flag) Then 
         LHol(1) = 0
         Call lstDel(nE, L1E, nL2, L2E, iCntl, LHol, qE, iE)
         Call lstAdd(nE, L1E, nL2, L2E, iCntl, LHol, qE, qiE, iE)
      Else
         qE(iE) = qiE
      End if

      Return  
      End



C ======================================================================
      Subroutine lstDel(nE, L1E, nL2, L2E, iCntl, LHol, qE, iE)
C ======================================================================
C Deletes element iE from the list {qE, L2E,L2E}. We ellimate the
C corresponding interval if it becomes small but do not check overflow 
C of the large interval.
C ======================================================================
      Implicit none
      include 'list.fd'
C ======================================================================
      Integer L1E(2,*), L2E(2,*), LHol(*), nE, nL2, iE, iCntl(4)
      Real*8  qE(*)

      Integer i, i1,i2,i3, n,nmin, iEt,iEs, findIntervalHard

C ======================================================================
      If(L1E(1,iE) + L1E(2,iE).EQ.0) Then
         Call errMesMBA(5103, 'lstDel', 'no such element in the list')
      End if

      i1 = findIntervalHard(nL2, L1E, L2E, qE, iE)

c ... remove iE from the list
      iEs = L1E(1,iE)
      iEt = L1E(2,iE)

      If(iEt.NE.0) L1E(1,iEt) = iEs
      If(iEs.NE.0) L1E(2,iEs) = iEt
      nE = nE - 1

c ... update the ordered list 
      L2E(2,i1) = L2E(2,i1) - 1
      n = L2E(2,i1)
      nmin = max(AniMinInterval, iCntl(1) / AniMinFraction)

      If(n.LT.nmin) Then
         If(nL2.EQ.1) Return

         If(i1.LT.nL2) Then
            i2 = i1 + 1
            L2E(2,i1) = L2E(2,i1) + L2E(2,i2)
            i3 = i2
         Else
            i2 = i1 - 1
            L2E(2,i2) = L2E(2,i2) + L2E(2,i1)
            i3 = i1 
         End if

         nL2 = nL2 - 1
         Do i = i3, nL2
            L2E(1,i) = L2E(1,i+1)
            L2E(2,i) = L2E(2,i+1)
         End do
      End if

c ... check that iE is not a pointer
      If(iE.EQ.L2E(1,i1)) L2E(1,i1) = L1E(2,iE)

c ... make a new hole
      L1E(1,iE) = 0
      L1E(2,iE) = 0
      LHol(1) = LHol(1) + 1
      LHol(LHol(1)+1) = iE

      Return  
      End



C ======================================================================
      Subroutine lstAdd(nE, L1E, nL2, L2E, iCntl, LHol, qE, qiE, iE)
C ======================================================================
C Adds element iE to the list {qE, L2E,L2E}. We split the corresponding 
C interval into two if it becomes small.
C ======================================================================
      Implicit none
      include 'list.fd'

      Integer L1E(2,*), L2E(2,*), LHol(*), nE, nL2, iE, iCntl(4)
      Real*8  qE(*), qiE

      Integer i, i1,i2, m,n,nmax, iEt,iEs, findInterval 
      Logical flagFP
  
C ======================================================================
      Call fpcheck(qiE, flagFP)
      If(flagFP) Then
         Call errMesMBA(5113, 'lstAdd', 'bad input value(NAN/INF)')
      End if

      i1 = findInterval(nL2, L2E, qE, qiE)

c ... fill in a hole if any
      nE = nE + 1
      If(LHol(1).eq.0) Then
         iE = nE
      Else
         iE = LHol(LHol(1)+1)
         LHol(1) = LHol(1) - 1
      End if
      qE(iE) = qiE

c ... zero interval
      If(i1.EQ.0) Then
         iEt = L2E(1,1)
         L2E(1,1) = iE 
         i1 = 1

         L1E(1,iE)  = 0
         L1E(2,iE)  = iEt
         L1E(1,iEt) = iE

         Goto 100
      End if

c ... middle interval
      iEt = L2E(1,i1)
      n   = L2E(2,i1)

      Do i = 1, n
         iEs = L1E(2,iEt)
         If(iEs.EQ.0) Then
            L1E(1,iE) = iEt
            L1E(2,iE) = 0
            L1E(2,iEt) = iE 
            Goto 100
         Else If(qiE.LE.qE(iEs)) Then
            L1E(1,iEs) = iE
            L1E(2,iEt) = iE

            L1E(1,iE) = iEt
            L1E(2,iE) = iEs
            Goto 100
         End if
         iEt = iEs
      End do
      Write(*,*) 'oops...no slot for', qiE, iEs, L2E(1,i1+1), qE(iEs)

c ... update the ordered list
 100  Continue
      L2E(2,i1) = L2E(2,i1) + 1
      n = L2E(2,i1)
      nmax = AniMaxFraction * iCntl(1)

      If(n.GT.nmax) Then
         nL2 = nL2 + 1
         i2 = i1 + 1
         Do i = nL2, i2, -1
            L2E(1,i) = L2E(1,i-1)
            L2E(2,i) = L2E(2,i-1)
         End do

         m = n / 2
         L2E(2,i1) = m
         L2E(2,i2) = n - m 

         iEt = L2E(1,i1)
         Do i = 1, m
            iEt = L1E(2,iEt)
         End do
         L2E(1,i2) = iEt
      End if

      Return
      End



C ======================================================================
      Integer Function findInterval(nL2, L2E, qE, qiE)
C ======================================================================
C We return interval that contains qiE
C ======================================================================
      Implicit none
      Integer nL2, L2E(2,*)
      Real*8  qE(*), qiE

      Integer i, i1,i2,i3, iEt

C ======================================================================
      If(qiE.LT.qE(L2E(1,1))) Then
         findInterval = 0
         Return
      End If

      i1 = 1
      i2 = nL2
      Do While(i1.LT.i2-1)
         i3 = (i1+i2) / 2
         If(qiE.LE.qE(L2E(1,i3))) Then
            i2 = i3
         Else
            i1 = i3
         End If
      End do

      If(qiE.GT.qE(L2E(1,i2))) i1 = i2
      findInterval = i1

      Return
      End



C ======================================================================
      Integer Function findIntervalHard(nL2, L1E, L2E, qE, iE)
C ======================================================================
C We return interval that contains existing iE
C ======================================================================
      Implicit none
      Integer nL2, L1E(2,*), L2E(2,*), iE
      Real*8  qE(*)

      Integer i, i1,i2,i3, iEt
      Real*8  qiE

C ======================================================================
      qiE = qE(iE)

      i1 = 1
      i2 = nL2
      Do While(i1.LT.i2-1)
         i3 = (i1+i2) / 2
         If(qiE.LE.qE(L2E(1,i3))) Then
            i2 = i3
         Else
            i1 = i3
         End If
      End do

c ... a difficult case when all qE's are the same
      iEt = L2E(1,i1)
      Do While(i1.LT.nL2)
         i2 = i1 + 1

         findIntervalHard = i1
         If(qiE.LT.qE(L2E(1,i2))) Return 

         Do i = 1, L2E(2,i1)
            If(iEt.EQ.iE) Return
            iEt = L1E(2,iEt)
         End do
         i1 = i1 + 1
      End do
      findIntervalHard = nL2

      Return
      End



C ======================================================================
      Subroutine weirdCopy(nE, iOrd, L2E)
C ======================================================================
C Allows to overload unknowns, introduced for consistency
C ======================================================================
      Integer iOrd(*), L2E(2,*)

      Do n = nE, 1, -1
         L2E(1,n) = iOrd(n)
      End do

      Return
      End

      

C ======================================================================
      Integer Function PrevL2(nL2, L2E, qE, qiE)
C ======================================================================
C Compatability (temporary/obsolete) routine
C ======================================================================
      Implicit none
      Integer L2E(2,*), nL2, findInterval
      Real*8  qE(*), qiE

      PrevL2 = findInterval(nL2, L2E, qE, qiE)

      Return
      End



C ======================================================================
       Integer Function PrevL1(L1E, L2E, iPrevL2, qE, qiE)
C ======================================================================
C Compatability (temporary ?) routine
C ======================================================================
      Implicit none
      Integer L1E(2,*), L2E(2,*), iPrevL2
      Real*8  qE(*),qiE

      Integer j,k

C ======================================================================
      If(iPrevL2.le.0) Then
         Call errMesMBA(5005, 'PrevL1', 'wrong iPrevL2')
      End if

      j = L2E(1,iPrevL2)

      Do while (.true.)
         If (qiE.lt.qE(j)) Then
            PrevL1 = L1E(1,j)
            Return
         End If
         k = L1E(2,j)
         If (k.eq.0) Then
            PrevL1 = j
            Return
         End If
         j = k
      End do

      Return
      End

