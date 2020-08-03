C ======================================================================
      subroutine LocalCoarse(
     &        nP, nPmax, nF, nFmax, nE, nEmax,
     &           XYP, IPE, IPF, lbF, 
     &           CoarseRule, ilevel, maxlevel, history,
     &           MaxWi, iW, iERR)
C ======================================================================
      implicit none

c ... user procedure
      external  CoarseRule

c ... standard mesh arrays 
c ... number of points, triangles and boundary edges
      Integer  nP, nPmax, nE, nEmax, nF, nFmax

c ... maxlevel - maximum number of level bisection
      Integer maxlevel

c ... coordinates of mesh points 
      Double precision   XYP(2,nPmax)

c ... connectivity table for triangles and triangle labels
      Integer  IPE(3,nEmax)

c ... connectivity table for boundary edges, and edge labels
      Integer  IPF(2,nFmax), lbF(nFmax)

c ... level of bisection
      Integer  ilevel

c ... history of bisection
      Logical history(maxlevel, nPmax)
      
c ... work memory
      Integer   MaxWi, iERR
      Integer   iW(MaxWi)
      
c ... local point
      Integer   iP_pbisE, iP_level, iP_zyk
      Integer   ip_IEF, ip_verf, ip_pconfE
      Integer   ip_realP, ip_realE
      Integer   ip_liW, MaxWil

C ======================================================================
c ... ip_IEF, ip_verf, ip_pbisE, ip_pconfE
      iP_pbisE   =                 1
      iP_level   = iP_pbisE        + nEmax
      iP_zyk     = iP_level        + nEmax
      ip_IEF     = iP_zyk          + 6
      ip_verf    = ip_IEF          + 3*nEmax
      ip_pconfE  = ip_verf         + nEmax
      ip_realE   = ip_pconfE       + 3*nEmax
      ip_realP   = ip_realE        + nEmax
      ip_liW     = ip_realP        + nPmax
      
      If(MaxWi .le. 11*nEmax + nPmax + 7) then
         iERR = 1010
         Call errMesRCB(iERR, 'LocalCoarse','MaxWi  is too small')
      End if
      MaxWil = MaxWi - 11*nEmax - nPmax - 7

c ... define ip_IEF, ip_pbisE
      call OrientandMakeNeighbors(
     &           nP, nF, nE,
     &           XYP, IPE, iW(ip_IEF), IPF, lbF, iW(ip_pbisE),
     &           iW(ip_zyk), MaxWil, iW(ip_liW), iERR)
     
c ... define ip_verf
      call CoarseRule(nE, IPE, XYP, iW(ip_verf), ilevel)

c ... bisection
      call Untria(nP, nE, IPE, iW(ip_IEF), XYP, 
     &            iW(ip_pbisE), iW(ip_verf), iW(ip_pconfE), 
     &      iW(ip_level), maxlevel, history, 
     &            iW(ip_realE), iW(ip_realP), iW(ip_zyk))

c ... define IPF
      call RestoreBndEdges(nF, nFmax, nE, IPE, iW(ip_IEF),
     &                     IPF, lbF, iW(ip_zyk), iERR)
                  
      Return
      End



C  =====================================================================
c  UNTRIA   produces (local) coarsing of  a given set of triangles 
c  the algorithm is a modified version of that due to Baensch
C  =====================================================================
       Subroutine Untria(
     &            nP, nE, 
     &           IPE, IEF, XYP, 
     &            pbisE, verf, pconfE,
     &           level, maxlevel, history,
     &            realE, realP, zyk)
C  =====================================================================
c
c                          3
c                         /\
c                        /  \
c                       /    \
c                      /      \
c                     /        \       
c                    /          \
c                   /            \
c                  /              \
c                 /                \
c               1 ------------------  2
c
c
c
c  XYP(l, i)    : l-th coordinate of the i-th node
c  verf(k)      : flag whether to bisect triangle k or not
c  IPE(l, k)    : No. of l-th vertex of triangle k (anticlockwise)
c  IEF(l, k)    : neighbour of triangle k opposite to vertex l
c  pconfE(l, k) : flag whether there's a "conflict" on edge l
c  pbisE(k)     : No. of the refinement edge of triangle k
C  ================================================================ 
      implicit none

      Double precision  XYP(2,*)
      Integer           nP, nE
      Integer           IPE(3,*), IEF(3,*)
      Integer           pbisE(*), verf(*), pconfE(3,*)
      Integer           level(*), maxlevel
      Integer           realE(*), realP(*)
      Logical           history(maxlevel, *)
                      
      Integer           zyk(6)
      Integer           i, j, k, maxlevelmesh

C  =====================================================================
c ... initialization:

      Do i = 1, nE
        realE(i) = 1
      End do  
      
      Do i = 1, nP
        realP(i) = 1
      End do  

      Do i = 1, nE
        Do j = 1, 3
          pconfE(j, i) = 0
        End do
      End do

      Do i = 1, nE
       IF (verf(i) .ge. 1) goto 1
      End do

      Return

c ... coarse elements until there is  no conflict:
1     Continue

      maxlevelmesh=0
      Do i = 1, nE
       IF ((realE(i) .eq. 1) .and.  
     &      (verf(i) .ge. 1) .and.
     &      (level(i) .gt. maxlevelmesh))
     &   maxlevelmesh = level(i)
      End do

      Do i = 1, nE
        If((verf(i).ge.1) .and. (realE(i).eq.1) 
     &                    .and. (level(i).eq.maxlevelmesh)) 
     &    Call Unteil(i, IPE, IEF, 
     &                zyk, pbisE, verf, pconfE,
     &                level, maxlevel, history, realE, realP)
      End do
         
      Do i = 1, nE
        If ((verf(i) .ge. 1) .and. (realE(i) .eq. 1))  go to 1
      End do

c ... delete point
      k = 0
      Do i = 1, nP 
        If (realP(i) .eq. 1) then
          k = k + 1
          realP(i) = k
          Do j = 1, 2
            XYP(j,k) = XYP(j,i)
          End do 
        End if
      End do
      nP = k

c ... delete triangles
      k = 0
      do i = 1, nE
        if (realE(i) .eq. 1) then
          k = k + 1
          realE(i) = k
        end if
      end do

      Do i = 1, nE
       If (realE(i) .ge. 1) then
         Do j = 1, 3
           IPE(j, realE(i)) = realP(IPE(j, i))
           if (IEF(j, i).ge.1) then
             IEF(j,realE(i)) = realE(IEF(j, i))
            else
             IEF(j,realE(i)) = IEF(j,i)
           end if  
         end do
         pbisE(realE(i)) = pbisE(i)
         level(realE(i)) = level(i)
         Do j=1, level(i)
           history(j, realE(i)) = history(j, i)
         End do  
       End if
      End do 
      nE = k

      Return
      End



C  =====================================================================
c  UNTEIL coarse the triangle i from i and neightbor
c
c  INPUT:      i - no. of triangle 
C  =====================================================================
      Subroutine Unteil(
     &           i, IPE, IEF, 
     &           zyk, pbisE, verf, pconfE,
     &           level, maxlevel, history, realE, realP)
C  =====================================================================
      implicit none

      Integer           i
      Integer           IPE(3,*), IEF(3,*)
      Integer           zyk(6)
      Integer  pbisE(*), verf(*), pconfE(3,*)
      Integer  level(*), maxlevel
      Integer           realE(*), realP(*)
      Logical           history(maxlevel, *)
     
      Integer           pcorE(4,2)
      Integer           ii, j, jj,  k, kk, v
      Integer           deltr, keeptr, delin, keepin
      Integer           jcikl, lcikl, ncikl 
      Logical  flag

C  =====================================================================
c
c                   /|\                                   /\          
c                  / | \                                 /  \        
c                 /  |  \                               /    \       
c                /   |   \                             /      \      
c      pbisE(i) /    |    \  pbisE(ii)   ----->       /        \     
c              /     |     \                         /          \    
c             /      |      \                       /     i      \   
c            /   i   |  ii   \                     /              \  
c           /        |        \                   /                \ 
c           -------------------                   ------------------ 
c                                                     pbisE(i)
c
C  =====================================================================
      If (level(i) .eq. 0) Then
         verf(i) = 0
         Return
      End if   

c ... define pcorE for triangle i
      pcorE(1, 1) = i
      pcorE(2, 1) = pbisE(i)
      If (history(level(i),i)) then
        kk = 2
       Else
        kk = 1
      End if
      pcorE(3, 1) = zyk(pbisE(i)+kk)
      Do kk = 1, 3
        If ((kk .ne. pcorE(2, 1)) .and. (kk .ne. pcorE(3, 1))) then
           pcorE(4, 1) = kk
           go to 1
        End if
      End do
1     Continue      

c ... define triangle ii
      ii = IEF (pcorE(4, 1),i)
      
c ... define pcorE for triangle i 
      pcorE(1, 2) = ii
      Do kk = 1, 3
        If (IPE(kk, ii) .eq. IPE(pcorE(2, 1), i)) pcorE(2, 2) = kk    
        If (IPE(kk, ii) .eq. IPE(pcorE(3, 1), i)) pcorE(3, 2) = kk    
      End do  
      Do kk = 1, 3
        If ((kk .ne. pcorE(2, 2)) .and. (kk .ne. pcorE(3, 2))) then
           pcorE(4, 2) = kk
           go to 2
        End if
      End do
2     Continue      
  
      verf(ii) = max(verf(ii), 1)

c ... define conflict nodes
      v = IPE(pcorE(2, 1), i)
      flag = .false.

c ... control of correspondence of level, pbisE and history
c     of triangles i and ii 
      If (level(i) .ne. level(ii)) then
         If (pconfE(pcorE(2, 1), i) .eq. 1) Return
         goto 10  
      End if
    
c ... define flag
      If (pconfE(pbisE(i),i) .eq. 0) then
        realP(v) = 0
        flag = .true.
      End if
      
c ... define delete triangle
      IF (ii .gt. i) then
        deltr  = ii
        keeptr = i
        delin  = 2
        keepin = 1
       Else
        deltr  = i
        keeptr = ii
        delin  = 1
        keepin = 2
      End if

      
c ... delete triangle deltr
      realE(deltr) = 0     

c ... define IPE
      IPE(pcorE(2, keepin), keeptr) = IPE(pcorE(4, delin), deltr)

c ... define IEF
      IEF(pcorE(4, keepin), keeptr) = IEF(pcorE(2, delin), deltr)
       
      j = IEF(pcorE(2, delin), deltr)
      If (j .gt. 0) then
        Do kk = 1, 3
          If (IEF(kk, j) .eq. deltr) then
            IEF(kk, j) = keeptr
            go to 3
          End if  
        End do
3       Continue
      End if 
      
 
      j = IEF(pcorE(3, delin), deltr)
      If (j .gt. 0) then
        Do kk = 1, 3
          If (IEF(kk, j) .eq. deltr) then
            IEF(kk, j) = keeptr
            go to 4
          End if
        End do
4       Continue
      End if 

c ... define pbisE
      pbisE(keeptr) = pcorE(3, keepin)

c ... define pconfE
      pconfE(pcorE(2, keepin), keeptr) = pconfE(pcorE(4, delin), deltr)

c ... define verf
      verf(keeptr) = max(max(verf(keeptr),verf(deltr)) - 1,
     &                   pconfE(1,keeptr), 
     &                   pconfE(2,keeptr), 
     &                   pconfE(3,keeptr)) 

c ... define level 
      level(keeptr) = level(keeptr) - 1

      if (.not. flag) Return

10    continue
     
      Do k = 1, 2
        jcikl = pcorE(1, k)
        lcikl = pcorE(3, k)
        ncikl = pcorE(4, k)

20      continue
        jj = IEF(lcikl, jcikl)
        If (jj .le. 0) go to 30
        If ((jj .eq. i) .or. (jj .eq. ii)) go to 40
        If (verf(jj) .lt. 1) verf(jj) = 1
        If (flag) Then
          Do kk = 1, 3
            If ((IPE(kk, jj) .eq. v) .and. 
     &          (pconfE(kk,jj) .ne. 1)) then
                    pconfE(kk,jj) = 1
                    go to 21
            End if        
          End do
21        Continue          
        End if
        
        Do kk = 1, 3
          If ((IPE(kk, jj) .eq. IPE(ncikl, jcikl))) then
            lcikl = kk    
            go to 22
          End if  
        End do
22      Continue              

        Do kk = 1, 3
          If ((IPE(kk, jj) .ne. v) .and. (kk. ne. lcikl)) then
            ncikl = kk         
            goto 23
          End if  
        End do
23      Continue

        jcikl = jj
        go to 20
        
30      continue
      End do

40    continue      
        
      Return
      End
        
