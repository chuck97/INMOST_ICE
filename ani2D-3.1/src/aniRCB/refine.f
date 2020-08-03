C ======================================================================
      Subroutine LocalRefine(
     &        nP, nPmax, nF, nFmax, nE, nEmax,
     &           XYP, IPE, IPF, lbF, lbE,
     &        RefineRule, ilevel, maxlevel, history,
     &           MaxWi, iW, iERR)
C ======================================================================
      implicit none

c ... user procedure
      EXTERNAL  RefineRule

c ... standard mesh arrays 
c ... number of points, triangles and boundary edges
      Integer  nP, nPmax, nE, nEmax, nF, nFmax

c ... maxlevel - maximum number of level bisection
      Integer maxlevel

c ... coordinates of mesh points 
      Double precision   XYP(2,nPmax)

c ... connectivity table for triangles and triangle labels
      Integer  IPE(3,nEmax), lbE(nEmax)

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
      Integer   ip_IEF, ip_verf, ip_fconfE
      Integer   ip_liW, MaxWil
C ================================================================    


c ... ip_IEF, ip_verf, ip_pbisE, ip_fconfE
      iP_pbisE   = 1
      iP_level   = iP_pbisE  + nEmax
      iP_zyk     = iP_level  + nEmax
      ip_IEF     = iP_zyk    + 6
      ip_verf    = ip_IEF    + 3*nEmax
      ip_fconfE  = ip_verf   + 3*nEmax
      ip_liW     = ip_fconfE + 3*nEmax
      
      If (MaxWi .le. 11*nEmax+7) then
          iERR = 1010
          call errMesRCB(iERR, 'LocalRefine','MaxWi  is too small')
      End if
      MaxWil = MaxWi - 11*nEmax-7

c ... define ip_IEF, ip_pbisE
      call OrientandMakeNeighbors(
     &           nP, nF, nE,
     &           XYP, IPE, iW(ip_IEF), IPF, lbF, iW(ip_pbisE),
     &           iW(ip_zyk), MaxWil, iW(ip_liW), iERR)
     
c ... define ip_verf
      call RefineRule(nE, IPE, XYP, iW(ip_verf), ilevel)

c ... bisection
      call Tria(nP, nPmax, nE, nEmax, IPE, iW(ip_IEF), XYP, 
     &      iW(ip_pbisE),  lbE, iW(ip_verf), iW(ip_fconfE), 
     &      iW(ip_level), maxlevel, history, iW(ip_zyk), iERR)

c ... define IPF
      call RestoreBndEdges(nF, nFmax, nE, IPE, iW(ip_IEF),
     &                     IPF, lbF, iW(ip_zyk), iERR)
                  
      Return
      End



C  =====================================================================
c  TRIA   produces (local) refinement of  a given set of triangles 
c  the algorithm is a modified version of that due to Baensch
C  =====================================================================
       Subroutine Tria (
     &            nP, nPmax, nE, nEmax,
     &            IPE, IEF, XYP, 
     &            pbisE,  lbE, verf, fconfE,
     &            level, maxlevel, history, zyk, iERR)
C  =====================================================================
c
c                          3
c                         /\
c                        /  \
c                       /    \
c                      /      \
c              l = 2  /        \  l = 1
c                    /          \
c                   /            \
c                  /              \
c                 /                \
c               1 ------------------  2
c                        l = 3
c
c
c  XYP(l, i)    : l-th coordinate of the i-th node
c  verf(k)      : flag whether to bisect triangle k or not
c  IPE(l, k)    : No. of l-th vertex of triangle k (anticlockwise)
c  IEF(l, k)    : neighbour of triangle k opposite to vertex l
c  fconfE(l, k) : flag whether there's a "conflict" on edge l
c  pbisE(k)     : No. of the refinement edge of triangle k
C  =====================================================================
      implicit none

      Double precision XYP(2,*)
      Integer          nP, nPmax, nE, nEmax
      Integer          maxlevel
      Integer          IPE(3,*), IEF(3,*)
      Integer          pbisE(*)
      Integer          lbE(*)
      Integer          verf(*)
      Integer          fconfE(3,*)
      Integer          level(*)
      Integer          iERR
      Logical          history(maxlevel,*)
                      
      Integer          zyk(6)
      Integer          i, j
C  ================================================================ 
c ... initialization:
      Do i = 1, nE
        Do j = 1, 3
          fconfE(j, i) = 0
        End do
      End do

c ... bisect elemenes until there is  no conflict:
1     Continue

      Do i = 1, nE
        If (verf(i) .ne. 0)  
     &    call Teil(i, nP, nPmax, nE, nEmax, IPE, IEF, XYP,
     &           zyk, pbisE,  lbE, verf, fconfE,
     &           level, maxlevel, history, iERR)
      End do
         
      Do i = 1, nE
        If (verf(i) .ne. 0)  go to 1
      End do
c
      Return
      End



C  ================================================================ 
c  TEIL divides the triangle i into two new ones
c
c  INPUT:      i - no. of triangle 
C  ================================================================ 
      Subroutine Teil(
     &           i, nP, nPmax, nE, nEmax, IPE, IEF, XYP, 
     &           zyk, pbisE,  lbE, verf, fconfE,
     &           level, maxlevel, history, iERR)
C  ================================================================ 
      implicit none

      Integer           i
      Integer           nP, nPmax, nE, nEmax
      Integer           maxlevel
      Integer           IPE(3,*), IEF(3,*)
      Double precision  XYP(2,*) 
      Integer           zyk(6)
      Integer           pbisE(*)
      Integer           lbE(*)
      Integer           verf(*)
      Integer           fconfE(3,*)
      Integer           level(*)
      Integer           iERR
      Logical           history(maxlevel, *)
      
      Integer           ii, l, j
      logical           flag
C  ================================================================ 
c
c                 /\                                 |
c                /  \                               /| \
c               /    \                             / |  \
c              /      \                           /  |   \
c             /        \     ----->     pbisE(i) /   |    \  pbisE(ii)
c            /          \                       /    |     \
c           /     i      \                     /     |      \
c          /              \                   /   i  |  ii   \
c         /                \                 /       |        \
c         ------------------                 ------------------
c             l = pbisE(i)
c
C  ================================================================ 

c ... l - refinement edge
      l = pbisE(i)
      If (fconfE(l, i) .eq. 0) Then
         flag = .false.
       else
         flag = .true.
      End if   
      
c ... one more triangle:  ii - No. of new triangle
      If (nE .ge. nEmax) Then
        iERR = 1006
        call errMesRCB(iERR, 'teil','nEmax is too small')
      End if  

c ... introduce new triangle
      nE = nE+1
      verf(nE) = max(verf(i)-1,0)
      verf(i)  = verf(nE)
      lbE(nE) = lbE(i)
      ii = nE

c ... define level 
      level(i)  = level(i) + 1
      level(ii) = level(i) 
      
c ... define history
      history(level(i), i)   = .false.
      history(level(ii), ii) = .true.
      
c ... define the value of  some variables for the
c     new triangle and change values for the old one
      fconfE(zyk(l+1), ii) = fconfE(zyk(l+1), i)
      fconfE(zyk(l+1), i)  = 0
      fconfE(l, i)         = 0
      fconfE(l, ii)        = 0
      fconfE(zyk(l+2), ii) = 0

      IEF(zyk(l+1), ii) = IEF(zyk(l+1), i)
      IEF(zyk(l+1), i)  = ii
      IEF(zyk(l+2), ii) = i

      pbisE(i)  = zyk(l+2)
      pbisE(ii) = zyk(l+1)

      IPE(zyk(l+2), ii) = IPE(zyk(l+2), i)
      IPE(l, ii)        = IPE(l, i)

c ...is there a "conflict" for i or ii:
      verf(i)  = max (verf(i),  fconfE(zyk(l+2), i))
      verf(ii) = max (verf(ii), fconfE(zyk(l+1), ii))
      
      j = IEF(l, i)
      If (j .gt. 0) Then
       
c ...adjust all variables in the elements meeting
c    at the refinement edge:
        Call bisectpoint(i, j, l, ii, flag, 
     &       XYP, nP, nPmax, zyk, IPE, IEF, fconfE, verf, iERR)
     
       else
c ...refinement edge is a part of the boundary:
        call addpoint (i, j, l, ii, XYP, nP, nPmax,
     &                 zyk, IPE, IEF, iERR)
           
      End if

      j = IEF(zyk(l+1), ii)
      If (fconfE(zyk(l+1), ii) .eq. 1) Then
        flag = .true.
       else
        flag = .false.
      End if  
c ...adjust all variables in the "right" triangle
c     If l+1 is not a boundary edge:
      If (j .gt. 0) Then
        call sidepoint (i, j, ii, flag, zyk, IEF)
      End if
c
      Return
      End



C  ================================================================ 
c  BISECTPOINT  deals with the triangle meeting
c        at the refinement edge:
c
c  i    : no. of the triangle
c  j    : no. of the neighbor
c  flag : = true if fconfE(l, i) = 1
C  =====================================================================
      Subroutine bisectpoint(
     &           i, j, l, ii, flag, XYP, 
     &           nP, nPmax, zyk, IPE, IEF, fconfE, verf, iERR)
C  =====================================================================
      implicit none

      Double precision  XYP(2,*)
      Integer           nP, nPmax, zyk(6), IPE(3,*), IEF(3,*)
      Integer           fconfE(3,*)
      Integer           verf(*)
      Integer           iERR
     
      Integer           i, j, l, ii, ll
      Integer           k, kk, k2, j2, lauf
      Integer           j3, jj
      logical           flag
C  ================================================================ 

       k = 0
       Do kk = 1, 3
         If (IEF(kk, j) .eq. i) then
           k = kk
           go to 1
         End if  
       End do
1      Continue       

       If(flag) Then
       
c ... flag  =  .true. : j was divided before
         IPE(zyk(l+2), i)  = IPE(zyk(k+1), j)
         IPE(zyk(l+1), ii) = IPE(zyk(k+1), j)

         j2 = j
         k2 = k
         Do lauf = 1, 3
             j3 = IEF(zyk(k2+2), j2)
             Do ll = 1, 3
                 If(IEF(ll, j3).eq.j2) goto 4
             End do  
4            k2 = ll
             j2 = j3
             If(IEF(zyk(k2+2), j2).eq.i) Then
                 jj = j2
                 kk = zyk(k2+2)
                 goto 5
             End if
         End do  

5        IEF(kk, jj) = ii
         IEF(l,  ii) = jj

         verf(j)  = max (verf(j),  fconfE(zyk(k+1),  j))
         verf(jj) = max (verf(jj), fconfE(zyk(kk+2), jj))
       else
       
c  ...flag  =  .false.  : j hasn't been bisected so far we have to create a new node
         verf(j) = max (verf(j), 1) 
         fconfE(k, j) = 1

         IEF(l, ii) = j
         IEF(k, j)  = ii
c add a new node:
         If (nP.ge.nPmax) Then
            iERR = 1003
            Call errMesRCB(iERR, 'teil','nPmax is too small')
         End if   

         nP = nP+1
         IPE(zyk(l+2), i)  = nP
         IPE(zyk(l+1), ii) = nP

         XYP(1, nP) = ( XYP(1, IPE(zyk(l+1), i))
     *                + XYP(1, IPE(zyk(l+2), ii)) )/2.d0
         XYP(2, nP) = ( XYP(2, IPE(zyk(l+1), i))
     *                + XYP(2, IPE(zyk(l+2), ii)) )/2.d0
       End if

       Return
       End



C  =====================================================================
c  SIDEPOINT  gives the information about triangle ii to the "right"
c            neighbor fconfE(l+1, i)
C  ===================================================================== 
      subroutine sidepoint(i, j, ii, flag, zyk, IEF)
C  =====================================================================
      implicit none

      Integer           zyk(6), IEF(3,*)
      Integer           i, j, ii
      Integer           k, kk, k2, jj, j2, j3, lauf, ll
      logical           flag
C  =====================================================================

      k = 0
      Do kk = 1, 3
        If (IEF(kk, j) .eq. i) then
          k = kk
          go to 1
        End if  
      End do  
1     Continue

      If (k .ne. 0) Then
       
        IEF(k, j) = ii
       
        If (flag) Then
            j2 = j
            k2 = k
            Do lauf = 1, 3
              j3 = IEF(zyk(k2+2), j2)
              
              Do ll = 1, 3
                If (IEF(ll, j3).eq.j2) then
                  k2 = ll
                  go to 2
                End if  
              End do
2             Continue              

              j2 = j3
              If (IEF(zyk(k2+2), j2).eq.i) Then
                jj = j2
                kk = zyk(k2+2)
                goto 3
              End if
              
            End do
3           IEF(kk, jj) = ii
         End if
         
      End if
        
      Return
      End



C  ===================================================================== 
c   ADDPOINT the refinement edge is now a boundary edge
C  ===================================================================== 
      subroutine addpoint (
     &       i, j, l, ii, XYP, nP, nPmax, zyk, IPE, IEF, iERR)
C  =====================================================================
      implicit none
c
      Double precision  XYP(2,*)
      Integer           nP, nPmax, zyk(6), IPE(3,*), IEF(3,*)
      Integer           iERR
      
      Integer           i, j, l, ii
C  ===================================================================== 
c ...add a new node:
       If (nP.ge.nPmax) Then
          iERR = 1003
          call errMesRCB(iERR, 'teil','nPmax is too small')
       End if  

       nP = nP+1
       IPE(zyk(l+2), i)  = nP
       IPE(zyk(l+1), ii) = nP

c ...in case of a curved boundary here's the right place to put in the information:
       XYP(1, nP) = ( XYP(1, IPE(zyk(l+1), i))
     *              + XYP(1, IPE(zyk(l+2), ii)) )/2.d0
       XYP(2, nP) = ( XYP(2, IPE(zyk(l+1), i))
     *              + XYP(2, IPE(zyk(l+2), ii)) )/2.d0

       IEF(l, i) = j
       IEF(l, ii) = j

       Return
       End

