c ======================================================================
c Initialization routine
c ======================================================================
      Subroutine InitializeRCB(nE, nEmax, XYP, IPE, MaxWi, iW, iERR)
c ======================================================================
      implicit none

      Integer             nE, nEmax, MaxWi, iERR
      Double precision    XYP(2,*)
      Integer             IPE(3,*)
      Integer             iW(*) 

      Integer             iP_pbisE, iP_level 
      Integer             iP_zyk
      
c ======================================================================
      if (MaxWi.lt.11*nEmax+6) then
          iERR=1 ! LocalRefine needs at least 11*nEmax+6
          return
      end if

      iP_pbisE   =              1
      iP_level   = iP_pbisE   + nEmax
      iP_zyk     = iP_level   + nEmax
  
      call InitializeMeshData (nE, XYP, IPE, iW(iP_pbisE),
     &                 iW(iP_level), iW(iP_zyk))

      Return
      End



c ======================================================================
      Subroutine InitializeMeshData(nE, XYP, IPE, pbisE, level, zyk)
c ======================================================================
      implicit none

      Integer             nE
      Double precision    XYP(2,*)
      Integer             IPE(3,*)
      Integer             pbisE(*)
      Integer             level(*)
      Integer             zyk(6)
      
      Integer             i, i1, i2, j, j1
      Double precision    s, s1

c ======================================================================
c
c Define  zyk 
c
      zyk(1) = 1 
      zyk(2) = 2 
      zyk(3) = 3 
      zyk(4) = 1 
      zyk(5) = 2
      zyk(6) = 3 
c
c Define pbisE
c
      Do i = 1, nE
      
        s  = 0
        j1 = 0
        
        Do j = 1, 3
        
          i1 = IPE(zyk(j + 1),i)
          i2 = IPE(zyk(j + 2),i)
          
          call Dist(XYP(1,i1),XYP(1,i2),s1)
          
          If (s1 .gt. s) then
            s  = s1
            j1 =j
          End if
          
        End do
        
        pbisE(i)=j1 ! Initial split of element by a median emanating from vertex j1.
                    ! The applied rule is to split the longest edge of triangle. 
                    ! No conformity should be cared of!
                    ! The user is free to change the rule.
        
      End do
        
c
c Define level
c
       Do i = 1, nE
         level(i)=0
       End do

      Return
      End



c ======================================================================
      Subroutine Dist(pt1,pt2,res)
c ======================================================================
      implicit none
      
      Double precision  pt1(2), pt2(2), res
c ======================================================================
       res = dsqrt((pt1(1) - pt2(1)) **2 +
     &             (pt1(2) - pt2(2)) **2)
                         
      Return
      End
      


c ======================================================================
c Addition auxiliary procedure.
c ======================================================================
      subroutine OrientandMakeNeighbors(
     &           nP, nF, nE,
     &           XYP, IPE, IEF, IPF, lbF, pbisE, 
     &           zyk, MaxWi, iW, iERR)
c ======================================================================
      implicit none

      Integer             nP, nF, nE
      Integer             IPE(3,*), IEF(3,*)
      Integer             IPF(2,*), lbF(*)
      Double precision    XYP(2,*)
      Integer             pbisE(*)
      Integer             zyk(6), MaxWi, iW(*)
      Integer             iERR

      Integer             maxmetIPF
      Integer             i,j,i1,i2,i3,i4
      Integer             k,l,m,k1,k2,m1,m2
      Double precision    s,ax,ay,bx,by

c ======================================================================
c
c  Check of the orientation of triangul 
c
       Do i = 1, nE
   
c ... forming the vector from 1 and 2 vertices of triangle  
         ax = XYP(1,IPE(2,i)) - XYP(1,IPE(1,i))
         ay = XYP(2,IPE(2,i)) - XYP(2,IPE(1,i))
         s  = dsqrt(ax**2 + ay**2)
         ax = ax / s
         ay = ay / s
   
c ... forming the vector from 2 and 3 vertices of triangle  
         bx = XYP(1,IPE(3,i)) - XYP(1,IPE(2,i))
         by = XYP(2,IPE(3,i)) - XYP(2,IPE(2,i))
         s  = dsqrt(ax**2 + ay**2)
         bx = bx / s
         by = by / s

c ... compute the vector multiply and it sign
         s = ax * by - ay * bx
         If (s .lt. 0) then
           i1 = IPE(1,i)
           IPE(1,i) = IPE(3,i)
           IPE(3,i) = i1
           IF (pbisE(i) .ne. 2) pbisE(i) = 4 - pbisE(i)
         End if
         
       End do  
c
c ... generate superelement for each node
c
       Do i = 1, nP+1
         iW(i) = 0
       End do
       
       Do i = 1, nE
         Do j = 1, 3
           k = IPE(j,i)
           iW(k+1) = iW(k+1) + 1
         End do
       End do
       
       iW(1)=1
       Do i = 2, nP+1
         iW(i) = iW(i-1) + iW(i) +1
         If (iW(i) + nP + 1  .gt. MaxWi) then 
           iERR = 1010
           call errMesRCB(iERR, 'Init_mesh_gen','MaxWi is too small')
         End if
       End do
       
       Do i = 1, nP
         j = nP + 1 + iW(i)
         iW(j)=0
       End do
       
       Do i = 1, nE
         Do j =1, 3
           k = IPE(j,i)
           m = nP + 1 + iW(k)
           iW(m) = iW(m) + 1
           l = iW(m)
           iW(m + l) = i
         End do
       End do
       
       Do i = 1, nP
         iW(i) = iW(i) + nP + 1
       End do
c        
c ... form IEF
c
       maxmetIPF=0
       Do i = 1, nF
         j = lbF(i)
         if (j .gt. maxmetIPF) maxmetIPF = j
       End do
       maxmetIPF = maxmetIPF +1

       Do i = 1,nE
         Do j = 1,3
           IEF(j,i)=-maxmetIPF
         End do
       End do
c
       Do i = 1, nE
         Do j = 1, 3
           If (IEF(j,i) .eq. -maxmetIPF) then

             i1 = IPE(zyk(j + 1),i)
             i2 = IPE(zyk(j + 2),i)
          
             i3 = iW(i1)
             i4 = iW(i2)
             Do l = i3 + 1, i3 + iW(i3)
               If (iW(l) .ne. i) then
                 Do  k = i4 + 1, i4 + iW(i4)
                   If (iW(k) .ne. i) then
                     If (iW(k) .eq. iW(l)) then
                     
                       IEF(j,i) = iW(k)
                       go to 10
                       
                     End if   
                   End if    
                 End do
               End if
             End do
             
             k1 = max(i1,i2)
             k2 = min(i1,i2)
             Do l = 1, nF
               m1 = max(IPF(1,l), IPF(2,l))  
               m2 = min(IPF(1,l), IPF(2,l))  
               If ((k1 .eq. m1) .and. (k2 .eq. m2)) then
                 IEF(j,i) = -lbF(l)
                 go to 10
               End if
             End do
             iERR = 1011
             call errMesRCB(iERR, 'Init_mesh_gen',
     &                'Mistake in the bound date')
           End if
10         Continue

         End do
       End do 
c
      Return
      End
      


c ======================================================================
c Recover boundary edges.
c ======================================================================
      Subroutine RestoreBndEdges(nF, nFmax, nE, 
     &                           IPE, IEF, IPF, lbF, zyk, iERR)
c ======================================================================
      implicit none

      Integer             nF, nFmax, nE
      Integer             IPE(3,*), IEF(3,*)
      Integer       IPF(2,*), lbF(*), zyk(6), iERR

      Integer             i,j,k

c ======================================================================

      nF = 0

      Do i = 1, nE
        Do j =1, 3
        
          k = IEF(j,i)
          If (k .le. 0) Then
            nF = nF + 1
            If(nF .gt. nFmax) then 
               iERR = 1004
               Call errMesRCB(iERR, 'Post_procesing', 
     &                        'nFmax is too small')
            End if
            IPF(1, nF) = IPE(zyk(j + 1), i)
            IPF(2, nF) = IPE(zyk(j + 2), i)
            lbF(nF)    = - k
          End if        
          
        End do
      End do
        
      Return
      End
