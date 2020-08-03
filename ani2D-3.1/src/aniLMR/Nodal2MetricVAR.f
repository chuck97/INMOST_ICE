C=========================================================================
C  The routine generates the metric "metric" from a given continuous
C  piecewise linear nodal function "u" and the mesh "tri,vrt,bnd".
C
C  Limitations to the mesh: each boundary node has to have at least
C  one adjacent interior node or one adjacent boundary node which is 
C  adjacent to an interior one.
C=========================================================================
      Subroutine Nodal2MetricVAR(u, 
     &                           vrt,nv, tri,nt, bnd,nb, metric, 
     &                           Nrmem,rmem, Nimem,imem)
C=========================================================================
C  Input: u    - function defined at mesh vertices
C         nv   - the number of vertices
C         vrt  - coordinates of these vertices
C         nt   - the number of triangles
C         tri  - the connecticity table
C         nb   - the number of boundary edges
C         bnd  - the list of boundary edges
C
C  Output: metric - tensor metric
C
C  Work arrays: rmem - real*8 array of length Nrmem
C               imem - integer array of length Nimem
C=========================================================================
      implicit none

      Real*8            u(*)

C Mesh
      Integer           nv,nt,nb
      double precision  vrt(2,*)
      integer           tri(3,*), bnd(2,*)

C Metric
      double precision  metric(3,*)

C Work arrays
      integer           Nrmem,Nimem
      integer           imem(*)
      double precision  rmem(*)

C Local variables
      double precision  DX12,DY12,DX13,DY13,DET,Dxu,Dyu,DxNi,DyNi
      integer           i,j,k,m,it,iv,L1,L2,L3, Ni(3), ki,km 
      double precision  meas(3),weight(1000)
      integer           lbuf, kbuf(1000), kv
C  weight,kbuf: number of triangles sharing a node is a priori less than 1000
      logical           flag

C=========================================================================
c ... imemory check
      if (Nrmem.lt.nv.or.Nimem.lt.2*nv+3*nt) then
         write(*,*) 'Increase Nrmem to', nv
         write(*,*) 'Increase Nimem to',2*nv+3*nt
         Call errMesg(1000, 'Nodal2MetricVAR', 'see above')
      end if

c ... initialize BadNode,Nhosttri,Hosttri
      do iv = 1, nv
         imem(iv+nv) = 0
      End do

c ... generate Nhosttri
      do it = 1, nt
         do i = 1, 3
            imem(tri(i,it)+nv) = imem(tri(i,it)+nv) + 1
         End do
      End do

c ... pointers to Hosttri
      imem( 1 ) = 2*nv + 1
      do iv = 2, nv
         imem( iv ) = imem( iv-1 ) + imem( nv + iv-1 )
      End do

      do iv = 1, nv
         imem( nv + iv ) = 0
      End do

      do it = 1, nt
         do i = 1, 3
            iv = tri(i,it)
            k = imem( iv ) + imem( nv + iv )
            imem( k ) = it
            imem( nv + iv ) = imem( nv + iv ) + 1     
         End do
      End do

c ... restore BadNode, Nhosttri..
      do iv = 1, nv
         imem(iv) = 0
         imem(iv+nv) = 0
      End do

      do it = 1, nt
         do i = 1, 3
            imem(tri(i,it)+nv) = imem(tri(i,it)+nv) + 1
         End do
      End do

c ... mark boundary nodes
      do i = 1, nb
         imem(bnd(1,i)+nv) = -abs(imem(bnd(1,i)+nv))
         imem(bnd(1,i)) = 1
         imem(bnd(2,i)+nv) = -abs(imem(bnd(2,i)+nv))
         imem(bnd(2,i)) = 1
      End do

c ... initialize Mass and metric
      do iv = 1, nv
c        Mass(iv) = 0d0
         rmem(iv) = 0d0
         metric(1,iv) = 0d0
         metric(2,iv) = 0d0
         metric(3,iv) = 0d0
      End do

c ... recover the Hessian: H_ij  = M^{-1} integral_Om du/dx_i dNi/dx_j d Om
      ki = 1
      do iv = 1, nv
      do i  = 1, iabs(imem(nv+iv))
         it = imem(2*nv+ki)
         ki = ki + 1
   
         L1 = tri(1,it)
         L2 = tri(2,it)
         L3 = tri(3,it)

         DX12 = VRT(1,L1)-VRT(1,L2)
         DX13 = VRT(1,L1)-VRT(1,L3)
         DY12 = VRT(2,L1)-VRT(2,L2)
         DY13 = VRT(2,L1)-VRT(2,L3)
         DET = DX12*DY13-DX13*DY12

        if (imem(nv+iv).lt.0) then
          do j = 1, 3
           if (imem(tri(j,it)).eq.0) rmem(iv) = rmem(iv) + DABS(DET)/24
          End do
        else

c  ... nodal basis function
         Ni(1) = 0
         Ni(2) = 0
         Ni(3) = 0
         if (iv.eq.L1) Ni(1) = 1d0
         if (iv.eq.L2) Ni(2) = 1d0
         if (iv.eq.L3) Ni(3) = 1d0

c        Mass(iv) = Mass(iv) + DABS(DET)/6d0
         rmem(iv) = rmem(iv) + DABS(DET)/6d0

         Dxu = ( (u(L1)-u(L2))*DY13 - (u(L1)-u(L3))*DY12)/DET
         Dyu = (-(u(L1)-u(L2))*DX13 + (u(L1)-u(L3))*DX12)/DET

         DxNi = ( (Ni(1)-Ni(2))*DY13 - (Ni(1)-Ni(3))*DY12)/DET
         DyNi = (-(Ni(1)-Ni(2))*DX13 + (Ni(1)-Ni(3))*DX12)/DET

         metric(1,iv) = metric(1,iv) - DABS(DET)*Dxu*DxNi/2
         metric(2,iv) = metric(2,iv) - DABS(DET)*Dyu*DyNi/2
         metric(3,iv) = metric(3,iv) - DABS(DET)*(Dxu*DyNi+Dyu*DxNi)/4
       end if
      End do
      End do

C     Inner nodes..
      do iv = 1, nv
         if (imem(nv+iv).gt.0) then
            metric(1,iv) = metric(1,iv)/rmem(iv)
            metric(2,iv) = metric(2,iv)/rmem(iv)
            metric(3,iv) = metric(3,iv)/rmem(iv)
         end if
      End do

c ... extrapolate the metric from interior nodes to boundary and avearage 
c     the extrapolations. Exlude respective nodes from the set of bad ones.
      ki = 1
      do iv = 1, nv
       if (imem(nv+iv).gt.0) then
          do i = 1, imem(nv+iv)
             ki = ki + 1
          End do
       else
          lbuf = 0
          do i  = 1, -imem(nv+iv)
             it = imem(2*nv+ki)
             ki = ki + 1
             L1 = tri(1,it)
             L2 = tri(2,it)
             L3 = tri(3,it)

             DX12 = VRT(1,L1)-VRT(1,L2)
             DX13 = VRT(1,L1)-VRT(1,L3)
             DY12 = VRT(2,L1)-VRT(2,L2)
             DY13 = VRT(2,L1)-VRT(2,L3)

             DET = DX12*DY13-DX13*DY12

             do j = 1, 3
                kv = tri(j,it)
                flag = .true.
                do k = 1, lbuf
                   if (kbuf(k).eq.kv) then
                      flag = .false.
                      weight(k) = weight(k) + DABS(DET)/24
                   end if
                End do
                if (flag.and.imem(kv).eq.0) then
                   lbuf = lbuf + 1
                   kbuf(lbuf) = kv
                   weight(lbuf) = DABS(DET)/24
                end if
             End do
          End do
          do k = 1, lbuf
             kv = kbuf(k)
             do j = 1, 3
                metric(j,iv) = metric(j,iv) 
     &               + metric(j,kv)*weight(k)/rmem(iv)
             End do
             imem(nv+iv) = iabs(imem(nv+iv))
          End do
       end if
      End do
     
      do iv = 1, nv
         if (imem(iv).ne.0.and.imem(nv+iv).gt.0) then
            imem(nv+iv) = - imem(nv+iv)
            imem(iv)      = 0
         end if
      End do
                
c ... treating the BadNodes by secondary extrapolating
      km = 1
      do iv = 1, nv
c       if (BadNode(iv).eq.1) then
        if (imem(iv).eq.1) then
           i = 0
           do k = 1, 3
              meas(k) = 0d0
           End do
           do m = 1, -imem(nv+iv)
c             it = Hosttri(m,iv)
              it = imem(2*nv+km)
              km = km + 1
              do j = 1, 3
c                if (BadNode(tri(j,it)).eq.0) then
                 if (imem(tri(j,it)).eq.0) then
                    i = i + 1
                    do k = 1, 3
                       meas(k) = meas(k) + metric(k,tri(j,it))  
                    End do
                 end if
              End do
           End do

           If (i.eq.0) then
              Call errMesg(1000, 'Nodal2MetricVAR', 
     &                     'Unexpected logical path')
           Else
              do k = 1, 3
                 metric(k,iv) = meas(k)/i
              End do
           End if
        Else
           Do m = 1, iabs(imem(nv+iv))
              km = km + 1
           End do
        End if
      End do

c ... make the Hessian to be elliptic using SpectralModule from libanimba.a
      Do iv = 1, nv
         Call SpectralModule(metric(1,iv), det)
      End do

      Return
      End 


