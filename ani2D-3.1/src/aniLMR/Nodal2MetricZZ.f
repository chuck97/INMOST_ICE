C=========================================================================
C  The routine generates the metric "measure" from a given continuous
C  piecewise linear function "u" and the mesh "tri, vrt, bnd".
C  The algorithm is based on the ZZ recovery method for interior nodes.
C=========================================================================
      Subroutine Nodal2MetricZZ(u, 
     &                          Vrt, nP, Tri,nE, measure, 
     &                          Nrmem,rmem, Nimem,imem )
C=========================================================================
C  Input: u    - function defined at mesh vertices
C         nP  - the number of vertices
C         Vrt  - coordinates of these vertices
C         nE  - the number of triangles
C         Tri  - the connecticity table
C
C  Output: measure - tensor metric
C
C  Work arrays: rmem - real*8 array of length Nrmem
C               imem - integer array of length Nimem
C=========================================================================
      implicit None
      Double precision  u(*)

C Mesh
      Integer           nP, nE, Tri(3,*)
      Double precision  vrt(2,*)

C Metric
      Double precision  measure(3,*)

C Work arrays
      Integer           Nrmem, Nimem
      Integer           imem(*)
      Double precision  rmem(*)

C Local variables
      double precision  DX12,DY12,DX13,DY13, DET, Dx,Dy
      integer           it,iv,L1,L2,L3, ip1,ip2,ip3,ip4, iERR

C=========================================================================
c ... memory check
      If(Nrmem.lt.4*nP+2*nE .OR. Nimem.lt.3*nP+3*nE) Then
         Write(*,'(A)') 'Required integer memory (Nimem):', 3*nP + 3*nE
         Write(*,'(A)') 'Required real*8  memory (Nrmem):', 4*nP + 2*nE

         Call errMesLMR(1001, 
     &       'Nodal2MetricZZ', 'Not enough integer or real*8 memory')
      End if

c ... pointers to rmem
      ip1 = 2*nP 
      ip2 = ip1 + nE 
      ip3 = ip2 + nE 
      ip4 = ip3 + nP 

c ... Recover the constant gradients
      do it = 1, nE
         L1 = tri(1,it)
         L2 = tri(2,it)
         L3 = tri(3,it)

         DX12 = Vrt(1,L1)-Vrt(1,L2)
         DX13 = Vrt(1,L1)-Vrt(1,L3)
         DY12 = Vrt(2,L1)-Vrt(2,L2)
         DY13 = Vrt(2,L1)-Vrt(2,L3)

         DET = DX12*DY13-DX13*DY12

         Dx = ((u(L1)-u(L2))*DY13 - (u(L1)-u(L3))*DY12)/DET
         Dy =(-(u(L1)-u(L2))*DX13 + (u(L1)-u(L3))*DX12)/DET

         rmem(ip1+it) = Dx
         rmem(ip2+it) = Dy
      end do

c ... recover the linear gradients
c     the allocated memory (checked above) is sufficient
      call P02P1(nP, nE, Vrt, Tri, 
     &           rmem(ip1+1), rmem(ip3+1), Nimem, imem, iERR)

      call P02P1(nP, nE, Vrt, Tri, 
     &           rmem(ip2+1), rmem(ip4+1), Nimem, imem, iERR)


c     the allocated memory (checked above) is sufficient
c ... recover the constant gradients
      do it = 1, nE
         L1 = tri(1,it)
         L2 = tri(2,it)
         L3 = tri(3,it)

         DX12 = Vrt(1,L1)-Vrt(1,L2)
         DX13 = Vrt(1,L1)-Vrt(1,L3)
         DY12 = Vrt(2,L1)-Vrt(2,L2)
         DY13 = Vrt(2,L1)-Vrt(2,L3)

         DET = DX12*DY13-DX13*DY12

         Dx = ((rmem(ip3+L1)-rmem(ip3+L2))*DY13 
     &       - (rmem(ip3+L1)-rmem(ip3+L3))*DY12)/DET
         Dy =(-(rmem(ip3+L1)-rmem(ip3+L2))*DX13 
     &       + (rmem(ip3+L1)-rmem(ip3+L3))*DX12)/DET

         rmem(ip1+it) = Dx
         rmem(ip2+it) = Dy
      end do


c ... recover the linear gradients
c     the allocated memore checked above is sufficient
      call P02P1(nP, nE, Vrt, Tri, 
     &           rmem(ip1+1), rmem(ip3+1), Nimem, imem, iERR)

      do iv = 1, nP
         measure(1,iv) = rmem(ip3+iv)
      end do

      call P02P1(nP, nE, Vrt, Tri, 
     &           rmem(ip2+1), rmem(ip3+1), Nimem, imem, iERR)

      do iv = 1, nP
         measure(3,iv) = 0.5d0 * rmem(ip3+iv)
      end do


c ... recover the constant gradients
      do it = 1, nE
         L1 = tri(1,it)
         L2 = tri(2,it)
         L3 = tri(3,it)

         DX12 = Vrt(1,L1)-Vrt(1,L2)
         DX13 = Vrt(1,L1)-Vrt(1,L3)
         DY12 = Vrt(2,L1)-Vrt(2,L2)
         DY13 = Vrt(2,L1)-Vrt(2,L3)

         DET = DX12*DY13-DX13*DY12

         Dx = ((rmem(ip4+L1)-rmem(ip4+L2))*DY13 
     &       - (rmem(ip4+L1)-rmem(ip4+L3))*DY12)/DET
         Dy =(-(rmem(ip4+L1)-rmem(ip4+L2))*DX13 
     &       + (rmem(ip4+L1)-rmem(ip4+L3))*DX12)/DET

         rmem(ip1+it) = Dx
         rmem(ip2+it) = Dy
      end do


c ... recover the linear gradients
c     the allocated memory checked above is suffisient
      call P02P1(nP, nE, Vrt, Tri, 
     &           rmem(ip1+1), rmem(ip4+1), Nimem, imem, iERR)

      do iv = 1, nP
         measure(3,iv) = measure(3,iv) + 0.5d0 * rmem(ip4+iv)
      end do

      call P02P1(nP, nE, Vrt, Tri, 
     &           rmem(ip2+1), rmem(ip4+1), Nimem, imem, iERR)

      do iv = 1, nP
         measure(2,iv) = rmem(ip4+iv)
      end do


c ... Make the Hessian elliptic using SpectralModule() from libanimba.a 
      do iv = 1, nP
         Call SpectralModule(measure(1,iv), det)
      end do

      return
      end 


