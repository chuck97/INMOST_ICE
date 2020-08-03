C ======================================================================
      Integer Function DexactU(x, y, label, dDATA, iDATA, iSYS, Coef)
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'
C ======================================================================
      Real*8  dDATA(*), x, y, Coef(4, 4)
      Integer iDATA(*), iSYS(*), label

      Real*8  x0,y0, x1,y1, a, f, g, r, phi, r0, r1, pi, e

C ======================================================================
      iSYS(1) = 1
      iSYS(2) = 1

c ... test 1 (H-function)
      If( .FALSE. ) Then
         x0 = 9.0D0
         y0 = 9.0D0
         a  = 2.0D0

         x1 = (x - x0)/6
         y1 = (y - y0)/6

         f = dsin(6*y1) - 3*x1
         g = dsin(6*x1) - 3*y1

         Coef(1,1) = (x1**2 * y1 + y1**3) / 16**3 + dtanh(a*f*g)
      End if


c ... test 2 (isotropic singularity)
      If( .FALSE. ) Then
         x0 = 0.5
         y0 = -0.2
         Coef(1,1) = ((x-x0)**2 - (dsqrt(1d1)*y-y0)**2)/
     &               ((x-x0)**2 + (dsqrt(1d1)*y-y0)**2)**2
      End if


c ... test 3: spiral function
      If( .FALSE. ) Then
         pi = 4D0 * datan(1D0)

         x0 = 10.0
         y0 =  9.0
         a = 2.0D0
         e = 1.15D0

         x1 = x - x0
         y1 = y - y0
         r = dsqrt(x1 * x1 + y1 * y1)

         If( r.EQ.0D0 ) Then
            phi = 0D0
         Else If(x1.GE.0D0 ) Then
            phi = dasin(y1 / r)
         Else If(x1.LT.0D0 ) Then
            phi = pi - dasin(y1 / r)
         End if

         r0 = a * (e ** phi) 
         r1 = a * (e ** (phi + 2*pi)) 

         Coef(1,1) = (x1**2 * y1 + y1**3) / 16**3
     &             + dtanh(2*dabs(r-r0)) * dtanh(2*dabs(r-r1))
      End if


c ... test 4 (snake)
      If( .TRUE. ) Then
         x0 = 0.0D0
         y0 = 0.0D0
         a = 6D0

         x1 = (x - x0)
         y1 = (y - y0)

         f = dsin(5*y1) - 2*x1

         Coef(1,1) = (x1**2 * y1 + y1**3) + dtanh(a*f)
      End if


c ... test 5
c     Coef(1,1) = 2 * x * y

      DexactU = TENSOR_SCALAR

      Return
      End



C ======================================================================
      Integer Function DgradU(x, y, label, dDATA, iDATA, iSYS, Coef)
C ======================================================================
      implicit none
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(4, 4)
      Integer iDATA(*), iSYS(*), label

      Real*8  u, a, f, g, x0,y0, x1,y1, v, ux,vx, uy,vy

C ======================================================================
      iSYS(1) = 2
      iSYS(2) = 1

c ... test 1 (H-function)
      If( .FALSE. ) Then
         x0 = 9.0D0
         y0 = 9.0D0
         a  = 2.0D0

         x1 = (x - x0)/6
         y1 = (y - y0)/6

         f = dsin(6*y1) - 3*x1
         g = dsin(6*x1) - 3*y1

         u = (dcosh(a*f*g))**2
         Coef(1,1) = 2*x1*y1 / 6 / 16**3           
     &             + a*(-g/2 + f*dcos(6*x1))/u
         Coef(2,1) = (x1**2 + 3*y1**2) / 6 / 16**3 
     &             + a*(-f/2 + g*dcos(6*y1))/u
      End if


c ... test 2 (isotropiv singularity)
      If( .FALSE. ) Then
         x0 = 0.5
         y0 =-0.2

         u = (x-x0)**2 - (dsqrt(1d1)*y-y0)**2
         v = (x-x0)**2 + (dsqrt(1d1)*y-y0)**2

         ux = 2*(x-x0)
         vx = 2*(x-x0)

         uy =-2*dsqrt(1d1) * (dsqrt(1d1)*y-y0)
         vy = 2*dsqrt(1d1) * (dsqrt(1d1)*y-y0)

         Coef(1,1) =( ux*v**2 - 2*u*v*vx ) / v**4
         Coef(2,1) =( uy*v**2 - 2*u*v*vy ) / v**4
      End if


c ... test 4 (snake)
      If( .TRUE. ) Then
         x0 = 0.0D0
         y0 = 0.0D0
         a = 6D0

         x1 = x - x0
         y1 = y - y0

         f = dsin(5*y1) - 2*x1
         u = (dcosh(a*f))**2

         Coef(1,1) = (2*x1*y1)         + a*(-2)/u
         Coef(2,1) = (x1**2 + 3*y1**2) + a*(dcos(5*y1)*5)/u
      End if


c ... test 5
c     Coef(1,1) = 2*x
c     Coef(2,1) = 2*y

      DgradU = TENSOR_GENERAL

      Return
      End




