C ======================================================================
      Subroutine SpectralModule(HesP, detG)
C ======================================================================
      implicit none
      include 'magic.fd'
C ======================================================================
C Routines sets the minimal eigenvalue of HesP to some constant
C ======================================================================
      Real*8  HesP(3), detG
    
C LOCAL VARIABLES
      Real*8  A(2, 2), E(2), rW(10)
      Integer info

C ======================================================================
      A(1, 1) = HesP(1)
      A(2, 2) = HesP(2)
      A(1, 2) = HesP(3)

      Call dsyev('V', 'U', 2, A, 2, E, rW, 10, info)
      If(info.NE.0) Call errMesMBA(3011, 'SpectralModule', 
     &                            'Error in Lapack routine dsyev')

      E(1) = dabs(E(1))
      E(2) = dabs(E(2))

      E(1) = max(E(1), E(2) * AniRatio)

      If(E(2).EQ.0D0) Then
         E(1) = AniEigenvalue
         E(2) = AniEigenvalue
      End if

      HesP(1) = E(1) * A(1, 1) ** 2 + E(2) * A(1, 2) ** 2
      HesP(2) = E(1) * A(2, 1) ** 2 + E(2) * A(2, 2) ** 2
      HesP(3) = E(1) * A(1, 1) * A(2, 1) + E(2) * A(1, 2) * A(2, 2)

c     detG = HesP(1) * HesP(2) - HesP(3) ** 2
      detG = E(1) * E(2)

      Return
      End



C ======================================================================
      Subroutine SpectralModuleOpt(HesP, detG, ek)
C ======================================================================
      implicit none
      include 'magic.fd'
C ======================================================================
C Routines sets the minimal eigenvalue of HesP to some constant
C ======================================================================
      Real*8  HesP(3), detG, ek(2, 3)
    
C LOCAL VARIABLES
      Real*8  V(2, 2), E(2), rW(10), a, b, v1, v2, mu, lambda, s
      Integer i, k, info

C ======================================================================
      V(1, 1) = HesP(1)
      V(2, 2) = HesP(2)
      V(1, 2) = HesP(3)

      Call dsyev('V', 'U', 2, V, 2, E, rW, 10, info)
      If(info.NE.0) Call errMesMBA(3011, 'SpectralModuleOpt', 
     &                            'Error in Lapack routine dsyev')

      detG = E(1) * E(2)

      E(1) = dabs(E(1))
      E(2) = dabs(E(2))

      E(1) = max(E(1), E(2) * AniRatio)

      If(E(2).EQ.0D0) Then
         E(1) = AniEigenvalue
         E(2) = AniEigenvalue
      End if

      If(detG.GT.0) Then
c ... return original metric
         HesP(1) = E(1) * V(1,1) ** 2 + E(2) * V(1, 2) ** 2
         HesP(2) = E(1) * V(2,1) ** 2 + E(2) * V(2, 2) ** 2
         HesP(3) = E(1) * V(1,1) * V(2, 1) + E(2) * V(1, 2) * V(2, 2)
      Else
c ... normalize eigenvectors by \sqrt{eigenvalue}
         Do i = 1, 2
            V(1, i) = V(1, i) * dsqrt(E(i))
            V(2, i) = V(2, i) * dsqrt(E(i))
         End do

c ... calculate  a = mu + lambda,  b = mu - lambda
         a = 0D0
         b = 0D0
         Do k = 1, 3
            v1 = V(1,1) * ek(1, k) + V(2,1) * ek(2, k)
            v2 = V(1,2) * ek(1, k) + V(2,2) * ek(2, k)

            a = a + (v1 + v2) ** 2
            b = b + (v1 - v2) ** 2
         End do

c ... recover mu & lambda
         mu     = (a + b) / 2
         lambda = (a - b) / 2

         s = 1D0 / dsqrt(mu ** 2 - lambda ** 2)
         mu     = s * mu 
         lambda = s * lambda

c ... calculate the metric
         HesP(1) = mu * (V(1,1) ** 2 + V(1,2) ** 2) 
     &           - 2*lambda * V(1,1) * V(1,2) 
         HesP(2) = mu * (V(2,1) ** 2 + V(2,2) ** 2) 
     &           - 2*lambda * V(2,1) * V(2,2) 
         HesP(3) = mu     * (V(1,1) * V(2,1) + V(1,2) * V(2,2))
     &           - lambda * (V(1,1) * V(2,2) + V(1,2) * V(2,1))

         detG = -detG
      End if

      Return
      End
