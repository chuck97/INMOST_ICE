C ======================================================================
      Subroutine TriMetric(xy1, xy2, xy3, Error, H)
C ======================================================================
C Routine computes the metric H for a triangle given by its vertices.
C ======================================================================
      implicit none

      include 'magic.fd'
C ======================================================================
      Real*8   xy1(2), xy2(2), xy3(2)
      Real*8   Error(3), H(3)

c Local variables
      Integer  i, j, i1, i2, i3
      Real*8   B(3, 3), alpha(3), s, p, t

      Integer  iref(4)
      DATA     iref/1,2,3,1/

C ======================================================================
c ... calculate matrix B = b_i . b_j where b_i is the bubble function 
      s = 1D0 / 180
      t = 1D0 / 360

      Do i = 1, 3
         B(i, i) = s
         Do j = i + 1, 3
            B(i, j) = t 
            B(j, i) = t 
         End do
      End do

c ... define edge errors
      Do i1 = 1, 3
         i2 = iref(i1 + 1)
         i3 = iref(i2 + 1)

         alpha(i3) = Error(i1)
      End do

c ... fix alphas
      p = alpha(1) + alpha(2) + alpha(3)

      If(p.EQ.0D0) Then
         alpha(1) = 1D-8
         p = alpha(1)
      End if 

c ... calculate the norm of error
      s = 0D0
      Do i = 1, 3
          Do j = i, 3
            s = s + B(i, j) * alpha(i) * alpha(j)
          End do
      End do

      s = dsqrt(s) / p

      Do i = 1, 3
         H(i) = dabs(alpha(i)) * s
      End do

c ... find the SPD matrix H such that (H e_i, e_i) = alpha(i), 
      Call updateH(xy1, xy2, xy3, H)

      Return
      End



C=======================================================================
      Subroutine TriGradMetric(xy1, xy2, xy3, Error, H, errL2, errH1)
C=======================================================================
C The routine computes the gradient metric H for a triangle given by
C its vertices. 
C=======================================================================
      implicit none

      Real*8   xy1(2), xy2(2), xy3(2)
      Real*8   Error(3), H(3), errL2, errH1

c Local variables
      Integer  i, j, k, i1, i2, i3
      Real*8   GradLambda(2, 3), Lambda(3), area, s, p
      Real*8   B(3, 3), alpha(3)

      Integer  delta(3, 3), iref(4)
      DATA     delta/0,2,0, 0,0,2, 2,0,0/
      DATA     iref/1,2,3,1/

C=======================================================================
c ... calculate matrix B = \grad b_i . \grad b_j where b_i is the bubble
c     the bubble is 1/4 at the edge mid-point (see the articles)
      Call GradU(xy1, xy2, xy3, 1D0, 0D0, 0D0, GradLambda(1,1), area)
      Call GradU(xy1, xy2, xy3, 0D0, 1D0, 0D0, GradLambda(1,2), area)
      Call GradU(xy1, xy2, xy3, 0D0, 0D0, 1D0, GradLambda(1,3), area)

      Do i = 1, 3
         Lambda(i) = GradLambda(1,i)**2 + GradLambda(2,i)**2
      End do

      Do i = 1, 3
         Do j = i, 3
            s = 0D0
            Do k = 1, 3
               s = s + Lambda(k) * (1 - delta(i,k)) * (1 - delta(j,k))
            End do
            B(i, j) = s / 12
            B(j, i) = B(i, j)
         End do
      End do

c ... compute alpha's using the given discrete function. 
c ... note that bubble i3 is defined of edge [i1 i2] 
      Do i1 = 1, 3
         i2 = iref(i1 + 1)
         i3 = iref(i2 + 1)

         alpha(i3) = Error(i1)
      End do

c ... fix alphas
      p = alpha(1) + alpha(2) + alpha(3)

      If(p.EQ.0D0) Then
         alpha(1) = 1D-8
         p = alpha(1)
      End if 

c ... calculate the energy norm of error
      s = 0D0
      Do i = 1, 3
         Do j = 1, 3
            s = s + B(i, j) * alpha(i) * alpha(j)
c           s = s + B(i, j) * Error(i) * Error(j)
         End do
      End do

      errH1 = s * area
      errL2 = (alpha(1)**2 + alpha(2)**2 + alpha(3)**2) * area / 180

      s = s / p

c ... gamma's (H) are scaled by alphas
      Do i = 1, 3
         H(i) = s * alpha(i) 
      End do

c ... find the SPD matrix H such that (H e_i, e_i) = alpha(i), 
      Call updateH(xy1, xy2, xy3, H)

      Return
      End



C=======================================================================
      Subroutine TriGradMetricVector(xy1, xy2, xy3, Error1, Error2, H)
C=======================================================================
C The routine computes the gradient metric H for a triangle given by
C its vertices. The metric is common for TWO functions (e.g.velocity 
C components).
C=======================================================================
      implicit none

      Real*8   xy1(2), xy2(2), xy3(2)
      Real*8   Error1(3), Error2(3), H(3)

c Local variables
      Integer  i, j, k, i1, i2, i3
      Real*8   GradLambda(2, 3), Lambda(3), area, s, p
      Real*8   B(3, 3), alpha1(3), alpha2(3)

      Integer  delta(3, 3), iref(4)
      DATA     delta/0,0,2, 2,0,0, 0,2,0/
      DATA     iref/1,2,3,1/

C=======================================================================
c ... calculate matrix B = \grad b_i . \grad b_j where b_i is the bubble
      Call GradU(xy1, xy2, xy3, 1D0, 0D0, 0D0, GradLambda(1,1), area)
      Call GradU(xy1, xy2, xy3, 0D0, 1D0, 0D0, GradLambda(1,2), area)
      Call GradU(xy1, xy2, xy3, 0D0, 0D0, 1D0, GradLambda(1,3), area)

      Do i = 1, 3
         Lambda(i) = GradLambda(1,i)**2 + GradLambda(2,i)**2
      End do

      Do i = 1, 3
         Do j = 1, 3
            s = 0D0
            Do k = 1, 3
               s = s + Lambda(k) * (1 - delta(i,k)) * (1 - delta(j,k))
            End do
            B(i, j) = s / 12
         End do
      End do


c ... compute alpha's using the given discrete function. 
c ... note that bubble i3 is defined of edge [i1 i2] 
      Do i1 = 1, 3
         i2 = iref(i1 + 1)
         i3 = iref(i2 + 1)

         alpha1(i3) = Error1(i1)
         alpha2(i3) = Error2(i1)
      End do

c ... fix alphas
      p = alpha1(1) + alpha2(1) 
     &  + alpha1(2) + alpha2(2) 
     &  + alpha1(3) + alpha2(3) 
c     p = max( alpha1(1) , alpha2(1) )
c    &  + max( alpha1(2) , alpha2(2) )
c    &  + max( alpha1(3) , alpha2(3) )

      If(p.EQ.0D0) Then
         alpha1(1) = 1D-8
         p = alpha1(1)
      End if 

c ... calculate the energy norm of error
      s = 0D0
      Do i = 1, 3
          Do j = 1, 3
             s = s + B(i, j) * alpha1(i) * alpha1(j)
     &             + B(i, j) * alpha2(i) * alpha2(j)
          End do
      End do

      s = s / p

c ... gamma's (H) are scaled by alphas
      Do i = 1, 3
         H(i) = s * (alpha1(i) + alpha2(i)) 
c        H(i) = s * max(alpha1(i) , alpha2(i)) 
      End do

c ... find the SPD matrix H such that (H e_i, e_i) = alpha(i), 
      Call updateH(xy1, xy2, xy3, H)

      Return
      End



C=======================================================================
      Subroutine updateH(xy1, xy2, xy3, H)
C=======================================================================
C  Routine finds the SPD matrix H such that (M e_i, e_i) = H(i)
C  and conditions of Lemma 1 are satisfied. The result is returned in H.
C
C  If H is close to singular, we modify it as described in Lemma 1. This
C  must be verified. If H is indefinite, we take its spectral module. 
C=======================================================================
      implicit none

      include 'magic.fd'
      include 'fem2Dtri.fd'

      Real*8   xy1(*), xy2(*), xy3(*), H(3)

      Integer  i, k, iLoop, info, IPIV(3)
      Real*8   ax(2, 3), B(3, 3), H0(3), s, det

C=======================================================================
      Do i = 1, 2
         ax(i, 1) = xy3(i) - xy2(i)
         ax(i, 2) = xy1(i) - xy3(i)
         ax(i, 3) = xy2(i) - xy1(i)
      End do

      Do i = 1, 3
         H0(i) = H(i)
      End do

      Do iLoop = 1, 2
         Do i = 1, 3
            B(i, 1) = ax(1, i) ** 2
            B(i, 2) = ax(2, i) ** 2
            B(i, 3) = 2 * ax(1, i) * ax(2, i)

            H(i) = H0(i)
         End do

         Call dgesv(3, 1, B, 3, IPIV, H, 3, info)
         If(info.NE.0) Call errMesLMR(3011, 'updateH',
     &                      'Error in the LAPACK routine dgesv')

         s   = H(1) * H(2)
         det = s - H(3) ** 2

         If(dabs(det).LT.1D-8*dabs(s)) Then
            If(iLoop.EQ.2) Call errMesLMR(6007, 'updateH',
     &                          'System error: H must be nonsingular')

            k = 1
            Do i = 2, 3
               If(H0(i).GT.H0(k)) k = i
            End do
 
            If(H0(k).GT.0D0) Then
               H0(k) = H0(k) * 1.0001
            Else
               H(1) = AniEigenvalue
               H(2) = AniEigenvalue
               H(3) = 0d0

               Goto 100
            End if
         Else
            Call SpectralModule(H, det)
c ... a metric for indefinite Hessians
c           Call SpectralModuleOpt(H, det, ax)
            Goto 100
         End if
      End do

 100  Continue

      Return
      End

