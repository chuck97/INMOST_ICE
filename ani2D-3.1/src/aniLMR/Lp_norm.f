C ==========================================================
      Subroutine Lp_norm(nP, Lp, Metric)
C ==========================================================
C  Routine computes the metric for L_p norm using the metric
C  generated for the maximum norm.
C
C     Lp    - norm for which the metric is to be adjusted:
C             Lp > 0  means  L_p     norm
C             Lp = 0  means  maximum norm (L_infinity)
C ==========================================================
      Real*8  Metric(3, *), Lp, det
C ==========================================================
      If (Lp.eq.0d0) Return
      If (Lp.lt.0d0) then
         Call errMesg(1000, 'Lp_norm', 'negative Lp')
      End if

      Do n = 1, nP
         det = Metric(1, n) * Metric(2, n) - Metric(3, n) ** 2

         det = det ** (-1D0 / (2 * Lp + 2))

         Do i = 1, 3
            Metric(i, n) = Metric(i, n) * det
         End do
      End do

      Return
      End



C ==========================================================
      Subroutine Lp_gradnorm(nP, Lp, Metric)
C ==========================================================
C  Routine computes the metric for L_p norm of the gradient
C  using the metric generated for the maximum norm.
C
C     Lp    - norm for which the metric is to be adjusted:
C             Lp > 0  means  L_p     norm
C             Lp = 0  means  maximum norm (L_infinity)
C ==========================================================
      Real*8  Metric(3, *), Lp, det
C ==========================================================
      If (Lp.eq.0d0) Return
      If (Lp.lt.0d0) Then
         Call errMesg(1000, 'Lp_gradnorm', 'negative Lp')
      End if

      Do n = 1, nP
         det = Metric(1, n) * Metric(2, n) - Metric(3, n) ** 2

         det = det ** (-1D0 / (2 + Lp))

         Do i = 1, 3
            Metric(i, n) = Metric(i, n) * det
         End do
      End do

      Return
      End


