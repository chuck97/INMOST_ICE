C ======================================================================
C Parameterization of the unit circle centered at the origin
C ======================================================================
      Subroutine CrvFunction_user(tc, xyc, iFnc)
C ======================================================================
C  The routine computes the Cartesian coordinates of point
C  xyc from its parametric coordinate tc.
C
C  tc     - the given parametric coordinate of point
C  xyc(2) - the Cartesian coordinate of the same point
C  iFnc   - the function number for computing
C
C  On input :  tc, iFnc
C  On output:  xyc(2)
C
C *** Remarks:
C        1. This is the default routine to link library libani2D-2.x.a
C ======================================================================
      Real*8  tc, xyc(2)
      Real*8  PI, R

C ======================================================================
      PI = 4*datan(1D0)
      R  = 0.5D0

      If(iFnc.EQ.1) Then
         xyc(1) = -2.0 + 4.0 * tc
         xyc(2) =  1.0 + 0.5 * datan(8*tc-4) / datan(4D0)

      Else If(iFnc.EQ.2) Then
         xyc(1) = -2.0 + 4.0 * tc
         xyc(2) = -1.0 - 0.5 * datan(8*tc-4) / datan(4D0)

      Else If(iFnc.EQ.3) Then
         xyc(1) = R * dcos( PI*tc )
         xyc(2) = R * dsin( PI*tc )

      Else If(iFnc.EQ.4) Then
         xyc(1) = -R * dcos( PI*tc )
         xyc(2) = -R * dsin( PI*tc )

      Else
         Write(*, '(A,I5)') 'Undefined function =', iFnc
         Stop
      End if

      Return
      End



