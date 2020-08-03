C ==========================================================
C Parameterization if the unit circle centered at the origin
C ==========================================================
      Subroutine CrvFunction_user(tc, xyc, iFnc)
C ==========================================================
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
C        1. This is the default routine to link library libani2D.a
C ==========================================================
      Real*8  tc, xyc(2)
      Real*8  PI

C ==========================================================
      If(iFnc.EQ.1) Then
         PI=4*datan(1d0)
         xyc(1) = dcos( PI*tc )
         xyc(2) = dsin( PI*tc )
      Else
         Write(*, '(A,I5)') 'Undefined function =', iFnc
         Stop
      End if

      Return
      End



