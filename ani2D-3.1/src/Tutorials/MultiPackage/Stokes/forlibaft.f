c ======================================================================
      Subroutine userboundary(i, t, x, y)
c ======================================================================
      Integer i
      Real*8  t, x(*), y(*), PI, R
c ======================================================================
      If(i.EQ.1) Then
        x(1) = -2D0 + 4.0D0 * t  ! 0 <= t <= 1.0 
        y(1) =  1D0 + 0.5D0 * datan(8*t-4) / datan(4.0D0)

      Else If(i.EQ.2) Then
        x(1) = -2D0 + 4.0D0 * t  ! 0 <= t <= 1.0 
        y(1) = -1D0 - 0.5D0 * datan(8*t-4) / datan(4.0D0)

      Else If(i.EQ.3) Then
        PI = 4 * datan(1D0)
        R  = 0.5D0
 
        x(1) = R * dcos(PI * t)  ! 0 <= t <= 1.0 
        y(1) = R * dsin(PI * t)

      Else If(i.EQ.4) Then
        PI = 4 * datan(1D0)
        R  = 0.5D0
 
        x(1) = -R * dcos(PI * t)  ! 0 <= t <= 1.0
        y(1) = -R * dsin(PI * t)
      End if

      Return
      End

