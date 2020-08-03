C ======================================================================
      Integer Function ANI_MetricFunction(x, y, Metric)
C ======================================================================
C  This routine defines a metric at point (x,y). The
C  metric is a 2x2 positive definite symmetric tensor:
C
C                M11   M12
C      Metric =     
C                M12   M22
C
C  Only the upper triangular part of Metric must be defined.
C ======================================================================
      Real*8  x, y, Metric(2, 2)

      Metric(1,1) = 1D0
      Metric(2,2) = 1D0

      Metric(1,2) = 0D0

      ANI_MetricFunction = 0

      Return
      End



C ======================================================================
      Subroutine ANI_CrvFunction(tc, xyc, iFnc)
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
C  *** Remarks:
C         1. This is the dummy routine.
C
C ======================================================================
      Real*8  tc, xyc(2)

      Return
      End



