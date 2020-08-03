      SUBROUTINE slpcg(
     >  prevec, IPREVEC, iW,rW, 
     >  matvec, IMATVEC, IA,JA,A,
     >  WORK, MW, NW,
     >  N, RHS, SOL,
     >  ITER, RESID,
     >  INFO, NUNIT )
c---->
      IMPLICIT NONE
c---->------------------------------------------------------------------<
c  Argument types:
c
      EXTERNAL  matvec, prevec
      INTEGER   IMATVEC(*), IPREVEC(*)
      INTEGER   N, MW, NW, ITER, INFO, NUNIT
      REAL*8    RESID
      REAL*8    RHS(*), SOL(*), WORK(MW,NW)
      REAL*8    A(*),  rW(*)
      INTEGER   IA(*),JA(*),iW(*)

c---->
c  Argument Descriptions:
c
c  prevec   : extern : Precondition-vector routine
c  IPREVEC  : input  : Configuration data for 'prevec'
c  matvec   : extern : Matrix-vector multiply routine
c  IMATVEC  : input  : Configuration data for 'matvec'
c
c  WORK     : work   : Workspace (MW,NW)
c  MW       : input  : leading  dimension of workspace >= N
c  NW       : input  : trailing dimension of workspace >= 4
c
c  N        : input  : Length of vectors
c  RHS      : input  : RHS vector
c  SOL      : in/out : Initial guess / iterated solution
c  ITER     : in/out : Maximum iterations / actual iterations
c  RESID    : in/out : Convergence target / Norm of final residual
c  INFO     : output : = 0, converged
c                    : > 0, did not converge
c                    : < 0, error with input
c---->
c  External routine specifications:
c
c    matvec( IMATVEC, A, X, B, Y )  <=>  Y = A * Mat * X + B * Y
c    prevec( IPREVEC, i, X, Y )  <=>  Y = (MatP_{i})^{-1} * X
c      where MatP is the approximation of Mat
c---->------------------------------------------------------------------<
c  Local Parameters
c
      REAL*8    ZERO,ONE
      PARAMETER ( ZERO = 0.0 , ONE = 1.0 )
c---->------------------------------------------------------------------<
c  Local Variables:
c
      INTEGER MAXIT
      INTEGER JR, JP, JQ, JZ
      REAL*8  RHO, RHOPREV, ALPHA, BETA, TOL, TMP2
c
c---->------------------------------------------------------------------<
c  External BLAS, etc.:
c
      EXTERNAL  ddot,daxpy,dcopy,dscal
      REAL*8    ddot
      INTRINSIC sqrt
c---->------------------------------------------------------------------<
c
c    Test the input parameters.
c
      INFO = 0
c
      if ( N .eq. 0 ) then
         return
      else if ( N .lt. 0 ) then
         INFO = -10
      else if ( MW .lt. N ) then
         INFO = -20
      else if ( NW .lt. 4 ) then
         INFO = -30
      else if ( ITER .le. 0 ) then
         INFO = -40
      endif
c
      if ( INFO .ne. 0 ) return
c---->------------------------------------------------------------------<
c  Save input iteration limit and convergence tolerance
c
      MAXIT = ITER
      TOL   = RESID
c---->
c  Alias workspace columns.
c
      JR  = 1
      JP  = JR + 1
      JQ  = JP + 1
      JZ  = JQ + 1
c---->
c  Set initial residual
c
      call dcopy( N, RHS, 1, WORK(1,JR), 1 )
c
      TMP2 = ddot( N, SOL, 1, SOL, 1 )
      if ( TMP2 .ne. ZERO ) then
        call matvec( IMATVEC, -ONE, SOL, ONE, WORK(1,JR) ,
     &               IA,JA,A )
      endif
c---->
      TMP2 = ddot( N, WORK(1,JR), 1, WORK(1,JR), 1 )
      RESID = sqrt( TMP2 )
c---->
      ITER = 0
      if ( RESID .lt. TOL ) GOTO 20
c---->------------------------------------------------------------------<
c  PCG  iteration point
c---->--
   10   continue
c
          ITER = ITER + 1
c---->----
          call prevec( IPREVEC, 0, WORK(1,JR), WORK(1,JZ), iW,rW )

          RHOPREV = RHO
          RHO = ddot( N, WORK(1,JR), 1, WORK(1,JZ), 1 )

          IF (ITER.eq.1) THEN
             call dcopy( N, WORK(1,JZ), 1, WORK(1,JP), 1 )
          ELSE
             BETA = RHO / RHOPREV
             call dscal( N, BETA, WORK(1,JP), 1 )
             call daxpy( N, ONE, WORK(1,JZ), 1, WORK(1,JP), 1 )
          END IF

          call matvec( IMATVEC, ONE, WORK(1,JP), ZERO, WORK(1,JQ) ,
     &               IA,JA,A )

          TMP2 = ddot( N, WORK(1,JP), 1, WORK(1,JQ), 1 )
          ALPHA = RHO / TMP2

          call daxpy( N, ALPHA, WORK(1,JP), 1, SOL, 1 )
          call daxpy( N,-ALPHA, WORK(1,JQ), 1, WORK(1,JR), 1 )
c---->----
c  Check convergence
          TMP2 = ddot( N, WORK(1,JR), 1, WORK(1,JR), 1 )
          RESID = sqrt( TMP2 )
c---->------------------------------------------------------------------<
c  Continue PCG loop while:
c    1)  Less than maximum iteration
c    2)  Have not converged
c         print*,'pcg: ',ITER, RESID
       if (NUNIT.gt.0)    write(NUNIT,*) ITER, RESID , '[A'
c
          if ( ITER .lt. MAXIT .and. RESID .ge. TOL ) go to 10
c---->--
c
c  Convergence failure?
c
        if ( ITER .ge. MAXIT .and. RESID .ge. TOL ) INFO = 1
c---->------------------------------------------------------------------<
c  Output
c
  20    continue
        TMP2 = ddot( N , SOL, 1, SOL, 1 )
        TMP2 = sqrt( TMP2 )
        if ( NUNIT .gt. 0 ) then
          WRITE(NUNIT,9000) ITER,RESID,TMP2
 9000     FORMAT(3x,'SLPCG ',I4,' : ',E16.10,' (SOL ',E16.10,')')
        end if
c---->------------------------------------------------------------------<
      return
      end

