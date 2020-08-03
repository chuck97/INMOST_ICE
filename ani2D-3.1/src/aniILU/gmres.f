c---->------------------------------------------------------------------<
c  GMRES with "right" preconditioning
c  Templates for the solution if linear Systems...
c  http://www.netlib.org
c---->------------------------------------------------------------------<
      SUBROUTINE slgmres(
     >  prevec, IPREVEC, iW,rW, 
     >  matvec, IMATVEC, IA,JA,A,
     >  WORK, MW, NW, H, MH, NH,
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
      INTEGER   N, MW, NW, MH, NH, ITER, INFO, NUNIT
      REAL*8    RESID
      REAL*8    RHS(*), SOL(*), WORK(MW,NW), H(MH,NH)
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
c  NW       : input  : trailing dimension of workspace >= IRESTART+3
c  H        : work   : Workspace (MH,NH)
c  MH       : input  : >= IRESTART+1
c  NH       : input  : >= IRESTART+6
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
      PARAMETER ( ZERO = 0d0 , ONE = 1d0 )
c---->------------------------------------------------------------------<
c  Local Variables:
c
      INTEGER I, K, MAXIT, LH, IRESTART, ICR
      INTEGER JAV, JCS, JSN, JR, JS, JV, JW, JY, JHREF
      REAL*8  AA, BB, RNORM, TOL, TMP, TMP2
c
c---->------------------------------------------------------------------<
c  External BLAS, etc.:
c
      EXTERNAL  ddot,daxpy,dcopy,drot,drotg,dscal
      REAL*8    ddot
      INTRINSIC dsqrt, min, dabs
c---->------------------------------------------------------------------<
c
c    Test the input parameters.
c
      INFO = 0
c
      IRESTART = min( NW - 3 , MH - 1 , NH - 6 )
c
      if ( N .eq. 0 ) then
         return
      else if ( N .lt. 0 ) then
         INFO = -10
      else if ( MW .lt. N ) then
         INFO = -20
      else if ( IRESTART .lt. 1 ) then
         INFO = -30
      else if ( ITER .le. 0 ) then
         INFO = -40
      endif
c
      if ( INFO .ne. 0 ) return
c
      LH = IRESTART + 1
c---->------------------------------------------------------------------<
c  Save input iteration limit and convergence tolerance
c
      MAXIT = ITER
      TOL   = RESID
c---->
c  Alias workspace columns.
c
      JR  = 1
      JW  = JR + 1
      JAV = JW
      JV  = JAV + 1
c---->
c  Set up extra columns in H.
c
      JCS   = IRESTART + 1
      JSN   = JCS + 1
      JS    = JSN + 1
      JY    = JS + 1
      JHREF = JY + 1
c---->
c  Set initial residual (av is temporary workspace here).
c
      call dcopy( N, RHS, 1, WORK(1,JR), 1 )
c
      TMP = ddot( N, SOL, 1, SOL, 1 )
      if ( TMP .ne. ZERO ) then
        call matvec( IMATVEC, -ONE, SOL, ONE, WORK(1,JR) ,
     &               IA,JA,A )
      endif
c---->
      TMP = ddot( N, WORK(1,JR), 1, WORK(1,JR), 1 )
      RESID = dsqrt( TMP )
c---->
      ITER = 0
      ICR  = 0
      if ( RESID .lt. TOL ) GOTO 30
c---->------------------------------------------------------------------<
c  Restart point
c---->
   10 continue
c
c Increment the restart counter to possibly shift the preconditioner
c Reset the iteration counter
c
        ICR = ICR + 1
        I = 0
c
c  Construct the first column of V (gmres basis).
c
        call dcopy( N, WORK(1,JR), 1, WORK(1,JV), 1)
c---->--
        TMP = ddot( N, WORK(1,JV), 1, WORK(1,JV), 1 )
        TMP = dsqrt( TMP )
c---->--
        call dscal( N, ONE / TMP, work(1,JV), 1 )
c---->--
c  Initialize s to the elementary vector e1 scaled by 'beta'.
c
        H(1,JS) = TMP
        do K = 2, LH
           H(K,JS) = ZERO
        end do
c---->------------------------------------------------------------------<
c  GMRES normal iteration point
c---->--
   20   continue
c
          I = I + 1
          ITER = ITER + 1
c---->----
          call prevec( IPREVEC, ICR, WORK(1,JV+I-1), WORK(1,JR), iW,rW )
          call matvec( IMATVEC, ONE, WORK(1,JR), ZERO, WORK(1,JW) ,
     &               IA,JA,A )
c---->----
c  Construct I-th column of h orthnormal to the previous I-1 columns.
c
          call slgmresbasis(
     >      I, N, H(1,I), WORK(1,JV), MW,
     >      WORK(1,JW), H(1,JHREF))
c
c  Apply Givens rotations to the I-th column of H. This
c  "updating" of the QR factorization effectively reduces
c  the Hessenberg matrix to upper triangular form during
c  the restart iterations.
c
          do K = 1,I-1
            call drot(1,H(K,I),MH,H(K+1,I),MH,H(K,JCS),H(K,JSN))
          end do
c---->----
c  Construct the I-th rotation matrix, and apply it to h so that
c  H(i+1,i) = 0.
c
          aa = H(I,I)
          bb = H(I+1,I)
          call drotg( aa, bb, H(I,JCS), H(I,JSN) )
          call drot(1,H(I,I),MH,H(I+1,I),MH,H(I,JCS),h(I,JSN))
c---->----
c  Apply the i-th rotation matrix to [ s(i), s(i+1) ]'. this
c  gives an approximation of the residual norm.
c
          call drot(1,H(I,JS),MH,H(I+1,JS),MH,H(I,JCS),H(I,JSN))
          RESID = dabs( H(I+1,JS) )
c---->------------------------------------------------------------------<
c  Continue GMRES loop while:
c    1)  Less than maximum iteration
c    2)  Less than restart iteration
c    3)  Have not converged (estimated)
c         print*, ITER, RESID
          if (NUNIT.gt.0)    write(NUNIT,*) ITER, RESID , '[A'
c
          if ( ITER .lt. MAXIT .and.
     >         I .lt. IRESTART .and. RESID .ge. TOL ) go to 20
c---->------------------------------------------------------------------<
c  slgmresupdate generates  dX~ = M dX  =>  X = X + M^{-1} dX~
c
        call dcopy( N, ZERO, 0, WORK(1,JAV), 1 )
        call slgmresupdate(
     >         I,N,   WORK(1,JAV),
     >         H,MH,H(1,JY),H(1,JS),WORK(1,JV),MW)
       call prevec(IPREVEC,ICR,WORK(1,JAV),WORK(1,JR), iW,rW)
       call daxpy( N, ONE, WORK(1,JR), 1, SOL, 1 )

       if ( RESID .le. TOL ) go to 30 
c---->------
c  Compute the residual vector JR, its norm, and check for convergence
c  (JAV is temporary)
c
        call dcopy( N, RHS, 1, WORK(1,JR), 1 )
        call matvec( IMATVEC, -ONE, SOL, ONE, WORK(1,JR) ,
     &               IA,JA,A )
c---->--
        TMP = ddot( N, WORK(1,JR), 1, WORK(1,JR), 1 )
        H(I+1,JS) = dsqrt( TMP )
c
        RESID = H(I+1,JS)
c
c  Convergence failure?
c
        if ( ITER .ge. MAXIT .and. RESID .ge. TOL ) INFO = 1
c---->------------------------------------------------------------------<
c  Output
c
 30     continue
        TMP = ddot( N , SOL, 1, SOL, 1 )
        TMP = dsqrt( TMP )
        if ( NUNIT .gt. 0 ) then
c         WRITE(NUNIT,9000) ITER,RESID
c9000     FORMAT(3x,'SLGMRES ',I4,' : ',E16.10)
          WRITE(NUNIT,9000) ITER,RESID,TMP
 9000     FORMAT(3x,'SLGMRES ',I4,' : ',E16.10,' (SOL ',E16.10,')')
        end if
c---->------------------------------------------------------------------<
c Continue RESTART loop while:
c   1)  Not failed
c   2)  A restart iteration
c   3)  Not converged
c
        if ( INFO .eq. 0 .and.
     >       I .ge. IRESTART .and. RESID .ge. TOL ) go to 10
c---->------------------------------------------------------------------<
      return
      end

c---->------------------------------------------------------------------<
      SUBROUTINE slgmresupdate( I, N, X, H, LDH, Y, S, V, LDV )
c---->
      IMPLICIT NONE
c---->
      INTEGER   N, I, LDH, LDV
      REAL*8    X(*), Y(*), S(*), H(LDH,*), V(LDV,*)
      EXTERNAL  daxpy, dcopy, dtrsv
c---->
      REAL*8    ONE
      PARAMETER ( ONE = 1d0 )
c---->
c  Local Variables
c
c---->
c  Solve H * y = s for upper triangular H
c
      call dcopy( I, S, 1, Y, 1 )
      call dtrsv( 'U', 'N', 'N', I, H, LDH, Y, 1 )
c
c  Update current solution vector X = X + V * Y
c
      call dgemv( 'N', N, I, ONE, V, LDV, Y, 1, ONE, X, 1 )
c
      return
      end
c---->------------------------------------------------------------------<
c---->------------------------------------------------------------------<
      SUBROUTINE slgmresbasis( I, N, H, V, LDV, W, HREF )
      IMPLICIT NONE
c---->
      INTEGER   I, N, LDV
      REAL*8    H(*), W(*), HREF(*), V(LDV,*)
c---->
      REAL*8  ONE, ZERO, EPS, COS45, ETA, SQRTHALF
      PARAMETER ( ONE = 1d0 , ZERO = 0d0 )
      PARAMETER ( EPS = 1.0E-15 , SQRTHALF = 0.70710678 )
      PARAMETER ( ETA = EPS / SQRTHALF )
c---->
c  Local variables
c
      INTEGER K
      REAL*8  NORMW, TMP, TMP2
c---->
c  Extern BLAS, etc.
c
      REAL*8    ddot
      EXTERNAL  daxpy, ddot
      INTRINSIC dsqrt
c---->
c---->
      TMP = ddot( N, W, 1, W, 1 )
      NORMW = dsqrt( TMP )
c---->
c  Modified Gram-Schmidt
c
      do K = 1, I
        TMP = ddot( N, W, 1, V(1,K), 1 )
        H(K) = TMP
        call daxpy( N, -H(K), V(1,K), 1, W, 1 )
      end do
c---->
      TMP = ddot( N, W, 1, W, 1 )
      H(I+1) = dsqrt( TMP )
      call dcopy( N, W, 1, V(1,I+1), 1 )
      call dscal( N, ONE/H(i+1), V(1,I+1), 1 )
c---->
c  Iterative Refinement Process - Kahan Test for Refinement
c                                 See Parlett's Book, p. 107
c
      if ( H(I+1) .lt. ETA * NORMW ) then
        call dgemv('T',N,I,ONE,V,LDV,W,1,ZERO,HREF,1)
        call dgemv('N',N,I,-ONE,V,LDV,HREF,1,ONE,W,1)
        call daxpy(I,ONE,HREF,1,H,1)
c---->--
        TMP = ddot(N,W,1,W,1)
        TMP = dsqrt( TMP )
c---->--
        if ( TMP .lt. ETA * H(I+1) ) then
          do K = 1,N
            W(K) = 0.0
          end do
          TMP = 0.0
        end if
        H(I+1) = TMP
c---->--
        call dcopy(N,W,1,V(1,I+1),1)
c---->--
        if ( H(I+1) .ne. 0.0 ) then
          call dscal(N, ONE/H(I+1), V(1,I+1), 1 )
        end if
      end if
c---->
      return
      end

