c ======================================================================
c The program reads data for the Bratu problem, initializes ILU2 
c preconditioner, and solves the Bratu problem using the inexact Newton 
c method with inner linear BiCGstab iterations (assuming that the Jacobian 
c is not known). The precondtioner for BiCGstab is an ILU2 preconditioner 
c for the mesh Laplacian operator.
c ======================================================================
c
c  Delta u + d  u'_x + lambda exp(u) = 0   in  domain
c                                 u  = 0   on  boundary
c
c  Parameters:
c      n       - dimension of the discrete problem 
c      x       - vector of the initial guess and the solution
c
c      nx,ny   - number of grid nodes in X and Y directions, n = nx*ny
c      lamdba  - Bratu problem parameter (e.g.1)
c      d       - Bratu problem parameter (e.g.1)
c ======================================================================
      Program  mainBRATU
      Implicit None
        
c ... maximum size of the discrete problem
      Integer   maxn
      Parameter(maxn = 250 000)

c ... work arrays
c  11*maxn - real*8  space for inexact Newton method with inner BiCGstab solver
c  20*maxn - real*8  space for preconditioner and matrix
c  30*maxn - integer space for preconditioner and matrix
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 11*maxn + 20*maxn, MaxWi = 30*maxn)

      Double precision rW(MaxWr)
      Integer          iW(MaxWi)

c ... Bratu problem data
      Integer          nx, ny
      Double precision h, cl,cr, h2l, lambda, d
     
c ... pointers for matrix kept in CSR format
      Integer  iPia, iPja, rPa

c ... inexact Newton data
      EXTERNAL fBratu, ddot, prevec2, laplace_csr 

c  fbratu is the vector function of the nonlinear system
c  prevec2 is the preconditioner evaluation for the jacobian
c  make_csr generates a matrix for which the preconditioner is initialized
      Integer          IPREVEC, N, iParBratu(2), INFO(5), rPwork,LenrW
      Double Precision RESID, STPTOL, SOL(maxn), rParBratu(3) 

c
c  INFO(5) is the control array for driver routine slInexactNewton() in the 
c  libinb-2.x.a library. On  input, it is as follows:
c       INFO(1) <= initial value for successful termination flag (=0)
c       INFO(2) <= maximal number of linear iterations per Newton step 
c       INFO(3) <= maximal number of nonlinear iterations 
c       INFO(4) <= maximal number of backtracks   
c       INFO(5) <= printing level (0 kepps mum)
c
c  On output, it is as follows:
c       INFO(1) - termination flag; its value has the folowing meanings:
c              0 => normal termination: ||F||.le.RESID or ||step||.le.STPTOL.
c              1 => maximum nonlinear iterations reached without success.
c              2 => failure to evaluate F.
c              3 => in jvnewton, J*v failure.
c              4 => in jvnewton, P(inverse)*v failure.
c              5 => in drvnewton, insufficient initial model norm reduction
c                   for adequate progress. NOTE: This can occur for several
c                   reasons; examine itrmks on return from the Krylov
c                   solver for further information. 
c              6 => in btnewton, failure to reach an acceptable step through
c                   backtracking.
c              7 => insufficient work memory for slInexactNewton (should be at least 11*N)
c       INFO(2) - number of linear iterations
c       INFO(3) - number of nonlinear iterations
c       INFO(4) - number of backtracks
c       INFO(5) - number of function evaluations
 

c ... for aniilu-2.x library
      Double Precision tau1, tau2, partlur, partlurout
      Integer          verb, ierr, UsedWr, UsedWi, iPilu, rPilu
      Integer          MaxWrILU, MaxWiILU


c ... local variables
      Integer          i


C ======================================================================
C  STAGE 1: read and set parameters of Bratu problem
C ======================================================================
      nx = 20
      d  = 1
      lambda = 1

      Write(*,'(/,A,I4,2(A,E12.4))') 
     &   ' Parameters for the Bratu problem: nx =', nx, 
     &   ', d =', d, ', lambda =', lambda
        
      ny = nx
      iParBratu(1) = nx
      iParBratu(2) = ny
      N = nx*ny
      h = 1.d0/dfloat(nx + 1) 

      cl  = 1.d0 - h*d/2.d0 
      cr  = 1.d0 + h*d/2.d0 
      h2l = h*h*lambda

      rParBratu(1) = cl
      rParBratu(2) = cr
      rParBratu(3) = h2l

      If(N.gt.maxn) Then
         Write(*,'(A)') 'increase maxn'
         Stop 911
      End if


c=======================================================================
c  STAGE 2: The matrix of discrete Laplacian in the CSR format
c           It will be used to define the preconditioner        
c ======================================================================
      iPia = 1
      iPja = iPia + N+1
      rPa  = 1

      Call laplace_csr(nx-1, iW(iPia), iW(iPja), rW(rPa))


C ======================================================================
C  STAGE 3: Initialization of the ILU2 preconditioner
C ======================================================================
      verb = 0
      tau1 = 5d-2
      tau2 = 5d-3
      partlur = 0.5d0
      ierr = 0

      iPilu = iPja + 5*N
      rPilu = rPa + 5*N

      MaxWrILU = MaxWr - rPilu
      MaxWiILU = MaxWi - iPilu

      Call iluoo(n, iW(iPia), iW(iPja), rW(rPa), tau1, tau2, verb,
     &           rW(rPilu), iW(iPilu), MaxWrILU, MaxWiILU,
     &           partlur,partlurout, UsedWr, UsedWi, ierr)

      If(ierr.ne.0) Then
         Write(*,'(A)') 'Initialization of iluoo failed, ierr=', ierr
         Stop 911
      End if

c Shift the preconditioner data to the beginning of work array 
      Do i = 1, UsedWi
         iW(i) = iW(i + iPilu - 1)
      End do

      Do i = 1, UsedWr
         rW(i) = rW(i + rPilu - 1)
      End do


C ======================================================================
C  STAGE 4: iterative solution by slInexactNewton
C ======================================================================
      IPREVEC = N
      RESID   = 1d-10
      STPTOL  = 1d-7

      INFO(1) = 0    ! initializing successful termination flag for Newton
      INFO(2) = 1000 ! maximal number of linear iterations 
      INFO(3) = 100  ! maximal number of nonlinear iterations 
      INFO(4) = 10   ! maximal number of backtracks   
      INFO(5) = 0    ! print level (0 keeps mum)

      rPwork = rPa + UsedWr
      LenrW  = MaxWr - rPwork + 1

c Initial guess for nonlinear iterations
      Do i = 1, N
         SOL(i) = 0.d0
      End do
        
      Call slInexactNewton(prevec2, IPREVEC, iW, rW,   
     &                     fBratu, rParBratu, iParBratu,
     &                     N, SOL,
     &                     RESID, STPTOL, 
     &                     rW(rPwork), LenrW,
     &                     INFO)
        
      If(INFO(1).ne.0) Then
         Write(*,*) 'Failed to solve the problem INFO=', INFO(1)
         Stop 911
      Else
         Write(*,*)
         Write(*,*) '||F(SOL)|| = ', RESID
         Write(*,*) 'Number of linear iterations:   ', INFO(2)
         Write(*,*) 'Number of nonlinear iterations:', INFO(3)
         Write(*,*) 'Number of backtracks:          ', INFO(4)
         Write(*,*) 'Number of function evaluations:', INFO(5)
      End if
        
      Stop
      End



C ======================================================================
c This user routine evaluates a nonlinear residual for the Bratu problem.
c
c INPUT:
c    n       dimension of vectors
c    xcur    current vector 
c    rpar    double precision user-supplied parameters
c    ipar    integer user-supplied parameters
c
c OUTPUT:
c    fcur    nonlinear residual vector (zero for the solution)
c    itrmf   flag for successful termination of the routine
C ======================================================================
      Subroutine fBratu(n, xcur, fcur, rpar, ipar, itrmf)
C ======================================================================
      Implicit none
      Integer          i, itrmf, j, j1, j2, n, ipar(*), nx, ny 
      Double precision cl, cr, h2l, xcur(n), fcur(n), rpar(*)

C ======================================================================
c  Set local variables from user-supplied parameters.
      nx  = ipar(1)
      ny  = ipar(2)
      cl  = rpar(1) 
      cr  = rpar(2)
      h2l = rpar(3)

c  Evaluate nonlinear residual.
      Do 100 j = 1, ny 
         j1 = (j - 1)*nx + 2 
         j2 = j*nx - 1

         Do 110 i = j1, j2
            fcur(i) = cr*xcur(i+1) + cl*xcur(i-1) - 4.d0*xcur(i) 
     $              + h2l*dexp(xcur(i)) 
 110     Continue

         j1 = j1 - 1
         fcur(j1) = cr*xcur(j1+1) - 4.d0*xcur(j1) + h2l*dexp(xcur(j1)) 

         j2 = j2 + 1 
         fcur(j2) = cl*xcur(j2-1) - 4.d0*xcur(j2) + h2l*dexp(xcur(j2)) 

         If(j.ne.1) then 
            Do 120 i = j1, j2 
               fcur(i) = fcur(i) + xcur(i-nx)
 120        Continue
         Endif

         If(j.ne.ny) then 
            Do 130 i = j1, j2 
               fcur(i) = fcur(i) + xcur(i+nx)
 130        Continue
         Endif
 100  Continue 

c  Set the termination flag to success.

      itrmf = 0

      Return
      End



C ======================================================================
      Subroutine laplace_csr(nx,ia,ja,a)
C ======================================================================
c This routine generate the 5-point Laplacian matrix and stores it
c in the compressed sparce row format.
C ======================================================================
      Implicit none
      Integer   nx,N,str
      Real *8   a(5*(nx+1)**2 - 2*(nx+1))
      Integer   ia((nx+1)**2 + 1), ja(5*(nx+1)**2 - 2*(nx+1))

C ======================================================================
      ia(1) = 1
      ia(2) = 4
      ja(1) = 1
      ja(2) = 2
      ja(3) = 2 + nx

      a(1) =  4.d0
      a(2) = -1.d0
      a(3) = -1.d0

      N = (nx+1)**2
      Do str = 2, N-1
         If(str.le.nx+1) Then
            ia(str+1) = ia(str) + 4
                        
            ja(ia(str) + 0) = str - 1
             a(ia(str) + 0) = -1.d0 

            ja(ia(str) + 1) = str
             a(ia(str) + 1) = 4.d0

            ja(ia(str) + 2) = str + 1 
             a(ia(str) + 2) = -1.d0

            ja(ia(str) + 3) = str + 1+nx 
             a(ia(str) + 3) = -1.d0

         Else If(str.le.N - (nx+1)) Then
            ia(str + 1) = ia(str) + 5

            ja(ia(str) + 0) = str - (nx+1)
             a(ia(str) + 0) = -1.d0

            ja(ia(str) + 1) = str - 1
             a(ia(str) + 1) = -1.d0

            ja(ia(str) + 2) = str
             a(ia(str) + 2) = 4.d0

            ja(ia(str) + 3) = str +1
             a(ia(str) + 3) = -1.d0

            ja(ia(str) + 4) = str +  nx+1
             a(ia(str) + 4) = -1.d0

         Else
            ia(str+1) = ia(str) + 4

            ja(ia(str) + 0) = str - (nx+1)
             a(ia(str) + 0) = -1.d0

            ja(ia(str) + 1) = str - 1
             a(ia(str) + 1) = -1.d0 

            ja(ia(str) + 2) = str
             a(ia(str) + 2) = 4.d0

            ja(ia(str) + 3) = str + 1 
             a(ia(str) + 3) = -1.d0
         End if
      End do

      Do str = nx+1,N - (nx+1), nx+1
         a(ia(str) + 3) = 0.d0
      End do

      Do str = nx+2, N - nx , nx+1
         a(ia(str) + 1) = 0.d0
      End do

      ia(N+1) = ia(N) + 3

      ja(ia(N)) = N - (nx+1)
       a(ia(N)) = -1.d0

      ja(ia(N) + 1) = N - 1
       a(ia(N) + 1) = -1.d0

      ja(ia(N) + 2) = N
       a(ia(N) + 2) = 4.d0

       a(ia(nx+1) + 2) = 0.d0
       a(ia(N-nx) + 1) = 0.d0 

      Return
      End
                                                     

