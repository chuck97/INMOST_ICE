C ======================================================================
C The program reads a system in CSR format, initialize ILU0 preconditioner,
C and solve the system by preconditioned GMRES(IRESTART) method.
C IRESTART is the restart parameter limiting the size of Krylov subspace.
C ======================================================================
      Program  mainGmresILU0
      implicit none

C ======================================================================
      Integer   IRESTART
      Parameter (IRESTART = 20) 
C ======================================================================


C Maximum size of matrix and the maximum number of non-zero entries
      Integer   maxn, maxnz 
      Parameter(maxn = 100 000, maxnz = 1 000 000)

C Arrays for the matrix sparse row (CSR) format
      Integer  ia(maxn+1), ja(maxnz)
      Real*8    a(maxnz), f(maxn), u(maxn)

C Work arrays to keep ILU factors and IRESTART+3 vectors of Krylov subspace 
      Integer   mKryl, MaxWr, MaxWi
      Parameter(mKryl = IRESTART+3,
     &          MaxWr = maxnz + mKryl*maxn, MaxWi = maxnz + 2*maxn+1)
      Real*8    rW(MaxWr)
      Integer   iW(MaxWi)

C GMRES data
      External  matvec, prevec0
      Integer   ITER, INFO, NUNIT
      Real*8    RESID
      Integer   MH, NH
      Parameter ( MH = IRESTART + 1, NH = IRESTART + 6 )
      Real*8    H(MH,NH)


C ILU0 data
      Integer   ierr, ipaLU, ipjLU, ipjU, ipiw

C External routines from BLAS library
      Real*8    ddot
      External  ddot, dcopy    

C Local variables
      Integer   imatvec(1), iprevec(1), ipKrylov
      Real*8    resinit
      Integer   n,i,j,nz

C ======================================================================
C  STAGE 1: read system in CSR format
C ======================================================================
      Open(10,file='../src/Tutorials/PackageILU/CSRsystem',
     &        status='OLD',ERR=1000)
        Read(10,*,ERR=1001) n
        If(n.gt.maxn) Then
           Write(*,*) 'Increase maxn to', n
           Stop 911
        End if

        Read(10,*,ERR=1001) (ia(i), i=1,n+1)
        nz = ia(n+1)-1
        If(nz.gt.maxnz) then
           Write(*,*) 'Increase maxnz to', nz
           Stop 911
        End if

        Read(10,*,ERR=1001) (ja(j), j=1,nz)
        Read(10,*,ERR=1001) (a(j),  j=1,nz)
        Read(10,*,ERR=1001) (f(i),  i=1,n)
      Close(10)


C ======================================================================
C  STAGE 2: initialization of the preconditioner
C ======================================================================
      ipaLU = 1
      ipKrylov = ipaLU+nz
      ipjU  = 1
      ipjLU = ipjU+n+1
      ipiw  = ipjLU+nz !work array of length n

      Call ilu0(n,a,ja,ia, rW(ipaLU),iW(ipjLU),iW(ipjU),iW(ipiw),ierr)

      If(ierr.ne.0) Then
         Write(*,'(A,I7)') 'initialization of ilu0 failed, ierr =', ierr
         Goto 1002
      End if
c Keep data in rW and iW up to rW(nz) and iW(nz+n+1) !


C ======================================================================
C  STAGE 3: set initial guess and compute initial residual
C ======================================================================
c  set initial guess to 0
      Call dcopy(n,0d0,0,u,1)

c  compute initial residual norm
      resinit = ddot(n,f,1,f,1)
      resinit = dsqrt(resinit)
      If(resinit.eq.0d0) Then
         Write(*,'(A)') 'rhs=0, nothing to solve!'
         Goto 1002
      End if

C ======================================================================
C  STAGE 4: iterative solution
C ======================================================================
      ITER = 1000             ! max number of iterations
      RESID = 1d-8 * resinit  ! threshold for \|RESID\|
      INFO  = 0               ! no troubles on input
      NUNIT = 6               ! output channel
      iprevec(1) = n          ! single entry required: system size 
      imatvec(1) = n          ! single entry required: system size 

      call slgmres(
     >     prevec0, iprevec, iW,rW,
     >     matvec,  imatvec, ia,ja,a,
     >     rW(ipKrylov), n, mKryl, H, MH, NH,
     >     n, f, u,
     >     ITER, RESID,
     >     INFO, NUNIT)
      If(INFO.ne.0) Then
         Write(*,'(A)') 'GMRES failed'
         Goto 1002
      End if

      Stop

 1000 Continue
      Write(*,'(A)') 'Cannot open file CSRsystem'
      Stop 911
 1001 Continue
      Write(*,'(A)') 'Corrupted data in CSRsystem'
      Stop 911
 1002 Continue
      Stop 911

      End

