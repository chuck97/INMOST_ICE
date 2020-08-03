C ======================================================================
C The program reads a system in CSR format, initialize ILU2 preconditioner,
C and solve the system by preconditioned BiCGstab method
C ======================================================================
      Program  mainILU2
      implicit none

C Maximum size of matrix and the maximum number of non-zero entries
      Integer   maxn, maxnz 
      Parameter(maxn = 100 000, maxnz = 1 000 000)

C Arrays for the matrix in sparse row (CSR) format
      Integer  ia(maxn+1), ja(maxnz)
      Real*8    a(maxnz), f(maxn), u(maxn)

C Work arrays
      Integer   MaxWr, MaxWi
      Parameter(MaxWr = 5*maxnz, MaxWi = 6*maxnz)
      Real*8    rW(MaxWr)
      Integer   iW(MaxWi)

C BiCGStab data
      External  matvec, prevec2
      Integer   ITER, INFO, NUNIT
      Real*8    RESID

C ILU data
      Real*8    tau1,tau2,partlur,partlurout
      Integer   verb, ierr, UsedWr, UsedWi

C External routines from the BLAS library
      Real*8    ddot
      External  ddot, dcopy  

C Local variables
      Integer   n,i,j,nz, ipBCG
      Integer   imatvec(1), iprevec(1)
      Real*8    resinit

C ======================================================================
C  STAGE 1: read system in CSR format
C ======================================================================
      Open(10,file='../src/Tutorials/PackageILU/CSRsystem',
     &        status='OLD',ERR=1000)
        Read(10,*,ERR=1001) n
        If(n.gt.maxn) then
           Write(*,'(A,I6)') 'Increase maxn to', n
           Stop 911
        End if

        Read(10,*,ERR=1001) (ia(i), i=1,n+1)
        nz = ia(n+1)-1
        If(nz.gt.maxnz) then
           Write(*,'(A,I6)') 'Increase maxnz to', nz
           Stop 911
        End if

        Read(10,*,ERR=1001) (ja(j), j=1,nz)
        Read(10,*,ERR=1001) (a(j),  j=1,nz)
        Read(10,*,ERR=1001) (f(i),  i=1,n)
      Close(10)


C ======================================================================
C  STAGE 2: initialization of the preconditioner
C ======================================================================
      verb = 0 ! verbose no
      tau1 = 1d-2
      tau2 = 1d-3
      partlur = 0.5
      ierr = 0
 
      Call iluoo(n, ia, ja, a, tau1, tau2, verb,
     &           rW, iW, MaxWr, MaxWi, partlur, partlurout,
     &           UsedWr, UsedWi, ierr)
      If(ierr.ne.0) Then
         Write(*,'(A,I7)') 'Initialization of iluoo failed, ierr=', ierr
         Stop 911
      End if
      Write(*,'(A)') ' Recommended memory ILU2 + BiCGs:'
      Write(*,'(A,F12.2)') ' partlur =', partlurout
      Write(*,'(A,I12)') ' MaxWr   =', UsedWr+8*n
      Write(*,'(A,I12)') ' MaxWi   =', UsedWi

      If(UsedWr+8*n.gt.MaxWr) Then
         Write(*,'(A,I7)') 'Increase MaxWr to ', UsedWr+8*n
         Stop 911
      End if
      ipBCG = UsedWr + 1


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
         Stop 911
      End if


C ======================================================================
C  STAGE 4: iterative solution
C ======================================================================
      ITER = 1000             ! max number of iterations
      RESID = 1d-8 * resinit  ! threshold for \|RESID\|
      INFO  = 0               ! no troubles on imput
      NUNIT = 6               ! output channel
      iprevec(1) = n          ! single entry required: system size 
      imatvec(1) = n          ! single entry required: system size 

      call slpbcgs(
     >     prevec2, iprevec, iW,rW,
     >     matvec,  imatvec, ia,ja,a,
     >     rW(ipBCG), n, 8,
     >     n, f, u,
     >     ITER, RESID,
     >     INFO, NUNIT)

      If(INFO.ne.0) Then
         Write(*,'(A)') 'BiCGStab failed'
         Stop 911
      End if

      Stop

 1000 Continue
      Write(*,'(A)') 'Cannot open file CSRsystem'
      Stop 911
 1001 Continue
      Write(*,'(A)') 'Corrupted data in CSRsystem'
      Stop 911

      End

