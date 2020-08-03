c-------------------------------------------------------------
c This is the simplest test program which reads a sparce matrix 
c in the CSR format and solves a linear problem with the UMFPACK 
c package. 

c The file 'CSRsystem' has to be created by mainGenFemSys.f prior 
c running this program!
c
c Since UMFPACK uses sorted the CSC format, the matrix is transformed
c from the CSR iformat to the sorted CSC format by routines in 
c aniLU/aux_umf.f.
c 
c Other useful operations may be found in UMFPACK/Demo/umf4hb.fx.
c
c   link:  -liblu.a -lblas 
c   where:  liblu.a = UMFPACK/Demo/umf4_f77wrapper.o + libumfpack.a + libamd.a
c-------------------------------------------------------------
      PROGRAM TESTUMFPACK
      IMPLICIT NONE

c from: umf4_f77wrapper.f
      EXTERNAL  umf4def, umf4sym, umf4num, umf4sol, umf4fsym, umf4fnum
c from: aux_umf.f
      EXTERNAL  IDSORT, TRANSP

      INTEGER   NMAX,  NZMAX
c maximum order of the system
      PARAMETER ( NMAX = 50000 )
c maximum number of non-zero entries in the matrix
      PARAMETER ( NZMAX= NMAX * 7 )
c system data
      integer   IA(NMAX+1),JA(NZMAX), iwk(NZMAX) ! iwk is needed ONLY to transform CSR to CSC
      real*8    A(NZMAX), RHS(NMAX), SOL(NMAX)
c umfpack data
      integer symbolic(2), numeric(2), sys
      real*8  control(20), info(90)
c local
      integer i,j,N,NZ,ierr

c read CSR system
      Open(10,file='../src/Tutorials/PackageLU/CSRsystem',ERR=1001)
        Read(10,*) n,nz
        If(n.GT.NMAX) Then
           Write(*,*)'increase NMAX to', n
           Stop 911
        End if

        If(nz.gt.NZMAX) Then
           Write(*,*) 'increase NZMAX to', nz
           Stop 911
        End if

        Read(10,*) (ia(j),j=1,n+1)
        Read(10,*) (ja(j),j=1,nz)
        Read(10,*) (a(j),j=1,nz)
        Read(10,*) (RHS(j),j=1,n)
      Close(10)

      If(nz.NE.ia(n+1)-1) Then
         Write(*,*) 'Inconsistency of data'
         Stop 911
      End if

c order within each row and transform system to CSC  and make it 0-bazed
      Do i=1,n
         j=ia(i+1)-ia(i)
         Call IDSORT (JA(ia(i)), A(ia(i)), j, 2)
      End do

      Call TRANSP(n, n, a,ja,ia, iwk, ierr)

      If(ierr.NE.0) Then
         Write(*,*) 'transp failed'
         Stop 911
      End if

      Do i = 1, n+1
         IA(i) = IA(i) - 1
      End do
      Do j = 1, nz
         JA(j) = JA(j) - 1
      End do
c now the matrix is in 0-bazed CSC format with sorted row indexes.

c set the default control parameters in the Control array
      Call umf4def( control )
c print error messages only
      control(1) = 1

c pre-order and symbolic analysis
      Call umf4sym( n,n, ia,ja,a, symbolic,control,info)

      If(info(1).LT.0) Then
         Write(*,*) 'Error occurred in umf4sym: ', info(1)
         Stop 911
      End if
      Write(*,*) 'time of symbolic analysis:', info(16)

c numeric factorization
      Call umf4num(ia,ja,a, symbolic,numeric,control,info)

      If(info(1) .LT. 0) Then
         Write(*,*) 'Error occurred in umf4num: ', info(1)
         Stop 911
      End if

c print statistics for the numeric factorization
c call umf4pinf (control, info) could also be done.
      Write(*,90) info (66),
     $    (info (41) * info (4)) / 2**20,
     $    (info (42) * info (4)) / 2**20,
     $    info (43), info (44), info (45)
90    Format ('numeric factorization:',/,
     $    '   time:    ', e10.2, /,
     $    '   actual numeric LU statistics:', /,
     $    '   size of LU:    ', f10.2, ' (MB)', /,
     $    '   memory needed: ', f10.2, ' (MB)', /,
     $    '   flop count:    ', e10.2, /
     $    '   nnz (L):       ', f10.0, /
     $    '   nnz (U):       ', f10.0)

c free the symbolic analysis data
      Call umf4fsym (symbolic)

c solve Ax=b, without iterative refinement
      sys = 0
      Call umf4sol (sys, SOL, RHS, numeric, control, info)

      If(info(1) .LT. 0) Then
         Write(*,*) 'Error occurred in umf4sol: ', info(1)
         Stop 911
      End if


c free the numeric factorization data
      Call umf4fnum (numeric)

      Stop
 
 1001 Continue
      Write(*,*) 'no file CSRsystem'
      Stop 911

      End

      
