C ======================================================================
      Subroutine mulAcsr(N, IA, JA, A, U, AU)
C ======================================================================
c  Computing  AU = A * U, where A is in compressed sparse row format
C ======================================================================
      implicit none

      Integer  N, IA(*), JA(*)
      Real*8   A(*), U(*), AU(*)

      Integer  i, k
C ======================================================================
      Do k = 1, N
         AU(k) = 0D0
      End do

      Do k = 1, N
         Do i = IA(k), IA(k + 1) - 1
            AU(k) = AU(k) + A(i) * U(JA(i))
         End do
      End do
      Return
      End



C ======================================================================
      Subroutine mulAcsc(N, IA, JA, A, U, AU)
C ======================================================================
c  Computing  AU = A * U, where A is in compressed sparse column format
C ======================================================================
      implicit none
      Integer N, IA(*), JA(*)
      Real*8  A(*), U(*), AU(*)

      Integer i, j, k
C ======================================================================
      Do k = 1, N
         AU(k) = 0D0
      End do

      Do k = 1, N
         Do i = IA(k), IA(k + 1) - 1
            j = JA(i) 
            AU(j) = AU(j) + A(i) * U(k)
         End do
      End do
      Return
      End



C ======================================================================
      Subroutine mulAcsc0(N, IA, JA, A, U, AU)
C ======================================================================
c  Computing  AU = A * U, where A is in compressed sparse zero-based 
c  column format
C ======================================================================
      implicit none
      Integer N, IA(*), JA(*)
      Real*8  A(*), U(*), AU(*)

      Integer i, j, k
C ======================================================================
      Do k = 1, N
         AU(k) = 0D0
      End do

      Do k = 1, N
         Do i = IA(k) + 1, IA(k + 1)
            j = JA(i) + 1
            AU(j) = AU(j) + A(i) * U(k)
         End do
      End do
      Return
      End



C ======================================================================
      Subroutine CSR2CSC(NA, IA, JA, A, NB, IB, JB, B)
C ======================================================================
c  Calculating the transpose matrix in AMG format. Basically, we calculate 
c  the column based sparse format.
C
C  *** Note: 1. size(IB) = max(JA) + 1
C            2. size(JB) = size(JA)
C            3. size(B) = size(A) 
C ======================================================================
      Integer NA, IA(*), JA(*)
      Integer NB, IB(*), JB(*)
      Real*8  A(*), B(*)

      M = IA(NA + 1) - 1

c ... compute the number of rows
      NB = 0
      Do i = 1, M
         NB = max(NB, JA(i))
      End do


c ... compute the number of non-zero elements in each row of B
      Do i = 1, NB + 1
         IB(i) = 0
      End do

      Do i = 1, M
         j = JA(i)
         IB(j) = IB(j) + 1 
      End do 

      Do i = 1, NB 
         IB(i + 1) = IB(i + 1) + IB(i)
      End do

c ... copy entries of matrix A to matrix B
      Do i = NA, 1, -1
         Do j = IA(i), IA(i + 1) - 1
            k = JA(j)

            l = IB(k)
            IB(k) = l - 1

            JB(l) = i
             B(l) = A(j)
         End do
      End do

      Do i = 1, NB + 1
         IB(i) = IB(i) + 1
      End do

      Return
      End



C ======================================================================
      Subroutine AMG2CSC(NA, IA, JA, A, NB, IB, JB, B)
C ======================================================================
c  Converting the AMG format (CSR) to the compressed 0-based 
C  column format (CSC).
C
C  *** Note: 1. size(IB) = max(JA) + 1
C            2. size(JB) = size(JA)
C            3. size(B) = size(A) 
C ======================================================================
      Integer NA, IA(*), JA(*)
      Integer NB, IB(*), JB(*)
      Real*8  A(*), B(*)

      M = IA(NA + 1) - 1

c ... compute the number of rows
      NB = 0
      Do i = 1, M
         NB = max(NB, JA(i))
      End do


c ... compute the number of non-zero elements in each row
      Do i = 1, NB + 1
         IB(i) = 0
      End do

      Do i = 1, M
         j = JA(i)
         IB(j) = IB(j) + 1 
      End do 

      Do i = 1, NB
         IB(i + 1) = IB(i + 1) + IB(i)
      End do


c ... copy matrix A to matrix B
      Do i = NA, 1, -1
         Do j = IA(i), IA(i + 1) - 1
            k = JA(j)

            l = IB(k)
            IB(k) = l - 1

            JB(l) = i
             B(l) = A(j)
         End do
      End do

      Do i = 1, NB + 1
         IB(i) = IB(i) + 1
      End do

      Return
      End



C ======================================================================
      Subroutine CSC2AMG(NA, IA, JA, A, NB, IB, JB, B)
C ======================================================================
c  Converting the compressed 0-based column format (CSC) to the
C  AMG (row-wice) format.
C
C  *** Note: 1. size(IB) = max(JA) + 1
C            2. size(JB) = size(JA)
C            3. size(B) = size(A) 
C ======================================================================
      Integer NA, IA(*), JA(*)
      Integer NB, IB(*), JB(*)
      Real*8  A(*), B(*)

      Call AMG2CSC(NA, IA, JA, A, NB, IB, JB, B)

      Do 10 i = 1, NB
         j1 = IB(i)
         j2 = IB(i + 1) - 1

         Do j = j1, j2
            k = JB(j)
            If(k.EQ.i) Then
               Call swapii(JB(j), JB(j1)) 
               Call swapdd( B(j),  B(j1))
               Goto 10
            End if
         End do
 10   Continue

      Return
      End



C ======================================================================
      Subroutine CSC2CSC0(nB, IB, JB)
C ======================================================================
      Integer nB, IB(*), JB(*)

      Do i = 1, nB + 1
         IB(i) = IB(i) - 1
      End do

      Do i = 1, IB(nB + 1)
         JB(i) = JB(i) - 1
      End do

      Return
      End



C ======================================================================
      Subroutine mulABgen(nA, IA, JA, A, 
     &                    nB, IB, JB, B, 
     &                    nC, IC, JC, C, MaxC, iW)
C ======================================================================
c  Computing  C = A * B where all matrices are in the sparse row format
c  The size of working array is the number of columns of B. The AMG
c  format is enforced for matrix C.
C ======================================================================
      Integer IA(*), JA(*), IB(*), JB(*), IC(*), JC(*)
      Real*8  A(*), B(*), C(*)

      Integer iW(*)
C ======================================================================
      nC = nA

      mC = 0
      iC(1) = 1

c ... the number of columns of B
      nColB = 0
      Do n = 1, IB(nB+1) - 1
         nColB = max(nColB, JB(n)) 
      End do

      Do n = 1, nColB
         iW(n) = 0
      End do
      
      Do i = 1, nA
         If(mC + nColB.GT.MaxC) Call errMesFEM(1004, 
     &        'algebra', 'Not enough memory for the matrix')

         If(i.LE.nColB) Then
            mC = mC + 1
            JC(mC) = i
            C(mC)  = 0D0
            iW(i)  = mC
         End if

         Do n = IA(i), IA(i + 1) - 1
            k = JA(n)
            Do l = IB(k), IB(k + 1) - 1
               j = JB(l)
               m = iW(j)

               If(m.NE.0) Then 
                  C(m) = C(m) + A(n) * B(l)
               Else
                  mC = mC +  1
                  JC(mC) = j
                  C(mC)  = A(n) * B(l)
                  iW(j)  = mC
               End if
            End do
         End do
         IC(i + 1) = mC + 1

         Do k = IC(i), mC 
            iW(JC(k)) = 0
         End do
      End do

      Return
      End



C ======================================================================
      Subroutine addBgen(nA, IA, JA, A, nB, IB, JB, B)
C ======================================================================
c  Computing  A = A + B. The sparcity structure of A must be
c  bigger than that of B.
C ======================================================================
      Integer IA(*), JA(*), IB(*), JB(*)
      Real*8  A(*), B(*)

C ======================================================================
      Do i = 1, nB
         k1 = IA(i)   
         k2 = IA(i + 1) + 1

         j1 = IB(i)
         j2 = IB(i + 1) - 1

         Do j = j1, j2
            Call findSE(k2 - k1, JA(k1), JB(j), k)

            If(k.EQ.0) Call errMesFEM(2011, 'algebra',
     &                 'Structure of B bigger than structure of A.')

            k = k1 + k - 1
            A(k) = A(k) + B(j)
         End do
      End do
 
      Return
      End



C ======================================================================
      Real*8 Function DotProd(N, A, B)
C ======================================================================
C Routine computes scalar product of two vectors
C ======================================================================
      Real*8 A(*), B(*)

      DotProd = 0D0
      Do i = 1, N
         DotProd = DotProd + A(i) * B(i)
      End do

      Return
      End



C ======================================================================
      Subroutine diagonal(nA, IA, JA, A, DA)
C ======================================================================
C Routines retrs diagonal of martix A.
C ======================================================================
      Integer IA(*), JA(*)
      Real*8   A(*), DA(*)

C ======================================================================
      Do 100 i = 1, nA
         Do j = IA(i), IA(i+1) - 1
            k = JA(j)
            If(k.EQ.i) Then
               DA(i) = A(k)
               Goto 100
            End if
         End do
 100  Continue

      Return
      End



C ======================================================================
      Subroutine gershgorin(N, IA, JA, A, center, radius)
C ======================================================================
C  Calculates the center and radius of the biggest Gershgorin circle.
C ======================================================================
      implicit none

      Integer  N, IA(*), JA(*)
      Real*8   A(*), center, radius

      Integer  i, k
      Real*8   sc, sr
C ======================================================================
      radius = 0D0

      Do k = 1, N
         sr = 0D0
         Do i = IA(k), IA(k + 1) - 1
            sr = sr + dabs(A(i))
            If(k.EQ.JA(i)) sc = A(i)
         End do

         If(radius.LT.sr) Then
            radius = sr
            center = sc
         End if
      End do
      Return
      End



C ======================================================================
      Subroutine minmax(A, N1, N2)
C ======================================================================
      Integer N1, N2
      Real*8  A(*), amin, amax, aavg
C ======================================================================
      amin = A(N1)
      amax = A(N1)
      aavg = A(N1)

      Do i = N1 + 1, N2
         amin = min(amin, A(i))
         amax = max(amax, A(i))
         aavg = aavg + A(i)
      End do

      aavg = aavg / (N2 - N1 + 1)

      Write(*,5000) N1, N2, amin, amax, aavg

 5000 Format('A[',I5,':',I5,']  Bounds:', 3E14.6)

      Return
      End


C ======================================================================
      Subroutine printA(nA, IA, JA, A)
C ======================================================================
      Integer IA(*), JA(*)
      Real*8  A(*)

C ======================================================================
      Do i = 1, nA
         j1 = IA(i)   
         j2 = IA(i + 1) - 1

         Write(*, '(A,I2,A,100I6)') 'Row[', i, ']:', (JA(j), j=j1,j2)
         Write(*, '(A,100F7.3)')    '         ', (A(j), j=j1,j2)
      End do
 
      Return
      End



C ======================================================================
      Subroutine printAfull(nA, IA, JA, A, rW)
C ======================================================================
      Integer IA(*), JA(*)
      Real*8  A(*), rW(*)

C ======================================================================
c ... the number of columns of A
      nColA = 0
      Do n = 1, IA(nA+1) - 1
         nColA = max(nColA, JA(n)) 
      End do

      Do i = 1, nA
         j1 = IA(i)   
         j2 = IA(i + 1) - 1

         Do j = 1, nColA
            rW(j) = 0D0
         End do

         Do j = j1, j2 
            rW(JA(j)) = A(j)
         End do 

c        Write(*, n'(3000F11.5)') (rW(j), j=1,nColA)
         Write(*, *) (rW(j), j=1,nColA)
      End do
 
      Return
      End

