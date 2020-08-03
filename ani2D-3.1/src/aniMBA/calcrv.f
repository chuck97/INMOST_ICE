C ======================================================================
      Subroutine aniCrv(tc, XYP, iFnc, CrvFunction)
C ======================================================================
C The function is used as a buffer between the real function and
C the code. It rescales a physical point to the computational domain.
C ======================================================================
      Real*8   tc, XYP(2)
      Integer  iFnc
      EXTERNAL CrvFunction

      Real*8  refXYP(2), scaXYP(2)
      Integer isON
      Common /rescale/refXYP, scaXYP, isON

C ======================================================================
      Call CrvFunction(tc, XYP, iFnc)

      If(isON.EQ.1) Then
         Do i = 1, 2
            XYP(i) = (XYP(i) - refXYP(i)) * scaXYP(i)
         End do
      End if

      Return
      End



C ======================================================================
      Subroutine scale2Square(nP, nF, XYP, lbC, flag)
C ======================================================================
C Routine scales the model to the square [0.1, 0.9]^2. We allow
C 10% freedom for curved edges. No scaling will be performed if
C there are no curved faces and the mesh is already inside the
C unit square.
C ======================================================================
      Real*8  XYP(2, *)
      Integer lbC(*)
      Logical flag

C ======================================================================
      Real*8  minXYP(2), maxXYP(2), scale, size

      Real*8  refXYP(2), scaXYP(2)
      Integer isON
      Common /rescale/refXYP, scaXYP, isON

C ======================================================================
      If(flag) Then
         Do i = 1, 2
            minXYP(i) = XYP(i, 1)
            maxXYP(i) = XYP(i, 1)
         End do

         Do n = 2, nP
            Do i = 1, 2
               minXYP(i) = min(minXYP(i), XYP(i, n))
               maxXYP(i) = max(maxXYP(i), XYP(i, n))
            End do
         End do

c  ...  add 10% for the bounding box
         isON = 0
         Do i = 1, 2
            If(minXYP(i).LT.0D0 .OR. maxXYP(i).GT.1D0) isON = 1
         End do

         If(isON.EQ.0) Then
            Do n = 1, nF
               If(lbC(n).GT.0) Then
                  isON = 1
                  Goto 100
               End if
            End do
         End if

 100     Continue
         If(isON.EQ.1) Then
            Do i = 1, 2
               size = (maxXYP(i) - minXYP(i)) / 10
               minXYP(i) = minXYP(i) - size
               maxXYP(i) = maxXYP(i) + size
            End do

            Do i = 1, 2
               refXYP(i) = minXYP(i)
               scaXYP(i) = 1D0 / (maxXYP(i) - minXYP(i))
            End do

            scale = min(scaXYP(1), scaXYP(2))

            Do i = 1, 2
               scaXYP(i) = scale
            End do

            Do n = 1, nP
               Do i = 1, 2
                  XYP(i, n) = (XYP(i, n) - refXYP(i)) * scaXYP(i)
               End do
            End do
         End if
      Else If(isON.EQ.1) Then 
         Do i = 1, 2
            scaXYP(i) = 1D0 / scaXYP(i)
         End do

         Do n = 1, nP
            Do i = 1, 2
               XYP(i, n) = refXYP(i) + XYP(i, n) * scaXYP(i)
            End do
         End do
      End if

      Return
      End

 

C ======================================================================
      Subroutine scaleBack(XYPi, XYPo)
C ======================================================================
C  Routine computes physical coordinates of point XYPi 
C ======================================================================
      Real*8  XYPi(2), XYPo(2)

      Real*8  refXYP(2), scaXYP(2)
      Integer isON
      Common /rescale/refXYP, scaXYP, isON

C ======================================================================
      If(isON.EQ.1) Then
         Do i = 1, 2
            XYPo(i) = refXYP(i) + XYPi(i) / scaXYP(i)
         End do
      Else
         Do i = 1, 2
            XYPo(i) = XYPi(i)
         End do
      End if

      Return
      End



C ======================================================================
      Subroutine scale2Meshes(nP, XYP, nP2, XYP2, flag)
C ======================================================================
C Routine applies scale2Square() algorithm for flag = .TRUE.
C ======================================================================
      Real*8  XYP(2, *), XYP2(2, *)
      Logical flag

C ======================================================================
      Real*8  minXYP(2), maxXYP(2)

      Real*8  refXYP(2), scaXYP(2)
      Integer isON
      Common /rescale/refXYP, scaXYP, isON

C ======================================================================
      If(flag) Then
         Do i = 1, 2
            minXYP(i) = XYP(i, 1)
            maxXYP(i) = XYP(i, 1)
         End do

         Do n = 2, nP
            Do i = 1, 2
               minXYP(i) = min(minXYP(i), XYP(i, n))
               maxXYP(i) = max(maxXYP(i), XYP(i, n))
            End do
         End do

         Do n = 1, nP2
            Do i = 1, 2
               minXYP(i) = min(minXYP(i), XYP2(i, n))
               maxXYP(i) = max(maxXYP(i), XYP2(i, n))
            End do
         End do

c  ...  add 10% for the bounding box
         isON = 0
         Do i = 1, 2
            If(minXYP(i).LT.0D0 .OR. maxXYP(i).GT.1D0) isON = 1
         End do

         If(isON.EQ.1) Then
            Do i = 1, 2
               refXYP(i) = minXYP(i)
               scaXYP(i) = 1D0 / (maxXYP(i) - minXYP(i))
            End do

            Do n = 1, nP
               Do i = 1, 2
                  XYP(i, n) = (XYP(i, n) - refXYP(i)) * scaXYP(i)
               End do
            End do

            Do n = 1, nP2
               Do i = 1, 2
                  XYP2(i, n) = (XYP2(i, n) - refXYP(i)) * scaXYP(i)
               End do
            End do
         End if

      Else If(isON.EQ.1) Then 
         Do n = 1, nP
            Do i = 1, 2
               XYP(i, n) = refXYP(i) + XYP(i, n) / scaXYP(i)
            End do
         End do

         Do n = 1, nP2
            Do i = 1, 2
               XYP2(i, n) = refXYP(i) + XYP2(i, n) / scaXYP(i)
            End do
         End do
      End if

      Return
      End

 





