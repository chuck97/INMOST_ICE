C ======================================================================
      Subroutine mbaNodal(
C ======================================================================
c group (M)
     &      nP, nPfix, MaxP, XYP, lbP, fixP,
     &      nF, nFfix, MaxF, IPF, lbF, fixF,
     &      nC,              Crv, lbC, CrvFunction,
     &      nE, nEfix, MaxE, IPE, lbE, fixE,
c group (CONTROL)
     &      nEStar, Quality, control, Metric,
c group (W)
     &      MaxWr, MaxWi, rW, iW, iERR)
C ======================================================================
      implicit none
      include 'lintrp.fd'
      include 'status.fd'
C ======================================================================
C Routine adapts the mesh to the discrete metric given at mesh nodes.
C ======================================================================
C ======= MESH DESCRIPTION =======
C  
C     nP    - the number of points (P)
C     nPfix - the number of fixed points
C     MaxP  - the maximal number of points
C
C     XYP(2, MaxP) - the Cartesian coordinates of mesh points
C     lbP(MaxP)    - point indentificator (a non-negative number)
C     fixP(MaxP)   - list of fixed points
C
C     nF    - the number of boundary and interface edges (F)
C     nFfix - the number of fixed edges
C     MaxF  - the maximal number of edges
C
C     IPF(2, MaxF) - connectivity list of boundary and interface edges
C     lbF(MaxF)    - boundary edge identificator (Example: unit square 
C           has 4 boundaries which may have or not different ids. To
C           automatically recognize corners of the square, boundaries 
C           must have different colors. This is not required if the 
C           corner points are in the list of fixed points. 
C           See colors.fd for more detail.
C
C     fixF(MaxF) - list of fixed edges
C
C     nC - the number of curved edges
C
C     Crv(2, MaxF) - parametrizations of curvilinear edges 
C                       column 1 - parameter for the starting point
C                       column 2 - parameter for the terminal point
C
C           Parameters for the inner points are computed by linear 
C           interpolation between two given numbers. The Cartesian 
C           coordinates are computed using the user-given formulas 
C           in routine CrvFunction.
C
C     lbC(MaxF) - zero or positive function number for computing the 
C                 Cartesian coordinates
C
C     CrvFunction - user-created routine that computes the Cartesian 
C           coordinates of point xyc from its parametric coordinate tc:
C
C           Subroutine CrvFunction(tc, xyc, lbC)
C              tc     - parametric coordinate of point xyc
C              xyc(2) - Cartesian coordinate of the same point
C              lbC    - the function number associated with curved edge
C     
C     nE    - the number of triangles (E)
C     nEfix - the number of fixed triangles
C     MaxE  - the maximal number of triangles
C
C     IPE(3, MaxE) - connectivity list of triangles
C     lbE(MaxE)    - triangle indentificator (a positive number)
C     fixE(MaxE)   - list of fixed triangles
C
C ======= CONTROL ========
C
C     nEstar  - the desired number of triangles
C  
C     Quality - On input, the requested quality for the final grid
C               (a positive number between 0 and 1)
C               On output, the reached quality for the final grid
C
C     control(6) - integer array with control parameters:
C         MaxSkipE = control(1)
C         MaxQItr  = control(2)
C         status   = control(3)
C         flagAuto = control(4).GT.0
C         iPrint   = control(5)
C         iErrMesg = control(6)
C
C     MaxSkipE  - the maximal number of skipped triangles
C     MaxQItr   - the maximal number of local grid modifications
C         status   - advanced control of mesh generation (see below)
C         flagAuto - flag controling the mesh generation:
C                    TRUE  - automatic recover of missing mesh elements
C                    FALSE - rigorous checking of input data 
C         iPrint   - the level of output information (between 0 and 9)
C         iErrMesg - flag controling interior code termination.
C                    0 - unrecoverable code error
C                    * - user will decided based on the error message  
C
C     Metric(3, nP) - Real*8 array containing the metric defined at
C             mesh points. The metric is a 2x2 positive definite 
C             symmetric tensor:
C
C                            M11   M12   
C                   Metric = 
C                            M12   M22
C
C             Each column of this array stores the upper triangular 
C             entries in the following order: M11, M22, and M12.
C
C ======= WORKING MEMORY ========
C
C     MaxWr     - the maximal space for real*8 arrays
C     MaxWi     - the maximal space for integer arrays
C
C     rW(MaxWr) - the real*8  working array. More details are below.
C     iW(MaxWi) - the integer working array. More details are below.
C     iERR - the error code:
C                  0 - correct termination
C                 >0 - error message as described in error.f 
C
C ======================================================================
C Note about status = control(3)
C     0 or negative - no additional control
C     positive      - enforce some of the mesh properties. The detailed 
C          description of available properties is in file status.fd. 
C          Variable status is equal to the sum of positive numbers 
C          corresponding to the desired properties. Here is the list of 
C          user-controled properties:
C            1  - the final mesh will not have boundary elements 
C            4  - the algorithm will not change boundary edges 
C            8  - missing material interfaces and boundary edges created 
C                 by the code will be removed from the final mesh
C            16 - the algorithm will not change boundary points.
C
C          Example: If only first two properties are required, set 
C              status = 5. Do not use the number since status is the 
C              input/output parameter.
C
C ======================================================================
C Note about working arrays. They return some additional statistics.
C  rW(1)= hStar - average size of triangles (with respect to the given 
C                 metric)
C
C  iW(1) - the minimal required integer memory allocation
C  iW(2) - the minimal required real*8  memory allocation
C  iW(3) - the number of performed iterations
C  iW(4:nP+3) - colors of mesh points as described in file color.fd
C
C ======================================================================
C This package uses non-standard notations for the mesh. These notation
C are more convient for various additional connectivity arrays. This
C table shows correspondance between our notations and the standard one:
C 
C  nP = nv,   nPfix = nvfix   and so on...
C  XYP(2,*) = vrt(2,*) 
C
C  nF = nb,   nFfix = nbfix   and so on...
C  IPF(2,*) = bnd(2,*)
C
C  nE = nt,   nEfix = ntfix   and so on...
C  IPE(3,*) = tri(3,*)   
C ======================================================================
C  Note:
C       Input parameters:  MaxP, MaxF, MaxE, nPv,
C                          IPV, IFV, IEV, lbE, flagAuto, 
C                          nEStar, MaxSkipE, MaxQItr,
C                          Metric, Quality, MaxWr, MaxWi, iPrint
C
C       Input / output:    nP, nF, nE, XYP, IPF, IPE,
C                          ParCrv, iFnc, Sol, status, rW, iW
C
C       Output parameters: rQuality, iERR
C
C ======================================================================
C Note on a possible choice for input parameters. 
C       MaxP > nP
C       MaxF > nF
C       MaxE > nE
C
C       nPfix >= 0
C       nFfix >= 0
C       nEfix >= 0
C
C       MaxSkipE = 100
C       MaxQItr = 15000
C       Quality = 0.6
C
C       MaxWr > 7 * MaxP + 10 * nP + MaxF + MaxE (approximate estimate)
C       MaxWi > 6 * MaxP + 10 * nP + 19 * MaxF + 11 * MaxE + 12 * nE
C       (approximately)
C
C ======================================================================
C *** Authors: K. Lipnikov (lipnikov@gmail.com)
C              Y. Vassilevski (yuri.vasilevski@gmail.com)
C ======================================================================
C group (M)
      Integer  nP, nPfix, MaxP, lbP(*), fixP(*)
      Real*8  XYP(2, *)

      Integer  nF, nFfix, MaxF, IPF(2, *), lbF(*), fixF(*)

      Integer  nC, lbC(*)
      Real*8   Crv(2, *)
      EXTERNAL CrvFunction

      Integer  nE, nEfix, MaxE, IPE(3, *), lbE(*), fixE(*)

C group (CONTROL)
      Integer  nEStar, control(*)
      Real*8   Metric(3, *), Quality

C group (W)
      Integer  MaxWr, MaxWi, iW(*), iERR
      Real*8  rW(*)


C LOCAL VARIABLEs
      Integer  ANI_MetricFunction
      EXTERNAL ANI_MetricFunction

      Integer  MaxSkipE, MaxQItr, status, nQItr, iPrint, iErrMesg
      Integer  flagFixShape
      Real*8   hStar, rQuality
      Logical  flagAnalytic, flagAuto

      Integer  iICP, iIEP, iIFE, iIEE, iXYPw, iIPEw, iIEPw, inEPw
      Integer  iIHolP, iIHolF, iIHolE
      Integer  iL1E, iL2E, iL1Et, iL2Et, inL2t, inStept
      Integer  iqE, iLFnc, iILt, inEt, itE, irSE, iiSE, iiW
      Integer  iHesP, iHesPw, idG

      Integer  i, k, n, nWr, nWi, mrLINTRP, miLINTRP

C ======================================================================
C group (Common blocks)
      Integer iDomBnd, iMatBnd
      Common /aniBND/ iDomBnd, iMatBnd
 
      Real*8  refXYP(2), scaXYP(2)
      Integer isON
      Common /rescale/refXYP, scaXYP, isON

C ======================================================================
      iERR = 0

c ... unpack control
      MaxSkipE = control(1)
      MaxQItr  = control(2)
      status   = control(3)
      flagAuto = control(4).GT.0
      iPrint   = control(5)
      iErrMesg = control(6)

      If(MaxSkipE.LE.0) MaxSkipE = max(100,  nE / 20)
      If(MaxQItr .LE.0) MaxQItr  = max(1000, nE * 3)


c ... default scaling (no scaling)
      isON = 0
      Do i = 1, 2
         refXYP(i) = 0D0
         scaXYP(i) = 1D0
      End do


c ... copy the input metric to working file
      If(3*MaxP.GT.MaxWr) Then
         iERR = 1002
         Goto 1000
      End if
      
      iHesP = 1
      k = 0
      Do n = 1, nP
         Do i = 1, 3
            k = k + 1
            rW(k) = Metric(i, n)
         End do
      End do


c ... refine initial mesh when nE is very small
c ... it increases robustness of the code 
      Do while(nE.LE.nEStar / 15 .AND. nE.LE.500 .AND. nEfix+nFfix.EQ.0)
         iIFE = 1
         iiW  = iIFE + 3 * nE
         nWi  = iiW  + 3 * nE + nP 
         If(nWi.GT.MaxWi) Goto 100

         If(iPrint.GE.1) Write(*,5001) nP, nE

         Call uniformRefinement(
     &        nP, MaxP, nF, MaxF, nE, MaxE,
     &        XYP, IPF, lbF, IPE, lbE,
     &        CrvFunction, Crv, lbC, iW(iIFE),
     &        rW(iHesP), 3, iW(iiW), MaxWi)
      End do 


 100  miLINTRP = 10 * nP + 3 * nE + 6
      mrLINTRP =  4 * nP + MaxH + 4

      inEt = 1
      inStept = inEt + MaxF
      inL2t = inStept + 4 * MaxF
      iLFnc = inL2t + MaxF
      iILt  = iLFnc + MaxF
      iL1Et = iILt + MaxF
      iL2Et = iL1Et + 2 * MaxF
      iIHolP = iL2Et + 2 * MaxF
      iIHolF = iIHolP + MaxP
      iIHolE = iIHolF + MaxF
      iICP = iIHolE + MaxE
      iIEP = iICP + MaxP
      iIFE = iIEP + MaxP
      iIEE = iIFE + 3 * MaxE
      iL1E = iIEE + 3 * MaxE
      iL2E = iL1E + 2 * MaxE
      iIPEw = iL2E + 2 * MaxE
      iiSE  = iIPEw + 3 * nE
      iIEPw = iiSE + miLINTRP
c ... we need twice less memory for backReferences
      inEPw = iIEPw + max(6 * nE, 4 * MaxF)
      nWi   = inEPw + max(3 * MaxP, 2 * MaxF)


      itE = iHesP + 3 * MaxP
      idG = itE + MaxF
      iqE = idG + MaxP
      iHesPw = iqE + MaxE
      iXYPw = iHesPw + 3 * nP
      irSE = iXYPw + 2 * nP
      nWr  = irSE + max(mrLINTRP, max(nE, MaxF))


      iW(1) = nWi
      iW(2) = nWr
      If(nWi.GT.MaxWi) Then
         iERR = 1001
         Goto 1000
      End if

      If(nWr.GT.MaxWr) Then
         iERR = 1002
         Goto 1000
      End if


      Do n = itE, nWr
         rW(n) = 0D0
      End do

      Do n = 1, nWi
         iW(n) = 0
      End do


c ... scale geometry to unit cube
      Call scale2Square(nP, nF, XYP, lbC, .TRUE.)


c ... print Ani2D header
      If(iPrint.GE.1) Write(*, 5004) Quality, nEStar, MaxSkipE, MaxQItr


c ... set up default status
      Call setStatus(flagAuto, status, iPrint)


c ... call the main module
      flagAnalytic = .FALSE.
      flagFixShape = 0

      Call ani2(
c group (M)
     &      nP, nPfix, MaxP, XYP, lbP, fixP,
     &      nF, nFfix, MaxF, IPF, lbF, fixF,
     &      nC,              Crv, lbC, CrvFunction,
     &      nE, nEfix, MaxE, IPE, lbE, fixE,
c group (M-EXT)
     &      iW(iICP), iW(iIEP), iW(iIFE), iW(iIEE),
     &      iW(iIHolP), iW(iIHolF), iW(iIHolE),
     &      rW(iXYPw), iW(iIPEw), rW(iHesPw), 
     &      iW(iIEPw), iW(inEPw), iW(iiSE), rW(irSE),
c group (CRV)
     &      iW(iL1Et), iW(iL2Et), rW(itE),
     &      iW(inL2t), iW(inStept), iW(inEt),
     &      iW(iLFnc), iW(iILt),
c group (Q)
     &      nEStar, hStar, 
     &      Quality, rQuality, rW(iHesP), rW(idG),
     &      iW(iL1E), iW(iL2E), rW(iqE),
c group (CONTROL)
     &      flagFixShape, flagAuto, status,
     &      ANI_MetricFunction, flagAnalytic,
     &      MaxSkipE, MaxQItr, nQItr,
     &      iPrint, iERR)


c ... rescale geometry back
      Call scale2Square(nP, nF, XYP, lbC, .FALSE.)


c ... returning sadditional information
      iW(3) = nQItr
      Do n = 1, nP
         iW(n + 3) = iW(iICP + n - 1)
      End do

      rW(1) = hStar
      Quality = rQuality


 1000 If(iERR.EQ.0 .OR. iERR.EQ.1000) Return
      Call errMesMBA(iERR, 'mbaNodal', 'See error.f for more detail')

      Return

 5001 Format('     Auto mesh refinement:', I6, ' points and', 
     &                                     I7, ' triangles')

 5004 Format('MBA: STONE FLOWER! (1997-2014), version 3.1', /,
     &       5X,'Target Quality=', F4.2, ' (nEStar:', I7, 
     &       ', SkipE:', I6, ', maxITR:', I8,')') 
      End


