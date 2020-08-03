C ======================================================================
C Examples of scalar coeffients:
C Dnull       = coefficient is one (makes assembling faster)
C D1x1_one    = coefficient is one
C D1x1_scalar = coefficient equals dDATA
C ======================================================================
      Integer Function ANI_Dnull(x, y, label, dDATA, iDATA, iSYS, Coef)
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      ANI_Dnull = TENSOR_NULL

      Return
      End


C ======================================================================
      Integer Function ANI_D1x1_one(x, y, label, dDATA,iDATA,iSYS, Coef)
C ======================================================================
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = 1D0

      ANI_D1x1_one = TENSOR_SCALAR

      Return
      End



C ======================================================================
      Integer Function ANI_D1x1_scalar(x, y, label, 
     &                                 dDATA, iDATA, iSYS, Coef)
C ======================================================================
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = dDATA(1)

      ANI_D1x1_scalar = TENSOR_SCALAR

      Return
      End



C ======================================================================
C Examples of tensors
C ANI_D2x2_one    = identity tensor
C ANI_D2x2_scalar = scalar tensor equal to dDATA 
C ANI_D2x2_diag   = diagonal tensor equal to diag{dDATA(1), dDATA(2)} 
C ======================================================================
      Integer Function ANI_D2x2_one(x, y, label, dDATA,iDATA,iSYS, Coef)
C ======================================================================
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 2

      Do i = 1, 2
         Do j = 1, 2
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = 1D0
      End do

      ANI_D2x2_one = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
      Integer Function ANI_D2x2_scalar(x, y, label, 
     &                                 dDATA, iDATA, iSYS, Coef)
C ======================================================================
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 2

      Do i = 1, 2
         Do j = 1, 2
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = dDATA(1)
      End do

      ANI_D2x2_scalar = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
      Integer Function ANI_D2x2_diag(x, y, label, 
     &                               dDATA, iDATA, iSYS, Coef)
C ======================================================================
      include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 2
      iSYS(2) = 2

      Do i = 1, 2
         Do j = 1, 2
            Coef(i, j) = 0D0
         End do
         Coef(i, i) = dDATA(i)
      End do

      ANI_D2x2_diag = TENSOR_SYMMETRIC

      Return
      End



C ======================================================================
C Examples of boundary conditions
C ANI_BC_zero   = zero Dirichlet boundary condition
C ANI_BC_scalar = Dirichlet boundary condition equals to dDATA
C ======================================================================
      Integer Function ANI_BC_zero(x, y, label, dDATA,iDATA, iSYS, Coef)
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = 0D0

      ANI_BC_zero = BC_DIRICHLET

      Return
      End



C ======================================================================
      Integer Function ANI_BC_scalar(x, y, label, 
     &                               dDATA, iDATA, iSYS, Coef)
C ======================================================================
      Include 'fem2Dtri.fd'

      Real*8  dDATA(*), x, y, Coef(MaxTensorSize, *)
      Integer iDATA(*), label, iSYS(*)

      iSYS(1) = 1
      iSYS(2) = 1

      Coef(1, 1) = dDATA(1)

      ANI_BC_scalar = BC_DIRICHLET

      Return
      End




