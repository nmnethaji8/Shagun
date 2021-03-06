!!----------------------------COMMONMOD----------------------------!!
MODULE COMMONMOD
IMPLICIT NONE  
  
  ! PROTECTED ensures LNODE can only be modified inside COMMONMOD
  INTEGER(KIND=4),PROTECTED::LNODE
  INTEGER(KIND=4),PARAMETER::NLMAXN=1100,NFSB = 850

  INTEGER(KIND=4),ALLOCATABLE::NODEID(:),NWALLID(:,:)
  INTEGER(KIND=4)::ICELLX,ICELLY,ICELLZ
  INTEGER(KIND=4)::IXMAX,IYMAX,IZMAX
  INTEGER(KIND=4)::KW, MBAS, I_WM
  INTEGER(KIND=4),PARAMETER::MABSM=10

  REAL(KIND=8),ALLOCATABLE::COORX(:,:),COORY(:,:),COORZ(:,:)
  REAL(KIND=8),ALLOCATABLE::SNX(:),SNY(:),SNZ(:)
  REAL(KIND=8),ALLOCATABLE::SMX(:),SMY(:),SMZ(:)
  REAL(KIND=8),ALLOCATABLE::SSX(:),SSY(:),SSZ(:)
  !REAL(KIND=8),ALLOCATABLE::DN1(:)
  REAL(KIND=8),ALLOCATABLE::TANK_A(:,:)
  REAL(KIND=8)::h, SLOPE
  REAL(KIND=8)::DDL, DT, TOTAL_TIME, H0
  !REAL(KIND=8)::DN00, DN01

  REAL(KIND=8),PARAMETER::GRA=-9.81D0, PI=ATAN(1D0)*4D0
  REAL(KIND=8)::R1=2.1D0!,R2=3.1D0,R3=2.1D0
  REAL(KIND=8)::RCELL=3.2D0

CONTAINS

  SUBROUTINE INITCOMMONMOD
  IMPLICIT NONE

    ALLOCATE(NODEID(-7:LNODE), NWALLID(LNODE,4))

    ALLOCATE(COORX(LNODE,3), COORY(LNODE,3), COORZ(LNODE,3))
    ALLOCATE(SNX(LNODE), SNY(LNODE), SNZ(LNODE))
    ALLOCATE(SMX(LNODE), SMY(LNODE), SMZ(LNODE))
    ALLOCATE(SSX(LNODE), SSY(LNODE), SSZ(LNODE))
    !ALLOCATE(DN1(LNODE))
    ALLOCATE(TANK_A(LNODE,3))

  END SUBROUTINE INITCOMMONMOD

  

  SUBROUTINE SETLNODE(LNODEIN)
  IMPLICIT NONE

    INTEGER(KIND=4),INTENT(IN)::LNODEIN

    LNODE = LNODEIN

  END SUBROUTINE SETLNODE


END MODULE COMMONMOD
!!--------------------------END COMMONMOD--------------------------!!