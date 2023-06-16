SUBROUTINE TRMTOS(YDGEOMETRY,LDNHDYN,LDNHX,LDTRANSPOSE,PSPVOR,PSPDIV,PSPT,PSPSPD,PSPSVD,PSPSNHX, &
 & PSPGFL,PSPSP,PSPAUX,PSPSEL3D,PSPSEL2D,&
 & PSPVORG,PSPDIVG,PSPTG,PSPSPDG,PSPSVDG,PSPSNHXG,&
 & PSPGFLG,PSPSPG,PSPAUXG,PSPSEL3DG,PSPSEL2DG,&
 & LDSELECT3D,LDSELECT2D,LDFULLM)  

!**** *trmtos * - Transposition from horizontal to vertical spectral
!                 coefficients

!     Purpose.
!     --------
!              Transpose data from SPXX arrays to SPxxG arrays.
!              The semi-implicit calculations require vertical
!              columns of PSPT, PSPDIV, and PSPSP.
!              This is done during the spectral space calculations.
!              Also used in Jb calculations.
!              This routine is the inverse of TRSTOM.


!**   Interface.
!     ----------
!        *call* *trmtos

!        Explicit arguments :
!        --------------------
!          PSPVOR     -  distributed vorticity
!          PSPDIV     -  distributed divergence
!          PSPT       -  distributed temperature
!          PSPSPD     -  distributed NH pressure departure variable
!          PSPSVD     -  distributed NH vertical divergence variable
!          PSPSNHX    -  distributed NH "X" part divergence variable
!          PSPGFL     -  distributed spectral GFL
!          PSPSP      -  distributed surface pressure
!          PSPAUX     -  distributed auxiliary field
!          PSPSEL3D   -  3d fields for selection
!          PSPSEL2D   -  2d fields for selection
!          PSPVORG    -  vorticity columns
!          PSPDIVG    -  divergence columns
!          PSPTG      -  temperature columns
!          PSPSPDG    -  NH pressure departure variable columns
!          PSPSVDG    -  NH vertical divergence variable columns
!          PSPSNHXG   -  NH "X" part divergence variable columns
!          PSPGFLG    -  GFL columns
!          PSPSPG     -  surface pressure
!          PSPAUX     -  auxiliary field columns
!          PSPSEL3DG  -  selected 3d fields
!          PSPSEL2DG  -  selected 2d fields
!          LDSELECT3D - .T.=> select this 3d field
!          LDSELECT2D - .T.=> select this 2d field
!          LDFULLM    - .T. if full m-columns are used

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      MPP Group *ECMWF*
!      Original       : 95-10-01 

!     Modifications.
!     --------------
!      M.Hamrud       : 03-10-01 CY28 Cleaning
!      M.Hamrud       : 03-12-01 CY28R1 Cleaning
!      K. Yessad      : 05-02-08 NH variables + cleaning
!      L. Isaksen     : 04-09-01 Optional arguments, more flexible
!      M. Fisher      : 05-08-22 Even more flexible (selected 2d,3d fields)
!      G. Mozdzynski  : 08-01-01 Cleanup
!      F. Vana + NEC  : 08-09-09 optimization
!      T. Wilhelmsson : 09-09-22 Add auxiliary field
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib 02-Jul-2015 Optimization (partial)
!      R. El Khatib 10-Dec-2020 Optimization by overlaping packs/unpacks with comms
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE EXCHANGE_MS_MOD, ONLY : FIELDLIST, ADD3DF, ADD3DFL, ADD2DF, ADD2DFL, NEXCHANGE_MTOS, EXCHANGE_MS,TERMINATE_LIST

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,          INTENT(IN)    :: YDGEOMETRY
LOGICAL           ,          INTENT(IN)    :: LDNHDYN
LOGICAL           ,          INTENT(IN)    :: LDNHX
LOGICAL           ,          INTENT(IN)    :: LDTRANSPOSE
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPVOR(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPDIV(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPT(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSPD(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSVD(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSNHX(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPGFL(:,:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSP(:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPAUX(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSEL3D(:,:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSEL2D(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPVORG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPDIVG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPTG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPSPDG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPSVDG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPSNHXG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(INOUT) :: PSPGFLG(:,:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPSPG(:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPAUXG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(INOUT) :: PSPSEL3DG(:,:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(INOUT) :: PSPSEL2DG(:,:) 
LOGICAL           ,OPTIONAL, INTENT(IN)    :: LDSELECT3D(:)
LOGICAL           ,OPTIONAL, INTENT(IN)    :: LDSELECT2D(:)
LOGICAL           ,OPTIONAL, INTENT(IN)    :: LDFULLM 

TYPE (FIELDLIST) :: YLLIST

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TRMTOS',0,ZHOOK_HANDLE)

!$ACC ENTER DATA COPYIN(YLLIST)
CALL ADD3DF (YLLIST, PSPVOR, PSPVORG, "PSPVOR")
CALL ADD3DF (YLLIST, PSPDIV, PSPDIVG, "PSPDIV")
CALL ADD3DF (YLLIST, PSPT  , PSPTG,   "PSPT"  )
CALL ADD3DF (YLLIST, PSPAUX, PSPAUXG, "PSPAUX")
IF (LdNHDYN) CALL ADD3DF (YLLIST, PSPSPD,  PSPSPDG,  "PSPSPD")
IF (LdNHDYN) CALL ADD3DF (YLLIST, PSPSVD,  PSPSVDG,  "PSPSVD")
IF (LdNHX)   CALL ADD3DF (YLLIST, PSPSNHX, PSPSNHXG, "PSPSNHX")
CALL ADD3DFL (YLLIST, PSPGFL,   PSPGFLG,   "PSPGFL")
CALL ADD3DFL (YLLIST, PSPSEL3D, PSPSEL3DG, "PSPSEL3D", LDSELECT3D)
CALL ADD2DF (YLLIST, PSPSP, PSPSPG, "PSPSP")
CALL ADD2DFL (YLLIST, PSPSEL2D, PSPSEL2DG, "PSPSEL2D", LDSELECT2D)

CALL EXCHANGE_MS (YDGEOMETRY,LDTRANSPOSE,YLLIST, KDIR=NEXCHANGE_MTOS, LDFULLM=LDFULLM)

#if defined(_OPENACC)
CALL TERMINATE_LIST(YLLIST)
#endif

IF (LHOOK) CALL DR_HOOK('TRMTOS',1,ZHOOK_HANDLE)

END SUBROUTINE TRMTOS
