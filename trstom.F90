SUBROUTINE TRSTOM(YDGEOMETRY,LDNHDYN,LDNHX,LDTRANSPOSE,PSPVORG,PSPDIVG,PSPTG,PSPSPDG,PSPSVDG,PSPSNHXG, &
 & PSPGFLG,PSPSPG,PSPAUXG,PSPSEL3DG,PSPSEL2DG,&
 & PSPVOR,PSPDIV,PSPT,PSPSPD,PSPSVD,PSPSNHX,&
 & PSPGFL,PSPSP,PSPAUX,PSPSEL3D,PSPSEL2D,&
 & LDSELECT3D,LDSELECT2D,LDFULLM,LDNEEDPS)  

!**** *trstom * - Transposition from vertical to horizontal spectral
!                 coefficients

!     Purpose.
!     --------
!              Transpose data from SPXXG arrays to SPxx arrays.
!              The semi-implicit calculations require vertical
!              columns of PSPT, PSPDIV, and restructured PSPSP.
!              This routine is called during the spectral computations.
!              This routine is the inverse of TRMTOS.

!**   Interface.
!     ----------
!        *call* *trstom(...)

!        Explicit arguments :
!        --------------------
!         PSPVORG     - vorticity columns
!         PSPDIVG     - divergence columns
!         PSPTG       - temperature columns
!         PSPSPDG     - NH pressure departure variable columns
!         PSPSVDG     - NH vertical divergence variable columns
!         PSPSNHXG    - NH "X" part divergence variable columns
!         PSPGFLG     - GFL columns
!         PSPSPG      - surface pressure
!         PSPAUX      - auxiliary field columns
!         PSPSEL3DG   - 3d fields
!         PSPSEL2DG   - 2d fields
!         PSPVOR      - distributed vorticity
!         PSPDIV      - distributed divergence
!         PSPT        - distributed temperature
!         PSPSPD      - distributed NH pressure departure variable
!         PSPSVD      - distributed NH vertical divergence variable
!         PSPSNHX     - distributed NH "X" part divergence variable
!         PSPGFL      - distributed humidity
!         PSPSP       - distributed surface pressure
!         PSPAUX      - distributed auxiliary field
!         PSPSEL3D    - 3d fields (only those indicated
!                       by LSELECT3D are updated)
!         PSPSEL2D    - 2d fields (only those indicated
!                       by LSELECT2D are updated)
!         LDSELECT3D  - .T.=> select this 3d field
!         LDSELECT2D  - .T.=> select this 2d field
!         LDFULLM     - .T. if full m-columns are used
!         LDNEEDPS    - .T. if all PE's recv Ps

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
!      D.Salmond      : 01-11-23 LIMP_NOOLAP Option for non-overlapping
!                                message passing and buffer packing
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
USE EXCHANGE_MS_MOD, ONLY : FIELDLIST, ADD3DF, ADD3DFL, ADD2DF, ADD2DFL, NEXCHANGE_STOM, EXCHANGE_MS,TERMINATE_LIST

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),INTENT(IN)    :: YDGEOMETRY
LOGICAL       ,INTENT(IN)    :: LDNHDYN 
LOGICAL       ,INTENT(IN)    :: LDNHX
LOGICAL       ,INTENT(IN)    :: LDTRANSPOSE
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPVORG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPDIVG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPTG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSPDG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSVDG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSNHXG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPGFLG(:,:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSPG(:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPAUXG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSEL3DG(:,:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(IN)    :: PSPSEL2DG(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPVOR(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPDIV(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPT(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPSPD(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPSVD(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPSNHX(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(INOUT) :: PSPGFL(:,:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPSP(:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(OUT)   :: PSPAUX(:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(INOUT) :: PSPSEL3D(:,:,:) 
REAL(KIND=JPRB)   ,OPTIONAL, INTENT(INOUT) :: PSPSEL2D(:,:) 
LOGICAL           ,OPTIONAL, INTENT(IN)    :: LDSELECT3D(:)
LOGICAL           ,OPTIONAL, INTENT(IN)    :: LDSELECT2D(:)
LOGICAL           ,OPTIONAL, INTENT(IN)    :: LDFULLM 
LOGICAL           ,OPTIONAL, INTENT(IN)    :: LDNEEDPS 

TYPE (FIELDLIST) :: YLLIST

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TRSTOM',0,ZHOOK_HANDLE)

!$ACC ENTER DATA COPYIN(YLLIST)
CALL ADD3DF (YLLIST, PSPVOR, PSPVORG, "PSPVOR")
CALL ADD3DF (YLLIST, PSPDIV, PSPDIVG, "PSPDIV")
CALL ADD3DF (YLLIST, PSPT  , PSPTG,   "PSPT"  )
CALL ADD3DF (YLLIST, PSPAUX, PSPAUXG, "PSPAUX")
IF (LDNHDYN) CALL ADD3DF (YLLIST, PSPSPD,  PSPSPDG,  "PSPSPD")
IF (LDNHDYN) CALL ADD3DF (YLLIST, PSPSVD,  PSPSVDG,  "PSPSVD")
IF (LDNHX)   CALL ADD3DF (YLLIST, PSPSNHX, PSPSNHXG, "PSPSNHX")
CALL ADD3DFL (YLLIST, PSPGFL,   PSPGFLG,   "PSPGFL")
CALL ADD3DFL (YLLIST, PSPSEL3D, PSPSEL3DG, "PSPSEL3D", LDSELECT3D)
CALL ADD2DF (YLLIST, PSPSP, PSPSPG, "PSPSP", LDBCAST=LDNEEDPS)
CALL ADD2DFL (YLLIST, PSPSEL2D, PSPSEL2DG, "PSPSEL2D", LDSELECT2D)

CALL EXCHANGE_MS (YDGEOMETRY,LDTRANSPOSE,YLLIST, KDIR=NEXCHANGE_STOM, LDFULLM=LDFULLM)

#if defined(_OPENACC)
CALL TERMINATE_LIST(YLLIST)
#endif

IF (LHOOK) CALL DR_HOOK('TRSTOM',1,ZHOOK_HANDLE)

END SUBROUTINE TRSTOM
