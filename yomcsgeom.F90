MODULE YOMCSGEOM

USE PARKIND1, ONLY : JPIM, JPRB, JPRD

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    * Computational sphere horizontal geometry

!     RCOLON(NGPTOT) cosine of longitude on transformed sphere
!     RSILON(NGPTOT)   sine        "             "         "
!     RINDX (NGPTOT) Longitude index
!     RINDY (NGPTOT) Latitude index
!     RATATH(NGPTOT) RA*TAN(THETA) on real sphere
!     RATATX(NGPTOT) Curvature term for LAM (for u eq.)

TYPE TCSGEOM
  REAL(KIND=JPRD), POINTER :: RCOLON(:) => NULL()
  REAL(KIND=JPRD), POINTER :: RSILON(:) => NULL()
  REAL(KIND=JPRB), POINTER :: RINDX (:) => NULL()
  REAL(KIND=JPRB), POINTER :: RINDY (:) => NULL()
  REAL(KIND=JPRB), POINTER :: RATATH(:) => NULL()
  REAL(KIND=JPRB), POINTER :: RATATX(:) => NULL()
END TYPE TCSGEOM

TYPE TCSGEOM_BLOCKED
  REAL(KIND=JPRD), POINTER :: RCOLON(:,:) => NULL()
  REAL(KIND=JPRD), POINTER :: RSILON(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: RINDX (:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: RINDY (:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: RATATH(:,:) => NULL()
  REAL(KIND=JPRB), POINTER :: RATATX(:,:) => NULL()
END TYPE TCSGEOM_BLOCKED

END
