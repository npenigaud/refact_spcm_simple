#if defined(_OPENACC)
SUBROUTINE MXTURS(KLX,KVX,KVXS,KIX,KIXS,tnsmax,PA,PB,PC,PY,PX,tbloc,pas,pbs,pcs,pys,pxs)
!$acc routine vector
#else
SUBROUTINE MXTURS(KLX,KVX,KVXS,KIX,tnsmax,PA,PB,PC,PY,PX)
#endif

!**** *MXTURS*   - Resolution of a set of pentadiagonal symmetric systems.

!     Purpose.    Resolution of a set of pentadiagonal symmetric systems,
!     --------    once known the LU decomposition of these symmetric matrixes.

!    Example, for KLX=4:

!    PY is known, PX is unknown. We want to find PX, solution of the
!    following system, for l=1 to KVX, i=1 to KIX.
!                            S(l)*PX(l,i)=PY(l,i)
!    where S(l) are symmetric pentadiagonal matrixes, the LU decomposition
!    of which yields:

!           I PA(l,1)    0       0       0    I
!    L(l) = I PB(l,1) PA(l,2)    0       0    I
!           I PC(l,1) PB(l,2) PA(l,3)    0    I
!           I    0    PC(l,2) PB(l,3) PA(l,4) I

!           I    1    ZB(l,1) ZC(l,1)    0    I
!    U(l) = I    0       1    ZB(l,2) ZC(l,2) I
!           I    0       0       1    ZB(l,3) I
!           I    0       0       0       1    I

!    where ZB(l,j)=PB(l,j)/PA(l,j); ZC(l,j)=PC(l,j)/PA(l,j) .

!    S(l)=L(l)*U(l)

!    We call routine MXTURE to inverse L, then U.

!**   Interface.
!     ----------
!        *CALL* *MXTURS(KLX,KVX,KVXS,KIX,PA,PB,PC,PY,PX)

!        Explicit arguments :
!        --------------------
!         KLX:        - Dimension of the system.                    (input)
!         KVX:        - Number of variables (second                 (input)
!                       dimension of PX and PY).
!         KVXS:       - Surdimension corresponding to KVX.          (input)
!         KIX:        - Number of variables multiplied by the same
!                       matrix.                                     (input)
!         PA,PB,PC:   - non-zero diagonals of the triangular        (input)
!                       matrixes L (see figures above).
!                       Caution! All PA coefficients must be non 0.
!         PY:         - known vector.                               (input)
!         PX:         - unknown vector.                             (output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Calls MXTURE.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD: MARCH 1994.

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!     ------------------------------------------------------------------

IMPLICIT NONE

#if defined(_OPENACC)
INTEGER(KIND=JPIM),INTENT(IN),value :: KLX
INTEGER(KIND=JPIM),INTENT(IN),value :: KVX
INTEGER(KIND=JPIM),INTENT(IN),value :: KVXS
INTEGER(KIND=JPIM),INTENT(IN),value :: KIX
INTEGER(KIND=JPIM),INTENT(IN),value :: KIXS
INTEGER(KIND=JPIM),INTENT(IN),value :: tnsmax
REAL(KIND=JPRB) ,INTENT(IN) :: PA(klx,kvx)
REAL(KIND=JPRB) ,INTENT(IN) :: PB(klx,kvx)
REAL(KIND=JPRB) ,INTENT(IN) :: PC(klx,kvx)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PY(tnsmax+1,kixs,kvx)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PX(tnsmax+1,kixs,kvx)
integer(kind=jpim),intent(in),value :: tbloc
REAL(KIND=JPRB) ,INTENT(INout) :: PAs(tbloc+3,kvxs)
REAL(KIND=JPRB) ,INTENT(INout) :: PBs(tbloc+3,kvxs)
REAL(KIND=JPRB) ,INTENT(INout) :: PCs(tbloc+3,kvxs)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PYs(tbloc+3,kvxs,kixs)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PXs(tbloc+3,kvxs,kixs)
#else
INTEGER(KIND=JPIM),INTENT(IN) :: KLX
INTEGER(KIND=JPIM),INTENT(IN) :: KVX
INTEGER(KIND=JPIM),INTENT(IN) :: KVXS
INTEGER(KIND=JPIM),INTENT(IN) :: KIX
integer(kind=JPIM),intent(in) :: tnsmax
REAL(KIND=JPRB) ,INTENT(IN) :: PA(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PB(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PC(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PY(KVXS,tnsmax+1,KIX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PX(KVXS,tnsmax+1,KIX)
#endif


!     ------------------------------------------------------------------


#include "mxture.h"

!     ------------------------------------------------------------------
!!IF (LHOOK) CALL DR_HOOK('MXTURS',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    INVERSION OF THE TWO TRIANGULAR TRIDIAGONAL MATRIXES.
!              -----------------------------------------------------
#if defined(_OPENACC)
CALL MXTURE(KLX,KVX,KVXS,KIX,KIXS,tnsmax,-2,.TRUE. ,PA,PB,PC,PY,PX,tbloc,pas,pbs,pcs,pys,pxs)
CALL MXTURE(KLX,KVX,KVXS,KIX,KIXS,tnsmax, 1,.FALSE.,PA,PB,PC,PY,PX,tbloc,pas,pbs,pcs,pys,pxs)
#else
CALL MXTURE(KLX,KVX,KVXS,KIX,tnsmax,-2,.TRUE. ,PA,PB,PC,PY,PX)
CALL MXTURE(KLX,KVX,KVXS,KIX, tnsmax,1,.FALSE.,PA,PB,PC,PY,PX)
#endif
!     ------------------------------------------------------------------

!!IF (LHOOK) CALL DR_HOOK('MXTURS',1,ZHOOK_HANDLE)
END SUBROUTINE MXTURS

