FUNCTION MYSENDSET(KSETS,KMYSET,KSET)


!**** *MYSENDSET* RETURNS SET NUMBER TO SEND TO

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        ISENDSET = MYSENDSET(KSETS,KMYSET,KSET)

!        Explicit arguments :  
!        --------------------
!                  input:   KSETS

!        Implicit arguments :  NONE
!        --------------------
!     Method.
!     -------

!     Externals.
!     ----------
!         NONE

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE
INTEGER(KIND=JPIM) :: MYSENDSET
INTEGER(KIND=JPIM),INTENT(IN)  :: KSETS,KMYSET,KSET
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

!*       1.    Check input argument for validity 
!              ---------------------------------

IF (LHOOK) CALL DR_HOOK('MYSENDSET',0,ZHOOK_HANDLE)
IF(KSETS < 1 .OR. KMYSET > KSETS .OR. KSET > KSETS-1) THEN

  CALL ABOR1(' MYSENDSET: INVALID ARGUMENT ')

ELSE

!*       2.    Compute output parameters
!              -------------------------

  MYSENDSET = MOD(KMYSET+KSET-1,KSETS)+1

ENDIF

IF (LHOOK) CALL DR_HOOK('MYSENDSET',1,ZHOOK_HANDLE)
END FUNCTION MYSENDSET
