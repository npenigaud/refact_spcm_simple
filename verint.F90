SUBROUTINE VERINT(KPROMA,KSTART,KPROF,KLEVIN,KLEVOUT,PINTE,PIN,POUT,KTYPE,KCHUNK)

!**** *VERINT*   Vertical integral

!     Purpose.
!     --------
!          This subroutine computes the vertical integral (with respect
!          to eta) of a function given at full model
!          levels using a general scheme
!          The integral is given either from the top (KTYPE=0) down
!          (or from the bottom (KTYPE=1) up) to each full model level

!**   Interface.
!     ----------
!        *CALL* *VERINT(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA    - horizontal dimension.
!          KSTART    - first element of work.
!          KPROF     - depth of work.
!          KLEVIN    - vertical dimension for array PIN.
!          KLEVOUT   - vertical dimension for array POUT.
!          PINTE     - matrix operator used to perform vertical integrals.
!          PIN       - Input field
!          KTYPE     - starting point of the integral (0=top, 1=bottom)
!          KCHUNK    - chunking size (to maintain bit reproducibility)

!        OUTPUT:
!                    eta
!                     _
!                    |
!          POUT :    | PIN deta   at each half model level
!                   _|
!                 KTYPE

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none

!     Reference.

!     Author.
!     -------
!         M. Hortal (ECMWF)

!     Modifications.
!     --------------
!        Original : MAY 2000
!        D.SALMOND : APRIL 2002 FORTRAN Matrix multiply replace by BLAS routine
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        Y.Seity : MAY 2009 bf for B-Level parallelisation
!        J.Hague : OCT 2012: Parallelise call to DGEMM if not in parallel region
!        P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!        F. Vana  05-Mar-2015  Support for single precision
!        F. Vana  14-Jan-2020  Exclusive usage of double precision
!        P. Gillies & F. Vana  22-Jan-2020  Bit reproducible chunking
!        H. Petithomme (September 2020): full rewrite for optimisation and new options
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRD, JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULERR
USE OML_MOD  , ONLY : OML_IN_PARALLEL

#if defined(_OPENACC)
USE CUBLAS
#endif

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA,KLEVIN,KLEVOUT,KSTART,KPROF,KTYPE
REAL(KIND=JPRD)   ,INTENT(IN)    :: PINTE(KLEVOUT,KLEVIN) 
REAL(KIND=JPRB),TARGET,INTENT(IN)    :: PIN(KPROMA,KLEVIN)
REAL(KIND=JPRB),TARGET,INTENT(OUT):: POUT(KPROMA,KLEVOUT) 
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KCHUNK

INTEGER(KIND=JPIM) :: JLEN
LOGICAL,PARAMETER :: LLSINGLE=JPRB/=JPRD
LOGICAL :: LPAR
INTEGER(KIND=JPIM) ::  JLEV, JROF
REAL(KIND=JPRD),CONTIGUOUS,POINTER :: ZIN(:,:)
REAL(KIND=JPRD),CONTIGUOUS,POINTER :: ZOUT(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE,ZHOOK_HANDLE_XGEMM,ZHOOK_HANDLE2

#include "abor1.intfb.h"

IF (KSTART > KPROF) RETURN

IF (LHOOK) CALL DR_HOOK('VERINT',0,ZHOOK_HANDLE)

#ifdef PARKIND1_SINGLE        
  ALLOCATE(ZOUT(KPROMA,KLEVOUT))
  ALLOCATE(ZIN(KPROMA,KLEVIN))
#if defined(_OPENACC)
  !$ACC DATA CREATE(ZIN,ZOUT) PRESENT(PIN)
  !$ACC PARALLEL PRIVATE(JLEV,JROF) DEFAULT(NONE)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO JLEV=1,KLEVIN
    DO JROF=KSTART,KPROF
      ZIN(JROF,JLEV) = PIN(JROF,JLEV)
    ENDDO
  ENDDO
  !$ACC END PARALLEL
#else
  ZIN(KSTART:KPROF,:) = PIN(KSTART:KPROF,:)
#endif
#else
  ZOUT => POUT
  ZIN => PIN
#endif

LPAR = OML_IN_PARALLEL()

!$ACC DATA PRESENT(ZIN,ZOUT,POUT,PINTE)

IF (LPAR) THEN
  IF (LHOOK) CALL DR_HOOK('VERINT_DGEMM_1',0,ZHOOK_HANDLE_XGEMM)

  CALL DGEMM('N','T',KPROF-KSTART+1,KLEVOUT,KLEVIN, &
       & 1.0_JPRD,ZIN,KPROMA,PINTE,KLEVOUT,0.0_JPRD,ZOUT,KPROMA)  

  IF (LHOOK) CALL DR_HOOK('VERINT_DGEMM_1',1,ZHOOK_HANDLE_XGEMM)
ELSE
  IF (LHOOK) CALL DR_HOOK('VERINT_DGEMM_2',0,ZHOOK_HANDLE_XGEMM)

  if (.true.) then
#if defined(_OPENACC)
  !$ACC HOST_DATA USE_DEVICE(ZIN,ZOUT,PINTE)
    CALL CUBLASDGEMM('N','T',KPROMA,KLEVOUT,KLEVIN,&
         & 1.0_JPRD,ZIN,KPROMA,PINTE,KLEVOUT,0.0_JPRD,ZOUT,KPROMA)
  !$ACC END HOST_DATA
  !$ACC WAIT
#else
   ! Chunking across KPROMA
!$OMP PARALLEL DO PRIVATE(JROF,JLEN)
    DO JROF=KSTART,KPROF,KCHUNK
      JLEN=MIN(KCHUNK,KPROF-JROF+1)
      CALL DGEMM('N','T',JLEN,KLEVOUT,KLEVIN, &
           & 1.0_JPRD,ZIN(JROF,1),KPROMA,PINTE,KLEVOUT,0.0_JPRD,ZOUT(JROF,1),KPROMA)
    ENDDO
!$OMP END PARALLEL DO
#endif
  else
    ! Chunking across KLEVOUT
!$OMP PARALLEL DO PRIVATE(JLEV,JLEN)
    DO JLEV=1,KLEVOUT,KCHUNK
      JLEN=MIN(KCHUNK,KLEVOUT-JLEV+1)
      CALL DGEMM('N','T',KPROF-KSTART+1,JLEN,KLEVIN, &
           & 1.0_JPRD,ZIN,KPROMA,PINTE(JLEV,1),KLEVOUT,0.0_JPRD,ZOUT(1,JLEV),KPROMA)
    ENDDO
!$OMP END PARALLEL DO
  end if

  IF (LHOOK) CALL DR_HOOK('VERINT_DGEMM_2',1,ZHOOK_HANDLE_XGEMM)
ENDIF

IF(KTYPE == 1) THEN
  ! warning: dependence on last level in OMP case, last level is done separately
IF (LHOOK) CALL DR_HOOK('VERINT_calcul',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
  !$ACC PARALLEL PRIVATE(JLEV,JROF) if (.not.lpar) DEFAULT(NONE)
  !$ACC LOOP GANG VECTOR COLLAPSE(2) 
  DO JLEV=1,KLEVOUT-1
    DO JROF=KSTART,KPROF
      POUT(JROF,JLEV)=ZOUT(JROF,JLEV)-ZOUT(JROF,KLEVOUT) 
    ENDDO
  ENDDO
  !$ACC END PARALLEL 

  ! last level substraction summarizes to zeroing
  !$ACC PARALLEL private(JROF) DEFAULT(NONE)
  !$ACC LOOP GANG VECTOR
  do JROF=KSTART,KPROF
    POUT(JROF,KLEVOUT)=0._JPRB
  enddo
  !$ACC END PARALLEL
#else

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JROF) if (.not.lpar)
  DO JLEV=1,KLEVOUT-1
    DO JROF=KSTART,KPROF
      POUT(JROF,JLEV)=ZOUT(JROF,JLEV)-ZOUT(JROF,KLEVOUT)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  ! last level substraction summarizes to zeroing
  POUT(KSTART:KPROF,KLEVOUT)=0._JPRB
#endif
IF (LHOOK) CALL DR_HOOK('VERINT_calcul',1,ZHOOK_HANDLE2)

ELSEIF (KTYPE /= 0) THEN
  WRITE(NULERR,*) ' INVALID KTYPE IN VERINT =',KTYPE
  CALL ABOR1(' VERINT: ABOR1 CALLED')
ELSE IF (LLSINGLE) THEN
#if defined(_OPENACC)
!$ACC PARALLEL PRIVATE(JROF,JLEV) DEFAULT(NONE)
!$ACC LOOP GANG VECTOR COLLAPSE(2)
DO JLEV=1,KLEVOUT
  DO JROF=KSTART,KPROF 
    POUT(JROF,JLEV)=ZOUT(JROF,JLEV)
  ENDDO
ENDDO
!$ACC END PARALLEL
#else
!$OMP PARALLEL DO SCHEDULE(STATIC) if (.not.lpar)
  DO JLEV=1,KLEVOUT
    POUT(KSTART:KPROF,JLEV) = ZOUT(KSTART:KPROF,JLEV)
  ENDDO
!$OMP END PARALLEL DO
#endif
ENDIF

#ifdef PARKIND1_SINGLE 
!$ACC END DATA
#endif

!$ACC END DATA
IF (LLSINGLE) THEN
  DEALLOCATE(ZOUT)
  DEALLOCATE(ZIN)
ENDIF

IF (LHOOK) CALL DR_HOOK('VERINT',1,ZHOOK_HANDLE)
END SUBROUTINE VERINT
