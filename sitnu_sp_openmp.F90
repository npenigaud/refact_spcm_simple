SUBROUTINE SITNU_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP)

!**** *SITNU_SP_OPENMP*   - Continuity equation for semi-implicit.

!     Purpose.
!     --------
!           Evaluate operators Tau and Nu in semi-implicit.

!**   Interface.
!     ----------
!        *CALL* *SITNU_SP_OPENMP(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME LEVEL

!           TYPICAL VALUES ARE  NDLSUR,1  FOR GRID POINT ARRAY
!                               1,NFLSUR  FOR SPECTRAL ARRAY

!        PD    : DIVERGENCE
!        PT    : TEMPERATURE
!        PSP   : SURFACE PRESSURE

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      Modified : 09-Oct-2007 by K. YESSAD: possibility to have a specific
!                 value of LVERTFE in the SI NH linear model.
!      F. Vana + NEC 28-Apr-2009: OpenMP + optimization
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      G. Mozdzynski Oct 2012: OpenMP optimization
!      K. Yessad (Dec 2016): Prune obsolete options.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      R.Brozkova + NEC: Mar 2021: Optimization for vector (NEC)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : TCST
USE YOMDYN       , ONLY : TDYN


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC

REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KSPEC,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT(KSPEC,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSP(KSPEC) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSDIV(KSPEC,0:KLEV+1)
REAL(KIND=JPRB) :: ZOUT(KSPEC,0:KLEV)

REAL(KIND=JPRB) :: ZSDIVX(0:KLEV, KSPEC)
INTEGER(KIND=JPIM) :: JLEV, JSPEC
REAL(KIND=JPRB) :: ZREC
REAL(KIND=JPRB) :: ZDETAH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE,ZHOOK_HANDLE2

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SITNU_SP_OPENMP',0,ZHOOK_HANDLE)

ASSOCIATE( YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDCVER=>YDGEOMETRY%YRCVER) 
ASSOCIATE(SIALPH=>YDDYN%SIALPH, SIDELP=>YDDYN%SIDELP, SILNPR=>YDDYN%SILNPR, SIRDEL=>YDDYN%SIRDEL, &
 & SIRPRN=>YDDYN%SIRPRN, SITLAF=>YDDYN%SITLAF, SITR=>YDDYN%SITR, YDDIMV=>YDGEOMETRY%YRDIMV)

!     ------------------------------------------------------------------

!*       1.    SUM DIVERGENCE AND COMPUTES TEMPERATURE.
!              ----------------------------------------

IF(YDCVER%LVERTFE) THEN
!$ACC DATA PRESENT(PD,PSP,PT,YDDYN%SIDELP,YDDYN%SITLAF,YDGEOMETRY%YRVETA,YDCST,YDDYN%SIRPRN)
!$ACC DATA CREATE(ZSDIV,ZOUT) PRESENT(YDDYN,YDGEOMETRY,YDGEOMETRY%YRVETA%VFE_RDETAH) 

IF (LHOOK) CALL DR_HOOK('SITNU_transpose1',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
!$ACC PARALLEL PRIVATE(JLEV,JSPEC,ZDETAH) DEFAULT(NONE)
!$ACC LOOP GANG 
#else
!$OMP PARALLEL PRIVATE(JLEV,JSPEC,ZDETAH)
!$OMP DO SCHEDULE(STATIC) 
#endif
  DO JLEV=1,KLEV
    ZDETAH=YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)*YDDYN%SIDELP(JLEV)
    !$ACC LOOP VECTOR
    DO JSPEC=1,KSPEC
      ZSDIV(JSPEC,JLEV)=PD(JSPEC,JLEV)*ZDETAH
    ENDDO
  ENDDO
#if defined(_OPENACC)
!$ACC END PARALLEL
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

IF (LHOOK) CALL DR_HOOK('SITNU_transpose1',1,ZHOOK_HANDLE2)

  IF (KSPEC>=1) THEN
!DEC$ IVDEP
IF (LHOOK) CALL DR_HOOK('SITNU_cond_lim',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
    !$ACC PARALLEL PRIVATE(JSPEC) DEFAULT(NONE)
    !$ACC LOOP GANG
    do JSPEC=1,KSPEC
      ZSDIV(JSPEC,0)=0.0_JPRB
      ZSDIV(JSPEC,KLEV+1)=0.0_JPRB
    enddo
    !$ACC END PARALLEL
#else
    ZSDIV(1:KSPEC,0)=0.0_JPRB
    ZSDIV(1:KSPEC,KLEV+1)=0.0_JPRB
#endif
IF (LHOOK) CALL DR_HOOK('SITNU_cond_lim',1,ZHOOK_HANDLE2)

    CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',KSPEC,1,KSPEC,KLEV,ZSDIV,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)
  ENDIF

IF (LHOOK) CALL DR_HOOK('SITNU_transpose2',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
!$ACC PARALLEL PRIVATE(JLEV,JSPEC,ZREC) DEFAULT(NONE)
!$ACC LOOP GANG
#else
!$OMP PARALLEL PRIVATE(JLEV,JSPEC,ZREC)
!$OMP DO SCHEDULE(STATIC) 
#endif
  DO JLEV=1,KLEV
    ZREC=(1.0_JPRB/YDDYN%SITLAF(JLEV))*YDCST%RKAPPA*YDDYN%SITR
    !$ACC LOOP VECTOR
    DO JSPEC=1,KSPEC
      PT(JSPEC,JLEV)=ZOUT(JSPEC,JLEV-1)*ZREC
    ENDDO
  ENDDO
#if defined(_OPENACC)
!$ACC END PARALLEL
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

IF (LHOOK) CALL DR_HOOK('SITNU_transpose2',1,ZHOOK_HANDLE2)

IF (LHOOK) CALL DR_HOOK('SITNU_calcul1',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
  !$ACC PARALLEL LOOP PRIVATE(JSPEC) DEFAULT(NONE)
  DO JSPEC=1,KSPEC
    PSP(JSPEC)=ZOUT(JSPEC,KLEV)*YDDYN%SIRPRN
  ENDDO
  !$ACC END PARALLEL
#else
  DO JSPEC=1,KSPEC
    PSP(JSPEC)=ZOUT(JSPEC,KLEV)*SIRPRN
  ENDDO
#endif
IF (LHOOK) CALL DR_HOOK('SITNU_calcul1',1,ZHOOK_HANDLE2)

!$ACC END DATA
!$ACC END DATA

ELSE

  ZSDIVX(0, :)=0.0_JPRB

!$OMP PARALLEL PRIVATE(JLEV,JSPEC)
!$OMP DO SCHEDULE(STATIC)
  DO JSPEC=1,KSPEC
    DO JLEV=1,KLEV
      ZSDIVX(JLEV, JSPEC)=ZSDIVX(JLEV-1, JSPEC)+PD(JLEV,JSPEC)*SIDELP(JLEV)
      PT(JLEV,JSPEC)=YDCST%RKAPPA*SITR*(SIRDEL(JLEV)*SILNPR(JLEV)*ZSDIVX(JLEV-1, JSPEC)&
       & +SIALPH(JLEV)*PD(JLEV,JSPEC))
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  DO JSPEC=1,KSPEC
    PSP(JSPEC)=ZSDIVX(KLEV, JSPEC)*SIRPRN
  ENDDO

ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SITNU_SP_OPENMP',1,ZHOOK_HANDLE)

END SUBROUTINE SITNU_SP_OPENMP

