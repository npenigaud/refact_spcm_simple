SUBROUTINE SITNU_GP (YDGEOMETRY, YDCST, YDDYN, KLON, KLEV, PD, PT, PSP)

!**** *SITNU_GP*   - Continuity equation for semi-implicit.

!     Purpose.
!     --------
!           Evaluate operators Tau and Nu in semi-implicit.

!**   Interface.
!     ----------
!        *CALL* *SITNU_GP(...)

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
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSP(KLON) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSDIV(KLON,0:KLEV+1)
REAL(KIND=JPRB) :: ZOUT(KLON,0:KLEV)


REAL(KIND=JPRB) :: ZSDIVX(KLON, 0:KLEV)
INTEGER(KIND=JPIM) :: JLEV, JLON
REAL(KIND=JPRB) :: ZREC
REAL(KIND=JPRB) :: ZDETAH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SITNU_GP',0,ZHOOK_HANDLE)

ASSOCIATE( YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDCVER=>YDGEOMETRY%YRCVER) 
ASSOCIATE(SIALPH=>YDDYN%SIALPH, SIDELP=>YDDYN%SIDELP, SILNPR=>YDDYN%SILNPR, SIRDEL=>YDDYN%SIRDEL, &
 & SIRPRN=>YDDYN%SIRPRN, SITLAF=>YDDYN%SITLAF, SITR=>YDDYN%SITR, YDDIMV=>YDGEOMETRY%YRDIMV)

!     ------------------------------------------------------------------

!*       1.    SUM DIVERGENCE AND COMPUTES TEMPERATURE.
!              ----------------------------------------

IF(YDCVER%LVERTFE) THEN

  DO JLEV=1,KLEV
    ZDETAH=YDVETA%VFE_RDETAH(JLEV)
    DO JLON=1,KLON
      ZSDIV(JLON,JLEV)=PD(JLON,JLEV)*SIDELP(JLEV)*ZDETAH
    ENDDO
  ENDDO

  IF (KLON>=1) THEN
!DEC$ IVDEP
    ZSDIV(1:KLON,0)=0.0_JPRB
    ZSDIV(1:KLON,KLEV+1)=0.0_JPRB
    CALL VERDISINT(YDVFE,YDCVER,'ITOP','11',KLON,1,KLON,KLEV,ZSDIV,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)
  ENDIF

  DO JLEV=1,KLEV
    DO JLON=1,KLON
      ZREC=1.0_JPRB/SITLAF(JLEV)
      PT(JLON,JLEV)=YDCST%RKAPPA*SITR*ZOUT(JLON,JLEV-1)*ZREC
    ENDDO
  ENDDO
  DO JLON=1,KLON
    PSP(JLON)=ZOUT(JLON,KLEV)*SIRPRN
  ENDDO

ELSE

  ZSDIVX(:, 0)=0.0_JPRB

  DO JLEV=1,KLEV
    DO JLON=1,KLON
      ZSDIVX(JLON, JLEV)=ZSDIVX(JLON, JLEV-1)+PD(JLON,JLEV)*SIDELP(JLEV)
      PT(JLON,JLEV)=YDCST%RKAPPA*SITR*(SIRDEL(JLEV)*SILNPR(JLEV)*ZSDIVX(JLON, JLEV-1)&
       & +SIALPH(JLEV)*PD(JLON,JLEV))
    ENDDO
  ENDDO

  DO JLON=1,KLON
    PSP(JLON)=ZSDIVX(JLON, KLEV)*SIRPRN
  ENDDO

ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SITNU_GP',1,ZHOOK_HANDLE)

END SUBROUTINE SITNU_GP

