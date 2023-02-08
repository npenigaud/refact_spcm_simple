
SUBROUTINE SIGAM_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP)

!**** *SIGAM_SP_OPENMP* - Solve hydrostatic operator in semi-implicit

!     Purpose.
!     --------
!           Operator gamma to compute p.

!**   Interface.
!     ----------
!        *CALL* *SIGAM_SP_OPENMP(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME LEVEL

!           TYPICAL VALUES ARE  NDLSUR,1  FOR GRID POINT ARRAY
!                               1,NFLSUR  FOR SPECTRAL ARRAY

!        PD    : DIVERGENCE       (output)
!        PT    : TEMPERATURE      (input)
!        PSP   : SURFACE PRESSURE (input)
!        KNLON : NUMBER OF VERTICAL COLUMNS TREATED
!        KFLEVG: NUMBER OF ELEMENTS IN A VERTICAL COLUMN

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
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      Modified : 09-Oct-2007 by K. YESSAD: possibility to have a specific
!                 value of LVERTFE in the SI linear model.
!      F. Vana + NEC 28-Apr-2009: OpenMP
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      G. Mozdzynski Oct 2012: OpenMP optimization
!      K. Yessad (Dec 2016): Prune obsolete options.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      R.Brozkova + NEC 03-Mar-2021: Optimization for vector (NEC)
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
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KSPEC)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSPHI(KSPEC,0:KLEV+1)
REAL(KIND=JPRB) :: ZOUT(KSPEC,0:KLEV)


REAL(KIND=JPRB) :: ZSPHIX(0:KLEV, KSPEC)
INTEGER(KIND=JPIM) :: JLEV, JSPEC
CHARACTER(LEN=4):: CLOPER
REAL(KIND=JPRB) :: ZDETAH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIGAM_SP_OPENMP',0,ZHOOK_HANDLE)

ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE, YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(SIALPH=>YDDYN%SIALPH, SILNPR=>YDDYN%SILNPR, SIRPRG=>YDDYN%SIRPRG)
!     ------------------------------------------------------------------

!*       1.    SUM GEOPOTENTIAL, COMPUTES P AND PUT IT IN PD.
!              ----------------------------------------------

CLOPER='IBOT'

IF (YDCVER%LVERTFE) THEN

  IF (YDCVER%LVFE_COMPATIBLE) CLOPER='INTG'

!$OMP PARALLEL PRIVATE(JLEV,JSPEC,ZDETAH)
!$OMP DO SCHEDULE(STATIC) 
  DO JLEV=1,KLEV
    ZDETAH=YDVETA%VFE_RDETAH(JLEV)
    DO JSPEC=1,KSPEC
      ZSPHI(JSPEC,JLEV)=-YDCST%RD*PT(JLEV,JSPEC)*SILNPR(JLEV)*ZDETAH
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ZSPHI(:,0)=0.0_JPRB
  ZSPHI(:,KLEV+1)=0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,CLOPER,'11',KSPEC,1,KSPEC,KLEV,ZSPHI,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)

!$OMP PARALLEL PRIVATE(JLEV,JSPEC)
!$OMP DO SCHEDULE(STATIC) 
  DO JLEV=1,KLEV
    DO JSPEC=1,KSPEC
      PD(JLEV,JSPEC)=ZOUT(JSPEC,JLEV-1)+PSP(JSPEC)*SIRPRG
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

ELSE

  ZSPHIX(KLEV, :)=0.0_JPRB

!$OMP PARALLEL PRIVATE(JSPEC,JLEV)
!$OMP DO SCHEDULE(STATIC)
  DO JSPEC=1,KSPEC,1
    DO JLEV=KLEV,1,-1
      ZSPHIX(JLEV-1, JSPEC)=ZSPHIX(JLEV, JSPEC)+YDCST%RD*PT(JLEV,JSPEC)*SILNPR(JLEV)
      PD(JLEV,JSPEC)=ZSPHIX(JLEV, JSPEC)+SIALPH(JLEV)*YDCST%RD*PT(JLEV,JSPEC)+PSP(JSPEC)*SIRPRG
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

ENDIF

!      -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SIGAM_SP_OPENMP',1,ZHOOK_HANDLE)

END SUBROUTINE SIGAM_SP_OPENMP
