#if defined(_OPENACC)
SUBROUTINE SIGAM_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT,PSP,zsphi,zout)
#else
SUBROUTINE SIGAM_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP)
#endif

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

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD(KSPEC,klev)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KSPEC,klev)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KSPEC)

!     ------------------------------------------------------------------
#if defined(_OPENACC)
REAL(KIND=JPRB),intent(inout) :: ZSPHI(KSPEC,0:KLEV+1)
REAL(KIND=JPRB),intent(inout) :: ZOUT(KSPEC,0:KLEV)
#else
REAL(KIND=JPRB) :: ZSPHI(KSPEC,0:KLEV+1)
REAL(KIND=JPRB) :: ZOUT(KSPEC,0:KLEV)
#endif


REAL(KIND=JPRB) :: ZSPHIX(0:KLEV, KSPEC)
INTEGER(KIND=JPIM) :: JLEV, JSPEC
CHARACTER(LEN=4):: CLOPER
REAL(KIND=JPRB) :: ZDETAH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE,zhook_handle2

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIGAM_SP_OPENMP',0,ZHOOK_HANDLE)

ASSOCIATE( YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(SIALPH=>YDDYN%SIALPH, SILNPR=>YDDYN%SILNPR, SIRPRG=>YDDYN%SIRPRG)
!     ------------------------------------------------------------------

!*       1.    SUM GEOPOTENTIAL, COMPUTES P AND PUT IT IN PD.
!              ----------------------------------------------
CLOPER='IBOT'

IF (YDCVER%LVERTFE) THEN

  IF (YDCVER%LVFE_COMPATIBLE) CLOPER='INTG'
!$acc data present(pt,psp,pd,klev)
!$acc data present(ydgeometry%YRVETA%vfe_rdetah,yddyn,YDDYN%SILNPR,YDDYN%SIALPH,YDDYN%SIRPRG,ydgeometry%YRVETA,ydcst)
!$acc data present(zsphi,zout)

if (lhook) call dr_hook('SIGAM_transpose1',0,zhook_handle2)

#if defined(_OPENACC)
!$ACC PARALLEL PRIVATE(JLEV,JSPEC,ZDETAH) vector_length(128) DEFAULT(NONE)
!$ACC LOOP GANG VECTOR collapse(2)
#else
!$OMP PARALLEL PRIVATE(JLEV,JSPEC,ZDETAH)
!$OMP DO SCHEDULE(STATIC) 
#endif
DO JLEV=1,KLEV

!      ZDETAH=-YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)*YDCST%RD*YDDYN%SILNPR(JLEV)
!  !$acc loop vector
  DO JSPEC=1,KSPEC
!    DO JLEV=1,KLEV
      !ZDETAH=-YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)*YDCST%RD*YDDYN%SILNPR(JLEV)
      ZSPHI(JSPEC,JLEV)=-PT(JSPEC,JLEV)*YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)*YDCST%RD*YDDYN%SILNPR(JLEV)!!ZDETAH
    ENDDO
  ENDDO
#if defined(_OPENACC)
!$ACC END PARALLEL
#else
!$OMP END DO
!$OMP END PARALLEL
#endif
if (lhook) call dr_hook('SIGAM_transpose1',1,zhook_handle2)

if (lhook) call dr_hook('SIGAM_cond_lim',0,zhook_handle2)
#if defined(_OPENACC)
!$acc parallel private(jspec) default(none)
!$acc loop gang vector
#else
!$omp parallel do private(jspec)  !!dans le code initial fait par :
#endif
do jspec=1,kspec
  ZSPHI(jspec,0)=0.0_JPRB
  ZSPHI(jspec,KLEV+1)=0.0_JPRB
enddo
#if defined(_OPENACC)
!$acc end parallel
#else
!$omp end parallel do
#endif
if (lhook) call dr_hook('SIGAM_cond_lim',1,zhook_handle2)

  CALL VERDISINT(ydgeometry%YrVFE,YDCVER,CLOPER,'11',KSPEC,1,KSPEC,KLEV,ZSPHI,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)
if (lhook) call dr_hook('SIGAM_transpose2',0,zhook_handle2)
#if defined(_OPENACC)
!$ACC PARALLEL PRIVATE(JLEV,JSPEC) vector_length(128) DEFAULT(NONE)
!$ACC LOOP gang vector collapse(2)
#else
!$OMP PARALLEL PRIVATE(JLEV,JSPEC)
!$OMP DO SCHEDULE(STATIC) 
#endif
!do jspec=1,kspec
!  zdetah=PSP(JSPEC)*YDDYN%SIRPRG
!  !$acc loop vector

  DO JLEV=1,KLEV
!    !$acc loop vector
    DO JSPEC=1,KSPEC
      PD(JSPEC,JLEV)=ZOUT(JSPEC,JLEV-1)+PSP(JSPEC)*YDDYN%SIRPRG
    ENDDO
  ENDDO
#if defined(_OPENACC)
!$ACC END PARALLEL
#else
!$OMP END DO
!$OMP END PARALLEL
#endif

if (lhook) call dr_hook('SIGAM_transpose2',1,zhook_handle2)
!$acc end data
!$acc end data
!$acc end data

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
