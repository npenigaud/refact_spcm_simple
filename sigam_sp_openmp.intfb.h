INTERFACE

#if defined(_OPENACC)
SUBROUTINE SIGAM_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP,zsphi,zout)
#else
SUBROUTINE SIGAM_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP)
#endif
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : TCST
USE YOMDYN       , ONLY : TDYN
IMPLICIT NONE
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC
#if defined(_OPENACC)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD(KSPEC,klev)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KSPEC,klev)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KSPEC)
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: zsphi(KSPEC,0:klev+1)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: zout(KSPEC,0:klev)
#else
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KSPEC)
#endif
END SUBROUTINE SIGAM_SP_OPENMP

END INTERFACE
