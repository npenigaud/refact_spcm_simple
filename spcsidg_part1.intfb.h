INTERFACE
#if defined(_OPENACC)
SUBROUTINE SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, KMLOC, PSDIVP, PSPDIVP,zsdivpl,zspdivpl)
!$acc routine vector
#else
SUBROUTINE SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, KMLOC, PSDIVP, PSPDIVP)
#endif
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN       , ONLY : TDYN
IMPLICIT NONE
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN),value    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN),value    :: KMLOC
#if defined(_OPENACC)
REAL(KIND=JPRB),   INTENT(IN)    :: PSDIVP (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSPDIVP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
real(kind=JPRB),   intent(inout) :: zsdivpl(ydgeometry%yrdimv%nflevg,ydgeometry%yrdim%nsmax+1,2)
real(kind=JPRB),   intent(inout) :: zspdivpl(ydgeometry%yrdimv%nflevg,ydgeometry%yrdim%nsmax+1,2)
#else
REAL(KIND=JPRB),   INTENT(IN)    :: PSDIVP (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSPDIVP(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
#endif
END SUBROUTINE SPCSIDG_PART1

END INTERFACE
