INTERFACE
#if defined(_OPENACC)
SUBROUTINE SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, PSDIVP, PSPDIVP,taillec,&
  &pas,pbs,pcs,entree,sortie,param_mxture,kmlocsta,kmlocend)
!$acc routine vector
#else
SUBROUTINE SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, PSDIVP, PSPDIVP,kmlocsta,kmlocend)
#endif
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN       , ONLY : TDYN
IMPLICIT NONE
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCsta
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCend

REAL(KIND=JPRB),   INTENT(IN)    :: PSDIVP (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSPDIVP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
#if defined(_OPENACC)
integer(kind=jpim)               :: taillec
real(kind=jprb),   intent(in)    :: param_mxture(:,:,:)
real(kind=jprb),   intent(inout) :: pas(2001)
real(kind=jprb),   intent(inout) :: pbs(2001)
real(kind=jprb),   intent(inout) :: pcs(2001)
real(kind=jprb),   intent(inout) :: entree(2001)
real(kind=jprb),   intent(inout) :: sortie(2001)
#endif
END SUBROUTINE SPCSIDG_PART1

END INTERFACE
