INTERFACE
#if defined(_OPENACC)
SUBROUTINE MXTURE(KLX,KVX,KVXS,KIX,tnsmax,KT,LDMT,PA,PB,PC,PY,PX)
!$acc routine vector
#else
SUBROUTINE MXTURE(KLX,KVX,KVXS,KIX,KT,LDMT,PA,PB,PC,PY,PX)
#endif
USE PARKIND1 ,ONLY : JPIM ,JPRB
INTEGER(KIND=JPIM),INTENT(IN) :: KLX
INTEGER(KIND=JPIM),INTENT(IN) :: KVX
INTEGER(KIND=JPIM),INTENT(IN) :: KVXS
INTEGER(KIND=JPIM),INTENT(IN) :: KIX
#if defined(_OPENACC)
integer(kind=jpim),intent(in) :: tnsmax
#endif
INTEGER(KIND=JPIM),INTENT(IN) :: KT
LOGICAL ,INTENT(IN) :: LDMT
REAL(KIND=JPRB) ,INTENT(IN) :: PA(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PB(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PC(KVX,KLX)
#if defined(_OPENACC)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PY(tnsmax+1,kvxs,KIX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PX(tnsmax+1,kvxs,KIX)
#else
REAL(KIND=JPRB) ,INTENT(INOUT) :: PY(KVXS,KLX,KIX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PX(KVXS,KLX,KIX)
#endif
END SUBROUTINE MXTURE
END INTERFACE
