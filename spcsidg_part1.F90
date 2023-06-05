#if defined(_OPENACC)
SUBROUTINE SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, KMLOC, PSDIVP,PSPDIVP,zsdivpl,zspdivpl)
!$acc routine vector
#else
SUBROUTINE SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, KMLOC, PSDIVP, PSPDIVP)
#endif
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN       , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
#if defined(_OPENACC)
INTEGER(KIND=JPIM),INTENT(IN),value    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN),value    :: KMLOC
#else
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOC
#endif

REAL(KIND=JPRB),   INTENT(IN)    :: PSDIVP (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSPDIVP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
#if defined(_OPENACC)
REAL(KIND=JPRB),   intent(inout) :: ZSDIVPL (1:YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2)
REAL(KIND=JPRB),   intent(inout) :: ZSPDIVPL(1:YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2)
#endif


!     ------------------------------------------------------------------
#if defined(_OPENACC)
!!REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRDIMV%NFLEVG,1:YDGEOMETRY%YRDIM%NSMAX+1,2)
!!REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRDIMV%NFLEVG,1:YDGEOMETRY%YRDIM%NSMAX+1,2)
#else
REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,2)
#endif

INTEGER(KIND=JPIM) :: II, IS0, IS02, ISE, JN,compteur
INTEGER(KIND=JPIM) :: IM, ISTA, IEND

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxture.h"
#include "mxturs.h"

!     ------------------------------------------------------------------

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART1',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX,NFLEVG=>YDDIMV%NFLEVG,SIHEG=>YDDYN%SIHEG,SIHEG2=>YDDYN%SIHEG2,NSPSTAF=>YDMP%NSPSTAF)

!             Inversion of two tridiagonal systems (Helmholtz equation)
!                --> (SIMI*DIVprim(t+dt)).

!             Reorganisation of divergence

!$acc data present(psdivp,pspdivp,zsdivpl,zspdivpl)
!$acc data present(YDLAP,YDLAP%MYMS,NSPSTAF,SIHEG)
IM=YDLAP%MYMS(KMLOC)
ISTA=NSPSTAF(IM)
IEND=ISTA+2*(NSMAX+1-IM)-1

IS0=YDLAP%NSE0L(KMLOC)
IS02=0
II=MIN(IM,1)+1
ZSDIVPL(:,:,:)=0.0_JPRB
ZSPDIVPL(:,:,:)=0.0_JPRB

#if defined(_OPENACC)
!$acc loop vector private(ise,compteur,jn) collapse(2)
do compteur=1,nflevg
  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    ZSDIVPL(JN-im+1,compteur,1)=PSDIVP(ISE,compteur)
    ZSDIVPL(JN-im+1,compteur,2)=PSDIVP(ISE+1,compteur)
  enddo
ENDDO
#else
!$omp parallel do private(compteur,jn,ise) !!pas de parallelisation dans code initial, pas d inversion pour le moment
DO JN=IM,NSMAX
  do compteur=1,nflevg
    ISE=ISTA+2*(JN-IM)
    ZSDIVPL(compteur,JN,1:2)=PSDIVP(ISE:ISE+1,compteur)
  enddo
ENDDO
!$omp end parallel do
#endif

IF (IM > 0) THEN

  !               Inversion of a symmetric matrix.
#if defined(_OPENACC)
  CALL MXTURS(NSMAX+1-IM,NFLEVG,NFLEVG,II,nsmax,&
   & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
   & ZSDIVPL,ZSPDIVPL) 
#else
  CALL MXTURS(NSMAX+1-IM,NFLEVG,NFLEVG,II,&
   & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
   & ZSDIVPL,ZSPDIVPL)  
#endif
ELSE

  !               Inversion of a non-symmetric matrix.
#if defined(_OPENACC)
  CALL MXTURE(NSMAX+1-IM,NFLEVG,NFLEVG,II,nsmax,-2,.TRUE.,&
   & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
   & ZSDIVPL,ZSPDIVPL)  
  CALL MXTURE(NSMAX+1-IM,NFLEVG,NFLEVG,II,nsmax,3,.FALSE.,&
   & SIHEG(1,IS0+1,1),SIHEG2(1,IS02+1,2),&
   & SIHEG2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)
#else
  CALL MXTURE(NSMAX+1-IM,NFLEVG,NFLEVG,II,-2,.TRUE.,&
   & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
   & ZSDIVPL,ZSPDIVPL)  
  CALL MXTURE(NSMAX+1-IM,NFLEVG,NFLEVG,II,3,.FALSE.,&
   & SIHEG(1,IS0+1,1),SIHEG2(1,IS02+1,2),&
   & SIHEG2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)
#endif  
ENDIF

#if defined(_OPENACC)
!$acc loop vector private(ise,compteur,jn) collapse(2)
do compteur=1,nflevg
  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    PSPDIVP(ISE,compteur)=ZSPDIVPL(JN-im+1,compteur,1)
    PSPDIVP(ISE+1,compteur)=ZSPDIVPL(JN-im+1,compteur,2)
  enddo
ENDDO
#else
!$omp parallel do private(jn,compteur,ise) collapse(2)
DO JN=IM,NSMAX
  do compteur=1,nflevg
    ISE=ISTA+2*(JN-IM)
    PSPDIVP(ISE:ISE+1,compteur)=ZSPDIVPL(compteur,JN,1:2)
  enddo
ENDDO
!$omp end parallel do
#endif
!$acc end data
!$acc end data
END ASSOCIATE
END ASSOCIATE

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART1',1,ZHOOK_HANDLE)
END SUBROUTINE SPCSIDG_PART1

