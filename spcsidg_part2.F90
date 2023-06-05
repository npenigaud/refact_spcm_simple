#if defined(_OPENACC)
SUBROUTINE SPCSIDG_PART2(YDGEOMETRY,KSPEC2V,KMLOC,PSPDIVG,PHELP,zsdivpl,zspdivpl)
!$acc routine vector
#else
SUBROUTINE SPCSIDG_PART2(YDGEOMETRY,KSPEC2V,KMLOC,PSPDIVG,PHELP)
#endif
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
#if defined(_OPENACC)
INTEGER(KIND=JPIM),INTENT(IN),value    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN),value    :: KMLOC
#else
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOC
#endif

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDIVG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHELP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
#if defined(_OPENACC)
real(kind=JPRB)   ,intent(inout) :: zsdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)
real(kind=JPRB)   ,intent(inout) :: zspdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)
#endif

#if defined(_OPENACC)
!!REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRDIMV%NFLEVG,1:YDGEOMETRY%YRDIM%NSMAX+1,2)
!!REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRDIMV%NFLEVG,1:YDGEOMETRY%YRDIM%NSMAX+1,2)
#else
REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,ydgeometry%yrdimv%nflevg,2)
REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,ydgeometry%yrdimv%nflevg,2)
#endif

INTEGER(KIND=JPIM) :: II, IS0, ISE, JN,compteur
INTEGER(KIND=JPIM) :: IM, ISTA, IEND

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxptma.h"

!     ------------------------------------------------------------------
!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART2',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDSPGEOM=>YDGEOMETRY%YSPGEOM,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NFLEVG=>YDDIMV%NFLEVG, SCGMAP=>YDSPGEOM%SCGMAP,NSPSTAF=>YDMP%NSPSTAF)
!     ------------------------------------------------------------------

!$acc data present(zsdivpl,zspdivpl)
!$acc data present(YDLAP,ydlap%MYMS,ydlap%nse0l,nspstaf,scgmap,phelp,pspdivg)

IM=YDLAP%MYMS(KMLOC)
ISTA=NSPSTAF(IM)
IEND=ISTA+2*(NSMAX+1-IM)-1

IS0=YDLAP%NSE0L(KMLOC)
II=MIN(IM,1)+1
!           ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GM**2 * DIVprim(t+dt)) .

ZSDIVPL(:,:,:)=0.0_JPRB
ZSPDIVPL(:,:,:)=0.0_JPRB

!           Reorganisation of ZSDIVP (Back to the USSR)
#if defined(_OPENACC)
!$acc loop vector private(ISE,compteur,jn) collapse(2)
do compteur=1,nflevg
  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    ZSDIVPL(JN-im+1,compteur,1)=PSPDIVG(ISE,compteur)
    ZSDIVPL(JN-im+1,compteur,2)=PSPDIVG(ISE+1,compteur)
  enddo
ENDDO
#else
!$omp parallel do private(compteur,jn,ise) !!pas de parallelisation dans code initial
do compteur=1,nflevg
  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    ZSDIVPL(JN,compteur,1:2)=PSPDIVG(ISE:ISE+1,compteur)
  ENDDO
enddo
!$omp end parallel do
#endif

!        ZSPDIV=(DIVprim(t+dt)) --> ZPSPDIVG=(GMBAR**2 * DIVprim(t+dt)).

#if defined(_OPENACC)
CALL MXPTMA(NSMAX+1-IM,NFLEVG,NFLEVG,II,ydgeometry%yrdim%nsmax,SCGMAP(IS0+1,1),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & ZSDIVPL,ZSPDIVPL)  
#else
CALL MXPTMA(NSMAX+1-IM,NFLEVG,NFLEVG,II,SCGMAP(IS0+1,1),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & ZSDIVPL,ZSPDIVPL)  
#endif

!           Reorganisation of ZSPDIVPL

#if defined(_OPENACC)
!$acc loop vector private(ISE,compteur,jn) collapse(2)
do compteur=1,nflevg
  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    PHELP(ISE,compteur)=ZSPDIVPL(JN-im+1,compteur,1)
    PHELP(ISE+1,compteur)=ZSPDIVPL(JN-im+1,compteur,2)
  enddo
ENDDO
#else
!$omp parallel do private(compteur,jn,ise)  !!pas de parallelisation dans code initial
do compteur=1,nflevg
  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    PHELP(ISE:ISE+1,compteur)=ZSPDIVPL(JN,compteur,1:2)
  ENDDO
enddo
!$omp end parallel do
#endif

!$acc end data
!$acc end data
END ASSOCIATE
END ASSOCIATE

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART2',1,ZHOOK_HANDLE)

END SUBROUTINE SPCSIDG_PART2

