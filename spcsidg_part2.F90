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
INTEGER(KIND=JPIM),INTENT(IN),value    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN),value    :: KMLOC
#if defined(_OPENACC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDIVG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHELP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
real(kind=JPRB)   ,intent(inout) :: zsdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)
real(kind=JPRB)   ,intent(inout) :: zspdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)
#else
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHELP(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
#endif

#if defined(_OPENACC)
!!REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRDIMV%NFLEVG,1:YDGEOMETRY%YRDIM%NSMAX+1,2)
!!REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRDIMV%NFLEVG,1:YDGEOMETRY%YRDIM%NSMAX+1,2)
#else
REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,2)
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
!$acc loop vector private(ISE,compteur)
DO JN=IM,NSMAX
  ISE=ISTA+2*(JN-IM)
  do compteur=1,nflevg
    ZSDIVPL(JN-im+1,compteur,1:2)=PSPDIVG(ISE:ISE+1,compteur)
  enddo
ENDDO
#else
DO JN=IM,NSMAX
  ISE=ISTA+2*(JN-IM)
  ZSDIVPL(:,JN,1:2)=PSPDIVG(:,ISE:ISE+1)
ENDDO
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
!$acc loop vector private(ISE,compteur)
DO JN=IM,NSMAX
  ISE=ISTA+2*(JN-IM)
  do compteur=1,nflevg
    PHELP(ISE:ISE+1,compteur)=ZSPDIVPL(JN-im+1,compteur,1:2)
  enddo
ENDDO
#else
DO JN=IM,NSMAX
  ISE=ISTA+2*(JN-IM)
  PHELP(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)
ENDDO
#endif

!$acc end data
!$acc end data
END ASSOCIATE
END ASSOCIATE

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART2',1,ZHOOK_HANDLE)

END SUBROUTINE SPCSIDG_PART2

