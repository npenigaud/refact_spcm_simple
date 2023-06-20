#if defined(_OPENACC)
SUBROUTINE SPCSIDG_PART2(YDGEOMETRY,KSPEC2V,PSPDIVG,PHELP,ZSDIVPL,ZSPDIVPL,KMLOCSTA,KMLOCEND)
#else
SUBROUTINE SPCSIDG_PART2(YDGEOMETRY,KSPEC2V,PSPDIVG,PHELP,KMLOCSTA,KMLOCEND)
#endif
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCSTA
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCEND

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDIVG(KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHELP(KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
#if defined(_OPENACC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: ZSDIVPL (YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2,500)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: ZSPDIVPL(YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2,500)
#else
REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2)
REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2)
#endif

INTEGER(KIND=JPIM) :: II, IS0, ISE, JN,JCNTV,JMLOC,ji,klx
INTEGER(KIND=JPIM) :: IM, ISTA, IEND

!     ------------------------------------------------------------------

#include "mxptma.h"

!     ------------------------------------------------------------------
!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART2',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDSPGEOM=>YDGEOMETRY%YSPGEOM,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NFLEVG=>YDDIMV%NFLEVG, SCGMAP=>YDSPGEOM%SCGMAP,NSPSTAF=>YDMP%NSPSTAF)
!     ------------------------------------------------------------------

#if defined(_OPENACC)

!$ACC DATA PRESENT(zsdivpl,zspdivpl,nsmax,nflevg,ydgeometry,ydgeometry%yspgeom,ydgeometry%yrmp,ydspgeom,ydmp) 
!$ACC DATA PRESENT(YDLAP,ydlap%MYMS,ydlap%nse0l,nspstaf,scgmap,phelp,pspdivg)

!$ACC PARALLEL DEFAULT(NONE) PRIVATE(IM,ISTA,IEND,IS0,II,KLX)
!$ACC LOOP GANG COLLAPSE(3)
DO JMLOC=KMLOCSTA,KMLOCEND
 DO JCNTV=1,NFLEVG
  DO JI=1,2
   IM=YDLAP%MYMS(JMLOC)
   ISTA=NSPSTAF(IM)
   IEND=ISTA+2*(NSMAX+1-IM)-1
   KLX=NSMAX+1-IM
   IS0=YDLAP%NSE0L(JMLOC)
   II=MIN(IM,1)+1
!           ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GM**2 * DIVprim(t+dt)) .

!           Reorganisation of ZSDIVP (Back to the USSR)
    !$ACC LOOP VECTOR PRIVATE(ISE)
    DO JN=IM,NSMAX
      ISE=ISTA+2*(JN-IM)
      ZSDIVPL(JN-im+1,JCNTV,ji,JMLOC)=PSPDIVG(ISE+ji-1,JCNTV)
    ENDDO

!        ZSPDIV=(DIVprim(t+dt)) --> ZPSPDIVG=(GMBAR**2 * DIVprim(t+dt)).
if (.true.) then  !!fonction routine vector copié dans routine principale
if (ii==1 .and. ji==2) then
       !$acc loop vector
       do jn=1,klx
         zspdivpl(jn,JCNTV,2,JMLOC)=zsdivpl(jn,JCNTV,2,JMLOC)
       enddo
elseIF (KLX >= 4) THEN
      zspdivpl(1,JCNTV,JI,JMLOC) = scgmap (is0+1,1)*zsdivpl(1,JCNTV,JI,JMLOC)+scgmap(is0+1,2)*zsdivpl(2,JCNTV,JI,JMLOC)+scgmap(is0+1,3)*zsdivpl(3,JCNTV,JI,JMLOC)
      zspdivpl(2,JCNTV,JI,JMLOC) = scgmap(is0+1,2)*zsdivpl(1,JCNTV,JI,JMLOC)&
       & +scgmap(is0+2,1)*zsdivpl(2,JCNTV,JI,JMLOC)&
       & +scgmap(is0+2,2)*zsdivpl(3,JCNTV,JI,JMLOC)&
       & +scgmap(is0+2,3)*zsdivpl(4,JCNTV,JI,JMLOC)  

  !$acc loop vector 
  do jn=3,klx-2 !!jl=jn
        zspdivpl(jn,JCNTV,JI,JMLOC) = scgmap(is0+jn-2,3)*zsdivpl(jn-2,JCNTV,JI,JMLOC)&
         & +scgmap(is0+jn-1,2)*zsdivpl(jn-1,JCNTV,JI,JMLOC)&
         & +scgmap(is0+jn,1  )*zsdivpl(jn,JCNTV  ,JI,JMLOC)&
         & +scgmap(is0+jn,2  )*zsdivpl(jn+1,JCNTV,JI,JMLOC)&
         & +scgmap(is0+jn,3  )*zsdivpl(jn+2,JCNTV,JI,JMLOC)  
  ENDDO
      zspdivpl(KLX-1,JCNTV,JI,JMLOC) = scgmap(is0+KLX-3,3)*zsdivpl(KLX-3,JCNTV,JI,JMLOC)&
       & +scgmap(is0+KLX-2,2)*zsdivpl(KLX-2,JCNTV,JI,JMLOC)&
       & +scgmap (is0+KLX-1,1)*zsdivpl(KLX-1,JCNTV,JI,JMLOC)&
       & +scgmap(is0+KLX-1,2)*zsdivpl(KLX,JCNTV  ,JI,JMLOC)  
      zspdivpl(KLX,JCNTV,JI,JMLOC) = scgmap(is0+KLX-2,3)*zsdivpl(KLX-2,JCNTV,JI,JMLOC)&
       & +scgmap(is0+KLX-1,2)*zsdivpl(KLX-1,JCNTV,JI,JMLOC)&
       & +scgmap (is0+KLX,1  )*zsdivpl(KLX,JCNTV  ,JI,JMLOC)  

ELSEIF (KLX == 3) THEN
      zspdivpl(1,JCNTV,JI,JMLOC) = scgmap(is0+ 1,1)*zsdivpl(1,JCNTV,JI,JMLOC)+scgmap(is0+1,2)*zsdivpl(2,JCNTV,JI,JMLOC)+scgmap(is0+1,3)*zsdivpl(3,JCNTV,JI,JMLOC)
      zspdivpl(2,JCNTV,JI,JMLOC) = scgmap(is0+1,2)*zsdivpl(1,JCNTV,JI,JMLOC)+scgmap(is0+2,1)*zsdivpl(2,JCNTV,JI,JMLOC)+scgmap(is0+2,2)*zsdivpl(3,JCNTV,JI,JMLOC)
      zspdivpl(3,JCNTV,JI,JMLOC) = scgmap(is0+1,3)*zsdivpl(1,JCNTV,JI,JMLOC)+scgmap(is0+2,2)*zsdivpl(2,JCNTV,JI,JMLOC)+scgmap (is0+3,1)*zsdivpl(3,JCNTV,JI,JMLOC)

ELSEIF (KLX == 2) THEN
      zspdivpl(1,JCNTV,JI,JMLOC) = scgmap(is0+1,1)*zsdivpl(1,JCNTV,JI,JMLOC)+scgmap(is0+1,2)*zsdivpl(2,JCNTV,JI,JMLOC)
      zspdivpl(2,JCNTV,JI,JMLOC) = scgmap(is0+1,2)*zsdivpl(1,JCNTV,JI,JMLOC)+scgmap(is0+2,1)*zsdivpl(2,JCNTV,JI,JMLOC)

ELSEIF (KLX == 1) THEN
      zspdivpl(1,JCNTV,JI,JMLOC) = scgmap(is0+1,1)*zsdivpl(1,JCNTV,JI,JMLOC)

ENDIF
else   !!fonction dans une routine vector
IF (IM>0 .OR. JI==1) THEN
CALL MXPTMA(NSMAX+1-IM,1,1,1,NSMAX,SCGMAP(IS0+1,1),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & ZSDIVPL(1,JCNTV,JI,JMLOC),ZSPDIVPL(1,JCNTV,JI,JMLOC))  

ELSE
       !$ACC LOOP VECTOR
       DO JN=1,KLX
         ZSPDIVPL(JN,JCNTV,2,JMLOC)=ZSDIVPL(JN,JCNTV,2,JMLOC)  !!a changer pour 0
       ENDDO

ENDIF
endif   !!manière d'appeler mxptma
!           Reorganisation of ZSPDIVPL

    !$acc loop vector private(ISE)
    DO JN=IM,NSMAX
      ISE=ISTA+2*(JN-IM)
      PHELP(ISE+JI-1,JCNTV)=ZSPDIVPL(JN-IM+1,JCNTV,JI,JMLOC)
    ENDDO
   ENDDO !!ji
  ENDDO !!JCNTV
ENDDO !!JMLOC
!$ACC END PARALLEL
!$ACC END DATA
!$ACC END DATA

#else

!$OMP PARALLEL DO PRIVATE(JMLOC,IM,ISTA,IEND,IS0,II,JN,ISE,ZSDIVPL,ZSPDIVPL)
do JMLOC=KMLOCSTA,KMLOCEND
  IM=YDLAP%MYMS(JMLOC)
  ISTA=NSPSTAF(IM)
  IEND=ISTA+2*(NSMAX+1-IM)-1

  IS0=YDLAP%NSE0L(JMLOC)
  II=MIN(IM,1)+1
!           ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GM**2 * DIVprim(t+dt)) .

ZSDIVPL(:,:,:)=0.0_JPRB
ZSPDIVPL(:,:,:)=0.0_JPRB

!           Reorganisation of ZSDIVP (Back to the USSR)
  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    ZSDIVPL(JN-im+1,:,1)=PSPDIVG(ISE,:)
    ZSDIVPL(JN-im+1,:,2)=PSPDIVG(ISE+1,:)
  ENDDO

!        ZSPDIV=(DIVprim(t+dt)) --> ZPSPDIVG=(GMBAR**2 * DIVprim(t+dt)).

CALL MXPTMA(NSMAX+1-IM,NFLEVG,NFLEVG,II,nsmax,SCGMAP(IS0+1,1),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & ZSDIVPL,ZSPDIVPL)  

!           Reorganisation of ZSPDIVPL

  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
     PHELP(ISE,:)=ZSPDIVPL(JN-im+1,:,1)
     PHELP(ISE+1,:)=ZSPDIVPL(JN-im+1,:,2)
  ENDDO

ENDDO !!JMLOC
!$OMP END PARALLEL DO
#endif

END ASSOCIATE
END ASSOCIATE

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART2',1,ZHOOK_HANDLE)

END SUBROUTINE SPCSIDG_PART2

