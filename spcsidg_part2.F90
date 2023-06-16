#if defined(_OPENACC)
SUBROUTINE SPCSIDG_PART2(YDGEOMETRY,KSPEC2V,PSPDIVG,PHELP,zsdivpl,zspdivpl,kmlocsta,kmlocend)
#else
SUBROUTINE SPCSIDG_PART2(YDGEOMETRY,KSPEC2V,PSPDIVG,PHELP,kmlocsta,kmlocend)
#endif
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN)    :: kmlocsta
INTEGER(KIND=JPIM),INTENT(IN)    :: kmlocend

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDIVG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHELP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
#if defined(_OPENACC)
real(kind=JPRB)   ,intent(inout) :: zsdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2,130)
real(kind=JPRB)   ,intent(inout) :: zspdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2,130)
#else
REAL(KIND=JPRB) :: ZSDIVPL (ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)
REAL(KIND=JPRB) :: ZSPDIVPL(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)
#endif

INTEGER(KIND=JPIM) :: II, IS0, ISE, JN,compteur,jmloc,ji,klx
INTEGER(KIND=JPIM) :: IM, ISTA, IEND

!     ------------------------------------------------------------------

#include "mxptma.h"

!     ------------------------------------------------------------------
!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART2',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDSPGEOM=>YDGEOMETRY%YSPGEOM,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NFLEVG=>YDDIMV%NFLEVG, SCGMAP=>YDSPGEOM%SCGMAP,NSPSTAF=>YDMP%NSPSTAF)
!     ------------------------------------------------------------------

#if defined(_OPENACC)

!$acc data present(zsdivpl,zspdivpl,nsmax,nflevg,ydgeometry,ydgeometry%yspgeom,ydgeometry%yrmp,ydspgeom,ydmp) 
!$acc data present(YDLAP,ydlap%MYMS,ydlap%nse0l,nspstaf,scgmap,phelp,pspdivg)

!$acc parallel default(none) private(im,ista,iend,is0,ii,klx)
!$acc loop gang collapse(3)
do jmloc=kmlocsta,kmlocend
 do compteur=1,nflevg
  do ji=1,2
   IM=YDLAP%MYMS(jmloc)
   ISTA=NSPSTAF(IM)
   IEND=ISTA+2*(NSMAX+1-IM)-1
   klx=nsmax+1-im
   IS0=YDLAP%NSE0L(jmloc)
   II=MIN(IM,1)+1
!           ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GM**2 * DIVprim(t+dt)) .

!           Reorganisation of ZSDIVP (Back to the USSR)
    !$acc loop vector private(ISE)
    DO JN=IM,NSMAX
      ise=ista+2*(jn-im)
      ZSDIVPL(JN-im+1,compteur,ji,jmloc)=PSPDIVG(ISE+ji-1,compteur)
    ENDDO

!        ZSPDIV=(DIVprim(t+dt)) --> ZPSPDIVG=(GMBAR**2 * DIVprim(t+dt)).
if (.true.) then  !!fonction routine vector copié dans routine principale
if (ii==1 .and. ji==2) then
       !$acc loop vector
       do jn=1,klx
         zspdivpl(jn,compteur,2,jmloc)=zsdivpl(jn,compteur,2,jmloc)
       enddo
elseIF (KLX >= 4) THEN
      zspdivpl(1,compteur,JI,jmloc) = scgmap (is0+1,1)*zsdivpl(1,compteur,JI,jmloc)+scgmap(is0+1,2)*zsdivpl(2,compteur,JI,jmloc)+scgmap(is0+1,3)*zsdivpl(3,compteur,JI,jmloc)
      zspdivpl(2,compteur,JI,jmloc) = scgmap(is0+1,2)*zsdivpl(1,compteur,JI,jmloc)&
       & +scgmap(is0+2,1)*zsdivpl(2,compteur,JI,jmloc)&
       & +scgmap(is0+2,2)*zsdivpl(3,compteur,JI,jmloc)&
       & +scgmap(is0+2,3)*zsdivpl(4,compteur,JI,jmloc)  

  !$acc loop vector 
  do jn=3,klx-2 !!jl=jn
        zspdivpl(jn,compteur,JI,jmloc) = scgmap(is0+jn-2,3)*zsdivpl(jn-2,compteur,JI,jmloc)&
         & +scgmap(is0+jn-1,2)*zsdivpl(jn-1,compteur,JI,jmloc)&
         & +scgmap(is0+jn,1  )*zsdivpl(jn,compteur  ,JI,jmloc)&
         & +scgmap(is0+jn,2  )*zsdivpl(jn+1,compteur,JI,jmloc)&
         & +scgmap(is0+jn,3  )*zsdivpl(jn+2,compteur,JI,jmloc)  
  ENDDO
      zspdivpl(KLX-1,compteur,JI,jmloc) = scgmap(is0+KLX-3,3)*zsdivpl(KLX-3,compteur,JI,jmloc)&
       & +scgmap(is0+KLX-2,2)*zsdivpl(KLX-2,compteur,JI,jmloc)&
       & +scgmap (is0+KLX-1,1)*zsdivpl(KLX-1,compteur,JI,jmloc)&
       & +scgmap(is0+KLX-1,2)*zsdivpl(KLX,compteur  ,JI,jmloc)  
      zspdivpl(KLX,compteur,JI,jmloc) = scgmap(is0+KLX-2,3)*zsdivpl(KLX-2,compteur,JI,jmloc)&
       & +scgmap(is0+KLX-1,2)*zsdivpl(KLX-1,compteur,JI,jmloc)&
       & +scgmap (is0+KLX,1  )*zsdivpl(KLX,compteur  ,JI,jmloc)  

ELSEIF (KLX == 3) THEN
      zspdivpl(1,compteur,JI,jmloc) = scgmap(is0+ 1,1)*zsdivpl(1,compteur,JI,jmloc)+scgmap(is0+1,2)*zsdivpl(2,compteur,JI,jmloc)+scgmap(is0+1,3)*zsdivpl(3,compteur,JI,jmloc)
      zspdivpl(2,compteur,JI,jmloc) = scgmap(is0+1,2)*zsdivpl(1,compteur,JI,jmloc)+scgmap(is0+2,1)*zsdivpl(2,compteur,JI,jmloc)+scgmap(is0+2,2)*zsdivpl(3,compteur,JI,jmloc)
      zspdivpl(3,compteur,JI,jmloc) = scgmap(is0+1,3)*zsdivpl(1,compteur,JI,jmloc)+scgmap(is0+2,2)*zsdivpl(2,compteur,JI,jmloc)+scgmap (is0+3,1)*zsdivpl(3,compteur,JI,jmloc)

ELSEIF (KLX == 2) THEN
      zspdivpl(1,compteur,JI,jmloc) = scgmap(is0+1,1)*zsdivpl(1,compteur,JI,jmloc)+scgmap(is0+1,2)*zsdivpl(2,compteur,JI,jmloc)
      zspdivpl(2,compteur,JI,jmloc) = scgmap(is0+1,2)*zsdivpl(1,compteur,JI,jmloc)+scgmap(is0+2,1)*zsdivpl(2,compteur,JI,jmloc)

ELSEIF (KLX == 1) THEN
      zspdivpl(1,compteur,JI,jmloc) = scgmap(is0+1,1)*zsdivpl(1,compteur,JI,jmloc)

ENDIF
else   !!fonction dans une routine vector
if (im>0 .or. ji==1) then
CALL MXPTMA(NSMAX+1-IM,1,1,1,nsmax,SCGMAP(IS0+1,1),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
 & ZSDIVPL(1,compteur,ji,jmloc),ZSPDIVPL(1,compteur,ji,jmloc))  

else
       !$acc loop vector
       do jn=1,klx
         zspdivpl(jn,compteur,2,jmloc)=zsdivpl(jn,compteur,2,jmloc)
       enddo

endif
endif   !!manière d'appeler mxptma
!           Reorganisation of ZSPDIVPL

    !$acc loop vector private(ISE)
    DO JN=IM,NSMAX
      ise=ista+2*(jn-im)
      PHELP(ISE+ji-1,compteur)=ZSPDIVPL(JN-im+1,compteur,ji,jmloc)
    ENDDO
   enddo !!ji
  enddo !!compteur
enddo !!jmloc
!$acc end parallel
!$acc end data
!$acc end data

#else

!$omp parallel do private(jmloc,im,ista,iend,is0,ii,jn,ise,zsdivpl,zspdivpl)
do jmloc=kmlocsta,kmlocend
  IM=YDLAP%MYMS(jmloc)
  ISTA=NSPSTAF(IM)
  IEND=ISTA+2*(NSMAX+1-IM)-1

  IS0=YDLAP%NSE0L(jmloc)
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

enddo !!jmloc
!$omp end parallel do
#endif

END ASSOCIATE
END ASSOCIATE

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART2',1,ZHOOK_HANDLE)

END SUBROUTINE SPCSIDG_PART2

