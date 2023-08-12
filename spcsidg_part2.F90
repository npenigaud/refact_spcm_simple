SUBROUTINE SPCSIDG_PART2(YDGEOMETRY,KSPEC2V,PSPDIVG,PHELP,kmlocsta,kmlocend)
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

#else
REAL(KIND=JPRB) :: ZSDIVPL (ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)
REAL(KIND=JPRB) :: ZSPDIVPL(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)
#endif

INTEGER(KIND=JPIM) :: II, IS0, ISE, JN,compteur,jmloc,ji,klx,longueur
INTEGER(KIND=JPIM) :: IM, ISTA, IEND,jv,jl,kvx

!     ------------------------------------------------------------------

#include "mxptma.h"

!     ------------------------------------------------------------------
!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART2',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDSPGEOM=>YDGEOMETRY%YSPGEOM,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NFLEVG=>YDDIMV%NFLEVG,NSPSTAF=>YDMP%NSPSTAF)
!     ------------------------------------------------------------------

#if defined(_OPENACC)

!$acc data present(ydgeometry%yrdim%nsmax,ydgeometry%yrdimv%nflevg,ydgeometry,ydgeometry%yspgeom,ydgeometry%yrmp) 
!$acc data present(ydgeometry%yrlap%MYMS,ydgeometry%yrlap%nse0l,ydgeometry%yrmp%nspstaf,ydgeometry%yspgeom%scgmap,phelp,pspdivg)
is0=ydgeometry%yrlap%nse0L(kmlocsta)
longueur=ydgeometry%yrlap%nse0L(kmlocend)-is0+ydgeometry%yrdim%nsmax-ydgeometry%yrlap%myms(kmlocend)+1
if (.FALSE.) then
call mxptma(longueur,nflevg,nflevg,ydgeometry%yspgeom%scgmap(is0+1,1),&
&ydgeometry%yspgeom%scgmap(is0+1,2),ydgeometry%yspgeom%scgmap(is0+1,3),&
&ydgeometry%yspgeom%scgmap(is0+1,2),ydgeometry%yspgeom%scgmap(is0+1,3),pspdivg,phelp)   
else
klx=longueur
kvx=nflevg
  !$acc parallel private(jv,ji) default(none) async(1)
  !$acc loop gang
  do jv=1,kvx 
    !$acc loop vector
    do ji=-1,0   
      phelp(2+ji,jv) = ydgeometry%yspgeom%scgmap(is0+1,1)*pspdivg(2+ji,jv)&
       &+ydgeometry%yspgeom%scgmap(is0+1,2)*pspdivg(4+ji,jv)&
       &+ydgeometry%yspgeom%scgmap(is0+1,3)*pspdivg(6+ji,jv)
      phelp(4+ji,jv) = ydgeometry%yspgeom%scgmap(is0+1,2)*pspdivg(2+ji,jv)&
       & +ydgeometry%yspgeom%scgmap(is0+2,1)*pspdivg(4+ji,jv)&
       & +ydgeometry%yspgeom%scgmap(is0+2,2)*pspdivg(6+ji,jv)&
       & +ydgeometry%yspgeom%scgmap(is0+2,3)*pspdivg(8+ji,jv)  
    ENDDO
  enddo
  !$acc end parallel 

  !$acc parallel private(jl,jv,ji) default(none) async(2)
  !$acc loop gang vector collapse(3) private(jl,jv,ji)
  do jv=1,kvx
    do jl=3,klx-2        
     do ji=-1,0
!    !$acc cache(px(jl-2:jl+2,:),pbi(jl-2),pci(jl-1),pa(jl),pbs(jl),pcs(jl))
        phelp(2*JL+ji,jv) = ydgeometry%yspgeom%scgmap(is0+jl-2,3)*pspdivg(2*(JL-2)+ji,jv)&
         & +ydgeometry%yspgeom%scgmap(is0+jl-1,2)*pspdivg(2*(JL-1)+ji,jv)&
         & +ydgeometry%yspgeom%scgmap(is0+jl,1)*pspdivg(2*JL+ji,jv  )&
         & +ydgeometry%yspgeom%scgmap(is0+jl,2)*pspdivg(2*(JL+1)+ji,jv)&
         & +ydgeometry%yspgeom%scgmap(is0+jl,3)*pspdivg(2*(JL+2)+ji,jv)
!        if (jv==10) then
!          print *,"PX indice",pspdivg(2*jl+ji,jv),2*jl+ji
!          print *,"PY indice",phelp(2*jl+ji,jv),2*jl+ji
!          print *,"coeffs",ydgeometry%yspgeom%scgmap(is0+jl-2,3),ydgeometry%yspgeom%scgmap(is0+jl-2,2),ydgeometry%yspgeom%scgmap(is0+jl,1),&
 !          &ydgeometry%yspgeom%scgmap(is0+jl,3),ydgeometry%yspgeom%scgmap(is0+jl,3)
  !      endif

 !        phelp(2*JL-1:2*jl,jv) = ydgeometry%yspgeom%scgmap(is0+jl-2,3)*pspdivg(2*(JL-2)-1:2*(jl-2),jv)&
 !        & +ydgeometry%yspgeom%scgmap(is0+jl-1,2)*pspdivg(2*(JL-1)-1:2*(jl-1),jv)&
 !        & +ydgeometry%yspgeom%scgmap(is0+jl,1)*pspdivg(2*JL-1:2*jl,jv  )&
 !        & +ydgeometry%yspgeom%scgmap(is0+jl,2)*pspdivg(2*(JL+1)-1:2*(jl+1),jv)&
 !        & +ydgeometry%yspgeom%scgmap(is0+jl,3)*pspdivg(2*(JL+2)-1:2*(jl+2),jv)
 !       if (jv==10) then
 !         print *,"PX",pspdivg(2*jl-1,jv),pspdivg(2*jl,jv)
 !         print *,"PY",phelp(2*jl-1,jv),phelp(2*jl,jv)
 !         print *,"coeffs",ydgeometry%yspgeom%scgmap(is0+jl-2,3),ydgeometry%yspgeom%scgmap(is0+jl-2,2),ydgeometry%yspgeom%scgmap(is0+jl,1),&
 !          &ydgeometry%yspgeom%scgmap(is0+jl,3),ydgeometry%yspgeom%scgmap(is0+jl,3)
 !       endif
       
       enddo

     ENDDO
  ENDDO
  !$acc end parallel
  !$acc parallel private(jv,ji) default(none) async(3)
  !$acc loop gang
  do jv=1,kvx
    !$acc loop vector  
    do ji=-1,0
      phelp(2*(KLX-1)+ji,jv) = ydgeometry%yspgeom%scgmap(is0+klx-3,3)*pspdivg(2*(KLX-3)+ji,jv)&
       & +ydgeometry%yspgeom%scgmap(is0+klx-2,2)*pspdivg(2*(KLX-2)+ji,jv)&
       & +ydgeometry%yspgeom%scgmap(is0+klx-1,1)*pspdivg(2*(KLX-1)+ji,jv)&
       & +ydgeometry%yspgeom%scgmap(is0+klx-1,2)*pspdivg(2*KLX+ji,jv)
      phelp(2*KLX+ji,jv) = ydgeometry%yspgeom%scgmap(is0+klx-2,3)*pspdivg(2*(KLX-2)+ji,jv)&
       & +ydgeometry%yspgeom%scgmap(is0+klx-1,2)*pspdivg(2*(KLX-1)+ji,jv)&
       & +ydgeometry%yspgeom%scgmap(is0+klx,1)*pspdivg(2*KLX+ji,jv  )  
    enddo
   ENDDO
  !$acc end parallel
  !$acc wait

endif

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

CALL MXPTMA(NSMAX+1-IM,NFLEVG,NFLEVG,II,nsmax,ydgeometry%yspgeom%scgmap(IS0+1,1),&
 & ydgeometry%yspgeom%scgmap(IS0+1,2),ydgeometry%yspgeom%scgmap(IS0+1,3),&
 & ydgeometry%yspgeom%scgmap(IS0+1,2),ydgeometry%yspgeom%scgmap(IS0+1,3),&
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
if (.FALSE.) then
!$acc update host(phelp)
do jmloc=kmlocsta,kmlocend
  IM=YDLAP%MYMS(jmloc)
  ISTA=NSPSTAF(IM)
  IEND=ISTA+2*(NSMAX+1-IM)-1

  IS0=YDLAP%NSE0L(jmloc)
  II=MIN(IM,1)+1
  print *,"VALEUR DE JMLOC",jmloc
  print *,"plage spectrale pour niveau 1"
  print *,phelp(ista:iend,1)
  print *,"plage spectrale pour niveau 10"
  print *,phelp(ista:iend,10)
  print *,"plage spectrale pour niveau 105"
  print *,phelp(ista:iend,105)
enddo
endif
!$acc end data
!$acc end data
END ASSOCIATE
END ASSOCIATE
!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART2',1,ZHOOK_HANDLE)

END SUBROUTINE SPCSIDG_PART2

