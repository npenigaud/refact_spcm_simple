#if defined(_OPENACC)
SUBROUTINE SPCSIDG_PART2(YDGEOMETRY,KSPEC2V,PSPDIVG,PHELP,zsdivpl,zspdivpl,kmlocsta,kmlocend)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY

INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN)    :: kmlocsta
integer(kind=jpim),intent(in)    :: kmlocend

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDIVG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHELP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
real(kind=JPRB)   ,intent(inout) :: zsdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)
real(kind=JPRB)   ,intent(inout) :: zspdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2)

INTEGER(KIND=JPIM) :: II, IS0, ISE, JN,compteur
INTEGER(KIND=JPIM) :: IM, ISTA, IEND,kmloc

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDSPGEOM=>YDGEOMETRY%YSPGEOM,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NFLEVG=>YDDIMV%NFLEVG, SCGMAP=>YDSPGEOM%SCGMAP,NSPSTAF=>YDMP%NSPSTAF)

!$acc data present(zsdivpl,zspdivpl)
!$acc data present(YDLAP,ydlap%MYMS,ydlap%nse0l,nspstaf,scgmap,phelp,pspdivg)

kix=2
klx=nsmax
kvx=nflevl
kvxs=nflevl
!$acc parallel default(none)
!$acc loop gang collapse(3) private( )
do kmloc=kmlocsta,kmlocend
  do compteur=1,nflevl
   do ji=1,kix

     IM=YDLAP%MYMS(KMLOC) !!! traiter le cas km=0
     ISTA=NSPSTAF(IM)
     IEND=ISTA+2*(NSMAX+1-IM)-1

     IS0=YDLAP%NSE0L(KMLOC)
     II=MIN(IM,1)+1

     !! decalage=()* !!parametres : scgmap(is0+1,1 à 3)

     !$acc loop vector private(ise)
     do jn=1,iend
       ise=ista+2*(jn-1)
       entree3(jn)=pspdivg(ise+ji,compteur)
     enddo

     







enddo !!kmloc

!$acc end data
!$acc end data
END ASSOCIATE
END ASSOCIATE

END SUBROUTINE SPCSIDG_PART2




#else

SUBROUTINE SPCSIDG_PART2(YDGEOMETRY,KSPEC2V,PSPDIVG,PHELP,kmlocsta,kmlocend)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
integer(kind=jpim),intent(in)    :: kmlocsta
integer(kind=jpim),intent(in)    :: kmlocend

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDIVG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PHELP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
!REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,ydgeometry%yrdimv%nflevg,2)
!REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,ydgeometry%yrdimv%nflevg,2)
REAL(KIND=JPRB),allocatable :: ZSDIVPL (:,:,:)
REAL(KIND=JPRB),allocatable :: ZSPDIVPL(:,:,:)


INTEGER(KIND=JPIM) :: II, IS0, ISE, JN,compteur
INTEGER(KIND=JPIM) :: IM, ISTA, IEND,kmloc

!     ------------------------------------------------------------------

#include "mxptma.h"
!     ------------------------------------------------------------------

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDSPGEOM=>YDGEOMETRY%YSPGEOM,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NFLEVG=>YDDIMV%NFLEVG, SCGMAP=>YDSPGEOM%SCGMAP,NSPSTAF=>YDMP%NSPSTAF)
!     ------------------------------------------------------------------
print *,"taille scgmap",size(scgmap,1)
print *,"taille scgmap",size(scgmap,2)
do kmloc=kmlocsta,kmlocend
allocate(zsdivpl(YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,ydgeometry%yrdimv%nflevg,2))
allocate(ZSPDIVPL(YDGEOMETRY%YRLAP%MYMS(KMLOC):YDGEOMETRY%YRDIM%NSMAX,ydgeometry%yrdimv%nflevg,2))



 IM=YDLAP%MYMS(KMLOC)
 ISTA=NSPSTAF(IM)
 IEND=ISTA+2*(NSMAX+1-IM)-1

 IS0=YDLAP%NSE0L(KMLOC)
 print *,"is0",is0
 II=MIN(IM,1)+1
 !           ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GM**2 * DIVprim(t+dt)) .

 ZSDIVPL(:,:,:)=0.0_JPRB
 ZSPDIVPL(:,:,:)=0.0_JPRB

 !           Reorganisation of ZSDIVP (Back to the USSR)
 !$omp parallel do private(compteur,jn,ise) !!pas de parallelisation dans code initial
 do compteur=1,nflevg
  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    ZSDIVPL(JN,compteur,1:2)=PSPDIVG(ISE:ISE+1,compteur)
  ENDDO
 enddo
 !$omp end parallel do

 !        ZSPDIV=(DIVprim(t+dt)) --> ZPSPDIVG=(GMBAR**2 * DIVprim(t+dt)).

 CALL MXPTMA(NSMAX+1-IM,NFLEVG,NFLEVG,II,SCGMAP(IS0+1,1),&
  & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
  & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
  & ZSDIVPL,ZSPDIVPL)  

 !           Reorganisation of ZSPDIVPL

 !$omp parallel do private(compteur,jn,ise)  !!pas de parallelisation dans code initial
 do compteur=1,nflevg
  DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    PHELP(ISE:ISE+1,compteur)=ZSPDIVPL(JN,compteur,1:2)
  ENDDO
 enddo
 !$omp end parallel do
deallocate(zsdivpl)
deallocate(zspdivpl)
enddo !!kmloc

END ASSOCIATE
END ASSOCIATE

END SUBROUTINE SPCSIDG_PART2
#endif
