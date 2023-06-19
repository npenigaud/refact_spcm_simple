#if defined(_OPENACC)
SUBROUTINE SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, PSDIVP,PSPDIVP,&
  &zsdivpl,zspdivpl,param_mxture, kmlocsta, kmlocend)

#else
SUBROUTINE SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, PSDIVP, PSPDIVP ,kmlocsta, kmlocend)
#endif
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN       , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCSTA
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCEND


REAL(KIND=JPRB),   INTENT(IN)    :: PSDIVP (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSPDIVP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
#if defined(_OPENACC)
INTEGER(KIND=JPIM)               :: taillec
INTEGER(KIND=JPIM), PARAMETER    :: tbloc=62
INTEGER(KIND=JPIM), PARAMETER     :: bloclev=8
REAL(KIND=JPRB),   INTENT(IN)    :: param_mxture(:,:,:)
REAL(KIND=JPRB) :: pas(tbloc+3,bloclev)
REAL(KIND=JPRB) :: pbs(tbloc+3,bloclev)
REAL(KIND=JPRB) :: pcs(tbloc+3,bloclev)
REAL(KIND=JPRB) :: entree(tbloc+3,bloclev,2)
REAL(KIND=JPRB) :: sortie(tbloc+3,bloclev,2)
REAL(KIND=JPRB),INTENT(INOUT)    :: ZSDIVPL (1:YDGEOMETRY%YRDIM%NSMAX+1,2,YDGEOMETRY%YRDIMV%NFLEVG,500)
REAL(KIND=JPRB),INTENT(INOUT)    :: ZSPDIVPL(1:YDGEOMETRY%YRDIM%NSMAX+1,2,YDGEOMETRY%YRDIMV%NFLEVG,500)

#else
REAL(KIND=JPRB) :: ZSDIVPL  (YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NSMAX+1,2)
REAL(KIND=JPRB) :: ZSPDIVPL (YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NSMAX+1,2)
#endif

INTEGER(KIND=JPIM) :: II, IS0, IS02, ISE, JN,compteur,jmloc,ji,decalage1,klx
INTEGER(KIND=JPIM) :: IM, ISTA, IEND,jl,jlb,decalage,reste,compteurc,compteurb

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxture.h"
#include "mxturs.h"

!     ------------------------------------------------------------------

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX,NFLEVG=>YDDIMV%NFLEVG,SIHEG=>YDDYN%SIHEG,SIHEG2=>YDDYN%SIHEG2,NSPSTAF=>YDMP%NSPSTAF)

!             Inversion of two tridiagonal systems (Helmholtz equation)
!                --> (SIMI*DIVprim(t+dt)).
!$acc data present(psdivp,pspdivp,nsmax,nflevg) create(pas,pbs,pcs,entree,sortie)
!$acc data present(YDLAP,YDLAP%MYMS,NSPSTAF,SIHEG,siheg2,param_mxture) create(zsdivpl,zspdivpl)

#if defined(_OPENACC)

!$acc parallel private(ii,im,ista,klx,pas,pbs,pcs,entree,sortie,jlb,decalage,reste,tbloc,compteur,decalage1) default(none)
!$acc cache(pas(1:tbloc+3,1:bloclev),pbs(1:tbloc+3,1:bloclev),pcs(1:tbloc+3,1:bloclev),entree(1:tbloc+3,1:bloclev,1:2),sortie(1:tbloc+3,1:bloclev,1:2))
!$acc loop gang collapse(2) 
do jmloc=kmlocsta,kmlocend
  do compteurb=1,(nflevg-1)/bloclev+1

      IM=YDLAP%MYMS(jmloc)
      ISTA=NSPSTAF(IM)
      klx=nsmax+1-im
      II=min(IM,1)+1
      !!zsdivpl(:,:,(compteurb-1)*bloclev+1:min(nflevg,compteurb*bloclev),jmloc)=0.0_JPRB
      !!zspdivpl(:,:,(compteurb-1)*bloclev+1:min(nflevg,compteurb*bloclev),jmloc)=0.0_JPRB
     
      !$acc loop vector private(ise,compteurc,compteur,ji)      
      DO JN=IM,NSMAX
        ISE=ISTA+2*(JN-IM)
        do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
          compteur=compteurc+(compteurb-1)*bloclev
          do ji=1,2
            ZSDIVPL(JN-im+1,ji,compteur,jmloc)=PSDIVP(ISE+ji-1,compteur)
          enddo
        enddo
      ENDDO 

      compteur=(compteurb-1)*bloclev
      decalage1=compteur*(nsmax+1-im)
     
      IF (IM > 0) THEN

        !               Inversion of a symmetric matrix.

        CALL MXTURS(NSMAX+1-IM,min(bloclev,nflevg-compteur),bloclev,II,2,NSMAX,&
          & param_mxture(decalage1+1,jmloc,1),param_mxture(decalage1+1,jmloc,2),&
          & param_mxture(decalage1+1,jmloc,3),&
          & ZSDIVPL(1,1,1+compteur,jmloc),ZSPDIVPL(1,1,1+compteur,jmloc),&
          & tbloc,pas,pbs,pcs,entree,sortie)    

     ELSE
      do ji=1,2
        !               Inversion of a non-symmetric matrix.
        if (ji==1) then

          CALL MXTURE(NSMAX+1-IM,min(bloclev,nflevg-compteur),bloclev,II,2,NSMAX,-2,.TRUE.,&
            & param_mxture(decalage1+1,jmloc,1),param_mxture(decalage1+1,jmloc,2),&
            & param_mxture(decalage1+1,jmloc,3),&
            & ZSDIVPL(1,1,1+compteur,jmloc),ZSPDIVPL(1,1,1+compteur,jmloc),&
            & tbloc,pas,pbs,pcs,entree,sortie)

          CALL MXTURE(NSMAX+1-IM,min(bloclev,nflevg-compteur),bloclev,II,2,NSMAX,3,.FALSE.,&
            & param_mxture(decalage1+1,jmloc,1),param_mxture(decalage1+1,jmloc,4),&
            & param_mxture(decalage1+1,jmloc,5),&
            & ZSDIVPL(1,1,1+compteur,jmloc),ZSPDIVPL(1,1,1+compteur,jmloc),&
            & tbloc,pas,pbs,pcs,entree,sortie) 

        else
          !$acc loop vector private(compteur,compteurc)
          do jn=1,nsmax+1-im
            do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
              compteur=compteurc+(compteurb-1)*bloclev
              zspdivpl(jn,ji,compteur,jmloc)=0.0_JPRB !!zsdivpl(jn,ji,compteur,jmloc)
            enddo
          enddo
        endif
       enddo !!ji im 0
      ENDIF

      !$acc loop vector private(ise,compteurc,compteur,ji) 
      DO JN=IM,NSMAX
        ISE=ISTA+2*(JN-IM)
        do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)          
          compteur=compteurc+(compteurb-1)*bloclev
          do ji=1,2
            PSPDIVP(ISE+ji-1,compteur)=ZSPDIVPL(JN-im+1,ji,compteur,jmloc)
          enddo
        enddo
      enddo
  enddo  !!compteur
ENDDO    !!jmloc
!$acc end parallel

!$acc end data
!$acc end data

#else

!$omp parallel do private(jmloc,im,ista,iend,is0,is02,ii,jn,ise,zsdivpl,zspdivpl)
do jmloc=kmlocsta,kmlocend

IM=YDLAP%MYMS(jmloc)
ISTA=NSPSTAF(IM)
IEND=ISTA+2*(NSMAX+1-IM)-1

IS0=YDLAP%NSE0L(jmloc)
IS02=0
II=MIN(IM,1)+1
ZSDIVPL(:,:,:)=0.0_JPRB
ZSPDIVPL(:,:,:)=0.0_JPRB

DO JN=IM,NSMAX
    ISE=ISTA+2*(JN-IM)
    ZSDIVPL(:,JN-im+1,1)=PSDIVP(ISE,:)
    ZSDIVPL(:,JN-im+1,2)=PSDIVP(ISE+1,:)
ENDDO

IF (IM > 0) THEN

  !               Inversion of a symmetric matrix.
  CALL MXTURS(NSMAX+1-IM,NFLEVG,NFLEVG,II,nsmax,&
   & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
   & ZSDIVPL,ZSPDIVPL)  
ELSE

  !               Inversion of a non-symmetric matrix.
  CALL MXTURE(NSMAX+1-IM,NFLEVG,NFLEVG,II,nsmax,-2,.TRUE.,&
   & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
   & ZSDIVPL,ZSPDIVPL)  
  CALL MXTURE(NSMAX+1-IM,NFLEVG,NFLEVG,II,nsmax,3,.FALSE.,&
   & SIHEG(1,IS0+1,1),SIHEG2(1,IS02+1,2),&
   & SIHEG2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)
ENDIF

DO JN=IM,NSMAX
!!  do compteur=1,nflevg
    ISE=ISTA+2*(JN-IM)
    PSPDIVP(ISE,:)=ZSPDIVPL(:,JN-im+1,1)
    PSPDIVP(ISE+1,:)=ZSPDIVPL(:,JN-im+1,2)
!!  enddo
ENDDO

enddo !!jmloc
!$omp end parallel do

#endif

!!$acc update host(pspdivp)
!write (0,*) __FILE__,';',__LINE__
!write (0,*), pspdivp(:,:)
!call flush(0)

END ASSOCIATE
END ASSOCIATE

END SUBROUTINE SPCSIDG_PART1

