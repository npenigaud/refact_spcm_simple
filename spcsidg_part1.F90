#if defined(_OPENACC)
SUBROUTINE SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, PSDIVP,PSPDIVP,taillec,&
  &zsdivpl,zspdivpl ,pas,pbs,pcs,entree,sortie,param_mxture, kmlocsta, kmlocend)

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
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCsta
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOCend


REAL(KIND=JPRB),   INTENT(IN)    :: PSDIVP (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSPDIVP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
#if defined(_OPENACC)
INTEGER(KIND=JPIM)               :: taillec
INTEGER(KIND=JPIM), PARAMETER    :: tbloc=62
INTEGER(KIND=JPIM), PARAMETER     :: bloclev=8
REAL(KIND=JPRB),   INTENT(INOUT) :: ZSDIVPL (1:YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2,499)
REAL(KIND=JPRB),   INTENT(INOUT) :: ZSPDIVPL(1:YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2,499)
REAL(KIND=JPRB),   INTENT(IN)    :: param_mxture(:,:,:)
REAL(KIND=JPRB),   INTENT(INOUT) :: pas(tbloc+3,bloclev)
REAL(KIND=JPRB),   INTENT(INOUT) :: pbs(tbloc+3,bloclev)
REAL(KIND=JPRB),   INTENT(INOUT) :: pcs(tbloc+3,bloclev)
REAL(KIND=JPRB),   INTENT(INOUT) :: entree(tbloc+3,bloclev,2)
REAL(KIND=JPRB),   INTENT(INOUT) :: sortie(tbloc+3,bloclev,2)

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
!$acc data present(psdivp,pspdivp,zsdivpl,zspdivpl,nsmax,nflevg,pas,pbs,pcs,entree,sortie)
!$acc data present(YDLAP,YDLAP%MYMS,NSPSTAF,SIHEG,siheg2,param_mxture)

#if defined(_OPENACC)

!$acc parallel private(im,ista,klx,pas,pbs,pcs,entree,sortie,jlb,decalage,reste,tbloc) default(none)
!$acc cache(pas(1:tbloc+3,bloclev),pbs(1:tbloc+3,bloclev),pcs(1:tbloc+3,bloclev),entree(1:tbloc+3,1:bloclev,1:2),sortie(1:tbloc+3,1:bloclev,1:2))
!$acc loop gang collapse(2) 
do jmloc=kmlocsta,kmlocend
  do compteurb=1,(nflevg-1)/bloclev+1

      IM=YDLAP%MYMS(jmloc)
      ISTA=NSPSTAF(IM)
      klx=nsmax+1-im
      !$acc loop vector private(ise,compteurc,compteur,ji)      
      DO JN=IM,NSMAX
        ISE=ISTA+2*(JN-IM)
        do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
          compteur=compteurc+(compteurb-1)*bloclev
          do ji=1,2
            ZSDIVPL(JN-im+1,compteur,ji,jmloc)=PSDIVP(ISE+ji-1,compteur)
          enddo
        enddo
      ENDDO 
      
      IF (IM > 0) THEN

        !               Inversion of a symmetric matrix.

        IF (KLX >= 3) THEN
          !$acc loop seq
          do jlb=1,(klx-3)/tbloc+1  !!klx-3+1 elements, de 3 a klx
            decalage=(jlb-1)*tbloc
            !$acc loop vector private(compteurc,compteur,decalage1)
            do jl=decalage+1,min(decalage+tbloc+2,klx)
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                compteur=compteurc+(compteurb-1)*bloclev
                decalage1=(compteur-1)*(nsmax+1-im)
                pas(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,1)
                pbs(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,2)
                pcs(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,3)
                do ji=1,2
                  entree(jl-decalage,compteurc,ji)=zsdivpl(jl,compteur,ji,jmloc)
                enddo
              enddo
            enddo
            if (jlb==1) then
              !$acc loop vector private(compteur) collapse(2)
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                do ji=1,2
                  sortie(1,compteurc,ji)=entree(1,compteurc,ji)/pas(1,compteurc)
                  sortie(2,compteurc,ji)=(entree(2,compteurc,ji)-pbs(1,compteurc)*sortie(1,compteurc,ji))/pas(2,compteurc)
                enddo
              enddo
            else
              !$acc loop vector private(compteur) collapse(2)
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                do ji=1,2
                  compteur=compteurc+(compteurb-1)*bloclev
                  sortie(1,compteurc,ji)=zspdivpl(decalage+1,compteur,ji,jmloc)
                  sortie(2,compteurc,ji)=zspdivpl(decalage+2,compteur,ji,jmloc)
                enddo
              enddo
            endif
            !$acc loop vector private(jl) collapse(2)
            do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
              do ji=1,2
                do jl=3,min(decalage+tbloc+2,klx)-decalage
                  sortie(JL,compteurc,ji)=(entree(JL,compteurc,ji)-pcs(JL-2,compteurc)*sortie(JL-2,compteurc,ji)&
                   & -pbs(JL-1,compteurc)*sortie(JL-1,compteurc,ji))/pas(JL,compteurc)  
                enddo
              enddo
            enddo
            !$acc loop vector private(compteurc,ji)
            do jl=decalage+1,min(decalage+tbloc+2,klx)
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                do ji=1,2
                  zspdivpl(jl,compteurc+(compteurb-1)*bloclev,ji,jmloc)=sortie(jl-decalage,compteurc,ji)
                enddo
              enddo
            enddo
          enddo
        ELSE
          !$acc loop vector private(compteur,decalage1) collapse(2)
          do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
            do ji=1,2
              compteur=compteurc+(compteurb-1)*bloclev
              decalage1=(compteur-1)*(nsmax+1-im)
              zspdivpl(1,compteur,ji,jmloc)=zsdivpl(1,compteur,ji,jmloc)/param_mxture(decalage1+1,jmloc,1)
              IF (KLX >= 2) THEN
                zspdivpl(2,compteur,ji,jmloc)=(zsdivpl(2,compteur,ji,jmloc)&
                 &-param_mxture(decalage1+1,jmloc,2)*zspdivpl(1,compteur,ji,jmloc))/param_mxture(decalage1+2,jmloc,1)
              ENDIF
            enddo
          enddo
        ENDIF

       !$acc loop vector private(compteurc,compteur,ji)
       do jl=1,klx
         do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
           compteur=compteurc+(compteurb-1)*bloclev
           do ji=1,2
             zsdivpl(JL,compteur,ji,jmloc)=zspdivpl(JL,compteur,ji,jmloc)
           enddo
         enddo
       ENDDO

       IF (KLX >= 3) THEN
         reste=mod(klx-3,tbloc)+1
         !$acc loop seq
         do jlb=(klx-3)/tbloc+1,1,-1  !!klx-3+1 elements, de klx-2 a 1
           decalage=(jlb-2)*tbloc+reste
           !$acc loop vector private(compteur,compteurc,decalage1,ji)
           do jl=decalage+tbloc+2,max(decalage+1,1),-1
             do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
               compteur=compteurc+(compteurb-1)*bloclev
               decalage1=(compteur-1)*(nsmax+1-im)
               pas(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,1)
               pbs(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,2)
               pcs(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,3)
               do ji=1,2
                 entree(jl-decalage,compteurc,ji)=zsdivpl(jl,compteur,ji,jmloc)
               enddo
             enddo
           enddo
           if (jlb==(klx-3)/tbloc+1) then
              !$acc loop vector collapse(2)
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                do ji=1,2
                  sortie(klx-decalage,compteurc,ji)=entree(klx-decalage,compteurc,ji)
                  sortie(klx-1-decalage,compteurc,ji)=entree(klx-1-decalage,compteurc,ji)&
                   &-pbs(klx-1-decalage,compteurc)/pas(klx-1-decalage,compteurc)*sortie(klx-decalage,compteurc,ji)
                enddo
              enddo
           else
              !$acc loop vector private(compteur) collapse(2)
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                do ji=1,2
                  compteur=compteurc+(compteurb-1)*bloclev
                  sortie(tbloc+1,compteurc,ji)=zspdivpl(decalage+tbloc+1,compteur,ji,jmloc)
                  sortie(tbloc+2,compteurc,ji)=zspdivpl(decalage+tbloc+2,compteur,ji,jmloc)
                enddo
              enddo
           endif
           !$acc loop vector private(jl) collapse(2)
           do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
             do ji=1,2
             do jl=tbloc,max(decalage+1,1)-decalage,-1
                 sortie(JL,compteurc,ji)=entree(JL,compteurc,ji)-pbs(jl,compteurc)/pas(jl,compteurc)*sortie(JL+1,compteurc,ji)&
                   &-pcs(jl,compteurc)/pas(jl,compteurc)*sortie(JL+2,compteurc,ji)     
               enddo
             enddo
           enddo
           !$acc loop vector private(compteurc,ji)
           do jl=decalage+tbloc+2,max(decalage+1,1),-1
             do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
               do ji=1,2
                 zspdivpl(jl,compteurc+(compteurb-1)*bloclev,ji,jmloc)=sortie(jl-decalage,compteurc,ji)
               enddo
             enddo
           enddo
         enddo
       ELSE
         !$acc loop vector private(compteur,decalage1) collapse(2)
         do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
           do ji=1,2
             compteur=compteurc+(compteurb-1)*bloclev
             decalage1=(compteur-1)*(nsmax+1-im)
             zspdivpl(KLX,compteur,ji,jmloc)=zsdivpl(KLX,compteur,ji,jmloc)
             IF (KLX >= 2) THEN
               zspdivpl(KLX-1,compteur,ji,jmloc)=zsdivpl(KLX-1,compteur,ji,jmloc)&
                &-param_mxture(decalage1+klx-1,jmloc,2)/param_mxture(decalage1+klx-1,jmloc,1)*zspdivpl(KLX,compteur,ji,jmloc)
             ENDIF
           enddo
         enddo
       ENDIF


     ELSE
      do ji=1,2
        !               Inversion of a non-symmetric matrix.
       if (ji==1) then

        IF (KLX >= 3) THEN
          !$acc loop seq
          do jlb=1,(klx-3)/tbloc+1  !!klx-3+1 elements, de 3 a klx
            decalage=(jlb-1)*tbloc
            !$acc loop vector private(compteurc,compteur,decalage1)
            do jl=decalage+1,min(decalage+tbloc+2,klx)
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                compteur=compteurc+(compteurb-1)*bloclev
                decalage1=(compteur-1)*(nsmax+1-im)
                pas(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,1)
                pbs(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,2)
                pcs(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,3)
                entree(jl-decalage,compteurc,ji)=zsdivpl(jl,compteur,ji,jmloc)
              enddo
            enddo
            if (jlb==1) then
              !$acc loop vector
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                sortie(1,compteurc,ji)=entree(1,compteurc,ji)/pas(1,compteurc)
                sortie(2,compteurc,ji)=(entree(2,compteurc,ji)-pbs(1,compteurc)*sortie(1,compteurc,ji))/pas(2,compteurc)
              enddo
            else
              !$acc loop vector private(compteur)
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                compteur=compteurc+(compteurb-1)*bloclev
                sortie(1,compteurc,ji)=zspdivpl(decalage+1,compteur,ji,jmloc)
                sortie(2,compteurc,ji)=zspdivpl(decalage+2,compteur,ji,jmloc)
              enddo
            endif
            !$acc loop vector private(jl)
            do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
              do jl=3,min(decalage+tbloc+2,klx)-decalage
                sortie(JL,compteurc,ji)=(entree(JL,compteurc,ji)-pcs(JL-2,compteurc)*sortie(JL-2,compteurc,ji)&
                 & -pbs(JL-1,compteurc)*sortie(JL-1,compteurc,ji))/pas(JL,compteurc)  
              enddo
            enddo
            !$acc loop vector private(compteurc)
            do jl=decalage+1,min(decalage+tbloc+2,klx)
              do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
                zspdivpl(jl,compteurc+(compteurb-1)*bloclev,ji,jmloc)=sortie(jl-decalage,compteurc,ji)
              enddo
            enddo
          enddo
        ELSE
          !$acc loop vector private(compteur,decalage1)
          do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
            compteur=compteurc+(compteurb-1)*bloclev
            decalage1=(compteur-1)*(nsmax+1-im)
            zspdivpl(1,compteur,ji,jmloc)=zsdivpl(1,compteur,ji,jmloc)/param_mxture(decalage1+1,jmloc,1)
            IF (KLX >= 2) THEN
              zspdivpl(2,compteur,ji,jmloc)=(zsdivpl(2,compteur,ji,jmloc)&
               &-param_mxture(decalage1+1,jmloc,2)*zspdivpl(1,compteur,ji,jmloc))/param_mxture(decalage1+2,jmloc,1)
            ENDIF
          enddo
        ENDIF

       !$acc loop vector private(compteurc,compteur)
       do jl=1,klx
         do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
           compteur=compteurc+(compteurb-1)*bloclev
           zsdivpl(JL,compteur,ji,jmloc)=zspdivpl(JL,compteur,ji,jmloc)
         enddo
       ENDDO

     IF (KLX >= 3) THEN
       reste=mod(klx-3,tbloc)+1
       !$acc loop seq
       do jlb=(klx-3)/tbloc+1,1,-1  !!klx-3+1 elements, de klx-2 a 1
         decalage=(jlb-2)*tbloc+reste
         !$acc loop vector private(compteurc,compteur,decalage1)
         do jl=decalage+tbloc+2,max(decalage+1,1),-1
           do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
             compteur=compteurc+(compteurb-1)*bloclev
             decalage1=(compteur-1)*(nsmax+1-im)
             pbs(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,4)
             pcs(jl-decalage,compteurc)=param_mxture(decalage1+jl,jmloc,5)
             entree(jl-decalage,compteurc,ji)=zsdivpl(jl,compteur,ji,jmloc)
           enddo
         enddo
         if (jlb==(klx-3)/tbloc+1) then
           !$acc loop vector
           do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
             sortie(klx-decalage,compteurc,ji)=entree(klx-decalage,compteurc,ji)
             sortie(klx-1-decalage,compteurc,ji)=entree(klx-1-decalage,compteurc,ji)&
               &-pbs(klx-1-decalage,compteurc)*sortie(klx-decalage,compteurc,ji)
           enddo
         else
           !$acc loop vector private(compteur)
           do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
             compteur=compteurc+(compteurb-1)*bloclev
             sortie(tbloc+1,compteurc,ji)=zspdivpl(decalage+tbloc+1,compteur,ji,jmloc)
             sortie(tbloc+2,compteurc,ji)=zspdivpl(decalage+tbloc+2,compteur,ji,jmloc)
           enddo
         endif
         !$acc loop vector private(jl)
         do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
           !$acc loop seq
           do jl=tbloc,max(decalage+1,1)-decalage,-1
             sortie(JL,compteurc,ji)=entree(JL,compteurc,ji)-pbs(jl,compteurc)*sortie(JL+1,compteurc,ji)&
               &-pcs(jl,compteurc)*sortie(JL+2,compteurc,ji)
           enddo
         enddo
         !$acc loop vector private(compteurc)
         do jl=decalage+tbloc+2,max(decalage+1,1),-1
           do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
             zspdivpl(jl,compteurc+(compteurb-1)*bloclev,ji,jmloc)=sortie(jl-decalage,compteurc,ji)
           enddo
         enddo
       enddo
     ELSE
       !$acc loop vector private(compteur,decalage1)
       do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
         compteur=compteurc+(compteurb-1)*bloclev
         decalage1=(compteur-1)*(nsmax+1-im)
         zspdivpl(KLX,compteur,ji,jmloc)=zsdivpl(KLX,compteur,ji,jmloc)
         IF (KLX >= 2) THEN
           zspdivpl(KLX-1,compteur,ji,jmloc)=zsdivpl(KLX-1,compteur,ji,jmloc)&
             &-param_mxture(decalage1+klx-1,jmloc,4)*zspdivpl(KLX,compteur,ji,jmloc)
         ENDIF
       enddo
     ENDIF


        else
          !$acc loop vector private(compteur,compteurc)
          do jn=1,nsmax+1-im
            do compteurc=1,min(bloclev,nflevg-(compteurb-1)*bloclev)
              compteur=compteurc+(compteurb-1)*bloclev
              zspdivpl(jn,compteur,ji,jmloc)=zsdivpl(jn,compteur,ji,jmloc)
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
            PSPDIVP(ISE+ji-1,compteur)=ZSPDIVPL(JN-im+1,compteur,ji,jmloc)
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

