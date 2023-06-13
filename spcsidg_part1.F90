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
integer(kind=jpim)               :: taillec
integer(kind=jpim), parameter    :: tbloc=254!!97
REAL(KIND=JPRB),   intent(inout) :: ZSDIVPL (1:YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2,499)
REAL(KIND=JPRB),   intent(inout) :: ZSPDIVPL(1:YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2,499)
real(kind=jprb),   intent(in)    :: param_mxture(:,:,:)
real(kind=jprb),   intent(inout) :: pas(tbloc+3)
real(kind=jprb),   intent(inout) :: pbs(tbloc+3)
real(kind=jprb),   intent(inout) :: pcs(tbloc+3)
real(kind=jprb),   intent(inout) :: entree(tbloc+3)
real(kind=jprb),   intent(inout) :: sortie(tbloc+3)

#else
REAL(KIND=JPRB) :: ZSDIVPL  (YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NSMAX+1,2)
REAL(KIND=JPRB) :: ZSPDIVPL (YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NSMAX+1,2)
#endif

INTEGER(KIND=JPIM) :: II, IS0, IS02, ISE, JN,compteur,jmloc,ji,decalage1,klx
INTEGER(KIND=JPIM) :: IM, ISTA, IEND,jl,jlb,decalage,reste

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

!$acc parallel private(im,ista,decalage1,klx,pas,pbs,pcs,entree,sortie,jlb,jl,decalage,reste,tbloc) default(none)
!$acc cache(pas(1:tbloc+3),pbs(1:tbloc+3),pcs(1:tbloc+3),entree(1:tbloc+3),sortie(1:tbloc+3))
!$acc loop gang collapse(3) 
do jmloc=kmlocsta,kmlocend
  do compteur=1,nflevg
    do ji=1,2

      IM=YDLAP%MYMS(jmloc)
      ISTA=NSPSTAF(IM)
      decalage1=(compteur-1)*(nsmax+1-im)
      klx=nsmax+1-im
      !$acc loop vector private(ise)      
      DO JN=IM,NSMAX
        ISE=ISTA+2*(JN-IM)
        ZSDIVPL(JN-im+1,compteur,ji,jmloc)=PSDIVP(ISE+ji-1,compteur)
       ENDDO 
      
      IF (IM > 0) THEN

        !               Inversion of a symmetric matrix.

        IF (KLX >= 3) THEN
          !$acc loop seq
          do jlb=1,(klx-3)/tbloc+1  !!klx-3+1 elements, de 3 a klx
            decalage=(jlb-1)*tbloc
            !$acc loop vector
            do jl=decalage+1,min(decalage+tbloc+2,klx)
              pas(jl-decalage)=param_mxture(decalage1+jl,jmloc,1)
              pbs(jl-decalage)=param_mxture(decalage1+jl,jmloc,2)
              pcs(jl-decalage)=param_mxture(decalage1+jl,jmloc,3)
              entree(jl-decalage)=zsdivpl(jl,compteur,ji,jmloc)
            enddo
            if (jlb==1) then
              sortie(1)=entree(1)/pas(1)
              sortie(2)=(entree(2)-pbs(1)*sortie(1))/pas(2)
            else
              sortie(1)=zspdivpl(decalage+1,compteur,ji,jmloc)
              sortie(2)=zspdivpl(decalage+2,compteur,ji,jmloc)
            endif
            !$acc loop seq
            do jl=3,min(decalage+tbloc+2,klx)-decalage
              sortie(JL)=(entree(JL)-pcs(JL-2)*sortie(JL-2)&
               & -pbs(JL-1)*sortie(JL-1))/pas(JL)  
            enddo
            !$acc loop vector
            do jl=decalage+1,min(decalage+tbloc+2,klx)
              zspdivpl(jl,compteur,ji,jmloc)=sortie(jl-decalage)
            enddo
          enddo
        ELSE
          zspdivpl(1,compteur,ji,jmloc)=zsdivpl(1,compteur,ji,jmloc)/param_mxture(decalage1+1,jmloc,1)
          IF (KLX >= 2) THEN
              zspdivpl(2,compteur,ji,jmloc)=(zsdivpl(2,compteur,ji,jmloc)&
               &-param_mxture(decalage1+1,jmloc,2)*zspdivpl(1,compteur,ji,jmloc))/param_mxture(decalage1+2,jmloc,1)
          ENDIF
        ENDIF

       !$acc loop vector
       do jl=1,klx
         zsdivpl(JL,compteur,ji,jmloc)=zspdivpl(JL,compteur,ji,jmloc)
       ENDDO

       IF (KLX >= 3) THEN
         reste=mod(klx-3,tbloc)+1
         !$acc loop seq
         do jlb=(klx-3)/tbloc+1,1,-1  !!klx-3+1 elements, de klx-2 a 1
           decalage=(jlb-2)*tbloc+reste
           !$acc loop vector
           do jl=decalage+tbloc+2,max(decalage+1,1),-1
             pas(jl-decalage)=param_mxture(decalage1+jl,jmloc,1)
             pbs(jl-decalage)=param_mxture(decalage1+jl,jmloc,2)
             pcs(jl-decalage)=param_mxture(decalage1+jl,jmloc,3)
             entree(jl-decalage)=zsdivpl(jl,compteur,ji,jmloc)
           enddo
           if (jlb==(klx-3)/tbloc+1) then
              !!pxs(klx-decalage)=pys(klx-decalage)
              !!pxs(klx-1-decalage)=pys(klx-1-decalage)-pbs(klx-1-decalage)/pas(klx-1-decalage)*pxs(klx-decalage)
              sortie(klx-decalage)=entree(klx-decalage)
              sortie(klx-1-decalage)=entree(klx-1-decalage)-pbs(klx-1-decalage)/pas(klx-1-decalage)*sortie(klx-decalage)
           else
              sortie(tbloc+1)=zspdivpl(decalage+tbloc+1,compteur,ji,jmloc)
              sortie(tbloc+2)=zspdivpl(decalage+tbloc+2,compteur,ji,jmloc)
           endif
           !$acc loop seq
           do jl=tbloc,max(decalage+1,1)-decalage,-1
             sortie(JL)=entree(JL)-pbs(jl)/pas(jl)*sortie(JL+1)-pcs(jl)/pas(jl)*sortie(JL+2)     
           enddo
           !$acc loop vector
           do jl=decalage+tbloc+2,max(decalage+1,1),-1
             zspdivpl(jl,compteur,ji,jmloc)=sortie(jl-decalage)
           enddo
         enddo
       ELSE
         zspdivpl(KLX,compteur,ji,jmloc)=zsdivpl(KLX,compteur,ji,jmloc)
         IF (KLX >= 2) THEN
           zspdivpl(KLX-1,compteur,ji,jmloc)=zsdivpl(KLX-1,compteur,ji,jmloc)&
            &-param_mxture(decalage1+klx-1,jmloc,2)/param_mxture(decalage1+klx-1,jmloc,1)*zspdivpl(KLX,compteur,ji,jmloc)
         ENDIF
       ENDIF


     ELSE

        !               Inversion of a non-symmetric matrix.
       if (ji==1) then

        IF (KLX >= 3) THEN
          !$acc loop seq
          do jlb=1,(klx-3)/tbloc+1  !!klx-3+1 elements, de 3 a klx
          decalage=(jlb-1)*tbloc
          !$acc loop vector
          do jl=decalage+1,min(decalage+tbloc+2,klx)
            pas(jl-decalage)=param_mxture(decalage1+jl,jmloc,1)
            pbs(jl-decalage)=param_mxture(decalage1+jl,jmloc,2)
            pcs(jl-decalage)=param_mxture(decalage1+jl,jmloc,3)
            entree(jl-decalage)=zsdivpl(jl,compteur,ji,jmloc)
          enddo
          if (jlb==1) then
            sortie(1)=entree(1)/pas(1)
            sortie(2)=(entree(2)-pbs(1)*sortie(1))/pas(2)
          else
            sortie(1)=zspdivpl(decalage+1,compteur,ji,jmloc)
            sortie(2)=zspdivpl(decalage+2,compteur,ji,jmloc)
          endif
          !$acc loop seq
          do jl=3,min(decalage+tbloc+2,klx)-decalage
            sortie(JL)=(entree(JL)-pcs(JL-2)*sortie(JL-2)&
             & -pbs(JL-1)*sortie(JL-1))/pas(JL)  
          enddo
          !$acc loop vector
          do jl=decalage+1,min(decalage+tbloc+2,klx)
            zspdivpl(jl,compteur,ji,jmloc)=sortie(jl-decalage)
          enddo
        enddo
      ELSE
        zspdivpl(1,compteur,ji,jmloc)=zsdivpl(1,compteur,ji,jmloc)/param_mxture(decalage1+1,jmloc,1)
        IF (KLX >= 2) THEN
            zspdivpl(2,compteur,ji,jmloc)=(zsdivpl(2,compteur,ji,jmloc)&
             &-param_mxture(decalage1+1,jmloc,2)*zspdivpl(1,compteur,ji,jmloc))/param_mxture(decalage1+2,jmloc,1)
        ENDIF
      ENDIF

     !$acc loop vector
     do jl=1,klx
       zsdivpl(JL,compteur,ji,jmloc)=zspdivpl(JL,compteur,ji,jmloc)
     ENDDO

     IF (KLX >= 3) THEN
       reste=mod(klx-3,tbloc)+1
       !$acc loop seq
       do jlb=(klx-3)/tbloc+1,1,-1  !!klx-3+1 elements, de klx-2 a 1
         decalage=(jlb-2)*tbloc+reste
         !$acc loop vector
         do jl=decalage+tbloc+2,max(decalage+1,1),-1
           pbs(jl-decalage)=param_mxture(decalage1+jl,jmloc,4)
           pcs(jl-decalage)=param_mxture(decalage1+jl,jmloc,5)
           entree(jl-decalage)=zsdivpl(jl,compteur,ji,jmloc)
         enddo
         if (jlb==(klx-3)/tbloc+1) then
           sortie(klx-decalage)=entree(klx-decalage)
           sortie(klx-1-decalage)=entree(klx-1-decalage)-pbs(klx-1-decalage)*sortie(klx-decalage)
         else
           sortie(tbloc+1)=zspdivpl(decalage+tbloc+1,compteur,ji,jmloc)
           sortie(tbloc+2)=zspdivpl(decalage+tbloc+2,compteur,ji,jmloc)
         endif
         !$acc loop seq
         do jl=tbloc,max(decalage+1,1)-decalage,-1
           sortie(JL)=entree(JL)-pbs(jl)*sortie(JL+1)-pcs(jl)*sortie(JL+2)
         enddo
         !$acc loop vector
         do jl=decalage+tbloc+2,max(decalage+1,1),-1
           zspdivpl(jl,compteur,ji,jmloc)=sortie(jl-decalage)
         enddo
       enddo
     ELSE
       zspdivpl(KLX,compteur,ji,jmloc)=zsdivpl(KLX,compteur,ji,jmloc)
       IF (KLX >= 2) THEN
         zspdivpl(KLX-1,compteur,ji,jmloc)=zsdivpl(KLX-1,compteur,ji,jmloc)&
           &-param_mxture(decalage1+klx-1,jmloc,4)*zspdivpl(KLX,compteur,ji,jmloc)
       ENDIF
     ENDIF


        else
          !$acc loop vector
          do jn=1,nsmax+1-im
            zspdivpl(jn,compteur,ji,jmloc)=zsdivpl(jn,compteur,ji,jmloc)
          enddo
        endif
      ENDIF

      !$acc loop vector private(ise) 
      DO JN=IM,NSMAX
        ISE=ISTA+2*(JN-IM)
        PSPDIVP(ISE+ji-1,compteur)=ZSPDIVPL(JN-im+1,compteur,ji,jmloc)
      enddo
    enddo!!ji
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

