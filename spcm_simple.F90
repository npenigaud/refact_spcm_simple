SUBROUTINE SPCM_SIMPLE (YDGEOMETRY,YDMODEL,PSPSP,PSPVOR,PSPDIV,PSPT,PSPSPD,PSPSVD)

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
use YOMMP0             , only : MYSETN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)         ,INTENT(INOUT) :: YDMODEL
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPSP (YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPVOR(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPDIV(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPT  (YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPSPD(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPSVD(YDGEOMETRY%YRDIMV%NFLEVL,YDGEOMETRY%YRDIM%NSPEC2)

#include "spcsi_cor.intfb.h"
#include "spcsi_str.intfb.h"
#include "trmtos.intfb.h"
#include "trstom.intfb.h"
#include "spcimpfinit.intfb.h"
#include "spcimpfpost.intfb.h"

#if defined(_OPENACC)
REAL(KIND=JPRB),ALLOCATABLE  :: param_mxture(:,:,:)
#else
REAL(KIND=JPRB)  :: SIMIT(YDGEOMETRY%yrdimv%nflevg,YDGEOMETRY%yrdimv%nflevg)
REAL(KIND=JPRB)  :: SIMOT(YDGEOMETRY%yrdimv%nflevg,YDGEOMETRY%yrdimv%nflevg)
#endif

REAL(KIND=JPRB), ALLOCATABLE :: ZSPVORG2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPDIVG2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPTG2  (:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPSPDG2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPSVDG2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPSPG2 (:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPVORG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPDIVG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPTG  (:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPSPDG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPSVDG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPSPG (:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPTNDSI_VORG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPTNDSI_DIVG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPTNDSI_TG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPAUX (:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPAUXG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: PSPVOR2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: PSPDIV2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: PSPT2  (:,:)
REAL(KIND=JPRB), ALLOCATABLE :: PSPSPD2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: PSPSVD2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: PSPSP2 (:)

LOGICAL :: LLONEM, LDTRANSPOSE

INTEGER(KIND=JPIM) :: IM, ISPEC2V
INTEGER (KIND=JPIM) :: JMLOC, ISTA, IEND

REAL(KIND=JPHOOK) ::  ZHOOK_HANDLE,zhook_handle2,zhook_handle3
INTEGER(KIND=JPIM)::  icnt1,icnt2,is0,is02,taille

IF (LHOOK) CALL DR_HOOK('SPCM_SIMPLE',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,  YDMP=>YDGEOMETRY%YRMP, YDDYN=>YDMODEL%YRML_DYN%YRDYN,   &
&  YDDYNA=>YDMODEL%YRML_DYN%YRDYNA, YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YDSPDDH=>YDMODEL%YRML_DIAG%YRSPDDH, &
&  YDLAP=>YDGEOMETRY%YRLAP)

ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, NSPEC2V=>YDMP%NSPEC2V, NSPEC2VF=>YDMP%NSPEC2VF, LIMPF=>YDDYN%LIMPF, &
& NFLSUR=>YDDIMV%NFLSUR, NSPEC2=>YDDIM%NSPEC2, NUMP=>YDDIM%NUMP, NSMAX=>YDDIM%NSMAX,nflevl=>yddimv%nflevl, &
& simi=>YDDYN%SIMI,simo=>YDDYN%SIMO,siheg=>yddyn%siheg,siheg2=>yddyn%siheg2, &
& nptrmf=>ydmp%nptrmf,LSIDG=>YDDYN%LSIDG)

#if defined(_OPENACC)
LDTRANSPOSE=.FALSE.
#else
LDTRANSPOSE=.TRUE.
#endif

#if defined(_OPENACC)
if (LSIDG) then
  !!transposition temporaire des parametres mxture
  !!calcul de la taille max d une plage de params
  taille=0
  do jmloc=nptrmf(mysetn),nptrmf(mysetn+1)-1
    im=ydlap%myms(jmloc)  
    taille=max(taille,nflevg*(nsmax+1-im))
  enddo

!!  print *,"taille",taille
!!  print *,"nptrmf(mysetn)",nptrmf(mysetn)
!!  print *,"nptrmf(mysetn+1)-1",nptrmf(mysetn+1)-1
!!  call flush(0)

  !!allocation, passage sur carte et initialisation
  allocate(param_mxture(taille,nptrmf(mysetn+1)-1,5))
  !$acc enter data create(param_mxture)
!!!!!!  !$acc data create(param_mxture)
  param_mxture(:,:,:)=0.0_JPRB

  !!charge et transp données sur carte
  do jmloc=nptrmf(mysetn),nptrmf(mysetn+1)-1
    im=ydlap%myms(jmloc)  
    is0=ydlap%nse0l(jmloc)
    is02=0
    if (.true.) then
      !$acc parallel loop private(icnt2)
      do icnt1=is0+1,is0+nsmax+1-im
        do icnt2=1,nflevg
          param_mxture(icnt1-is0-1+1+(nsmax+1-im)*(icnt2-1),jmloc,1)=siheg(icnt2,icnt1,1) 
          param_mxture(icnt1-is0-1+1+(nsmax+1-im)*(icnt2-1),jmloc,2)=siheg(icnt2,icnt1,2) 
          param_mxture(icnt1-is0-1+1+(nsmax+1-im)*(icnt2-1),jmloc,3)=siheg(icnt2,icnt1,3) 
        enddo
      enddo
      !$acc end parallel

      !$acc parallel loop private(icnt2)
      do icnt1=is02+1,is02+nsmax+1-im
        do icnt2=1,nflevg
          param_mxture(icnt1-is02-1+1+(nsmax+1-im)*(icnt2-1),jmloc,4)=siheg2(icnt2,icnt1,2) 
          param_mxture(icnt1-is02-1+1+(nsmax+1-im)*(icnt2-1),jmloc,5)=siheg2(icnt2,icnt1,3) 
        enddo
      enddo
      !$acc end parallel
    else
      !$acc parallel loop private(icnt2)
      do icnt1=is0+1,is0+nsmax+1-im
        do icnt2=1,nflevg
          param_mxture((icnt1-is0-1)*nflevg+1+(icnt2-1),jmloc,1)=siheg(icnt2,icnt1,1) 
          param_mxture((icnt1-is0-1)*nflevg+1+(icnt2-1),jmloc,2)=siheg(icnt2,icnt1,2) 
          param_mxture((icnt1-is0-1)*nflevg+1+(icnt2-1),jmloc,3)=siheg(icnt2,icnt1,3) 
        enddo
      enddo
      !$acc end parallel

      !$acc parallel loop private(icnt2)
      do icnt1=is02+1,is02+nsmax+1-im
        do icnt2=1,nflevg
          param_mxture((icnt1-is02-1)*nflevg+1+(icnt2-1),jmloc,4)=siheg2(icnt2,icnt1,2) 
          param_mxture((icnt1-is02-1)*nflevg+1+(icnt2-1),jmloc,5)=siheg2(icnt2,icnt1,3) 
        enddo
      enddo
      !$acc end parallel

    endif !!choix du sens de transposition

  enddo   !!jmloc

else      !! non lsidg lsidg 
  allocate(param_mxture(1,1,1))
  !$acc enter data create(param_mxture)
endif


#else
!!transposition de SIMIT SIMOT
!$omp parallel do private(icnt1,icnt2)
do icnt1=1,nflevg
  do icnt2=1,nflevg
    SIMIT(icnt1,icnt2)=simi(icnt2,icnt1)
  enddo
enddo
!$omp end parallel do

!$omp parallel do private(icnt1,icnt2)
do icnt1=1,nflevg
  do icnt2=1,nflevg
    SIMOT(icnt1,icnt2)=simo(icnt2,icnt1)
  enddo
enddo
!$omp end parallel do
#endif

LLONEM = .NOT. LIMPF

IF (LLONEM) THEN
  ISPEC2V=NSPEC2VF
ELSE
  ISPEC2V=NSPEC2V
ENDIF

ALLOCATE(ZSPVORG2(NFLEVG,ISPEC2V))
ALLOCATE(ZSPDIVG2(NFLEVG,ISPEC2V))
ALLOCATE(ZSPTG2  (NFLEVG,ISPEC2V))
ALLOCATE(ZSPSPDG2(NFLEVG,ISPEC2V))
ALLOCATE(ZSPSVDG2(NFLEVG,ISPEC2V))
ALLOCATE(ZSPSPG2 (ISPEC2V))

ALLOCATE(PSPVOR2(nspec2,nflevl))
ALLOCATE(PSPDIV2(nspec2,nflevl))
ALLOCATE(PSPT2  (nspec2,nflevl))
ALLOCATE(PSPSPD2(nspec2,nflevl))
ALLOCATE(PSPSVD2(nspec2,nflevl))
ALLOCATE(PSPSP2 (nspec2))

ALLOCATE(ZSPVORG(ISPEC2V,nflevg))
ALLOCATE(ZSPDIVG(ISPEC2V,nflevg))
ALLOCATE(ZSPTG  (ISPEC2V,nflevg))
ALLOCATE(ZSPSPDG(ISPEC2V,nflevg))
ALLOCATE(ZSPSVDG(ISPEC2V,nflevg))
ALLOCATE(ZSPSPG (ISPEC2V))

ALLOCATE(ZSPTNDSI_VORG(1,1))
ALLOCATE(ZSPTNDSI_DIVG(1,1))
ALLOCATE(ZSPTNDSI_TG(1,1))

IF (LIMPF) THEN
  ALLOCATE(ZSPAUX(NFLSUR,NSPEC2))
  ALLOCATE(ZSPAUXG(NFLEVG,ISPEC2V))

!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JMLOC,IM,ISTA,IEND)
  DO JMLOC=1,NUMP
    IM=YDLAP%MYMS(JMLOC)
    ISTA=YDLAP%NASM0(IM)
    IEND=ISTA+2*(NSMAX+1-IM)-1
    CALL SPCIMPFINIT(YDMODEL%YRCST,YDGEOMETRY,YDMODEL%YRML_GCONF%YRRIP,YDDYN,IM,ISTA,IEND,PSPVOR,ZSPAUX)
  ENDDO
!$OMP END PARALLEL DO

  !! LDTRANSPOSE to .false. as no data transposition implemented yet
  CALL TRMTOS(YDGEOMETRY,YDDYNA%LNHDYN,YDDYNA%LNHX,.FALSE.,&
    & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPT=PSPT,PSPSPD=PSPSPD,&
    & PSPSVD=PSPSVD,PSPSP=PSPSP,PSPAUX=ZSPAUX,& 
    & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,PSPSPDG=ZSPSPDG,&
    & PSPSVDG=ZSPSVDG,PSPSPG=ZSPSPG,PSPAUXG=ZSPAUXG,&
    & LDFULLM=LLONEM)

  CALL SPCSI_COR(YDGEOMETRY, YDMODEL%YRCST, YDLDDH, YDMODEL%YRML_GCONF%YRRIP, YDDYN, ISPEC2V, &
  & LLONEM, ZSPVORG, ZSPDIVG, ZSPTG, ZSPSPG, ZSPTNDSI_VORG, ZSPTNDSI_DIVG, ZSPTNDSI_TG,       &
  & PSPAUXG=ZSPAUXG)

  CALL TRSTOM(&
    & YDGEOMETRY,YDDYNA%LNHDYN,YDDYNA%LNHX,.FALSE.,&
    & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,PSPSPDG=ZSPSPDG,&
    & PSPSVDG=ZSPSVDG,PSPSPG=ZSPSPG,&
    & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPT=PSPT,PSPSPD=PSPSPD,&
    & PSPSVD=PSPSVD,PSPSP=PSPSP,&
    & LDFULLM=LLONEM,LDNEEDPS=.TRUE.)  

!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JMLOC,IM,ISTA,IEND)
    DO JMLOC=1,NUMP
      IM=YDLAP%MYMS(JMLOC)
      ISTA=YDLAP%NASM0(IM)
      IEND=ISTA+2*(NSMAX+1-IM)-1
      CALL SPCIMPFPOST(YDMODEL%YRCST,YDGEOMETRY,YDMODEL%YRML_GCONF%YRRIP,YDDYN,IM,ISTA,IEND,PSPDIV,PSPVOR)
    ENDDO
!$OMP END PARALLEL DO

ELSE

#if defined(_OPENACC)
    pspsp2(:)=pspsp(:)
    do icnt1=1,nspec2
      do icnt2=1,nflevl
        pspvor2(icnt1,icnt2)=pspvor(icnt2,icnt1)
        pspdiv2(icnt1,icnt2)=pspdiv(icnt2,icnt1)
        pspt2(icnt1,icnt2)=pspt(icnt2,icnt1)
        pspspd2(icnt1,icnt2)=pspspd(icnt2,icnt1)
        pspsvd2(icnt1,icnt2)=pspsvd(icnt2,icnt1)
      enddo
    enddo
#endif

    if (lhook) call dr_hook('SPCM_SIMPLE_transferts1a',0,zhook_handle2)
    !$acc data create(zspvorg,zspdivg,zsptg,zspspdg,zspsvdg,zspspg) 
    if (lhook) call dr_hook('SPCM_SIMPLE_transferts1a',1,zhook_handle2)
    if (lhook) call dr_hook('SPCM_SIMPLE_transferts1b',0,zhook_handle2)
    !$acc data copy(pspvor2,pspdiv2,pspt2,pspspd2,pspsvd2,pspsp2)
    if (lhook) call dr_hook('SPCM_SIMPLE_transferts1b',1,zhook_handle2)

  IF (LHOOK) CALL DR_HOOK('SPCM_SIMPLE_utile',0,ZHOOK_HANDLE3)

#if defined(_OPENACC)

  CALL TRMTOS(YDGEOMETRY,YDDYNA%LNHDYN,YDDYNA%LNHX,LDTRANSPOSE,&
    & PSPVOR=PSPVOR2,PSPDIV=PSPDIV2,PSPT=PSPT2,PSPSPD=PSPSPD2,&
    & PSPSVD=PSPSVD2,PSPSP=PSPSP2,&
    & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,PSPSPDG=ZSPSPDG,&
    & PSPSVDG=ZSPSVDG,PSPSPG=ZSPSPG,&
    & LDFULLM=LLONEM)

  CALL SPCSI_STR(YDGEOMETRY, YDMODEL%YRCST, YDLDDH, YDMODEL%YRML_GCONF%YRRIP, YDDYN, ISPEC2V, &
  & ZSPVORG, ZSPDIVG, ZSPTG, ZSPSPG, ZSPTNDSI_VORG, ZSPTNDSI_DIVG, ZSPTNDSI_TG,&
  & param_mxture)

  CALL TRSTOM(&
    & YDGEOMETRY,YDDYNA%LNHDYN,YDDYNA%LNHX,LDTRANSPOSE,&
    & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,PSPSPDG=ZSPSPDG,&
    & PSPSVDG=ZSPSVDG,PSPSPG=ZSPSPG,&
    & PSPVOR=PSPVOR2,PSPDIV=PSPDIV2,PSPT=PSPT2,PSPSPD=PSPSPD2,&
    & PSPSVD=PSPSVD2,PSPSP=PSPSP2,&
    & LDFULLM=LLONEM,LDNEEDPS=.TRUE.)  

#else 

if (LDTRANSPOSE) then

  CALL TRMTOS(YDGEOMETRY,YDDYNA%LNHDYN,YDDYNA%LNHX,LDTRANSPOSE,&
    & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPT=PSPT,PSPSPD=PSPSPD,&
    & PSPSVD=PSPSVD,PSPSP=PSPSP,&
    & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,PSPSPDG=ZSPSPDG,&
    & PSPSVDG=ZSPSVDG,PSPSPG=ZSPSPG,&
    & LDFULLM=LLONEM)

  CALL SPCSI_STR(YDGEOMETRY, YDMODEL%YRCST, YDLDDH, YDMODEL%YRML_GCONF%YRRIP, YDDYN, ISPEC2V, &
  & ZSPVORG, ZSPDIVG, ZSPTG, ZSPSPG, ZSPTNDSI_VORG, ZSPTNDSI_DIVG, ZSPTNDSI_TG,&
  & SIMIT,SIMOT)

  CALL TRSTOM(&
    & YDGEOMETRY,YDDYNA%LNHDYN,YDDYNA%LNHX,LDTRANSPOSE,&
    & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,PSPSPDG=ZSPSPDG,&
    & PSPSVDG=ZSPSVDG,PSPSPG=ZSPSPG,&
    & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPT=PSPT,PSPSPD=PSPSPD,&
    & PSPSVD=PSPSVD,PSPSP=PSPSP,&
    & LDFULLM=LLONEM,LDNEEDPS=.TRUE.)  

else

    pspsp2(:)=pspsp(:)
    do icnt1=1,nspec2
      do icnt2=1,nflevl
        pspvor2(icnt1,icnt2)=pspvor(icnt2,icnt1)
        pspdiv2(icnt1,icnt2)=pspdiv(icnt2,icnt1)
        pspt2(icnt1,icnt2)=pspt(icnt2,icnt1)
        pspspd2(icnt1,icnt2)=pspspd(icnt2,icnt1)
        pspsvd2(icnt1,icnt2)=pspsvd(icnt2,icnt1)
      enddo
    enddo

  CALL TRMTOS(YDGEOMETRY,YDDYNA%LNHDYN,YDDYNA%LNHX,LDTRANSPOSE,&
    & PSPVOR=PSPVOR2,PSPDIV=PSPDIV2,PSPT=PSPT2,PSPSPD=PSPSPD2,&
    & PSPSVD=PSPSVD2,PSPSP=PSPSP2,&
    & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,PSPSPDG=ZSPSPDG,&
    & PSPSVDG=ZSPSVDG,PSPSPG=ZSPSPG,&
    & LDFULLM=LLONEM)

  CALL SPCSI_STR(YDGEOMETRY, YDMODEL%YRCST, YDLDDH, YDMODEL%YRML_GCONF%YRRIP, YDDYN, ISPEC2V, &
  & ZSPVORG, ZSPDIVG, ZSPTG, ZSPSPG, ZSPTNDSI_VORG, ZSPTNDSI_DIVG, ZSPTNDSI_TG,&
  & SIMIT,SIMOT)

  CALL TRSTOM(&
    & YDGEOMETRY,YDDYNA%LNHDYN,YDDYNA%LNHX,LDTRANSPOSE,&
    & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,PSPSPDG=ZSPSPDG,&
    & PSPSVDG=ZSPSVDG,PSPSPG=ZSPSPG,&
    & PSPVOR=PSPVOR2,PSPDIV=PSPDIV2,PSPT=PSPT2,PSPSPD=PSPSPD2,&
    & PSPSVD=PSPSVD2,PSPSP=PSPSP2,&
    & LDFULLM=LLONEM,LDNEEDPS=.TRUE.) 

    pspsp(:)=pspsp2(:)
    do icnt2=1,nspec2
      do icnt1=1,nflevl
        pspvor(icnt1,icnt2)=pspvor2(icnt2,icnt1)
        pspdiv(icnt1,icnt2)=pspdiv2(icnt2,icnt1)
        pspt(icnt1,icnt2)=pspt2(icnt2,icnt1)
        pspspd(icnt1,icnt2)=pspspd2(icnt2,icnt1)
        pspsvd(icnt1,icnt2)=pspsvd2(icnt2,icnt1)
      enddo
    enddo


endif

#endif

    IF (LHOOK) CALL DR_HOOK('SPCM_SIMPLE_utile',1,ZHOOK_HANDLE3)

    if (lhook) call dr_hook('SPCM_SIMPLE_transferts2b',0,zhook_handle2)
    !$acc end data
    if (lhook) call dr_hook('SPCM_SIMPLE_transferts2b',1,zhook_handle2)
    if (lhook) call dr_hook('SPCM_SIMPLE_transferts2a',0,zhook_handle2)
    !$acc end data
    if (lhook) call dr_hook('SPCM_SIMPLE_transferts2a',1,zhook_handle2)

#if defined(_OPENACC)
    pspsp(:)=pspsp2(:)
    do icnt2=1,nspec2
      do icnt1=1,nflevl
        pspvor(icnt1,icnt2)=pspvor2(icnt2,icnt1)
        pspdiv(icnt1,icnt2)=pspdiv2(icnt2,icnt1)
        pspt(icnt1,icnt2)=pspt2(icnt2,icnt1)
        pspspd(icnt1,icnt2)=pspspd2(icnt2,icnt1)
        pspsvd(icnt1,icnt2)=pspsvd2(icnt2,icnt1)
      enddo
    enddo

#endif

ENDIF

IF (ALLOCATED(ZSPVORG)) DEALLOCATE(ZSPVORG)
IF (ALLOCATED(ZSPDIVG)) DEALLOCATE(ZSPDIVG)
IF (ALLOCATED(ZSPTG))   DEALLOCATE(ZSPTG)
IF (ALLOCATED(ZSPSPDG)) DEALLOCATE(ZSPSPDG)
IF (ALLOCATED(ZSPSVDG)) DEALLOCATE(ZSPSVDG)
IF (ALLOCATED(ZSPSPG))  DEALLOCATE(ZSPSPG)

IF (ALLOCATED(ZSPVORG2)) DEALLOCATE(ZSPVORG2)
IF (ALLOCATED(ZSPDIVG2)) DEALLOCATE(ZSPDIVG2)
IF (ALLOCATED(ZSPTG2))   DEALLOCATE(ZSPTG2)
IF (ALLOCATED(ZSPSPDG2)) DEALLOCATE(ZSPSPDG2)
IF (ALLOCATED(ZSPSVDG2)) DEALLOCATE(ZSPSVDG2)
IF (ALLOCATED(ZSPSPG2))  DEALLOCATE(ZSPSPG2)

if (allocated(pspvor2)) DEALLOCATE(PSPVOR2)
if (allocated(pspdiv2)) DEALLOCATE(PSPDIV2)
if (allocated(pspt2)) DEALLOCATE(PSPT2)
if (allocated(pspspd2)) DEALLOCATE(PSPSPD2)
if (allocated(pspsvd2)) DEALLOCATE(PSPSVD2)
if (allocated(pspsp2)) DEALLOCATE(PSPSP2)

IF (ALLOCATED(ZSPTNDSI_VORG)) DEALLOCATE(ZSPTNDSI_VORG)
IF (ALLOCATED(ZSPTNDSI_DIVG)) DEALLOCATE(ZSPTNDSI_DIVG)
IF (ALLOCATED(ZSPTNDSI_TG))   DEALLOCATE(ZSPTNDSI_TG)
#if defined(_OPENACC)
  !$acc exit data delete(param_mxture)
!!!!!!  !$acc end data
  if (allocated(param_mxture)) deallocate(param_mxture)
#endif

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SPCM_SIMPLE',1,ZHOOK_HANDLE)

END SUBROUTINE SPCM_SIMPLE

