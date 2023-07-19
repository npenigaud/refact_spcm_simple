#if defined(_OPENACC)
SUBROUTINE SPCSI_STR(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KSPEC2V,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
 & taillec,param_mxture)
#else
SUBROUTINE SPCSI_STR(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KSPEC2V,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
 & simit,simot)
#endif

!**** *SPCSI* - SPECTRAL SPACE SEMI-IMPLICIT COMPUTATIONS FOR HYD MODEL.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPCSI(..)

!        Explicit arguments :
!        --------------------  KM      - Zonal wavenumber 
!                              KMLOC   - Zonal wavenumber (DM-local numbering)
!                              KSTA    - First column processed
!                              KEND    - Last column processed
!                              PSPVORG - Vorticity columns
!                              PSPDIVG - Divergence columns
!                              PSPTG   - Temperature columns
!                              PSPSPG  - Surface Pressure
!                              PSPTNDSI_VORG - [D vor/Dt]_SI
!                              PSPTNDSI_DIVG - [D div/Dt]_SI
!                              PSPTNDSI_TG   - [D T/Dt]_SI

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-11-24 (before 1997 spcsi.F was part of spc.F)

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      N.Wedi        08-Mar-2005 remove mass correction      
!      K.Yessad 09-Dec-2004: move mass correction in SPCMASCOR + cleanings.
!      K. Yessad 15-May-2006: memory optimisations for stretched geometry
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      K. Yessad (Aug 2009): remove LSITRIC option.
!      F. Voitus: add DDH diagnostics.
!      T.Wilhelmsson 09-09-25: Remove LFULLM requirement for LIMPF
!      K. Yessad (Feb 2012): tests on LL3D, LLDOSI in the caller, simplifications.
!      P. Marguinaud (Nov 2012): Fix unallocated array arguments
!      P. Marguinaud (Sep 2012) : Make PSPAUXG optional
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      O. Marsden (May 2016): Removed redundant geometry arguments
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0       , ONLY : MYSETW, MYSETV, MYSETN
USE YOMDYN       , ONLY : TDYN
USE YOMLDDH      , ONLY : TLDDH
USE YOMRIP       , ONLY : TRIP
USE YOMCST       , ONLY : TCST

#if defined(_OPENACC)
use cublas
#endif

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(kspec2V,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
#if defined(_OPENACC)
integer(kind=jpim)               :: taillec
real(kind=jprb)   ,intent(in)    :: param_mxture(:,:,:)
real(kind=jprb)                  :: pa(taillec)
real(kind=jprb)                  :: pb(taillec)
real(kind=jprb)                  :: pc(taillec)
real(kind=jprb)                  :: entree(taillec)
real(kind=jprb)                  :: sortie(taillec)
real(kind=JPRB)                  :: zsdivpl(ydgeometry%yrdim%NSMAX+1,ydgeometry%yrdimv%NFLEVG,2,ydgeometry%yrdim%nump)
real(kind=JPRB)                  :: zspdivpl(ydgeometry%yrdim%NSMAX+1,ydgeometry%yrdimv%NFLEVG,2,ydgeometry%yrdim%nump)
#else
real(kind=jprb)   ,intent(in)    :: simit(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
real(kind=jprb)   ,intent(in)    :: simot(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
#endif


REAL(KIND=JPRB) :: ZSDIVP (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSPDIVP(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZSDIV  (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZHELP  (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZST    (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSP    (KSPEC2V)


INTEGER(KIND=JPIM) :: IN, IOFF, JLEV, JSP  
INTEGER(KIND=JPIM) :: JMLOC, IM, ISTA, IEND
REAL(KIND=JPRB) :: ZBDT, ZBDT2

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE,zhook_handle2
real(kind=jprb)    :: zbdtaux

!     ------------------------------------------------------------------

#include "mxmaop.h"
#include "mxptma.h"
#include "mxture.h"
#include "mxturs.h"
#include "abor1.intfb.h"
#include "sigam_sp_openmp.intfb.h"
#include "sitnu_sp_openmp.intfb.h"
#include "spcsidg_part0.intfb.h"
#include "spcsidg_part1.intfb.h"
#include "spcsidg_part2.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPCSI_STR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,   &
&  YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE( LSIDG=>YDDYN%LSIDG,                                 &
&   LRSIDDH=>YDLDDH%LRSIDDH,SIMI=>YDDYN%SIMI,SIMO=>YDDYN%SIMO)
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------

!*       1.    MEMORY TRANSFER.
!              ----------------
IF (LRSIDDH) THEN
  ! DDH memory transfer
  PSPTNDSI_DIVG=-PSPDIVG
  PSPTNDSI_TG  =-PSPTG  
  !the case of surface pressure has not been treated yet
ENDIF

call flush(0)

!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

!*        2.1  Preliminary initialisations.

IOFF=ydgeometry%yrmp%NPTRSVF(MYSETV)-1

ZBDT=yddyn%RBTS2*ydrip%TDT
ZBDT2=(ZBDT*ydgeometry%yrgem%RSTRET)**2

!*        2.2  OpenACC memory
#if defined(_OPENACC)
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts1',0,ZHOOK_HANDLE2)
!$acc data create(zsdivp,zspdivp,zsdiv,zhelp,zst,zsp,zsdivpl,zspdivpl,pa,pb,pc,entree,sortie)
!$acc data present(YDGEOMETRY,YDGEOMETRY%YRLAP,YDGEOMETRY%YRLAP%NVALUE,YDGEOMETRY%YRLAP%RLAPIN,YDGEOMETRY%YRLAP%RLAPDI,ydgeometry%yrdimv%NFLEVG,ydgeometry%yrdim%NSMAX,YDDYN,YDDYN%SIVP)
!$acc data present(pspdivg,psptg,pspspg,YDRIP,ydgeometry%yrmp%NPTRMF,ydgeometry%YRLAP,ydgeometry%yrmp%NSPSTAF,ydgeometry%yrgem%RSTRET)
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts1',1,ZHOOK_HANDLE2)
#endif


!*        2.3  Computes right-hand side of Helmholtz equation.
#if defined(_OPENACC)
CALL SIGAM_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,ydgeometry%yrdimv%NFLEVG,KSPEC2V,ZSDIV,PSPTG(:,:),PSPSPG(1:kspec2v),zsdivpl,zspdivpl) 
#else
CALL SIGAM_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,ydgeometry%yrdimv%NFLEVG,KSPEC2V,ZSDIV,PSPTG(:,:),PSPSPG(1:kspec2v)) !!!ispcol remplacé par kspec2V
#endif
IF (LSIDG) THEN

if (lhook) CALL DR_HOOK('SPCSI_sidg0',0,zhook_handle2)
#if defined(_OPENACC)
  !$acc parallel private(im,ista,iend) default(none)
  !$acc loop gang collapse(2) 
  DO JMLOC=ydgeometry%yrmp%NPTRMF(MYSETN), ydgeometry%yrmp%NPTRMF(MYSETN+1)-1
    do jlev=1,ydgeometry%yrdimv%NFLEVG
  
      IM=ydgeometry%YRLAP%MYMS(jMLOC)
      ISTA=ydgeometry%yrmp%NSPSTAF(IM)
      IEND=ISTA+2*(ydgeometry%yrdim%NSMAX+1-IM)-1
      IF (IM > 0) THEN
        !$acc loop vector private(in) 
        DO JSP=ISTA,IEND
          IN=ydgeometry%YRLAP%NVALUE(JSP+IOFF)
          ZSDIV(jsp,JLEV)=ydgeometry%YRLAP%RLAPIN(IN)*PSPDIVG(jsp,JLEV)-ZBDT*ZSDIV(jsp,JLEV)
        ENDDO

      ELSE

        !$acc loop vector private(in)
        DO JSP=ISTA,IEND
          IN=ydgeometry%YRLAP%NVALUE(JSP+IOFF)
          ZSDIV(jsp,JLEV)=PSPDIVG(jsp,JLEV)-ZBDT*ydgeometry%YRLAP%RLAPDI(IN)*ZSDIV(jsp,JLEV)
        ENDDO

      ENDIF
    enddo 
  ENDDO
  !$acc end parallel
#else
  DO JMLOC=ydgeometry%yrmp%NPTRMF(MYSETN), ydgeometry%yrmp%NPTRMF(MYSETN+1)-1
    CALL SPCSIDG_PART0 (YDGEOMETRY, YDDYN, YDRIP, KSPEC2V, JMLOC, ZSDIV, PSPDIVG)
  ENDDO
#endif
if (lhook) CALL DR_HOOK('SPCSI_sidg0',1,zhook_handle2) 
ELSE

  ! Case of No Stretching

if (lhook) call dr_hook('SPCSI_boucle1',0,zhook_handle2)
#if defined(_OPENACC)
!$acc PARALLEL loop collapse(2) PRIVATE(JSP,JLEV,IN,zbdtaux) default(none)
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN,zbdtaux)
#endif
  DO JLEV=1,ydgeometry%yrdimv%NFLEVG
    DO JSP=1,kspec2v
      IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
      zbdtaux=ZBDT*YDGEOMETRY%YRLAP%RLAPDI(IN)
      ZSDIV(JSP,jlev)=PSPDIVG(JSP,jlev)-zbdtaux*ZSDIV(JSP,jlev)
    ENDDO
  ENDDO
#if defined(_OPENACC)
!$acc END PARALLEL
#else
!$OMP END PARALLEL DO
#endif
if (lhook) call dr_hook('SPCSI_boucle1',1,zhook_handle2)
ENDIF


!*        2.4  Solve Helmholtz equation

!           Current space --> vertical eigenmodes space.

IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop1',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
   !$acc host_data use_device(SIMI,ZSDIV,ZSDIVP)
     CALL cublasDgemm('N','T',kspec2v,ydgeometry%yrdimv%NFLEVG,ydgeometry%yrdimv%NFLEVG,1.0_JPRB,&
      &zsdiv,kspec2v,SIMI,ydgeometry%yrdimv%NFLEVG,0.0_JPRB,ZSDIVP(1,1),kspec2v)  !!!!ispcol remplacé par kspec2V, ksta par 1, suppression de (1,1)?
   !$acc end host_data
   !$acc wait
!!           pa       kad       pb         kbd    pc       kca      kar    kac   kbc
!!!!!$acc host_data use_device(yddyn%yddyn%SIMI,ydgeometry%yrdimv%NFLEVG)
!!!!call mxmaop(zsdiv,1,kspec2v,yddyn%yddyn%SIMI,1,ydgeometry%yrdimv%NFLEVG,zsdivp,1,kspec2v,kspec2v,ydgeometry%yrdimv%NFLEVG,ydgeometry%yrdimv%NFLEVG)!! mxmaop fait une transposition en _OPENACC
!!!!!$acc end host_data
#else
CALL MXMAOP(zsdiv,1,kspec2v,simit,1,ydgeometry%yrdimv%NFLEVG,ZSDIVP,1,kspec2v,kspec2v,ydgeometry%yrdimv%NFLEVG,ydgeometry%yrdimv%NFLEVG)  
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop1',1,ZHOOK_HANDLE2)
IF (LSIDG) THEN

if (lhook) CALL DR_HOOK('SPCSI_sidg1',0,zhook_handle2)
#if defined(_OPENACC) 
    CALL SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V,ZSDIVP,ZSPDIVP,taillec,zsdivpl,zspdivpl,pa,pb,pc,entree,sortie,param_mxture,ydgeometry%yrmp%NPTRMF(mysetn),ydgeometry%yrmp%NPTRMF(mysetn+1)-1)
#else
    CALL SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, ZSDIVP,ZSPDIVP,ydgeometry%yrmp%NPTRMF(mysetn),ydgeometry%yrmp%NPTRMF(mysetn+1)-1)
#endif
if (lhook) CALL DR_HOOK('SPCSI_sidg1',1,zhook_handle2) 

ELSE
  !                 Inversion of a diagonal system (Helmholtz equation)
  !                 --> (yddyn%SIMI*DIVprim(t+dt)).

if (lhook) call dr_hook('SPCSI_boucle2',0,zhook_handle2)
#if defined(_OPENACC)
    !$acc parallel private(JSP,JLEV,zbdtaux) default(none)
    !$acc loop gang
#else
    !$omp parallel do private(jsp,jlev,zbdtaux) !!pas de parallelisation code initial
#endif
    DO JLEV=1,ydgeometry%yrdimv%NFLEVG
      zbdtaux=ZBDT2*YDDYN%SIVP(JLEV)
      !$acc loop vector
      DO JSP=1,KSPEC2V
        ZSPDIVP(JSP,jlev)=ZSDIVP(JSP,jlev)&
         & /(1.0_JPRB-zbdtaux*YDGEOMETRY%YRLAP%RLAPDI(YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)))  
      ENDDO
    ENDDO
#if defined(_OPENACC)
    !$acc end parallel
#else
    !$omp end parallel do !!pas de parallelisation de cette boucle dans le code initial
#endif

if (lhook) call dr_hook('SPCSI_boucle2',1,zhook_handle2)

ENDIF

!           Vertical eigenmodes space --> current space.

if (lhook) call dr_hook('SPCSI_mxmaop2',0,zhook_handle2)
#if defined(_OPENACC)
!$acc host_data use_device(SIMO,ZSPDIVP,PSPDIVG)
CALL cublasDgemm('N','T',kspec2v,ydgeometry%yrdimv%NFLEVG,ydgeometry%yrdimv%NFLEVG,1.0_JPRB,&    !!ispcol remplacé par kspec2v
&ZSPDIVP(1,1),kspec2v,SIMO,ydgeometry%yrdimv%NFLEVG,0.0_JPRB,PSPDIVG(1,1),kspec2v) !!2ksta remplacés par 1,pourrait partir
!$acc end host_data
!$acc wait
!!           pa        kad      pb            kbd    pc       kca     kar    kac   kbc
!!!$acc host_data use_device(yddyn%simo,ydgeometry%yrdimv%NFLEVG)
!!call mxmaop(zspdivp,1,kspec2v,yddyn%simo,1,ydgeometry%yrdimv%NFLEVG,pspdivg,1,kspec2v,kspec2v,ydgeometry%yrdimv%NFLEVG,ydgeometry%yrdimv%NFLEVG) !!il y avait une erreur ici
!!!$acc end host_data
#else
CALL MXMAOP(zspdivp,1,kspec2v,simot,1,ydgeometry%yrdimv%NFLEVG,PSPDIVG,1,kspec2V,kspec2v,ydgeometry%yrdimv%NFLEVG,ydgeometry%yrdimv%NFLEVG)  
#endif
if (lhook) call dr_hook('SPCSI_mxmaop2',1,zhook_handle2)

IF (LSIDG) THEN

if (lhook) CALL DR_HOOK('SPCSI_sidg2',0,zhook_handle2)
#if defined(_OPENACC) 
    CALL SPCSIDG_PART2 (YDGEOMETRY, KSPEC2V,PSPDIVG,ZHELP,zsdivpl,zspdivpl,ydgeometry%yrmp%NPTRMF(mysetn),ydgeometry%yrmp%NPTRMF(mysetn+1)-1)
#else
    CALL SPCSIDG_PART2 (YDGEOMETRY, KSPEC2V, PSPDIVG, ZHELP,ydgeometry%yrmp%NPTRMF(mysetn),ydgeometry%yrmp%NPTRMF(mysetn+1)-1)
#endif
if (lhook) CALL DR_HOOK('SPCSI_sidg2',1,zhook_handle2) 

ELSE

  !       ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GMBAR**2 * DIVprim(t+dt)) .
if (lhook) call dr_hook('SPCSI_boucle3',0,zhook_handle2)

#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JSP,JLEV) default(none)
!$acc loop gang
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
#endif
  DO JLEV=1,ydgeometry%yrdimv%NFLEVG
    !$acc loop vector
    DO JSP=1,kspec2v
      ZHELP(JSP,jlev)=PSPDIVG(JSP,jlev)*ydgeometry%yrgem%RSTRET*ydgeometry%yrgem%RSTRET
    ENDDO
  ENDDO
#if defined(_OPENACC)
!$acc END PARALLEL
#else
!$OMP END PARALLEL DO
#endif
if (lhook) call dr_hook('SPCSI_boucle3',1,zhook_handle2)

ENDIF

!       If LSIDG:
!         (GM**2 * DIVprim(t+dt)) --> [ tau * (GM**2 * DIVprim(t+dt)) ]
!                                 and [  nu * (GM**2 * DIVprim(t+dt)) ]
!       or if not LSIDG:
!         (GMBAR**2 * DIVprim(t+dt)) --> [ tau * (GMBAR**2 * DIVprim(t+dt)) ]
!                                    and [  nu * (GMBAR**2 * DIVprim(t+dt)) ]

#if defined(_OPENACC)
CALL SITNU_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,ydgeometry%yrdimv%NFLEVG,KSPEC2V,ZHELP,ZST,ZSP,zsdivpl,zspdivpl)
#else
CALL SITNU_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,ydgeometry%yrdimv%NFLEVG,KSPEC2V,ZHELP,ZST,ZSP)
#endif

!*       2.5  Increment Temperature and surface pressure

if (lhook) call dr_hook('SPCSI_boucle4',0,zhook_handle2)
#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JSP,JLEV) default(none)
!$acc loop gang vector collapse(2)
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
#endif
DO JLEV=1,ydgeometry%yrdimv%NFLEVG
  DO JSP=1,kspec2v
    PSPTG(JSP,jlev)=PSPTG(JSP,jlev)-ZBDT*ZST(JSP,jlev)
  ENDDO
ENDDO
#if defined(_OPENACC)
!$acc end parallel
#else
!$OMP END PARALLEL DO
#endif

if (lhook) call dr_hook('SPCSI_boucle4',1,zhook_handle2)
if (lhook) call dr_hook('SPCSI_boucle4_b',0,zhook_handle2)


#if defined(_OPENACC)
!$acc parallel private(jsp) default(none)
!$acc loop gang vector
#else
!$OMP PARALLEL DO PRIVATE(JSP)
#endif
do jsp=1,kspec2v
  pspspg(jsp)=pspspg(jsp)-zbdt*zsp(jsp)
enddo
#if defined(_OPENACC)
!$acc end parallel
#else
!$OMP END PARALLEL DO
#endif

if (lhook) call dr_hook('SPCSI_boucle4_b',1,zhook_handle2)

#if defined(_OPENACC)
if (lhook) call dr_hook('SPCSI_transferts2',0,zhook_handle2)
!$acc end data
!$acc end data 
!$acc end data


if (lhook) call dr_hook('SPCSI_transferts2',1,zhook_handle2)
#endif

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF SI TERM AT t+dt FOR DDH.
!              ---------------------------------------

IF (LRSIDDH) THEN
  PSPTNDSI_DIVG=PSPTNDSI_DIVG + PSPDIVG
  PSPTNDSI_TG=PSPTNDSI_TG + PSPTG
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPCSI_STR',1,ZHOOK_HANDLE)
END SUBROUTINE SPCSI_STR

