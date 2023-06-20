#if defined(_OPENACC)
SUBROUTINE SPCSI_STR(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KSPEC2V,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
 & param_mxture)
#else
SUBROUTINE SPCSI_STR(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KSPEC2V,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
 & SIMIT,SIMOT)
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
!      K. Yessad (Feb 2012): tests on LL3D, LLDOSI in the CALLer, simplifications.
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
USE CUBLAS
#endif

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
#if defined(_OPENACC)
REAL(KIND=JPRB)                  :: ZSDIVPL(YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2,500)
REAL(KIND=JPRB)                  :: ZSPDIVPL(YDGEOMETRY%YRDIM%NSMAX+1,YDGEOMETRY%YRDIMV%NFLEVG,2,500)
REAL(KIND=jprb)   ,INTENT(IN)    :: param_mxture(:,:,:)
#else
REAL(KIND=jprb)   ,INTENT(IN)    :: SIMIT(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=jprb)   ,INTENT(IN)    :: SIMOT(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
#endif


REAL(KIND=JPRB) :: ZSDIVP (KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSPDIVP(KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZSDIV  (KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZHELP  (KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZST    (KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSP    (KSPEC2V)


INTEGER(KIND=JPIM) :: IN, IOFF, JLEV, JSP  
INTEGER(KIND=JPIM) :: JMLOC, IM, ISTA, IEND
REAL(KIND=JPRB) :: ZBDT, ZBDT2

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE,ZHOOK_HANDLE2
REAL(KIND=JPRB)    :: ZBDTAUX

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
& YDLAP=>YDGEOMETRY%YRLAP, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, LSIDG=>YDDYN%LSIDG, RBTS2=>YDDYN%RBTS2,                                &
& SIMI=>YDDYN%SIMI, SIMO=>YDDYN%SIMO, SIVP=>YDDYN%SIVP, RSTRET=>YDGEM%RSTRET, LRSIDDH=>YDLDDH%LRSIDDH,  &
& NPTRSV=>YDMP%NPTRSV, NPTRSVF=>YDMP%NPTRSVF, NSPEC2V=>YDMP%NSPEC2V, NSPEC2VF=>YDMP%NSPEC2VF,           &
& TDT=>YDRIP%TDT,NPTRMF=>YDMP%NPTRMF,NSPSTAF=>YDMP%NSPSTAF,NSMAX=>YDDIM%NSMAX)
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

!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

!*        2.1  Preliminary initialisations.

IOFF=NPTRSVF(MYSETV)-1

ZBDT=RBTS2*TDT
ZBDT2=(ZBDT*RSTRET)**2

!*        2.2  OpenACC memory
#if defined(_OPENACC)
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts1',0,ZHOOK_HANDLE2)
!$ACC DATA CREATE(ZSDIVP,ZSPDIVP,ZSDIV,ZHELP,ZST,ZSP,ZSDIVPL,ZSPDIVPL)
!$ACC DATA PRESENT(YDGEOMETRY,YDGEOMETRY%YRLAP,YDGEOMETRY%YRLAP%NVALUE,YDGEOMETRY%YRLAP%RLAPIN,YDGEOMETRY%YRLAP%RLAPDI,NFLEVG,NSMAX,YDDYN,YDDYN%SIVP,RSTRET)
!$ACC DATA PRESENT(PSPDIVG,PSPTG,PSPSPG,YDRIP,NPTRMF,MYSETN)
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts1',1,ZHOOK_HANDLE2)
#endif


!*        2.3  Computes right-hand side of Helmholtz equation.

CALL SIGAM_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,KSPEC2V,ZSDIV,PSPTG(:,:),PSPSPG(1:KSPEC2V)) 

IF (LSIDG) THEN

IF (LHOOK) CALL DR_HOOK('SPCSI_sidg0',0,ZHOOK_HANDLE2)

CALL SPCSIDG_PART0(YDGEOMETRY, YDDYN, YDRIP, KSPEC2V, ZSDIV,&
  &PSPDIVG,NPTRMF(MYSETN),NPTRMF(MYSETN+1)-1)

IF (LHOOK) CALL DR_HOOK('SPCSI_sidg0',1,ZHOOK_HANDLE2) 

ELSE

  ! Case of No Stretching

IF (LHOOK) CALL DR_HOOK('SPCSI_boucle1',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(JSP,JLEV,IN,ZBDTAUX) DEFAULT(NONE)
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN,ZBDTAUX)
#endif
  DO JLEV=1,NFLEVG
    DO JSP=1,KSPEC2V
      IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
      ZBDTAUX=ZBDT*YDGEOMETRY%YRLAP%RLAPDI(IN)
      ZSDIV(JSP,JLEV)=PSPDIVG(JSP,JLEV)-ZBDTAUX*ZSDIV(JSP,JLEV)
    ENDDO
  ENDDO
#if defined(_OPENACC)
!$ACC END PARALLEL
#else
!$OMP END PARALLEL DO
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle1',1,ZHOOK_HANDLE2)
ENDIF


!*        2.4  Solve Helmholtz equation

!           Current space --> vertical eigenmodes space.

IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop1',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
CALL MXMAOP(ZSDIV,1,KSPEC2V,SIMI,1,NFLEVG,ZSDIVP,1,KSPEC2V,KSPEC2V,NFLEVG,NFLEVG)
#else
CALL MXMAOP(ZSDIV,1,KSPEC2V,SIMIT,1,NFLEVG,ZSDIVP,1,KSPEC2V,KSPEC2V,NFLEVG,NFLEVG)  
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop1',1,ZHOOK_HANDLE2)

IF (LSIDG) THEN

IF (LHOOK) CALL DR_HOOK('SPCSI_sidg1',0,ZHOOK_HANDLE2)
#if defined(_OPENACC) 
    CALL SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V,ZSDIVP,ZSPDIVP,ZSDIVPL,ZSPDIVPL,param_mxture,NPTRMF(MYSETN),NPTRMF(MYSETN+1)-1)
#else
    CALL SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, ZSDIVP,ZSPDIVP,NPTRMF(MYSETN),NPTRMF(MYSETN+1)-1)
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_sidg1',1,ZHOOK_HANDLE2) 

ELSE
  !                 Inversion of a diagonal system (Helmholtz equation)
  !                 --> (SIMI*DIVprim(t+dt)).

IF (LHOOK) CALL DR_HOOK('SPCSI_boucle2',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
    !$ACC PARALLEL PRIVATE(JSP,JLEV,ZBDTAUX) DEFAULT(NONE)
    !$ACC LOOP GANG
#else
    !$OMP PARALLEL DO PRIVATE(JSP,JLEV,ZBDTAUX) 
#endif
    DO JLEV=1,NFLEVG
      ZBDTAUX=ZBDT2*YDDYN%SIVP(JLEV)
      !$ACC LOOP VECTOR
      DO JSP=1,KSPEC2V
        ZSPDIVP(JSP,JLEV)=ZSDIVP(JSP,JLEV)&
         & /(1.0_JPRB-ZBDTAUX*YDGEOMETRY%YRLAP%RLAPDI(YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)))  
      ENDDO
    ENDDO
#if defined(_OPENACC)
    !$ACC END PARALLEL
#else
    !$OMP END PARALLEL DO
#endif

IF (LHOOK) CALL DR_HOOK('SPCSI_boucle2',1,ZHOOK_HANDLE2)

ENDIF

!           Vertical eigenmodes space --> current space.

IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop2',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
CALL MXMAOP(ZSPDIVP,1,KSPEC2V,YDDYN%simo,1,NFLEVG,PSPDIVG,1,KSPEC2V,KSPEC2V,NFLEVG,NFLEVG)
#else
CALL MXMAOP(ZSPDIVP,1,KSPEC2V,SIMOT,1,NFLEVG,PSPDIVG,1,KSPEC2V,KSPEC2V,NFLEVG,NFLEVG)  
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop2',1,ZHOOK_HANDLE2)

IF (LSIDG) THEN

IF (LHOOK) CALL DR_HOOK('SPCSI_sidg2',0,ZHOOK_HANDLE2)
#if defined(_OPENACC) 
    CALL SPCSIDG_PART2 (YDGEOMETRY, KSPEC2V,PSPDIVG,ZHELP,ZSDIVPL,ZSPDIVPL,NPTRMF(MYSETN),NPTRMF(MYSETN+1)-1)
#else
    CALL SPCSIDG_PART2 (YDGEOMETRY, KSPEC2V, PSPDIVG, ZHELP,NPTRMF(MYSETN),NPTRMF(MYSETN+1)-1)
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_sidg2',1,ZHOOK_HANDLE2) 

ELSE

  !       ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GMBAR**2 * DIVprim(t+dt)) .
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle3',0,ZHOOK_HANDLE2)

#if defined(_OPENACC)
!$ACC PARALLEL PRIVATE(JSP,JLEV) DEFAULT(NONE)
!$ACC LOOP GANG
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
#endif
  DO JLEV=1,NFLEVG
    !$ACC LOOP VECTOR
    DO JSP=1,KSPEC2V
      ZHELP(JSP,JLEV)=PSPDIVG(JSP,JLEV)*RSTRET*RSTRET
    ENDDO
  ENDDO
#if defined(_OPENACC)
!$ACC END PARALLEL
#else
!$OMP END PARALLEL DO
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle3',1,ZHOOK_HANDLE2)

ENDIF

!       IF LSIDG:
!         (GM**2 * DIVprim(t+dt)) --> [ tau * (GM**2 * DIVprim(t+dt)) ]
!                                 and [  nu * (GM**2 * DIVprim(t+dt)) ]
!       or IF not LSIDG:
!         (GMBAR**2 * DIVprim(t+dt)) --> [ tau * (GMBAR**2 * DIVprim(t+dt)) ]
!                                    and [  nu * (GMBAR**2 * DIVprim(t+dt)) ]


CALL SITNU_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,KSPEC2V,ZHELP,ZST,ZSP)

!*       2.5  Increment Temperature and surface pressure

IF (LHOOK) CALL DR_HOOK('SPCSI_boucle4',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
!$ACC PARALLEL PRIVATE(JSP,JLEV) DEFAULT(NONE)
!$ACC LOOP COLLAPSE(2)
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
#endif
DO JLEV=1,NFLEVG
  DO JSP=1,KSPEC2V
    PSPTG(JSP,JLEV)=PSPTG(JSP,JLEV)-ZBDT*ZST(JSP,JLEV)
  ENDDO
ENDDO
#if defined(_OPENACC)
!$ACC END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

#if defined(_OPENACC)
!$ACC PARALLEL LOOP PRIVATE(JSP) DEFAULT(NONE)
#else
!$OMP PARALLEL DO PRIVATE(JSP)
#endif
do JSP=1,KSPEC2V
  PSPSPG(JSP)=PSPSPG(JSP)-ZBDT*ZSP(JSP)
ENDdo
#if defined(_OPENACC)
!$ACC END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

IF (LHOOK) CALL DR_HOOK('SPCSI_boucle4',1,ZHOOK_HANDLE2)

#if defined(_OPENACC)
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts2',0,ZHOOK_HANDLE2)
!$ACC END DATA
!$ACC END DATA 
!$ACC END DATA


IF (LHOOK) CALL DR_HOOK('SPCSI_transferts2',1,ZHOOK_HANDLE2)
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

