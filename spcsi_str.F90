#if defined(_OPENACC)
SUBROUTINE SPCSI_STR(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KSPEC2V,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
 & zsdiv,zhelp,zsp,zst,zsdivp,zspdivp,zsphi,zout)
#else
SUBROUTINE SPCSI_STR(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KSPEC2V,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG)
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
#if defined(_OPENACC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(kspec2V,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(1,1)!!chgt ici
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(1,1)!!là 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(1,1)  !!et là
REAL(KIND=JPRB)   ,intent(inout) :: ZSDIV  (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,intent(inout) :: ZHELP  (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,intent(inout) :: ZST    (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,intent(inout) :: ZSP    (KSPEC2V)
REAL(KIND=JPRB)   ,intent(inout) :: ZSDIVP (max(kspec2v,YDGEOMETRY%YRMP%NSPEC2VF),YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,intent(inout) :: ZSPDIVP(max(kspec2v,YDGEOMETRY%YRMP%NSPEC2VF),YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,intent(inout) :: ZSPHI  (kspec2v,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB)   ,intent(inout) :: ZOUT  (kspec2v,0:YDGEOMETRY%YRDIMV%NFLEVG)
#else
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)


!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSDIVP (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB) :: ZSPDIVP(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)

REAL(KIND=JPRB) :: ZSDIV  (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB) :: ZHELP  (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB) :: ZST    (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB) :: ZSP    (KSPEC2V)
#endif

INTEGER(KIND=JPIM) :: IN, IOFF, JLEV, JSP  
INTEGER(KIND=JPIM) :: JMLOC, IM, ISTA, IEND
REAL(KIND=JPRB) :: ZBDT, ZBDT2

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE,zhook_handle2
real(kind=jprb)    :: intermediaire

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
!$acc data present(YDGEOMETRY,YDGEOMETRY%YRLAP,YDGEOMETRY%YRLAP%NVALUE,YDGEOMETRY%YRLAP%RLAPIN,YDGEOMETRY%YRLAP%RLAPDI,nflevg,nsmax,YDDYN,YDDYN%SIVP,rstret)
!$acc data present(pspdivg,psptg,pspspg)
!$acc data present(zsdiv,zhelp,zsp,zst,zsdivp,zspdivp,zsphi,zout)
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts1',1,ZHOOK_HANDLE2)
#endif


!*        2.3  Computes right-hand side of Helmholtz equation.

#if defined(_OPENACC)
CALL SIGAM_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,KSPEC2V,ZSDIV,PSPTG(:,:),PSPSPG(1:kspec2v),zsphi,zout) !!!ispcol remplacé par kspec2V
#else
CALL SIGAM_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,KSPEC2V,ZSDIV,PSPTG,PSPSPG)
#endif

IF (LSIDG) THEN

  DO JMLOC=NPTRMF(MYSETN), NPTRMF(MYSETN+1)-1
    CALL SPCSIDG_PART0 (YDGEOMETRY, YDDYN, YDRIP, KSPEC2V, JMLOC, ZSDIV, PSPDIVG)
  ENDDO

ELSE

  ! Case of No Stretching

if (lhook) call dr_hook('SPCSI_boucle1',0,zhook_handle2)
#if defined(_OPENACC)
!$acc PARALLEL loop collapse(2) PRIVATE(JSP,JLEV,IN,intermediaire) default(none)
  DO JLEV=1,nflevg
    DO JSP=1,kspec2v
      IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
      intermediaire=ZBDT*YDGEOMETRY%YRLAP%RLAPDI(IN)
      ZSDIV(JSP,jlev)=PSPDIVG(JSP,jlev)-intermediaire*ZSDIV(JSP,jlev)
    ENDDO
  ENDDO
!$acc END PARALLEL
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
  DO JSP=1,KSPEC2V
    DO JLEV=1,NFLEVG
      IN=YDLAP%NVALUE(JSP+IOFF)
      ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-ZBDT*YDLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
#endif
ENDIF

if (lhook) call dr_hook('SPCSI_boucle1',1,zhook_handle2)

!*        2.4  Solve Helmholtz equation

!           Current space --> vertical eigenmodes space.

IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop1',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
   !$acc host_data use_device(SIMI,ZSDIV,ZSDIVP)
     CALL cublasDgemm('N','T',kspec2v,nflevg,nflevg,1.0_JPRB,&
      &zsdiv,kspec2v,simi,nflevg,0.0_JPRB,ZSDIVP(1,1),kspec2v)  !!!!ispcol remplacé par kspec2V, ksta par 1, suppression de (1,1)?
   !$acc end host_data
   !$acc wait
#else
CALL MXMAOP(SIMI,1,NFLEVG,ZSDIV,1,NFLEVG,ZSDIVP,1,NFLEVG,NFLEVG,NFLEVG,KSPEC2V)  
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop1',1,ZHOOK_HANDLE2)

IF (LSIDG) THEN

  DO JMLOC=NPTRMF(MYSETN), NPTRMF(MYSETN+1)-1
    CALL SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, JMLOC, ZSDIVP, ZSPDIVP)
  ENDDO

ELSE
  !                 Inversion of a diagonal system (Helmholtz equation)
  !                 --> (SIMI*DIVprim(t+dt)).

if (lhook) call dr_hook('SPCSI_boucle2',0,zhook_handle2)
#if defined(_OPENACC)
    !$acc parallel private(JSP,JLEV,intermediaire) default(none)
    !$acc loop gang
    DO JLEV=1,NFLEVG
      intermediaire=ZBDT2*YDDYN%SIVP(JLEV)
      !$acc loop vector
      DO JSP=1,KSPEC2V
        ZSPDIVP(JSP,jlev)=ZSDIVP(JSP,jlev)&
         & /(1.0_JPRB-intermediaire*YDGEOMETRY%YRLAP%RLAPDI(YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)))  
      ENDDO
    ENDDO
    !$acc end parallel
#else
  DO JSP=1,KSPEC2V
    DO JLEV=1,NFLEVG
      ZSPDIVP(JLEV,JSP)=ZSDIVP(JLEV,JSP)&
       & /(1.0_JPRB-ZBDT2*SIVP(JLEV)*YDLAP%RLAPDI(YDLAP%NVALUE(JSP+IOFF)))  
    ENDDO
  ENDDO
#endif
if (lhook) call dr_hook('SPCSI_boucle2',1,zhook_handle2)

ENDIF

!           Vertical eigenmodes space --> current space.

if (lhook) call dr_hook('SPCSI_mxmaop2',0,zhook_handle2)
#if defined(_OPENACC)
!$acc host_data use_device(SIMO,ZSPDIVP,PSPDIVG)
CALL cublasDgemm('N','T',kspec2v,nflevg,nflevg,1.0_JPRB,&    !!ispcol remplacé par kspec2v
&ZSPDIVP(1,1),kspec2v,SIMO,NFLEVG,0.0_JPRB,PSPDIVG(1,1),kspec2v) !!2ksta remplacés par 1,pourrait partir
!$acc end host_data
!$acc wait
#else
CALL MXMAOP(SIMO,1,NFLEVG,ZSPDIVP,1,NFLEVG,PSPDIVG,1,NFLEVG,NFLEVG,NFLEVG,KSPEC2V)  
#endif
if (lhook) call dr_hook('SPCSI_mxmaop2',1,zhook_handle2)

IF (LSIDG) THEN

  DO JMLOC=NPTRMF(MYSETN), NPTRMF(MYSETN+1)-1
    CALL SPCSIDG_PART2 (YDGEOMETRY, KSPEC2V, JMLOC, PSPDIVG, ZHELP)
  ENDDO

ELSE

  !       ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GMBAR**2 * DIVprim(t+dt)) .
if (lhook) call dr_hook('SPCSI_boucle3',0,zhook_handle2)

#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JSP,JLEV) default(none)
!$acc loop gang
  DO JLEV=1,NFLEVG
    !$acc loop vector
    DO JSP=1,kspec2v
      ZHELP(JSP,jlev)=PSPDIVG(JSP,jlev)*RSTRET*RSTRET
    ENDDO
  ENDDO
!$acc END PARALLEL
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
  DO JSP=1,KSPEC2V
    DO JLEV=1,NFLEVG
      ZHELP(JLEV,JSP)=PSPDIVG(JLEV,JSP)*RSTRET*RSTRET
    ENDDO
  ENDDO
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
call SITNU_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,KSPEC2V,ZHELP,ZST,ZSP,ZSPHI,ZOUT)
#else
CALL SITNU_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,KSPEC2V,ZHELP,ZST,ZSP)
#endif

!*       2.5  Increment Temperature and surface pressure

if (lhook) call dr_hook('SPCSI_boucle4',0,zhook_handle2)
#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JSP,JLEV) default(none)
!$acc loop collapse(2)
DO JLEV=1,nflevg
  DO JSP=1,kspec2v
    PSPTG(JSP,jlev)=PSPTG(JSP,jlev)-ZBDT*ZST(JSP,jlev)
  ENDDO
ENDDO
!$acc end parallel

!$acc parallel loop private(jsp) default(none)
do jsp=1,kspec2v
  pspspg(jsp)=pspspg(jsp)-zbdt*zsp(jsp)
enddo
!$acc end parallel

#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
DO JSP=1,KSPEC2V
  DO JLEV=1,NFLEVG
    PSPTG(JLEV,JSP)=PSPTG(JLEV,JSP)-ZBDT*ZST(JLEV,JSP)
  ENDDO
  PSPSPG(JSP)=PSPSPG(JSP)-ZBDT*ZSP(JSP)
ENDDO
!$OMP END PARALLEL DO
#endif
if (lhook) call dr_hook('SPCSI_boucle4',1,zhook_handle2)

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

