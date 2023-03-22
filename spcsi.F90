SUBROUTINE SPCSI(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KSPEC2V,LDONEM,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
 ! --- INPUT OPTIONAL --------------------------------------------------------
 & PSPAUXG)

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
!                              LDONEM  - T if only one m if processed
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

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC2V
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
LOGICAL           ,INTENT(IN)    :: LDONEM 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(KSPEC2V) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PSPAUXG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSDIVP (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB) :: ZSPDIVP(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)

REAL(KIND=JPRB) :: ZSDIV  (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB) :: ZHELP  (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB) :: ZST    (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB) :: ZSP    (KSPEC2V)

INTEGER(KIND=JPIM) :: IN, IOFF, JLEV, JSP  
INTEGER(KIND=JPIM) :: JMLOC, IM, ISTA, IEND
REAL(KIND=JPRB) :: ZBDT, ZBDT2

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"
#include "mxptma.h"
#include "mxture.h"
#include "mxturs.h"
#include "abor1.intfb.h"
#include "sigam_sp_openmp.intfb.h"
#include "spcimpfsolve.intfb.h"
#include "sitnu_sp_openmp.intfb.h"
#include "spcsidg_part1.intfb.h"
#include "spcsidg_part2.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPCSI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,   &
& YDLAP=>YDGEOMETRY%YRLAP, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, LIMPF=>YDDYN%LIMPF, LSIDG=>YDDYN%LSIDG, RBTS2=>YDDYN%RBTS2,            &
& SIMI=>YDDYN%SIMI, SIMO=>YDDYN%SIMO, SIVP=>YDDYN%SIVP, RSTRET=>YDGEM%RSTRET, LRSIDDH=>YDLDDH%LRSIDDH,  &
& NPTRSV=>YDMP%NPTRSV, NPTRSVF=>YDMP%NPTRSVF, NSPEC2V=>YDMP%NSPEC2V, NSPEC2VF=>YDMP%NSPEC2VF,           &
& TDT=>YDRIP%TDT,NPTRMF=>YDMP%NPTRMF,NSPSTAF=>YDMP%NSPSTAF,NSMAX=>YDDIM%NSMAX)
!     ------------------------------------------------------------------

!*       0.    TESTINGS.
!              ---------

IF (LIMPF .AND. .NOT.PRESENT(PSPAUXG)) THEN
  CALL ABOR1(' SPCSI: If LIMPF=T, argument PSPAUXG must be present!')
ENDIF

!     ------------------------------------------------------------------

!*       1.    MEMORY TRANSFER.
!              ----------------

IF (LRSIDDH) THEN
  ! DDH memory transfer
  IF (LIMPF) PSPTNDSI_VORG=-PSPVORG
  PSPTNDSI_DIVG=-PSPDIVG
  PSPTNDSI_TG  =-PSPTG  
  !the case of surface pressure has not been treated yet
ENDIF

!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

!*        2.1  Preliminary initialisations.

IF (LDONEM) THEN
  IOFF=NPTRSVF(MYSETV)-1
ELSE
  IOFF=NPTRSV(MYSETV)-1
ENDIF

ZBDT=RBTS2*TDT
ZBDT2=(ZBDT*RSTRET)**2


!*        2.3  Computes right-hand side of Helmholtz equation.

IF( .NOT.LDONEM ) CALL GSTATS(1655,0) ! Main routines and loops in SIGAM chain are parallel
CALL SIGAM_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,KSPEC2V,ZSDIV,PSPTG,PSPSPG)
IF( .NOT.LDONEM ) CALL GSTATS(1655,1)

IF( .NOT.LDONEM ) CALL GSTATS(1656,0)

IF (LSIDG) THEN

  DO JMLOC=NPTRMF(MYSETN), NPTRMF(MYSETN+1)-1
    IM=YDLAP%MYMS(JMLOC)
    ISTA=NSPSTAF(IM)
    IEND=ISTA+2*(NSMAX+1-IM)-1
  
    IF (IM > 0) THEN
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
      DO JSP=ISTA,IEND
        DO JLEV=1,NFLEVG
          IN=YDLAP%NVALUE(JSP+IOFF)
          ZSDIV(JLEV,JSP)=YDLAP%RLAPIN(IN)*PSPDIVG(JLEV,JSP)-ZBDT*ZSDIV(JLEV,JSP)      
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ELSE
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
      DO JSP=ISTA,IEND
        DO JLEV=1,NFLEVG
          IN=YDLAP%NVALUE(JSP+IOFF)
          ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-ZBDT*YDLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)  
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF
  
  ENDDO

ELSE

  ! Case of No Stretching
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
  DO JSP=1,KSPEC2V
    DO JLEV=1,NFLEVG
      IN=YDLAP%NVALUE(JSP+IOFF)
      ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-ZBDT*YDLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!        Add [F] * result to rhs of Helmholtz equation

IF (LIMPF) THEN
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
  DO JSP=1,KSPEC2V
    DO JLEV=1,NFLEVG
      ZSDIV(JLEV,JSP)=ZSDIV(JLEV,JSP) + PSPAUXG(JLEV,JSP)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF
IF( .NOT.LDONEM ) CALL GSTATS(1656,1)

!*        2.4  Solve Helmholtz equation

!           Current space --> vertical eigenmodes space.

IF( .NOT.LDONEM ) CALL GSTATS(1660,0) ! MXMAOP Call to SGEMMX Parallelised
CALL MXMAOP(SIMI,1,NFLEVG,ZSDIV,1,NFLEVG,ZSDIVP,1,NFLEVG,NFLEVG,NFLEVG,KSPEC2V)  
IF( .NOT.LDONEM ) CALL GSTATS(1660,1)

IF (LSIDG) THEN

  DO JMLOC=NPTRMF(MYSETN), NPTRMF(MYSETN+1)-1
    IM=YDLAP%MYMS(JMLOC)
    ISTA=NSPSTAF(IM)
    IEND=ISTA+2*(NSMAX+1-IM)-1
    CALL SPCSIDG_PART1 (YDGEOMETRY, YDDYN, KSPEC2V, IM, JMLOC, ISTA, IEND, ZSDIVP, ZSPDIVP)
  ENDDO

ELSE

  !             Case with NO Stretching :

  IF (LIMPF) THEN

    !               Solve complex pentadiagonal system

    CALL SPCIMPFSOLVE(YDGEOMETRY,YDCST,YDRIP,YDDYN,.FALSE.,.FALSE.,LDONEM,ZSDIVP,ZSPDIVP)

  ELSE

    !                 Inversion of a diagonal system (Helmholtz equation)
    !                 --> (SIMI*DIVprim(t+dt)).

    DO JSP=1,KSPEC2V
      DO JLEV=1,NFLEVG
        ZSPDIVP(JLEV,JSP)=ZSDIVP(JLEV,JSP)&
         & /(1.0_JPRB-ZBDT2*SIVP(JLEV)*YDLAP%RLAPDI(YDLAP%NVALUE(JSP+IOFF)))  
      ENDDO
    ENDDO
  ENDIF
ENDIF

!           Vertical eigenmodes space --> current space.

IF( .NOT.LDONEM ) CALL GSTATS(1660,0) ! MXMAOP Calls SGEMMX in parallel region
CALL MXMAOP(SIMO,1,NFLEVG,ZSPDIVP,1,NFLEVG,PSPDIVG,1,NFLEVG,NFLEVG,NFLEVG,KSPEC2V)  
IF( .NOT.LDONEM ) CALL GSTATS(1660,1)

IF (LSIDG) THEN

  DO JMLOC=NPTRMF(MYSETN), NPTRMF(MYSETN+1)-1
    IM=YDLAP%MYMS(JMLOC)
    ISTA=NSPSTAF(IM)
    IEND=ISTA+2*(NSMAX+1-IM)-1
    CALL SPCSIDG_PART2 (YDGEOMETRY, KSPEC2V, IM, JMLOC, ISTA, IEND, PSPDIVG, ZHELP)
  ENDDO

ELSE

  !       ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GMBAR**2 * DIVprim(t+dt)) .

  IF( .NOT.LDONEM ) CALL GSTATS(1656,0)
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
  DO JSP=1,KSPEC2V
    DO JLEV=1,NFLEVG
      ZHELP(JLEV,JSP)=PSPDIVG(JLEV,JSP)*RSTRET*RSTRET
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  IF( .NOT.LDONEM ) CALL GSTATS(1656,1)

ENDIF

!       If LSIDG:
!         (GM**2 * DIVprim(t+dt)) --> [ tau * (GM**2 * DIVprim(t+dt)) ]
!                                 and [  nu * (GM**2 * DIVprim(t+dt)) ]
!       or if not LSIDG:
!         (GMBAR**2 * DIVprim(t+dt)) --> [ tau * (GMBAR**2 * DIVprim(t+dt)) ]
!                                    and [  nu * (GMBAR**2 * DIVprim(t+dt)) ]

IF( .NOT.LDONEM ) CALL GSTATS(1657,0)  ! Main routines and loops in SITNU chain are parallel
CALL SITNU_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,KSPEC2V,ZHELP,ZST,ZSP)
IF( .NOT.LDONEM ) CALL GSTATS(1657,1)

!*       2.5  Increment Temperature and surface pressure

IF( .NOT.LDONEM ) CALL GSTATS(1656,0)
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
DO JSP=1,KSPEC2V
  DO JLEV=1,NFLEVG
    PSPTG(JLEV,JSP)=PSPTG(JLEV,JSP)-ZBDT*ZST(JLEV,JSP)
  ENDDO
  PSPSPG(JSP)=PSPSPG(JSP)-ZBDT*ZSP(JSP)
ENDDO
!$OMP END PARALLEL DO
IF( .NOT.LDONEM ) CALL GSTATS(1656,1)

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF SI TERM AT t+dt FOR DDH.
!              ---------------------------------------

IF (LRSIDDH) THEN
  IF (LIMPF) PSPTNDSI_VORG=PSPTNDSI_VORG + PSPVORG
  PSPTNDSI_DIVG=PSPTNDSI_DIVG + PSPDIVG
  PSPTNDSI_TG=PSPTNDSI_TG + PSPTG
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPCSI',1,ZHOOK_HANDLE)
END SUBROUTINE SPCSI

