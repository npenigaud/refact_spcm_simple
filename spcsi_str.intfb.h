INTERFACE
#if defined(_OPENACC)
SUBROUTINE SPCSI_STR(&
& YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KSPEC2V,&
& PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
& PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
& taillec,zsdivpl,zspdivpl,pa,pb,pc,entree,sortie,param_mxture)

#else
SUBROUTINE SPCSI_STR(&
& YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KSPEC2V,&
& PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
& PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
& simit,simot)
#endif
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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(KSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
#if defined(_OPENACC)
integer(kind=jpim)               :: taillec
real(kind=JPRB)   ,intent(inout) :: zsdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2,64)
real(kind=JPRB)   ,intent(inout) :: zspdivpl(ydgeometry%yrdim%nsmax+1,ydgeometry%yrdimv%nflevg,2,64)
real(kind=jprb)   ,intent(in) :: param_mxture(:,:,:)
real(kind=jprb)   ,intent(inout)    :: pa(taillec)
real(kind=jprb)   ,intent(inout)    :: pb(taillec)
real(kind=jprb)   ,intent(inout)    :: pc(taillec)
real(kind=jprb)   ,intent(inout)    :: entree(taillec)
real(kind=jprb)   ,intent(inout)    :: sortie(taillec)
#else
real(kind=jprb)   ,intent(in)    :: simit(ydgeometry%yrdimv%nflevg,ydgeometry%yrdimv%nflevg)
real(kind=jprb)   ,intent(in)    :: simot(ydgeometry%yrdimv%nflevg,ydgeometry%yrdimv%nflevg)
#endif
END SUBROUTINE SPCSI_STR

END INTERFACE
