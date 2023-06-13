SUBROUTINE SPCSIDG_PART0 (YDGEOMETRY, YDDYN, YDRIP, KSPEC2V, KMLOC, ZSDIV, PSPDIVG)
!$acc routine vector
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN       , ONLY : TDYN
USE YOMRIP       , ONLY : TRIP
USE YOMMP0       , ONLY : MYSETV

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN),value    :: KSPEC2V
INTEGER(KIND=JPIM),INTENT(IN),value    :: KMLOC
REAL(KIND=JPRB),   INTENT(INOUT) :: ZSDIV (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSPDIVG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB)  :: ZBDT

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM):: IN,IM,ISTA,IEND,JSP,JLEV,IOFF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -----------------------

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART0',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX,NFLEVG=>YDDIMV%NFLEVG,SIHEG=>YDDYN%SIHEG,SIHEG2=>YDDYN%SIHEG2,&
& NSPSTAF=>YDMP%NSPSTAF,TDT=>YDRIP%TDT,RBTS2=>YDDYN%RBTS2,NPTRSVF=>YDMP%NPTRSVF)
IM=YDLAP%MYMS(KMLOC)

ISTA=NSPSTAF(IM)
IEND=ISTA+2*(NSMAX+1-IM)-1
ZBDT=RBTS2*TDT
IOFF=NPTRSVF(MYSETV)-1
!$acc data present(zsdiv,ydlap,ydlap%rlapin,pspdivg)
    IF (IM > 0) THEN
#if defined(_OPENACC)
      !!$acc loop vector private(jsp,jlev,in) !!collapse(2) private(jsp,jlev,in)
      !$acc loop vector private(jsp,jlev,in) collapse(2)
#else
      !$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
#endif
      do jlev=1,nflevg
       DO JSP=ISTA,IEND
          IN=YDLAP%NVALUE(JSP+IOFF)
 !!         do jlev=1,nflevg
          ZSDIV(jsp,JLEV)=YDLAP%RLAPIN(IN)*PSPDIVG(jsp,JLEV)-ZBDT*ZSDIV(jsp,JLEV)
       ENDDO
      ENDDO
#if defined(_OPENACC)

#else
     !$OMP END PARALLEL DO
#endif

    ELSE

#if defined(_OPENACC)
      !!$acc loop vector private(jsp,in,jlev) !!collapse(2) private(jsp,in,jlev)
      !$acc loop vector private(jsp,in,jlev) collapse(2)
#else
      !$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
#endif
      do jlev=1,nflevg
        DO JSP=ISTA,IEND
          IN=YDLAP%NVALUE(JSP+IOFF)
  !!        do jlev=1,nflevg
          ZSDIV(jsp,JLEV)=PSPDIVG(jsp,JLEV)-ZBDT*YDLAP%RLAPDI(IN)*ZSDIV(jsp,JLEV)
        ENDDO
      ENDDO
#if defined(_OPENACC)

#else
      !$OMP END PARALLEL DO
#endif

    ENDIF
!$acc end data

END ASSOCIATE
END ASSOCIATE

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART0',1,ZHOOK_HANDLE)
END SUBROUTINE SPCSIDG_PART0

