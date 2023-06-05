SUBROUTINE SPCSIDG_PART0 (YDGEOMETRY, YDDYN, YDRIP, KSPEC2V, KMLOC, ZSDIV, PSPDIVG)
!!$acc routine gang
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
#if defined(_OPENACC)
REAL(KIND=JPRB),   INTENT(INOUT) :: ZSDIV (kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSPDIVG(kspec2v,YDGEOMETRY%YRDIMV%NFLEVG)
#else
REAL(KIND=JPRB),   INTENT(INOUT) :: ZSDIV (YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
REAL(KIND=JPRB),   INTENT(INOUT) :: PSPDIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSPEC2V)
#endif

REAL(KIND=JPRB)  :: ZBDT
real(kind=jprb)  :: rlapinin,rlapdiin
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM):: IN,IM,ISTA,IEND,JSP,JLEV,IOFF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -----------------------

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART0',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP,YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX,NFLEVG=>YDDIMV%NFLEVG,SIHEG=>YDDYN%SIHEG,SIHEG2=>YDDYN%SIHEG2,&
& NSPSTAF=>YDMP%NSPSTAF,TDT=>YDRIP%TDT,RBTS2=>YDDYN%RBTS2,NPTRSVF=>YDMP%NPTRSVF)
!$acc data present(zsdiv,ydlap,ydlap%rlapin,ydlap%rlapdi,pspdivg,ydlap%myms,ydmp,ydmp%nspstaf,ydmp%nptrsvf)
!$acc data present(nflevg,nsmax,rbts2,tdt,mysetv,ydgeometry%yrmp%nptrsvf,nptrsvf,nspstaf)
!$acc parallel num_gangs(16) num_workers(1) vector_length(64) default(none)
IM=YDLAP%MYMS(KMLOC)

ISTA=NSPSTAF(IM)
IEND=ISTA+2*(NSMAX+1-IM)-1
ZBDT=RBTS2*TDT
IOFF=NPTRSVF(MYSETV)-1
    IF (IM > 0) THEN
#if defined(_OPENACC)
!$acc loop gang private(jsp,jlev,in,rlapinin)
do jlev=1,nflevg
!$acc loop vector !!PRIVATE(JSP,IN,rlapinin,jlev)
      DO JSP=ISTA,IEND
        !!do jlev=1,nflevg
          IN=YDLAP%NVALUE(JSP+IOFF)
          rlapinin=ydlap%rlapin(in)
          ZSDIV(jsp,JLEV)=rlapinin*PSPDIVG(jsp,JLEV)-ZBDT*ZSDIV(jsp,JLEV)!!rlapinin*pspdivg      
       ENDDO
      ENDDO

#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
      DO JSP=ISTA,IEND
        DO JLEV=1,NFLEVG
          IN=YDLAP%NVALUE(JSP+IOFF)
          ZSDIV(JLEV,JSP)=YDLAP%RLAPIN(IN)*PSPDIVG(JLEV,JSP)-ZBDT*ZSDIV(JLEV,JSP)      
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
#endif
    ELSE
#if defined(_OPENACC)
!$acc loop gang private(jsp,in,rlapdiin,jlev)
do jlev=1,nflevg
!$acc loop vector !!PRIVATE(JSP,IN,rlapdiin,jlev)
      DO JSP=ISTA,IEND
        !!do jlev=1,nflevg
          IN=YDLAP%NVALUE(JSP+IOFF)
          rlapdiin=ydlap%rlapdi(in)*zbdt
          ZSDIV(jsp,JLEV)=PSPDIVG(jsp,JLEV)-ZSDIV(jsp,JLEV)*rlapdiin  
        ENDDO
      ENDDO

#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
      DO JSP=ISTA,IEND
        DO JLEV=1,NFLEVG
          IN=YDLAP%NVALUE(JSP+IOFF)
          ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-ZBDT*YDLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)  
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
#endif
    ENDIF
!$acc end parallel
!$acc end data 
!$acc end data

END ASSOCIATE
END ASSOCIATE

!!IF (LHOOK) CALL DR_HOOK('SPCSIDG_PART0',1,ZHOOK_HANDLE)
END SUBROUTINE SPCSIDG_PART0

