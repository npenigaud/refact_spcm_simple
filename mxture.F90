#if defined(_OPENACC)
SUBROUTINE MXTURE(KLX,KVX,KVXS,KIX,KIXS,tnsmax,KT,LDMT,PA,PB,PC,PY,PX,tbloc,pas,pbs,pcs,pys,pxs)
!$acc routine vector
#else
SUBROUTINE MXTURE(KLX,KVX,KVXS,KIX,tnsmax,KT,LDMT,PA,PB,PC,PY,PX)
#endif

!**** *MXTURE*   - Resolution of a set of triangular tridiagonal systems.

!     Purpose.    Resolution of a set of triangular tridiagonal systems.
!     --------

!    Example, for KLX=4:

!    PY is known, PX is unknown. We want to find PX, solution of the
!     following system, for l=1 to KVX, i=1 to KIX.

!    * KT = -2:

!    I PA(l,1)    0       0       0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I PB(l,1) PA(l,2)    0       0    I   I PX(l,2,i) I   I PY(l,2,i) I
!    I PC(l,1) PB(l,2) PA(l,3)    0    I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0    PC(l,2) PB(l,3) PA(l,4) I   I PX(l,4,i) I   I PY(l,4,i) I

!    * KT = -3:

!    I    1       0       0       0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I PB(l,1)    1       0       0    I   I PX(l,2,i) I   I PY(l,2,i) I
!    I PC(l,1) PB(l,2)    1       0    I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0    PC(l,2) PB(l,3)    1    I   I PX(l,4,i) I   I PY(l,4,i) I

!    Dummy array PA is not used in this case.

!    * KT = 1:

!    I    1    ZB(l,1) ZC(l,1)    0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I    0       1    ZB(l,2) ZC(l,2) I   I PX(l,2,i) I   I PY(l,2,i) I
!    I    0       0       1    ZB(l,3) I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0       0       0       1    I   I PX(l,4,i) I   I PY(l,4,i) I

!    where ZB(l,j)=PB(l,j)/PA(l,j); ZC(l,j)=PC(l,j)/PA(l,j) .

!    * KT = 2:

!    I PA(l,1) PB(l,1) PC(l,1)    0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I    0    PA(l,2) PB(l,2) PC(l,2) I   I PX(l,2,i) I   I PY(l,2,i) I
!    I    0       0    PA(l,3) PB(l,3) I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0       0       0    PA(l,4) I   I PX(l,4,i) I   I PY(l,4,i) I

!    * KT = 3:

!    I    1    PB(l,1) PC(l,1)    0    I   I PX(l,1,i) I   I PY(l,1,i) I
!    I    0       1    PB(l,2) PC(l,2) I   I PX(l,2,i) I   I PY(l,2,i) I
!    I    0       0       1    PB(l,3) I * I PX(l,3,i) I = I PY(l,3,i) I
!    I    0       0       0       1    I   I PX(l,4,i) I   I PY(l,4,i) I

!    Dummy array PA is not used in this case.

!**   Interface.
!     ----------
!        *CALL* *MXTURE(KLX,KVX,KVXS,KIX,KT,LDMT,PA,PB,PC,PY,PX)

!        Explicit arguments :
!        --------------------
!         KLX:        - Dimension of the system.                    (input)
!         KVX:        - Number of variables (second                 (input)
!                       dimension of PX and PY).
!         KVXS:       - Surdimension corresponding to KVX.          (input)
!         KIX:        - Number of variables multiplied by the same
!                       matrix.                                     (input)
!         KT:         - Type of matrix (see figures above).         (input)
!         LDMT:       - .T.: copy of PX on PY at the end of subr.   (input)
!                       .F.: no copy.
!         PA,PB,PC:   - non-zero diagonals of the triangular        (input)
!                       matrixes (see figures above).
!                       Caution! All PA coefficients must be non 0.
!         PY:         - known vector.                               (input)
!         PX:         - unknown vector.                             (output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None.
!        Called by SPC,SPCAD,SPC2,SPC2AD,MXTURS.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD: OCTOBER 1993.

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!     ------------------------------------------------------------------

IMPLICIT NONE

#if defined(_OPENACC)
INTEGER(KIND=JPIM),INTENT(IN),value :: KLX
INTEGER(KIND=JPIM),INTENT(IN),value :: KVX
INTEGER(KIND=JPIM),INTENT(IN),value :: KVXS
INTEGER(KIND=JPIM),INTENT(IN),value :: KIX
INTEGER(KIND=JPIM),INTENT(IN),value :: KIXS
INTEGER(KIND=JPIM),INTENT(IN),value :: tnsmax
INTEGER(KIND=JPIM),INTENT(IN),value :: KT
LOGICAL ,INTENT(IN),value :: LDMT
REAL(KIND=JPRB) ,INTENT(IN) :: PA(klx,kvx)
REAL(KIND=JPRB) ,INTENT(IN) :: PB(klx,kvx)
REAL(KIND=JPRB) ,INTENT(IN) :: PC(klx,kvx)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PY(tnsmax+1,kixs,kvx)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PX(tnsmax+1,kixs,kvx)
integer(kind=jpim), intent(in),value :: tbloc
REAL(KIND=JPRB) ,INTENT(INout) :: PAs(tbloc+3,kvxs)
REAL(KIND=JPRB) ,INTENT(INout) :: PBs(tbloc+3,kvxs)
REAL(KIND=JPRB) ,INTENT(INout) :: PCs(tbloc+3,kvxs)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PYs(tbloc+3,kvxs,kixs)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PXs(tbloc+3,kvxs,kixs)
#else
INTEGER(KIND=JPIM),INTENT(IN) :: KLX
INTEGER(KIND=JPIM),INTENT(IN) :: KVX
INTEGER(KIND=JPIM),INTENT(IN) :: KVXS
INTEGER(KIND=JPIM),INTENT(IN) :: KIX
integer(kind=jpim),intent(in) :: tnsmax
INTEGER(KIND=JPIM),INTENT(IN) :: KT
LOGICAL ,INTENT(IN) :: LDMT
REAL(KIND=JPRB) ,INTENT(IN) :: PA(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PB(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(IN) :: PC(KVX,KLX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PY(KVXS,tnsmax+1,KIX)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PX(KVXS,tnsmax+1,KIX)
#endif


!     ------------------------------------------------------------------

#if defined(_OPENACC)
INTEGER(KIND=JPIM) :: IT, JL,jlb,reste,decalage,jv,ji
#else
INTEGER(KIND=JPIM) :: IIX, IT, JI, JL, JV
REAL(KIND=JPRB) :: ZBB, ZCC
#endif

!     ------------------------------------------------------------------
!!IF (LHOOK) CALL DR_HOOK('MXTURE',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

!$acc data present(pa,pb,pc,py,px)

#if defined(_OPENACC)
IT=MIN(3,MAX(-3,KT))
IF (KT == 0.OR.KT == -1) IT=-2
#else
IIX=KIX-MOD(KIX,2)
IT=MIN(3,MAX(-3,KT))
IF (KT == 0.OR.KT == -1) IT=-2
#endif

!      ----------------------------------------------------------------

!*       2.    COMPUTATION OF PX.
!              ------------------

#if defined(_OPENACC)
IF (IT == -3) THEN



ELSEIF (IT == -2) THEN
!bloclev => kvxs min(nflevg-(compteurb-1)*bloclev => kvx 
!pys pys pxs pxs 
!param_mxture => PA PB PC à passer avec le bon décalage 
!decalage1 => COMPTEURC*kvx si param_mxture passé avec decalage de compteurb
!COMPTEURC => jv
        IF (KLX >= 3) THEN
          !$acc loop seq
          do jlb=1,(klx-3)/tbloc+1  !!klx-3+1 elements, de 3 a klx
            decalage=(jlb-1)*tbloc
            !$acc loop vector private(jv,ji)
            do jl=decalage+1,min(decalage+tbloc+2,klx)
              do jv=1,kvx
                pas(jl-decalage,jv)=pa(jl,jv)
                pbs(jl-decalage,jv)=pb(jl,jv)
                pcs(jl-decalage,jv)=pc(jl,jv)
                do ji=1,kix
                  pys(jl-decalage,jv,ji)=PY(jl,ji,jv)
                enddo
              enddo
            enddo
            if (jlb==1) then
              !$acc loop vector collapse(2)
              do jv=1,kvx
                do ji=1,kix
                  pxs(1,jv,ji)=pys(1,jv,ji)/pas(1,jv)
                  pxs(2,jv,ji)=(pys(2,jv,ji)-pbs(1,jv)*pxs(1,jv,ji))/pas(2,jv)
                enddo
              enddo
            else
              !$acc loop vector collapse(2)
              do jv=1,kvx
                do ji=1,kix
                  pxs(1,jv,ji)=PX(decalage+1,ji,jv)
                  pxs(2,jv,ji)=PX(decalage+2,ji,jv)
                enddo
              enddo
            endif
            !$acc loop vector private(jl) collapse(2)
            do jv=1,kvx
              do ji=1,kix
                do jl=3,min(decalage+tbloc+2,klx)-decalage
                  pxs(JL,jv,ji)=(pys(JL,jv,ji)-pcs(JL-2,jv)*pxs(JL-2,jv,ji)&
                   & -pbs(JL-1,jv)*pxs(JL-1,jv,ji))/pas(JL,jv)  
                enddo
              enddo
            enddo
            !$acc loop vector private(jv,ji)
            do jl=decalage+1,min(decalage+tbloc+2,klx)
              do jv=1,kvx
                do ji=1,kix
                  PX(jl,ji,jv)=pxs(jl-decalage,jv,ji)
                enddo
              enddo
            enddo
          enddo
        ELSE
          !$acc loop vector collapse(2)
          do jv=1,kvx
            do ji=1,kix
              PX(1,ji,jv)=PY(1,ji,jv)/pa(1,jv) 
              IF (KLX >= 2) THEN
                PX(2,ji,jv)=(PY(2,ji,jv)&
                 &-pb(1,jv)*PX(1,ji,jv))/pa(2,jv)
              ENDIF
            enddo
          enddo
        ENDIF

ELSEIF (IT == 1) THEN

       IF (KLX >= 3) THEN
         reste=mod(klx-3,tbloc)+1
         !$acc loop seq
         do jlb=(klx-3)/tbloc+1,1,-1  !!klx-3+1 elements, de klx-2 a 1
           decalage=(jlb-2)*tbloc+reste
           !$acc loop vector private(jv,ji)
           do jl=decalage+tbloc+2,max(decalage+1,1),-1
             do jv=1,kvx
               pas(jl-decalage,jv)=pa(jl,jv)
               pbs(jl-decalage,jv)=pb(jl,jv)
               pcs(jl-decalage,jv)=pc(jl,jv)
               do ji=1,kix
                 pys(jl-decalage,jv,ji)=PY(jl,ji,jv)
               enddo
             enddo
           enddo
           if (jlb==(klx-3)/tbloc+1) then
              !$acc loop vector collapse(2)
              do jv=1,kvx
                do ji=1,kix
                  pxs(klx-decalage,jv,ji)=pys(klx-decalage,jv,ji)
                  pxs(klx-1-decalage,jv,ji)=pys(klx-1-decalage,jv,ji)&
                   &-pbs(klx-1-decalage,jv)/pas(klx-1-decalage,jv)*pxs(klx-decalage,jv,ji)
                enddo
              enddo
           else
              !$acc loop vector collapse(2)
              do jv=1,kvx
                do ji=1,kix
                  pxs(tbloc+1,jv,ji)=PX(decalage+tbloc+1,ji,jv)
                  pxs(tbloc+2,jv,ji)=PX(decalage+tbloc+2,ji,jv)
                enddo
              enddo
           endif
           !$acc loop vector private(jl) collapse(2)
           do jv=1,kvx
             do ji=1,kix
               !$acc loop seq
               do jl=tbloc,max(decalage+1,1)-decalage,-1
                 pxs(JL,jv,ji)=pys(JL,jv,ji)-pbs(jl,jv)/pas(jl,jv)*pxs(JL+1,jv,ji)&
                   &-pcs(jl,jv)/pas(jl,jv)*pxs(JL+2,jv,ji)     
               enddo
             enddo
           enddo
           !$acc loop vector private(jv,ji)
           do jl=decalage+tbloc+2,max(decalage+1,1),-1
             do jv=1,kvx
               do ji=1,kix
                 PX(jl,ji,jv)=pxs(jl-decalage,jv,ji)
               enddo
             enddo
           enddo
         enddo
       ELSE
         !$acc loop vector collapse(2)
         do jv=1,kvx
           do ji=1,kix
             PX(KLX,ji,jv)=PY(KLX,ji,jv)
             IF (KLX >= 2) THEN
               PX(KLX-1,ji,jv)=PY(KLX-1,ji,jv)&
                &-pb(klx-1,jv)/pa(klx-1,jv)*PX(KLX,ji,jv)
             ENDIF
           enddo
         enddo
       ENDIF


ELSEIF (IT == 2) THEN
 
!  IF (KLX >= 3) THEN
!    reste=mod(klx-3,tbloc)+1
!    !$acc loop seq
!    do jlb=(klx-3)/tbloc+1,1,-1  !!klx-3+1 elements, de klx-2 a 1
!      decalage=(jlb-2)*tbloc+reste
!      !$acc loop vector
!      do jl=decalage+tbloc+2,max(decalage+1,1),-1
!        pas(jl-decalage)=pa(jl)
!        pbs(jl-decalage)=pb(jl)
!        pcs(jl-decalage)=pc(jl)
!        pys(jl-decalage)=py(jl)
!      enddo
!      if (jlb==(klx-3)/tbloc+1) then
!        pxs(klx-decalage)=pys(klx-decalage)
!        pxs(klx-1-decalage)=(pys(klx-1-decalage)-pbs(klx-1-decalage)*pxs(klx-decalage))/pas(klx-1-decalage)
!      else
!        pxs(tbloc+1:tbloc+2)=px(decalage+tbloc+1:decalage+tbloc+2)
!      endif
!      !$acc loop seq
!      do jl=tbloc,max(decalage+1,1)-decalage,-1
!        PXs(JL)=(PYs(JL)-PBs(JL)*PXs(JL+1)&
!         & -PCs(JL)*PXs(JL+2))/PAs(JL)  
!      enddo
!      !$acc loop vector
!      do jl=decalage+tbloc+2,max(decalage+1,1),-1
!        px(jl)=pxs(jl-decalage)
!      enddo
!    enddo
!  ELSE
!    PX(KLX)=PY(KLX)/PA(KLX)
!    IF (KLX >= 2) THEN
!      PX(KLX-1)=&
!       & (PY(KLX-1)-PB(KLX-1)*PX(KLX))/PA(KLX-1)  
!    ENDIF
!  ENDIF

ELSEIF (IT == 3) THEN

     IF (KLX >= 3) THEN
       reste=mod(klx-3,tbloc)+1
       !$acc loop seq
       do jlb=(klx-3)/tbloc+1,1,-1  !!klx-3+1 elements, de klx-2 a 1
         decalage=(jlb-2)*tbloc+reste
         !$acc loop vector private(jv,ji)
         do jl=decalage+tbloc+2,max(decalage+1,1),-1
           do jv=1,kvx
             pbs(jl-decalage,jv)=pb(jl,jv)
             pcs(jl-decalage,jv)=pc(jl,jv)
             do ji=1,kix
               pys(jl-decalage,jv,ji)=PY(jl,ji,jv)
             enddo
           enddo
         enddo
         if (jlb==(klx-3)/tbloc+1) then
           !$acc loop vector collapse(2)
           do jv=1,kvx
             do ji=1,kix
               pxs(klx-decalage,jv,ji)=pys(klx-decalage,jv,ji)
               pxs(klx-1-decalage,jv,ji)=pys(klx-1-decalage,jv,ji)&
                 &-pbs(klx-1-decalage,jv)*pxs(klx-decalage,jv,ji)
             enddo
           enddo
         else
           !$acc loop vector collapse(2)
           do jv=1,kvx
             do ji=1,kix
               pxs(tbloc+1,jv,ji)=PX(decalage+tbloc+1,ji,jv)
               pxs(tbloc+2,jv,ji)=PX(decalage+tbloc+2,ji,jv)
             enddo
           enddo
         endif
         !$acc loop vector private(jl) collapse(2)
         do jv=1,kvx
           do ji=1,kix
             !$acc loop seq
             do jl=tbloc,max(decalage+1,1)-decalage,-1
               pxs(JL,jv,ji)=pys(JL,jv,ji)-pbs(jl,jv)*pxs(JL+1,jv,ji)&
                 &-pcs(jl,jv)*pxs(JL+2,jv,ji)
             enddo
           enddo
         enddo
         !$acc loop vector private(jv,ji)
         do jl=decalage+tbloc+2,max(decalage+1,1),-1
           do jv=1,kvx
             do ji=1,kix
               PX(jl,ji,jv)=pxs(jl-decalage,jv,ji)
             enddo
           enddo
         enddo
       enddo
     ELSE
       !$acc loop vector collapse(2)
       do jv=1,kvx
         do ji=1,kix
           PX(KLX,ji,jv)=PY(KLX,ji,jv)
           IF (KLX >= 2) THEN
             PX(KLX-1,ji,jv)=PY(KLX-1,ji,jv)&
               &-pb(klx-1,jv)*PX(KLX,ji,jv)
           ENDIF
         enddo
       enddo
     ENDIF


ENDIF

#else

IF (IT == -3) THEN

  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,1,JI)=PY(JV,1,JI)
      PX(JV,1,JI+1)=PY(JV,1,JI+1)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,1,JI)=PY(JV,1,JI)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        ZBB=PB(JV,1)
        PX(JV,2,JI)=PY(JV,2,JI)-ZBB*PX(JV,1,JI)
        PX(JV,2,JI+1)=PY(JV,2,JI+1)-ZBB*PX(JV,1,JI+1)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        ZBB=PB(JV,1)
        PX(JV,2,JI)=PY(JV,2,JI)-ZBB*PX(JV,1,JI)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=3,KLX
        DO JV=1,KVX
          ZBB=PB(JV,JL-1)
          ZCC=PC(JV,JL-2)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL-1,JI)-ZCC*PX(JV,JL-2,JI)
          PX(JV,JL,JI+1)=PY(JV,JL,JI+1)&
           & -ZBB*PX(JV,JL-1,JI+1)-ZCC*PX(JV,JL-2,JI+1)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=3,KLX
        DO JV=1,KVX
          ZBB=PB(JV,JL-1)
          ZCC=PC(JV,JL-2)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL-1,JI)-ZCC*PX(JV,JL-2,JI)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT == -2) THEN

  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,1,JI)=PY(JV,1,JI)/PA(JV,1)
      PX(JV,1,JI+1)=PY(JV,1,JI+1)/PA(JV,1)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,1,JI)=PY(JV,1,JI)/PA(JV,1)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        PX(JV,2,JI)=(PY(JV,2,JI)-PB(JV,1)*PX(JV,1,JI))/PA(JV,2)
        PX(JV,2,JI+1)=(PY(JV,2,JI+1)-PB(JV,1)*PX(JV,1,JI+1))/PA(JV,2)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        PX(JV,2,JI)=(PY(JV,2,JI)-PB(JV,1)*PX(JV,1,JI))/PA(JV,2)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=3,KLX
        DO JV=1,KVX
          PX(JV,JL,JI)=(PY(JV,JL,JI)-PC(JV,JL-2)*PX(JV,JL-2,JI)&
           & -PB(JV,JL-1)*PX(JV,JL-1,JI))/PA(JV,JL)  
          PX(JV,JL,JI+1)=(PY(JV,JL,JI+1)-PC(JV,JL-2)*PX(JV,JL-2,&
           & JI+1)&
           & -PB(JV,JL-1)*PX(JV,JL-1,JI+1))/PA(JV,JL)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=3,KLX
        DO JV=1,KVX
          PX(JV,JL,JI)=(PY(JV,JL,JI)-PC(JV,JL-2)*PX(JV,JL-2,JI)&
           & -PB(JV,JL-1)*PX(JV,JL-1,JI))/PA(JV,JL)  
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT == 1) THEN

  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)
      PX(JV,KLX,JI+1)=PY(JV,KLX,JI+1)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        ZBB=PB(JV,KLX-1)/PA(JV,KLX-1)
        PX(JV,KLX-1,JI)=PY(JV,KLX-1,JI)-ZBB*PX(JV,KLX,JI)
        PX(JV,KLX-1,JI+1)=PY(JV,KLX-1,JI+1)-ZBB*PX(JV,KLX,JI+1)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        ZBB=PB(JV,KLX-1)/PA(JV,KLX-1)
        PX(JV,KLX-1,JI)=PY(JV,KLX-1,JI)-ZBB*PX(JV,KLX,JI)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          ZBB=PB(JV,JL)/PA(JV,JL)
          ZCC=PC(JV,JL)/PA(JV,JL)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL+1,JI)-ZCC*PX(JV,JL+2,JI)
          PX(JV,JL,JI+1)=PY(JV,JL,JI+1)&
           & -ZBB*PX(JV,JL+1,JI+1)-ZCC*PX(JV,JL+2,JI+1)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          ZBB=PB(JV,JL)/PA(JV,JL)
          ZCC=PC(JV,JL)/PA(JV,JL)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL+1,JI)-ZCC*PX(JV,JL+2,JI)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT == 2) THEN
  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)/PA(JV,KLX)
      PX(JV,KLX,JI+1)=PY(JV,KLX,JI+1)/PA(JV,KLX)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)/PA(JV,KLX)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        PX(JV,KLX-1,JI)=&
         & (PY(JV,KLX-1,JI)-PB(JV,KLX-1)*PX(JV,KLX,JI))/PA(JV,KLX-1)  
        PX(KLX-1,jv,JI+1)=&
         & (PY(JV,KLX-1,JI+1)-PB(JV,KLX-1)*PX(JV,KLX,JI+1))/PA(JV,&
         & KLX-1)  
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        PX(JV,KLX-1,JI)=&
         & (PY(JV,KLX-1,JI)-PB(JV,KLX-1)*PX(JV,KLX,JI))/PA(JV,KLX-1)  
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          PX(JV,JL,JI)=(PY(JV,JL,JI)-PB(JV,JL)*PX(JV,JL+1,JI)&
           & -PC(JV,JL)*PX(JV,JL+2,JI))/PA(JV,JL)  
          PX(JV,JL,JI+1)=(PY(JV,JL,JI+1)-PB(JV,JL)*PX(JV,JL+1,JI+&
           & 1)&
           & -PC(JV,JL)*PX(JV,JL+2,JI+1))/PA(JV,JL)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          PX(JV,JL,JI)=(PY(JV,JL,JI)-PB(JV,JL)*PX(JV,JL+1,JI)&
           & -PC(JV,JL)*PX(JV,JL+2,JI))/PA(JV,JL)  
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ELSEIF (IT == 3) THEN

  DO JI=1,IIX,2
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)
      PX(JV,KLX,JI+1)=PY(JV,KLX,JI+1)
    ENDDO
  ENDDO
  DO JI=IIX+1,KIX
    DO JV=1,KVX
      PX(JV,KLX,JI)=PY(JV,KLX,JI)
    ENDDO
  ENDDO

  IF (KLX >= 2) THEN
    DO JI=1,IIX,2
      DO JV=1,KVX
        ZBB=PB(JV,KLX-1)
        PX(JV,KLX-1,JI)=PY(JV,KLX-1,JI)-ZBB*PX(JV,KLX,JI)
        PX(JV,KLX-1,JI+1)=PY(JV,KLX-1,JI+1)-ZBB*PX(JV,KLX,JI+1)
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JV=1,KVX
        ZBB=PB(JV,KLX-1)
        PX(JV,KLX-1,JI)=PY(JV,KLX-1,JI)-ZBB*PX(JV,KLX,JI)
      ENDDO
    ENDDO
  ENDIF

  IF (KLX >= 3) THEN
    DO JI=1,IIX,2
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          ZBB=PB(JV,JL)
          ZCC=PC(JV,JL)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL+1,JI)-ZCC*PX(JV,JL+2,JI)
          PX(JV,JL,JI+1)=PY(JV,JL,JI+1)&
           & -ZBB*PX(JV,JL+1,JI+1)-ZCC*PX(JV,JL+2,JI+1)  
        ENDDO
      ENDDO
    ENDDO
    DO JI=IIX+1,KIX
      DO JL=KLX-2,1,-1
        DO JV=1,KVX
          ZBB=PB(JV,JL)
          ZCC=PC(JV,JL)
          PX(JV,JL,JI)=PY(JV,JL,JI)-ZBB*PX(JV,JL+1,JI)-ZCC*PX(JV,JL+2,JI)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

ENDIF
#endif
!      ----------------------------------------------------------------

!*       3.    FINAL MEMORY TRANSFER.
!              ----------------------

IF (LDMT) THEN
#if defined(_OPENACC)

 !$acc loop vector private(ji,jv)
  do jl=1,klx
    do ji=1,kix
      do jv=1,kvx    
        PY(JL,ji,jv)=PX(JL,ji,jv)
      enddo
    enddo 
  ENDDO

#else

  DO JI=1,KIX
    DO JL=1,KLX
      DO JV=1,KVX
        PY(JV,JL,JI)=PX(JV,JL,JI)
      ENDDO
    ENDDO
  ENDDO

#endif

ENDIF

!$acc end data

!     ------------------------------------------------------------------

!!IF (LHOOK) CALL DR_HOOK('MXTURE',1,ZHOOK_HANDLE)
END SUBROUTINE MXTURE

