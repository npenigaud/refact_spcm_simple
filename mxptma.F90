#if defined(_OPENACC)
SUBROUTINE MXPTMA(KLX,KVX,KVXS,PA,PBI,PCI,PBS,PCS,PX,PY)
#else
SUBROUTINE MXPTMA(KLX,KVX,KVXS,KIX,tnsmax,PA,PBI,PCI,PBS,PCS,PX,PY)
#endif

!**** *MXPTMA*   - Multiplication of a pentadiagonal matrix by a matrix.

!     Purpose.    Multiplication of a pentadiagonal matrix by a matrix.
!     --------

!    Example, for KLX=5,KVX=2,KIX=1.

!    I PA (1) PBS(1) PCS(1)   0      0    I   I PX(1,1) PX(2,1) I
!    I PBI(1) PA (2) PBS(2) PCS(2)   0    I   I PX(1,2) PX(2,2) I
!    I PCI(1) PBI(2) PA (3) PBS(3) PCS(3) I * I PX(1,3) PX(2,3) I
!    I   0    PCI(2) PBI(3) PA (4) PBS(4) I   I PX(1,4) PX(2,4) I
!    I   0      0    PCI(3) PBI(4) PA (5) I   I PX(1,5) PX(2,5) I

!      I PY(1,1) PY(2,1) I
!      I PY(1,2) PY(2,2) I
!    = I PY(1,3) PY(2,3) I
!      I PY(1,4) PY(2,4) I
!      I PY(1,5) PY(2,5) I

!**   Interface.
!     ----------
!        *CALL* *MXPTMA(KLX,KVX,KVXS,KIX,PA,PBI,PCI,PBS,PCS,PX,PY)

!        Explicit arguments :
!        --------------------
!         KLX:         - Dimension of the matrix.                   (input)
!         KVX,KIX:     - Number of vectors to be multiplied is      (input)
!                        KVX*KIX.
!         KVXS:        - Surdimension corresponding to KVX.         (input)
!         PA:          - Diagonal of the matrix.                    (input)
!         PBI,PCI:     - Lower diagonals of the matrix.             (input)
!         PBS,PCS:     - Upper diagonals of the matrix.             (input)
!         PX:          - Initial vector:                            (input)
!         PY:          - Final vector:                              (output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None.

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
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVXS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX
#if defined(_OPENACC)

#else
INTEGER(KIND=JPIM),INTENT(IN)    :: KIX 
integer(kind=jpim),intent(in) :: tnsmax
#endif 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBI(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBS(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCS(KLX)
#if defined(_OPENACC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX(2*KLX,kvxs) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PY(2*KLX,kvxs) 
#else 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX(tnsmax+1,KVXS,KIX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PY(tnsmax+1,KVXS,KIX) 
#endif

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JI, JL, JV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
!!IF (LHOOK) CALL DR_HOOK('MXPTMA',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF PY.
!              ------------------
  !$acc data present(px,py,pa,pbi,pci,pbs,pcs)

#if defined(_OPENACC)
IF (KLX >= 4) THEN
  !$acc parallel private(ji,jv) default(none) async(1)
  !$acc loop gang
  do jv=1,kvx    
    !$acc loop vector
    do ji=-1,0
      PY(2+ji,jv) = PA (1)*PX(2+ji,jv)+PBS(1)*PX(4+ji,jv)+PCS(1)*PX(6+ji,jv)
      PY(4+ji,jv) = PBI(1)*PX(2+ji,jv)&
       & +PA (2)*PX(4+ji,jv)&
       & +PBS(2)*PX(6+ji,jv)&
       & +PCS(2)*PX(8+ji,jv)  
    enddo
  ENDDO
  !$acc end parallel

  !$acc parallel private(jl,jv) default(none) async(2)
  !$acc loop gang vector collapse(2)
  do jv=1,kvx
    do jl=3,klx-2        
!    !$acc cache(px(jl-2:jl+2,:),pbi(jl-2),pci(jl-1),pa(jl),pbs(jl),pcs(jl))
        PY(2*JL-1:2*jl,jv) = PCI(JL-2)*PX(2*(JL-2)-1:2*(jl-2),jv)&
         & +PBI(JL-1)*PX(2*(JL-1)-1:2*(jl-1),jv)&
         & +PA (JL  )*PX(2*JL-1:2*jl,jv  )&
         & +PBS(JL  )*PX(2*(JL+1)-1:2*(jl+1),jv)&
         & +PCS(JL  )*PX(2*(JL+2)-1:2*(jl+2),jv) 
!        PY(jl,jv) = PCI(JL-2)*PX((jl-2),jv)&
!         & +PBI(JL-1)*PX((jl-1),jv)&
!         & +PA (JL  )*PX(jl,jv  )&
!         & +PBS(JL  )*PX((jl+1),jv)&
!         & +PCS(JL  )*PX((jl+2),jv) 
      !  if (jv==10) then
      !    print *,"PX",PX(2*jl-1,jv),PX(2*jl,jv)
      !    print *,"PY",PY(2*jl-1,jv),PY(2*jl,jv)
      !    print *,"coeffs",PCI(jl-2),PBI(JL-1),PA(jl),PBS(jl),pcs(jl)
      !  endif
     ENDDO
  ENDDO
  !$acc end parallel
  
  !$acc parallel private(ji,jv) async(3)
  !$acc loop gang
  do jv=1,kvx
    !$acc loop vector
    do ji=-1,0
      PY(2*(KLX-1)+ji,jv) = PCI(KLX-3)*PX(2*(KLX-3)+ji,jv)&
       & +PBI(KLX-2)*PX(2*(KLX-2)+ji,jv)&
       & +PA (KLX-1)*PX(2*(KLX-1)+ji,jv)&
       & +PBS(KLX-1)*PX(2*KLX+ji,jv)
      PY(2*KLX+ji,jv) = PCI(KLX-2)*PX(2*(KLX-2)+ji,jv)&
       & +PBI(KLX-1)*PX(2*(KLX-1)+ji,jv)&
       & +PA (KLX  )*PX(2*KLX+ji,jv  )  
    enddo
  ENDDO
  !$acc end parallel
  !$acc wait
ELSEIF (KLX == 3) THEN
  !$acc parallel private(jv)
  !$acc loop vector
  do jv=1,kvx
      PY(1:2,jv) = PA (1)*PX(1:2,jv)+PBS(1)*PX(3:4,jv)+PCS(1)*PX(5:6,jv)
      PY(3:4,jv) = PBI(1)*PX(1:2,jv)+PA (2)*PX(3:4,jv)+PBS(2)*PX(5:6,jv)
      PY(5:6,jv) = PCI(1)*PX(1:2,jv)+PBI(2)*PX(3:4,jv)+PA (3)*PX(5:6,jv)
  ENDDO
  !$acc end parallel
ELSEIF (KLX == 2) THEN
  !$acc parallel private(jv)
  !$acc loop vector
  do jv=1,kvx
      PY(1:2,jv) = PA (1)*PX(1:2,jv)+PBS(1)*PX(3:4,jv)
      PY(3:4,jv) = PBI(1)*PX(1:2,jv)+PA (2)*PX(3:4,jv)
  ENDDO
  !$acc end parallel
ELSEIF (KLX == 1) THEN
  !$acc parallel private(jv)
  !$acc loop vector
  do jv=1,kvx
      PY(1:2,jv) = PA (1)*PX(1:2,jv)
  enddo
  !$acc end parallel
ENDIF
  !$acc end data

#else
IF (KLX >= 4) THEN
  do jv=1,kvx
    DO JI=1,KIX
      PY(1,jv,JI) = PA (1)*PX(1,jv,JI)+PBS(1)*PX(2,jv,JI)+PCS(1)*PX(3,jv,JI)
      PY(2,jv,JI) = PBI(1)*PX(1,jv,JI)&
       & +PA (2)*PX(2,jv,JI)&
       & +PBS(2)*PX(3,jv,JI)&
       & +PCS(2)*PX(4,jv,JI)  
    ENDDO
  ENDDO
  do jl=3,klx-2 
    DO JI=1,KIX
      DO JV=1,kvx
        PY(JL,jv,JI) = PCI(JL-2)*PX(JL-2,jv,JI)&
         & +PBI(JL-1)*PX(JL-1,jv,JI)&
         & +PA (JL  )*PX(JL,jv  ,JI)&
         & +PBS(JL  )*PX(JL+1,jv,JI)&
         & +PCS(JL  )*PX(JL+2,jv,JI)  
      ENDDO
    ENDDO
  ENDDO
  do jv=1,kvx
    DO JI=1,KIX
      PY(KLX-1,jv,JI) = PCI(KLX-3)*PX(KLX-3,jv,JI)&
       & +PBI(KLX-2)*PX(KLX-2,jv,JI)&
       & +PA (KLX-1)*PX(KLX-1,jv,JI)&
       & +PBS(KLX-1)*PX(KLX,jv  ,JI)  
      PY(KLX,jv,JI) = PCI(KLX-2)*PX(KLX-2,jv,JI)&
       & +PBI(KLX-1)*PX(KLX-1,jv,JI)&
       & +PA (KLX  )*PX(KLX,jv  ,JI)  
    ENDDO
  ENDDO

ELSEIF (KLX == 3) THEN
  do jv=1,kvx
    DO JI=1,KIX
      PY(1,jv,JI) = PA (1)*PX(1,jv,JI)+PBS(1)*PX(2,jv,JI)+PCS(1)*PX(3,jv,JI)
      PY(2,jv,JI) = PBI(1)*PX(1,jv,JI)+PA (2)*PX(2,jv,JI)+PBS(2)*PX(3,jv,JI)
      PY(3,jv,JI) = PCI(1)*PX(1,jv,JI)+PBI(2)*PX(2,jv,JI)+PA (3)*PX(3,jv,JI)
    ENDDO
  ENDDO

ELSEIF (KLX == 2) THEN
  do jv=1,kvx
    DO JI=1,KIX
      PY(1,jv,JI) = PA (1)*PX(1,jv,JI)+PBS(1)*PX(2,jv,JI)
      PY(2,jv,JI) = PBI(1)*PX(1,jv,JI)+PA (2)*PX(2,jv,JI)
    ENDDO
  ENDDO

ELSEIF (KLX == 1) THEN
  do jv=1,kvx
    DO JI=1,KIX
      PY(1,jv,JI) = PA (1)*PX(1,jv,JI)
    ENDDO
  ENDDO

ENDIF
#endif
!     ------------------------------------------------------------------

!!IF (LHOOK) CALL DR_HOOK('MXPTMA',1,ZHOOK_HANDLE)
END SUBROUTINE MXPTMA

