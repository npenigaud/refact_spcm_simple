#if defined(_OPENACC)
SUBROUTINE MXPTMA(KLX,KVX,KVXS,KIX,KTNSMAX,PA,PBI,PCI,PBS,PCS,PX,PY)
!$ACC ROUTINE VECTOR
#else
SUBROUTINE MXPTMA(KLX,KVX,KVXS,KIX,KTNSMAX,PA,PBI,PCI,PBS,PCS,PX,PY)
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

!     ------------------------------------------------------------------

IMPLICIT NONE

#if defined(_OPENACC)
INTEGER(KIND=JPIM),INTENT(IN),VALUE    :: KLX 
INTEGER(KIND=JPIM),INTENT(IN),VALUE    :: KVXS 
INTEGER(KIND=JPIM),INTENT(IN),VALUE    :: KIX 
INTEGER(KIND=JPIM),INTENT(IN),VALUE    :: KVX
INTEGER(KIND=JPIM),INTENT(IN),VALUE    :: KTNSMAX
#else
INTEGER(KIND=JPIM),INTENT(IN)    :: KLX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVXS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX
INTEGER(KIND=JPIM),INTENT(IN)    :: KTNSMAX
#endif 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBI(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBS(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCS(KLX)
#if defined(_OPENACC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX(KTNSMAX+1,KVXS,KIX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PY(KTNSMAX+1,KVXS,KIX) 
#else 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX(KTNSMAX+1,KVXS,KIX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PY(KTNSMAX+1,KVXS,KIX) 
#endif

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JI, JL, JV


!*       1.    COMPUTATION OF PY.
!              ------------------

IF (KLX >= 4) THEN
  DO JV=1,KVX
    DO JI=1,KIX
      PY(1,JV,JI) = PA (1)*PX(1,JV,JI)+PBS(1)*PX(2,JV,JI)+PCS(1)*PX(3,JV,JI)
      PY(2,JV,JI) = PBI(1)*PX(1,JV,JI)&
       & +PA (2)*PX(2,JV,JI)&
       & +PBS(2)*PX(3,JV,JI)&
       & +PCS(2)*PX(4,JV,JI)  
    ENDDO
  ENDDO
  !$ACC LOOP VECTOR 
  DO JL=3,KLX-2 
    DO JI=1,KIX
      DO JV=1,KVX
        PY(JL,JV,JI) = PCI(JL-2)*PX(JL-2,JV,JI)&
         & +PBI(JL-1)*PX(JL-1,JV,JI)&
         & +PA (JL  )*PX(JL,JV  ,JI)&
         & +PBS(JL  )*PX(JL+1,JV,JI)&
         & +PCS(JL  )*PX(JL+2,JV,JI)  
      ENDDO
    ENDDO
  ENDDO
  DO JV=1,KVX
    DO JI=1,KIX
      PY(KLX-1,JV,JI) = PCI(KLX-3)*PX(KLX-3,JV,JI)&
       & +PBI(KLX-2)*PX(KLX-2,JV,JI)&
       & +PA (KLX-1)*PX(KLX-1,JV,JI)&
       & +PBS(KLX-1)*PX(KLX,JV  ,JI)  
      PY(KLX,JV,JI) = PCI(KLX-2)*PX(KLX-2,JV,JI)&
       & +PBI(KLX-1)*PX(KLX-1,JV,JI)&
       & +PA (KLX  )*PX(KLX,JV  ,JI)  
    ENDDO
  ENDDO

ELSEIF (KLX == 3) THEN
  DO JV=1,KVX
    DO JI=1,KIX
      PY(1,JV,JI) = PA (1)*PX(1,JV,JI)+PBS(1)*PX(2,JV,JI)+PCS(1)*PX(3,JV,JI)
      PY(2,JV,JI) = PBI(1)*PX(1,JV,JI)+PA (2)*PX(2,JV,JI)+PBS(2)*PX(3,JV,JI)
      PY(3,JV,JI) = PCI(1)*PX(1,JV,JI)+PBI(2)*PX(2,JV,JI)+PA (3)*PX(3,JV,JI)
    ENDDO
  ENDDO

ELSEIF (KLX == 2) THEN
  DO JV=1,KVX
    DO JI=1,KIX
      PY(1,JV,JI) = PA (1)*PX(1,JV,JI)+PBS(1)*PX(2,JV,JI)
      PY(2,JV,JI) = PBI(1)*PX(1,JV,JI)+PA (2)*PX(2,JV,JI)
    ENDDO
  ENDDO

ELSEIF (KLX == 1) THEN
  DO JV=1,KVX
    DO JI=1,KIX
      PY(1,JV,JI) = PA (1)*PX(1,JV,JI)
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE MXPTMA

