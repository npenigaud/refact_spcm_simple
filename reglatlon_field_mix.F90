MODULE REGLATLON_FIELD_MIX

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMLUN   , ONLY : NULOUT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
IMPLICIT NONE
SAVE

!#include "abor1.intfb.h

! Define a regular lat-lon field. 
! Assume that:
!       latitude and longitude spacing is regular
!       That the grid spans the entire longitude circle
!
!     Hans Hersbach  9/12/2010
!     M. Fisher   7-March-2012 Use DEALLOCATE_IF_ASSOCIATED + remove semicolons


TYPE REGLATLON_FIELD
  INTEGER(KIND=JPIM)       :: NLAT               ! # of latitudes   ; >=1
  INTEGER(KIND=JPIM)       :: NLON               ! # of longitudes  ; >=1
  REAL(KIND=JPRB)          :: DLAT               ! Latitude  increment (degrees), can be positve/negative
  REAL(KIND=JPRB)          :: DLON               ! Longitude increment (degrees), is assumed to be positve
  REAL(KIND=JPRB) ,POINTER :: PFLD(:,:)=>NULL()  ! Field contents as (lat,lon)
  REAL(KIND=JPRB) ,POINTER :: PLAT(:)  =>NULL()  ! List of NLAT  Latitudes (degrees)
  REAL(KIND=JPRB) ,POINTER :: PSIN(:)  =>NULL()  ! List of NLAT  SIN(Latitudes)
  REAL(KIND=JPRB) ,POINTER :: PLON(:)  =>NULL()  ! List of NLON Longitudes (degrees)
END TYPE REGLATLON_FIELD

END
