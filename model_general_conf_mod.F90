! (C) Copyright 2017- ECMWF.
! (C) Copyright 2017- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE MODEL_GENERAL_CONF_MOD
  USE YOMRIP       , ONLY : TRIP
  IMPLICIT NONE

  TYPE MODEL_GENERAL_CONF_TYPE

    TYPE(TRIP)              :: YRRIP                   !! TEMPORARY TREATMENT OF TIME, SHOULD CHANGE AT CY45

  END TYPE MODEL_GENERAL_CONF_TYPE

END MODULE MODEL_GENERAL_CONF_MOD
