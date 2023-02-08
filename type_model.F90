! (C) Copyright 2011- ECMWF.
! (C) Copyright 2011- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!> Handle model configuration for the IFS model

MODULE TYPE_MODEL

USE PARKIND1                      , ONLY : JPIM
USE MODEL_GENERAL_CONF_MOD        , ONLY : MODEL_GENERAL_CONF_TYPE
USE MODEL_DYNAMICS_MOD            , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_DIAGNOSTICS_MOD         , ONLY : MODEL_DIAGNOSTICS_TYPE
USE YOMCST                        , ONLY : TCST

IMPLICIT NONE

TYPE, PUBLIC :: MODEL
!!  PRIVATE
  TYPE(MODEL_GENERAL_CONF_TYPE)         :: YRML_GCONF
  TYPE(MODEL_DYNAMICS_TYPE)             :: YRML_DYN
  TYPE(MODEL_DIAGNOSTICS_TYPE)          :: YRML_DIAG
  TYPE(TCST), POINTER                   :: YRCST => NULL()
END TYPE MODEL

END MODULE TYPE_MODEL
