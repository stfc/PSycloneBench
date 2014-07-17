!> Module for physical/mathematical parameters
MODULE physical_params_mod
  USE kind_params_mod
  IMPLICIT none

  PUBLIC

  !> Pi
  REAL(wp), PARAMETER :: pi    = 3.1415926535897932_wp  
  !> Acceleration due to gravity (ms^-2)
  REAL(wp), PARAMETER :: g     = 9.80665_wp
  !> earth rotation speed (s^(-1))
  REAL(wp), PARAMETER :: omega = 7.292116e-05_wp
  !> degree to radian
  REAL(wp), PARAMETER :: d2r   = pi / 180._wp                       

END MODULE physical_params_mod
