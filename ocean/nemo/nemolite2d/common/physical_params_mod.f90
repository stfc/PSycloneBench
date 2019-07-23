!> Module for physical/mathematical parameters
MODULE physical_params_mod
  USE kind_params_mod
  IMPLICIT none

  PUBLIC

  !> Pi
  REAL(go_wp), PARAMETER :: pi    = 3.1415926535897932_go_wp  
  !> Acceleration due to gravity (ms^-2)
  REAL(go_wp), PARAMETER :: g     = 9.80665_go_wp
  !> earth rotation speed (s^(-1))
  REAL(go_wp), PARAMETER :: omega = 7.292116e-05_go_wp
  !> degree to radian
  REAL(go_wp), PARAMETER :: d2r   = pi / 180._go_wp                       

END MODULE physical_params_mod
