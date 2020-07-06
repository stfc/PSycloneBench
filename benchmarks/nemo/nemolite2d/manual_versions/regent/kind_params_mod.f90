!> Module holding basic KIND parameters
!> This module imitates the naming of dl_esm_inf
!> It uses C interoperable types for safety with Regent.
module kind_params_mod
  use iso_c_binding
  implicit none

  public

  !> Douple precision kind parameter
  integer, parameter :: GO_WP = c_double

  ! Kind type for double precision
  integer, parameter :: GO_DP = c_double

end module kind_params_mod
