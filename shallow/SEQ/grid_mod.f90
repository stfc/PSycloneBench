module grid_mod
  use kind_params_mod
  implicit none

  private

  integer, parameter :: ARAKAWA_C = 0
  integer, parameter :: ARAKAWA_B = 1

  type, public :: grid_type
     !> The type of grid this is (e.g. Arakawa C Grid)
     integer :: grid_name
     !> Total number of grid points
     integer :: npts
     !> Extent of grid in x
     integer :: nx
     !> Extent of grid in y
     integer :: ny
     !> Grid spacing in x
     real(wp) :: dx
     !> Grid spacing in y
     real(wp) :: dy
  end type grid_type

  interface grid_type
     module procedure grid_constructor
  end interface grid_type

contains

  function grid_constructor() result(self)
    implicit none
    type(grid_type), target :: self
    write(*,*) "grid_constructor needs to read namelist file!"

    self%grid_name = ARAKAWA_C

  end function grid_constructor

end module grid_mod
