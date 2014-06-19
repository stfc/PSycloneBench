module grid_mod
  use kind_params_mod
  implicit none

  private

  integer, public, parameter :: ARAKAWA_C = 0
  integer, public, parameter :: ARAKAWA_B = 1

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

  function grid_constructor(grid_name) result(self)
    implicit none
    integer, intent(in) :: grid_name
    type(grid_type), target :: self

    write(*,*) "grid_constructor needs to read namelist file!"

    ! This case statement is mainly to check that the caller
    ! has specified a valid value for grid_name.
    select case(grid_name)

    case(ARAKAWA_C)
       self%grid_name = ARAKAWA_C
    case(ARAKAWA_B)
       self%grid_name = ARAKAWA_B
    case default
       write(*,*) 'grid_constructor: ERROR: unsupported grid type: ',grid_name
       stop
    end select

  end function grid_constructor

end module grid_mod
