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

  public grid_init

contains

  function grid_constructor(grid_name) result(self)
    implicit none
    integer, intent(in) :: grid_name
    type(grid_type), target :: self

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

  !============================================

  !> Initialise the supplied grid object for a 2D model
  !! consisting of m x n points.
  !! @param[inout] grid The object to initialise
  !! @param[in] m Extent of model in x dimension
  !! @param[in] n Extent of model in y dimension
  !! @param[in] dxarg Grid spacing in x dimension
  !! @param[in] dyarg Grid spacing in y dimension
  subroutine grid_init(grid, m, n, dxarg, dyarg)
    implicit none
    type(grid_type), intent(inout) :: grid
    integer, intent(in) :: m, n
    real(wp), intent(in) :: dxarg, dyarg

    select case(grid%grid_name)

    case(ARAKAWA_C)
       ! For an Arakawa C grid we need a grid that has extent 
       ! 1 greater than the extent of the model fields in
       ! order to allow for variable staggering 
       grid%nx = m+1
       grid%ny = n+1
    case default
       stop 'grid_init: ERROR: only Arakawa C grid implemented!'
    end select

    ! For a regular, orthogonal mesh the spatial resolution is constant
    grid%dx = dxarg
    grid%dy = dyarg

  end subroutine grid_init

end module grid_mod
