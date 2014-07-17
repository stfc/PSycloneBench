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

     !> Horizontal scale factors
     real(wp), allocatable :: e1t(:,:), e2t(:,:)
     real(wp), allocatable :: e1u(:,:), e2u(:,:)
     real(wp), allocatable :: e1v(:,:), e2v(:,:)  
     real(wp), allocatable :: e1f(:,:), e2f(:,:)
     real(wp), allocatable :: e12t(:,:), e12u(:,:), e12v(:,:)

     !> Latitude of grid points (?)
     real(wp), allocatable :: gphiu(:,:), gphiv(:,:), gphif(:,:)

     !> Coordinates of grid points in horizontal plane
     real(wp), allocatable :: xt(:,:), yt(:,:)

  end type grid_type

  interface grid_type
     module procedure grid_constructor
  end interface grid_type

  public grid_init

contains

  !============================================

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
    ! Locals
    integer :: ierr(5)
    integer :: ji, jj

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

    allocate(grid%e1t(m,n), grid%e2t(m,n), grid%e1u(0:m,n), grid%e2u(0:m,n), STAT=ierr(1))
    allocate(grid%e1f(0:m,0:n), grid%e2f(0:m,0:n), grid%e1v(m,0:n), grid%e2v(m,0:n), STAT=ierr(2)) 
    allocate(grid%e12t(m,n), grid%e12u(0:m,n), grid%e12v(m,0:n), STAT=ierr(3))
    allocate(grid%gphiu(0:m,n), grid%gphiv(m,0:n), grid%gphif(0:m,0:n), STAT=ierr(4))
    allocate(grid%xt(m,n), grid%yt(m,n), STAT=ierr(5))
    if( any(ierr /= 0, 1) )then
       stop 'grid_init: failed to allocate arrays'
    end if

    grid%e1t(1:jpi, 1:jpj)   = grid%dx
    grid%e2t(1:jpi, 1:jpj)   = grid%dy
    grid%e1u(0:jpi, 1:jpj)   = grid%dx
    grid%e2u(0:jpi, 1:jpj)   = grid%dy
    grid%e1v(1:jpi, 0:jpj)   = grid%dx
    grid%e2v(1:jpi, 0:jpj)   = grid%dy
    grid%e1f(0:jpi, 0:jpj)   = grid%dx
    grid%e2f(0:jpi, 0:jpj)   = grid%dy

     
    ! -here is an f-plane testing case
    grid%gphiu(0:jpi, 1:jpj) = 50._wp
    grid%gphiv(1:jpi, 0:jpj) = 50._wp
    grid%gphif(0:jpi, 0:jpj) = 50._wp

    grid%xt(1,1) = 0.0_wp + 0.5_wp * grid%e1t(1,1)
    grid%yt(1,1) = 0.0_wp + 0.5_wp * grid%e2t(1,1)

    DO ji = 2, jpi
      grid%xt(ji,1:jpj) = grid%xt(ji-1, 1:jpj) + grid%dx
    END DO
            
    DO jj = 2, jpj
      grid%yt(1:jpi,jj) = grid%yt(1:jpi, jj-1) + grid%dy
    END DO

  end subroutine grid_init

end module grid_mod
