module grid_mod
  use kind_params_mod
  use region_mod
  implicit none

  private

  ! Enumeration of possible grid types (we only actually
  ! support ARAKAWA_C at the moment)
  integer, public, parameter :: ARAKAWA_C = 0
  integer, public, parameter :: ARAKAWA_B = 1

  ! Enumeration of the four possible ways of staggering
  ! the fields.
  !> Points to North and East of T point have same
  !! i,j index (e.g. NEMO code).
  integer, public, parameter :: STAGGER_NE = 0
  integer, public, parameter :: STAGGER_NW = 1
  integer, public, parameter :: STAGGER_SE = 2
  !> Points to South and West of T point have same 
  !! i,j index (e.g. 'shallow' code)
  integer, public, parameter :: STAGGER_SW = 3

! If our model grid has total extent (i.e. including pts
! necessary to define the boundaries of the domain) of nx x ny.
! The type of the the points around the edges of the grid
! are determined according to the way in which the fields are
! staggered:
!
! e.g. for a NE stagger:
!
!     1       2       3.....          nx
!     V   F                           V    F
! ny  T   U   T       T       T   U   T    U
!     V   F   V   F   V   F   V   F   V    F
!     T       T   U   T   U   T   U   T    U
!             V   F   V   F   V
! 3   T       T   U   T   U   T       T
!             V   F   V   F   V
! 2   T   U   T   U   T   U   T       T
!     V   F
! 1   T   U   T       T       T       T
!
! That way, we have nx x ny grids of every grid-point type.

  type, public :: grid_type
     !> The type of grid this is (e.g. Arakawa C Grid)
     integer :: name
     !> Specifies the convention by which grid-point
     !! types are indexed.
     integer :: stagger
     !> Total number of grid points
     integer :: npts
     !> Extent of T-point grid in x. Note that this is the whole grid,
     !! not just the region that is simulated.
     integer :: nx
     !> Extent of T-point grid in y. Note that this is the whole grid,
     !! not just the region that is simulated.
     integer :: ny
     !> Grid spacing in x (m)
     real(wp) :: dx
     !> Grid spacing in y (m)
     real(wp) :: dy

     !> Nature of each T point: 1 == wet inside simulated region
     !!                         0 == land
     !!                        -1 == wet outside simulated region
     !! This is the key quantity that determines the region that
     !! is actually simulated.
     integer, allocatable :: tmask(:,:)

     !> Where on the grid our simulated domain sits.
     !! \todo Decide whether this is useful.
     type(region_type) :: simulation_domain

     !> Horizontal scale factors at t point (m)
     real(wp), allocatable :: e1t(:,:), e2t(:,:)
     !> Horizontal scale factors at u point (m)
     real(wp), allocatable :: e1u(:,:), e2u(:,:)
     !> Horizontal scale factors at v point (m)
     real(wp), allocatable :: e1v(:,:), e2v(:,:)  
     !> Horizontal scale factors at f point (m)
     real(wp), allocatable :: e1f(:,:), e2f(:,:)
     !> Unknown \todo Name these fields!
     real(wp), allocatable :: e12t(:,:), e12u(:,:), e12v(:,:)
     !> Latitude of u points
     real(wp), allocatable :: gphiu(:,:)
     !> Latitude of v points
     real(wp), allocatable :: gphiv(:,:)
     !> Latitude of f points
     real(wp), allocatable :: gphif(:,:)

     !> Coordinates of grid points in horizontal plane
     real(wp), allocatable :: xt(:,:), yt(:,:)

  end type grid_type

  interface grid_type
     module procedure grid_constructor
  end interface grid_type

  public grid_init

contains

  !============================================

  function grid_constructor(grid_name, grid_stagger) result(self)
    implicit none
    integer, intent(in) :: grid_name
    integer, intent(in) :: grid_stagger
    type(grid_type), target :: self

    ! This case statement is mainly to check that the caller
    ! has specified a valid value for grid_name.
    select case(grid_name)

    case(ARAKAWA_C)
       self%name = ARAKAWA_C
    case(ARAKAWA_B)
       self%name = ARAKAWA_B
    case default
       write(*,*) 'grid_constructor: ERROR: unsupported grid type: ', &
                  grid_name
       stop
    end select

    ! Ditto for the choice of grid staggering
    select case(grid_stagger)

    case(STAGGER_NE)
       self%stagger = STAGGER_NE
    case(STAGGER_NW)
       self%stagger = STAGGER_NW
    case(STAGGER_SE)
       self%stagger = STAGGER_SE
    case(STAGGER_SW)
       self%stagger = STAGGER_SW
    case default
       write(*,*) 'grid_constructor: ERROR: unsupported grid stagger: ', &
                  grid_stagger
       stop
    end select

  end function grid_constructor

  !============================================

  !> Initialise the supplied grid object for a 2D model
  !! consisting of m x n points. Ultimately, this routine should be
  !! general purpose but it is currently evolving towards that...
  !! @param[inout] grid The object to initialise
  !! @param[in] m Extent in x of domain for which we have information
  !! @param[in] n Extent in y of domain for which we have information
  !! @param[in] dxarg Grid spacing in x dimension
  !! @param[in] dyarg Grid spacing in y dimension
  !! @param[in] tmask Array holding the T-point mask which defines
  !!                  the contents of the domain.
  subroutine grid_init(grid, m, n, dxarg, dyarg, tmask)
    implicit none
    type(grid_type), intent(inout) :: grid
    integer, intent(in) :: m, n
    real(wp), intent(in) :: dxarg, dyarg
    integer, dimension(:,:), intent(in) :: tmask
    ! Locals
    integer :: ierr(5)
    integer :: ji, jj
    integer :: xstart, ystart ! Start of internal region of T-pts
    integer :: xstop, ystop ! End of internal region of T-pts

    select case(grid%name)

    case(ARAKAWA_C)
       grid%nx = m
       grid%ny = n
    case default
       stop 'grid_init: ERROR: only Arakawa C grid implemented!'
    end select

    ! For a regular, orthogonal mesh the spatial resolution is constant
    grid%dx = dxarg
    grid%dy = dyarg

    allocate(grid%e1t(m,n), grid%e2t(m,n), &
             grid%e1u(m,n), grid%e2u(m,n), &
             stat=ierr(1))
    allocate(grid%e1f(m,n), grid%e2f(m,n), &
             grid%e1v(m,n), grid%e2v(m,n),   &
             stat=ierr(2)) 
    allocate(grid%e12t(m,n), grid%e12u(m,n), grid%e12v(m,n), &
             stat=ierr(3))
    allocate(grid%gphiu(m,n), grid%gphiv(m,n), grid%gphif(m,n), &
             stat=ierr(4))
    allocate(grid%xt(m,n), grid%yt(m,n), grid%tmask(m,n), stat=ierr(5))

    if( any(ierr /= 0, 1) )then
       stop 'grid_init: failed to allocate arrays'
    end if

    ! Copy-in the externally-supplied T-mask
    grid%tmask(:,:) = tmask(:,:)

    ! Use the T mask to determine the dimensions of the
    ! internal, simulated region of the grid.
    call compute_internal_region(grid)

    ! Initialise the horizontal scale factors for a regular,
    ! orthogonal mesh. (Constant spatial resolution.)
    grid%e1t(:, :)   = grid%dx
    grid%e2t(:, :)   = grid%dy

    grid%e1u(:, :)   = grid%dx
    grid%e2u(:, :)   = grid%dy

    grid%e1v(:, :)   = grid%dx
    grid%e2v(:, :)   = grid%dy

    grid%e1f(:, :)   = grid%dx
    grid%e2f(:, :)   = grid%dy

    ! calculate t,u,v cell area
    do jj = 1, n
       do ji = 1, m
          grid%e12t(ji,jj) = grid%e1t(ji,jj) * grid%e2t(ji,jj)
       end do
    end do
  
    DO jj = 1, n
       DO ji = 1, m
          grid%e12u(ji,jj) = grid%e1u(ji,jj) * grid%e2u(ji,jj)
       END DO
    END DO

    DO jj = 1, n
       DO ji = 1, m
          grid%e12v(ji,jj) = grid%e1v(ji,jj) * grid%e2v(ji,jj)
       END DO
    END DO

    ! -here is an f-plane testing case
    grid%gphiu(:, :) = 50._wp
    grid%gphiv(:, :) = 50._wp
    grid%gphif(:, :) = 50._wp

    ! Co-ordinates of the T points
    xstart = grid%simulation_domain%xstart
    xstop  = grid%simulation_domain%xstop
    ystart = grid%simulation_domain%ystart
    ystop  = grid%simulation_domain%ystop
    grid%xt(xstart, :) = 0.0_wp + 0.5_wp * grid%e1t(xstart,:)
    grid%yt(:,ystart) = 0.0_wp + 0.5_wp * grid%e2t(:,ystart)

    DO ji = xstart+1, xstop
      grid%xt(ji,ystart:ystop) = grid%xt(ji-1, ystart:ystop) + grid%dx
    END DO
            
    DO jj = ystart+1, ystop
      grid%yt(xstart:xstop,jj) = grid%yt(xstart:xstop, jj-1) + grid%dy
    END DO


  end subroutine grid_init

  !================================================

  !> Use the T-mask to deduce the inner or simulated region
  !! of the supplied grid.
  subroutine compute_internal_region(grid)
    implicit none
    type(grid_type), intent(inout) :: grid

    ! Here we will loop over the grid points, looking for
    ! the first occurrence of wet points.
    ! However, for the moment we just hardwire the routine
    ! to return results appropriate for a T mask that has
    ! a shell of unit depth of boundary/external points:

    ! i= 1           nx 
    !    b   b   b   b   ny
    !    b   x   x   b  
    !    b   x   x   b  
    !    b   x   x   b   
    !    b   b   b   b   1
    !                    j

    ! The actual part of this domain that is simulated. The outer-most 
    ! rows and columns of T points are not in the domain but are needed
    ! to specify the type of boundary (whether hard or open).
    !> \todo Generate the bounds of the simulation domain by
    !! examining the T-point mask.
    !> \todo Isn't this just the same as the internal region of
    !! any T-point field??
    grid%simulation_domain%xstart = 2
    grid%simulation_domain%xstop  = grid%nx - 1
    grid%simulation_domain%ystart = 2
    grid%simulation_domain%ystop  = grid%ny - 1
    grid%simulation_domain%nx = grid%nx - 2
    grid%simulation_domain%ny = grid%ny - 2

  end subroutine compute_internal_region

  !================================================

end module grid_mod
