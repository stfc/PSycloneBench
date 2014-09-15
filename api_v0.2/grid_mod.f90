module grid_mod
  use kind_params_mod
  use region_mod
  use gocean_mod
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

  ! Enumeration of the four possible choices for 
  ! offsetting the grid-point types relative to the T point.
  !> Points to South and West of T point have same 
  !! i,j index (e.g. 'shallow' code)
  integer, public, parameter :: OFFSET_SW = 0
  integer, public, parameter :: OFFSET_SE = 1
  integer, public, parameter :: OFFSET_NW = 2
  !> Points to North and East of T point have same
  !! i,j index (e.g. NEMO code).
  integer, public, parameter :: OFFSET_NE = 3

  ! Enumeration of boundary-condition types
  !> Grid (model domain) has periodic boundary condition
  integer, public, parameter :: BC_PERIODIC = 0
  !> Grid (model domain) has external boundary conditions. This is a
  !! placeholder really as this is a complex area.
  integer, public, parameter :: BC_EXTERNAL = 1
  !> Grid (model domain) has no boundary conditions
  integer, public, parameter :: BC_NONE = 2

  !> The width of the halos we set up for implementing PBCs
  integer, parameter :: HALO_WIDTH_X = 1
  integer, parameter :: HALO_WIDTH_Y = 1

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
     !! is actually simulated. However, we also support the
     !! specification of a model consisting entirely of wet points
     !! with periodic boundary conditions. Since this does not
     !! require a T-mask, we do not allocate this array for that
     !! case.
     integer, allocatable :: tmask(:,:)

     !> The type of boundary conditions applied to the model domain
     !! in the x, y and z dimensions. Note that at this stage
     !! this is really only required for Periodic Boundary
     !! conditions.
     integer, dimension(3) :: boundary_conditions

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

     !> Coordinates of grid (T) points in horizontal plane
     real(wp), allocatable :: xt(:,:), yt(:,:)

  end type grid_type

  interface grid_type
     module procedure grid_constructor
  end interface grid_type

  public grid_init

contains

  !============================================

  !> Basic constructor for the grid type. Full details
  !! are fleshed-out by the grid_init() routine.
  function grid_constructor(grid_name, grid_stagger, &
                            boundary_conditions) result(self)
    implicit none
    integer, intent(in) :: grid_name
    integer, intent(in) :: grid_stagger
    !> The boundary conditions that this field is subject to
    integer, dimension(3), intent(in) :: boundary_conditions
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

    ! Store the boundary conditions that the model domain is
    ! subject to.
    self%boundary_conditions(1:3) = boundary_conditions(1:3)

  end function grid_constructor

  !============================================

  !> Initialise the supplied grid object for a 2D model
  !! consisting of m x n points. Ultimately, this routine should be
  !! general purpose but it is not there yet.
  !! N.B. the definition of m and n (the grid extents) depends on
  !! the type of boundary conditions that the model is subject to.
  !! For periodic boundary conditions they specify the extent of the
  !! simulated region (since we don't require the user to specify 
  !! the halos required to *implement* the PBCs). However, when a
  !! T-mask is used to define the model domain, m and n give the
  !! extents of that mask/grid. This, of necessity, includes boundary
  !! points. Therefore, the actual simulated region has an extent
  !! which is less than m x n.
  !! @param[inout] grid The object to initialise
  !! @param[in] m Extent in x of domain for which we have information
  !! @param[in] n Extent in y of domain for which we have information
  !! @param[in] dxarg Grid spacing in x dimension
  !! @param[in] dyarg Grid spacing in y dimension
  !! @param[in] tmask Array holding the T-point mask which defines
  !!                  the contents of the domain. Need not be
  !!                  supplied if domain is all wet and has PBCs.
  subroutine grid_init(grid, m, n, dxarg, dyarg, tmask)
    implicit none
    type(grid_type), intent(inout) :: grid
    integer,         intent(in)    :: m, n
    real(wp),        intent(in)    :: dxarg, dyarg
    integer, dimension(m,n), intent(in), optional :: tmask
    ! Locals
    integer :: ierr(5)
    integer :: ji, jj
    integer :: xstart, ystart ! Start of internal region of T-pts
    integer :: xstop, ystop ! End of internal region of T-pts

    ! Store the global dimensions of the grid.
    if( present(tmask) )then
       ! A T-mask has been supplied and that tells us everything
       ! about the extent of this model.
       grid%nx = m
       grid%ny = n
    else
       ! No T-mask has been supplied so we assume we're implementing
       ! periodic boundary conditions and allow for halos of width
       ! HALO_WIDTH_{X,Y} here.  Currently we put a halo on all four
       ! sides of our rectangular domain. This is actually unnecessary
       ! - depending on the variable staggering used only one of the
       ! E/W halos and one of the N/S halos are required. However,
       ! this is an optimisation and this framework must be developed
       ! in such a way that that optimisation is supported.
       grid%nx = m + 2*HALO_WIDTH_X
       grid%ny = n + 2*HALO_WIDTH_Y
    end if

    ! For a regular, orthogonal mesh the spatial resolution is constant
    grid%dx = dxarg
    grid%dy = dyarg

    allocate(grid%e1t(grid%nx,grid%ny), grid%e2t(grid%nx,grid%ny), &
             grid%e1u(grid%nx,grid%ny), grid%e2u(grid%nx,grid%ny), &
             stat=ierr(1))
    allocate(grid%e1f(grid%nx,grid%ny), grid%e2f(grid%nx,grid%ny), &
             grid%e1v(grid%nx,grid%ny), grid%e2v(grid%nx,grid%ny), &
             stat=ierr(2)) 
    allocate(grid%e12t(grid%nx,grid%ny), grid%e12u(grid%nx,grid%ny), &
             grid%e12v(grid%nx,grid%ny), stat=ierr(3))
    allocate(grid%gphiu(grid%nx,grid%ny), grid%gphiv(grid%nx,grid%ny), &
             grid%gphif(grid%nx,grid%ny), stat=ierr(4))
    allocate(grid%xt(grid%nx,grid%ny), grid%yt(grid%nx,grid%ny), stat=ierr(5))

    if( any(ierr /= 0, 1) )then
       call gocean_stop('grid_init: failed to allocate arrays')
    end if

    ! Copy-in the externally-supplied T-mask, if any
    if( present(tmask) )then
       allocate(grid%tmask(grid%nx,grid%ny), stat=ierr(1))
       if( ierr(1) /= 0 )then
          call gocean_stop('grid_init: failed to allocate array for T mask')
       end if
       grid%tmask(:,:) = tmask(:,:)
    else
       ! No T-mask supplied. Check that grid has PBCs in both
       ! x and y dimensions otherwise we won't know what to do.
       if( .not. ( (grid%boundary_conditions(1) == BC_PERIODIC) .and. &
                   (grid%boundary_conditions(2) == BC_PERIODIC) ) )then
          call gocean_stop('grid_init: ERROR: No T-mask supplied and '// &
                           'grid does not have periodic boundary conditions!')
       end if
    end if ! T-mask supplied

    ! Use the T mask to determine the dimensions of the
    ! internal, simulated region of the grid.
    ! This call sets grid%simulation_domain.
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
    do jj = 1, grid%ny
       do ji = 1, grid%nx
          grid%e12t(ji,jj) = grid%e1t(ji,jj) * grid%e2t(ji,jj)
       end do
    end do
  
    DO jj = 1, grid%ny
       DO ji = 1, grid%nx
          grid%e12u(ji,jj) = grid%e1u(ji,jj) * grid%e2u(ji,jj)
       END DO
    END DO

    DO jj = 1, grid%ny
       DO ji = 1, grid%nx
          grid%e12v(ji,jj) = grid%e1v(ji,jj) * grid%e2v(ji,jj)
       END DO
    END DO

    ! -here is an f-plane testing case
    ! i.e. the Coriolis parameter is set to a constant value.
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


    if( allocated(grid%tmask) )then

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
       ! This defines the internal region of any T-point field.
       grid%simulation_domain%xstart = 2
       grid%simulation_domain%xstop  = grid%nx - 1
       grid%simulation_domain%ystart = 2
       grid%simulation_domain%ystop  = grid%ny - 1

    else

       ! We don't have a T mask so we must have PBCs in both x and y
       ! dimensions. In this case, the grid dimensions stored in grid%{nx,ny}
       ! have already been adjusted in grid_init() such that they
       ! include the halos required to implement the PBCs.
       grid%simulation_domain%xstart = 2
       grid%simulation_domain%xstop  = grid%nx - 1
       grid%simulation_domain%ystart = 2
       grid%simulation_domain%ystop  = grid%ny - 1

    end if

    grid%simulation_domain%nx =  grid%simulation_domain%xstop -  &
                                 grid%simulation_domain%xstart + 1
    grid%simulation_domain%ny =  grid%simulation_domain%ystop -  &
                                 grid%simulation_domain%ystart + 1

  end subroutine compute_internal_region

  !================================================

end module grid_mod
