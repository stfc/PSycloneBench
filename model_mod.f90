MODULE model_mod
  use kind_params_mod
  use grid_mod
  use field_mod
  use timing_mod, ONLY: timer_init, timer_report
  implicit none

  public

  !> Problem size, read from namelist
  !! \todo Will live with the grid object
  integer, save  :: jpi, jpj
  !> start-end and record time steps
  integer, save  :: nit000, nitend, irecord 
  !> type (source?) of grid
  integer, save  :: jphgr_msh

  REAL(wp), save :: rdt                       !< time step
  REAL(wp), save :: cbfr                      !< bottom friction coefficient
  REAL(wp), save :: visc                      !< backgroud/constant viscosity 

  REAL(wp), save :: dep_const                 !< constant depth

  
CONTAINS

  !================================================

  subroutine model_init(grid)
    use gocean2d_io_mod
    implicit none
    type(grid_type), intent(inout) :: grid

    !> Problem size, read from namelist
    integer :: jpiglo, jpjglo
    real(wp) :: dx, dy

    ! Initialise timing system
    call timer_init()

    ! Read model configuration from namelist
    call read_namelist(jpiglo, jpjglo, dx, dy, &
                       nit000, nitend, irecord, &
                       jphgr_msh, dep_const, rdt, cbfr, visc)

    ! Set up mesh parameters
    call grid_init(grid, jpiglo, jpjglo, dx, dy)

    ! Store the grid dimensions in module variables - this is a temporary
    ! fix prior to carrying everything around in the grid object.
    jpi = grid%nx
    jpj = grid%ny

    ! Initialise model IO 'system'
    call model_write_init(jpiglo, jpjglo, irecord)

  end subroutine model_init

  !================================================

  SUBROUTINE model_finalise()
    use gocean2d_io_mod, only: model_write_finalise
    IMPLICIT none

    CALL model_write_finalise()

    CALL timer_report()

    CALL model_dealloc()
  
  end subroutine model_finalise

  !================================================

  subroutine model_dealloc()
    implicit none

  END SUBROUTINE model_dealloc

END MODULE model_mod
