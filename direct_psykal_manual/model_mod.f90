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

  REAL(wp), ALLOCATABLE, save :: ht(:,:), hu(:,:), hv(:,:), hf(:,:) 
  
  REAL(wp), ALLOCATABLE, save :: sshb(:,:), sshb_u(:,:), sshb_v(:,:)
  !> Sea-surface height at T, U and V points
  REAL(wp), ALLOCATABLE, save :: sshn(:,:), sshn_u(:,:), sshn_v(:,:)
  REAL(wp), ALLOCATABLE, save :: ssha(:,:), ssha_u(:,:), ssha_v(:,:)
  
  REAL(wp), ALLOCATABLE, save :: un(:,:),  vn(:,:)
  REAL(wp), ALLOCATABLE, save :: ua(:,:),  va(:,:)

CONTAINS

  !================================================

  subroutine model_init(grid)
    use gocean2d_io_mod
    implicit none
    type(grid_type), intent(inout) :: grid

    !> Problem size, read from namelist
    integer :: jpiglo, jpjglo
    real(wp) :: dx, dy

    integer :: ierr(5)

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

    ! Allocate all our fields with the same extents as the T-mask defining
    ! the whole model domain (i.e. including points that are outside the
    ! simulated region).
    ALLOCATE(ht(jpi,jpj), hu(jpi,jpj), hv(jpi,jpj), hf(jpi,jpj), STAT=ierr(1))

    ALLOCATE(sshb(jpi,jpj), sshb_u(jpi,jpj), sshb_v(jpi,jpj), STAT=ierr(2))
    ALLOCATE(sshn(jpi,jpj), sshn_u(jpi,jpj), sshn_v(jpi,jpj), STAT=ierr(3))
    ALLOCATE(ssha(jpi,jpj), ssha_u(jpi,jpj), ssha_v(jpi,jpj), STAT=ierr(4))

    ALLOCATE(un(jpi,jpj), vn(jpi,jpj), ua(jpi,jpj), va(jpi,jpj), STAT=ierr(5))

    IF(ANY(ierr /= 0, 1)) STOP "Failed to allocate solution arrays"

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

    DEALLOCATE(ht, hu, hv, hf)
    DEALLOCATE(sshb, sshb_u, sshb_v)
    DEALLOCATE(sshn, sshn_u, sshn_v)
    DEALLOCATE(ssha, ssha_u, ssha_v)
    DEALLOCATE(un, vn, ua, va)

  END SUBROUTINE model_dealloc

END MODULE model_mod
