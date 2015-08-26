MODULE model_mod
  use kind_params_mod
  use grid_mod
  use field_mod
  use timing_mod, ONLY: timer_init, timer_report
  implicit none

  public

  !> start-end and record time steps
  integer, save  :: nit000, nitend, irecord 
  !> type (source?) of grid
  integer, save  :: jphgr_msh

  real(wp), save :: rdt             !< time step
  REAL(wp), save :: cbfr            !< bottom friction coefficient
  REAL(wp), save :: visc            !< backgroud/constant viscosity 

  REAL(wp), save :: dep_const       !< constant depth

  
CONTAINS

  !================================================

  subroutine model_init(grid)
    use gocean2d_io_mod
    use gocean_mod, only: gocean_init
    implicit none
    type(grid_type), intent(inout) :: grid

    !> Problem size, read from namelist
    integer :: jpiglo, jpjglo
    real(wp) :: dx, dy
    integer  :: ierr
    integer, dimension(:,:), allocatable :: tmask

    ! Initialise timing system
    call timer_init()

    ! Read model configuration from namelist
    call read_namelist(jpiglo, jpjglo, dx, dy, &
                       nit000, nitend, irecord, &
                       jphgr_msh, dep_const, rdt, cbfr, visc)

    ! Set-up the T mask. This defines the model domain.
    allocate(tmask(jpiglo,jpjglo), stat=ierr)
    if(ierr /= 0)then
       stop 'model_init: failed to allocate T mask'
    end if

    call setup_tpoints_mask(jpiglo, jpjglo, tmask)

    ! Having specified the T points mask, we can set up mesh parameters
    call grid_init(grid, jpiglo, jpjglo, dx, dy, tmask)

    ! Initialise model IO 'system'
    call model_write_init(jpiglo, jpjglo, irecord)

    ! Clean-up. T-mask is now a part of the grid object.
    deallocate(tmask)

    ! Any other GOcean initialisation
    call gocean_init()

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

  !================================================

  subroutine setup_tpoints_mask(jpi, jpj, tmask)
    implicit none
    integer, intent(in) :: jpi, jpj
    integer, dimension(:,:), intent(inout) :: tmask
    ! Locals
    integer :: ji, jj

    ! jphgr_msh is read from the namelist file
    !jphgr_msh = 0    ! read in grid definition from file
    !jphgr_msh = 1    ! define manually 

    select case( jphgr_msh)

    case(0) ! read in grid from a coordinate file

       ! to be added
       stop "Reading the grid definition from file is not yet supported!"
       ! add reading data part here
       ! the following variables/arrays are needed:
       ! jpi, jpj, 
       ! pt(0:jpi+1)
       ! e1t(jpi,    jpj), e2t(jpi,    jpj)
       ! e1u(0:jpi,  jpj), e2u(0:jpi,  jpj)
       ! e1v(jpi,  0:jpj), e2v(jpi,  0:jpj)
       ! e1f(0:jpi,0:jpj), e2f(0:jpi,0:jpj)
       ! gphiu(0:jpi,jpj), gphiv(jpi,0:jpj), gphif(0:jpi,0:jpj) 
       ! xt(jpi,jpj), yt(jpi,jpj)
       ! ht(jpi,jpj), hu(jpi,jpj), hv(jpi,jpj)

    case(1)

       !##### a manually defined grid

       ! -size of each grid cell
       ! -depth on each T points
       ! -grid dimension

       ! Define Model solid/open Boundaries via the properties of t-cells

       tmask(:,:) = 1                 ! all inner cells

       ! -define solid/open boundaries
       DO jj = 1, jpj
          tmask(  1, jj) = 0              ! west solid boundary
          tmask(jpi, jj) = 0              ! east solid boundary
       END Do

       DO ji = 1, jpi
          tmask(ji,jpj) = 0               ! north solid boundary
       END Do

       DO ji = 1, jpi
          tmask(ji,1) = -1                    ! south open boundary
       END Do

    CASE DEFAULT
       ! undefined jphgr_msh value
       ! add interrupt here
       STOP "Wrong grid definition type, check your setup !!!!"
    END SELECT

  end subroutine setup_tpoints_mask

END MODULE model_mod
