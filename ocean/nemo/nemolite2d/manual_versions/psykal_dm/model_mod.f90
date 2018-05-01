MODULE model_mod
  use kind_params_mod
  use grid_mod
  use field_mod
  use dl_timer, ONLY: timer_init, timer_report
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
    use subdomain_mod, only: subdomain_type, decompose
    use parallel_mod, only: get_rank
    implicit none
    type(grid_type), intent(inout) :: grid

    !> Global problem size, read from namelist
    integer :: jpiglo, jpjglo
    real(wp) :: dx, dy
    integer  :: ierr
    integer, dimension(:,:), allocatable :: tmask
    type(subdomain_type), allocatable :: tile_list(:)
    integer :: rank

    ! Initialise timing system
    call timer_init()

    ! Read model configuration from namelist
    call read_namelist(jpiglo, jpjglo, dx, dy, &
                       nit000, nitend, irecord, &
                       jphgr_msh, dep_const, rdt, cbfr, visc)

    ! Work out the decomposition (uses the number of MPI ranks by
    ! default)
    tile_list = decompose(jpiglo, jpjglo)
    rank = get_rank()
    ! Take a copy of this process' subdomain definition
    grid%subdomain = tile_list(rank)

    ! Set-up the T mask. This defines the model domain on this process.
    call setup_tpoints_mask(grid%subdomain, jpiglo, jpjglo, tmask)

    ! Having specified the T points mask, we can set up mesh parameters
    call grid_init(grid, jpiglo, jpjglo, dx, dy, tmask)

    call dump_tmask(tmask, grid)

    ! Initialise model IO 'system'
    call model_write_init(jpiglo, jpjglo, irecord)

    ! Clean-up. T-mask is now a part of the grid object.
    deallocate(tmask)

  end subroutine model_init

  !================================================

  SUBROUTINE model_finalise()
    use parallel_mod, only: parallel_finalise
    use gocean2d_io_mod, only: model_write_finalise
    IMPLICIT none

    call model_write_finalise()

    call timer_report()

    call model_dealloc()

    call parallel_finalise()

  end subroutine model_finalise

  !================================================

  subroutine model_dealloc()
    implicit none

  END SUBROUTINE model_dealloc

  !================================================

  subroutine setup_tpoints_mask(subdomain, jpi, jpj, tmask)
    !> Allocate and set-up the T-points mask that defines the
    !! simulation domain on the current process.
    use gocean_mod, only: gocean_stop
    use subdomain_mod, only: subdomain_type
    implicit none
    ! The (sub)domain dealt with by the calling process
    type(subdomain_type), intent(in) :: subdomain
    ! Dimensions of global domain
    integer, intent(in) :: jpi, jpj
    integer, dimension(:,:), allocatable, intent(inout) :: tmask
    ! Locals
    integer :: ji, jj, ierr

    allocate(tmask(subdomain%nx, subdomain%ny), stat=ierr)
    if(ierr /= 0)then
       call gocean_stop('model_init: failed to allocate T mask')
    end if

    ! jphgr_msh is read from the namelist file
    !jphgr_msh = 0    ! read in grid definition from file
    !jphgr_msh = 1    ! define manually 

    select case( jphgr_msh)

    case(0) ! read in grid from a coordinate file

       ! to be added
       call gocean_stop("Reading the grid definition from file is "// &
                        "not yet supported!")
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

       tmask(:,:) = -1 ! Default all cells to being outside the domain

       ! Mark all inner cells
       tmask(subdomain%internal%xstart:subdomain%internal%xstop, &
             subdomain%internal%ystart:subdomain%internal%ystop) = 1
             
       ! -define solid/open boundaries
       if(subdomain%xstart == 1)then
          do jj = 1, subdomain%internal%ny
             tmask(subdomain%internal%xstart, jj) = 0 ! west solid boundary
          end do
       end if

       if(subdomain%xstop == jpi)then
          do jj = 1, subdomain%internal%ny
             tmask(subdomain%internal%xstop, jj) = 0 ! east solid boundary
          end do
       end if

       if(subdomain%ystop == jpj)then
          do ji = 1, subdomain%internal%nx
             tmask(ji, subdomain%internal%ystop) = 0 ! north solid boundary
          end do
       end if

       if(subdomain%ystart == 1)then
          do ji = 1, subdomain%internal%nx
             tmask(ji, subdomain%internal%ystart) = -1 ! south open boundary
          end do
       end if

    case DEFAULT
       ! undefined jphgr_msh value
       ! add interrupt here
       call gocean_stop("Wrong grid definition type (jphgr_msh must be 0 "// &
                        "or 1), check your namelist file.")
    end select

  end subroutine setup_tpoints_mask

  subroutine dump_tmask(tmask, grid)
    use parallel_mod, only: get_rank
    use grid_mod, only: grid_type
    implicit none
    type(grid_type), intent(in) :: grid
    integer, dimension(:,:), intent(in) :: tmask
    ! Locals
    integer :: ji, jj
    character(len=20) :: fname

    write(fname, '(I5.5)') get_rank()
    open(21, file='tmask_'//TRIM(fname)//'.dat', STATUS='UNKNOWN', &
         action='write')

    ! Loop over 'internal' T points
    DO jj = grid%subdomain%internal%ystart, grid%subdomain%internal%ystop, 1
       DO ji = grid%subdomain%internal%xstart, grid%subdomain%internal%xstop, 1

          write(21,'(2e16.7,1x,I3)') grid%xt(ji,jj), grid%yt(ji,jj), &
               tmask(ji,jj)
       END DO
       WRITE(21,*)
    END DO
          
    close(21)
  end subroutine dump_tmask

END MODULE model_mod