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

  real(go_wp), save :: rdt             !< time step
  REAL(go_wp), save :: cbfr            !< bottom friction coefficient
  REAL(go_wp), save :: visc            !< backgroud/constant viscosity 

  REAL(go_wp), save :: dep_const       !< constant depth

  
CONTAINS

  !================================================

  subroutine model_init(grid)
    use gocean2d_io_mod
    use decomposition_mod, only: decomposition_type
    use subdomain_mod, only: decompose
    use parallel_mod, only: get_rank, set_proc_grid, map_comms
    implicit none
    type(grid_type), intent(inout) :: grid

    !> Global problem size, read from namelist
    integer :: jpiglo, jpjglo
    real(go_wp) :: dx, dy
    integer  :: ierr
    integer, dimension(:,:), allocatable :: tmask
    type(decomposition_type) :: decomp
    integer :: rank

    ! Initialise timing system
    call timer_init()

    ! Read model configuration from namelist
    call read_namelist(jpiglo, jpjglo, dx, dy, &
                       nit000, nitend, irecord, &
                       jphgr_msh, dep_const, rdt, cbfr, visc)

    ! Work out the decomposition of the global domain (uses the number
    ! of MPI ranks by default)
    decomp = decompose(jpiglo, jpjglo)
    call set_proc_grid(decomp%nx, decomp%ny)

    ! Set-up the T mask on this process. This defines the model domain.
    rank = get_rank()
    call setup_tpoints_mask(decomp%subdomains(rank), jpiglo, jpjglo, tmask)

    ! Having specified the T points mask, we can set up mesh parameters
    call grid_init(grid, decomp, dx, dy, tmask)

    ! Write generated, decomposed T-mask out to file
    call dump_tmask(grid, raw=.false.)

    call map_comms(decomp, tmask, .false., ierr)

    ! Initialise model IO 'system'
    call model_write_init(jpiglo, jpjglo, irecord)

    ! Clean-up. T-mask is now a part of the grid object.
    deallocate(tmask)

  end subroutine model_init

  !================================================

  SUBROUTINE model_finalise()
    use gocean2d_io_mod, only: model_write_finalise
    IMPLICIT none

    call model_write_finalise()

    call timer_report()

    call model_dealloc()

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
    use decomposition_mod, only: subdomain_type
    implicit none
    ! The (sub)domain dealt with by the calling process
    type(subdomain_type), intent(in) :: subdomain
    ! Dimensions of global domain
    integer, intent(in) :: jpi, jpj
    integer, dimension(:,:), allocatable, intent(inout) :: tmask
    ! Locals
    integer :: ji, jj, ierr

    allocate(tmask(subdomain%global%nx, subdomain%global%ny), stat=ierr)
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

       ! A manually defined T mask
       ! Define Model solid/open Boundaries via the properties of t-cells

       tmask(:,:) = 0 ! Default all cells to being dry/outside the domain

       ! Mark all inner (wet) cells
       tmask(subdomain%internal%xstart:subdomain%internal%xstop, &
             subdomain%internal%ystart:subdomain%internal%ystop) = 1
             
       ! -define solid/open boundaries
       if(subdomain%global%xstart == 1)then
          ! West solid boundary
          do jj = subdomain%internal%ystart, subdomain%internal%ystop
             tmask(1:subdomain%internal%xstart-1, jj) = 0
          end do
       else
          ! Subdomain is not at the westernmost extent of the global
          ! domain so boundary is wet.
          do jj = subdomain%internal%ystart, subdomain%internal%ystop
             tmask(1:subdomain%internal%xstart-1, jj) = 1
          end do          
       end if

       if(subdomain%global%xstop == jpi)then
           ! East solid boundary
          do jj = subdomain%internal%ystart, subdomain%internal%ystop
             tmask(subdomain%internal%xstop+1:, jj) = 0
          end do
       else
          ! Subdomain is not at the easternmost extent of the global
          ! domain so boundary is wet.
          do jj = subdomain%internal%ystart, subdomain%internal%ystop
             tmask(subdomain%internal%xstop+1:, jj) = 1
          end do
       end if

       if(subdomain%global%ystop == jpj)then
          ! North solid boundary
          do ji = subdomain%internal%xstart, subdomain%internal%xstop
             tmask(ji, subdomain%internal%ystop+1) = 0
          end do
       else
          ! Subdomain is not at the northernmost extent of the domain so
          ! does not have a solid boundary
          tmask(subdomain%internal%xstart+1:subdomain%internal%xstop-1, &
                subdomain%internal%ystop:) = 1
       end if

       if(subdomain%global%ystart == 1)then
          ! South open boundary
          do ji = subdomain%internal%xstart, subdomain%internal%xstop
             ! Make the open boundary outside the computational domain
             ! (i.e. ystart - 1)
             tmask(ji, subdomain%internal%ystart-1) = -1
          end do
       else
          ! Subdomain is not at the southernmost extent of the domain so
          ! does not have an open boundary
          tmask(subdomain%internal%xstart+1:subdomain%internal%xstop-1, &
                1:subdomain%internal%ystart-1) = 1
       end if

    case DEFAULT
       ! undefined jphgr_msh value
       ! add interrupt here
       call gocean_stop("Wrong grid definition type (jphgr_msh must be 0 "// &
                        "or 1), check your namelist file.")
    end select

  end subroutine setup_tpoints_mask

  subroutine dump_tmask(grid, raw)
    use parallel_mod, only: get_rank
    use grid_mod, only: grid_type
    implicit none
    type(grid_type), intent(in) :: grid
    logical, optional, intent(in) :: raw
    ! Locals
    logical :: lraw
    integer :: ji, jj
    character(len=20) :: fname

    write(fname, '(I5.5)') get_rank()
    open(21, file='tmask_'//TRIM(fname)//'.dat', STATUS='UNKNOWN', &
         action='write')

    if(present(raw))then
       lraw = raw
    else
       lraw = .false.
    end if
    
    if(lraw)then
       do jj = 1, grid%subdomain%global%ny
          do ji = 1, grid%subdomain%global%nx
             write(21,'(3(I4,1x))') ji, jj, grid%tmask(ji,jj)
          end do
          write(21,*)
       end do
    else
       ! Loop over 'internal' T points
       do jj = grid%subdomain%internal%ystart, grid%subdomain%internal%ystop, 1
          do ji = grid%subdomain%internal%xstart, grid%subdomain%internal%xstop, 1

             write(21,'(2e16.7,1x,I3)') grid%xt(ji,jj), grid%yt(ji,jj), &
                  grid%tmask(ji,jj)
          end do
          write(21,*)
       end do
    end if
    close(21)
  end subroutine dump_tmask

END MODULE model_mod
