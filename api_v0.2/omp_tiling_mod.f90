module omp_tiling_mod
  use kind_params_mod
  use region_mod
  private

  ! For coarse-grained OpenMP tiling
  TYPE :: tile_type
     type(region_type) :: internal
     type(region_type) :: whole
     real(wp), dimension(:,:), allocatable :: data
  END TYPE tile_type

  INTEGER, SAVE                                    :: ntiles
  INTEGER, SAVE                                    :: ntilex, ntiley
  !> Whether or not to automatically compute dimensions of tiling grid
  LOGICAL, SAVE                                    :: auto_tile
  ! Extent of this 1D array will be ntiles which will default to
  ! the no. of threads
  TYPE(tile_type), ALLOCATABLE, SAVE, DIMENSION(:) :: tile
  INTEGER, SAVE                                    :: max_tile_width
  INTEGER, SAVE                                    :: max_tile_height

  public openmp_grid_init
  public tile_type, ntiles, tile

contains

  !================================================

  SUBROUTINE openmp_grid_init(xlen, ylen, nx_arg, ny_arg)
    use omp_lib, only: omp_get_max_threads
    implicit none
    !> Dimensions of the model mesh
    integer, intent(in) :: xlen, ylen
    !> Optional specification of the dimensions of the tiling grid
    integer, intent(in), optional :: nx_arg, ny_arg
    integer :: nx, ny
    INTEGER :: idx, idy, ival, jval ! For tile extent calculation
    integer :: internal_width, internal_height
    INTEGER :: ierr, nwidth
    INTEGER :: ji,jj, ith
    INTEGER :: nthreads       ! No. of OpenMP threads being used
    INTEGER :: jover, junder, idytmp, jover_per_row, junder_per_row
    INTEGER :: iover, iunder, idxtmp, iover_per_col, iunder_per_col
    LOGICAL :: nested_par     ! Whether OpenMP supports nested parallelism
    ! For doing stats on tile sizes
    INTEGER :: nvects, nvects_sum, nvects_min, nvects_max 
    LOGICAL, PARAMETER :: print_tiles = .TRUE.

    ! Set-up regular grid of tiles
    auto_tile = .TRUE.

    ! Dimensions of the grid of tiles. 
    if(get_grid_dims(nx,ny) )then
       ntilex = nx
       ntiley = ny
       auto_tile = .FALSE.
    else if( present(nx_arg) .and. present(ny_arg) )then
       ntilex = nx_arg
       ntiley = ny_arg
       auto_tile = .FALSE.
    else
       ntilex = 1
       ntiley = 1
    end if

    
    ntiles = ntilex * ntiley

    nthreads = 1
!$  nthreads = omp_get_max_threads()
    WRITE (*,"(/'Have ',I3,' OpenMP threads available.')") nthreads

    ! If we've not manually specified a grid of tiles then use the no. of
    ! threads
    IF(ntiles == 1 .AND. auto_tile)ntiles = nthreads

    ALLOCATE(tile(ntiles), Stat=ierr)
    IF(ierr /= 0 )THEN
       STOP 'Harness: ERROR: failed to allocate tiling structures'
    END IF

    IF(auto_tile)THEN

       ntilex = INT( SQRT(REAL(ntiles)) )
       DO WHILE(MOD(ntiles,ntilex) /= 0)
          ntilex = ntilex - 1
       END DO
       ntiley = ntiles/ntilex

       ! Match longest dimension of MPI domain to longest dimension of 
       ! thread grid
       IF(xlen > ylen)THEN
          IF( ntilex < ntiley )THEN
             ierr   = ntiley
             ntiley = ntilex
             ntilex = ierr
          END IF
       ELSE
          ! N >= M so want nthready >= nthreadx
          IF( ntiley < ntilex )THEN
             ierr   = ntiley
             ntiley = ntilex
             ntilex = ierr
          END IF
       END IF

    END IF ! automatic determination of tiling grid

    WRITE (*,"('OpenMP thread tiling using grid of ',I3,'x',I3)") ntilex,ntiley

    ! Tiles at left and right of domain only have single
    ! overlap. Every other tile has two overlaps. So: 
    ! xlen = (ntilex-2)*(idx-2) + 2*(idx-1)
    !      = ntilex.idx - 2.ntilex - 2.idx + 4 + 2.idx - 2
    !      = ntilex.idx + 2 - 2.ntilex
    !=> idx = (xlen - 2 + 2.ntilex)/ntilex
    ! where idx is the whole width of a tile.
    !idx = NINT(REAL(xlen - 2 + 2*ntilex)/REAL(ntilex))
    !idy = NINT(REAL(ylen - 2 + 2*ntiley)/REAL(ntiley))
    ! Alternatively, if we think about the internal regions of the tiles,
    ! then they should share the domain between them:
    internal_width = NINT(REAL(xlen) / REAL(ntilex))
    internal_height = NINT(REAL(ylen) / REAL(ntiley))

    ! Integer arithmetic means that ntiley tiles of height idy might
    ! actually span a height greater or less than N. If so, we try and
    ! reduce the height of each row by just one until we've accounted
    ! for the <jover> extra rows.
    !nwidth = (ntiley-2)*(idy-2) + 2*(idy-1)
    nwidth = ntiley * internal_height
    IF(nwidth > ylen)THEN
       jover  = nwidth - ylen
       junder = 0
    ELSE IF(nwidth < ylen)THEN
       jover  = 0
       junder = ylen - nwidth
    ELSE
       jover  = 0
       junder = 0
    END IF
    ! Ditto for x dimension
    !nwidth = (ntilex-2)*(idx-2) + 2*(idx-1)
    nwidth = ntilex * internal_width
    IF(nwidth > xlen)THEN
       iover  = nwidth - xlen
       iunder = 0
    ELSE IF(nwidth < xlen)THEN
       iover  = 0
       iunder = xlen - nwidth
    ELSE
       iover  = 0
       iunder = 0
    END IF

    ! For AVX (256-bit vector) instructions, I think we want
    ! MOD(idx,4) == 0 idx = idx + (4 - MOD(idx,4))

    WRITE(*,"('Tile width = ',I4,', tile height = ',I4)") &
         internal_width, internal_height
    WRITE(*,"('iover = ',I3,', iunder = ',I3)") iover, iunder
    WRITE(*,"('jover = ',I3,', junder = ',I3)") jover, junder

    ith = 1
    jval = 1

    nvects_max = 0
    nvects_min = 1000000
    nvects_sum = 0
    max_tile_width  = 0
    max_tile_height = 0

    IF(print_tiles)WRITE(*,"(/'Tile dimensions:')")

    DO jj = 1, ntiley, 1

       ! If necessary, correct the height of this tile row
       IF(jover > 0)THEN
          idytmp = internal_height - 1
          jover = jover - 1
       ELSE IF(junder > 0)THEN
          idytmp = internal_height + 1
          junder = junder - 1
       ELSE
          idytmp = internal_height
       END IF

       ival = 1

       DO ji = 1, ntilex, 1
         
          ! If necessary, correct the width of this tile column
          IF(iover > 0)THEN
             idxtmp = internal_width - 1
             iover = iover - 1
          ELSE IF(iunder > 0)THEN
             idxtmp = internal_width + 1
             iunder = iunder - 1
          ELSE
             idxtmp = internal_width
          END IF

          if(ji == 1)then
             tile(ith)%whole%xstart    = ival
             tile(ith)%internal%xstart = ival
          else
             tile(ith)%whole%xstart    = ival - 1
             tile(ith)%internal%xstart = ival
          end if
          
          IF(ji == ntilex)THEN
             tile(ith)%internal%xstop = xlen
             tile(ith)%whole%xstop = tile(ith)%internal%xstop
          ELSE
             tile(ith)%internal%xstop =  MIN(xlen-1, &
                                      tile(ith)%internal%xstart + idxtmp - 1)
             tile(ith)%whole%xstop = tile(ith)%internal%xstop + 1
          END IF
          
          if(jj == 1)then
             tile(ith)%whole%ystart    = jval
             tile(ith)%internal%ystart = jval
          else
             tile(ith)%whole%ystart    = jval - 1
             tile(ith)%internal%ystart = jval
          end if

          IF(jj /= ntiley)THEN
             tile(ith)%internal%ystop =  MIN(tile(ith)%internal%ystart+idytmp-1, &
                                             ylen-1)
             tile(ith)%whole%ystop = tile(ith)%internal%ystop + 1
          ELSE
             tile(ith)%internal%ystop = ylen
             tile(ith)%whole%ystop = tile(ith)%internal%ystop
          END IF

          IF(print_tiles)THEN
             WRITE(*,"('tile[',I4,'](',I4,':',I4,')(',I4,':',I4,'), interior:(',I4,':',I4,')(',I4,':',I4,') ')")                                       &
                  ith,                                                 &
                  tile(ith)%whole%xstart, tile(ith)%whole%xstop,       &
                  tile(ith)%whole%ystart, tile(ith)%whole%ystop,       &
                  tile(ith)%internal%xstart, tile(ith)%internal%xstop, &
                  tile(ith)%internal%ystart, tile(ith)%internal%ystop
          END IF

          ! Collect some data on the distribution of tile sizes for 
          ! loadbalance info
          nvects = (tile(ith)%internal%xstop - tile(ith)%internal%xstart + 1) &
                  * (tile(ith)%internal%ystop - tile(ith)%internal%ystart + 1)
          nvects_sum = nvects_sum + nvects
          nvects_min = MIN(nvects_min, nvects)
          nvects_max = MAX(nvects_max, nvects)

          ! For use when allocating tile-'private' work arrays
          max_tile_width  = MAX(max_tile_width, &
                          (tile(ith)%whole%xstop - tile(ith)%whole%xstart + 1) )
          max_tile_height = MAX(max_tile_height, &
                          (tile(ith)%whole%ystop - tile(ith)%whole%ystart + 1) )

          ival = tile(ith)%whole%xstop
          ith = ith + 1
       END DO
       jval = tile(ith-1)%whole%ystop
    END DO

    ! Print tile-size statistics
    WRITE(*,"(/'Mean tile size = ',F8.1,' pts = ',F7.1,' KB')") &
                                 REAL(nvects_sum)/REAL(ntiles), &
                                 REAL(8*nvects_sum)/REAL(ntiles*1024)
    WRITE(*,"('Min,max tile size (pts) = ',I6,',',I6)") nvects_min,nvects_max
    WRITE(*,"('Tile load imbalance (%) =',F6.2)") &
                           100.0*(nvects_max-nvects_min)/REAL(nvects_min)
    WRITE (*,"('Max tile dims are ',I4,'x',I4/)") max_tile_width, &
                                                  max_tile_height

  END SUBROUTINE openmp_grid_init

  !==============================================

  function get_grid_dims(nx, ny) result(success)
    implicit none
    integer, intent(inout) :: nx, ny
    logical :: success
    character(len=20) :: lstr
    integer :: idx, ierr

    success = .FALSE.

    call get_environment_variable(NAME='GOCEAN_OMP_GRID', VALUE=lstr, &
                                  STATUS=ierr)

    if(ierr /= 0)return

    ! We expect the string to have the format 'AxB' where A and B are
    ! integers.
    idx = index(lstr, 'x')
    if(idx == 0)then
       write (*,"(/'shallow_omp_mod::get_grid_dims: failed to parse ' &
                 &  'GOCEAN_OMP_GRID string: ',(A))") TRIM(lstr)
       write (*,"('   -  will use defaults for dimensions of tiling grid')")
       return
    endif

    read(lstr(1:idx-1),*,iostat=ierr) nx
    if(ierr /= 0)return

    read(lstr(idx+1:),*,iostat=ierr) ny
    if(ierr == 0)success = .TRUE.

  end function get_grid_dims

  !==============================================

end module omp_tiling_mod
