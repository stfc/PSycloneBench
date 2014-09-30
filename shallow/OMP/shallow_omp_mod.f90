module shallow_omp_mod

  private

  ! For coarse-grained OpenMP tiling
  TYPE :: tile_type
     INTEGER :: ibeg
     INTEGER :: iend
     INTEGER :: jbeg
     INTEGER :: jend
     ! For loop limits that are 2 instead of 1 and jp{i,j}-1 instead
     ! of jp{i,j}
     INTEGER :: ibegp1
     INTEGER :: iendm1
     INTEGER :: jbegp1
     INTEGER :: jendm1
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

contains

  !================================================

  SUBROUTINE openmp_grid_init(nx, ny)
    use topology_mod, only: M, N
    integer, intent(in), optional :: nx, ny
    INTEGER :: idx, idy, ival, jval ! For tile extent calculation
    INTEGER :: ierr, nwidth
    INTEGER :: ji,jj, ith
    INTEGER :: nthreads       ! No. of OpenMP threads being used
    INTEGER :: jover, junder, idytmp
    INTEGER :: iover, iunder, idxtmp
    LOGICAL :: nested_par     ! Whether OpenMP supports nested parallelism
    ! For doing stats on tile sizes
    INTEGER :: nvects, nvects_sum, nvects_min, nvects_max 
    LOGICAL, PARAMETER :: print_tiles = .TRUE.

    ! Set-up regular grid of tiles
    auto_tile = .TRUE.

    ! Dimensions of the grid of tiles. 
    if(present(nx))then
       ntilex = nx
       auto_tile = .FALSE.
    else
       ntilex = 1
    endif
    if(present(ny))then
       ntiley = ny
       auto_tile = .FALSE.
    else
       ntiley = 1
    end if

    ntiles = ntilex * ntiley

     nthreads = 1
!$   nthreads = omp_get_max_threads()
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
        IF(M > N)THEN
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
     ! M = (ntilex-2)*(idx-2) + 2*(idx-1)
     ! Rearranging this gives the following expressions for idx and idy:
     idx = (M + 6)/ntilex + 2
     idy = (N + 6)/ntiley + 2

     ! Integer arithmetic means that ntiley tiles of height idy might
     ! actually span a height greater or less than N. If so, we try and
     ! reduce the height of each row by just one until we've accounted
     ! for the <jover> extra rows.
     nwidth = (ntiley-2)*(idy-2) + 2*(idy-1)
     IF(nwidth > N)THEN
        jover  = nwidth - N
        junder = 0
     ELSE IF(nwidth < N)THEN
        jover  = 0
        junder = N - nwidth
     ELSE
        jover  = 0
        junder = 0
     END IF
     ! Ditto for x dimension
     nwidth = (ntilex-2)*(idx-2) + 2*(idx-1)
     IF(nwidth > M)THEN
        iover  = nwidth - M
        iunder = 0
     ELSE IF(nwidth < M)THEN
        iover  = 0
        iunder = M - nwidth
     ELSE
        iover  = 0
        iunder = 0
     END IF

     ! For AVX instructions, I think we want MOD(idx,4) == 0
     !idx = idx + (4 - MOD(idx,4))

     WRITE(*,"('Tile width = ',I4,', tile height = ',I4)") idx, idy
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
           idytmp = idy - 1
           jover = jover - 1
        ELSE IF(junder > 0)THEN
           idytmp = idy + 1
           junder = junder - 1
        ELSE
           idytmp = idy
        END IF

        ival = 1

        DO ji = 1, ntilex, 1
         
           ! If necessary, correct the width of this tile column
           IF(iover > 0)THEN
              idxtmp = idx - 1
              iover = iover - 1
           ELSE IF(iunder > 0)THEN
              idxtmp = idx + 1
              iunder = iunder - 1
           ELSE
              idxtmp = idx
           END IF

           tile(ith)%ibeg = ival
           tile(ith)%ibegp1 = tile(ith)%ibeg + 1

           IF(ji == ntilex)THEN
              tile(ith)%iend = M
              tile(ith)%iendm1 = M-1
           ELSE
              tile(ith)%iend = MIN(ival + idxtmp - 1, M)
              tile(ith)%iendm1 = tile(ith)%iend - 1
           END IF

           tile(ith)%jbeg = jval
           tile(ith)%jbegp1 = tile(ith)%jbeg + 1

           IF(jj == ntiley)THEN
              tile(ith)%jend = N
              tile(ith)%jendm1 = N-1
           ELSE
              tile(ith)%jend = MIN(jval + idytmp - 1, N)
              tile(ith)%jendm1 = tile(ith)%jend - 1
           END IF

           IF(print_tiles)THEN
              WRITE(*,"('tile[',I4,'](',I4,':',I4,')(',I4,':',I4,'), interior:(',I4,':',I4,')(',I4,':',I4,') ')") &
                 ith,                                &
                 tile(ith)%ibeg, tile(ith)%iend,     &
                 tile(ith)%jbeg, tile(ith)%jend,     &
                 tile(ith)%ibegp1, tile(ith)%iendm1, &
                 tile(ith)%jbegp1, tile(ith)%jendm1
           END IF

           ! Collect some data on the distribution of tile sizes for loadbalance info
           nvects = (tile(ith)%iendm1 - tile(ith)%ibegp1 + 1) * &
                    (tile(ith)%jendm1 - tile(ith)%jbegp1 + 1)
           nvects_sum = nvects_sum + nvects
           nvects_min = MIN(nvects_min, nvects)
           nvects_max = MAX(nvects_max, nvects)

           ! For use when allocating tile-'private' work arrays
           max_tile_width  = MAX(max_tile_width, (tile(ith)%iend - tile(ith)%ibeg + 1) )
           max_tile_height = MAX(max_tile_height, (tile(ith)%jend - tile(ith)%jbeg + 1) )

           ival = ival + idxtmp - 2
           ith = ith + 1
        END DO
        jval = jval + idytmp - 2
     END DO

     ! Print tile-size statistics
     WRITE(*,"(/'Mean tile size = ',F6.1,' cols = ',F7.1,' KB')") &
                                   REAL(nvects_sum)/REAL(ntiles), &
                                   REAL(8*nvects_sum)/REAL(ntiles*1024)
     WRITE(*,"('Min,max tile size = ',I4,',',I4)") nvects_min,nvects_max
     WRITE(*,"('Tile load imbalance (%) =',F5.1)") &
                              100.0*(nvects_max-nvects_min)/REAL(nvects_min)
     WRITE (*,"('Max tile dims are ',I4,'x',I4)") max_tile_width, max_tile_height

     ! Check to see whether nested parallelism is supported
     nested_par = .FALSE.
!$   CALL omp_set_nested(.TRUE.)
!$   nested_par = omp_get_nested()

     IF( nested_par )THEN
        WRITE (*,"(/'OpenMP implementation SUPPORTS nested parallelism.'/)")
     ELSE
!$      WRITE (*,"(/'OpenMP implementation DOES NOT support nested parallelism.'/)")        
     END IF

!$   CALL omp_set_max_active_levels(2)

  END SUBROUTINE openmp_grid_init

end module shallow_omp_mod
