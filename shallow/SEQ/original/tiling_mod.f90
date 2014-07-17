MODULE tiling_mod

  public

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

   ! Extent of this 1D array will be ntiles which, by default will
   ! be set to the no. of threads
   INTEGER, SAVE                                    :: ntiles
   INTEGER, SAVE                                    :: ntilex, ntiley
   TYPE(tile_type), ALLOCATABLE, SAVE, DIMENSION(:) :: tile
   INTEGER, SAVE                                    :: max_tile_width
   INTEGER, SAVE                                    :: max_tile_height

CONTAINS

  SUBROUTINE tiling_init(ntx, nty, jpi, jpj)
!$   USE omp_lib,            ONLY: omp_get_nested, omp_set_max_active_levels
    integer, intent(in) :: jpi, jpj
    integer, intent(in) :: ntx, nty
    ! Locals
    INTEGER :: idx, idy, ival, jval ! For tile extent calculation
    INTEGER :: ierr, nwidth
    INTEGER :: ji,jj, ith
    INTEGER :: nthreads       ! No. of OpenMP threads being used
    INTEGER :: jover, junder, idytmp
    INTEGER :: iover, iunder, idxtmp
    LOGICAL :: nested_par     ! Whether OpenMP supports nested parallelism
    integer :: jpim1, jpjm1
    ! For doing stats on tile sizes
    INTEGER :: nvects, nvects_sum, nvects_min, nvects_max 
    LOGICAL, PARAMETER :: print_tiles = .TRUE.
    ! Whether or not to automatically compute dimensions of tiling grid
    LOGICAL, PARAMETER :: auto_tile = .FALSE.

    jpim1 = jpi - 1
    jpjm1 = jpj - 1

    ntilex = ntx
    ntiley = nty

    ! Set-up regular grid of tiles

    ! Dimensions of the grid of tiles. Done here to save coding
    ! to parse command line for the moment.
    !ntilex = 2
    !ntiley = 16
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
        IF(jpi > jpj)THEN
           IF( ntilex < ntiley )THEN
              ierr   = ntiley
              ntiley = ntilex
              ntilex = ierr
           END IF
        ELSE
           ! jpj >= jpi so want nthready >= nthreadx
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
     ! jpi = (ntilex-2)*(idx-2) + 2*(idx-1)
     ! Rearranging this gives the following expressions for idx and idy:
     idx = (jpi + 6)/ntilex + 2
     idy = (jpj + 6)/ntiley + 2

     ! Integer arithmetic means that ntiley tiles of height idy might
     ! actually span a height greater or less than jpj. If so, we try and
     ! reduce the height of each row by just one until we've accounted
     ! for the <jover> extra rows.
     nwidth = (ntiley-2)*(idy-2) + 2*(idy-1)
     IF(nwidth > jpj)THEN
        jover  = nwidth - jpj
        junder = 0
     ELSE IF(nwidth < jpj)THEN
        jover  = 0
        junder = jpj - nwidth
     ELSE
        jover  = 0
        junder = 0
     END IF
     ! Ditto for x dimension
     nwidth = (ntilex-2)*(idx-2) + 2*(idx-1)
     IF(nwidth > jpi)THEN
        iover  = nwidth - jpi
        iunder = 0
     ELSE IF(nwidth < jpi)THEN
        iover  = 0
        iunder = jpi - nwidth
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
              tile(ith)%iend = jpi
              tile(ith)%iendm1 = jpim1
           ELSE
              tile(ith)%iend = MIN(ival + idxtmp - 1, jpi)
              tile(ith)%iendm1 = tile(ith)%iend - 1
           END IF

           tile(ith)%jbeg = jval
           tile(ith)%jbegp1 = tile(ith)%jbeg + 1

           IF(jj == ntiley)THEN
              tile(ith)%jend = jpj
              tile(ith)%jendm1 = jpjm1
           ELSE
              tile(ith)%jend = MIN(jval + idytmp - 1, jpj)
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
     WRITE(*,"(/'Mean tile size = ',F6.1,' elements or ',F7.1,' KB')") &
                                   REAL(nvects_sum)/REAL(ntiles), &
                                   REAL(8*nvects_sum)/REAL(ntiles*1024)
     WRITE(*,"('Min,max tile size = ',I4,',',I4)") nvects_min,nvects_max
     WRITE(*,"('Tile load imbalance (%) =',F5.1)") &
                              100.0*(nvects_max-nvects_min)/REAL(nvects_min)
     WRITE (*,"('Max tile dims are ',I4,'x',I4/)") max_tile_width, max_tile_height

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

  END SUBROUTINE tiling_init

END MODULE tiling_mod

!!$
!!$   SUBROUTINE tra_adv_tvd ( kt, cdtype, p2dt, pun, pvn, pwn,      &
!!$      &                                       ptb, ptn, pta, kjpt, nloops )
!!$      !!----------------------------------------------------------------------
!!$      !!                  ***  ROUTINE tra_adv_tvd  ***
!!$      !! 
!!$      !! **  Purpose :   Compute the now trend due to total advection of 
!!$      !!       tracers and add it to the general trend of tracer equations
!!$      !!
!!$      !! **  Method  :   TVD scheme, i.e. 2nd order centered scheme with
!!$      !!       corrected flux (monotonic correction)
!!$      !!       note: - this advection scheme needs a leap-frog time scheme
!!$      !!
!!$      !! ** Action : - update (pta) with the now advective tracer trends
!!$      !!             - save the trends 
!!$      !!----------------------------------------------------------------------
!!$!$    USE omp_lib,          ONLY: omp_get_thread_num
!!$      USE TestHarnessTools, ONLY: timer_start, timer_stop
!!$      !!
!!$      INTEGER                              , INTENT(in   ) ::   kt              ! ocean time-step index
!!$      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
!!$      INTEGER                              , INTENT(in   ) ::   kjpt            ! number of tracers
!!$      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt            ! vertical profile of tracer time-step
!!$      REAL(wp), DIMENSION(jpi,jpj,jpk     ), INTENT(in   ) ::   pun, pvn, pwn   ! 3 ocean velocity components
!!$      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb, ptn        ! before and now tracer fields
!!$      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta             ! tracer trend 
!!$      INTEGER,                               INTENT(in) :: nloops ! For timing
!!$      !!
!!$      INTEGER  ::   ji, jj, jk, jn, it       ! dummy loop indices  
!!$      REAL(wp) ::   zbtr, ztra        ! local scalar
!!$      REAL(wp) ::   zfp_wk   !   -      -
!!$      REAL(wp) ::   zfm_wk   !   -      -
!!$
!!$      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE ::   ztrdx, ztrdy, ztrdz
!!$      ! What were work arrays
!!$      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE, SAVE :: zwx,zwy,zwi,zwz
!!$      INTEGER :: tid
!!$      INTEGER :: itime, itime1 ! tags for kernel timer
!!$      REAL(wp) :: azwx, azwy, azwz
!!$
!!$      !!----------------------------------------------------------------------
!!$
!!$
!!$      IF( .not. ALLOCATED(zwx) )  THEN
!!$         IF(lwp) WRITE(numout,*)
!!$         IF(lwp) WRITE(numout,*) 'tra_adv_tvd : TVD advection scheme on ', cdtype
!!$         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
!!$         !
!!$         l_trd = .FALSE.
!!$         ! What were work arrays must now be explicitly allocated
!!$         ALLOCATE( zwx(jpi,jpj,jpk),       zwi(jpi,jpj,jpk),       &
!!$                   zwy(jpi,jpj,jpk),       zwz(jpi,jpj,jpk))
!!$      ENDIF
!!$      !
!!$      !
!!$      IF( l_trd )  THEN
!!$        ALLOCATE( ztrdx(jpi,jpj,jpk) )      ;      ztrdx(:,:,:) = 0.e0_wp
!!$        ALLOCATE( ztrdy(jpi,jpj,jpk) )      ;      ztrdy(:,:,:) = 0.e0_wp
!!$        ALLOCATE( ztrdz(jpi,jpj,jpk) )      ;      ztrdz(:,:,:) = 0.e0_wp
!!$      END IF
!!$      !
!!$!$OMP PARALLEL default(none), shared(ntiles,nloops,l_trd,jpi,jpj,jpk, &
!!$!$OMP          jpkm1,jpjm1,jpim1,kjpt,zwx,zwy,zwz,zwi,     &
!!$!$OMP          pta,ptb,ptn,pun,pvn,pwn,ztrdx,ztrdy,ztrdz,  &
!!$!$OMP          p2dt,tmask,e1t,e2t,e3t,tile), &
!!$!$OMP          private(it,jk,jj,ji,itime,itime1,tid,zfm_wk,zfp_wk, &
!!$!$OMP                  zbtr,ztra,azwx,azwy,azwz)
!!$
!!$      timing: DO it = 1, nloops, 1
!!$
!!$         CALL timer_start('kernel', itime)
!!$
!!$!$OMP DO SCHEDULE(RUNTIME)
!!$!, COLLAPSE(2)
!!$         DO tid = 1, ntiles, 1
!!$         DO jk = 1, jpkm1
!!$            DO jj = tile(tid)%jbegp1, tile(tid)%jendm1 ! 1, jpjm1
!!$               DO ji = tile(tid)%ibegp1, tile(tid)%iendm1 ! 1, jpim1
!!$                  zwi(ji,jj,jk) = 0.e0_wp
!!$               END DO
!!$            END DO
!!$         END DO
!!$         END DO
!!$!$OMP END DO
!!$         !
!!$         !                                                       ! ===========
!!$         DO jn = 1, kjpt                                         ! tracer loop
!!$            !                                                    ! ===========
!!$
!!$            CALL timer_start('tvd_upstream',itime1)
!!$
!!$!$OMP DO SCHEDULE(RUNTIME)
!!$            DO tid = 1, ntiles, 1
!!$
!!$               ! upstream tracer flux in the k direction
!!$               ! Surface value
!!$               IF( lk_vvl ) THEN   
!!$                  DO jj = tile(tid)%jbegp1, tile(tid)%jendm1 ! 1, jpj
!!$                     DO ji = tile(tid)%ibegp1, tile(tid)%iendm1 ! 1, jpi
!!$                        !zwz(ji,jj, 1) = 0.e0                         ! volume variable
!!$
!!$                        zbtr = -1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,1) )
!!$
!!$                        azwx = 0.5 * ( (pun(ji  ,jj,1) + ABS(pun(ji  ,jj,1)) - &
!!$                             (pun(ji-1,jj,1) - ABS(pun(ji-1,jj,1))) )*ptb(ji ,jj,1,jn) + &
!!$                             (pun(ji  ,jj,1) - ABS(pun(ji  ,jj,1)) )*ptb(ji+1,jj,1,jn) - &
!!$                             (pun(ji-1,jj,1) + ABS(pun(ji-1,jj,1)) )*ptb(ji-1,jj,1,jn) )
!!$
!!$                        azwy = 0.5 * ( (pvn(ji,jj  ,1) + ABS(pvn(ji,jj  ,1)) - &
!!$                                    (pvn(ji,jj-1,1) - ABS(pvn(ji,jj-1,1))) )*ptb(ji  ,jj,1,jn) + &
!!$                                    (pvn(ji,jj  ,1) - ABS(pvn(ji,jj  ,1)) )*ptb(ji,jj+1,1,jn) - &
!!$                                    (pvn(ji,jj-1,1) + ABS(pvn(ji,jj-1,1)) )*ptb(ji,jj-1,1,jn) )
!!$
!!$                        ! = zwz(ji,jj,jk) - zwz(ji ,jj ,jk+1) = - zwz(ji ,jj ,jk+1)
!!$                        ! because zwz(:,:,1)=0
!!$                        azwz = -0.5 * ( (pwn(ji,jj,2) + ABS( pwn(ji,jj,2) )) * ptb(ji,jj,2,jn) + &
!!$                                        (pwn(ji,jj,2) - ABS( pwn(ji,jj,2) )) * ptb(ji,jj,1,jn) )
!!$
!!$                        ! total intermediate advective trends
!!$                        ztra =   zbtr * (  azwx + azwy + azwz )
!!$                        ! update and guess with monotonic sheme
!!$                        pta(ji,jj,1,jn) =   pta(ji,jj,1,jn)         + ztra
!!$                        zwi(ji,jj,1)= ( ptb(ji,jj,1,jn) + p2dt(1) * ztra ) * tmask(ji,jj,1)
!!$                        zwi(ji,jj,jpk) = 0.e0 ! Moved from earlier init loop
!!$
!!$                     END DO
!!$                  END DO
!!$               ELSE                
!!$                  DO jj = tile(tid)%jbegp1, tile(tid)%jendm1 ! 1, jpj
!!$                     DO ji = tile(tid)%ibegp1, tile(tid)%iendm1 ! 1, jpi
!!$                        !zwz(ji,jj, 1) = pwn(ji,jj,1) * ptb(ji,jj,1,jn)   ! linear free surface 
!!$
!!$                        zbtr = -1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,1) )
!!$
!!$                        azwx = 0.5 * ( (pun(ji  ,jj,1) + ABS(pun(ji  ,jj,1)) - &
!!$                             (pun(ji-1,jj,1) - ABS(pun(ji-1,jj,1))) )*ptb(ji ,jj,1,jn) + &
!!$                             (pun(ji  ,jj,1) - ABS(pun(ji  ,jj,1)) )*ptb(ji+1,jj,1,jn) - &
!!$                             (pun(ji-1,jj,1) + ABS(pun(ji-1,jj,1)) )*ptb(ji-1,jj,1,jn) )
!!$
!!$                        azwy = 0.5 * ( (pvn(ji,jj  ,1) + ABS(pvn(ji,jj  ,1)) - &
!!$                                    (pvn(ji,jj-1,1) - ABS(pvn(ji,jj-1,1))) )*ptb(ji  ,jj,1,jn) + &
!!$                                    (pvn(ji,jj  ,1) - ABS(pvn(ji,jj  ,1)) )*ptb(ji,jj+1,1,jn) - &
!!$                                    (pvn(ji,jj-1,1) + ABS(pvn(ji,jj-1,1)) )*ptb(ji,jj-1,1,jn) )
!!$
!!$                        ! = zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) 
!!$                        ! = pwn(ji,jj,1) * ptb(ji,jj,1,jn) - zwz(ji  ,jj  ,2)
!!$                        azwz = pwn(ji,jj,1) * ptb(ji,jj,1,jn) - &
!!$                               0.5 * ( (pwn(ji,jj,2) + ABS( pwn(ji,jj,2) )) * ptb(ji,jj,2,jn) + &
!!$                                       (pwn(ji,jj,2) - ABS( pwn(ji,jj,2) )) * ptb(ji,jj,1,jn) )
!!$
!!$                        ! total intermediate advective trends
!!$                        ztra =   zbtr * (  azwx + azwy + azwz )
!!$                        ! update and guess with monotonic sheme
!!$                        pta(ji,jj,1,jn) =   pta(ji,jj,1,jn)         + ztra
!!$                        zwi(ji,jj,1)= ( ptb(ji,jj,1,jn) + p2dt(1) * ztra ) * tmask(ji,jj,1)
!!$                        zwi(ji,jj,jpk) = 0.e0 ! Moved from earlier init loop
!!$
!!$                     END DO
!!$                  END DO
!!$               ENDIF
!!$
!!$               ! Interior value
!!$               DO jk = 2, jpkm1
!!$                  DO jj = tile(tid)%jbegp1, tile(tid)%jendm1 ! 2, jpjm1
!!$                     DO ji = tile(tid)%ibegp1, tile(tid)%iendm1 ! 2, jpim1
!!$
!!$                     zbtr = -1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
!!$
!!$                     azwx = 0.5 * ( (pun(ji  ,jj,jk) + ABS(pun(ji  ,jj,jk)) - &
!!$                                     pun(ji-1,jj,jk) + ABS(pun(ji-1,jj,jk)) )*ptb(ji  ,jj,jk,jn) + &
!!$                                    (pun(ji  ,jj,jk) - ABS(pun(ji  ,jj,jk)) )*ptb(ji+1,jj,jk,jn) - &
!!$                                    (pun(ji-1,jj,jk) + ABS(pun(ji-1,jj,jk)) )*ptb(ji-1,jj,jk,jn) )
!!$
!!$                     azwy = 0.5 * ( (pvn(ji,jj  ,jk) + ABS(pvn(ji,jj  ,jk)) - &
!!$                                     pvn(ji,jj-1,jk) + ABS(pvn(ji,jj-1,jk)) )*ptb(ji,jj  ,jk,jn) + &
!!$                                    (pvn(ji,jj  ,jk) - ABS(pvn(ji,jj  ,jk)) )*ptb(ji,jj+1,jk,jn) - &
!!$                                    (pvn(ji,jj-1,jk) + ABS(pvn(ji,jj-1,jk)) )*ptb(ji,jj-1,jk,jn) )
!!$
!!$                     azwz = 0.5 * ( (pwn(ji,jj,jk  ) + ABS(pwn(ji,jj,jk  )) - &
!!$                                     pwn(ji,jj,jk+1) + ABS(pwn(ji,jj,jk+1))) * ptb(ji,jj,jk,jn) + &
!!$                                    (pwn(ji,jj,jk  ) - ABS(pwn(ji,jj,jk  ))) * ptb(ji,jj,jk-1,jn) - &
!!$                                    (pwn(ji,jj,jk+1) + ABS(pwn(ji,jj,jk+1))) * ptb(ji,jj,jk+1,jn) )
!!$
!!$                     ! total intermediate advective trends
!!$                     ztra =   zbtr * (  azwx & !zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk )   &
!!$                                      + azwy & !zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk )   &
!!$                                      + azwz ) !zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) )
!!$                     ! update and guess with monotonic sheme
!!$                     pta(ji,jj,jk,jn) =   pta(ji,jj,jk,jn)         + ztra
!!$                     zwi(ji,jj,jk)= ( ptb(ji,jj,jk,jn) + p2dt(jk) * ztra ) * tmask(ji,jj,jk)
!!$                  END DO
!!$               END DO
!!$            END DO
!!$
!!$
!!$            ! 3. antidiffusive flux : high order minus low order
!!$            ! --------------------------------------------------
!!$            ! antidiffusive flux on i and j
!!$
!!$            DO jk = 1, jpkm1
!!$               DO jj = tile(tid)%jbegp1, tile(tid)%jendm1 ! 1, jpjm1
!!$                  DO ji = tile(tid)%ibegp1, tile(tid)%iendm1 ! 1, fs_jpim1   ! vector opt.
!!$                     zwx(ji,jj,jk) = 0.5 * ( pun(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji+1,jj,jk,jn) ) &
!!$                                          - (pun(ji,jj,jk) + ABS( pun(ji,jj,jk) ))*ptb(ji,jj,jk,jn)    &
!!$                                          - (pun(ji,jj,jk) - ABS( pun(ji,jj,jk) ))*ptb(ji+1,jj,jk,jn) )
!!$                     zwy(ji,jj,jk) = 0.5 * ( pvn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj+1,jk,jn) ) &
!!$                                          - (pvn(ji,jj,jk) + ABS( pvn(ji,jj,jk) ))*ptb(ji,jj,jk,jn)    &
!!$                                          - (pvn(ji,jj,jk) - ABS( pvn(ji,jj,jk) ))*ptb(ji,jj+1,jk,jn) )
!!$                  END DO
!!$               END DO
!!$            END DO
!!$
!!$            DO jj = tile(tid)%jbegp1, tile(tid)%jendm1 ! 1, jpj
!!$               DO ji = tile(tid)%ibegp1, tile(tid)%iendm1 ! 1, jpi
!!$                  zwz(ji,jj,1) = 0.e0         ! Surface value
!!$                  zwz(ji,jj,jpk) = 0.e0       ! Bottom value   
!!$                  zwx(ji,jj,jpk) = 0.e0
!!$                  zwy(ji,jj,jpk) = 0.e0    
!!$               END DO
!!$            END DO
!!$
!!$            DO jk = 2, jpkm1          ! Interior value
!!$               DO jj = tile(tid)%jbegp1, tile(tid)%jendm1 ! 1, jpj
!!$                  DO ji = tile(tid)%ibegp1, tile(tid)%iendm1 ! 1, jpi
!!$                     zwz(ji,jj,jk) = 0.5 * pwn(ji,jj,jk) * ( ptn(ji,jj,jk,jn) + ptn(ji,jj,jk-1,jn) ) &
!!$!                                         - zwz(ji,jj,jk)
!!$                                   - 0.5 * ( (pwn(ji,jj,jk) + ABS(pwn(ji,jj,jk)))*ptb(ji,jj,jk,jn) + &
!!$                                             (pwn(ji,jj,jk) - ABS(pwn(ji,jj,jk)))*ptb(ji,jj,jk-1,jn) )
!!$                  END DO
!!$               END DO
!!$            END DO
!!$
!!$         END DO ! tid
!!$!$OMP END DO NOWAIT
!!$            CALL timer_stop(itime1)
!!$!$OMP BARRIER
!!$
!!$            ! Only the master thread calls the MPI library
!!$            CALL timer_start('tvd_lbc',itime1)
!!$!$OMP MASTER
!!$            CALL lbc_lnk( zwi, 'T', 1._wp )  
!!$            CALL lbc_lnk( zwx, 'U', -1._wp )   
!!$            CALL lbc_lnk( zwy, 'V', -1._wp )         ! Lateral bondary conditions
!!$            CALL lbc_lnk( zwz, 'W',  1._wp )
!!$!$OMP END MASTER
!!$            CALL timer_stop(itime1)
!!$!$OMP BARRIER
!!$
!!$            ! 4. monotonicity algorithm
!!$            ! -------------------------
!!$            CALL nonosc( ptb(1,1,1,jn), zwx, zwy, zwz, zwi, p2dt )
!!$
!!$            ! 5. final trend with corrected fluxes
!!$            ! ------------------------------------
!!$            CALL timer_start('tvd_finaltr',itime1)
!!$!$OMP DO SCHEDULE(RUNTIME), COLLAPSE(2)
!!$            DO tid = 1, ntiles, 1
!!$            DO jk = 1, jpkm1
!!$               DO jj = tile(tid)%jbegp1, tile(tid)%jendm1 ! 2, jpjm1
!!$                  DO ji = tile(tid)%ibegp1, tile(tid)%iendm1 ! fs_2, fs_jpim1   ! vector opt.  
!!$                     zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
!!$                     ! total advective trends
!!$                     ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
!!$                          &           + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )   &
!!$                          &           + zwz(ji,jj,jk) - zwz(ji  ,jj  ,jk+1) )
!!$                     ! add them to the general tracer trends
!!$                     pta(ji,jj,jk,jn) = pta(ji,jj,jk,jn) + ztra
!!$                  END DO
!!$               END DO
!!$            END DO
!!$            END DO
!!$!$OMP END DO NOWAIT
!!$            CALL timer_stop(itime1)
!!$!$OMP BARRIER
!!$
!!$         !
!!$         ENDDO ! loop over tracers
!!$
!!$         CALL timer_stop(itime)
!!$
!!$      ENDDO timing
!!$
!!$!$OMP END PARALLEL
!!$
!!$      !
!!$      IF( l_trd )  THEN
!!$        DEALLOCATE( ztrdx )   ;   DEALLOCATE( ztrdy )  ;  DEALLOCATE( ztrdz )
!!$      END IF
!!$      !
!!$   END SUBROUTINE tra_adv_tvd
