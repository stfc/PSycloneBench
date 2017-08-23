module dl_timer_parallel
  !> Module containing all routines that involve calls to MPI routines.
  use mpi
  use dl_timer_constants_mod
  implicit none

contains

  function is_parallel()
    !> Returns .TRUE. to indicate that dl_timer is built with MPI support.
    !! Aborts if MPI_Init() has not yet been called.
    logical :: is_parallel
    integer :: ierr

    call MPI_INITIALIZED(is_parallel, ierr)

    if(.not. is_parallel)then
      write(*, &
           "('TIMING: ERROR: timer_init() must be called after MPI_Init()!')")
      stop
    end if

  end function is_parallel

  !=========================================================================

  function get_rank()
    !> Returns the rank of this process in MPI_COMM_WORLD
    integer :: get_rank
    integer :: ierr
    call MPI_COMM_RANK(MPI_COMM_WORLD, get_rank, ierr)
    return
  end function get_rank

  !=========================================================================

  function num_ranks()
    !> Returns the number of ranks in MPI_COMM_WORLD
    integer :: num_ranks
    integer :: ierr
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_ranks, ierr)
    return
  end function num_ranks

  !=========================================================================

  subroutine calc_dm_timer_stats(nThreads, ntimers, region_names, &
                                 visit_counts, times,             &
                                 max_times, min_times, sum_times, &
                                 sum_counts, npes_active)
    integer,                                  intent(in) :: nThreads, ntimers
    !> The name of each timer on this process. On return holds a list of
    !! all unique active timers aggregated over ranks.
    character(len=LABEL_LEN), dimension(ntimers,nThreads), intent(inout) :: region_names
    !> The total time spent in each region by each thread on this process
    real(r_def), dimension(ntimers,nThreads),    intent(in) :: times
    !> The number of times each thread on this process has visited each
    !! timed region. (Can be zero.)
    integer(i_def64), dimension(ntimers, nThreads), intent(in) :: visit_counts
    !> Max/Min time spent in region for a single visit together with
    !! rank of corresponding PE. Time is in (1,x,y), rank is in (2,x,y).
    real(r_def), dimension(2,ntimers,nThreads), intent(out) :: max_times, &
                                                               min_times
    !> The time spent in each region, summed over all active PEs and visits 
    real(r_def), dimension(ntimers,nThreads),   intent(out) :: sum_times
    !> The no. of times each region is visited, summed over all PEs
    integer(i_def64), dimension(ntimers, nThreads), intent(out) :: sum_counts
    !> On return holds the number of PES that have visited each unique
    !! timed region
    integer, dimension(ntimers), intent(out) :: npes_active
    !-----------------------------------------------------------------------
    ! Locals
    real(r_def), allocatable, dimension(:,:,:) :: times_ranks
    ! The names of all timers on all ranks packed into a 1D array. We can do
    ! this because we know that each label is LABEL_LEN chars long
    character(len=1), allocatable, dimension(:) :: labels_ranks
    ! The names of all timers on the local process packed into a 1D array.
    character(len=1), allocatable, dimension(:) :: labels_merged
    character(len=LABEL_LEN), allocatable, dimension(:,:) :: region_names_by_rank
    !> Array to store gather of the counts from all timers on all ranks
    integer(i_def64), allocatable :: all_counts(:)
    !> Array to store gather of the times from all timers on all ranks
    real(r_def), allocatable :: all_times(:)
    !> The number of uniquely-named timed-regions across all PEs
    integer :: unique_region_count
    !> The name of each of these unique regions
    character(len=LABEL_LEN), allocatable, dimension(:) :: unique_region_labels
    !> unique_timer_map(timer, rank) gives the index of the unique timed
    !! region on PE rank. If that PE does not have the region then we
    !! store a zero.
    integer, allocatable, dimension(:,:) :: unique_region_map
    !> Will store appropriate MPI type for a 64-bit integer
    integer :: int64_mpi_type
    ! No. of times a given PE has visited a region
    integer(i_def64) :: ncounts_on_pe
    logical :: is_unique
    
    integer :: ierr, myrank, nranks
    integer :: j, jt, itimer, irank, index, ichar
    real(r_def) :: time, time_per_visit

    myrank = get_rank()
    nranks = num_ranks()

    allocate(times_ranks(2,ntimers,nThreads),        &
             labels_merged(ntimers*LABEL_LEN),       &
             labels_ranks(ntimers*nranks*LABEL_LEN), &
             region_names_by_rank(nranks,ntimers),   &
             unique_region_labels(ntimers),          &
             unique_region_map(ntimers, nranks),     &
             all_counts(ntimers*nranks), all_times(ntimers*nranks), &
             Stat=ierr)
    if(ierr /= 0)then
       write (*,*) 'TIMING: calc_dm_timer_stats: failed to allocate memory'
       return
    end if
    
    labels_ranks(:) = ""
    region_names_by_rank(:,:) = ""

    ! Pack all the timed-region labels on this rank into one long
    ! array of chars
    DO itimer = 1, ntimers
       index = (itimer-1)*LABEL_LEN
       do ichar = 1, LABEL_LEN
          labels_merged(index+ichar) = region_names(itimer,1)(ichar:ichar)
       end do
    END DO

    ! We must pack the timing data into an array suitable for the
    ! reduction operations
    do jt = 1, nThreads, 1
       do itimer = 1, ntimers !itimerCount(jt)
          times_ranks(1,itimer,jt) = times(itimer,jt)
          times_ranks(2,itimer,jt) = myrank
       end do
    end do

    ! ARPDBG this is just for thread 1 currently
    call MPI_Gather(labels_merged, ntimers*LABEL_LEN, MPI_CHARACTER, &
                    labels_ranks, ntimers*LABEL_LEN, MPI_CHARACTER, 0,   &
                    MPI_COMM_WORLD, ierr)
    ! No. of times each region is visited on each PE
    ! We have to set-up a 64-bit integer type for MPI
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER, 8, int64_mpi_type, ierr)
    call MPI_Gather(visit_counts(:,1), ntimers, int64_mpi_type, &
                    all_counts, ntimers, int64_mpi_type, 0, &
                    MPI_COMM_WORLD, ierr)

    ! The timings themselves
    call MPI_Gather(times(:,1), ntimers, MPI_DOUBLE_PRECISION, &
                    all_times, ntimers, MPI_DOUBLE_PRECISION, 0, &
                    MPI_COMM_WORLD, ierr)

    if(myrank == 0)then

       unique_region_map(:,:) = 0

       ! Unpack the timed-region labels from arrays back into strings
       do irank = 1, nranks
          do itimer = 1, ntimers
             index = ((irank-1)*ntimers + itimer - 1)*LABEL_LEN + 1
             do ichar = 1, LABEL_LEN
                region_names_by_rank(irank,itimer)(ichar:ichar) = labels_ranks(index+ichar-1)
             end do
          end do
       end do

       ! Now we must work out how the different timed regions on different
       ! ranks are related.
       unique_region_count = 0
       do irank = 1, nranks
          do itimer = 1, ntimers
             ! Skip blank labels
             if(len_trim(region_names_by_rank(irank,itimer)) == 0)cycle
             is_unique = .TRUE.
             do j = 1, unique_region_count
                if (region_names_by_rank(irank,itimer) == &
                    unique_region_labels(j)) then
                   is_unique = .FALSE.
                   exit
                end if
             end do
             if (is_unique) then
                ! We haven't seen a region with this name before
                unique_region_count = unique_region_count + 1
                unique_region_labels(unique_region_count) = &
                                          region_names_by_rank(irank,itimer)
                ! Store the index of this region on this PE
                unique_region_map(unique_region_count,irank) = itimer
             else
                ! We have seen this region before. Store its index on this PE.
                unique_region_map(j, irank) = itimer
             end if
          end do
       end do

       ! Loop over the unique timers we've identified and collect the stats
       ! from the corresponding timer on each rank
       min_times(1,:,:) = 1.0E20
       max_times(1,:,:) = -1.0
       sum_times(:,:) = 0.0_r_def
       sum_counts(:,:) = 0_i_def64
       npes_active(:) = 0
       do itimer = 1, unique_region_count
          ! Write back the labels into the original list so that they are
          ! in the same order as the rest of the data
          region_names(itimer,1) = unique_region_labels(itimer)
          do irank = 1, nranks
             ! What index did this timer have on rank irank?
             index = unique_region_map(itimer, irank)
             ! Check whether this rank has used this timed region - if not
             ! then we skip it
             if(index == 0) cycle
             ! It's possible that a process might register a region but then
             ! never visit it so check for that
             ncounts_on_pe = all_counts((irank-1)*ntimers + index)
             if(ncounts_on_pe == 0_i_def64) cycle

             ! Increment the count of PEs for which this timed region is
             ! active
             npes_active(itimer) = npes_active(itimer) + 1
             ! Sum up the total no. of times this region has been visited over
             ! all PEs
             sum_counts(itimer, 1) = sum_counts(itimer, 1) + ncounts_on_pe

             ! Stored time is the sum over all visits to the region on this PE
             time = all_times((irank-1)*ntimers + index)
             time_per_visit = time / REAL(ncounts_on_pe, kind=r_def)

             ! Minimum time spent in this region by any rank
             if(min_times(1,itimer,1) > time_per_visit)then
                min_times(1,itimer,1) = time_per_visit
                min_times(2,itimer,1) = irank-1
             end if
             ! Maximum time spent in this region by any rank
             if(max_times(1,itimer,1) < time_per_visit)then
                max_times(1,itimer,1) = time_per_visit
                max_times(2,itimer,1) = irank-1
             end if
             ! Total time spent in this region summed over ranks
             sum_times(itimer,1) = sum_times(itimer,1) + time
          end do
       end do

    end if ! rank 0

    return
  end subroutine calc_dm_timer_stats

end module dl_timer_parallel
