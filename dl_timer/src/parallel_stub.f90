module dl_timer_parallel
  !> Module containing empty stubs of routines that, in an MPI
  !! build, would make calls to MPI routines. This allows us
  !! to compile without requiring MPI to be installed.
  !! When performing an MPI build we use parallel_mpi.f90
  !! instead.
  use dl_timer_constants_mod
  implicit none

contains

  function is_parallel()
    logical :: is_parallel
    is_parallel = .FALSE.
  end function is_parallel

  !=========================================================================

  function get_rank()
    !> No MPI support so we only have one process - rank 0
    integer :: get_rank
    get_rank = 0
  end function get_rank

  !=========================================================================

  function num_ranks()
    !> No MPI support so we only have a single process
    integer :: num_ranks
    num_ranks = 1
  end function num_ranks

  !=========================================================================

  subroutine calc_dm_timer_stats(nThreads, ntimers, region_names, &
                                 visit_counts, times,             &
                                 max_times, min_times, sum_times, &
                                 sum_counts, npes_active)
    !> Stub interface. Routine itself does nothing in absence of MPI.
    character(len=LABEL_LEN), intent(inout) :: region_names(ntimers,nThreads)
    integer,                     intent(in) :: nThreads, ntimers
    integer(i_def64),            intent(in) :: visit_counts(ntimers, nThreads)
    real(r_def),                 intent(in) :: times(ntimers, nThreads)
    real(r_def), dimension(2,ntimers,nThreads), intent(out) :: max_times, &
                                                               min_times
    real(r_def), dimension(ntimers,nThreads),   intent(out) :: sum_times
    integer(i_def64), dimension(ntimers,nThreads), intent(out) :: sum_counts
    integer, dimension(ntimers), intent(out) :: npes_active
    ! Pretend to use our arguments so that the compiler doesn't complain
    region_names(:,:) = region_names(1,1)
    max_times(1,:,:) = times(:,:)
    min_times(:,:,:) = 0.0
    sum_times(:,:) = 0.0
    sum_counts(:,:) = visit_counts(:,:)
    npes_active(:) = 0
  end subroutine calc_dm_timer_stats

end module dl_timer_parallel
