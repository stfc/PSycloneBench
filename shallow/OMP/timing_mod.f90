MODULE timing_mod

  USE intel_timer_mod
!$ USE omp_lib
  IMPLICIT none

   PRIVATE

   !: double precision (real 8)
   INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

   REAL(wp), PARAMETER :: TOL_ZERO  = 1.0E-10

   REAL(wp), PARAMETER :: REAL_SIZE = 8.0_wp              ! in bytes
   REAL(wp), PARAMETER :: INT_SIZE  = 4.0_wp              ! in bytes
   REAL(wp), PARAMETER :: MB_SIZE   = 1024.0_wp*1024.0_wp ! no. of bytes in 1 MB

   INTEGER, PARAMETER :: numout = 6   ! Unit for stdout

   !-------------------------------------------------------------------
   ! Parameters and types for the timing routines

   INTEGER :: iclk_rate ! Ticks per second of Fortran timer
   INTEGER :: iclk_max  ! Max value that Fortran timer can return

   INTEGER, PARAMETER :: LABEL_LEN  = 128
   INTEGER, PARAMETER :: MAX_TIMERS = 60

   TYPE :: timer_type
      CHARACTER (LABEL_LEN) :: label
      REAL       (KIND=wp)  :: istart
      REAL       (KIND=wp)  :: total
      INTEGER               :: count
   END TYPE timer_type

   type :: timer_stats_type
      CHARACTER (LABEL_LEN) :: label
      real(wp) :: maxTime
      real(wp) :: minTime
      real(wp) :: avgTime
      !> Total time spent in this timer region (sum over all
      !! threads that visit it)
      real(wp) :: totTime
      !> Total no. of trips of all threads to this region
      integer :: tripCount
      !> The no. of OpenMP threads which enter this timer region
      integer :: threadCount
   end type timer_stats_type

   INTEGER, SAVE :: nThreads ! No. of OMP threads being used (1 if no OMP)
                             ! Set in init_time().

   TYPE(timer_type), ALLOCATABLE, SAVE, DIMENSION(:,:) :: timer

   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: itimerCount

   !> Whether to use the Intel-specific rdtsc timer (reads the Time Stamp 
   !! Counter register). If false then the Fortran intrinsic SYSTEM_CLOCK 
   !! is used.
   LOGICAL, PARAMETER :: use_rdtsc_timer = .FALSE.

   !> Whether to output individual timings for each thread
   LOGICAL, PARAMETER :: PER_THREAD_OUTPUT = .FALSE.

   !-------------------------------------------------------------------
   ! Publicly-accessible routines

   PUBLIC timer_init, time_in_s, timer_start, timer_stop, timer_report

 CONTAINS

   !======================================================================

   SUBROUTINE timer_init()
      IMPLICIT none
      ! Set-up timing
      INTEGER :: ji, ith, ierr

! Check that init_time hasn't been called from within an OMP PARALLEL
! region.
!$      IF(omp_get_num_threads() > 1)THEN
!$OMP MASTER
!$         WRITE(*,"('init_time: ERROR: cannot be called from within OpenMP PARALLEL region.')")
!$OMP END MASTER
!$OMP BARRIER
!$         STOP
!$      END IF

      ! Initialise the timer structures

      IF(use_rdtsc_timer)THEN
         iclk_rate =1
         iclk_max = 1
      ELSE
         CALL SYSTEM_CLOCK(COUNT_RATE=iclk_rate, COUNT_MAX=iclk_max)
         WRITE (*,"('System clock, cycles/sec =',I7,', max count = ',I11)") &
                iclk_rate, iclk_max
      END IF

      nThreads = 1
!$    nThreads = omp_get_max_threads()

      WRITE (*,"('Allocating timer structures for ',I3,' threads.')") nThreads

      ALLOCATE(timer(MAX_TIMERS,nThreads), itimerCount(nThreads), &
               Stat=ierr)

      IF(ierr /= 0)THEN
         WRITE (*,*) 'init_time: ERROR: failed to allocate timer structures'
         RETURN
      END IF

!$OMP PARALLEL DO default(none), shared(nThreads,itimerCount,timer), &
!$OMP             private(ith, ji)
      DO ith = 1, nThreads, 1
         DO ji=1,MAX_TIMERS,1
            itimerCount(ith) = 0
            timer(ji,ith)%label  = ""
            timer(ji,ith)%istart = 0_int64
            timer(ji,ith)%total  = 0_wp
            timer(ji,ith)%count  = 0
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE timer_init

!============================================================================

   REAL(wp) FUNCTION time_in_s(clk0,clk1)
     IMPLICIT none
     REAL(wp),    INTENT(in) :: clk0
     REAL(wp), INTENT(inout) :: clk1
     ! This routine only actually returns time in seconds if the
     ! Fortran intrinsic timer (SYSTEM_CLOCK) is being used. Otherwise
     ! iclk_rate has been set to unity and this routine simply returns
     ! the difference between its arguments.

     IF(clk1 < clk0)THEN
        clk1 = clk1 + REAL(iclk_max,wp)
     END IF

     time_in_s =  (clk1 - clk0)/REAL(iclk_rate,wp)

   END FUNCTION time_in_s

!============================================================================

   SUBROUTINE timer_start(label, idx)
     USE intel_timer_mod
     IMPLICIT none
     CHARACTER (*), INTENT(in) :: label
     INTEGER, INTENT(out) :: idx
     INTEGER :: ji, ith, iclk

     IF(LEN_TRIM(label) > LABEL_LEN)THEN
        WRITE(*,"('timer_start: ERROR: length of label >>',(A),'<< exceeds ',I2,' chars')") &
             TRIM(label), LABEL_LEN
        idx = -1
        RETURN
     END IF

     ith = 1
!$   ith = 1 + omp_get_thread_num()

     ! Search for existing timer
     DO ji=1,itimerCount(ith),1
        ! Shorter string is padded with blanks so that lengths match
        IF(timer(ji,ith)%label == label)EXIT 
     END DO

     IF( ji > itimerCount(ith) )THEN
        ! Create a new timer
        itimerCount(ith) = itimerCount(ith) + 1
        IF(itimerCount(ith) > MAX_TIMERS)THEN
           WRITE(*,"('timer_start: ERROR: max. no. of timers exceeded!')")
           WRITE(*,"('timer_start: ERROR: thread = ',I3,'label = ',(A))") ith, label
           idx = -1
           itimerCount(ith) = itimerCount(ith) - 1
           RETURN
        END IF
        timer(itimerCount(ith),ith)%label = TRIM(ADJUSTL(label))
        ji = itimerCount(ith)
     END IF

     ! Increment the count of no. of times we've used this timer
     timer(ji,ith)%count = timer(ji,ith)%count + 1

     ! Return integer tag
     idx = ji

     ! And finally record the current timer value
     IF(use_rdtsc_timer)THEN
        timer(ji,ith)%istart = REAL(getticks(), wp)
     ELSE
        CALL SYSTEM_CLOCK(iclk)
        timer(ji,ith)%istart = REAL(iclk, wp)
     END IF

   END SUBROUTINE timer_start

!============================================================================

   SUBROUTINE timer_stop(itag)
     IMPLICIT none
     INTEGER, INTENT(in) :: itag ! Flag identifying the timer
     ! Stop the specified timer and record the elapsed number of ticks
     ! since it was started.
     INTEGER :: iclk, ith
     INTEGER (kind=int64) :: iclk64

     IF(use_rdtsc_timer)THEN
        iclk64 = getticks()
     ELSE
        CALL SYSTEM_CLOCK(iclk)
        iclk64 = INT(iclk, int64)
     END IF
     
     IF(itag < 1)RETURN

     ith = 1
!$   ith = 1 + omp_get_thread_num()

     IF( use_rdtsc_timer )THEN
        timer(itag,ith)%total = timer(itag,ith)%total + &
                          (REAL(iclk64,wp) - timer(itag,ith)%istart)

     ELSE
        IF( iclk < timer(itag,ith)%istart )THEN
           iclk64 = iclk64 + INT(iclk_max,int64)
        END IF

        timer(itag,ith)%total = timer(itag,ith)%total + &
                          (REAL(iclk64,wp) - timer(itag,ith)%istart)
     END IF

   END SUBROUTINE timer_stop

!============================================================================

   SUBROUTINE timer_report()
     implicit none
     integer :: ji, jt, ierr, jtimer
     real    :: wtime, wtimePerCall
     !> Array for tot compute time of each thread
     real, allocatable, dimension(:) :: totThrdTime
     real(wp) :: inv_thr_cnt
     !> Min and Max total compute time (summed over all timed
     !! regions) of any of the threads
     real(wp) :: minTime, maxTime
     integer  :: num_distinct_timers
     logical  :: print_summary, found
     type(timer_stats_type), allocatable, dimension(:) :: timer_stats

     allocate( timer_stats(MAX_TIMERS), &
               totThrdTime(nThreads),  &
               Stat=ierr)

     IF(ierr /= 0)THEN
        WRITE(*,"('timer_report: ERROR: allocate failed!')")
        RETURN
     END IF

     WRITE(*,"(/'====================== Timing report ==============================')")
     IF(use_rdtsc_timer)THEN
        WRITE(*," (' Timed using Intel Time Stamp Counter. Units are millions of cycles.')")
     ELSE
        WRITE(*," (' Timed using Fortran SYSTEM_CLOCK intrinsic. Units are seconds.')")
     END IF

      IF(PER_THREAD_OUTPUT)THEN
         WRITE(*," (67('-'))")
         WRITE(*," ('Region',26x,'Counts      Total         Average')")
         WRITE(*," (67('-'))")
      END IF

      ! Zero array used to cacluate tot times
      totThrdTime(:) = 0.0
      print_summary = .TRUE.
      num_distinct_timers = 0

      DO jt = 1, nThreads, 1

         IF(PER_THREAD_OUTPUT)THEN
            IF(jt > 1)THEN
               WRITE(*," (34('- '))")
            END IF

            IF(nThreads > 1)WRITE(*," ('Thread ',I3)") jt
         END IF

         ! Loop over all timed sections for this thread
         DO ji=1,itimerCount(jt),1

            ! Check whether we've seen a timer with this label before
            found = .FALSE.
            do jtimer = 1, num_distinct_timers
               if(timer_stats(jtimer)%label(1:32) == timer(ji,jt)%label(1:32))then
                  found = .TRUE.
                  exit
               end if
            end do

            if(.not. found)then
               ! This is a timer we've not seen before so set up a
               ! struct to store its stats
               num_distinct_timers = num_distinct_timers + 1
               jtimer = num_distinct_timers
               timer_stats(jtimer)%label = timer(ji,jt)%label
               timer_stats(jtimer)%maxTime = -1.0
               timer_stats(jtimer)%minTime = HUGE(1.0)
               timer_stats(jtimer)%avgTime = 0.0
               timer_stats(jtimer)%totTime = 0.0
               timer_stats(jtimer)%tripCount = 0
               timer_stats(jtimer)%threadCount = 1
            else
               ! Increment the count of threads that have visited this timed
               ! region
               timer_stats(jtimer)%threadCount = &
                                timer_stats(jtimer)%threadCount + 1
            end if

            IF(use_rdtsc_timer)THEN
               ! Units are MILLIONS of cycles (counts)
               wtime = 1.0e-6_wp * timer(ji,jt)%total
            ELSE
               wtime = time_in_s(0._wp,timer(ji,jt)%total)
            END IF

            ! Sum the trip-counts of all threads for this region
            timer_stats(jtimer)%tripCount = timer_stats(jtimer)%tripCount &
                                   + timer(ji,jt)%count

            ! Calculate the time per call in this region
            wtimePerCall = wtime/REAL(timer(ji,jt)%count)

            ! Total time spent in this section by all threads that
            ! visit it
            timer_stats(jtimer)%totTime = timer_stats(jtimer)%totTime + wtime

            ! Sum the timed sections for this thread
            totThrdTime(jt) = totThrdTime(jt) + wtime

            ! Add the time for this thread to the sum for the average
            timer_stats(jtimer)%avgTime = timer_stats(jtimer)%avgTime &
                                          + wtimePerCall

            ! Update the min and max values
            if(timer_stats(jtimer)%minTime > wtimePerCall)then
               timer_stats(jtimer)%minTime = wtimePerCall
            end if
            if(timer_stats(jtimer)%maxTime < wtimePerCall)then
               timer_stats(jtimer)%maxTime = wtimePerCall
            end if

            IF(PER_THREAD_OUTPUT)THEN
               ! Truncate the label to 32 chars for table-formatting purposes
               WRITE(*,"((A),1x,I4,2x,E13.6,2x,E13.6)") &
                            timer(ji,jt)%label(1:32), timer(ji,jt)%count, &
                            wtime, wtimePerCall
            END IF
         END DO
      END DO

      IF(print_summary)THEN
         WRITE(*," (91('='))")
         WRITE(*," ('                     Times averaged over ',I3,' threads')") nThreads
         WRITE(*," ('Region',29x,'Counts',7x,'Min',8x,'Average',7x,'Max',7x,'Total')")
         WRITE(*," (91('-'))")
         jt = 1
         DO ji=1,num_distinct_timers,1

            inv_thr_cnt = 1.0d0 / REAL(timer_stats(ji)%threadCount)

            WRITE(*,"((A),1x,E10.4,1x,E11.5,1x,E11.5,1x,E11.5,1x,E11.5)") &
                  timer_stats(ji)%label(1:32), &
                  timer_stats(ji)%tripCount*inv_thr_cnt, &
                  timer_stats(ji)%minTime, &
                  timer_stats(ji)%avgTime*inv_thr_cnt, &
                  timer_stats(ji)%maxTime,&
                  timer_stats(ji)%totTime*inv_thr_cnt
         END DO
         WRITE(*," (91('='))")
      END IF

      ! Calc. the min and max compute time over the set of threads
      minTime = HUGE(1.0)
      maxTime = -1.0

      IF(PER_THREAD_OUTPUT) WRITE(*," ('Thread  Tot. compute time')")
      DO jt = 1, nThreads, 1
         IF(PER_THREAD_OUTPUT)THEN
            WRITE(*,"(2x,I3,2x,F13.4)") jt, totThrdTime(jt)
         END IF
         IF(minTime > totThrdTime(jt))minTime = totThrdTime(jt)
         IF(maxTime < totThrdTime(jt))maxTime = totThrdTime(jt)
      END DO

      IF(PER_THREAD_OUTPUT) WRITE(*," (82('='))")

      WRITE(*,"('Tot /thread compute t, min=',E12.5,' max=',E12.5,', % diff=',F6.2)") &
            minTime, maxTime,   &
            100.0*(maxTime-minTime)/minTime
      WRITE(*," (82('='))")

   END SUBROUTINE timer_report

!============================================================================

END MODULE timing_mod
