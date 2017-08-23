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
      !> The name of this timed region
      CHARACTER (LABEL_LEN) :: label
      !> Time at which region was most recently entered
      REAL       (KIND=wp)  :: istart
      !> Total time spent in this timed region (accumulated over
      !! all visits).
      REAL       (KIND=wp)  :: total
      !> The no. of times this timed region has been executed.
      INTEGER               :: count
      !> The no. of repeated intervals within this timed region.
      !! Used in timer_report() to produce a mean time per repeat.
      !! Default value is 1. User can specify value in call to
      !! timer_start().
      INTEGER               :: nrepeat
   END TYPE timer_type

   INTEGER, SAVE :: nThreads ! No. of OMP threads being used (1 if no OMP)
                             ! Set in init_time().

   TYPE(timer_type), ALLOCATABLE, SAVE, DIMENSION(:,:) :: timer

   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) :: itimerCount

   ! Whether to use the Intel-specific rdtsc timer (reads the Time Stamp 
   ! Counter register). If false then the Fortran intrinsic SYSTEM_CLOCK 
   ! is used.
   LOGICAL, PARAMETER :: use_rdtsc_timer = .FALSE.

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

   SUBROUTINE timer_start(label, idx, nrepeat)
      USE intel_timer_mod
      IMPLICIT none
      CHARACTER (*), INTENT(in) :: label
      INTEGER, INTENT(out) :: idx
      !> The number of repeated intervals inside this timed region.
      !! Used to report a time per interval in the report generated
      !! by timer_report().
      INTEGER, INTENT(in), OPTIONAL :: nrepeat
      INTEGER :: ji, ith, iclk

      IF(LEN_TRIM(label) > LABEL_LEN)THEN
         WRITE(*,"('timer_start: ERROR: length of label >>',(A),'<< exceeds ',I2,' chars')") &
              TRIM(label), LABEL_LEN
         idx = -1
         RETURN
      END IF

      ith = 1
!$    ith = 1 + omp_get_thread_num()

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
            WRITE(*,"('timer_start: ERROR: thread = ',I3,'label = ',(A))") &
                  ith, label
            idx = -1
            itimerCount(ith) = itimerCount(ith) - 1
            RETURN
         END IF

         ji = itimerCount(ith)

         ! Initialise this new timer structure
         timer(ji, ith)%label = TRIM(ADJUSTL(label))
         if(present(nrepeat))then
            timer(ji, ith)%nrepeat = nrepeat
         else
            ! No repeat specified so default to a value of unity.
            timer(ji, ith)%nrepeat = 1
         end if

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
!$    ith = 1 + omp_get_thread_num()

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

   !==========================================================================

   SUBROUTINE timer_report()
      IMPLICIT none
      INTEGER       :: jt
      LOGICAL       :: have_repeats

      ! Check whether any of our timed regions have a non-unity
      ! no. of repeats
      have_repeats = .false.
      do jt = 1, nThreads, 1
         if( ANY( timer(1:itimerCount(jt),jt)%nrepeat > 1) )then
            have_repeats = .TRUE.
            exit
         end if
      end do

      if(have_repeats)then
         call timer_report_with_repeats()
      else
         call timer_report_no_repeats()
      end if

    end SUBROUTINE timer_report

    !==========================================================================

    subroutine timer_report_no_repeats()
      implicit none
      INTEGER       :: ji, jt
      REAL(KIND=wp) :: wtime

      WRITE(*,"(/'====================== Timing report ==============================')")
      IF(use_rdtsc_timer)THEN
         WRITE(*," ('    Timed using Intel Time Stamp Counter. Units are counts.')")
      ELSE
         WRITE(*," (' Timed using Fortran SYSTEM_CLOCK intrinsic. Units are seconds.')")
      END IF
      WRITE(*," ('-------------------------------------------------------------------')")
      WRITE(*," ('Region',26x,'Counts      Total         Average')")
      WRITE(*," ('-------------------------------------------------------------------')")
      DO jt = 1, nThreads, 1
         IF(jt > 1)THEN
            WRITE(*," ('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')")
         END IF

         IF(nThreads > 1)WRITE(*," ('Thread ',I3)") jt

         DO ji=1,itimerCount(jt),1

            IF(use_rdtsc_timer)THEN
               wtime = timer(ji,jt)%total
            ELSE
               wtime = time_in_s(0._wp,timer(ji,jt)%total)
            END IF

            ! Truncate the label to 32 chars for table-formatting purposes
            WRITE(*,"((A),1x,I4,2x,E13.6,2x,E13.6)") &
                            timer(ji,jt)%label(1:32), timer(ji,jt)%count, &
                            wtime, wtime/REAL(timer(ji,jt)%count)
         END DO
      END DO
      WRITE(*," ('===================================================================')")

   END SUBROUTINE timer_report_no_repeats

   !==========================================================================

   SUBROUTINE timer_report_with_repeats()
      IMPLICIT none
      INTEGER       :: ji, jt
      REAL(KIND=wp) :: wtime

      WRITE(*,"(/34('='),' Timing report ',34('='))")
      IF(use_rdtsc_timer)THEN
         WRITE(*," ('    Timed using Intel Time Stamp Counter. Units are counts.')")
      ELSE
         WRITE(*," (' Timed using Fortran SYSTEM_CLOCK intrinsic. Units are seconds.')")
      END IF
      WRITE(*,"(83('-'))")
      WRITE(*,"('Region',26x,'Counts      Total         Average    Average/repeat')")
      WRITE(*,"(83('-'))")
      DO jt = 1, nThreads, 1
         IF(jt > 1)THEN
            WRITE(*,"(41('- '))")
         END IF

         IF(nThreads > 1)WRITE(*," ('Thread ',I3)") jt

         DO ji=1,itimerCount(jt),1

            IF(use_rdtsc_timer)THEN
               wtime = timer(ji,jt)%total
            ELSE
               wtime = time_in_s(0._wp,timer(ji,jt)%total)
            END IF

            ! Truncate the label to 32 chars for table-formatting purposes
            WRITE(*,"((A),1x,I4,2x,E13.6,2x,E13.6,2x,E13.6)")          &
                            timer(ji,jt)%label(1:32), timer(ji,jt)%count, &
                            wtime, wtime/REAL(timer(ji,jt)%count),        &
                            wtime/REAL(timer(ji,jt)%count * timer(ji,jt)%nrepeat)
         END DO
      END DO
      WRITE(*,"(83('='))")

   END SUBROUTINE timer_report_with_repeats

!============================================================================

END MODULE timing_mod
