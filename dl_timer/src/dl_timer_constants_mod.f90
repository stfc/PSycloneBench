!> Module containing constant definitions
module dl_timer_constants_mod
   implicit none

   !------------------------------------------------------------------
   ! Type definitions
   !> double precision (real 8)
   INTEGER, PARAMETER :: r_def = SELECTED_REAL_KIND(12,307)
   !> 32-bit integer
   INTEGER, PARAMETER :: i_def32 = selected_int_kind(9)
   !> 64-bit integer
   INTEGER, PARAMETER :: i_def64 = selected_int_kind(12)
   !> Single precision
   INTEGER, PARAMETER :: sp = KIND(1.0)

   !> Tolerance below which we consider a number to be zero
   REAL(r_def), PARAMETER :: TOL_ZERO  = 1.0E-10

   !> Unit for stdout
   INTEGER, PARAMETER :: OUT_UNIT = 6
   INTEGER, PARAMETER :: ERR_UNIT = 6

   !> Maximum length of the label for a timed region
   INTEGER, PARAMETER :: LABEL_LEN  = 50
   !> Maximum number of distinct timed regions that an application
   !! may have 
   INTEGER, PARAMETER :: MAX_TIMERS = 30
   !> How many samples to keep when recording a time-line
   INTEGER, PARAMETER :: TIME_SERIES_LEN = 10000

 end module dl_timer_constants_mod
