!> Provides an interface to the Intel-specifc getticks routine which
!! queries the Time-Stamp Counter (rdtsc) register.
MODULE intel_timer_mod
   USE iso_c_binding
   IMPLICIT none

   PUBLIC

   ! Need 64-bit integers when using the Intel counter for timing
   INTEGER, PARAMETER :: int64 = SELECTED_INT_KIND(14)

   INTERFACE

      FUNCTION rdtsc_available()  bind(C, name="rdtsc_available")
        USE iso_c_binding
        IMPLICIT none
        INTEGER (C_INT) :: rdtsc_available
      END FUNCTION rdtsc_available

      FUNCTION getticks()  bind(C, name="getticks")
        USE iso_c_binding
        IMPLICIT none
        INTEGER (C_INT64_T) :: getticks
      END FUNCTION getticks

   END INTERFACE

END MODULE intel_timer_mod
