MODULE intel_timer_mod
   USE iso_c_binding
   IMPLICIT none

   PUBLIC

   ! Need 64-bit integers when using the Intel counter for timing
   INTEGER, PARAMETER :: int64 = SELECTED_INT_KIND(14)

   INTERFACE

      FUNCTION getticks()  bind(C,name="getticks")
        USE iso_c_binding
        IMPLICIT none
        INTEGER (C_INT64_T) :: getticks
      END FUNCTION getticks

   END INTERFACE

END MODULE intel_timer_mod
