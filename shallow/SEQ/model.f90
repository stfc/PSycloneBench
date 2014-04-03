MODULE model
  USE shallow_IO
  USE timing, ONLY: timer_init, timer_report
  IMPLICIT none

  INTEGER :: m, n      !< global domain size
  INTEGER :: mp1, np1  !< m+1 and n+1 == array extents

  INTEGER :: itmax   !< number of timesteps

CONTAINS

  !================================================

  SUBROUTINE model_init()
    IMPLICIT none

    CALL timer_init()

    CALL read_namelist(m,n,itmax)

    CALL model_write_init(m,n)

  END SUBROUTINE model_init

  !================================================

  SUBROUTINE model_finalise()
    IMPLICIT none

    CALL model_write_finalise()

    CALL timer_report()
  
  END SUBROUTINE model_finalise

END MODULE model
