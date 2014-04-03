MODULE model
  USE shallow_IO
  USE timing, ONLY: timer_init, timer_report
  IMPLICIT none

  INTEGER :: m, n      !< global domain size
  INTEGER :: mp1, np1  !< m+1 and n+1 == array extents

  INTEGER :: itmax   !< number of timesteps

  ! solution arrays
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) ::                & 
                             u, v, p, unew, vnew, pnew,       & 
                             uold, vold, pold, cu, cv, z, h, psi  

CONTAINS

  !================================================

  SUBROUTINE model_init()
    IMPLICIT none

    CALL timer_init()

    CALL read_namelist(m,n,itmax)

    !     Set up arrays
    MP1 = M+1
    NP1 = N+1

    CALL model_alloc(mp1, np1)

    CALL model_write_init(m,n)

  END SUBROUTINE model_init

  !================================================

  SUBROUTINE model_finalise()
    IMPLICIT none

    CALL model_write_finalise()

    CALL timer_report()

    CALL model_dealloc()
  
  END SUBROUTINE model_finalise

  !================================================

  SUBROUTINE model_alloc(idimx, idimy)
    IMPLICIT none
    INTEGER, INTENT(in) :: idimx, idimy

    ALLOCATE( u(idimx,idimy),    v(idimx,idimy),    p(idimx,idimy) ) 
    ALLOCATE( unew(idimx,idimy), vnew(idimx,idimy), pnew(idimx,idimy) ) 
    ALLOCATE( uold(idimx,idimy), vold(idimx,idimy), pold(idimx,idimy) )
    ALLOCATE( cu(idimx,idimy),   cv(idimx,idimy) ) 
    ALLOCATE( z(idimx,idimy),    h(idimx,idimy),    psi(idimx,idimy) ) 

  END SUBROUTINE model_alloc

  !================================================

  SUBROUTINE model_dealloc()
    IMPLICIT none

    !> Free memory \todo Move to model_finalise()
    DEALLOCATE( u, v, p, unew, vnew, pnew, uold, vold, pold )
    DEALLOCATE( cu, cv, z, h, psi ) 

  END SUBROUTINE model_dealloc

END MODULE model
