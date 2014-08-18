MODULE model_mod
  use kind_params_mod
  use grid_mod
  use field_mod
  use shallow_io_mod
  use timing_mod, ONLY: timer_init, timer_report
  implicit none

  INTEGER :: itmax   !< number of timesteps

  TYPE(scalar_field_type) :: dt  !< model timestep (seconds)
  TYPE(scalar_field_type) :: tdt !< 2xdt apart from first step when is just dt

  ! solution arrays
  ! Fields are allocated with extents (M+1,N+1).
  ! Presumably the extra row and column are needed for periodic BCs.

  ! The finite-difference variables are defined on the mesh as follows:
  !            u              u          
  !    -------->-----P,H------>-----P,H
  !    |              |              |   
  !    |        rot   |        rot   |   
  !v  /|\      o     /|\ v    o     /|\ v
  !    |              |              |   
  !    |              |              |   
  !    -------->-----P,H------>-----P,H
  !    |              |              |   
  !    |        rot   |        rot   |   
  !v  /|\      o     /|\ v    o     /|\ v
  !    |              |              |   
  !    |              |              |   
  !    |------->------|------->------|   
  !            u              u          

  ! Calling the point at which P,H etc. are defined 'T' and the
  ! point at which rotation/vorticity defined 'f' then:
  !
  !          f(i,j+1) v(i,j+1)  f(i+1,j+1)  v(i+1,j+1)  f(i+2,j+1)
  !
  !          u(i,j)   T(i,j)    u(i+1,j)    T(i+1,j)    u(i+2,j)
  !
  !          f(i,j)   v(i,j)    f(i+1,j)    v(i+1,j)    f...
  !
  !          u(i,j-1) T(i,j-1)  u(i+1,j-1)  T...        u...
  !
  !  So, T is dual of f.
  !
  ! This is consistent with the form of the original code:
  !     DO J=1,N
  !        DO I=1,M
  !           CU(I+1,J) = .5*(P(I+1,J)+P(I,J))*U(I+1,J)
  !           CV(I,J+1) = .5*(P(I,J+1)+P(I,J))*V(I,J+1)
  !           Z(I+1,J+1) = &
  !              (FSDX*(V(I+1,J+1)-V(I,J+1))-FSDY*(U(I+1,J+1)-U(I+1,J))) &
  !                          /(P(I,J)+P(I+1,J)+P(I+1,J+1)+P(I,J+1))
  !           H(I,J) = P(I,J)+.25*(U(I+1,J)*U(I+1,J)+U(I,J)*U(I,J)     & 
  !                +V(I,J+1)*V(I,J+1)+V(I,J)*V(I,J))
  !        END DO
  !     END DO


CONTAINS

  !================================================

  subroutine model_init(grid)
    use physical_params_mod, only: physical_params_init
    use time_smooth_mod, only: time_smooth_init
    IMPLICIT none
    type(grid_type), intent(inout) :: grid

    !> Grid spacings currently hard-wired, as in original
    !! version of code.
    REAL(KIND=wp), PARAMETER :: dxloc=1.0E5, dyloc=1.0E5
    !> Parameter for time smoothing
    REAL(KIND=wp), PARAMETER :: alpha_loc = .001
    !> Hardwired model time-step (seconds)
    REAL(KIND=wp), PARAMETER :: dt_loc = 90.
    !> Problem size, read from namelist
    integer :: m, n

    call timer_init()

    call physical_params_init()

    CALL read_namelist(m,n,itmax)

    ! Assume the namelist specifies the extent of the
    ! internal domain so allow for boundaries/halos used
    ! to implement periodic boundary conditions.

    ! Set up mesh parameters
    CALL grid_init(grid, m, n, dxloc, dyloc)

    ! Set model time-step
    CALL set(dt, dt_loc)

    ! Initialise time-smoothing module
    CALL time_smooth_init(alpha_loc)

    ! Initialise model IO 'system'
    CALL model_write_init(m,n)

    ! Log model parameters
    CALL print_initial_values(m,n,dxloc,dyloc, dt%data, alpha_loc)

  END SUBROUTINE model_init

  !================================================

  SUBROUTINE model_finalise()
    IMPLICIT none

    CALL model_write_finalise()

    CALL timer_report()

    CALL model_dealloc()
  
  END SUBROUTINE model_finalise

  !================================================

  SUBROUTINE model_dealloc()
    IMPLICIT none

    !> Free memory \todo Move to model_finalise()
    !DEALLOCATE( u, v, p, unew, vnew, pnew, uold, vold, pold )
    !DEALLOCATE( cu, cv, z, h, psi ) 

  END SUBROUTINE model_dealloc

END MODULE model_mod
