MODULE model_mod
  use kind_params_mod
  use field_mod
  use shallow_io_mod
  use dl_timer, ONLY: timer_init, timer_report
  implicit none

  integer :: itmax   !< number of timesteps

  real(go_wp) :: dt  !< model timestep (seconds)
  real(go_wp) :: tdt !< 2xdt apart from first step when is just dt

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
    use decomposition_mod, only: decomposition_type
    use subdomain_mod, only: decompose
    use parallel_mod, only: get_rank, set_proc_grid
    use time_smooth_mod, only: time_smooth_init
    use grid_mod, only: grid_type, grid_init
    implicit none
    type(grid_type), intent(inout) :: grid

    !> Grid spacings currently hard-wired, as in original
    !! version of code.
    REAL(KIND=go_wp), PARAMETER :: dxloc=1.0D5, dyloc=1.0D5
    !> Parameter for time smoothing
    REAL(KIND=go_wp), PARAMETER :: alpha_loc = .001d0
    !> Hardwired model time-step (seconds)
    REAL(KIND=go_wp), PARAMETER :: dt_loc = 90.0d0
    !> Problem size, read from namelist
    integer :: m, n
    !> The domain decomposition of the model
    type(decomposition_type) :: decomp

    call timer_init()

    call physical_params_init()

    call read_namelist(m,n,itmax)

    ! Work out the decomposition of the global domain (uses the number
    ! of MPI ranks by default)
    decomp = decompose(m, n)
    call set_proc_grid(decomp%nx, decomp%ny)

    ! Assume the namelist specifies the extent of the
    ! internal domain so allow for boundaries/halos used
    ! to implement periodic boundary conditions.

    ! Set up mesh parameters
    CALL grid_init(grid, decomp, dxloc, dyloc)

    ! Set model time-step
    dt = dt_loc

    ! Initialise time-smoothing module
    CALL time_smooth_init(alpha_loc)

    ! Initialise model IO 'system'
    CALL model_write_init(m,n)

    ! Log model parameters
    CALL print_initial_values(m,n,dxloc,dyloc, dt, alpha_loc)

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

  END SUBROUTINE model_dealloc

END MODULE model_mod
