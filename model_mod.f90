MODULE model_mod
  use kind_params_mod
  use grid_mod
  use field_mod
  use timing_mod, ONLY: timer_init, timer_report
  implicit none

  public

  !> Problem size, read from namelist
  !! \todo Will live with the grid object
  integer, save  :: jpi, jpj
  !> start-end and record time steps
  integer, save  :: nit000, nitend, irecord 
  !> type (source?) of grid
  integer, save  :: jphgr_msh

  REAL(wp), save :: rdt                                   !< time step
  REAL(wp), save :: cbfr                                  !< bottom friction coefficient
  REAL(wp), save :: visc                                  !< backgroud/constant viscosity 

  REAL(wp), save :: dep_const                             !< constant depth

  INTEGER, ALLOCATABLE, save :: pt(:,:)   ! properties of t-cells 
                                    ! 1: water cell within computational domain
                                    ! 0: land cell
                                    !-1: water cell outside computational domain
  
  REAL(wp), ALLOCATABLE, save :: ht(:,:), hu(:,:), hv(:,:), hf(:,:) 
  
  REAL(wp), ALLOCATABLE, save :: sshb(:,:), sshb_u(:,:), sshb_v(:,:)
  REAL(wp), ALLOCATABLE, save :: sshn(:,:), sshn_u(:,:), sshn_v(:,:)
  REAL(wp), ALLOCATABLE, save :: ssha(:,:), ssha_u(:,:), ssha_v(:,:)
  
  REAL(wp), ALLOCATABLE, save :: un(:,:),  vn(:,:)
  REAL(wp), ALLOCATABLE, save :: ua(:,:),  va(:,:)
         

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

CONTAINS

  !================================================

  subroutine model_init(grid)
    use gocean2d_io_mod
    implicit none
    type(grid_type), intent(inout) :: grid

    !> Problem size, read from namelist
    integer :: jpiglo, jpjglo
    real(wp) :: dx, dy

    integer :: ierr(6)
    integer :: ios

    !! Read in model setup parameters 
    NAMELIST/namctl/ jpiglo, jpjglo, jphgr_msh, &
                     dx    , dy    , dep_const, &
                     nit000, nitend, irecord  , &
                     rdt   , cbfr  , visc

    CALL timer_init()

    !! Default values

    jpiglo      =      50               !  number of columns of model grid
    jpjglo      =     100               !  number of rows of model grid
    jphgr_msh   =       1               !  type of grid (0: read in a data file; 1: setup with following parameters)
    dx          =   1000._wp            !  grid size in x direction (m)
    dy          =   1000._wp            !  grid size in y direction (m)
    dep_const   =    100._wp            !  constant depth (m)
    nit000      =       1               !  first time step
    nitend      =    1000               !  end time step
    irecord     =       1               !  intervals to save results
    rdt         =     10._wp            !  size of time step (second) 
    cbfr        =   0.001_wp            !  bottom friction coefficeint
    visc        =     100._wp           !  horiz. kinematic viscosity coeff. 
 
    OPEN(111, file='namelist', STATUS='OLD')
    REWIND(111)
    READ(111, NML=namctl, IOSTAT = ios, ERR = 901)
901 IF(ios /= 0) STOP "err found in reading namelist file"
    WRITE(*,NML=namctl)
    
    CLOSE(111)

    ! Set up mesh parameters
    CALL grid_init(grid, jpiglo, jpjglo, dx, dy)

    ! Store the grid dimensions in module variables - this is a temporary
    ! fix prior to carrying everything around in the grid object.
    jpi = jpiglo
    jpj = jpjglo

    ALLOCATE(ht(jpi,jpj), hu(0:jpi,jpj), hv(jpi,0:jpj), hf(0:jpi,0:jpj), STAT=ierr(1))

    ALLOCATE(sshb(jpi,jpj), sshb_u(0:jpi,jpj), sshb_v(jpi,0:jpj), STAT=ierr(2))
    ALLOCATE(sshn(jpi,jpj), sshn_u(0:jpi,jpj), sshn_v(jpi,0:jpj), STAT=ierr(3))
    ALLOCATE(ssha(jpi,jpj), ssha_u(0:jpi,jpj), ssha_v(jpi,0:jpj), STAT=ierr(4))

    ALLOCATE(un(0:jpi,jpj), vn(jpi,0:jpj), ua(0:jpi,jpj), va(jpi,0:jpj), STAT=ierr(5))

    ALLOCATE(pt(0:jpi+1,0:jpj+1), STAT=ierr(6))

    IF(ANY(ierr /= 0, 1)) STOP "Failed to allocate solution arrays"

    ! Initialise model IO 'system'
    call model_write_init(jpiglo, jpjglo, irecord)

  end subroutine model_init

  !================================================

  SUBROUTINE model_finalise()
    use gocean2d_io_mod, only: model_write_finalise
    IMPLICIT none

    CALL model_write_finalise()

    CALL timer_report()

    CALL model_dealloc()
  
  end subroutine model_finalise

  !================================================

  subroutine model_dealloc()
    implicit none

    DEALLOCATE(ht, hu, hv, hf)
    DEALLOCATE(sshb, sshb_u, sshb_v)
    DEALLOCATE(sshn, sshn_u, sshn_v)
    DEALLOCATE(ssha, ssha_u, ssha_v)
    DEALLOCATE(un, vn, ua, va)
    DEALLOCATE(pt)

  END SUBROUTINE model_dealloc

END MODULE model_mod
