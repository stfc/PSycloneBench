MODULE model_mod
  use kind_params_mod
  use grid_mod
  use field_mod
  use timing_mod, ONLY: timer_init, timer_report
  implicit none

  !> start-end and record time steps
  integer  :: nit000, nitend, irecord 
  !> type (source?) of grid
  INTEGER  :: jphgr_msh

  REAL(wp) :: rdt                                   !< time step
  REAL(wp) :: cbfr                                  !< bottom friction coefficient
  REAL(wp) :: visc                                  !< backgroud/constant viscosity 

  REAL(wp) :: dep_const                             !< constant depth

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
    implicit none
    type(grid_type), intent(inout) :: grid

    !> Problem size, read from namelist
    integer :: jpiglo, jpjglo
    real(wp) :: dx, dy

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
    visc        =     100._wp            !  horizontal kinematic viscosity coefficient 
 
 
    OPEN(111, file='namelist', STATUS='OLD')
    REWIND(111)
    READ(111, NML=namctl, IOSTAT = ios, ERR = 901)
901 IF(ios /= 0) STOP "err found in reading namelist file"
    WRITE(*,NML=namctl)
    
    CLOSE(111)

    ! Set up mesh parameters
    CALL grid_init(grid, jpiglo, jpjglo, dx, dy)

    ! Initialise model IO 'system'
    !CALL model_write_init(m,n)

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
