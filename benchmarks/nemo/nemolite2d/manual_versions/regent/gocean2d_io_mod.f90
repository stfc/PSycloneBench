!This is equivalent to the ../../common/gocean2d_io_mod.f90 read_namelist
!routine, however all other functions from that module are stripped away
!to avoid dependencies on other modules.
module gocean2d_io_mod
  use kind_params_mod
  !  use field_mod
  implicit none

  private

  logical, save :: l_out   !< Whether or not to produce output  
  integer, save :: mprint  !< frequency of output    
  !> Extents of arrays to write \todo Carry with field object.
  integer, save :: jpi, jpj 

  public read_namelist
!  public model_write_init, model_write, model_write_finalise

contains

  !===================================================

  !> Reads the namelist file for user-specified control 
  !! parameters
  SUBROUTINE read_namelist(jpiglo, jpjglo, dx, dy,   &
                           nit000, nitend, irecord,  &
                           jphgr_msh, dep_const, rdt, cbfr, visc) bind (c, name='read_namelist')
    use, intrinsic :: iso_c_binding
    implicit none
    !> Extent of the mask that describes the area that
    !! contains the simulation domain
    integer,  intent(out) :: jpiglo, jpjglo
    !> Grid spacing
    real(go_wp), intent(out) :: dx, dy
    !> Start time-step, stop time-step and output interval
    integer, intent(out) :: nit000, nitend, irecord
    !> Whether grid is to be read from file or not
    integer, intent(out) :: jphgr_msh
    !> (Constant) depth of water in domain
    real(go_wp), intent(out) :: dep_const
    !> Model time-step (seconds?)
    real(go_wp), intent(out) :: rdt
    !> Coefficient of bottom friction
    real(go_wp), intent(out) :: cbfr
    !> Background/constant viscosity
    real(go_wp), intent(out) :: visc
    ! namelist input 
    character (len=8) :: nml_name = "namelist" 
    integer :: input_unit = 99
    integer :: ios

    !! Read in model setup parameters 
    NAMELIST/namctl/ jpiglo, jpjglo, jphgr_msh, &
                     dx    , dy    , dep_const, &
                     nit000, nitend, irecord  , &
                     rdt   , cbfr  , visc

    !! Default values

    jpiglo      =      50               !  number of columns of model grid
    jpjglo      =     100               !  number of rows of model grid
    jphgr_msh   =       1               !  type of grid (0: read in a data file; 1: setup with following parameters)
    dx          =   1000._go_wp            !  grid size in x direction (m)
    dy          =   1000._go_wp            !  grid size in y direction (m)
    dep_const   =    100._go_wp            !  constant depth (m)
    nit000      =       1               !  first time step
    nitend      =    1000               !  end time step
    irecord     =       1               !  intervals to save results
    rdt         =     10._go_wp            !  size of time step (second) 
    cbfr        =   0.001_go_wp            !  bottom friction coefficeint
    visc        =     100._go_wp           !  horiz. kinematic viscosity coeff. 
 
    OPEN(input_unit, file=nml_name, STATUS='OLD')
    REWIND(input_unit)
    READ(input_unit, NML=namctl, IOSTAT = ios)
    IF(ios /= 0) STOP "err found in reading namelist file"
    WRITE(*,NML=namctl)
    
    CLOSE(unit=input_unit)

  END SUBROUTINE read_namelist

END MODULE gocean2d_io_mod
