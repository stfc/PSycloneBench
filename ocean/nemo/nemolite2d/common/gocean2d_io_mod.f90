module gocean2d_io_mod
  use kind_params_mod
  use field_mod
  implicit none

  private

  logical, save :: l_out   !< Whether or not to produce output  
  integer, save :: mprint  !< frequency of output    
  !> Extents of arrays to write \todo Carry with field object.
  integer, save :: jpi, jpj 

  public read_namelist
  public model_write_init, model_write, model_write_finalise

contains

  !===================================================

  !> Reads the namelist file for user-specified control 
  !! parameters
  SUBROUTINE read_namelist(jpiglo, jpjglo, dx, dy,   &
                           nit000, nitend, irecord,  &
                           jphgr_msh, dep_const, rdt, cbfr, visc)
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

  !===================================================
  
  SUBROUTINE model_write_init(m, n, interval)
    IMPLICIT none
    INTEGER, INTENT(in) :: m, n, interval

    jpi = m
    jpj = n
    mprint = interval
    l_out = .true.

  END SUBROUTINE model_write_init
 
  !===================================================

  !> Write data for the current time step
  subroutine model_write(grid, istp, ht, sshn, un, vn)
    use kind_params_mod
    use grid_mod
    use parallel_mod, only: get_rank
    implicit none
    type(grid_type), intent(in) :: grid
    integer, intent(in) :: istp
    type(r2d_field), intent(in) :: ht, sshn, un, vn
    ! Locals
    integer :: ji, jj
    real(go_wp) :: rtmp1, rtmp2
    character(len=11) :: fname
    
    if( l_out .and. (mod(istp, mprint) .eq. 0) ) then

       ! output model results
       write(fname, '(I5.5,"_",I5.5)') istp, get_rank()
       open(21, file='go2d_'//fname//'.dat', STATUS='UNKNOWN', &
            action='write')

       rewind(21)

       ! Loop over 'internal' T points
       DO jj = sshn%internal%ystart, sshn%internal%ystop, 1
          DO ji = sshn%internal%xstart, sshn%internal%xstop, 1

             rtmp1 = 0.5_go_wp * (un%data(ji-1,jj) + un%data(ji,jj))
             rtmp2 = 0.5_go_wp * (vn%data(ji,jj-1) + vn%data(ji,jj))

             ! write "x-coord, y-coord, depth, ssh, u-velocity,
             ! v-velocity" to ASCII files

              write(21,'(6e16.7)') grid%xt(ji,jj), grid%yt(ji,jj), &
                                   ht%data(ji,jj), sshn%data(ji,jj), &
                                   rtmp1, rtmp2 
          END DO
          WRITE(21,*)
       END DO
          
       close(21)

       write(fname, '("sshn_",I5.5)') get_rank()
       open(21, file=trim(fname)//'.dat', STATUS='UNKNOWN', action='write')

       ! ARPDBG remove this debug-output code
       DO jj = 1, sshn%whole%ny, 1
          DO ji = 1, sshn%whole%nx, 1
             !write(21,*) grid%xt(ji,jj), grid%yt(ji,jj), sshn%data(ji,jj)
             write(21,*) ji, jj, sshn%data(ji,jj)
          end do
          write(21,*)
       end do

       close(21)
       
    end if ! we're doing output

  end subroutine model_write

  !===================================================

  SUBROUTINE model_write_finalise()
    IMPLICIT none
    
    ! Nothing to do here as we open and close ASCII file for
    ! each write.

  END SUBROUTINE model_write_finalise

  !===================================================

  ! Check error code
!!$  subroutine check(status, text)
!!$    implicit none
!!$      
!!$    integer, intent(in) :: status
!!$    character (len=*)   :: text
!!$    
!!$    if (status /= 0) then
!!$       write(6,*) "error ", status
!!$       write(6,*) text
!!$       stop 2
!!$    endif
!!$
!!$  end subroutine check

  !===================================================

END MODULE gocean2d_io_mod
