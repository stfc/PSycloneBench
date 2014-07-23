module gocean2d_io_mod
  use kind_params_mod
  use field_mod
  implicit none

  private

  logical, save :: l_out   !< Whether or not to produce output  
  integer, save :: mprint  !< frequency of output    
  !> Extents of arrays to write \todo Carry with field object.
  integer, save :: jpi, jpj 

  !> Interface to logging routines
  interface model_write_log
     module procedure write_log_ir, write_log_r
  end interface

  public read_namelist
  public model_write_init, model_write, model_write_finalise
  public model_write_log

contains

  !===================================================

  !> Reads the namelist file for user-specified control 
  !! parameters
  SUBROUTINE read_namelist(m, n, itmax)
    IMPLICIT none
    INTEGER, INTENT(out) :: m, n
    INTEGER, INTENT(out) :: itmax
    ! namelist input 
    CHARACTER (LEN=8) :: nml_name = "namelist" 
    INTEGER :: input_unit = 99
    INTEGER :: ierr

    NAMELIST/global_domain/ m, n, itmax, mprint
    NAMELIST/io_control/ l_out

    ! Initialise these vars to problem values as defense
    ! against failure to read them properly and also
    ! to squelch incorrect compiler warnings about them
    ! not being assigned to.
    m = -1
    n = -1
    itmax = 0

    !     Read in namelist 
    OPEN(unit=input_unit, file=nml_name, status='old',iostat=ierr)

    CALL check(ierr, "open "//nml_name)
    READ(unit=input_unit, nml=global_domain, iostat=ierr)
    CALL check(ierr, "read "//nml_name)
    READ(unit=input_unit, nml=io_control, iostat=ierr)
    CALL check(ierr, "read "//nml_name)

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
    implicit none
    type(grid_type), intent(in) :: grid
    integer, intent(in) :: istp
    real(wp), dimension(:,:), intent(in) :: ht, sshn, un, vn
    ! Locals
    integer :: ji, jj
    real(wp) :: rtmp1, rtmp2
    character(len=5) :: fname
    
    if( l_out .and. (mod(istp, mprint) .eq. 0) ) then

       ! output model results
       write(fname, '(I5.5)') istp
       !OPEN(21, file='go2d_'//fname//'.dat', STATUS='UNKNOWN')
       open(21, file='go2d_'//fname//'.csv', STATUS='UNKNOWN')
       rewind(21)

       DO jj = 1, jpj
          DO ji = 1, jpi
             rtmp1 = 0.5_wp * (un(ji-1,jj) + un(ji,jj))
             rtmp2 = 0.5_wp * (vn(ji,jj-1) + vn(ji,jj))

             ! write "x-coord, y-coord, depth, ssh, u-velocity, v-velocity" to ASCII files

             !WRITE(1,'(2f20.3, 2f15.4, 2e18.3)')  &            
             WRITE(21,'(f20.3,'','',f20.3,'','',f15.4,'','',f15.4,'','',f18.3,'','',f18.3)') &
             & grid%xt(ji,jj), grid%yt(ji,jj), ht(ji,jj), sshn(ji,jj),rtmp1, rtmp2 
          END DO
       END DO
          
       close(21)

    end if ! we're doing output

  end subroutine model_write

  !===================================================

  !> Write log entry with one integer and one real arg
  SUBROUTINE write_log_ir(fmtstr, istep, fvar)
    IMPLICIT none
    CHARACTER(LEN=*), INTENT(in) :: fmtstr
    INTEGER,          INTENT(in) :: istep
    REAL(wp),         INTENT(in) :: fvar

    WRITE(6,FMT=fmtstr) istep, fvar

  END SUBROUTINE write_log_ir

  !===================================================

  !> Write log entry with one real arg
  SUBROUTINE write_log_r(fmtstr, fvar)
    IMPLICIT none
    CHARACTER(LEN=*), INTENT(in) :: fmtstr
    REAL(wp),         INTENT(in) :: fvar

    WRITE(6,FMT=fmtstr) fvar

  END SUBROUTINE write_log_r

  !===================================================

  SUBROUTINE model_write_finalise()
    IMPLICIT none
    
    ! Nothing to do here as we open and close ASCII file for
    ! each write.

  END SUBROUTINE model_write_finalise

  !===================================================

  ! Check error code
  subroutine check(status, text)
    implicit none
      
    integer, intent(in) :: status
    character (len=*)   :: text
    
    if (status /= 0) then
       write(6,*) "error ", status
       write(6,*) text
       stop 2
    endif

  end subroutine check

  !===================================================

END MODULE gocean2d_io_mod
