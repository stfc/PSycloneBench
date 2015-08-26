!> Module containing basic utilities
module gocean_mod
  use kind_params_mod
  implicit none

  !> Interface to logging routines
  interface model_write_log
     module procedure write_log_a, write_log_ir, &
                      write_log_i, write_log_r
  end interface

contains

  !===================================================

  !> Initialise the GOcean environment
  subroutine gocean_init()
#if _OPENACC
    use openacc
#endif
    implicit none

#if _OPENACC
    call acc_init(acc_device_nvidia)
#endif
  end subroutine gocean_init

  !===================================================

  !> Stop the model run. Currently simply does
  !! a Fortran STOP.
  !! @param[in] msg Message to print - reason we're stopping
  subroutine gocean_stop(msg)
    use iso_fortran_env, only : error_unit ! access computing environment
    implicit none
    character(len=*), intent(in) :: msg

    write(error_unit, *) msg
    stop 1

  end subroutine gocean_stop

  !===================================================

  !> Write log entry with one integer and one real arg
  SUBROUTINE write_log_ir(fmtstr, istep, fvar)
    use iso_fortran_env, only : output_unit ! access computing environment
    IMPLICIT none
    CHARACTER(LEN=*), INTENT(in) :: fmtstr
    INTEGER,          INTENT(in) :: istep
    REAL(wp),         INTENT(in) :: fvar

    WRITE(output_unit,FMT=fmtstr) istep, fvar

  END SUBROUTINE write_log_ir

  !===================================================

  !> Write log entry with one integer arg
  SUBROUTINE write_log_i(fmtstr, istep)
    use iso_fortran_env, only : output_unit ! access computing environment
    IMPLICIT none
    CHARACTER(LEN=*), INTENT(in) :: fmtstr
    INTEGER,          INTENT(in) :: istep

    WRITE(output_unit,FMT=fmtstr) istep

  END SUBROUTINE write_log_i

  !===================================================

  !> Write log entry with one real arg
  subroutine write_log_r(fmtstr, fvar)
    use iso_fortran_env, only : output_unit ! access computing environment
    implicit none
    character(len=*), intent(in) :: fmtstr
    real(wp),         intent(in) :: fvar

    write(output_unit,FMT=fmtstr) fvar

  end subroutine write_log_r

  !===================================================

  !> Write log entry with just a string
  subroutine write_log_a(fmtstr, msg)
    use iso_fortran_env, only : output_unit ! access computing environment
    implicit none
    character(len=*), intent(in) :: fmtstr
    character(len=*), intent(in) :: msg

    write(output_unit,FMT=fmtstr) msg

  end subroutine write_log_a

  !===================================================

end module gocean_mod
