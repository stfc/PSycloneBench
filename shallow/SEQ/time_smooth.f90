MODULE time_smooth
  IMPLICIT none

  PRIVATE

  !> Parameter for time smoothing
  REAL(KIND=8) :: alpha

  PUBLIC time_smooth_init, manual_invoke_time_smooth

CONTAINS

  !===================================================

  !> Initialise the time-smoothing module. Sets parameter
  !! alpha that is used in the time-smooth kernel.
  SUBROUTINE time_smooth_init(alpha_tmp)
    IMPLICIT none
    REAL(KIND=8), INTENT(in) :: alpha_tmp

    alpha = alpha_tmp

  END SUBROUTINE time_smooth_init

  !===================================================

  !> Manual implementation of code to invoke the time-smoothing
  !! kernel
  SUBROUTINE manual_invoke_time_smooth(field, field_new, field_old)
    IMPLICIT none
    REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: field
    REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: field_new
    REAL(KIND=8), INTENT(inout), DIMENSION(:,:) :: field_old
    ! Locals
    INTEGER :: i, j
    INTEGER :: idim1, idim2
    
    ! Here we will query what should be field objects to get at
    ! raw data.
    idim1 = SIZE(field, 1)
    idim2 = SIZE(field, 2)

    ! Loop over 'columns'
    DO J=1,idim2
      DO I=1,idim1
         CALL time_smooth_code(i,j,field,field_new,field_old)
      END DO
    END DO

  END SUBROUTINE manual_invoke_time_smooth

  !===================================================

  !> Kernel to smooth supplied field in time
  SUBROUTINE time_smooth_code(i, j, field, field_new, field_old)
    IMPLICIT none
    INTEGER,      INTENT(in)                    :: i, j
    REAL(KIND=8), INTENT(in),    DIMENSION(:,:) :: field
    REAL(KIND=8), INTENT(in),    DIMENSION(:,:) :: field_new
    REAL(KIND=8), INTENT(inout), DIMENSION(:,:) :: field_old

    field_old(i,j) = field(i,j) + &
         alpha*(field_new(i,j) - 2.*field(i,j) + field_old(i,j))

  END SUBROUTINE time_smooth_code

END MODULE time_smooth
