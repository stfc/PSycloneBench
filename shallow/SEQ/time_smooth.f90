MODULE time_smooth
  IMPLICIT none

  PRIVATE

  !> Parameter for time smoothing
  REAL(KIND=8) :: alpha

  PUBLIC time_smooth_init, manual_invoke_time_smooth

CONTAINS

  !===================================================

  SUBROUTINE time_smooth_init(alpha_tmp)
    IMPLICIT none
    REAL(KIND=8), INTENT(in) :: alpha_tmp

    alpha = alpha_tmp

  END SUBROUTINE time_smooth_init

  !===================================================

  SUBROUTINE manual_invoke_time_smooth(field, field_new, field_old)
    IMPLICIT none
    REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: field
    REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: field_new
    REAL(KIND=8), INTENT(inout), DIMENSION(:,:) :: field_old
    ! Locals
    INTEGER :: i, j
    INTEGER :: idim1, idim2
    
    idim1 = SIZE(field, 1)
    idim2 = SIZE(field, 2)

    DO J=1,idim2
      DO I=1,idim1
        field_old(I,J) = field(I,J)+ &
             ALPHA*(field_new(I,J)-2.*field(I,J)+field_old(I,J))
      END DO
    END DO

  END SUBROUTINE manual_invoke_time_smooth

END MODULE time_smooth
