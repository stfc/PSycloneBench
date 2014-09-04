MODULE time_smooth_mod
  USE kind_params_mod
  USE kernel_mod
  use argument_mod
  IMPLICIT none

  PRIVATE

  PUBLIC time_smooth_init, invoke_time_smooth
  PUBLIC time_smooth_type, time_smooth_code

  !> Parameter for time smoothing
  REAL(wp) :: alpha

  !> The time smoothing operates in time rather than space
  !! and therefore takes three fields defined on any one
  !! of the four grid point types (T, U, V or Q).
  !! Presumably FE should be FD for us and maybe CELLS 
  !! should be COLUMNS?
  TYPE, EXTENDS(kernel_type) :: time_smooth_type
     TYPE(arg), DIMENSION(3) :: meta_args = &
          (/ arg(READ, EVERY, POINTWISE),     &
             arg(READ, EVERY, POINTWISE),     &
             arg(READWRITE , EVERY, POINTWISE)      &
           /)
     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     INTEGER :: ITERATES_OVER = DOFS
  CONTAINS
    procedure, nopass :: code => time_smooth_code
  END type time_smooth_type

CONTAINS

  !===================================================

  !> Initialise the time-smoothing module. Sets parameter
  !! alpha that is used in the time-smooth kernel.
  SUBROUTINE time_smooth_init(alpha_tmp)
    IMPLICIT none
    REAL(wp), INTENT(in) :: alpha_tmp

    alpha = alpha_tmp

  END SUBROUTINE time_smooth_init

  !===================================================

  !> Manual implementation of code to invoke the time-smoothing
  !! kernel
  SUBROUTINE invoke_time_smooth(field1, field1_new, field1_old, &
                                field2, field2_new, field2_old, &
                                field3, field3_new, field3_old)
    use topology_mod, only: M, N
    IMPLICIT none
    REAL(wp), INTENT(in), DIMENSION(:,:) :: field1, field2, field3
    REAL(wp), INTENT(in), DIMENSION(:,:) :: field1_new, field2_new, field3_new
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: field1_old, field2_old, field3_old
    ! Locals
    INTEGER :: i, j

    ! The time-smoothing is applied to a field at *every* grid point
    
    ! Loop over 'columns'
    DO J=1,N+1 !idim2
      DO I=1,M+1 !idim1
        CALL time_smooth_code(i,j,field1,field1_new,field1_old)
      END DO
    END DO

    ! Loop over 'columns'
    DO J=1,N+1 ! idim2
      DO I=1,M+1 ! idim1
         CALL time_smooth_code(i,j,field2,field2_new,field2_old)
      END DO
    END DO

    ! Loop over 'columns'
    DO J=1,N+1 ! idim2
      DO I=1,M+1 ! idim1
         CALL time_smooth_code(i,j,field3,field3_new,field3_old)
      END DO
    END DO

  END SUBROUTINE invoke_time_smooth

  !===================================================

  !> Kernel to smooth supplied field in time
  SUBROUTINE time_smooth_code(i, j, field, field_new, field_old)
    IMPLICIT none
    INTEGER,  INTENT(in)                    :: i, j
    REAL(wp), INTENT(in),    DIMENSION(:,:) :: field
    REAL(wp), INTENT(in),    DIMENSION(:,:) :: field_new
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: field_old

    field_old(i,j) = field(i,j) + &
         alpha*(field_new(i,j) - 2.*field(i,j) + field_old(i,j))

  END SUBROUTINE time_smooth_code

END MODULE time_smooth_mod
