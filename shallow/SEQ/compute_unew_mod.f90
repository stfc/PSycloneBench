MODULE compute_unew_mod
  USE kind_params_mod
  USE kernel_mod
  use argument_mod
  use field_mod
  implicit none

  private

  public manual_invoke_compute_unew
  public compute_unew_type, compute_unew_code

  type, extends(kernel_type) :: compute_unew_type
     type(arg), dimension(6) :: meta_args =    &
          (/ arg(WRITE, CU, POINTWISE),        & ! unew
             arg(READ,  CU, POINTWISE),        & ! uold
             arg(READ,  CF, POINTWISE),        & ! z
             arg(READ,  CV, POINTWISE),        & ! cv
             arg(READ,  CT, POINTWISE),        & ! h
             arg(READ,  R,  POINTWISE)         & ! tdt
           /)
     !> This kernel operates on fields that live on an
     !! orthogonal, regular grid.
     integer :: GRID_TYPE = ORTHOGONAL_REGULAR

     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS
  contains
    procedure, nopass :: code => compute_unew_code
  end type compute_unew_type

contains

  !===================================================

  subroutine manual_invoke_compute_unew(unew, uold, z, cv, h, tdt)
    implicit none
    type(r2d_field_type), intent(inout) :: unew
    type(r2d_field_type), intent(in)    :: uold, z, cv, h
    real(wp), intent(in) :: tdt
    ! Locals
    integer  :: I, J
    real(wp) :: dx, dy

    ! Note that we do not loop over the full extent of the field.
    ! Fields are allocated with extents (M+1,N+1).
    ! Presumably the extra row and column are needed for periodic BCs.
    ! We are updating a quantity on CU.
    ! This loop writes to unew(2:M+1,1:N) so this looks like
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  x  x  x   j=1

    ! unew(i,j) depends upon:
    !   uold(i,j)
    !   z(i,j+1),  z(i,j)
    !  cv(i,j),   cv(i,j+1), cv(i-1,j+1), cv(i-1,j)
    !   h(i,j),    h(i-1,j)

    ! Swap indices, e.g. XX(i+1,j) => YY(i,j+1)
    ! Any field on U replaced with field on V
    ! => produces same code for the update of corresponding field on V.

    ! Original code looked like:
    !
    ! DO J=1,N
    !   DO I=1,M
    !     UNEW(I+1,J) = UOLD(I+1,J)+                                     &
    !         TDTS8*(Z(I+1,J+1)+Z(I+1,J))*(CV(I+1,J+1)+CV(I,J+1)+CV(I,J) &
    !        +CV(I+1,J))-TDTSDX*(H(I+1,J)-H(I,J))                       
    !   END DO
    ! END DO
    dx = unew%grid%dx
    dy = unew%grid%dy

    DO J=unew%internal%ystart, unew%internal%ystop, 1
       DO I=unew%internal%xstart, unew%internal%xstop, 1

          CALL compute_unew_code(i, j, dx, dy, &
                                 unew%data, uold%data, &
                                 z%data, cv%data, h%data, tdt)
       END DO
    END DO

  END SUBROUTINE manual_invoke_compute_unew

  !===================================================

  subroutine compute_unew_code(i, j, dx, dy, unew, uold, z, cv, h, tdt)
    implicit none
    integer,  intent(in) :: I, J
    real(wp), intent(in) :: dx, dy
    real(wp), intent(out), dimension(:,:) :: unew
    real(wp), intent(in),  dimension(:,:) :: uold, z, cv, h
    real(wp), intent(in) :: tdt
    ! Locals
    real(wp) :: tdts8, tdtsdx
   
    !> These quantities are computed here because tdt is not
    !! constant. (It is == dt for first time step, 2xdt for
    !! all remaining time steps.)
    tdts8 = tdt/8.0d0
    tdtsdx = tdt/dx

    UNEW(I,J) = UOLD(I,J) +                                 &
                TDTS8*(Z(I,J+1)+Z(I,J)) *                   &
                (CV(I,J+1)+CV(I-1,J+1)+CV(I-1,J)+CV(I,J)) - &
                TDTSDX*(H(I,J)-H(I-1,J))

  end subroutine compute_unew_code

end module compute_unew_mod