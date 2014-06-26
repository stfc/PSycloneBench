MODULE compute_vnew_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use field_mod
  IMPLICIT none

  PRIVATE

  PUBLIC manual_invoke_compute_vnew
  PUBLIC compute_vnew_type, compute_vnew_code

  TYPE, EXTENDS(kernel_type) :: compute_vnew_type
     TYPE(arg), DIMENSION(6) :: meta_args =    &
          (/ arg(WRITE, CV, POINTWISE),        & ! vnew
             arg(READ,  CV, POINTWISE),        & ! vold
             arg(READ,  CF, POINTWISE),        & ! z
             arg(READ,  CU, POINTWISE),        & ! cu
             arg(READ,  CT, POINTWISE),        & ! h
             arg(READ,  R,  POINTWISE)         & ! tdt
           /)
     !> This kernel operates on fields that live on an
     !! orthogonal, regular grid.
     integer :: GRID_TYPE = ORTHOGONAL_REGULAR

     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     INTEGER :: ITERATES_OVER = DOFS
  CONTAINS
    procedure, nopass :: code => compute_vnew_code
  END TYPE compute_vnew_type

CONTAINS

  !===================================================

  subroutine manual_invoke_compute_vnew(vnew, vold, z, cu, h, tdt)
    implicit none
    type(r2d_field_type), intent(out) :: vnew
    type(r2d_field_type), intent(in)  :: vold, z, cu, h
    real(wp), intent(in) :: tdt
    ! Locals
    integer :: I, J
    real(wp) :: dx, dy

    ! Note that we do not loop over the full extent of the field.
    ! Fields are allocated with extents (M+1,N+1).
    ! Presumably the extra row and column are needed for periodic BCs.
    ! We are updating a quantity on CU.
    ! This loop writes to vnew(1:M,2:N+1) so this looks like
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  x  x  x  o 
    !  x  x  x  o   j=N
    !  x  x  x  o
    !  o  o  o  o   j=1

    ! Original code looked like:
    !
    ! DO J=1,N
    !   DO I=1,M
    !     VNEW(I,J+1) = VOLD(I,J+1)-TDTS8*(Z(I+1,J+1)+Z(I,J+1)) & 
    !         *(CU(I+1,J+1)+CU(I,J+1)+CU(I,J)+CU(I+1,J))        & 
    !         -TDTSDY*(H(I,J+1)-H(I,J))
    !   END DO
    ! END DO

    ! vnew(i,j) depends upon:
    !   vold(i,j)                                : CV
    !   z(i+1,j), z(i,j)                         : CF
    !    => lateral CF neighbours of the CV pt being updated
    !   cu(i,j), cu(i+1,j),cu(i,j-1),cu(i+1,j-1) : CU
    !    => all CU neighbours of the CV pt being updated 
    !   h(i,j),   h(i,j-1)                       : CT
    !    =>  vertical CT neighbours of the CV pt being updated

    !   x-------x-------fi+1j+1
    !   |       |       |
    !   |       |       |
    !   uij-----Tij-----ui+1j
    !   |       |       |
    !   |       |       |
    !   fij-----vij-----fi+1j
    !   |       |       |
    !   |       |       |
    !   uij-1- -Tij-1---ui+1j-1
    !

    dx = vnew%grid%dx
    dy = vnew%grid%dy

    DO J=vnew%internal%ystart, vnew%internal%ystop, 1
       DO I=vnew%internal%xstart, vnew%internal%xstop, 1

          CALL compute_vnew_code(i, j, dx, dy, &
                                 vnew%data, vold%data, &
                                 z%data, cu%data, h%data, tdt)
       END DO
    END DO

  END SUBROUTINE manual_invoke_compute_vnew

  !===================================================

  subroutine compute_vnew_code(i, j, dx, dy, vnew, vold, z, cu, h, tdt)
    implicit none
    integer,  intent(in) :: I, J
    real(wp), intent(in) :: dx, dy
    REAL(wp), INTENT(out), DIMENSION(:,:) :: vnew
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: vold, z, cu, h
    REAL(wp), INTENT(in) :: tdt
    ! Locals
    REAL(wp) :: tdts8, tdtsdy
   
    !> These quantities are computed here because tdt is not
    !! constant. (It is == dt for first time step, 2xdt for
    !! all remaining time steps.)
    tdts8 = tdt/8.0d0
    tdtsdy = tdt/dy

    VNEW(I,J) = VOLD(I,J)-TDTS8*(Z(I+1,J)+Z(I,J))                 & 
                *(CU(I+1,J)+CU(I,J)+CU(I,J-1)+CU(I+1,J-1))        & 
                 -TDTSDY*(H(I,J)-H(I,J-1))

  END SUBROUTINE compute_vnew_code

END MODULE compute_vnew_mod
