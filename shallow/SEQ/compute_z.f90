MODULE compute_z
  USE kind_params
  USE kernel_mod
  use argument_mod
  IMPLICIT none

  PRIVATE

  PUBLIC manual_invoke_compute_z
  PUBLIC compute_z_type, compute_z_code

  TYPE, EXTENDS(kernel_type) :: compute_z_type
     TYPE(arg), DIMENSION(6) :: meta_args =    &
          (/ arg(WRITE, CT, POINTWISE),        & ! z
             arg(READ,  CT, POINTWISE),        & ! p
             arg(READ,  CU, POINTWISE),        & ! u
             arg(READ,  CV, POINTWISE)         & ! v
           /)
     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     INTEGER :: ITERATES_OVER = DOFS
  CONTAINS
    procedure, nopass :: code => compute_z_code
  END TYPE compute_z_type

CONTAINS

  !===================================================

  SUBROUTINE manual_invoke_compute_z(z, p, u, v)
    IMPLICIT none
    REAL(wp), INTENT(out), DIMENSION(:,:) :: z
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: p, u, v
    ! Locals
    INTEGER :: I, J

    ! Note that we do not loop over the full extent of the field.
    ! Fields are allocated with extents (M+1,N+1).
    ! Presumably the extra row and column are needed for periodic BCs.
    ! We are updating a quantity on CF.
    ! This loop writes to z(2:M+1,2:N+1) so this looks like
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  x  x  x 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  o  o  o   j=1

    ! Quantity Z is Potential Enstrophy

    ! Original code looked like:
    !
    !   DO J=1,N
    !      DO I=1,M
    !         Z(I+1,J+1) =(FSDX*(V(I+1,J+1)-V(I,J+1))-FSDY*(U(I+1,J+1) & 
    !                -U(I+1,J)))/(P(I,J)+P(I+1,J)+P(I+1,J+1)+P(I,J+1))
    !      END DO
    !    END DO

    ! z(i,j) depends upon:
    !   p(i-1,j-1), p(i,j-1), p(i-1,j), p(i,j) : CT
    !    => all CT neighbours of the CF pt being updated
    !   u(i,j), u(i,j-1)                      : CU
    !    => vertical CU neighbours of the CF pt being updated 
    !   v(i,j), v(i-1,j)                      : CV
    !    => lateral CV neighbours of the CF pt being updated

    !   vi-1j+1--fij+1---vij+1---fi+1j+1
    !   |        |       |       |
    !   |        |       |       |
    !   Ti-1j----uij-----Tij-----ui+1j
    !   |        |       |       |
    !   |        |       |       |
    !   vi-1j----fij-----vij-----fi+1j
    !   |        |       |       |
    !   |        |       |       |
    !   Ti-1j-1--uij-1---Tij-1---ui+1j-1
    !

    DO J=2, SIZE(z, 2)
       DO I=2, SIZE(z, 1)

          CALL compute_z_code(i, j, z, p, u, v)
       END DO
    END DO

  END SUBROUTINE manual_invoke_compute_z

  !===================================================

  SUBROUTINE compute_z_code(i, j, z, p, u, v)
    USE mesh, ONLY: fsdx, fsdy
    IMPLICIT none
    INTEGER, INTENT(in) :: I, J
    REAL(wp), INTENT(out), DIMENSION(:,:) :: z
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: p, u, v

    Z(I,J) =(FSDX*(V(I,J)-V(I-1,J))-FSDY*(U(I,J)-U(I,J-1)))/ &
                 (P(I-1,J-1)+P(I,J-1)+P(I,J)+P(I-1,J))

  END SUBROUTINE compute_z_code

END MODULE compute_z
