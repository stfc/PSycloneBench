!> \brief Compute the potential vorticity, z
!! \detail Given the current pressure and velocity fields,
!! computes the potential voriticity.
module compute_z_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use field_mod
  implicit none

  private

  public manual_invoke_compute_z
  public compute_z_type, compute_z_code

  type, extends(kernel_type) :: compute_z_type
     type(arg), dimension(4) :: meta_args =    &
          (/ arg(WRITE, CF, POINTWISE),        & ! z
             arg(READ,  CT, POINTWISE),        & ! p
             arg(READ,  CU, POINTWISE),        & ! u
             arg(READ,  CV, POINTWISE)         & ! v
           /)
     !> This kernel operates on fields that live on an
     !! orthogonal, regular grid.
     integer :: GRID_TYPE = ORTHOGONAL_REGULAR

     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS
  contains
    procedure, nopass :: code => compute_z_code
  end type compute_z_type

contains

  !===================================================

  !> Manual implementation of the code needed to invoke
  !! compute_z_code().
  subroutine manual_invoke_compute_z(zfld, pfld, ufld, vfld)
    implicit none
    type(r2d_field_type), intent(inout) :: zfld
    type(r2d_field_type), intent(in)    :: pfld, ufld, vfld
    ! Locals
    integer :: I, J
    real(wp) :: dx, dy
    integer :: uxshift, uyshift, vxshift, vyshift, txshift, tyshift
    real(wp), dimension(:,:), pointer :: u, v, p, z

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

    dx = zfld%grid%dx
    dy = zfld%grid%dy

    ! In original shallow, F field internal points began at (2,2), U internal
    ! points at (2,1), V internal points at (1,2) and T internal pts at 1,1.
    !                    F      U      V       T
    !                  (2,2)  (2, 1) (1, 2)  (1,1)
    ! Original shift          (0,-1) (-1,0)  (-1,-1)
    !             (fxstart,fystart) (uxstart,uystart)
    ! If we're looping over F internal pts and kernel assumes a shift
    ! of (0,-1) to U internal pts then we need to shift by:
    ! iu = if - (fxstart - uxstart - 0)
    ! ju = jf - (fystart - uystart - 1)
    ! For V pts the kernel assumes a shift of (-1,0) relative to F:
    ! iv = if - (fxstart - vxstart - 1)
    ! jv = jf - (fystart - vystart - 0)
    ! For T pts the kernel assumes a shift of (-1,-1) relative to T:
    ! it = if - (fxstart - txstart - 1)
    ! jt = jf - (fystart - tystart - 1)
    uxshift = ufld%internal%xstart - zfld%internal%xstart
    uyshift = ufld%internal%ystart - zfld%internal%ystart + 1
    vxshift = vfld%internal%xstart - zfld%internal%xstart + 1
    vyshift = vfld%internal%ystart - zfld%internal%ystart
    txshift = pfld%internal%xstart - zfld%internal%xstart + 1
    tyshift = pfld%internal%ystart - zfld%internal%ystart + 1

    do J=zfld%internal%ystart, zfld%internal%ystop, 1
       do I=zfld%internal%xstart, zfld%internal%xstop, 1

!          call compute_z_code(i, j, dx, dy,         &
!                              zfld%data, pfld%data, &
!                              ufld%data, vfld%data)
          Zfld%data(I,J) =((4.0/dx)*(Vfld%data(i+vxshift,J+vyshift)-Vfld%data(i+vxshift-1,J+vyshift)) - &
                   (4.0/dy)*(Ufld%data(i+uxshift,J+uyshift)-Ufld%data(i+uxshift,J-1+uyshift)))/ &
                 (Pfld%data(I-1+txshift,J-1+tyshift)+ &
                  Pfld%data(I+txshift,J-1+tyshift)+ &
                  Pfld%data(I+txshift,J+tyshift)+ &
                  Pfld%data(I-1+txshift,J+tyshift))

       end do
    end do

  end subroutine manual_invoke_compute_z

  !===================================================

  !> Compute the potential vorticity on the grid point (i,j)
  subroutine compute_z_code(i, j, dx, dy, z, p, u, v)
    implicit none
    integer,  intent(in) :: I, J
    real(wp), intent(in) :: dx, dy
    real(wp), intent(out), dimension(:,:) :: z
    real(wp), intent(in),  dimension(:,:) :: p, u, v

    Z(I,J) =((4.0/dx)*(V(I,J)-V(I-1,J))-(4.0/dy)*(U(I,J)-U(I,J-1)))/ &
                 (P(I-1,J-1)+P(I,J-1)+P(I,J)+P(I-1,J))

  end subroutine compute_z_code

end module compute_z_mod
