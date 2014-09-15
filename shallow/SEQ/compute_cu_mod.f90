!> \brief Compute the mass flux in the x direction, cu
!! \detail Given the current pressure and velocity fields,
!! computes the mass flux in the x direction.
module compute_cu_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use field_mod
  implicit none

  private

  public manual_invoke_compute_cu
  public compute_cu_type, compute_cu_code

  type, extends(kernel_type) :: compute_cu_type
     type(arg), dimension(3) :: meta_args =    &
          (/ arg(WRITE, CU, POINTWISE),        & ! cu
             arg(READ,  CT, POINTWISE),        & ! p
             arg(READ,  CU, POINTWISE)         & ! u
           /)
     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS

     !> This kernel is written assuming that the *arrays*
     !! containing the fields that it uses are indexed relative to the
     !! T point in the following way: If there is a T point at index
     !! (i,j) then the U point with the same grid index is at
     !! (i+ushift(1),j+ushift(2)).
     integer :: ushift(2) = (/ -1,  0 /)
     integer :: vshift(2) = (/  0, -1 /)
     integer :: fshift(2) = (/ -1, -1 /)

!     integer :: index_offset = SW_OFFSET

  contains
    procedure, nopass :: code => compute_cu_code
  end type compute_cu_type

contains

  !===================================================

  !> Manual implementation of the code needed to invoke
  !! compute_cu_code().
  subroutine manual_invoke_compute_cu(cufld, pfld, ufld)
    implicit none
    type(r2d_field_type), intent(inout) :: cufld
    type(r2d_field_type), intent(in)    :: pfld, ufld
    ! Locals
    integer :: I, J
    integer, dimension(2) :: kern_tshift
    integer, dimension(2) :: tshift
    type(compute_cu_type) :: this_kernel

    ! Note that we do not loop over the full extent of the field.
    ! Fields are allocated with extents (M+1,N+1).
    ! Presumably the extra row and column are needed for periodic BCs.
    ! We are updating a quantity on CU.
    ! This loop writes to cu(2:M+1,1:N) so this looks like
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  x  x  x   j=1

    ! Quantity CU is mass flux in x direction.

    ! Original code looked like:
    !
    !    DO J=1,N
    !      DO I=1,M
    !           CU(I+1,J) = .5*(P(I+1,J)+P(I,J))*U(I+1,J)
    !      END DO
    !    END DO

    ! cu(i,j) depends upon:
    !   p(i-1,j), p(i,j) : CT
    !    => lateral CT neighbours of the CU pt being updated
    !   u(i,j)           : CU
    !    => the horiz. vel. component at the CU pt being updated

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
    ! The shifts from U to xx pts assumed by this kernel
!    kern_tshift(1) = this_kernel%tstart(1) - this_kernel%ustart(1)
!    kern_tshift(2) = this_kernel%tstart(2) - this_kernel%ustart(2)

    ! The shifts we must apply to account for the fact that our
    ! fields may not have the same relative positioning as
    ! assumed by the kernel
!    tshift(1) = cufld%internal%xstart - pfld%internal%xstart - kern_tshift(1)
!    tshift(2) = cufld%internal%ystart - pfld%internal%ystart - kern_tshift(2)

    do J=cufld%internal%ystart, cufld%internal%ystop
       do I=cufld%internal%xstart, cufld%internal%xstop

          call compute_cu_code(i, j, cufld%data, pfld%data, ufld%data)
       end do
    end do

  end subroutine manual_invoke_compute_cu

  !===================================================

  !> Compute the mass flux in the x direction at point (i,j)
  subroutine compute_cu_code(i, j, cu, p, u)
    implicit none
    integer,  intent(in) :: I, J
    real(wp), intent(out), dimension(:,:) :: cu
    real(wp), intent(in),  dimension(:,:) :: p, u

    CU(I,J) = 0.5d0*(P(i+1,J)+P(I,J))*U(I,J)

  end subroutine compute_cu_code

end module compute_cu_mod
