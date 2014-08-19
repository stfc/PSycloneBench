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

     !> This kernel is written assuming that the internal
     !! regions for each grid-point type begin at the
     !! following locations on the grid:
     integer :: tstart(2) = (/ 1, 1 /)
     integer :: ustart(2) = (/ 2, 1 /)
     integer :: vstart(2) = (/ 1, 2 /)
     integer :: fstart(2) = (/ 2, 2 /)

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
    integer :: kern_uxshift, kern_uyshift, kern_vxshift, kern_vyshift, kern_txshift, kern_tyshift
    integer, dimension(2) :: ushift, vshift, tshift
    type(compute_z_type) :: this_kernel

    dx = zfld%grid%dx
    dy = zfld%grid%dy

    ! In original shallow, F field internal points began at (2,2), U internal
    ! points at (2,1), V internal points at (1,2) and T internal pts at 1,1.
    !                    F      U      V       T
    !                  (2,2)  (2, 1) (1, 2)  (1,1)
    ! Original shift          (0,-1) (-1,0)  (-1,-1)

    ! The shifts from F to xx pts assumed by this kernel
    kern_uxshift = this_kernel%ustart(1) - this_kernel%fstart(1)
    kern_uyshift = this_kernel%ustart(2) - this_kernel%fstart(2)
    kern_vxshift = this_kernel%vstart(1) - this_kernel%fstart(1)
    kern_vyshift = this_kernel%vstart(2) - this_kernel%fstart(2)
    kern_txshift = this_kernel%tstart(1) - this_kernel%fstart(1)
    kern_tyshift = this_kernel%tstart(2) - this_kernel%fstart(2)

    ! The shifts we must apply to account for the fact that our
    ! fields may not have the same relative positioning as
    ! assumed by the kernel
    ushift(1) = ufld%internal%xstart - zfld%internal%xstart - kern_uxshift
    ushift(2) = ufld%internal%ystart - zfld%internal%ystart - kern_uyshift
    vshift(1) = vfld%internal%xstart - zfld%internal%xstart - kern_vxshift
    vshift(2) = vfld%internal%ystart - zfld%internal%ystart - kern_vyshift
    tshift(1) = pfld%internal%xstart - zfld%internal%xstart - kern_txshift
    tshift(2) = pfld%internal%ystart - zfld%internal%ystart - kern_tyshift

    do J=zfld%internal%ystart, zfld%internal%ystop, 1
       do I=zfld%internal%xstart, zfld%internal%xstop, 1

          call compute_z_code(i, j, ushift, vshift, tshift, dx, dy, &
                              zfld%data, &
                              pfld%data, &
                              ufld%data, &
                              vfld%data)

       end do
    end do

  end subroutine manual_invoke_compute_z

  !===================================================

  !> Compute the potential vorticity on the grid point (i,j)
  subroutine compute_z_code(i, j, ushift, vshift, tshift, dx, dy, z, p, u, v)
    implicit none
    integer,  intent(in) :: I, J
    integer,  intent(in),  dimension(2) :: ushift, vshift, tshift
    real(wp), intent(in) :: dx, dy
    real(wp), intent(out), dimension(:,:) :: z
    real(wp), intent(in),  dimension(:,:) :: p, u, v

    Z(I,J) =((4.0/dx)*( &
            V(I+vshift(1),J+vshift(2))-V(I-1+vshift(1),J+vshift(2)))- &
            (4.0/dy)*(  &
            U(I+ushift(1),J+ushift(2))-U(I+ushift(1),J-1+ushift(2))))/ &
            (P(I-1+tshift(1),J-1+tshift(2))+P(I+tshift(1),J-1+tshift(2))+&
             P(I+tshift(1),J+tshift(2))+P(I-1+tshift(1),J+tshift(2)))

  end subroutine compute_z_code

end module compute_z_mod
