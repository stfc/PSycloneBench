module kernel_args_mod
  use opencl_utils_mod
  implicit none

  !> TODO use GOcean infrastructure?
  integer, parameter :: wp = kind(1.0d0)

contains
  
  subroutine set_continuity_args(kernel, &
                         nx,             &
                         ssha_device,    &
                         sshn_device,    &
                         sshn_u_device,  &
                         sshn_v_device,  &
                         hu_device,      &
                         hv_device,&
                         un_device,&
                         vn_device,&
                         rdt,&
                         e12t_device)    
    integer(c_intptr_t), target :: kernel
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: ssha_device
    integer(c_intptr_t), target :: sshn_device, sshn_u_device, sshn_v_device
    integer(c_intptr_t), target :: un_device, vn_device
    integer(c_intptr_t), target :: hu_device, hv_device
    integer(c_intptr_t), target :: e12t_device
    real(kind=wp), target :: rdt
    ! Locals
    integer :: arg_idx, ierr

    ! set the kernel arguments
    arg_idx = 0
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(ssha_device), &
         C_LOC(ssha_device))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(sshn_device), &
         C_LOC(sshn_device))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(sshn_u_device), &
         C_LOC(sshn_u_device))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(sshn_v_device), &
         C_LOC(sshn_v_device))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(hu_device), &
         C_LOC(hu_device))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(hv_device), &
         C_LOC(hv_device))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(un_device), &
         C_LOC(un_device))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(vn_device), &
         C_LOC(vn_device))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(rdt), C_LOC(rdt))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kernel, arg_idx, sizeof(e12t_device), &
         C_LOC(e12t_device))
    call check_status("clSetKernelArg", ierr)

  end subroutine set_continuity_args

end module kernel_args_mod
