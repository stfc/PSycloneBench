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

    write(*,"('Set ',I2,' arguments for Continuity kernel')") arg_idx

  end subroutine set_continuity_args

  subroutine set_momu_args(kern,      &
		           nx,        &
		           ua_device, &
		           un_device, &
		           vn_device, &
		           hu_device, &
		           hv_device, &
		           ht_device, &
		           ssha_u_device, &
		           sshn_device,   &
		           sshn_u_device, &
		           sshn_v_device, &
		           tmask_device, &
		           e1u_device, e1v_device, &
		           e1t_device, e2u_device, &
		           e2t_device, e12u_device,&
		           gphiu_device,           &
		           rdt, cbfr, visc)
    integer(c_intptr_t), target :: kern
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: ssha_device, ssha_u_device
    integer(c_intptr_t), target :: sshn_device, sshn_u_device, sshn_v_device
    integer(c_intptr_t), target :: ua_device, un_device, vn_device
    integer(c_intptr_t), target :: hu_device, hv_device, ht_device
    integer(c_intptr_t), target :: e1u_device, e1v_device, e1t_device
    integer(c_intptr_t), target :: e2u_device, e2t_device, e12t_device
    integer(c_intptr_t), target :: e12u_device, tmask_device, gphiu_device
    real(kind=wp), target :: rdt, cbfr, visc
    ! Locals
    integer :: arg_idx, ret

    arg_idx = 0
    ret = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(ua_device), C_LOC(ua_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(ua_device), C_LOC(un_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(vn_device), C_LOC(vn_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(hu_device), C_LOC(hu_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(hv_device), C_LOC(hv_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(ht_device), C_LOC(ht_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(ssha_u_device), &
                         C_LOC(ssha_u_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(sshn_device), &
		         C_LOC(sshn_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(sshn_u_device), &
                         C_LOC(sshn_u_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(sshn_v_device), &
		         C_LOC(sshn_v_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(tmask_device), &
		         C_LOC(tmask_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(e1u_device), &
		         C_LOC(e1u_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(e1v_device), &
		         C_LOC(e1v_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(e1t_device), &
		         C_LOC(e1t_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(e2u_device), &
		         C_LOC(e2u_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1    
    ret = clSetKernelArg(kern, arg_idx, sizeof(e2t_device), &
		         C_LOC(e2t_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1    
    ret = clSetKernelArg(kern, arg_idx, sizeof(e12u_device), &
		         C_LOC(e12u_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1    
    ret = clSetKernelArg(kern, arg_idx, sizeof(gphiu_device), &
		         C_LOC(gphiu_device))
    call check_status("clSetKernelArg", ret);
    arg_idx = arg_idx + 1    
    ret = clSetKernelArg(kern, arg_idx, sizeof(rdt), C_LOC(rdt))
    call check_status("clSetKernelArg", ret);
    arg_idx = arg_idx + 1    
    ret = clSetKernelArg(kern, arg_idx, sizeof(cbfr), C_LOC(cbfr))
    call check_status("clSetKernelArg", ret);
    arg_idx = arg_idx + 1    
    ret = clSetKernelArg(kern, arg_idx, sizeof(visc), C_LOC(visc))
    call check_status("clSetKernelArg", ret);
    write(*,"('Set ',I2,' arguments for Momentum-u kernel')") arg_idx

  end subroutine set_momu_args

end module kernel_args_mod
