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

  !> Set the arguments for the u-Momentum kernel
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
    integer(c_intptr_t), target :: e2u_device, e2t_device
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
    ret = clSetKernelArg(kern, arg_idx, sizeof(un_device), C_LOC(un_device))
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

  
  !> Set the arguments for the v-Momentum kernel
  subroutine set_momv_args(kern,      &
		           nx,        &
		           va_device, &
		           un_device, &
		           vn_device, &
		           hu_device, &
		           hv_device, &
		           ht_device, &
		           ssha_v_device, &
		           sshn_device,   &
		           sshn_u_device, &
		           sshn_v_device, &
		           tmask_device, &
		           e1v_device, e1t_device, &
		           e2u_device, e2v_device, &
		           e2t_device, e12v_device,&
		           gphiv_device,           &
		           rdt, cbfr, visc)
    integer(c_intptr_t), target :: kern
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: ssha_device, ssha_v_device
    integer(c_intptr_t), target :: sshn_device, sshn_u_device, sshn_v_device
    integer(c_intptr_t), target :: va_device, un_device, vn_device
    integer(c_intptr_t), target :: hu_device, hv_device, ht_device
    integer(c_intptr_t), target :: e1v_device, e1t_device, e12v_device
    integer(c_intptr_t), target :: e2u_device, e2v_device, e2t_device
    integer(c_intptr_t), target :: tmask_device, gphiv_device
    real(kind=wp), target :: rdt, cbfr, visc
    ! Locals
    integer :: arg_idx, ret

    arg_idx = 0
    ret = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(va_device), C_LOC(va_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1
    ret = clSetKernelArg(kern, arg_idx, sizeof(un_device), C_LOC(un_device))
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
    ret = clSetKernelArg(kern, arg_idx, sizeof(ssha_v_device), &
                         C_LOC(ssha_v_device))
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
    ret = clSetKernelArg(kern, arg_idx, sizeof(e2v_device), &
		         C_LOC(e2v_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1    
    ret = clSetKernelArg(kern, arg_idx, sizeof(e2t_device), &
		         C_LOC(e2t_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1    
    ret = clSetKernelArg(kern, arg_idx, sizeof(e12v_device), &
		         C_LOC(e12v_device))
    call check_status("clSetKernelArg", ret)
    arg_idx = arg_idx + 1    
    ret = clSetKernelArg(kern, arg_idx, sizeof(gphiv_device), &
		         C_LOC(gphiv_device))
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
    arg_idx = arg_idx + 1

    write(*,"('Set ',I2,' arguments for Momentum-v kernel')") arg_idx

  end subroutine set_momv_args

  subroutine set_bc_ssh_args(kern, nx, ssha, tmask, rdt)
    integer(c_intptr_t), target :: kern
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: ssha, tmask
    real(kind=wp), target :: rdt
    ! Locals
    integer :: arg_idx, ierr

    ! set the kernel arguments
    arg_idx = 0
    ierr = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ! This argument is actually the current time step so we will
    ! set it (again) during the time-stepping loop.
    ierr = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    
    ierr = clSetKernelArg(kern, arg_idx, sizeof(ssha), C_LOC(ssha))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(tmask), C_LOC(tmask))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(rdt), C_LOC(rdt))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1

    write(*,"('Set ',I2,' arguments for bc-ssh kernel')") arg_idx

  end subroutine set_bc_ssh_args

  subroutine set_bc_solid_u_args(kern, nx, ua, tmask)
    integer(c_intptr_t), target :: kern
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: ua, tmask
    ! Locals
    integer :: arg_idx, ierr

    ! set the kernel arguments
    arg_idx = 0
    ierr = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(ua), C_LOC(ua))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(tmask), C_LOC(tmask))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1

    write(*,"('Set ',I2,' arguments for bc-solid-u kernel')") arg_idx
    
  end subroutine set_bc_solid_u_args

  subroutine set_bc_solid_v_args(kern, nx, va, tmask)
    integer(c_intptr_t), target :: kern
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: va, tmask
    ! Locals
    integer :: arg_idx, ierr

    ! set the kernel arguments
    arg_idx = 0
    ierr = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(va), C_LOC(va))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(tmask), C_LOC(tmask))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1

    write(*,"('Set ',I2,' arguments for bc-solid-v kernel')") arg_idx
  end subroutine set_bc_solid_v_args

  subroutine set_bc_flather_u_args(kern, nx, ua, hu, sshn_u, tmask)
    integer(c_intptr_t), target :: kern
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: ua, hu, sshn_u, tmask
    ! Locals
    integer :: arg_idx, ierr

    ! set the kernel arguments
    arg_idx = 0
    ierr = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(ua), C_LOC(ua))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(hu), C_LOC(hu))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(sshn_u), C_LOC(sshn_u))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(tmask), C_LOC(tmask))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1

    write(*,"('Set ',I2,' arguments for bc-flather-u kernel')") arg_idx
  end subroutine set_bc_flather_u_args

  subroutine set_bc_flather_v_args(kern, nx, va, hv, sshn_v, tmask)
    integer(c_intptr_t), target :: kern
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: va, hv, sshn_v, tmask
    ! Locals
    integer :: arg_idx, ierr

    ! set the kernel arguments
    arg_idx = 0
    ierr = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(va), C_LOC(va))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(hv), C_LOC(hv))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(sshn_v), C_LOC(sshn_v))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(tmask), C_LOC(tmask))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1

    write(*,"('Set ',I2,' arguments for bc-flather-v kernel')") arg_idx
    
  end subroutine set_bc_flather_v_args

  subroutine set_next_sshu_args(kern, nx, sshn_u, sshn, tmask, e12t, e12u)
    integer(c_intptr_t), target :: kern
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: sshn, sshn_u, tmask, e12t, e12u
    ! Locals
    integer :: arg_idx, ierr

    ! set the kernel arguments
    arg_idx = 0
    ierr = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(sshn_u), C_LOC(sshn_u))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(sshn), C_LOC(sshn))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(tmask), C_LOC(tmask))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(e12t), C_LOC(e12t))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(e12u), C_LOC(e12u))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1

    write(*,"('Set ',I2,' arguments for next-sshu kernel')") arg_idx
  end subroutine set_next_sshu_args

  subroutine set_next_sshv_args(kern, nx, sshn_v, sshn, tmask, e12t, e12v)
    integer(c_intptr_t), target :: kern
    integer, target, intent(in) :: nx
    integer(c_intptr_t), target :: sshn, sshn_v, tmask, e12t, e12v
    ! Locals
    integer :: arg_idx, ierr

    ! set the kernel arguments
    arg_idx = 0
    ierr = clSetKernelArg(kern, arg_idx, sizeof(nx), C_LOC(nx))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1    
    ierr = clSetKernelArg(kern, arg_idx, sizeof(sshn_v), C_LOC(sshn_v))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(sshn), C_LOC(sshn))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(tmask), C_LOC(tmask))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(e12t), C_LOC(e12t))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1
    ierr = clSetKernelArg(kern, arg_idx, sizeof(e12v), C_LOC(e12v))
    call check_status("clSetKernelArg", ierr)
    arg_idx = arg_idx + 1

    write(*,"('Set ',I2,' arguments for next-sshv kernel')") arg_idx

  end subroutine set_next_sshv_args

end module kernel_args_mod
