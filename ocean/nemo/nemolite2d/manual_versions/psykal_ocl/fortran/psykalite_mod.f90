module psykalite_mod
  use iso_c_binding
  use clfortran
  use opencl_utils_mod, only: check_status
  implicit none

contains

  subroutine invoke_kernels(istp, fld1, fld2, ...)
    character(len=*), intent(in) :: kernel1, kernel2, ...
    integer, parameter :: NUM_KERNELS = 10
    ! Array of kernel objects used in this invoke
    integer(c_intptr_t), save, target :: kernels(NUM_KERNELS)
    integer, save :: num_cmd_queues
    ! Array of command queues - used to achieve concurrent execution
    integer(c_intptr_t), pointer, save :: cmd_queues(:)

    integer(c_int32_t) :: ierr
    integer(c_size_t),target :: globalsize(2), localsize(2)
    integer(c_intptr_t), target :: write_event
    logical, save :: first_time

    if(first_time)then
       context = get_cl_context() ! TODO infrastructure
       device = get_cl_device() ! TODO infrastructure
       num_cmd_queues = get_num_cmd_queues()
       cmd_queues = get_cmd_queues()

       ! Get pointers to the kernel objects used in this invoke
       kernels(1) = get_kernel(name="Name of kernel1 from meta-data")
       kernels(2) = get_kernel(name="Name of kernel2 from meta-data")

       ! Make sure field data is on the device (they might have been put
       ! there by a previous invoke)
       if(.not. ssha_t%data_on_device)then
          ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_t%device_ptr, &
                                      CL_TRUE, 0_8, &
                                      size_in_bytes, C_LOC(ssha_t%data), &
                                      0, C_NULL_PTR, &
                                      C_LOC(write_event))
          ssha_t%data_on_device = .true.
       end if
       ! Make sure any required grid properties are on the device

       ! Wait for the last write to complete
       ierr = clWaitForEvents(1, C_LOC(write_event));
       call check_status("clWaitForEvents", ierr)
    end if

    ! Set-up the kernel arguments. We have to do this every time an invoke
    ! is called in case different invokes call the same kernel(s) but with
    ! different arguments.
    call set_continuity_args(kernels(K_CONTINUITY), &
            model_grid%nx,         &
            ssha_device,           &
            sshn_device,           &
            sshn_u_device,         &
            sshn_v_device,         &
            hu_device,             &
            hv_device,             &
            un_device,             &
            vn_device,             &
            rdt,                   &
            e12t_device)

    call set_momu_args(kernels(K_MOM_U), &
            model_grid%nx,    &
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

    call set_momv_args(kernels(K_MOM_V), &
            model_grid%nx,    &
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
            tmask_device,  &
            e1v_device, e1t_device, e2u_device, &
            e2v_device, e2t_device, e12v_device,&
            gphiv_device, &
            rdt, cbfr, visc)

    call set_bc_ssh_args(kernels(K_BC_SSH), &
            model_grid%nx,     &
            ssha_device,       &
            tmask_device,      &
            rdt)
    call set_bc_solid_u_args(kernels(K_BC_SOLID_U), &
            model_grid%nx,         &
            ua_device, tmask_device)
    call set_bc_solid_v_args(kernels(K_BC_SOLID_V), &
            model_grid%nx,         &
            va_device, tmask_device)
    call set_bc_flather_u_args(kernels(K_BC_FLATHER_U), &
            model_grid%nx,           &
            ua_device, hu_device, sshn_u_device, tmask_device)
    call set_bc_flather_v_args(kernels(K_BC_FLATHER_V), &
            model_grid%nx,           &
            va_device, hv_device, sshn_v_device, tmask_device)
    call set_next_sshu_args(kernels(K_NEXT_SSH_U), &
            model_grid%nx,         &
            sshn_u_device, sshn_device, &
            tmask_device, e12t_device, e12u_device)
    call set_next_sshv_args(kernels(K_NEXT_SSH_V), &
            model_grid%nx,         &
            sshn_v_device, sshn_device, &
            tmask_device, e12t_device, e12v_device)

    globalsize(1) = fld1%grid%nx
    globalsize(2) = fld1%grid%ny


    ! Execute the Continuity kernel
    ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(1), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR)
     call check_status('clEnqueueNDRangeKernel', ierr)

     ! Execute the u-Momentum kernel
     ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernels(2), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR)
     call check_status('clEnqueueNDRangeKernel', ierr)

     ! Execute the v-Momentum kernel
     ierr = clEnqueueNDRangeKernel(cmd_queues(3), kernels(3), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR)
     call check_status('clEnqueueNDRangeKernel', ierr)

     ! Pass the current time-step index to the bc-ssh kernel
     ierr = clSetKernelArg(kernels(4), 1, sizeof(istp), C_LOC(istp))
     call check_status("clSetKernelArg", ierr)

     ! This kernel updates ssha_t and therefore has a dependence on the
     ! continuity kernel
     ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(4), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR)
     call check_status("clEnqueueNDRangeKernel(bc-ssh)", ierr)

     ! Apply boundary conditions
     ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernels(5), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR)
     call check_status("clEnqueueNDRangeKernel(bc-solid-u)", ierr)

     ierr = clEnqueueNDRangeKernel(cmd_queues(3), kernels(6), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR)
     call check_status("clEnqueueNDRangeKernel(bc-solid-v)", ierr)

     ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernels(7), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR);
     call check_status("clEnqueueNDRangeKernel(bc-flather-u)", ierr)

     ierr = clEnqueueNDRangeKernel(cmd_queues(3), kernels(8), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR);
     call check_status("clEnqueueNDRangeKernel(bc-flather-v)", ierr)

     ! Copy 'after' fields to 'now' fields. We could simply swap
     ! pointers around here but that's an optimisation.
     ierr = clEnqueueCopyBuffer(cmd_queues(2), ua%device_ptr, un%device_ptr, &
          0_8, 0_8, size_in_bytes, 0, C_NULL_PTR, C_NULL_PTR)
     call check_status("clEnqueueCopyBuffer", ierr)
     ierr = clEnqueueCopyBuffer(cmd_queues(3), va%device_ptr, vn%device_ptr, &
          0_8, 0_8, size_in_bytes, 0, C_NULL_PTR, C_NULL_PTR)
     call check_status("clEnqueueCopyBuffer", ierr)
     ierr = clEnqueueCopyBuffer(cmd_queues(1), ssha%device_ptr, sshn%device_ptr, &
          0_8, 0_8, size_in_bytes, 0, C_NULL_PTR, C_LOC(write_events(1)))
     call check_status("clEnqueueCopyBuffer", ierr)

     ! Update of sshu and sshv fields
     ! We do sshu in the same queue as used to perform the copy of ssha to
     ! sshn so correct ordering is guaranteed.
     ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(9), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_LOC(ssh_events(1)))
     call check_status("clEnqueueNDRangeKernel(next-sshu)", ierr)

     ! This can be done in parallel with next_sshu but we need to ensure
     ! it only happens after the copy of ssha to sshn.
     ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernels(10), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          1, C_LOC(write_events(1)), C_LOC(ssh_events(2)))
     call check_status("clEnqueueNDRangeKernel(next-sshv)", ierr)

     ! Wait for the ssh updates to complete as we can't execute the
     ! Continuity kernel (at the beginning of the next loop iteration)
     ! until they have.
     ierr = clWaitForEvents(2, C_LOC(ssh_events));
     call check_status("clWaitForEvents", ierr)

   end subroutine invoke_kernels

end module psykalite_mod
