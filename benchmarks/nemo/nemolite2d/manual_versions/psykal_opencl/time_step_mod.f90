module time_step_mod
    use field_mod
    use kind_params_mod
    implicit none
    private
    public invoke_time_step
contains

    SUBROUTINE invoke_time_step(ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, &
                                vn, ua, ht, ssha_u, va, ssha_v, istp)
        USE field_mod
        USE kind_params_mod
        USE fortcl, ONLY: create_rw_buffer
        USE fortcl, ONLY: get_num_cmd_queues, get_cmd_queues, &
                          get_kernel_by_name
        USE clfortran
        USE iso_c_binding
        USE physical_params_mod, ONLY: omega, d2r, g
        USE model_mod, ONLY: rdt, cbfr, visc
        TYPE(r2d_field), intent(inout), target :: ssha_t, sshn_t, sshn_u, &
            sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
        INTEGER, intent(inout) :: istp
        INTEGER(KIND=c_intptr_t), target :: write_event
        INTEGER(KIND=c_size_t) size_in_bytes
        INTEGER(KIND=c_size_t), target :: localsize(2)
        INTEGER(KIND=c_size_t), target :: globalsize(2)
        INTEGER(KIND=c_intptr_t), target, save :: kernel_next_sshv_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_next_sshu_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_field_copy_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_flather_v_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_flather_u_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_solid_v_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_solid_u_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_ssh_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_momentum_v_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_momentum_u_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_continuity_code
        LOGICAL, save :: first_time=.true.
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), pointer, save :: cmd_queues(:)
        INTEGER, save :: num_cmd_queues
        INTEGER istop, jstop

        ! Look-up loop bounds
        istop = ssha_t%grid%subdomain%internal%xstop
        jstop = ssha_t%grid%subdomain%internal%ystop

        IF (first_time) THEN
            first_time = .false.
            ! Ensure OpenCL run-time is initialised for this PSy-layer module
            CALL psy_init
            num_cmd_queues = get_num_cmd_queues()
            cmd_queues => get_cmd_queues()
            kernel_continuity_code = get_kernel_by_name("continuity_code")
            kernel_momentum_u_code = get_kernel_by_name("momentum_u_code")
            kernel_momentum_v_code = get_kernel_by_name("momentum_v_code")
            kernel_bc_ssh_code = get_kernel_by_name("bc_ssh_code")
            kernel_bc_solid_u_code = get_kernel_by_name("bc_solid_u_code")
            kernel_bc_solid_v_code = get_kernel_by_name("bc_solid_v_code")
            kernel_bc_flather_u_code = get_kernel_by_name("bc_flather_u_code")
            kernel_bc_flather_v_code = get_kernel_by_name("bc_flather_v_code")
            kernel_field_copy_code = get_kernel_by_name("field_copy_code")
            kernel_field_copy_code = get_kernel_by_name("field_copy_code")
            kernel_field_copy_code = get_kernel_by_name("field_copy_code")
            kernel_next_sshu_code = get_kernel_by_name("next_sshu_code")
            kernel_next_sshv_code = get_kernel_by_name("next_sshv_code")
        END IF

        ! Set up local and global sizes (for memory blocking)
        globalsize = (/sshn_t%grid%nx, sshn_t%grid%ny/)
        localsize = (/1, 1/)

        ! Ensure field data is on device
        IF (.NOT. ssha_t%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8) * &
                            c_sizeof(ssha_t%data(1,1))
            ! Create buffer on device
            ssha_t%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_t%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_t%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ssha_t%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_t%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8) * &
                            c_sizeof(sshn_t%data(1,1))
            ! Create buffer on device
            sshn_t%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_t%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_u%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8) * &
                c_sizeof(sshn_u%data(1,1))
            ! Create buffer on device
            sshn_u%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_u%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_u%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_v%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8) * &
                            c_sizeof(sshn_v%data(1,1))
            ! Create buffer on device
            sshn_v%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_v%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_v%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. hu%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8) * &
                            c_sizeof(hu%data(1,1))
            ! Create buffer on device
            hu%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), hu%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(hu%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            hu%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. hv%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8) * &
                            c_sizeof(hv%data(1,1))
            ! Create buffer on device
            hv%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), hv%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(hv%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            hv%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. un%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8) * &
                            c_sizeof(un%data(1,1))
            ! Create buffer on device
            un%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(un%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            un%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. vn%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8) * &
                            c_sizeof(vn%data(1,1))
            ! Create buffer on device
            vn%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), vn%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(vn%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            vn%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (sshn_t%grid%area_t_device == 0) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8) * &
                            c_sizeof(sshn_t%grid%area_t(1,1))
            ! Create buffer on device
            sshn_t%grid%area_t_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), &
                sshn_t%grid%area_t_device, CL_TRUE, 0_8, size_in_bytes, &
                &C_LOC(sshn_t%grid%area_t), 0, C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. ua%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8) * &
                            c_sizeof(ua%data(1,1))
            ! Create buffer on device
            ua%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            ua%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. ht%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8) * &
                            c_sizeof(ht%data(1,1))
            ! Create buffer on device
            ht%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ht%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(ht%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            ht%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. ssha_u%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8) * &
                            c_sizeof(ssha_u%data(1,1))
            ! Create buffer on device
            ssha_u%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_u%device_ptr, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_u%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ssha_u%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF

        ! Ensure scalar data is on the device
        IF (un%grid%tmask_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8) * &
                            c_sizeof(un%grid%tmask(1,1))
            ! Create buffer on device
            un%grid%tmask_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%tmask_device, &
                CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%tmask), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF

        CALL continuity_code_set_args(&
            kernel_continuity_code, ssha_t%device_ptr, sshn_t%device_ptr, &
            sshn_u%device_ptr, sshn_v%device_ptr, hu%device_ptr, &
            hv%device_ptr, un%device_ptr, vn%device_ptr, &
            sshn_t%grid%area_t_device, rdt)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_continuity_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        
        ! Ensure field data is on device
        IF (un%grid%dx_u_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_u(1,1))
            ! Create buffer on device
            un%grid%dx_u_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_u_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_u), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%dx_v_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_v(1,1))
            ! Create buffer on device
            un%grid%dx_v_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_v_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_v), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%dx_t_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_t(1,1))
            ! Create buffer on device
            un%grid%dx_t_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_t_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_t), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%dy_u_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_u(1,1))
            ! Create buffer on device
            un%grid%dy_u_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_u_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_u), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%dy_t_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_t(1,1))
            ! Create buffer on device
            un%grid%dy_t_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_t_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_t), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%area_u_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%area_u(1,1))
            ! Create buffer on device
            un%grid%area_u_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%area_u_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%area_u), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%gphiu_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%gphiu(1,1))
            ! Create buffer on device
            un%grid%gphiu_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%gphiu_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%gphiu), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL momentum_u_code_set_args(kernel_momentum_u_code, ua%device_ptr, un%device_ptr, vn%device_ptr, hu%device_ptr, &
            &hv%device_ptr, ht%device_ptr, ssha_u%device_ptr, sshn_t%device_ptr, sshn_u%device_ptr, sshn_v%device_ptr, un%grid%tmask_device, &
            &un%grid%dx_u_device, un%grid%dx_v_device, un%grid%dx_t_device, un%grid%dy_u_device, un%grid%dy_t_device, un%grid%area_u_device, &
            &un%grid%gphiu_device, omega, d2r, g, rdt, cbfr, visc)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_momentum_u_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. va%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(va%data(1,1))
            ! Create buffer on device
            va%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), va%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            va%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. un%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%data(1,1))
            ! Create buffer on device
            un%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            un%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. vn%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(vn%data(1,1))
            ! Create buffer on device
            vn%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), vn%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(vn%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            vn%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. hu%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(hu%data(1,1))
            ! Create buffer on device
            hu%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), hu%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hu%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            hu%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. hv%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(hv%data(1,1))
            ! Create buffer on device
            hv%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), hv%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hv%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            hv%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. ht%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(ht%data(1,1))
            ! Create buffer on device
            ht%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ht%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ht%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            ht%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. ssha_v%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(ssha_v%data(1,1))
            ! Create buffer on device
            ssha_v%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_v%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ssha_v%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_t%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
            ! Create buffer on device
            sshn_t%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_t%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_u%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(sshn_u%data(1,1))
            ! Create buffer on device
            sshn_u%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_u%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_u%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_v%data_on_device) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(sshn_v%data(1,1))
            ! Create buffer on device
            sshn_v%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_v%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_v%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%tmask_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%tmask(1,1))
            ! Create buffer on device
            un%grid%tmask_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%tmask), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%dx_v_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_v(1,1))
            ! Create buffer on device
            un%grid%dx_v_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_v_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_v), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%dx_t_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_t(1,1))
            ! Create buffer on device
            un%grid%dx_t_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_t_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_t), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%dy_u_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_u(1,1))
            ! Create buffer on device
            un%grid%dy_u_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_u_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_u), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%dy_v_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_v(1,1))
            ! Create buffer on device
            un%grid%dy_v_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_v_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_v), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%dy_t_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_t(1,1))
            ! Create buffer on device
            un%grid%dy_t_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_t_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_t), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%area_v_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%area_v(1,1))
            ! Create buffer on device
            un%grid%area_v_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%area_v_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%area_v), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (un%grid%gphiv_device == 0) THEN
            size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%gphiv(1,1))
            ! Create buffer on device
            un%grid%gphiv_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%gphiv_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%gphiv), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL momentum_v_code_set_args(kernel_momentum_v_code, va%device_ptr, un%device_ptr, vn%device_ptr, hu%device_ptr, &
            &hv%device_ptr, ht%device_ptr, ssha_v%device_ptr, sshn_t%device_ptr, sshn_u%device_ptr, sshn_v%device_ptr, un%grid%tmask_device, &
            &un%grid%dx_v_device, un%grid%dx_t_device, un%grid%dy_u_device, un%grid%dy_v_device, un%grid%dy_t_device, un%grid%area_v_device, &
            &un%grid%gphiv_device, omega, d2r, g, rdt, cbfr, visc)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_momentum_v_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. ssha_t%data_on_device) THEN
            size_in_bytes = int(ssha_t%grid%nx*ssha_t%grid%ny, 8)*c_sizeof(ssha_t%data(1,1))
            ! Create buffer on device
            ssha_t%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_t%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ssha_t%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (ssha_t%grid%tmask_device == 0) THEN
            size_in_bytes = int(ssha_t%grid%nx*ssha_t%grid%ny, 8)*c_sizeof(ssha_t%grid%tmask(1,1))
            ! Create buffer on device
            ssha_t%grid%tmask_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_t%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, &
                &C_LOC(ssha_t%grid%tmask), 0, C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL bc_ssh_code_set_args(kernel_bc_ssh_code, istp, ssha_t%device_ptr, ssha_t%grid%tmask_device, rdt)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_ssh_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. ua%data_on_device) THEN
            size_in_bytes = int(ua%grid%nx*ua%grid%ny, 8)*c_sizeof(ua%data(1,1))
            ! Create buffer on device
            ua%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            ua%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (ua%grid%tmask_device == 0) THEN
            size_in_bytes = int(ua%grid%nx*ua%grid%ny, 8)*c_sizeof(ua%grid%tmask(1,1))
            ! Create buffer on device
            ua%grid%tmask_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%grid%tmask), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL bc_solid_u_code_set_args(kernel_bc_solid_u_code, ua%device_ptr, ua%grid%tmask_device)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_solid_u_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. va%data_on_device) THEN
            size_in_bytes = int(va%grid%nx*va%grid%ny, 8)*c_sizeof(va%data(1,1))
            ! Create buffer on device
            va%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), va%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            va%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (va%grid%tmask_device == 0) THEN
            size_in_bytes = int(va%grid%nx*va%grid%ny, 8)*c_sizeof(va%grid%tmask(1,1))
            ! Create buffer on device
            va%grid%tmask_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), va%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%grid%tmask), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL bc_solid_v_code_set_args(kernel_bc_solid_v_code, va%device_ptr, va%grid%tmask_device)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_solid_v_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. ua%data_on_device) THEN
            size_in_bytes = int(hu%grid%nx*hu%grid%ny, 8)*c_sizeof(ua%data(1,1))
            ! Create buffer on device
            ua%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            ua%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. hu%data_on_device) THEN
            size_in_bytes = int(hu%grid%nx*hu%grid%ny, 8)*c_sizeof(hu%data(1,1))
            ! Create buffer on device
            hu%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), hu%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hu%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            hu%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_u%data_on_device) THEN
            size_in_bytes = int(hu%grid%nx*hu%grid%ny, 8)*c_sizeof(sshn_u%data(1,1))
            ! Create buffer on device
            sshn_u%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_u%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_u%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (hu%grid%tmask_device == 0) THEN
            size_in_bytes = int(hu%grid%nx*hu%grid%ny, 8)*c_sizeof(hu%grid%tmask(1,1))
            ! Create buffer on device
            hu%grid%tmask_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), hu%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(hu%grid%tmask), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL bc_flather_u_code_set_args(kernel_bc_flather_u_code, ua%device_ptr, hu%device_ptr, sshn_u%device_ptr, &
            &hu%grid%tmask_device, g)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_flather_u_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), &
            &0, C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. va%data_on_device) THEN
            size_in_bytes = int(hv%grid%nx*hv%grid%ny, 8)*c_sizeof(va%data(1,1))
            ! Create buffer on device
            va%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), va%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            va%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. hv%data_on_device) THEN
            size_in_bytes = int(hv%grid%nx*hv%grid%ny, 8)*c_sizeof(hv%data(1,1))
            ! Create buffer on device
            hv%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), hv%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hv%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            hv%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_v%data_on_device) THEN
            size_in_bytes = int(hv%grid%nx*hv%grid%ny, 8)*c_sizeof(sshn_v%data(1,1))
            ! Create buffer on device
            sshn_v%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_v%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_v%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (hv%grid%tmask_device == 0) THEN
            size_in_bytes = int(hv%grid%nx*hv%grid%ny, 8)*c_sizeof(hv%grid%tmask(1,1))
            ! Create buffer on device
            hv%grid%tmask_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), hv%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(hv%grid%tmask), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL bc_flather_v_code_set_args(kernel_bc_flather_v_code, va%device_ptr, hv%device_ptr, sshn_v%device_ptr, &
            &hv%grid%tmask_device, g)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_flather_v_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), &
            &0, C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. un%data_on_device) THEN
            size_in_bytes = int(ua%grid%nx*ua%grid%ny, 8)*c_sizeof(un%data(1,1))
            ! Create buffer on device
            un%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), un%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            un%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. ua%data_on_device) THEN
            size_in_bytes = int(ua%grid%nx*ua%grid%ny, 8)*c_sizeof(ua%data(1,1))
            ! Create buffer on device
            ua%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            ua%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL field_copy_code_set_args(kernel_field_copy_code, un%device_ptr, ua%device_ptr)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_field_copy_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. vn%data_on_device) THEN
            size_in_bytes = int(va%grid%nx*va%grid%ny, 8)*c_sizeof(vn%data(1,1))
            ! Create buffer on device
            vn%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), vn%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(vn%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            vn%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. va%data_on_device) THEN
            size_in_bytes = int(va%grid%nx*va%grid%ny, 8)*c_sizeof(va%data(1,1))
            ! Create buffer on device
            va%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), va%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%data), 0, C_NULL_PTR, &
                &C_LOC(write_event))
            va%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL field_copy_code_set_args(kernel_field_copy_code, vn%device_ptr, va%device_ptr)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_field_copy_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. sshn_t%data_on_device) THEN
            size_in_bytes = int(ssha_t%grid%nx*ssha_t%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
            ! Create buffer on device
            sshn_t%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_t%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. ssha_t%data_on_device) THEN
            size_in_bytes = int(ssha_t%grid%nx*ssha_t%grid%ny, 8)*c_sizeof(ssha_t%data(1,1))
            ! Create buffer on device
            ssha_t%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_t%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            ssha_t%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL field_copy_code_set_args(kernel_field_copy_code, sshn_t%device_ptr, ssha_t%device_ptr)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_field_copy_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. sshn_u%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_u%data(1,1))
            ! Create buffer on device
            sshn_u%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_u%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_u%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_t%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
            ! Create buffer on device
            sshn_t%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_t%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (sshn_t%grid%tmask_device == 0) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%tmask(1,1))
            ! Create buffer on device
            sshn_t%grid%tmask_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, &
                &C_LOC(sshn_t%grid%tmask), 0, C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (sshn_t%grid%area_t_device == 0) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%area_t(1,1))
            ! Create buffer on device
            sshn_t%grid%area_t_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%area_t_device, CL_TRUE, 0_8, size_in_bytes, &
                &C_LOC(sshn_t%grid%area_t), 0, C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (sshn_t%grid%area_u_device == 0) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%area_u(1,1))
            ! Create buffer on device
            sshn_t%grid%area_u_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%area_u_device, CL_TRUE, 0_8, size_in_bytes, &
                &C_LOC(sshn_t%grid%area_u), 0, C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL next_sshu_code_set_args(kernel_next_sshu_code, sshn_u%device_ptr, sshn_t%device_ptr, sshn_t%grid%tmask_device, &
            &sshn_t%grid%area_t_device, sshn_t%grid%area_u_device)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_next_sshu_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Ensure field data is on device
        IF (.NOT. sshn_v%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_v%data(1,1))
            ! Create buffer on device
            sshn_v%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_v%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_v%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (.NOT. sshn_t%data_on_device) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
            ! Create buffer on device
            sshn_t%device_ptr = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
                &C_NULL_PTR, C_LOC(write_event))
            sshn_t%data_on_device = .true.
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (sshn_t%grid%tmask_device == 0) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%tmask(1,1))
            ! Create buffer on device
            sshn_t%grid%tmask_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, &
                &C_LOC(sshn_t%grid%tmask), 0, C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (sshn_t%grid%area_t_device == 0) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%area_t(1,1))
            ! Create buffer on device
            sshn_t%grid%area_t_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%area_t_device, CL_TRUE, 0_8, size_in_bytes, &
                &C_LOC(sshn_t%grid%area_t), 0, C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        IF (sshn_t%grid%area_v_device == 0) THEN
            size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%area_v(1,1))
            ! Create buffer on device
            sshn_t%grid%area_v_device = create_rw_buffer(size_in_bytes)
            ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%area_v_device, CL_TRUE, 0_8, size_in_bytes, &
                &C_LOC(sshn_t%grid%area_v), 0, C_NULL_PTR, C_LOC(write_event))
            ! Block until data copies have finished
            ierr = clFinish(cmd_queues(1))
        END IF
        CALL next_sshv_code_set_args(kernel_next_sshv_code, sshn_v%device_ptr, sshn_t%device_ptr, sshn_t%grid%tmask_device, &
            &sshn_t%grid%area_t_device, sshn_t%grid%area_v_device)
        ! Launch the kernel
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_next_sshv_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            &C_NULL_PTR, C_NULL_PTR)
        !
        ! Block until all kernels have finished
        ierr = clFinish(cmd_queues(1))
    END SUBROUTINE invoke_time_step

    SUBROUTINE continuity_code_set_args(kernel_obj, ssha_t, sshn_t, &
            sshn_u, sshn_v, hu, hv, un, vn, area_t, rdt)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: ssha_t, sshn_t, &
            sshn_u, sshn_v, hu, hv, un, vn, area_t
        REAL(KIND=go_wp), intent(in), target :: rdt
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the continuity_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(ssha_t), C_LOC(ssha_t))
        CALL check_status('clSetKernelArg: arg 0 of continuity_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(sshn_t), C_LOC(sshn_t))
        CALL check_status('clSetKernelArg: arg 1 of continuity_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(sshn_u), C_LOC(sshn_u))
        CALL check_status('clSetKernelArg: arg 2 of continuity_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(sshn_v), C_LOC(sshn_v))
        CALL check_status('clSetKernelArg: arg 3 of continuity_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(hu), C_LOC(hu))
        CALL check_status('clSetKernelArg: arg 4 of continuity_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(hv), C_LOC(hv))
        CALL check_status('clSetKernelArg: arg 5 of continuity_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(un), C_LOC(un))
        CALL check_status('clSetKernelArg: arg 6 of continuity_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(vn), C_LOC(vn))
        CALL check_status('clSetKernelArg: arg 7 of continuity_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(area_t), C_LOC(area_t))
        CALL check_status('clSetKernelArg: arg 8 of continuity_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 9, C_SIZEOF(rdt), C_LOC(rdt))
        CALL check_status('clSetKernelArg: arg 9 of continuity_code', ierr)
    END SUBROUTINE continuity_code_set_args

    SUBROUTINE momentum_u_code_set_args(kernel_obj, ua, un, vn, hu, hv, ht, &
            ssha_u, sshn_t, sshn_u, sshn_v, tmask, dx_u, dx_v, &
            dx_t, dy_u, dy_t, area_u, gphiu, omega, d2r, g, rdt, cbfr, visc)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: ua, un, vn, hu, hv, &
            ht, ssha_u, sshn_t, sshn_u, sshn_v, tmask, dx_u, dx_v, &
            dx_t, dy_u, dy_t, area_u, gphiu
        REAL(KIND=go_wp), intent(in), target :: omega
        REAL(KIND=go_wp), intent(in), target :: d2r
        REAL(KIND=go_wp), intent(in), target :: g
        REAL(KIND=go_wp), intent(in), target :: rdt
        REAL(KIND=go_wp), intent(in), target :: cbfr
        REAL(KIND=go_wp), intent(in), target :: visc
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the momentum_u_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(ua), C_LOC(ua))
        CALL check_status('clSetKernelArg: arg 0 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(un), C_LOC(un))
        CALL check_status('clSetKernelArg: arg 1 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(vn), C_LOC(vn))
        CALL check_status('clSetKernelArg: arg 2 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(hu), C_LOC(hu))
        CALL check_status('clSetKernelArg: arg 3 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(hv), C_LOC(hv))
        CALL check_status('clSetKernelArg: arg 4 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(ht), C_LOC(ht))
        CALL check_status('clSetKernelArg: arg 5 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(ssha_u), C_LOC(ssha_u))
        CALL check_status('clSetKernelArg: arg 6 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(sshn_t), C_LOC(sshn_t))
        CALL check_status('clSetKernelArg: arg 7 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(sshn_u), C_LOC(sshn_u))
        CALL check_status('clSetKernelArg: arg 8 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 9, C_SIZEOF(sshn_v), C_LOC(sshn_v))
        CALL check_status('clSetKernelArg: arg 9 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 10, C_SIZEOF(tmask), C_LOC(tmask))
        CALL check_status('clSetKernelArg: arg 10 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 11, C_SIZEOF(dx_u), C_LOC(dx_u))
        CALL check_status('clSetKernelArg: arg 11 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 12, C_SIZEOF(dx_v), C_LOC(dx_v))
        CALL check_status('clSetKernelArg: arg 12 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 13, C_SIZEOF(dx_t), C_LOC(dx_t))
        CALL check_status('clSetKernelArg: arg 13 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 14, C_SIZEOF(dy_u), C_LOC(dy_u))
        CALL check_status('clSetKernelArg: arg 14 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 15, C_SIZEOF(dy_t), C_LOC(dy_t))
        CALL check_status('clSetKernelArg: arg 15 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 16, C_SIZEOF(area_u), C_LOC(area_u))
        CALL check_status('clSetKernelArg: arg 16 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 17, C_SIZEOF(gphiu), C_LOC(gphiu))
        CALL check_status('clSetKernelArg: arg 17 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 18, C_SIZEOF(omega), C_LOC(omega))
        CALL check_status('clSetKernelArg: arg 18 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 19, C_SIZEOF(d2r), C_LOC(d2r))
        CALL check_status('clSetKernelArg: arg 19 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 20, C_SIZEOF(g), C_LOC(g))
        CALL check_status('clSetKernelArg: arg 20 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 21, C_SIZEOF(rdt), C_LOC(rdt))
        CALL check_status('clSetKernelArg: arg 21 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 22, C_SIZEOF(cbfr), C_LOC(cbfr))
        CALL check_status('clSetKernelArg: arg 22 of momentum_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 23, C_SIZEOF(visc), C_LOC(visc))
        CALL check_status('clSetKernelArg: arg 23 of momentum_u_code', ierr)
    END SUBROUTINE momentum_u_code_set_args

    SUBROUTINE momentum_v_code_set_args(kernel_obj, va, un, vn, hu, hv, ht, &
            ssha_v, sshn_t, sshn_u, sshn_v, tmask, dx_v, dx_t, &
            dy_u, dy_v, dy_t, area_v, gphiv, omega, d2r, g, rdt, cbfr, visc)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: va, un, vn, hu, hv, &
            ht, ssha_v, sshn_t, sshn_u, sshn_v, tmask, dx_v, dx_t, &
            dy_u, dy_v, dy_t, area_v, gphiv
        REAL(KIND=go_wp), intent(in), target :: omega
        REAL(KIND=go_wp), intent(in), target :: d2r
        REAL(KIND=go_wp), intent(in), target :: g
        REAL(KIND=go_wp), intent(in), target :: rdt
        REAL(KIND=go_wp), intent(in), target :: cbfr
        REAL(KIND=go_wp), intent(in), target :: visc
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the momentum_v_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(va), C_LOC(va))
        CALL check_status('clSetKernelArg: arg 0 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(un), C_LOC(un))
        CALL check_status('clSetKernelArg: arg 1 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(vn), C_LOC(vn))
        CALL check_status('clSetKernelArg: arg 2 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(hu), C_LOC(hu))
        CALL check_status('clSetKernelArg: arg 3 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(hv), C_LOC(hv))
        CALL check_status('clSetKernelArg: arg 4 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(ht), C_LOC(ht))
        CALL check_status('clSetKernelArg: arg 5 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(ssha_v), C_LOC(ssha_v))
        CALL check_status('clSetKernelArg: arg 6 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(sshn_t), C_LOC(sshn_t))
        CALL check_status('clSetKernelArg: arg 7 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(sshn_u), C_LOC(sshn_u))
        CALL check_status('clSetKernelArg: arg 8 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 9, C_SIZEOF(sshn_v), C_LOC(sshn_v))
        CALL check_status('clSetKernelArg: arg 9 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 10, C_SIZEOF(tmask), C_LOC(tmask))
        CALL check_status('clSetKernelArg: arg 10 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 11, C_SIZEOF(dx_v), C_LOC(dx_v))
        CALL check_status('clSetKernelArg: arg 11 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 12, C_SIZEOF(dx_t), C_LOC(dx_t))
        CALL check_status('clSetKernelArg: arg 12 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 13, C_SIZEOF(dy_u), C_LOC(dy_u))
        CALL check_status('clSetKernelArg: arg 13 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 14, C_SIZEOF(dy_v), C_LOC(dy_v))
        CALL check_status('clSetKernelArg: arg 14 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 15, C_SIZEOF(dy_t), C_LOC(dy_t))
        CALL check_status('clSetKernelArg: arg 15 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 16, C_SIZEOF(area_v), C_LOC(area_v))
        CALL check_status('clSetKernelArg: arg 16 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 17, C_SIZEOF(gphiv), C_LOC(gphiv))
        CALL check_status('clSetKernelArg: arg 17 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 18, C_SIZEOF(omega), C_LOC(omega))
        CALL check_status('clSetKernelArg: arg 18 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 19, C_SIZEOF(d2r), C_LOC(d2r))
        CALL check_status('clSetKernelArg: arg 19 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 20, C_SIZEOF(g), C_LOC(g))
        CALL check_status('clSetKernelArg: arg 20 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 21, C_SIZEOF(rdt), C_LOC(rdt))
        CALL check_status('clSetKernelArg: arg 21 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 22, C_SIZEOF(cbfr), C_LOC(cbfr))
        CALL check_status('clSetKernelArg: arg 22 of momentum_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 23, C_SIZEOF(visc), C_LOC(visc))
        CALL check_status('clSetKernelArg: arg 23 of momentum_v_code', ierr)
    END SUBROUTINE momentum_v_code_set_args

    SUBROUTINE bc_ssh_code_set_args(kernel_obj, istp, ssha_t, tmask, rdt)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: ssha_t, tmask
        INTEGER, intent(in), target :: istp
        REAL(KIND=go_wp), intent(in), target :: rdt
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the bc_ssh_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(istp), C_LOC(istp))
        CALL check_status('clSetKernelArg: arg 0 of bc_ssh_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(ssha_t), C_LOC(ssha_t))
        CALL check_status('clSetKernelArg: arg 1 of bc_ssh_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(tmask), C_LOC(tmask))
        CALL check_status('clSetKernelArg: arg 2 of bc_ssh_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(rdt), C_LOC(rdt))
        CALL check_status('clSetKernelArg: arg 3 of bc_ssh_code', ierr)
    END SUBROUTINE bc_ssh_code_set_args

    SUBROUTINE bc_solid_u_code_set_args(kernel_obj, ua, tmask)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: ua, tmask
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the bc_solid_u_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(ua), C_LOC(ua))
        CALL check_status('clSetKernelArg: arg 0 of bc_solid_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(tmask), C_LOC(tmask))
        CALL check_status('clSetKernelArg: arg 1 of bc_solid_u_code', ierr)
    END SUBROUTINE bc_solid_u_code_set_args

    SUBROUTINE bc_solid_v_code_set_args(kernel_obj, va, tmask)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: va, tmask
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the bc_solid_v_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(va), C_LOC(va))
        CALL check_status('clSetKernelArg: arg 0 of bc_solid_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(tmask), C_LOC(tmask))
        CALL check_status('clSetKernelArg: arg 1 of bc_solid_v_code', ierr)
    END SUBROUTINE bc_solid_v_code_set_args

    SUBROUTINE bc_flather_u_code_set_args(kernel_obj, ua, hu, sshn_u, tmask, g)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: ua, hu, sshn_u, tmask
        REAL(KIND=go_wp), intent(in), target :: g
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the bc_flather_u_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(ua), C_LOC(ua))
        CALL check_status('clSetKernelArg: arg 0 of bc_flather_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(hu), C_LOC(hu))
        CALL check_status('clSetKernelArg: arg 1 of bc_flather_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(sshn_u), C_LOC(sshn_u))
        CALL check_status('clSetKernelArg: arg 2 of bc_flather_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(tmask), C_LOC(tmask))
        CALL check_status('clSetKernelArg: arg 3 of bc_flather_u_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(g), C_LOC(g))
        CALL check_status('clSetKernelArg: arg 4 of bc_flather_u_code', ierr)
    END SUBROUTINE bc_flather_u_code_set_args

    SUBROUTINE bc_flather_v_code_set_args(kernel_obj, va, hv, sshn_v, tmask, g)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: va, hv, sshn_v, tmask
        REAL(KIND=go_wp), intent(in), target :: g
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the bc_flather_v_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(va), C_LOC(va))
        CALL check_status('clSetKernelArg: arg 0 of bc_flather_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(hv), C_LOC(hv))
        CALL check_status('clSetKernelArg: arg 1 of bc_flather_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(sshn_v), C_LOC(sshn_v))
        CALL check_status('clSetKernelArg: arg 2 of bc_flather_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(tmask), C_LOC(tmask))
        CALL check_status('clSetKernelArg: arg 3 of bc_flather_v_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(g), C_LOC(g))
        CALL check_status('clSetKernelArg: arg 4 of bc_flather_v_code', ierr)
    END SUBROUTINE bc_flather_v_code_set_args

    SUBROUTINE field_copy_code_set_args(kernel_obj, un, ua)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: un, ua
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the field_copy_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(un), C_LOC(un))
        CALL check_status('clSetKernelArg: arg 0 of field_copy_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(ua), C_LOC(ua))
        CALL check_status('clSetKernelArg: arg 1 of field_copy_code', ierr)
    END SUBROUTINE field_copy_code_set_args

    SUBROUTINE next_sshu_code_set_args(kernel_obj, sshn_u, sshn_t, tmask, area_t, area_u)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: sshn_u, sshn_t, tmask, area_t, area_u
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the next_sshu_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(sshn_u), C_LOC(sshn_u))
        CALL check_status('clSetKernelArg: arg 0 of next_sshu_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(sshn_t), C_LOC(sshn_t))
        CALL check_status('clSetKernelArg: arg 1 of next_sshu_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(tmask), C_LOC(tmask))
        CALL check_status('clSetKernelArg: arg 2 of next_sshu_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(area_t), C_LOC(area_t))
        CALL check_status('clSetKernelArg: arg 3 of next_sshu_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(area_u), C_LOC(area_u))
        CALL check_status('clSetKernelArg: arg 4 of next_sshu_code', ierr)
    END SUBROUTINE next_sshu_code_set_args

    SUBROUTINE next_sshv_code_set_args(kernel_obj, sshn_v, sshn_t, tmask, area_t, area_v)
        USE clfortran, ONLY: clSetKernelArg
        USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
        USE ocl_utils_mod, ONLY: check_status
        INTEGER(KIND=c_intptr_t), intent(in), target :: sshn_v, sshn_t, tmask, area_t, area_v
        INTEGER ierr
        INTEGER(KIND=c_intptr_t), target :: kernel_obj

        ! Set the arguments for the next_sshv_code OpenCL Kernel
        ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(sshn_v), C_LOC(sshn_v))
        CALL check_status('clSetKernelArg: arg 0 of next_sshv_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(sshn_t), C_LOC(sshn_t))
        CALL check_status('clSetKernelArg: arg 1 of next_sshv_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(tmask), C_LOC(tmask))
        CALL check_status('clSetKernelArg: arg 2 of next_sshv_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(area_t), C_LOC(area_t))
        CALL check_status('clSetKernelArg: arg 3 of next_sshv_code', ierr)
        ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(area_v), C_LOC(area_v))
        CALL check_status('clSetKernelArg: arg 4 of next_sshv_code', ierr)
    END SUBROUTINE next_sshv_code_set_args

    SUBROUTINE psy_init()
        USE fortcl, ONLY: ocl_env_init, add_kernels
        CHARACTER(LEN=30) kernel_names(11)
        LOGICAL, save :: initialised=.False.

        ! Check to make sure we only execute this routine once
        IF (.not. initialised) THEN
            initialised = .True.
            ! Initialise the OpenCL environment/device
            CALL ocl_env_init(1)
            ! The kernels this PSy layer module requires
            kernel_names(1) = "continuity_code"
            kernel_names(2) = "momentum_u_code"
            kernel_names(3) = "momentum_v_code"
            kernel_names(4) = "bc_ssh_code"
            kernel_names(5) = "bc_solid_u_code"
            kernel_names(6) = "bc_solid_v_code"
            kernel_names(7) = "bc_flather_u_code"
            kernel_names(8) = "bc_flather_v_code"
            kernel_names(9) = "field_copy_code"
            kernel_names(10) = "next_sshu_code"
            kernel_names(11) = "next_sshv_code"
            ! Create the OpenCL kernel objects. Expects to find all of the compiled
            ! kernels in PSYCLONE_KERNELS_FILE.
            CALL add_kernels(11, kernel_names)
        END IF
    END SUBROUTINE psy_init

end module time_step_mod
