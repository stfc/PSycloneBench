module time_step_mod
    use field_mod
    use kind_params_mod
    implicit none
    private
    public invoke_time_step
contains

    ! This is a manual version of the OpenCL PSy-layer. It executes the
    ! following workflow:
    ! 1 ) The first time its called, initialize OpenCL and populate the
    !     necessary resources (e.g. queues, kernel functions).
    ! 2 ) Calls set_arguments to set up each OpenCL kernel argument list.
    ! 3 ) If the device does not contain up-to-date data of the fields and
    !     scalars used, it copies them into the device.
    ! 4 ) Enqueue functions for each of the kernels to execute.
    ! 5 ) Barrier at the end to wait for all enqueued tasks to finish.
    SUBROUTINE invoke_time_step(ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, &
                                vn, ua, ht, ssha_u, va, ssha_v, istp)
        USE field_mod
        USE kind_params_mod
        USE fortcl, ONLY: create_rw_buffer
        USE fortcl, ONLY: get_num_cmd_queues, get_cmd_queues, &
                          get_kernel_by_name
        USE clfortran
        USE ocl_utils_mod, ONLY: check_status
        USE iso_c_binding
        USE physical_params_mod, ONLY: omega, d2r, g
        USE model_mod, ONLY: rdt, cbfr, visc
        TYPE(r2d_field), intent(inout), target :: ssha_t, sshn_t, sshn_u, &
            sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
        INTEGER, intent(inout) :: istp
        INTEGER(KIND=c_size_t), target :: localsize(2)
        INTEGER(KIND=c_size_t), target :: globalsize(2)
        INTEGER(KIND=c_intptr_t), target, save :: kernel_next_sshv_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_next_sshu_code
        INTEGER(KIND=c_intptr_t), target, save :: kernel_field_copy_code1
        INTEGER(KIND=c_intptr_t), target, save :: kernel_field_copy_code2
        INTEGER(KIND=c_intptr_t), target, save :: kernel_field_copy_code3
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
        INTEGER(KIND=c_intptr_t) :: ssha_t_device, sshn_t_device, sshn_u_device, &
            sshn_v_device, hu_device, hv_device, un_device, vn_device, ua_device, &
            ht_device, ssha_u_device, va_device, ssha_v_device
        INTEGER, save :: num_cmd_queues
        type(c_ptr) :: swap
        logical, parameter :: no_copy_optimization = .False.
        logical, parameter :: tasks_optimization = .False.

        IF (first_time) THEN
            ! Ensure OpenCL run-time is initialised for this PSy-layer module
            call psy_init()
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
            kernel_field_copy_code1 = get_kernel_by_name("field_copy_code")
            kernel_field_copy_code2 = get_kernel_by_name("field_copy_code")
            kernel_field_copy_code3 = get_kernel_by_name("field_copy_code")
            kernel_next_sshu_code = get_kernel_by_name("next_sshu_code")
            kernel_next_sshv_code = get_kernel_by_name("next_sshv_code")

            ! Initialise device buffers if they don't already exist
            call initialise_device_buffer(ssha_t)
            call initialise_device_buffer(sshn_t)
            call initialise_device_buffer(sshn_u)
            call initialise_device_buffer(sshn_v)
            call initialise_device_buffer(hu)
            call initialise_device_buffer(hv)
            call initialise_device_buffer(un)
            call initialise_device_buffer(vn)
            call initialise_device_buffer(ua)
            call initialise_device_buffer(ht)
            call initialise_device_buffer(ssha_u)
            call initialise_device_buffer(va)
            call initialise_device_buffer(ssha_v)
            ! Grid is shared, it only needs to be called for one field
            call initialise_device_grid(ssha_t)
        endif

        ! Cast the C pointer to the cl_mem type expected by OpenCL
        ssha_t_device = transfer(ssha_t%device_ptr, ssha_t_device)
        sshn_t_device = transfer(sshn_t%device_ptr, sshn_t_device)
        sshn_u_device = transfer(sshn_u%device_ptr, sshn_u_device)
        sshn_v_device = transfer(sshn_v%device_ptr, sshn_v_device)
        hu_device = transfer(hu%device_ptr, hu_device)
        hv_device = transfer(hv%device_ptr, hv_device)
        un_device = transfer(un%device_ptr, un_device)
        vn_device = transfer(vn%device_ptr, vn_device)
        ua_device = transfer(ua%device_ptr, ua_device)
        ht_device = transfer(ht%device_ptr, ht_device)
        ssha_u_device = transfer(ssha_u%device_ptr, ssha_u_device)
        va_device = transfer(va%device_ptr, va_device)
        ssha_v_device = transfer(ssha_v%device_ptr, ssha_v_device)

        IF (first_time .or. no_copy_optimization) THEN
            ! Before writing into the device we need to execute the set_args
            ! subroutines because some architectures (e.g. Xilinx FPGA) use
            ! this information for the memory placement into different device
            ! memory banks.
            CALL continuity_code_set_args(kernel_continuity_code, &
                ssha_t%internal%xstart - 1, ssha_t%internal%xstop - 1, &
                ssha_t%internal%ystart - 1, ssha_t%internal%ystop - 1, &
                ssha_t_device, sshn_t_device, sshn_u_device, sshn_v_device, &
                hu_device, hv_device, un_device, vn_device, &
                sshn_t%grid%area_t_device, rdt)

            CALL momentum_u_code_set_args(kernel_momentum_u_code, &
                ua%internal%xstart - 1, ua%internal%xstop - 1, &
                ua%internal%ystart - 1, ua%internal%ystop - 1, &
                ua_device, un_device, vn_device, hu_device, hv_device, &
                ht_device, ssha_u_device, sshn_t_device, sshn_u_device, &
                sshn_v_device, un%grid%tmask_device, un%grid%dx_u_device, &
                un%grid%dx_v_device, un%grid%dx_t_device, un%grid%dy_u_device, &
                un%grid%dy_t_device, un%grid%area_u_device, &
                un%grid%gphiu_device, omega, d2r, g, rdt, cbfr, visc)

            CALL momentum_v_code_set_args(kernel_momentum_v_code, &
                va%internal%xstart - 1, va%internal%xstop - 1, &
                va%internal%ystart - 1, va%internal%ystop - 1, va_device, &
                un_device, vn_device, hu_device, hv_device, ht_device, &
                ssha_v_device, sshn_t_device, sshn_u_device, sshn_v_device, &
                un%grid%tmask_device, un%grid%dx_v_device, un%grid%dx_t_device, &
                un%grid%dy_u_device, un%grid%dy_v_device, un%grid%dy_t_device, &
                un%grid%area_v_device, un%grid%gphiv_device, omega, d2r, g, &
                rdt, cbfr, visc)

            CALL bc_ssh_code_set_args(kernel_bc_ssh_code, &
                ssha_t%internal%xstart - 1, ssha_t%internal%xstop - 1, &
                ssha_t%internal%ystart - 1, ssha_t%internal%ystop - 1, istp, &
                ssha_t_device, ssha_t%grid%tmask_device, rdt)

            CALL bc_solid_u_code_set_args(kernel_bc_solid_u_code, &
                ua%whole%xstart - 1, ua%whole%xstop - 1, ua%whole%ystart - 1, &
                ua%whole%ystop - 1, ua_device, ua%grid%tmask_device)

            CALL bc_solid_v_code_set_args(kernel_bc_solid_v_code, &
                va%whole%xstart - 1, va%whole%xstop - 1, va%whole%ystart - 1, &
                va%whole%ystop - 1, va_device, va%grid%tmask_device)

            CALL bc_flather_u_code_set_args(kernel_bc_flather_u_code, &
                ua%whole%xstart - 1, ua%whole%xstop - 1, ua%whole%ystart - 1, &
                ua%whole%ystop - 1, ua_device, hu_device, sshn_u_device, &
                hu%grid%tmask_device, g)

            CALL bc_flather_v_code_set_args(kernel_bc_flather_v_code, &
                va%whole%xstart - 1, va%whole%xstop - 1, va%whole%ystart - 1, &
                va%whole%ystop - 1, va_device, hv_device, sshn_v_device, &
                hv%grid%tmask_device, g)

            CALL next_sshu_code_set_args(kernel_next_sshu_code, &
                sshn_u%internal%xstart - 1, sshn_u%internal%xstop - 1, &
                sshn_u%internal%ystart - 1, sshn_u%internal%ystop - 1, &
                sshn_u_device, ssha_t_device, sshn_t%grid%tmask_device, &
                sshn_t%grid%area_t_device, sshn_t%grid%area_u_device)

            CALL next_sshv_code_set_args(kernel_next_sshv_code, &
                sshn_v%internal%xstart - 1, sshn_v%internal%xstop - 1, &
                sshn_v%internal%ystart - 1, sshn_v%internal%ystop - 1, &
                sshn_v_device, ssha_t_device, sshn_t%grid%tmask_device, &
                sshn_t%grid%area_t_device, sshn_t%grid%area_v_device)

            CALL field_copy_code_set_args(kernel_field_copy_code1, &
                0, SIZE(un%data, 1) - 1, 0, SIZE(un%data, 2) - 1, &
                un_device, ua_device)

            CALL field_copy_code_set_args(kernel_field_copy_code2, &
                0, SIZE(vn%data, 1) - 1, 0, SIZE(vn%data, 2) - 1, &
                vn_device, va_device)

            CALL field_copy_code_set_args(kernel_field_copy_code3, &
                0, SIZE(sshn_t%data, 1) - 1, 0, SIZE(sshn_t%data, 2) - 1, &
                sshn_t_device, ssha_t_device)
        ENDIF

        IF (first_time) THEN
            ! Write initial data to device
            call ssha_t%write_to_device()
            call sshn_t%write_to_device()
            call sshn_u%write_to_device()
            call sshn_v%write_to_device()
            call hu%write_to_device()
            call hv%write_to_device()
            call un%write_to_device()
            call vn%write_to_device()
            call ua%write_to_device()
            call ht%write_to_device()
            call ssha_u%write_to_device()
            call va%write_to_device()
            call ssha_v%write_to_device()
            ! Grid is shared, it only needs to be called for one field
            call write_device_grid(ssha_t)
        END IF

        ! Set up local and global sizes
        if(tasks_optimization) then
            ! Only one task enqueued if using tasks optimization
            globalsize = (/1, 1/)
            localsize = (/1, 1/)
        else
            ! Groups of 64 contiguous ids for the whole domain if using NDRanges
            globalsize = (/sshn_t%grid%nx, sshn_t%grid%ny/)
            localsize = (/64, 1/)
            if (mod(globalsize(1), localsize(1)) /= 0) then
                stop "OpenCL nx global size must be a multiple of 64, " // &
                     "change the problem size or use DL_ESM_ALIGNMENT=64"
            endif
        endif

        ! Launch the kernels
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_continuity_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            C_NULL_PTR, C_NULL_PTR)
        ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernel_momentum_u_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            C_NULL_PTR, C_NULL_PTR)
        ierr = clEnqueueNDRangeKernel(cmd_queues(3), kernel_momentum_v_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            C_NULL_PTR, C_NULL_PTR)
        ! Needs to be done every time because of the changing istp value
        CALL bc_ssh_code_set_args(kernel_bc_ssh_code, &
            ssha_t%internal%xstart - 1, ssha_t%internal%xstop - 1, &
            ssha_t%internal%ystart - 1, ssha_t%internal%ystop - 1, istp, &
            ssha_t_device, ssha_t%grid%tmask_device, rdt)
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_ssh_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            C_NULL_PTR, C_NULL_PTR)
        ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernel_bc_solid_u_code, 2, &
            C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            C_NULL_PTR, C_NULL_PTR)
        ierr = clEnqueueNDRangeKernel(cmd_queues(3), kernel_bc_solid_v_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            C_NULL_PTR, C_NULL_PTR)
        ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernel_bc_flather_u_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), &
            0, C_NULL_PTR, C_NULL_PTR)
        ierr = clEnqueueNDRangeKernel(cmd_queues(3), kernel_bc_flather_v_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), &
            0, C_NULL_PTR, C_NULL_PTR)

        ierr = clFinish(cmd_queues(1))
        ierr = clFinish(cmd_queues(2))
        ierr = clFinish(cmd_queues(3))

        ! Next kernels can be brought before the Copy kernels by switching
        ! appropriately sshn_t_device <-> ssha_t_device in the set_args
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_next_sshu_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            C_NULL_PTR, C_NULL_PTR)
        ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_next_sshv_code, &
            2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
            C_NULL_PTR, C_NULL_PTR)

        if (.not. no_copy_optimization) then
            ! The set_args subroutine need to be called every time because the
            ! field_copy_code kernel is reused with 3 different argument sets.
            ! Alternatively we could clone the kernel 3 times in allkernels.cl
            ! with different names.
            CALL field_copy_code_set_args(kernel_field_copy_code1, &
                0, SIZE(un%data, 1) - 1, 0, SIZE(un%data, 2) - 1, &
                un_device, ua_device)
            ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_field_copy_code1, &
                2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
                C_NULL_PTR, C_NULL_PTR)
            CALL field_copy_code_set_args(kernel_field_copy_code2, &
                0, SIZE(vn%data, 1) - 1, 0, SIZE(vn%data, 2) - 1, &
                vn_device, va_device)
            ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_field_copy_code2, &
                2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
                C_NULL_PTR, C_NULL_PTR)
            CALL field_copy_code_set_args(kernel_field_copy_code3, &
                0, SIZE(sshn_t%data, 1) - 1, 0, SIZE(sshn_t%data, 2) - 1, &
                sshn_t_device, ssha_t_device)
            ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_field_copy_code3, &
                2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
                C_NULL_PTR, C_NULL_PTR)
        else
            swap = ssha_t%device_ptr
            ssha_t%device_ptr = sshn_t%device_ptr
            sshn_t%device_ptr = swap

            swap = ua%device_ptr
            ua%device_ptr = un%device_ptr
            un%device_ptr = swap

            swap = va%device_ptr
            va%device_ptr = vn%device_ptr
            vn%device_ptr = swap
        endif

        ! Block until all kernels have finished
        ierr = clFinish(cmd_queues(1))
        CALL check_status('Barrier end of iteration', ierr)

        if (first_time) then
            first_time = .false.
        endif
    END SUBROUTINE invoke_time_step

    subroutine initialise_device_buffer(field)
        USE fortcl, ONLY: create_rw_buffer
        type(r2d_field), intent(inout), target :: field
        integer(kind=c_size_t) size_in_bytes
        IF (.NOT. field%data_on_device) THEN
            size_in_bytes = int(field%grid%nx*field%grid%ny, 8) * &
                            c_sizeof(field%data(1,1))
            ! Create buffer on device
            field%device_ptr = transfer(create_rw_buffer(size_in_bytes), &
                                         field%device_ptr)
            field%data_on_device = .true.
            field%read_from_device_f => read_opencl
            field%write_to_device_f => write_opencl
        END IF
    end subroutine initialise_device_buffer

    ! This OpenCL manual implementation only ever read/writes whole buffers
    ! at once, so the coarse-grain (full rows) read/write functions below
    ! already provide sufficient performance.
    subroutine read_opencl(from, to, startx, starty, nx, ny, blocking)
        use iso_c_binding, only: c_ptr, c_intptr_t, c_size_t, c_sizeof
        USE ocl_utils_mod, ONLY: check_status
        use kind_params_mod, only: go_wp
        USE clfortran
        USE fortcl, ONLY: get_cmd_queues
        type(c_ptr), intent(in) :: from
        real(go_wp), dimension(:,:), intent(inout), target :: to
        integer, intent(in) :: startx, starty, nx, ny
        logical, intent(in) :: blocking
        INTEGER(c_size_t) :: size_in_bytes, offset_in_bytes
        integer(c_intptr_t) :: cl_mem
        INTEGER(c_intptr_t), pointer :: cmd_queues(:)
        integer :: ierr
        ! Copy complete ny rows (regardless of nx)
        size_in_bytes = int(size(to, 1) * ny, 8) * c_sizeof(to(1,1))
        offset_in_bytes = int(size(to, 1) * (starty - 1), 8) * c_sizeof(to(1,1))
        cl_mem = transfer(from, cl_mem)
        cmd_queues => get_cmd_queues()
        ierr = clEnqueueReadBuffer(cmd_queues(1), cl_mem, &
            CL_TRUE, offset_in_bytes, size_in_bytes, C_LOC(to), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueReadBuffer', ierr)
    end subroutine read_opencl

    subroutine write_opencl(from, to, startx, starty, nx, ny, blocking)
        use iso_c_binding, only: c_ptr, c_intptr_t, c_size_t, c_sizeof
        USE ocl_utils_mod, ONLY: check_status
        use kind_params_mod, only: go_wp
        USE clfortran
        USE fortcl, ONLY: get_cmd_queues
        real(go_wp), dimension(:,:), intent(in), target :: from
        type(c_ptr), intent(in) :: to
        integer, intent(in) :: startx, starty, nx, ny
        logical, intent(in) :: blocking
        integer(c_intptr_t) :: cl_mem
        INTEGER(c_size_t) :: size_in_bytes, offset_in_bytes
        INTEGER(c_intptr_t), pointer :: cmd_queues(:)
        integer :: ierr
        ! Copy complete ny rows (regardless of nx)
        size_in_bytes = int(size(from, 1) * ny, 8) * c_sizeof(from(1,1))
        offset_in_bytes = int(size(from, 1) * (starty - 1)) * c_sizeof(from(1,1))
        cl_mem = transfer(to, cl_mem)
        cmd_queues => get_cmd_queues()
        ierr = clEnqueueWriteBuffer(cmd_queues(1), cl_mem, &
            CL_TRUE, offset_in_bytes, size_in_bytes, C_LOC(from), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer', ierr)
    end subroutine write_opencl

    subroutine initialise_device_grid(field)
        USE fortcl, ONLY: create_rw_buffer
        type(r2d_field), intent(inout), target :: field
        integer(kind=c_size_t) size_in_bytes
        IF (field%grid%tmask_device == 0) THEN
            ! Create integer grid fields
            size_in_bytes = int(field%grid%nx * field%grid%ny, 8) * &
                            c_sizeof(field%grid%tmask(1,1))
            field%grid%tmask_device = create_rw_buffer(size_in_bytes)

            ! Create real grid buffers
            size_in_bytes = int(field%grid%nx * field%grid%ny, 8) * &
                            c_sizeof(field%grid%area_t(1,1))
            field%grid%area_t_device = create_rw_buffer(size_in_bytes)
            field%grid%area_u_device = create_rw_buffer(size_in_bytes)
            field%grid%area_v_device = create_rw_buffer(size_in_bytes)
            field%grid%dx_u_device = create_rw_buffer(size_in_bytes)
            field%grid%dx_v_device = create_rw_buffer(size_in_bytes)
            field%grid%dx_t_device = create_rw_buffer(size_in_bytes)
            field%grid%dy_u_device = create_rw_buffer(size_in_bytes)
            field%grid%dy_v_device = create_rw_buffer(size_in_bytes)
            field%grid%dy_t_device = create_rw_buffer(size_in_bytes)
            field%grid%gphiu_device = create_rw_buffer(size_in_bytes)
            field%grid%gphiv_device = create_rw_buffer(size_in_bytes)
        END IF
    end subroutine initialise_device_grid

    subroutine write_device_grid(field)
        USE fortcl, ONLY: get_cmd_queues
        USE clfortran
        USE ocl_utils_mod, ONLY: check_status
        type(r2d_field), intent(inout), target :: field
        integer(kind=c_size_t) size_in_bytes
        INTEGER(c_intptr_t), pointer :: cmd_queues(:)
        integer :: ierr
        cmd_queues => get_cmd_queues()
        ! Integer grid buffers
        size_in_bytes = int(field%grid%nx * field%grid%ny, 8) * &
                            c_sizeof(field%grid%tmask(1,1))
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%tmask_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%tmask), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer tmask', ierr)

        ! Real grid buffers
        size_in_bytes = int(field%grid%nx * field%grid%ny, 8) * &
                            c_sizeof(field%grid%area_t(1,1))
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%area_t_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%area_t), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer area_t_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%area_u_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%area_u), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer area_u_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%area_v_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%area_v), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer area_v_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%dx_u_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%dx_u), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer dx_u_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%dx_v_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%dx_v), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer dx_v_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%dx_t_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%dx_t), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer dx_t_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%dy_u_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%dy_u), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer dy_u_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%dy_v_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%dy_v), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer dy_v_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%dy_t_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%dy_t), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer dy_t_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%gphiu_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%gphiu), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer gphiu_device', ierr)
        ierr = clEnqueueWriteBuffer(cmd_queues(1), field%grid%gphiv_device, &
            CL_TRUE, 0_8, size_in_bytes, C_LOC(field%grid%gphiv), 0, &
            C_NULL_PTR, C_NULL_PTR)
        CALL check_status('clEnqueueWriteBuffer gphiv_device', ierr)
    end subroutine write_device_grid

    SUBROUTINE continuity_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, &
&area_t, rdt)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, area_t
      REAL(KIND=go_wp), intent(in), target :: rdt
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the continuity_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(ssha_t), C_LOC(ssha_t))
      CALL check_status('clSetKernelArg: arg 4 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 5 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 6 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 7 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(hu), C_LOC(hu))
      CALL check_status('clSetKernelArg: arg 8 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 9, C_SIZEOF(hv), C_LOC(hv))
      CALL check_status('clSetKernelArg: arg 9 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 10, C_SIZEOF(un), C_LOC(un))
      CALL check_status('clSetKernelArg: arg 10 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 11, C_SIZEOF(vn), C_LOC(vn))
      CALL check_status('clSetKernelArg: arg 11 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 12, C_SIZEOF(area_t), C_LOC(area_t))
      CALL check_status('clSetKernelArg: arg 12 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 13, C_SIZEOF(rdt), C_LOC(rdt))
      CALL check_status('clSetKernelArg: arg 13 of continuity_code', ierr)
    END SUBROUTINE continuity_code_set_args
    SUBROUTINE momentum_u_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, ua, un, vn, hu, hv, ht, ssha_u, sshn_t, sshn_u, &
&sshn_v, tmask, dx_u, dx_v, dx_t, dy_u, dy_t, area_u, gphiu, omega, d2r, g, rdt, cbfr, visc)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ua, un, vn, hu, hv, ht, ssha_u, sshn_t, sshn_u, sshn_v, tmask, dx_u, dx_v, &
&dx_t, dy_u, dy_t, area_u, gphiu
      REAL(KIND=go_wp), intent(in), target :: omega
      REAL(KIND=go_wp), intent(in), target :: d2r
      REAL(KIND=go_wp), intent(in), target :: g
      REAL(KIND=go_wp), intent(in), target :: rdt
      REAL(KIND=go_wp), intent(in), target :: cbfr
      REAL(KIND=go_wp), intent(in), target :: visc
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the momentum_u_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(ua), C_LOC(ua))
      CALL check_status('clSetKernelArg: arg 4 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(un), C_LOC(un))
      CALL check_status('clSetKernelArg: arg 5 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(vn), C_LOC(vn))
      CALL check_status('clSetKernelArg: arg 6 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(hu), C_LOC(hu))
      CALL check_status('clSetKernelArg: arg 7 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(hv), C_LOC(hv))
      CALL check_status('clSetKernelArg: arg 8 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 9, C_SIZEOF(ht), C_LOC(ht))
      CALL check_status('clSetKernelArg: arg 9 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 10, C_SIZEOF(ssha_u), C_LOC(ssha_u))
      CALL check_status('clSetKernelArg: arg 10 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 11, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 11 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 12, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 12 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 13, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 13 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 14, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 14 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 15, C_SIZEOF(dx_u), C_LOC(dx_u))
      CALL check_status('clSetKernelArg: arg 15 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 16, C_SIZEOF(dx_v), C_LOC(dx_v))
      CALL check_status('clSetKernelArg: arg 16 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 17, C_SIZEOF(dx_t), C_LOC(dx_t))
      CALL check_status('clSetKernelArg: arg 17 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 18, C_SIZEOF(dy_u), C_LOC(dy_u))
      CALL check_status('clSetKernelArg: arg 18 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 19, C_SIZEOF(dy_t), C_LOC(dy_t))
      CALL check_status('clSetKernelArg: arg 19 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 20, C_SIZEOF(area_u), C_LOC(area_u))
      CALL check_status('clSetKernelArg: arg 20 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 21, C_SIZEOF(gphiu), C_LOC(gphiu))
      CALL check_status('clSetKernelArg: arg 21 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 22, C_SIZEOF(omega), C_LOC(omega))
      CALL check_status('clSetKernelArg: arg 22 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 23, C_SIZEOF(d2r), C_LOC(d2r))
      CALL check_status('clSetKernelArg: arg 23 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 24, C_SIZEOF(g), C_LOC(g))
      CALL check_status('clSetKernelArg: arg 24 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 25, C_SIZEOF(rdt), C_LOC(rdt))
      CALL check_status('clSetKernelArg: arg 25 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 26, C_SIZEOF(cbfr), C_LOC(cbfr))
      CALL check_status('clSetKernelArg: arg 26 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 27, C_SIZEOF(visc), C_LOC(visc))
      CALL check_status('clSetKernelArg: arg 27 of momentum_u_code', ierr)
    END SUBROUTINE momentum_u_code_set_args
    SUBROUTINE momentum_v_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, va, un, vn, hu, hv, ht, ssha_v, sshn_t, sshn_u, &
&sshn_v, tmask, dx_v, dx_t, dy_u, dy_v, dy_t, area_v, gphiv, omega, d2r, g, rdt, cbfr, visc)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: va, un, vn, hu, hv, ht, ssha_v, sshn_t, sshn_u, sshn_v, tmask, dx_v, dx_t, &
&dy_u, dy_v, dy_t, area_v, gphiv
      REAL(KIND=go_wp), intent(in), target :: omega
      REAL(KIND=go_wp), intent(in), target :: d2r
      REAL(KIND=go_wp), intent(in), target :: g
      REAL(KIND=go_wp), intent(in), target :: rdt
      REAL(KIND=go_wp), intent(in), target :: cbfr
      REAL(KIND=go_wp), intent(in), target :: visc
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the momentum_v_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(va), C_LOC(va))
      CALL check_status('clSetKernelArg: arg 4 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(un), C_LOC(un))
      CALL check_status('clSetKernelArg: arg 5 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(vn), C_LOC(vn))
      CALL check_status('clSetKernelArg: arg 6 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(hu), C_LOC(hu))
      CALL check_status('clSetKernelArg: arg 7 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(hv), C_LOC(hv))
      CALL check_status('clSetKernelArg: arg 8 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 9, C_SIZEOF(ht), C_LOC(ht))
      CALL check_status('clSetKernelArg: arg 9 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 10, C_SIZEOF(ssha_v), C_LOC(ssha_v))
      CALL check_status('clSetKernelArg: arg 10 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 11, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 11 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 12, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 12 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 13, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 13 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 14, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 14 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 15, C_SIZEOF(dx_v), C_LOC(dx_v))
      CALL check_status('clSetKernelArg: arg 15 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 16, C_SIZEOF(dx_t), C_LOC(dx_t))
      CALL check_status('clSetKernelArg: arg 16 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 17, C_SIZEOF(dy_u), C_LOC(dy_u))
      CALL check_status('clSetKernelArg: arg 17 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 18, C_SIZEOF(dy_v), C_LOC(dy_v))
      CALL check_status('clSetKernelArg: arg 18 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 19, C_SIZEOF(dy_t), C_LOC(dy_t))
      CALL check_status('clSetKernelArg: arg 19 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 20, C_SIZEOF(area_v), C_LOC(area_v))
      CALL check_status('clSetKernelArg: arg 20 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 21, C_SIZEOF(gphiv), C_LOC(gphiv))
      CALL check_status('clSetKernelArg: arg 21 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 22, C_SIZEOF(omega), C_LOC(omega))
      CALL check_status('clSetKernelArg: arg 22 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 23, C_SIZEOF(d2r), C_LOC(d2r))
      CALL check_status('clSetKernelArg: arg 23 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 24, C_SIZEOF(g), C_LOC(g))
      CALL check_status('clSetKernelArg: arg 24 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 25, C_SIZEOF(rdt), C_LOC(rdt))
      CALL check_status('clSetKernelArg: arg 25 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 26, C_SIZEOF(cbfr), C_LOC(cbfr))
      CALL check_status('clSetKernelArg: arg 26 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 27, C_SIZEOF(visc), C_LOC(visc))
      CALL check_status('clSetKernelArg: arg 27 of momentum_v_code', ierr)
    END SUBROUTINE momentum_v_code_set_args
    SUBROUTINE bc_ssh_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, istp, ssha_t, tmask, rdt)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ssha_t, tmask
      INTEGER, intent(in), target :: istp
      REAL(KIND=go_wp), intent(in), target :: rdt
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_ssh_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(istp), C_LOC(istp))
      CALL check_status('clSetKernelArg: arg 4 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(ssha_t), C_LOC(ssha_t))
      CALL check_status('clSetKernelArg: arg 5 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 6 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(rdt), C_LOC(rdt))
      CALL check_status('clSetKernelArg: arg 7 of bc_ssh_code', ierr)
    END SUBROUTINE bc_ssh_code_set_args
    SUBROUTINE bc_solid_u_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, ua, tmask)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ua, tmask
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_solid_u_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(ua), C_LOC(ua))
      CALL check_status('clSetKernelArg: arg 4 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 5 of bc_solid_u_code', ierr)
    END SUBROUTINE bc_solid_u_code_set_args
    SUBROUTINE bc_solid_v_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, va, tmask)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: va, tmask
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_solid_v_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(va), C_LOC(va))
      CALL check_status('clSetKernelArg: arg 4 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 5 of bc_solid_v_code', ierr)
    END SUBROUTINE bc_solid_v_code_set_args
    SUBROUTINE bc_flather_u_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, ua, hu, sshn_u, tmask, g)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ua, hu, sshn_u, tmask
      REAL(KIND=go_wp), intent(in), target :: g
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_flather_u_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(ua), C_LOC(ua))
      CALL check_status('clSetKernelArg: arg 4 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(hu), C_LOC(hu))
      CALL check_status('clSetKernelArg: arg 5 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 6 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 7 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(g), C_LOC(g))
      CALL check_status('clSetKernelArg: arg 8 of bc_flather_u_code', ierr)
    END SUBROUTINE bc_flather_u_code_set_args
    SUBROUTINE bc_flather_v_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, va, hv, sshn_v, tmask, g)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: va, hv, sshn_v, tmask
      REAL(KIND=go_wp), intent(in), target :: g
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_flather_v_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(va), C_LOC(va))
      CALL check_status('clSetKernelArg: arg 4 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(hv), C_LOC(hv))
      CALL check_status('clSetKernelArg: arg 5 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 6 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 7 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(g), C_LOC(g))
      CALL check_status('clSetKernelArg: arg 8 of bc_flather_v_code', ierr)
    END SUBROUTINE bc_flather_v_code_set_args
    SUBROUTINE field_copy_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, un, ua)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: un, ua
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the field_copy_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(un), C_LOC(un))
      CALL check_status('clSetKernelArg: arg 4 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(ua), C_LOC(ua))
      CALL check_status('clSetKernelArg: arg 5 of field_copy_code', ierr)
    END SUBROUTINE field_copy_code_set_args
    SUBROUTINE next_sshu_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, sshn_u, sshn_t, tmask, area_t, area_u)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: sshn_u, sshn_t, tmask, area_t, area_u
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the next_sshu_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 4 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 5 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 6 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(area_t), C_LOC(area_t))
      CALL check_status('clSetKernelArg: arg 7 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(area_u), C_LOC(area_u))
      CALL check_status('clSetKernelArg: arg 8 of next_sshu_code', ierr)
    END SUBROUTINE next_sshu_code_set_args
    SUBROUTINE next_sshv_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, sshn_v, sshn_t, tmask, area_t, area_v)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: sshn_v, sshn_t, tmask, area_t, area_v
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the next_sshv_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 4 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 5 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 6 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(area_t), C_LOC(area_t))
      CALL check_status('clSetKernelArg: arg 7 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(area_v), C_LOC(area_v))
      CALL check_status('clSetKernelArg: arg 8 of next_sshv_code', ierr)
    END SUBROUTINE next_sshv_code_set_args
 
    SUBROUTINE psy_init()
        USE fortcl, ONLY: ocl_env_init, add_kernels
        CHARACTER(LEN=30) kernel_names(11)
        LOGICAL, save :: initialised=.False.

        ! Check to make sure we only execute this routine once
        IF (.not. initialised) THEN
            initialised = .True.
            ! Initialise the OpenCL environment/device
            ! parameters are: number of command queues, device selection,
            !                 profiling, out_of_order
            CALL ocl_env_init(3, 1, .false., .false.)
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
            ! kernels in FORTCL_KERNELS_FILE.
            CALL add_kernels(11, kernel_names)
        END IF
    END SUBROUTINE psy_init

end module time_step_mod
