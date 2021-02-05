! First PSy-layer level with a wrapper to the C PSy-layer using C_ISO_BINDINGS
! of the appropriate data.
module time_step_mod
    use field_mod
    use kind_params_mod
    implicit none
    public invoke_time_step

    ! Fortran to C wrapper interface using iso_c_bindings for the main invoke
    ! function.
    interface
        subroutine wrapper_c_invoke_time_step( &
            ! Fields
            ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, &
            va, ssha_v, &
            ! Grid
            tmask, area_t, area_u, area_v, dx_u, dx_v, dx_t, dy_u, dy_v, &
            dy_t, gphiu, gphiv, &
            ! Device pointers
            ssha_t_dp, sshn_t_dp, sshn_u_dp, sshn_v_dp, hu_dp, hv_dp, un_dp, &
            vn_dp, ua_dp, ht_dp, ssha_u_dp, va_dp, ssha_v_dp, &
            tmask_dp, area_t_dp, area_u_dp, area_v_dp, dx_u_dp, dx_v_dp, &
            dx_t_dp, dy_u_dp, dy_v_dp, dy_t_dp, gphiu_dp, gphiv_dp, &
            ! Scalars
            istp, internal_xstart, internal_xstop, internal_ystart, &
            internal_ystop, width, total_size, rdt, cbfr, visc, omega, d2r, g &
        ) bind (C, name="c_invoke_time_step")
            use iso_c_binding
            real(kind=c_double), intent(inout), dimension(*) :: ssha_t, &
                sshn_t, sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, &
                ssha_v, area_t, area_u, area_v, dx_u, dx_v, dx_t, dy_u, dy_v, &
                dy_t, gphiu, gphiv
            type(c_ptr), intent(inout) :: ssha_t_dp, sshn_t_dp, &
                sshn_u_dp, sshn_v_dp, hu_dp, hv_dp, un_dp, vn_dp, ua_dp, &
                ht_dp, ssha_u_dp, va_dp, ssha_v_dp
            integer(c_intptr_t), intent(inout) :: tmask_dp, area_t_dp, &
                area_u_dp, area_v_dp, dx_u_dp, dx_v_dp, dx_t_dp, dy_u_dp, &
                dy_v_dp, dy_t_dp, gphiu_dp, gphiv_dp
            integer(kind=c_int), intent(inout), dimension(*) :: tmask
            integer(kind=c_int), intent(in), value :: istp, internal_xstart, &
                internal_xstop, internal_ystart, internal_ystop, width, total_size
            real(kind=c_double), intent(in), value :: rdt, cbfr, visc, omega, &
                d2r, g
        end subroutine wrapper_c_invoke_time_step
    end interface    

    ! Fortran to C wrapper interface using iso_c_bindings for the function to
    ! read data from the device location 'from' to the host location 'to'
    interface
        subroutine wrapper_read_from_device(from, to, offset, nx, ny, &
                                            stride_gap) &
                bind(C, name="c_read_from_device")
            use iso_c_binding, only: c_ptr, c_int
            type(c_ptr), intent(in), value :: from
            type(c_ptr), intent(in), value :: to
            integer(c_int), intent(in), value :: offset, nx, ny, stride_gap
        end subroutine wrapper_read_from_device
    end interface

    interface
        subroutine write_to_device_c_interface(from, to, offset, nx, ny, &
                                               stride_gap) &
                bind(C, name="c_write_to_device")
            use iso_c_binding, only: c_int, c_ptr
            type(c_ptr), intent(in), value :: from
            type(c_ptr), intent(in), value :: to
            integer(c_int), intent(in), value :: offset, nx, ny, stride_gap
        end subroutine write_to_device_c_interface
    end interface
contains

    ! This invoke_time_step needs to fetch all necessary arrays, global
    ! variables, and the additional scalar values needed during the PSy-layer
    ! implementation, and pass them all to the wrapper interface.
    SUBROUTINE invoke_time_step(istp, ssha_t, ssha_u, ssha_v, sshn_t, sshn_u, &
        sshn_v, hu, hv, ht, ua, va, un, vn)

        USE physical_params_mod, ONLY: omega, d2r, g
        USE model_mod, ONLY: rdt, cbfr, visc

        implicit none

        TYPE(r2d_field), intent(inout) :: ssha_t, sshn_t, sshn_u, &
            sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
        INTEGER, intent(in) :: istp
        LOGICAL, save :: first_time=.true.

        ! TODO: issue #35 - Should this use %get_data() instead?
        call wrapper_c_invoke_time_step( &
            ! Fields
            ssha_t%data, &
            sshn_t%data, &
            sshn_u%data, &
            sshn_v%data, &
            hu%data, &
            hv%data, &
            un%data, &
            vn%data, &
            ua%data, &
            ht%data, &
            ssha_u%data, &
            va%data, &
            ssha_v%data, &
            ! Grid
            sshn_t%grid%tmask, &
            sshn_t%grid%area_t, &
            sshn_t%grid%area_u, &
            sshn_t%grid%area_v, &
            sshn_t%grid%dx_u, &
            sshn_t%grid%dx_v, &
            sshn_t%grid%dx_t, &
            sshn_t%grid%dy_u, &
            sshn_t%grid%dx_v, &
            sshn_t%grid%dx_t, &
            sshn_t%grid%gphiu, &
            sshn_t%grid%gphiv, &
            ! Field device pointers
            ssha_t%device_ptr, &
            sshn_t%device_ptr, &
            sshn_u%device_ptr, &
            sshn_v%device_ptr, &
            hu%device_ptr, &
            hv%device_ptr, &
            un%device_ptr, &
            vn%device_ptr, &
            ua%device_ptr, &
            ht%device_ptr, &
            ssha_u%device_ptr, &
            va%device_ptr, &
            ssha_v%device_ptr, &
            ! Grid device pointers
            sshn_t%grid%tmask_device, &
            sshn_t%grid%area_t_device, &
            sshn_t%grid%area_u_device, &
            sshn_t%grid%area_v_device, &
            sshn_t%grid%dx_u_device, &
            sshn_t%grid%dx_v_device, &
            sshn_t%grid%dx_t_device, &
            sshn_t%grid%dy_u_device, &
            sshn_t%grid%dx_v_device, &
            sshn_t%grid%dx_t_device, &
            sshn_t%grid%gphiu_device, &
            sshn_t%grid%gphiv_device, & 
            ! Scalars
            istp, &
            ssha_t%grid%subdomain%internal%xstart - 1, & ! 1 -> 0 indexing
            ssha_t%grid%subdomain%internal%xstop - 1, & ! 1 -> 0 indexing
            ssha_t%grid%subdomain%internal%ystart - 1, & ! 1 -> 0 indexing
            ssha_t%grid%subdomain%internal%ystop - 1, & ! 1 -> 0 indexing
            size(ssha_t%data, 1), & ! Size of the contiguous dimension
            size(ssha_t%data), & ! Total array size
            rdt, &
            cbfr, &
            visc, &
            omega, &
            d2r, &
            g &
        )

        if (first_time) then
            first_time = .false.
            ! Mark data_on_device flags
            ssha_t%data_on_device = .true.
            sshn_t%data_on_device = .true.
            sshn_u%data_on_device = .true.
            sshn_v%data_on_device = .true.
            hu%data_on_device = .true.
            hv%data_on_device = .true.
            un%data_on_device = .true.
            vn%data_on_device = .true.
            ua%data_on_device = .true.
            ht%data_on_device = .true.
            ssha_u%data_on_device = .true.
            va%data_on_device = .true.
            ssha_v%data_on_device = .true.

            ! Specify device data retrieving methods
            ssha_t%read_from_device_c => wrapper_read_from_device
            sshn_t%read_from_device_c => wrapper_read_from_device
            sshn_u%read_from_device_c => wrapper_read_from_device
            sshn_v%read_from_device_c => wrapper_read_from_device
            hu%read_from_device_c => wrapper_read_from_device
            hv%read_from_device_c => wrapper_read_from_device
            un%read_from_device_c => wrapper_read_from_device
            vn%read_from_device_c => wrapper_read_from_device
            ua%read_from_device_c => wrapper_read_from_device
            ht%read_from_device_c => wrapper_read_from_device
            ssha_u%read_from_device_c => wrapper_read_from_device
            va%read_from_device_c => wrapper_read_from_device
            ssha_v%read_from_device_c => wrapper_read_from_device
        endif

    END SUBROUTINE invoke_time_step
end module time_step_mod
