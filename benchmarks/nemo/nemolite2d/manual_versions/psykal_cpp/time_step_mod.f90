module time_step_mod
    use field_mod
    use kind_params_mod
    implicit none
    public invoke_time_step

    ! Fortran to C wrapper interface using iso_c_bindings.
    interface
        subroutine wrapper_c_invoke_time_step( &
            ! Fields
            ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, &
            va, ssha_v, &
            ! Grid
            tmask, area_t, area_u, area_v, dx_u, dx_v, dx_t, dy_u, dy_v, &
            dy_t, gphiu, gphiv, &
            ! Scalars
            istp, internal_xstart, internal_xstop, internal_ystart, &
            internal_ystop, width, rdt, cbfr, visc, omega, d2r, g &
        ) bind (C, name="c_invoke_time_step")
            use iso_c_binding
            real(kind=c_double), intent(inout), dimension(*) :: ssha_t, &
                sshn_t, sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, &
                ssha_v, area_t, area_u, area_v, dx_u, dx_v, dx_t, dy_u, dy_v, &
                dy_t, gphiu, gphiv
            integer(kind=c_int), intent(inout), dimension(*) :: tmask
            integer(kind=c_int), intent(in), value :: istp, internal_xstart, &
                internal_xstop, internal_ystart, internal_ystop, width
            real(kind=c_double), intent(in), value :: rdt, cbfr, visc, omega, &
                d2r, g
        end subroutine wrapper_c_invoke_time_step
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
            ! Scalars
            istp, &
            ssha_t%grid%subdomain%internal%xstart - 1, & ! 1 -> 0 indexing
            ssha_t%grid%subdomain%internal%xstop - 1, & ! 1 -> 0 indexing
            ssha_t%grid%subdomain%internal%ystart - 1, & ! 1 -> 0 indexing
            ssha_t%grid%subdomain%internal%ystop - 1, & ! 1 -> 0 indexing
            size(ssha_t%data, 1), & ! Size of the contiguous dimension
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
            ssha_t%read_from_device_f => read_openmp
            sshn_t%read_from_device_f => read_openmp
            sshn_u%read_from_device_f => read_openmp
            sshn_v%read_from_device_f => read_openmp
            un%read_from_device_f => read_openmp
            vn%read_from_device_f => read_openmp
            ua%read_from_device_f => read_openmp
            ssha_u%read_from_device_f => read_openmp
            va%read_from_device_f => read_openmp
            ssha_v%read_from_device_f => read_openmp
        endif

    END SUBROUTINE invoke_time_step

  subroutine read_openmp(from, to, startx, starty, nx, ny, blocking)
    ! Function that specify how to retrieve the device data 'from' to a host
    ! location 'to', this function will be called by the infrastructure
    ! whenever the data is needed on the host.
    use iso_c_binding, only: c_ptr
    type(c_ptr), intent(in) :: from
    real(go_wp), dimension(:,:), intent(inout), target :: to
    integer, intent(in) :: startx, starty, nx, ny
    logical, intent(in) :: blocking

    ! Currently non-blocking reads are only requested by halo_exchanges which
    ! this manual version doesn't have, so it is safe to ignore. In the future
    ! we could use the 'async' clause if non-blocking synchronisations are
    ! needed.

    !$omp target update from(to)
  end subroutine read_openmp

end module time_step_mod
