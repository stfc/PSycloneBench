module time_step_mod
    use field_mod
    use kind_params_mod
    implicit none
    public invoke_time_step

    ! Signature of the C++ function we want to call
    ! extern "C" void c_invoke_time_step(double * c_ssha_t, double * c_sshn_t,
    !   double * c_sshn_u, double * c_sshn_v, double * c_hu, double * c_hv,
    !   double * c_hn, double * c_vn, double * c_ua, double * c_ht,
    !   double * c_ssha_u, double * c_va, double * c_ssha_v, int istp, int nsize)

    ! Fortran to C wrapper function as recommended in
    ! http://fortranwiki.org/fortran/show/Generating+C+Interfaces
    interface
        subroutine wrapper_c_invoke_time_step( &
            ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, &
            un, vn, ua, ht, ssha_u, va, ssha_v, istp, nx, ny &
        ) bind (C, name="c_invoke_time_step")
            use iso_c_binding
            real(kind=c_double), intent(inout), dimension(*) :: ssha_t, sshn_t, &
                sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
            integer(c_int), intent(in), value :: istp, nx, ny
        end subroutine wrapper_c_invoke_time_step
    end interface    

contains

    ! invoke_time_step needs to fetch all data as raw c arrays, the necessary
    ! global variables, and the additional values needed during the PSy-layer
    ! implementation, and pass them all to the C++ function using proper
    ! iso_c_bindings.
    SUBROUTINE invoke_time_step(istp, ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, &
            un, vn, ua, ht, ssha_u, va, ssha_v)
        USE iso_c_binding
        USE field_mod
        USE kind_params_mod
        USE physical_params_mod, ONLY: omega, d2r, g
        USE model_mod, ONLY: rdt, cbfr, visc

        implicit none

        TYPE(r2d_field), intent(inout), target :: ssha_t, sshn_t, sshn_u, &
            sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
        INTEGER, intent(in) :: istp
        INTEGER:: istop, jstop, nx, ny


        LOGICAL, save :: first_time=.true.

        ! Look-up loop bounds
        istop = ssha_t%grid%subdomain%internal%xstop
        jstop = ssha_t%grid%subdomain%internal%ystop
        nx = ssha_t%grid%nx
        ny = ssha_t%grid%ny

        ! Any steps that need to be done just the first time
        IF (first_time) THEN
            first_time = .false.
        END IF

        ! FIXME: Should we use %get_data() instead? The dl_esm_inf has some
        ! infrastructure for device_ptr and dirty data to manage when data is
        ! modified that may be handy here.
        call wrapper_c_invoke_time_step( &
            ssha_t%data, sshn_t%data, sshn_u%data, sshn_v%data, hu%data, &
            hv%data, un%data, vn%data, ua%data, ht%data, ssha_u%data, &
            va%data, ssha_v%data, istp, nx, ny)

    END SUBROUTINE invoke_time_step
end module time_step_mod
