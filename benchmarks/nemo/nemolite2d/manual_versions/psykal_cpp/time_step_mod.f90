module time_step_mod
    use field_mod
    use kind_params_mod
    implicit none
    public invoke_time_step

contains

    ! invoke_time_step needs to fetch all data as raw c arrays, the necessary
    ! global variables, and the additional values needed during the PSy-layer
    ! implementation, and pass them all to the C++ function using proper
    ! iso_c_bindings.
    SUBROUTINE invoke_time_step(ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, \
            vn, ua, ht, ssha_u, va, ssha_v, istp)
        USE field_mod
        USE kind_params_mod
        USE iso_c_binding
        USE physical_params_mod, ONLY: omega, d2r, g
        USE model_mod, ONLY: rdt, cbfr, visc

        implicit none

        TYPE(r2d_field), intent(inout), target :: ssha_t, sshn_t, sshn_u, \
            sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
        INTEGER, intent(inout) :: istp
        INTEGER:: istop, jstop


        LOGICAL, save :: first_time=.true.

        ! Look-up loop bounds
        istop = ssha_t%grid%subdomain%internal%xstop
        jstop = ssha_t%grid%subdomain%internal%ystop

        ! Any steps that need to be done just the first time
        IF (first_time) THEN
            first_time = .false.
        END IF

        call c_invoke_time_step()

    END SUBROUTINE invoke_time_step
end module time_step_mod
