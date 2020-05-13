!BEGINSOURCE 
  PROGRAM gocean2d
    USE dl_timer
    USE grid_mod
    USE field_mod
    USE initialisation_mod, ONLY: initialisation
    USE model_mod
    USE boundary_conditions_mod
    USE gocean2d_io_mod, ONLY: model_write
    USE gocean_mod, ONLY: model_write_log

    !> GOcean2d is a Horizontal 2D hydrodynamic ocean model initially developed
    !! by Hedong Liu, UK National Oceanography Centre (NOC), which:
    !!   1) uses structured grid
    !!   2) uses direct data addressing structures

    IMPLICIT NONE

    !> The grid on which our fields are defined
    TYPE(grid_type), target :: model_grid
    !> Current ('now') sea-surface height at different grid points
    TYPE(r2d_field) sshn_u_fld, sshn_v_fld, sshn_t_fld
    !> 'After' sea-surface height at different grid points
    TYPE(r2d_field) ssha_u_fld, ssha_v_fld, ssha_t_fld
    !> Distance from sea-bed to mean sea level at the different grid points.
    !! This is not time varying.
    TYPE(r2d_field) ht_fld, hu_fld, hv_fld
    !> Current ('now') velocity components
    TYPE(r2d_field) un_fld, vn_fld
    !> 'After' velocity components
    TYPE(r2d_field) ua_fld, va_fld

    ! time stepping index
    INTEGER istp
    INTEGER itimer0
    INTEGER(KIND=i_def64) nrepeat

    ! Create the model grid. We use a NE offset (i.e. the U, V and F
    ! points immediately to the North and East of a T point all have the
    ! same i,j index).  This is the same offset scheme as used by NEMO.
    model_grid = grid_type(arakawa_c,                          (/bc_external,bc_external,bc_none/),                          offset_ne)
    !  BC_PERIODIC, BC_NON_PERIODIC ??

    !! read in model parameters and configure the model grid
    CALL model_init(model_grid)

    ! Create fields on this grid

    ! Sea-surface height now (current time step)
    sshn_u_fld = r2d_field(model_grid, u_points)
    sshn_v_fld = r2d_field(model_grid, v_points)
    sshn_t_fld = r2d_field(model_grid, t_points)

    ! Sea-surface height 'after' (next time step)
    ssha_u_fld = r2d_field(model_grid, u_points)
    ssha_v_fld = r2d_field(model_grid, v_points)
    ssha_t_fld = r2d_field(model_grid, t_points)

    ! Distance from sea-bed to mean sea level
    hu_fld = r2d_field(model_grid, u_points)
    hv_fld = r2d_field(model_grid, v_points)
    ht_fld = r2d_field(model_grid, t_points)

    ! Velocity components now (current time step)
    un_fld = r2d_field(model_grid, u_points)
    vn_fld = r2d_field(model_grid, v_points)

    ! Velocity components 'after' (next time step)
    ua_fld = r2d_field(model_grid, u_points)
    va_fld = r2d_field(model_grid, v_points)

    !! setup model initial conditions
    CALL initialisation(ht_fld, hu_fld, hv_fld, sshn_u_fld, sshn_v_fld, sshn_t_fld, un_fld, vn_fld)

    CALL model_write(model_grid, 0, ht_fld, sshn_t_fld, un_fld, vn_fld)

    ! Start timer for time-stepping section
    nrepeat = nitend - nit000 + 1
    CALL timer_start(itimer0, label='Time-stepping', num_repeats=nrepeat)

    !! time stepping
    DO istp = nit000, nitend, 1

      !call model_write_log("('istp == ',I6)",istp)

      CALL step(istp, ua_fld, va_fld, un_fld, vn_fld, sshn_t_fld, sshn_u_fld, sshn_v_fld, ssha_t_fld, ssha_u_fld, ssha_v_fld, hu_fld, hv_fld, ht_fld)

      CALL model_write(model_grid, istp, ht_fld, sshn_t_fld, un_fld, vn_fld)

    END DO 

    ! Stop the timer for the time-stepping section
    CALL timer_stop(itimer0)

    ! Compute and output some checksums for error checking
    CALL model_write_log("('ua checksum = ',E16.8)", field_checksum(ua_fld))
    CALL model_write_log("('va checksum = ',E16.8)", field_checksum(va_fld))

    !! finalise the model run
    CALL model_finalise

    CALL model_write_log("((A))", 'Simulation finished!!')

  END PROGRAM gocean2d

  !+++++++++++++++++++++++++++++++++++

  SUBROUTINE step(istp, ua, va, un, vn, sshn_t, sshn_u, sshn_v, ssha_t, ssha_u, ssha_v, hu, hv, ht)
    USE psy_gocean2d, ONLY: invoke_0
    USE kind_params_mod
    USE grid_mod
    USE field_mod
    USE model_mod, ONLY: rdt
    ! The model time-step
    USE gocean2d_io_mod, ONLY: model_write
    USE continuity_mod, ONLY: continuity
    USE momentum_mod, ONLY: momentum_u, momentum_v
    USE boundary_conditions_mod, ONLY: bc_ssh, bc_solid_u, bc_solid_v, bc_flather_u, bc_flather_v
    USE time_update_mod, ONLY: next_sshu, next_sshv
    USE infrastructure_mod, ONLY: copy
    IMPLICIT NONE
    !> The current time step
    INTEGER, intent(inout) :: istp
    TYPE(r2d_field), intent(inout) :: un, vn, sshn_t, sshn_u, sshn_v
    TYPE(r2d_field), intent(inout) :: ua, va, ssha_t, ssha_u, ssha_v
    TYPE(r2d_field), intent(inout) :: hu, hv, ht

    CALL invoke_0(ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, rdt, ua, ht, ssha_u, va, ssha_v, istp)

  END SUBROUTINE step
