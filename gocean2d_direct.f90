program gocean2d
  use timing_mod
  use grid_mod
  use field_mod
  use initialisation_mod, only: initialisation
  use model_mod
  use boundary_conditions_mod
  use gocean2d_io_mod, only: model_write
  use gocean_mod,      only: model_write_log

  !> A Horizontal 2D hydrodynamic ocean model which
  !!   1) using structured grid
  !!   2) using direct data addressing structures

  implicit none

  !> The grid on which our fields are defined
  type(grid_type), target :: model_grid
  !> Current ('now') sea-surface height at different grid points
  type(r2d_field) :: sshn_u_fld, sshn_v_fld, sshn_t_fld
  !> 'After' sea-surface height at different grid points
  type(r2d_field) :: ssha_u_fld, ssha_v_fld, ssha_t_fld
  !> Distance from sea-bed to mean sea level at the different grid points.
  !! This is not time varying.
  type(r2d_field) :: ht_fld, hu_fld, hv_fld
  !> Current ('now') velocity components
  type(r2d_field) :: un_fld, vn_fld
  !> 'After' velocity components
  type(r2d_field) :: ua_fld, va_fld

  ! time stepping index
  integer :: istp   
  integer :: itimer0

  ! Create the model grid. We use a NE offset (i.e. the U, V and F
  ! points immediately to the North and East of a T point all have the
  ! same i,j index).  This is the same offset scheme as used by NEMO.
  model_grid = grid_type(ARAKAWA_C, &
  !  BC_PERIODIC, BC_NON_PERIODIC ??
                         (/BC_EXTERNAL,BC_EXTERNAL,BC_NONE/), &
                         OFFSET_NE)

  !! read in model parameters and configure the model grid 
  CALL model_init(model_grid)

  ! Create fields on this grid

  ! Sea-surface height now (current time step)
  sshn_u_fld = r2d_field(model_grid, U_POINTS)
  sshn_v_fld = r2d_field(model_grid, V_POINTS)
  sshn_t_fld = r2d_field(model_grid, T_POINTS)

  ! Sea-surface height 'after' (next time step)
  ssha_u_fld = r2d_field(model_grid, U_POINTS)
  ssha_v_fld = r2d_field(model_grid, V_POINTS)
  ssha_t_fld = r2d_field(model_grid, T_POINTS)

  ! Distance from sea-bed to mean sea level
  hu_fld = r2d_field(model_grid, U_POINTS)
  hv_fld = r2d_field(model_grid, V_POINTS)
  ht_fld = r2d_field(model_grid, T_POINTS)

  ! Velocity components now (current time step)
  un_fld = r2d_field(model_grid, U_POINTS)
  vn_fld = r2d_field(model_grid, V_POINTS)

  ! Velocity components 'after' (next time step)
  ua_fld = r2d_field(model_grid, U_POINTS)
  va_fld = r2d_field(model_grid, V_POINTS)

  !! setup model initial conditions
  call initialisation(ht_fld, hu_fld, hv_fld, &
                      sshn_u_fld, sshn_v_fld, sshn_t_fld, &
                      un_fld, vn_fld)

  call model_write(model_grid, 0, ht_fld, sshn_t_fld, un_fld, vn_fld)

  ! Start timer for time-stepping section
  CALL timer_start('Time-stepping',itimer0)

  !! time stepping 
  do istp = nit000, nitend, 1

     !call model_write_log("('istp == ',I6)",istp)

     call step(model_grid, istp, &
               ua_fld, va_fld, un_fld, vn_fld, &
               sshn_t_fld, sshn_u_fld, sshn_v_fld, &
               ssha_t_fld, ssha_u_fld, ssha_v_fld, &
               hu_fld, hv_fld, ht_fld)
  end do

  ! Stop the timer for the time-stepping section
  call timer_stop(itimer0)

  ! Compute and output some checksums for error checking
  call model_write_log("('ua checksum = ',E16.8)", &
                       field_checksum(ua_fld))
  call model_write_log("('va checksum = ',E16.8)", &
                       field_checksum(va_fld))

  !! finalise the model run
  call model_finalise()
  
  call model_write_log("((A))", 'Simulation finished!!')

end program gocean2d

!+++++++++++++++++++++++++++++++++++

subroutine step(grid, istp, &
                ua, va, un, vn, &
                sshn, sshn_u, sshn_v, ssha, ssha_u, ssha_v, &
                hu, hv, ht)
  use kind_params_mod
  use grid_mod
  use field_mod
  use time_step_mod, only: invoke_time_step
  use gocean2d_io_mod, only: model_write
  implicit none
  type(grid_type), intent(in) :: grid
  !> The current time step
  integer,         intent(in) :: istp
  type(r2d_field), intent(inout) :: un, vn, sshn, sshn_u, sshn_v
  type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
  type(r2d_field), intent(in) :: hu, hv, ht

  call invoke_time_step(istp, ssha, ssha_u, ssha_v, &
                        sshn, sshn_u, sshn_v, &
                        hu, hv, ht, ua, va, un, vn)

!  call invoke(                                               &
!              continuity(istp, ssha, sshn_t, sshn_u, sshn_v, &
!                         hu, hv, un, vn),                    &
!              momentum_u(ua, un, vn,                         &
!                         ssha_u, sshn_t, sshn_u, sshn_v),    &
!              momentum_v(va, un, vn, hu, hv, ht,             &
!                         ssha_v, sshn_t, sshn_u, sshn_v),    &
!              bc_ssh(istp, ssha),                            &
!              bc_solid_u(ua),                                &
!              bc_solid_v(va),                                &
!              bc_flather_u(ua, hu, sshn_u),                  &
!              bc_flather_v(va, hv, sshn_v),                  &
!              copy_field(ua, un),                            &
!              copy_field(va, vn),                            &
!              copy_field(ssha, sshn_t),                      &
!              next_sshu(sshn_u, sshn_t),                     &
!              next_sshv(sshn_v, sshn_t)                      &
!             )


  call model_write(grid, istp, ht, sshn, un, vn)

end subroutine step

