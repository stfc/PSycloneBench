program gocean2d
  use dl_timer, only: timer_start, timer_stop, timer_init, timer_report, i_def64
  use grid_mod
  use field_mod
  use initialisation_mod, only: initialisation
  use model_mod
  use gocean2d_io_mod, only: model_write
  use gocean_mod,      only: model_write_log, gocean_initialise, &
                             gocean_finalise

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
  integer     :: istp  
  real(go_wp) :: rstp 
  integer     :: itimer0

  ! Scratch space for logging messages
  character(len=160) :: log_str

  ! Initialise GOcean infrastructure
  call gocean_initialise()

  ! Create the model grid. We use a NE offset (i.e. the U, V and F
  ! points immediately to the North and East of a T point all have the
  ! same i,j index).  This is the same offset scheme as used by NEMO.
  model_grid = grid_type(GO_ARAKAWA_C, &
                         (/GO_BC_EXTERNAL,GO_BC_EXTERNAL,GO_BC_NONE/), &
                         GO_OFFSET_NE)

  !! read in model parameters and configure the model grid 
  CALL model_init(model_grid)

  ! Create fields on this grid

  ! Sea-surface height now (current time step)
  sshn_u_fld = r2d_field(model_grid, GO_U_POINTS)
  sshn_v_fld = r2d_field(model_grid, GO_V_POINTS)
  sshn_t_fld = r2d_field(model_grid, GO_T_POINTS)

  ! Sea-surface height 'after' (next time step)
  ssha_u_fld = r2d_field(model_grid, GO_U_POINTS)
  ssha_v_fld = r2d_field(model_grid, GO_V_POINTS)
  ssha_t_fld = r2d_field(model_grid, GO_T_POINTS)

  ! Distance from sea-bed to mean sea level
  hu_fld = r2d_field(model_grid, GO_U_POINTS)
  hv_fld = r2d_field(model_grid, GO_V_POINTS)
  ht_fld = r2d_field(model_grid, GO_T_POINTS)

  ! Velocity components now (current time step)
  un_fld = r2d_field(model_grid, GO_U_POINTS)
  vn_fld = r2d_field(model_grid, GO_V_POINTS)

  ! Velocity components 'after' (next time step)
  ua_fld = r2d_field(model_grid, GO_U_POINTS)
  va_fld = r2d_field(model_grid, GO_V_POINTS)

  !! setup model initial conditions
  call initialisation(ht_fld, hu_fld, hv_fld, &
                      sshn_u_fld, sshn_v_fld, sshn_t_fld, &
                      un_fld, vn_fld)

  call model_write(model_grid, 0, ht_fld, sshn_t_fld, un_fld, vn_fld)

  write(log_str, "('Simulation domain = (',I4,':',I4,',',I4,':',I4,')')") &
                       model_grid%subdomain%global%xstart, &
                       model_grid%subdomain%global%xstop,  &
                       model_grid%subdomain%global%ystart, &
                       model_grid%subdomain%global%ystop
  call model_write_log("((A))", TRIM(log_str))

  ! Start timer for time-stepping section
  CALL timer_start(itimer0, label='Time-stepping', &
                   num_repeats=INT(nitend-nit000+1,kind=i_def64) )

  !! time stepping 
  do istp = nit000, nitend, 1

     call step(istp,                               &
               ua_fld, va_fld, un_fld, vn_fld,     &
               sshn_t_fld, sshn_u_fld, sshn_v_fld, &
               ssha_t_fld, ssha_u_fld, ssha_v_fld, &
               hu_fld, hv_fld, ht_fld)

     call model_write(model_grid, istp,                &
                      ht_fld, sshn_t_fld, un_fld, vn_fld)

  end do

  ! Stop the timer for the time-stepping section
  call timer_stop(itimer0)

  ! Compute and output some checksums for error checking
  call model_write_log("('ua checksum = ', E16.8)", &
                       field_checksum(ua_fld))
  call model_write_log("('va checksum = ', E16.8)", &
                       field_checksum(va_fld))

  !! finalise the model run
  call model_finalise()

  call model_write_log("((A))", 'Simulation finished!!')

  call gocean_finalise()

end program gocean2d

!+++++++++++++++++++++++++++++++++++

subroutine step(istp,           &
                ua, va, un, vn, &
                sshn, sshn_u, sshn_v, ssha, ssha_u, ssha_v, &
                hu, hv, ht)
  use kind_params_mod
  use grid_mod
  use field_mod
  use time_step_mod, only: invoke_time_step
  use gocean2d_io_mod, only: model_write
  implicit none
  !> The current time step
  integer,         intent(in) :: istp
  type(r2d_field), intent(inout) :: un, vn, sshn, sshn_u, sshn_v
  type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
  type(r2d_field), intent(inout) :: hu, hv, ht

  call invoke_time_step(istp, ssha, ssha_u, ssha_v, &
                        sshn, sshn_u, sshn_v, &
                        hu, hv, ht, ua, va, un, vn)

end subroutine step

