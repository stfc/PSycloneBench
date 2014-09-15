program gocean2d
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
  type(r2d_field_type) :: sshn_u_fld, sshn_v_fld, sshn_t_fld
  !> 'After' sea-surface height at different grid points
  type(r2d_field_type) :: ssha_u_fld, ssha_v_fld, ssha_t_fld
  !> Distance from sea-bed to mean sea level at the different grid points.
  !! This is not time varying.
  type(r2d_field_type) :: ht_fld, hu_fld, hv_fld
  !> Current ('now') velocity components
  type(r2d_field_type) :: un_fld, vn_fld
  !> 'After' velocity components
  type(r2d_field_type) :: ua_fld, va_fld

  ! time stepping index
  integer  :: istp   

  ! Create the model grid. We use a NE staggering (i.e. the U, V and F 
  ! points to the North and East of a T point all have the same i,j index). 
  ! This is the same staggering as used by NEMO.
  model_grid = grid_type(ARAKAWA_C, STAGGER_NE, &
                         (/BC_EXTERNAL,BC_EXTERNAL,BC_NONE/))

  !! read in model parameters and configure the model grid 
  CALL model_init(model_grid)

  ! Create fields on this grid

  ! Sea-surface height now (current time step)
  sshn_u_fld = r2d_field_type(model_grid, U_POINTS)
  sshn_v_fld = r2d_field_type(model_grid, V_POINTS)
  sshn_t_fld = r2d_field_type(model_grid, T_POINTS)

  ! Sea-surface height 'after' (next time step)
  ssha_u_fld = r2d_field_type(model_grid, U_POINTS)
  ssha_v_fld = r2d_field_type(model_grid, V_POINTS)
  ssha_t_fld = r2d_field_type(model_grid, T_POINTS)

  ! Distance from sea-bed to mean sea level
  hu_fld = r2d_field_type(model_grid, U_POINTS)
  hv_fld = r2d_field_type(model_grid, V_POINTS)
  ht_fld = r2d_field_type(model_grid, T_POINTS)

  ! Velocity components now (current time step)
  un_fld = r2d_field_type(model_grid, U_POINTS)
  vn_fld = r2d_field_type(model_grid, V_POINTS)

  ! Velocity components 'after' (next time step)
  ua_fld = r2d_field_type(model_grid, U_POINTS)
  va_fld = r2d_field_type(model_grid, V_POINTS)

  !! setup model initial conditions
  call initialisation(ht_fld, hu_fld, hv_fld, &
                      sshn_u_fld, sshn_v_fld, sshn_t_fld, &
                      un_fld, vn_fld)

  call model_write(model_grid, 0, ht_fld, sshn_t_fld, un_fld, vn_fld)

  !! time stepping 
  do istp = nit000, nitend, 1

     call model_write_log("('istp == ',I6)",istp)

     call step(model_grid, istp, &
               ua_fld, va_fld, un_fld, vn_fld, &
               sshn_t_fld, sshn_u_fld, sshn_v_fld, &
               ssha_t_fld, ssha_u_fld, ssha_v_fld, &
               hu_fld, hv_fld, ht_fld)
  end do

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
  use model_mod, only: rdt
  use momentum_mod, only: invoke_momentum_u, invoke_momentum_v
  use continuity_mod, only: invoke_continuity
  use time_update_mod, only: invoke_next_sshu, invoke_next_sshv
  use boundary_conditions_mod
  use gocean2d_io_mod, only: model_write
  implicit none
  type(grid_type), intent(in) :: grid
  integer,         intent(in) :: istp
  type(r2d_field_type), intent(inout) :: un, vn, sshn, sshn_u, sshn_v
  type(r2d_field_type), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
  type(r2d_field_type), intent(in) :: hu, hv, ht
  ! Locals
  real(wp) :: rtime

  rtime = REAL(istp, wp) * rdt

  CALL invoke_continuity(ssha, sshn, sshn_u, sshn_v, &
                         hu, hv, un, vn)

  CALL invoke_momentum_u(ua, un, vn, &
                         hu, hv, ht, &
                         ssha_u, sshn, sshn_u, sshn_v)

  CALL invoke_momentum_v(va, un, vn, &
                         hu, hv, ht, &
                         ssha_v, sshn, sshn_u, sshn_v)

  ! Apply open and solid boundary conditions
  CALL invoke_bc_ssh(rtime, ssha)
  CALL invoke_bc_solid_u(ua)
  CALL bc_v_solid(va)
  CALL bc_u_flather(ua, hu, sshn_u)
  CALL bc_v_flather(va, hv, sshn_v)

  ! Time update of fields
  call copy_field(ua, un)
  call copy_field(va, vn)
  call copy_field(ssha, sshn)
  call invoke_next_sshu(sshn_u, sshn)
  call invoke_next_sshv(sshn_v, sshn)

  call model_write(grid, istp, ht, sshn, un, vn)

end subroutine step
