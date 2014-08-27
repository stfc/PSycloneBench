program gocean2d
  use grid_mod
  use field_mod
  use initialisation_mod, only: initialisation
  use model_mod
  use boundary_conditions_mod
  !!! A Horizontal 2D hydrodynamic ocean model which
  !!   1) using structured grid
  !!   2) using direct data addressing structures

  implicit none

  !> The grid on which our fields are defined
  type(grid_type), target :: model_grid
  !> Current ('now') sea-surface height at different grid points
  type(r2d_field_type) :: sshn_u_fld, sshn_v_fld, sshn_t_fld
  !> 'After' sea-surface height at different grid points
  type(r2d_field_type) :: ssha_u_fld, ssha_v_fld, ssha_t_fld
  !> Depth at the different grid points
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
  model_grid = grid_type(ARAKAWA_C, STAGGER_NE)

  !! read in model parameters and read in or setup model grid 
  CALL model_init(model_grid)

  call boundary_conditions_init(model_grid)

  ! Create fields on this grid
  sshn_u_fld = r2d_field_type(model_grid, &
                              U_POINTS,   &
                              (/BC_EXTERNAL,BC_EXTERNAL/))
  sshn_v_fld = r2d_field_type(model_grid, &
                              V_POINTS,   &
                              (/BC_EXTERNAL,BC_EXTERNAL/))
  sshn_t_fld = r2d_field_type(model_grid, &
                              T_POINTS,   &
                              (/BC_EXTERNAL,BC_EXTERNAL/))

  ssha_u_fld = r2d_field_type(model_grid, &
                              U_POINTS,   &
                              (/BC_EXTERNAL,BC_EXTERNAL/))
  ssha_v_fld = r2d_field_type(model_grid, &
                              V_POINTS,   &
                              (/BC_EXTERNAL,BC_EXTERNAL/))
  ssha_t_fld = r2d_field_type(model_grid, &
                              T_POINTS,   &
                              (/BC_EXTERNAL,BC_EXTERNAL/))

  hu_fld = r2d_field_type(model_grid, &
                          U_POINTS,   &
                          (/BC_EXTERNAL,BC_EXTERNAL/))
  hv_fld = r2d_field_type(model_grid, &
                          V_POINTS,   &
                          (/BC_EXTERNAL,BC_EXTERNAL/))
  ht_fld = r2d_field_type(model_grid, &
                          T_POINTS,   &
                          (/BC_EXTERNAL,BC_EXTERNAL/))

  un_fld = r2d_field_type(model_grid, &
                          U_POINTS,   &
                          (/BC_EXTERNAL,BC_EXTERNAL/))
  vn_fld = r2d_field_type(model_grid, &
                          V_POINTS,   &
                          (/BC_EXTERNAL,BC_EXTERNAL/))

  ua_fld = r2d_field_type(model_grid, &
                          U_POINTS,   &
                          (/BC_EXTERNAL,BC_EXTERNAL/))
  va_fld = r2d_field_type(model_grid, &
                          V_POINTS,   &
                          (/BC_EXTERNAL,BC_EXTERNAL/))

  !! setup model initial conditions
  call initialisation(ht_fld, hu_fld, hv_fld, &
                      sshn_u_fld, sshn_v_fld, sshn_t_fld, &
                      un_fld, vn_fld)

  !! time stepping 
  DO istp = nit000, nitend, 1
     print*, 'istp == ', istp

     CALL step(model_grid, istp, &
               ua_fld, va_fld, un_fld, vn_fld, &
               sshn_t_fld, sshn_u_fld, sshn_v_fld, &
               ssha_t_fld, ssha_u_fld, ssha_v_fld, &
               hu_fld, hv_fld, ht_fld)
  END DO

  !! finalise the model run
  call model_finalise()
  
  WRITE(*,*) 'Simulation finished!!'
end program gocean2d

!+++++++++++++++++++++++++++++++++++

subroutine step(grid, istp, &
                ua, va, un, vn, &
                sshn, sshn_u, sshn_v, ssha, ssha_u, ssha_v, &
                hu, hv, ht)
  use kind_params_mod
  use grid_mod
  use field_mod
  use model_mod, only: rdt, irecord
  use momentum_mod, only: invoke_momentum_u, invoke_momentum_v
  use continuity_mod, only: invoke_continuity
  use time_update_mod, only: next
  use boundary_conditions_mod
  use gocean2d_io_mod, only: model_write
  implicit none
  type(grid_type), intent(in) :: grid
  integer,         intent(in) :: istp
  type(r2d_field_type), intent(in) :: un, vn, sshn, sshn_u, sshn_v
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
  CALL bc_ssh(rtime, ssha)
  CALL bc_u_solid(ua)
  CALL bc_v_solid(va)
  CALL bc_u_flather(ua, hu, sshn_u)
  CALL bc_v_flather(va, hv, sshn_v)

  CALL next(grid)

  IF(MOD(istp, irecord) == 0)  CALL model_write(grid, istp, ht%data, &
                                                sshn%data, un%data, vn%data)

end subroutine step

!+++++++++++++++++++++++++++++++++++

SUBROUTINE output(grid, istp, ht, sshn, un, vn)
  use kind_params_mod
  use grid_mod
  use model_mod, only: jpi, jpj
  implicit none
  type(grid_type), intent(in) :: grid
  integer, intent(in) :: istp
  real(wp), dimension(:,:), intent(in) :: ht, sshn, un, vn
  ! Locals
  integer :: ji, jj
  real(wp) :: rtmp1, rtmp2

  ! output model results
  CHARACTER(len=5) :: fname
  WRITE(fname, '(I5.5)') istp
  !OPEN(21, file='go2d_'//fname//'.dat', STATUS='UNKNOWN')
  OPEN(21, file='go2d_'//fname//'.csv', STATUS='UNKNOWN')
  REWIND(21)

  DO jj = 1, jpj
     DO ji = 1, jpi
        rtmp1 = 0.5_wp * (un(ji-1,jj) + un(ji,jj))
        rtmp2 = 0.5_wp * (vn(ji,jj-1) + vn(ji,jj))

        ! write "x-coord, y-coord, depth, ssh, u-velocity, v-velocity" to ASCII files

        !WRITE(1,'(2f20.3, 2f15.4, 2e18.3)')  &            
        WRITE(1,'(f20.3,'','',f20.3,'','',f15.4,'','',f15.4,'','',f18.3,'','',f18.3)') &
             & grid%xt(ji,jj), grid%yt(ji,jj), ht(ji,jj), sshn(ji,jj),rtmp1, rtmp2 
     END DO
  END DO
          
  CLOSE(21)

END SUBROUTINE output

