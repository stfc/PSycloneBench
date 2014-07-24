program gocean2d
  use grid_mod
  use field_mod
  use model_mod
  use boundary_conditions_mod
  !!! A Horizontal 2D hydrodynamic ocean model which
  !!   1) using structured grid
  !!   2) using direct data addressig structures

  implicit none

  !> The grid on which our fields are defined
  type(grid_type), target :: model_grid
  type(r2d_field_type) :: blah

  ! time stepping index
  integer  :: istp   

  ! Create the model grid. We use a NE staggering (i.e. the U, V and F 
  ! points to the North and East of a T point all have the same i,j index). 
  ! This is the same staggering as used by NEMO.
  model_grid = grid_type(ARAKAWA_C, STAGGER_NE)

  !! read in model parameters and read in or setup model grid 
  CALL model_init(model_grid)

  call boundary_conditions_init(model_grid)

  !! setup model initial condition
  CALL initialisation(model_grid)

  !! time stepping 
  DO istp = nit000, nitend, 1
     print*, 'istp == ', istp
     CALL step(model_grid, istp)
  END DO

  !! finalise the model run
  call model_finalise()
  
  WRITE(*,*) 'Simulation finished!!'
end program gocean2d

!+++++++++++++++++++++++++++++++++++

subroutine initialisation(grid)
  use kind_params_mod
  use model_mod
  use grid_mod
  use boundary_conditions_mod, only: bc
  implicit none
  type(grid_type), intent(in) :: grid
  ! define (or read in) initil ssh and velocity fields
  !         ! split this part into ssh, sshu, sshv, u, v kernels 
  integer :: ji, jj
  integer :: itmp1, itmp2
  real(wp) :: rtmp1

  DO ji=1,jpi
     DO jj =1, jpj
        sshn(ji,jj) = 0.0_wp
     END DO
  END DO

  DO ji=0,jpi
     DO jj =1, jpj
        itmp1 = min(ji+1,jpi)
        itmp2 = max(ji,1)
        rtmp1 = grid%e12t(itmp1,jj) * sshn(itmp1,jj) + grid%e12t(itmp2,jj) * sshn(itmp2,jj)
        sshn_u(ji,jj) = 0.5_wp * rtmp1 / grid%e12u(ji,jj)
     END DO
  END DO

  DO ji=1,jpi
     DO jj =0, jpj
        itmp1 = min(jj+1,jpj)
        itmp2 = max(jj,1)
        rtmp1 = grid%e12t(ji,itmp1) * sshn(ji,itmp1) + grid%e12t(ji,itmp2) * sshn(ji,itmp2)
        sshn_v(ji,jj) = 0.5_wp * rtmp1 / grid%e12v(ji,jj)
     END DO
  END DO

  DO jj =1, jpj
     DO ji=0,jpi
        un(ji,jj) = 0._wp
     END DO
  END DO

  DO jj =0, jpj
     DO ji=1,jpi
        vn(ji,jj) = 0._wp
     END DO
  END DO
  
  CALL bc(0._wp, grid, sshn_u, sshn_v, ssha, ua, va, hu, hv)
  
end subroutine initialisation

!+++++++++++++++++++++++++++++++++++

subroutine step(grid, istp)
  use kind_params_mod
  use grid_mod
  use model_mod, only: jpi, jpj, ua, va, un, vn, &
                       ssha, sshn, sshn_u, sshn_v, ssha_u, ssha_v, &
                       hu, hv, ht, rdt, irecord
  use momentum_mod, only: momentum
  use continuity_mod, only: invoke_continuity
  use time_update_mod, only: next
  use boundary_conditions_mod, only: bc
  use gocean2d_io_mod, only: model_write
  implicit none
  type(grid_type), intent(in) :: grid
  integer,         intent(in) :: istp
  real(wp) :: rtime

  rtime = REAL(istp, wp) * rdt

  CALL invoke_continuity(grid, ssha, sshn, sshn_u, sshn_v, &
                         hu, hv, un, vn)

  CALL momentum(grid, jpi, jpj, ua, va, un, vn, &
                sshn, sshn_u, sshn_v, ssha_u, ssha_v, &
                hu, hv, ht)

  ! open and solid boundary condition
  CALL bc(rtime, grid, sshn_u, sshn_v, ssha, ua, va, hu, hv)

  CALL next(grid)

  IF(MOD(istp, irecord) == 0)  CALL model_write(grid, istp, ht, sshn, un, vn)

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

