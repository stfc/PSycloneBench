module initialisation_mod
  implicit none

contains

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

    sshn(:,:) = 0.0_wp

    DO ji=1,grid%nx
       DO jj =1, grid%ny
          itmp1 = min(ji+1,grid%nx)
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

end module initialisation_mod
