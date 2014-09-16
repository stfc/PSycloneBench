module initialisation_mod
  use field_mod
  implicit none

contains

  subroutine initialisation(ht_fld, hu_fld, hv_fld, &
                            sshn_t_fld, sshn_u_fld, sshn_v_fld, &
                            un_fld, vn_fld)
    use kind_params_mod
    use model_mod
    use grid_mod
    !use boundary_conditions_mod, only: bc
    implicit none
    type(r2d_field), intent(inout) :: ht_fld, hu_fld, hv_fld
    type(r2d_field), intent(inout) :: sshn_t_fld, sshn_u_fld, sshn_v_fld
    type(r2d_field), intent(inout) :: un_fld, vn_fld

    ! define (or read in) initial ssh and velocity fields
    !         ! split this part into ssh, sshu, sshv, u, v kernels 
    integer :: ji, jj
    integer :: itmp1, itmp2
    real(wp) :: rtmp1

    ! Depth at various grid points. This is constant here and is set to 
    ! the value read from the namelist file in the model_mod module.
    ht_fld%data(:,:) = dep_const 
    hu_fld%data(:,:) = dep_const 
    hv_fld%data(:,:) = dep_const 

    ! Sea surface height at T points ?
    sshn_t_fld%data(:,:) = 0.0_wp

    ! Sea-surface height at u points
    ! In original code this loop is over 0:jpi,1:jpj
    DO ji=1,sshn_u_fld%internal%xstop
       DO jj =1, sshn_u_fld%internal%ystop
          itmp1 = min(ji+1,sshn_u_fld%internal%nx)
          itmp2 = max(ji  ,1)
          rtmp1 = sshn_u_fld%grid%e12t(itmp1,jj) * sshn_t_fld%data(itmp1,jj) + &
                  sshn_u_fld%grid%e12t(itmp2,jj) * sshn_t_fld%data(itmp2,jj)
          sshn_u_fld%data(ji,jj) = 0.5_wp * rtmp1 / sshn_u_fld%grid%e12u(ji,jj)
       END DO
    END DO

    ! Sea-surface height at v points
    ! In original code this loop is over 1:jpi,0:jpj
    DO jj =1, sshn_v_fld%internal%ystop
      DO ji=1, sshn_v_fld%internal%xstop
        itmp1 = min(jj+1,sshn_v_fld%internal%ny)
        itmp2 = max(jj  ,1)
        rtmp1 = sshn_v_fld%grid%e12t(ji,itmp1) * sshn_t_fld%data(ji,itmp1) + &
                sshn_v_fld%grid%e12t(ji,itmp2) * sshn_t_fld%data(ji,itmp2)
        sshn_v_fld%data(ji,jj) = 0.5_wp * rtmp1 / sshn_v_fld%grid%e12v(ji,jj)
      END DO
    END DO

    ! Horizontal component of velocity (at U pts)
    un_fld%data(:,:) = 0._wp

    ! Vertical component of velocity (at V pts)
    vn_fld%data(:,:) = 0._wp
  
    ! Original bc() routine acts on ssha, ua and va. None of which
    ! this initialisation routine sets. Therefore, don't think it's
    ! needed here.
    !CALL bc(0._wp, sshn_u_fld, sshn_v_fld, ssha_t_fld, &
    !        ua_fld, va_fld, hu_fld, hv_fld)
  
  end subroutine initialisation

end module initialisation_mod
