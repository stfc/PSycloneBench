module ldfslp_mod

  IMPLICIT NONE

contains

  subroutine ldf_slp(uslp,wslpi,vslp,wslpj,...)
    !! ...
    !! ** Purpose :   Compute the slopes of neutral surface
    !! ...
    type(r3d_field), intent(in)  :: ...
    type(r3d_field), intent(out) :: uslp,wslpi,vslp,wslpj

    ! local fields
    type(r3d_field) :: zgru,zgrv,zwz,zww,zdzr

    zgru = r3d_field(model_grid, U_POINTS)
    zgrv = r3d_field(model_grid, V_POINTS)
    zwz  = r3d_field(model_grid, U_POINTS)
    zww  = r3d_field(model_grid, V_POINTS)
    zdzr = r3d_field(model_grid, T_POINTS)

    call invoke( horizontal_density_gradient_u(zgru,prd),                   &
                 horizontal_density_gradient_v(zgrv,prd) )
    if (ln_zps) then
      call invoke( partial_step_correction_u(zgru,gru),                     &
                   partial_step_correction_v(zgrv,grv) )
    end if
    call invoke( vertical_density_gradient_at_t(zdzr,prd,pn2),              &
                 slope_mixing_layer_u(uslpml,prd,pn2,zgru,zdzr),            &
                 slope_mixing_layer_v(vslpml,prd,pn2,zgrv,zdzr),            &
                 slope_mixing_layer_wi(wslpiml,prd,pn2,zgru),               &
                 slope_mixing_layer_wj(wslpjml,prd,pn2,zgrv),               &
                 vertical_density_gradient_at_u(zwz,zgru,zdzr,uslpml),      &
                 vertical_density_gradient_at_v(zww,zgrv,zdzr,vslpml),      &
                 horizontal_shapiro_filter(uslp,zwz),                       &
                 horizontal_shapiro_filter(vslp,zww),                       &
                 vertical_density_gradient_at_wi(zwz,zgru,prd,pn2,wslpiml), &
                 vertical_density_gradient_at_wj(zww,zgrv,prd,pn2,wslpjml), &
                 horizontal_shapiro_filter(wslpi,zwz),                      &
                 horizontal_shapiro_filter(wslpj,zww) )
  end subroutine ldf_slp

end module ldfslp_mod
