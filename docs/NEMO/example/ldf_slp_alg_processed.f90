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

    call invoke_0(zgru,prd,zgrv)
    if (ln_zps) then
      call invoke_1(zgru,gru,zgrv,grv)
    end if
    call invoke_2(zdzr,prd,pn2,uslpml,zgru,vslpml,zgrv,wslpiml,wslpjml,zwz,uslpml, &
                  zww,uslp,vslp,zwz,wslpi,wslpj)
  end subroutine ldf_slp

end module ldfslp_mod
