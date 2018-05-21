MODULE psy_gocean2d
  IMPLICIT NONE
CONTAINS
  SUBROUTINE invoke_0(ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, rdt, ua, ht, ssha_u, va, ssha_v, istp)

    USE field_mod, only : r2d_field
    USE kind_params_mod, only : wp
    USE grid_mod, only : grid_type
    
    TYPE(r2d_field), intent(inout) :: ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
    REAL(KIND=wp), intent(inout) :: rdt
    INTEGER, intent(inout) :: istp

    INTEGER istop, jstop
    type(grid_type), pointer :: grid
    
    grid=>ssha_t%grid

    ! Look-up loop bounds
    istop = grid%simulation_domain%xstop
    jstop = grid%simulation_domain%ystop
    !
    call invoke_0_deref(ssha_t%data, sshn_t%data, sshn_u%data, sshn_v%data, hu%data, hv%data, un%data, vn%data, rdt, ua%data, ht%data, ssha_u%data, va%data, ssha_v%data, istp, istop, jstop, sshn_t%data, size(sshn_t%data,1), size(sshn_t%data,2), grid%area_t, size(grid%area_t,1), size(grid%area_t,2), grid%area_u, size(grid%area_u,1), size(grid%area_u,2), grid%area_v, size(grid%area_v,1), size(grid%area_v,2), grid%tmask, size(grid%tmask,1), size(grid%tmask,2), grid%dx_u, size(grid%dx_u,1), size(grid%dx_u,2), grid%dx_v, size(grid%dx_v,1), size(grid%dx_v,2), grid%dx_t, size(grid%dx_t,1), size(grid%dx_t,2), grid%dy_u, size(grid%dy_u,1), size(grid%dy_u,2), grid%dy_v, size(grid%dy_v,1), size(grid%dy_v,2), grid%dy_t, size(grid%dy_t,1), size(grid%dy_t,2), grid%gphiu, size(grid%gphiu,1), size(grid%gphiu,2), grid%gphiv, size(grid%gphiv,1), size(grid%gphiv,2))
    !
  END SUBROUTINE invoke_0
END MODULE psy_gocean2d

SUBROUTINE invoke_0_deref(ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, rdt, ua, ht, ssha_u, va, ssha_v, istp, istop, jstop, n, m, grid_area_t, gat_dim1, gat_dim2, grid_area_u, gau_dim1, gau_dim2, grid_area_v, gav_dim1, gav_dim2, grid_tmask, gtm_dim1, gtm_dim2, grid_dx_u, gxu_dim1, gxu_dim2, grid_dx_v, gxv_dim1, gxv_dim2, grid_dx_t, gxt_dim1, gxt_dim2, grid_dy_u, gyu_dim1, gyu_dim2, grid_dy_v, gyv_dim1, gyv_dim2, grid_dy_t, gyt_dim1, gyt_dim2, grid_gphiu, gphiu_dim1, gphiu_dim2, grid_gphiv, gphiv_dim1, gphiv_dim2)
  
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

  integer, intent(in) :: gat_dim1, gat_dim2, gau_dim1, gau_dim2, gav_dim1, gav_dim2, gtm_dim1, gtm_dim2
  integer, intent(in) :: gxu_dim1, gxu_dim2, gxv_dim1, gxv_dim2, gxt_dim1, gxt_dim2
  integer, intent(in) :: gyu_dim1, gyu_dim2, gyv_dim1, gyv_dim2, gyt_dim1, gyt_dim2
  integer, intent(in) :: gphiu_dim1, gphiu_dim2, gphiv_dim1, gphiv_dim2
  REAL(KIND=wp), intent(in) :: grid_area_t(gat_dim1, gat_dim2)
  REAL(KIND=wp), intent(in) :: grid_area_u(gau_dim1, gau_dim2)
  REAL(KIND=wp), intent(in) :: grid_area_v(gav_dim1, gav_dim2)
  integer, intent(in) :: grid_tmask(gtm_dim1, gtm_dim2)
  REAL(KIND=wp), intent(in) :: grid_dx_u(gxu_dim1, gxu_dim2)
  REAL(KIND=wp), intent(in) :: grid_dx_v(gxv_dim1, gxv_dim2)
  REAL(KIND=wp), intent(in) :: grid_dx_t(gxt_dim1, gxt_dim2)
  REAL(KIND=wp), intent(in) :: grid_dy_u(gyu_dim1, gyu_dim2)
  REAL(KIND=wp), intent(in) :: grid_dy_v(gyv_dim1, gyv_dim2)
  REAL(KIND=wp), intent(in) :: grid_dy_t(gyt_dim1, gyt_dim2)
  REAL(KIND=wp), intent(in) :: grid_gphiu(gphiu_dim1, gphiu_dim2)
  REAL(KIND=wp), intent(in) :: grid_gphiv(gphiv_dim1, gphiv_dim2)
  
  REAL(KIND=wp), dimension(n,m), intent(inout)  :: ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
  REAL(KIND=wp), intent(in) :: rdt
  INTEGER, intent(in) :: istp
  INTEGER, intent(in) :: istop, jstop
  INTEGER, intent(in) :: n, m
  !
  INTEGER :: j, i
  !
  DO j=2,jstop
     DO i=2,istop
        CALL continuity_code_wrap(i, j, ssha_t, size(ssha_t,1), size(ssha_t,2), sshn_t, size(sshn_t,1), size(sshn_t,2), sshn_u, size(sshn_u,1), size(sshn_u,2), sshn_v, size(sshn_v,1), size(sshn_v,2), hu, size(hu,1), size(hu,2), hv, size(hv,1), size(hv,2), un, size(un,1), size(un,2), vn, size(vn,1), size(vn,2), rdt, grid_area_t, size(grid_area_t,1), size(grid_area_t,2))
     END DO
  END DO
  DO j=2,jstop
     DO i=2,istop-1
        CALL momentum_u_code_wrap(i, j, ua, size(ua,1), size(ua,2), un, size(un,1), size(un,2), vn, size(vn,1), size(vn,2), hu, size(hu,1), size(hu,2), hv, size(hv,1), size(hv,2), ht, size(ht,1), size(ht,2), ssha_u, size(ssha_u,1), size(ssha_u,2), sshn_t, size(sshn_t,1), size(sshn_t,2), sshn_u, size(sshn_u,1), size(sshn_u,2), sshn_v, size(sshn_v,1), size(sshn_v,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2), grid_dx_u, size(grid_dx_u,1), size(grid_dx_u,2), grid_dx_v, size(grid_dx_v,1), size(grid_dx_v,2), grid_dx_t, size(grid_dx_t,1), size(grid_dx_t,2), grid_dy_u, size(grid_dy_u,1), size(grid_dy_u,2), grid_dy_t, size(grid_dy_t,1), size(grid_dy_t,2), grid_area_u, size(grid_area_u,1), size(grid_area_u,2), grid_gphiu, size(grid_gphiu,1), size(grid_gphiu,2))
     END DO
  END DO
      DO j=2,jstop-1
        DO i=2,istop
          CALL momentum_v_code_wrap(i, j, va, size(va,1), size(va,2), un, size(un,1), size(un,2), vn, size(vn,1), size(vn,2), hu, size(hu,1), size(hu,2), hv, size(hv,1), size(hv,2), ht, size(ht,1), size(ht,2), ssha_v, size(ssha_v,1), size(ssha_v,2), sshn_t, size(sshn_t,1), size(sshn_t,2), sshn_u, size(sshn_u,1), size(sshn_u,2), sshn_v, size(sshn_v,1), size(sshn_v,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2), grid_dx_v, size(grid_dx_v,1), size(grid_dx_v,2), grid_dx_t, size(grid_dx_t,1), size(grid_dx_t,2), grid_dy_u, size(grid_dy_u,1), size(grid_dy_u,2), grid_dy_v, size(grid_dy_v,1), size(grid_dy_v,2), grid_dy_t, size(grid_dy_t,1), size(grid_dy_t,2), grid_area_v, size(grid_area_v,1), size(grid_area_v,2), grid_gphiv, size(grid_gphiv,1), size(grid_gphiv,2))
        END DO 
      END DO 
      DO j=2,jstop
        DO i=2,istop
          CALL bc_ssh_code_wrap(i, j, istp, ssha_t, size(ssha_t,1), size(ssha_t,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop
          CALL bc_solid_u_code_wrap(i, j, ua, size(ua,1), size(ua,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
        END DO 
      END DO 
      DO j=1,jstop
        DO i=1,istop+1
          CALL bc_solid_v_code_wrap(i, j, va, size(va,1), size(va,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop
          CALL bc_flather_u_code_wrap(i, j, ua, size(ua,1), size(ua,2), hu, size(hu,1), size(hu,2), sshn_u, size(sshn_u,1), size(sshn_u,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
        END DO 
      END DO 
      DO j=1,jstop
        DO i=1,istop+1
          CALL bc_flather_v_code_wrap(i, j, va, size(va,1), size(va,2), hv, size(hv,1), size(hv,2), sshn_v, size(sshn_v,1), size(sshn_v,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop+1
          CALL field_copy_code_wrap(i, j, un, size(un,1), size(un,2), ua, size(ua,1), size(ua,2))
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop+1
          CALL field_copy_code_wrap(i, j, vn, size(vn,1), size(vn,2), va, size(va,1), size(va,2))
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop+1
          CALL field_copy_code_wrap(i, j, sshn_t, size(sshn_t,1), size(sshn_t,2), ssha_t, size(ssha_t,1), size(ssha_t,2))
        END DO 
      END DO 
      DO j=2,jstop
        DO i=2,istop-1
          CALL next_sshu_code_wrap(i, j, sshn_u, size(sshn_u,1), size(sshn_u,2), sshn_t, size(sshn_t,1), size(sshn_t,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2), grid_area_t, size(grid_area_t,1), size(grid_area_t,2), grid_area_u, size(grid_area_u,1), size(grid_area_u,2))
        END DO 
      END DO 
      DO j=2,jstop-1
        DO i=2,istop
          CALL next_sshv_code_wrap(i, j, sshn_v, size(sshn_v,1), size(sshn_v,2), sshn_t, size(sshn_t,1), size(sshn_t,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2), grid_area_t, size(grid_area_t,1), size(grid_area_t,2), grid_area_v, size(grid_area_v,1), size(grid_area_v,2))
        END DO 
      END DO 
END SUBROUTINE invoke_0_deref

subroutine continuity_code_wrap(i, j, ssha_t, at_dim1, at_dim2, sshn_t, nt_dim1, nt_dim2, sshn_u, nu_dim1, nu_dim2, sshn_v, nv_dim1, nv_dim2, hu, hu_dim1, hu_dim2, hv, hv_dim1, hv_dim2, un, un_dim1, un_dim2, vn, vn_dim1, vn_dim2, rdt, grid_area_t, gat_dim1, gat_dim2)
  USE continuity_mod, ONLY: continuity_code
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: i, j
  integer, intent(in)     ::at_dim1, at_dim2, nt_dim1, nt_dim2, nu_dim1, nu_dim2, nv_dim1, nv_dim2, hu_dim1, hu_dim2, hv_dim1, hv_dim2, un_dim1, un_dim2, vn_dim1, vn_dim2, gat_dim1, gat_dim2
  real(wp), intent(inout) :: ssha_t(at_dim1, at_dim2)
  real(wp), intent(in)    :: sshn_t(nt_dim1, nt_dim2)
  real(wp), intent(in)    :: sshn_u(nu_dim1, nu_dim2)
  real(wp), intent(in)    :: sshn_v(nv_dim1, nv_dim2)
  real(wp), intent(in)    :: hu(hu_dim1, hu_dim2)
  real(wp), intent(in)    :: hv(hv_dim1, hv_dim2)
  real(wp), intent(in)    :: un(un_dim1, un_dim2)
  real(wp), intent(in)    :: vn(vn_dim1, vn_dim2)
  real(wp), intent(in)    :: rdt
  real(wp), intent(in)    :: grid_area_t(gat_dim1, gat_dim2)
  !
  CALL continuity_code(i, j, ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, rdt, grid_area_t)
  !
end subroutine continuity_code_wrap

subroutine momentum_u_code_wrap(i, j, ua, ua_dim1, ua_dim2, un, un_dim1, un_dim2, vn, vn_dim1, vn_dim2, hu, hu_dim1, hu_dim2, hv, hv_dim1, hv_dim2, ht, ht_dim1, ht_dim2, ssha_u, ssha_u_dim1, ssha_u_dim2, sshn_t, sshn_t_dim1, sshn_t_dim2, sshn_u, sshn_u_dim1, sshn_u_dim2, sshn_v, sshn_v_dim1, sshn_v_dim2, grid_tmask, g_dim1, g_dim2, grid_dx_u, dxu_dim1, dxu_dim2, grid_dx_v, dxv_dim1, dxv_dim2, grid_dx_t, dxt_dim1, dxt_dim2, grid_dy_u, dyu_dim1, dyu_dim2, grid_dy_t, dyt_dim1, dyt_dim2, grid_area_u, gau_dim1, gau_dim2, grid_gphiu, gphiu_dim1, gphiu_dim2)
  USE momentum_mod, ONLY: momentum_u_code
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: i, j
  integer, intent(in)     :: ua_dim1, ua_dim2, un_dim1, un_dim2, vn_dim1, vn_dim2, hu_dim1, hu_dim2, hv_dim1, hv_dim2, ht_dim1, ht_dim2, ssha_u_dim1, ssha_u_dim2, sshn_t_dim1, sshn_t_dim2, sshn_u_dim1, sshn_u_dim2, sshn_v_dim1, sshn_v_dim2, g_dim1, g_dim2, dxu_dim1, dxu_dim2, dxv_dim1, dxv_dim2, dxt_dim1, dxt_dim2, dyu_dim1, dyu_dim2, dyt_dim1, dyt_dim2, gau_dim1, gau_dim2, gphiu_dim1, gphiu_dim2
  real(wp), intent(inout) :: ua(ua_dim1, ua_dim2)
  real(wp), intent(in)    :: un(un_dim1, un_dim2)
  real(wp), intent(in)    :: vn(vn_dim1, vn_dim2)
  real(wp), intent(in)    :: hu(hu_dim1, hu_dim2)
  real(wp), intent(in)    :: hv(hv_dim1, hv_dim2)
  real(wp), intent(in)    :: ht(ht_dim1, ht_dim2)
  real(wp), intent(in)    :: ssha_u(ssha_u_dim1, ssha_u_dim2)
  real(wp), intent(in)    :: sshn_t(sshn_t_dim1, sshn_t_dim2)
  real(wp), intent(in)    :: sshn_u(sshn_u_dim1, sshn_u_dim2)
  real(wp), intent(in)    :: sshn_v(sshn_v_dim1, sshn_v_dim2)
  integer,  intent(in)    :: grid_tmask(g_dim1, g_dim2)
  real(wp), intent(in)    :: grid_dx_u(dxu_dim1, dxu_dim2)
  real(wp), intent(in)    :: grid_dx_v(dxv_dim1, dxv_dim2)
  real(wp), intent(in)    :: grid_dx_t(dxt_dim1, dxt_dim2)
  real(wp), intent(in)    :: grid_dy_u(dyu_dim1, dyu_dim2)
  real(wp), intent(in)    :: grid_dy_t(dyt_dim1, dyt_dim2)
  real(wp), intent(in)    :: grid_area_u(gau_dim1, gau_dim2)
  real(wp), intent(in)    :: grid_gphiu(gphiu_dim1, gphiu_dim2)
  !
  CALL momentum_u_code(i, j, ua, un, vn, hu, hv, ht, ssha_u, sshn_t, sshn_u, sshn_v, grid_tmask, grid_dx_u, grid_dx_v, grid_dx_t, grid_dy_u, grid_dy_t, grid_area_u, grid_gphiu)
  !
end subroutine momentum_u_code_wrap

subroutine momentum_v_code_wrap(i, j, va, va_dim1, va_dim2, un, un_dim1, un_dim2, vn, vn_dim1, vn_dim2, hu, hu_dim1, hu_dim2, hv, hv_dim1, hv_dim2, ht, ht_dim1, ht_dim2, ssha_v, ssha_v_dim1, ssha_v_dim2, sshn_t, sshn_t_dim1, sshn_t_dim2, sshn_u, sshn_u_dim1, sshn_u_dim2, sshn_v, sshn_v_dim1, sshn_v_dim2, grid_tmask, g_dim1, g_dim2, grid_dx_v, dxv_dim1, dxv_dim2, grid_dx_t, dxt_dim1, dxt_dim2, grid_dy_u, dyu_dim1, dyu_dim2, grid_dy_v, dyv_dim1, dyv_dim2, grid_dy_t, dyt_dim1, dyt_dim2, grid_area_v, gav_dim1, gav_dim2, grid_gphiv, gphiv_dim1, gphiv_dim2)
  USE momentum_mod, ONLY: momentum_v_code
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: i, j
  integer, intent(in)     :: va_dim1, va_dim2, un_dim1, un_dim2, vn_dim1, vn_dim2, hu_dim1, hu_dim2, hv_dim1, hv_dim2, ht_dim1, ht_dim2, ssha_v_dim1, ssha_v_dim2, sshn_t_dim1, sshn_t_dim2, sshn_u_dim1, sshn_u_dim2, sshn_v_dim1, sshn_v_dim2, g_dim1, g_dim2, dxv_dim1, dxv_dim2, dxt_dim1, dxt_dim2, dyu_dim1, dyu_dim2, dyv_dim1, dyv_dim2, dyt_dim1, dyt_dim2, gav_dim1, gav_dim2, gphiv_dim1, gphiv_dim2
  real(wp), intent(inout) :: va(va_dim1, va_dim2)
  real(wp), intent(in)    :: un(un_dim1, un_dim2)
  real(wp), intent(in)    :: vn(vn_dim1, vn_dim2)
  real(wp), intent(in)    :: hu(hu_dim1, hu_dim2)
  real(wp), intent(in)    :: hv(hv_dim1, hv_dim2)
  real(wp), intent(in)    :: ht(ht_dim1, ht_dim2)
  real(wp), intent(in)    :: ssha_v(ssha_v_dim1, ssha_v_dim2)
  real(wp), intent(in)    :: sshn_t(sshn_t_dim1, sshn_t_dim2)
  real(wp), intent(in)    :: sshn_u(sshn_u_dim1, sshn_u_dim2)
  real(wp), intent(in)    :: sshn_v(sshn_v_dim1, sshn_v_dim2)
  integer,  intent(in)    :: grid_tmask(g_dim1, g_dim2)
  real(wp), intent(in)    :: grid_dx_v(dxv_dim1, dxv_dim2)
  real(wp), intent(in)    :: grid_dx_t(dxt_dim1, dxt_dim2)
  real(wp), intent(in)    :: grid_dy_u(dyu_dim1, dyu_dim2)
  real(wp), intent(in)    :: grid_dy_v(dyv_dim1, dyv_dim2)
  real(wp), intent(in)    :: grid_dy_t(dyt_dim1, dyt_dim2)
  real(wp), intent(in)    :: grid_area_v(gav_dim1, gav_dim2)
  real(wp), intent(in)    :: grid_gphiv(gphiv_dim1, gphiv_dim2)
  !
  CALL momentum_v_code(i, j, va, un, vn, hu, hv, ht, ssha_v, sshn_t, sshn_u, sshn_v, grid_tmask, grid_dx_v, grid_dx_t, grid_dy_u, grid_dy_v, grid_dy_t, grid_area_v, grid_gphiv)
  !
end subroutine momentum_v_code_wrap

subroutine bc_ssh_code_wrap(i, j, istp, ssha_t, ssh_dim1, ssh_dim2, grid_tmask, g_dim1, g_dim2)
  USE boundary_conditions_mod, ONLY: bc_ssh_code
  implicit none
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: i, j, istp
  integer, intent(in)     :: ssh_dim1, ssh_dim2, g_dim1, g_dim2
  real(wp), intent(inout) :: ssha_t(ssh_dim1, ssh_dim2)
  integer, intent(in)     :: grid_tmask(g_dim1, g_dim2)
  !
  CALL bc_ssh_code(i, j, istp, ssha_t, grid_tmask)
  !
end subroutine bc_ssh_code_wrap

subroutine bc_solid_u_code_wrap(i, j, ua, u_dim1, u_dim2, grid_tmask, g_dim1, g_dim2)
  USE boundary_conditions_mod, ONLY: bc_solid_u_code
  implicit none
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: i, j
  integer, intent(in)     :: u_dim1, u_dim2, g_dim1, g_dim2
  real(wp), intent(inout) :: ua(u_dim1, u_dim2)
  integer, intent(in)     :: grid_tmask(g_dim1, g_dim2)
  !
  CALL bc_solid_u_code(i, j, ua, grid_tmask)
  !
end subroutine bc_solid_u_code_wrap


subroutine bc_solid_v_code_wrap(i, j, va, v_dim1, v_dim2, grid_tmask, g_dim1, g_dim2)
  USE boundary_conditions_mod, ONLY: bc_solid_v_code
  implicit none
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: i, j
  integer, intent(in)     :: v_dim1, v_dim2, g_dim1, g_dim2
  real(wp), intent(inout) :: va(v_dim1, v_dim2)
  integer, intent(in)     :: grid_tmask(g_dim1, g_dim2)
  !
  CALL bc_solid_v_code(i, j, va, grid_tmask)
  !
end subroutine bc_solid_v_code_wrap

subroutine bc_flather_u_code_wrap(i, j, ua, u_dim1, u_dim2, hu, hu_dim1, hu_dim2, sshn_u, ssh_dim1, ssh_dim2, grid_tmask, g_dim1, g_dim2)
  USE boundary_conditions_mod, ONLY: bc_flather_u_code
  implicit none
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: i, j
  integer, intent(in)     :: u_dim1, u_dim2, hu_dim1, hu_dim2, ssh_dim1, ssh_dim2, g_dim1, g_dim2
  real(wp), intent(inout) :: ua(u_dim1, u_dim2)
  real(wp), intent(in)    :: hu(hu_dim1, hu_dim2)
  real(wp), intent(in)    :: sshn_u(ssh_dim1, ssh_dim2)
  integer, intent(in)     :: grid_tmask(g_dim1, g_dim2)
  !
  CALL bc_flather_u_code(i, j, ua, hu, sshn_u, grid_tmask)
  !
end subroutine bc_flather_u_code_wrap

subroutine bc_flather_v_code_wrap(i, j, va, v_dim1, v_dim2, hv, hv_dim1, hv_dim2, sshn_v, ssh_dim1, ssh_dim2, grid_tmask, g_dim1, g_dim2)
  USE boundary_conditions_mod, ONLY: bc_flather_v_code
  implicit none
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: i, j
  integer, intent(in)     :: v_dim1, v_dim2, hv_dim1, hv_dim2, ssh_dim1, ssh_dim2, g_dim1, g_dim2
  real(wp), intent(inout) :: va(v_dim1, v_dim2)
  real(wp), intent(in)    :: hv(hv_dim1, hv_dim2)
  real(wp), intent(in)    :: sshn_v(ssh_dim1, ssh_dim2)
  integer, intent(in)     :: grid_tmask(g_dim1, g_dim2)
  !
  CALL bc_flather_v_code(i, j, va, hv, sshn_v, grid_tmask)
  !
end subroutine bc_flather_v_code_wrap

subroutine field_copy_code_wrap(i, j, field_out, out_dim1, out_dim2, field_in, in_dim1, in_dim2)
  USE infrastructure_mod, ONLY: field_copy_code
  implicit none
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: out_dim1, out_dim2, in_dim1, in_dim2
  integer, intent(in)     :: i, j
  real(wp), intent(inout) :: field_out(out_dim1, out_dim2)
  real(wp), intent(in)    :: field_in(in_dim1, in_dim2)
  !
  CALL field_copy_code(i, j, field_out, field_in)
  !
end subroutine field_copy_code_wrap

subroutine  next_sshu_code_wrap(i, j, sshn_u, u_dim1, u_dim2, sshn_t, t_dim1, t_dim2, grid_tmask, gtm_dim1, gtm_dim2, grid_area_t, gat_dim1, gat_dim2, grid_area_u, gau_dim1, gau_dim2)
  USE time_update_mod, ONLY: next_sshu_code
  implicit none
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer, intent(in)     :: i, j
  integer, intent(in)     :: u_dim1, u_dim2, t_dim1, t_dim2, gtm_dim1, gtm_dim2, gat_dim1, gat_dim2, gau_dim1, gau_dim2
  real(wp), intent(inout) :: sshn_u(u_dim1, u_dim2)
  real(wp), intent(in)    :: sshn_t(t_dim1, t_dim2)
  integer, intent(in)     :: grid_tmask(gtm_dim1, gtm_dim2)
  real(wp), intent(in)    :: grid_area_t(gat_dim1, gat_dim2)
  real(wp), intent(in)    :: grid_area_u(gau_dim1, gau_dim2)
  !
  CALL next_sshu_code(i, j, sshn_u, sshn_t, grid_tmask, grid_area_t, grid_area_u)
  !
end subroutine next_sshu_code_wrap

subroutine next_sshv_code_wrap(i, j, sshn_v, v_dim1, v_dim2, sshn_t, t_dim1, t_dim2, grid_tmask, gtm_dim1, gtm_dim2, grid_area_t, gat_dim1, gat_dim2, grid_area_v, gav_dim1, gav_dim2)
  USE time_update_mod, ONLY: next_sshv_code
  implicit none
  INTEGER, PARAMETER      :: wp = SELECTED_REAL_KIND(12,307)
  integer,  intent(in)    :: i, j
  integer,  intent(in)    :: v_dim1, v_dim2, t_dim1, t_dim2, gtm_dim1, gtm_dim2, gat_dim1, gat_dim2, gav_dim1, gav_dim2
  integer,  intent(in)    :: grid_tmask(gtm_dim1, gtm_dim2)
  real(wp), intent(in)    :: grid_area_t(gat_dim1, gat_dim2)
  real(wp), intent(in)    :: grid_area_v(gav_dim1, gav_dim2)
  real(wp), intent(inout) :: sshn_v(v_dim1, v_dim2)
  real(wp), intent(in)    :: sshn_t(t_dim1, t_dim2)
  !
  CALL next_sshv_code(i, j, sshn_v, sshn_t, grid_tmask, grid_area_t, grid_area_v)
  !
end subroutine next_sshv_code_wrap
