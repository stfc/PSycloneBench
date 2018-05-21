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
     !$omp task out(ssha_t), in(sshn_t, sshn_u, sshn_v, hu, hv, un, vn, grid_area_t)
     DO i=2,istop
        CALL continuity_code_wrap(i, j, ssha_t, size(ssha_t,1), size(ssha_t,2), sshn_t, size(sshn_t,1), size(sshn_t,2), sshn_u, size(sshn_u,1), size(sshn_u,2), sshn_v, size(sshn_v,1), size(sshn_v,2), hu, size(hu,1), size(hu,2), hv, size(hv,1), size(hv,2), un, size(un,1), size(un,2), vn, size(vn,1), size(vn,2), rdt, grid_area_t, size(grid_area_t,1), size(grid_area_t,2))
     END DO
     !$omp end task
  END DO
  DO j=2,jstop
     !$omp task inout(ua) in(un, vn, hu, hv, ht, ssha_u, sshn_t, sshn_u, sshn_v, grid_tmask, grid_dx_u, grid_dx_v, grid_dx_t, grid_dy_u, grid_dy_t, grid_area_u, grid_gphiu)
     DO i=2,istop-1
        CALL momentum_u_code_wrap(i, j, ua, size(ua,1), size(ua,2), un, size(un,1), size(un,2), vn, size(vn,1), size(vn,2), hu, size(hu,1), size(hu,2), hv, size(hv,1), size(hv,2), ht, size(ht,1), size(ht,2), ssha_u, size(ssha_u,1), size(ssha_u,2), sshn_t, size(sshn_t,1), size(sshn_t,2), sshn_u, size(sshn_u,1), size(sshn_u,2), sshn_v, size(sshn_v,1), size(sshn_v,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2), grid_dx_u, size(grid_dx_u,1), size(grid_dx_u,2), grid_dx_v, size(grid_dx_v,1), size(grid_dx_v,2), grid_dx_t, size(grid_dx_t,1), size(grid_dx_t,2), grid_dy_u, size(grid_dy_u,1), size(grid_dy_u,2), grid_dy_t, size(grid_dy_t,1), size(grid_dy_t,2), grid_area_u, size(grid_area_u,1), size(grid_area_u,2), grid_gphiu, size(grid_gphiu,1), size(grid_gphiu,2))
     END DO
     !$omp end task
  END DO
  DO j=2,jstop-1
     !$omp task inout(va) in(un, vn, hu, hv, ht, ssha_v, sshn_t, sshn_u, sshn_v, grid_tmask, grid_dx_v, grid_dx_t, grid_dy_u, grid_dy_v, grid_dy_t, grid_area_v, grid_gphiv)
     DO i=2,istop
        CALL momentum_v_code_wrap(i, j, va, size(va,1), size(va,2), un, size(un,1), size(un,2), vn, size(vn,1), size(vn,2), hu, size(hu,1), size(hu,2), hv, size(hv,1), size(hv,2), ht, size(ht,1), size(ht,2), ssha_v, size(ssha_v,1), size(ssha_v,2), sshn_t, size(sshn_t,1), size(sshn_t,2), sshn_u, size(sshn_u,1), size(sshn_u,2), sshn_v, size(sshn_v,1), size(sshn_v,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2), grid_dx_v, size(grid_dx_v,1), size(grid_dx_v,2), grid_dx_t, size(grid_dx_t,1), size(grid_dx_t,2), grid_dy_u, size(grid_dy_u,1), size(grid_dy_u,2), grid_dy_v, size(grid_dy_v,1), size(grid_dy_v,2), grid_dy_t, size(grid_dy_t,1), size(grid_dy_t,2), grid_area_v, size(grid_area_v,1), size(grid_area_v,2), grid_gphiv, size(grid_gphiv,1), size(grid_gphiv,2))
     END DO
     !$omp end task
  END DO
  DO j=2,jstop
     !$omp task inout(ssha_t) in(grid_tmask)
     DO i=2,istop
        CALL bc_ssh_code_wrap(i, j, istp, ssha_t, size(ssha_t,1), size(ssha_t,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
     END DO
     !$omp end task
  END DO
  DO j=1,jstop+1
     !$omp task inout(ua) in(grid_tmask)
     DO i=1,istop
        CALL bc_solid_u_code_wrap(i, j, ua, size(ua,1), size(ua,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
     END DO
     !$omp end task
  END DO
  DO j=1,jstop
     !$omp task inout(va) in(grid_tmask)
     DO i=1,istop+1
        CALL bc_solid_v_code_wrap(i, j, va, size(va,1), size(va,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
     END DO
     !$omp end task
  END DO
  DO j=1,jstop+1
     !$omp task inout(ua) in(hu, sshn_u, grid_tmask)
     DO i=1,istop
        CALL bc_flather_u_code_wrap(i, j, ua, size(ua,1), size(ua,2), hu, size(hu,1), size(hu,2), sshn_u, size(sshn_u,1), size(sshn_u,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
      END DO
      !$omp end task
  END DO
  DO j=1,jstop
     !$omp task inout(va) in(hv, sshn_v, grid_tmask)
     DO i=1,istop+1
        CALL bc_flather_v_code_wrap(i, j, va, size(va,1), size(va,2), hv, size(hv,1), size(hv,2), sshn_v, size(sshn_v,1), size(sshn_v,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2))
     END DO
     !$omp end task
  END DO
  DO j=1,jstop+1
     !$omp task out(un) in(ua)
     DO i=1,istop+1
        CALL field_copy_code_wrap(i, j, un, size(un,1), size(un,2), ua, size(ua,1), size(ua,2))
     END DO
     !$omp end task
  END DO
  DO j=1,jstop+1
     !$omp task out(vn) in(va)
     DO i=1,istop+1
        CALL field_copy_code_wrap(i, j, vn, size(vn,1), size(vn,2), va, size(va,1), size(va,2))
     END DO
     !$omp end task
  END DO
  DO j=1,jstop+1
     !$omp task out(sshn_t) in(ssha_t)
     DO i=1,istop+1
        CALL field_copy_code_wrap(i, j, sshn_t, size(sshn_t,1), size(sshn_t,2), ssha_t, size(ssha_t,1), size(ssha_t,2))
     END DO
     !$omp end task
  END DO
  DO j=2,jstop
     !$omp task inout(sshn_u) in(sshn_t, grid_tmask, grid_area_t, grid_area_u)
     DO i=2,istop-1
        CALL next_sshu_code_wrap(i, j, sshn_u, size(sshn_u,1), size(sshn_u,2), sshn_t, size(sshn_t,1), size(sshn_t,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2), grid_area_t, size(grid_area_t,1), size(grid_area_t,2), grid_area_u, size(grid_area_u,1), size(grid_area_u,2))
     END DO
     !$omp end task
  END DO
  DO j=2,jstop-1
     !$omp task inout(sshn_v) in(sshn_t, grid_tmask, grid_area_t, grid_area_v)
     DO i=2,istop
        CALL next_sshv_code_wrap(i, j, sshn_v, size(sshn_v,1), size(sshn_v,2), sshn_t, size(sshn_t,1), size(sshn_t,2), grid_tmask, size(grid_tmask,1), size(grid_tmask,2), grid_area_t, size(grid_area_t,1), size(grid_area_t,2), grid_area_v, size(grid_area_v,1), size(grid_area_v,2))
     END DO
     !$omp end task
  END DO
  
END SUBROUTINE invoke_0_deref
