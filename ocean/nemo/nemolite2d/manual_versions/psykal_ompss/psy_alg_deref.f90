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
