  MODULE psy_gocean2d
    USE field_mod
    USE kind_params_mod
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE invoke_0(ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, rdt, ua, ht, ssha_u, va, ssha_v, istp)
      USE time_update_mod, ONLY: next_sshv_code
      USE time_update_mod, ONLY: next_sshu_code
      USE infrastructure_mod, ONLY: field_copy_code
      USE boundary_conditions_mod, ONLY: bc_flather_v_code
      USE boundary_conditions_mod, ONLY: bc_flather_u_code
      USE boundary_conditions_mod, ONLY: bc_solid_v_code
      USE boundary_conditions_mod, ONLY: bc_solid_u_code
      USE boundary_conditions_mod, ONLY: bc_ssh_code
      USE momentum_mod, ONLY: momentum_v_code
      USE momentum_mod, ONLY: momentum_u_code
      USE continuity_mod, ONLY: continuity_code
      TYPE(r2d_field), intent(inout) :: ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
      REAL(KIND=wp), intent(inout) :: rdt
      INTEGER, intent(inout) :: istp
      INTEGER j
      INTEGER i
      INTEGER istop, jstop
      !
      ! Look-up loop bounds
      istop = ssha_t%grid%simulation_domain%xstop
      jstop = ssha_t%grid%simulation_domain%ystop
      !
      DO j=2,jstop
        DO i=2,istop
          CALL continuity_code(i, j, ssha_t%data, sshn_t%data, sshn_u%data, sshn_v%data, hu%data, hv%data, un%data, vn%data, rdt, sshn_t%grid%area_t)
        END DO 
      END DO 
      DO j=2,jstop
        DO i=2,istop-1
          CALL momentum_u_code(i, j, ua%data, un%data, vn%data, hu%data, hv%data, ht%data, ssha_u%data, sshn_t%data, sshn_u%data, sshn_v%data, un%grid%tmask, un%grid%dx_u, un%grid%dx_v, un%grid%dx_t, un%grid%dy_u, un%grid%dy_t, un%grid%area_u, un%grid%gphiu)
        END DO 
      END DO 
      DO j=2,jstop-1
        DO i=2,istop
          CALL momentum_v_code(i, j, va%data, un%data, vn%data, hu%data, hv%data, ht%data, ssha_v%data, sshn_t%data, sshn_u%data, sshn_v%data, un%grid%tmask, un%grid%dx_v, un%grid%dx_t, un%grid%dy_u, un%grid%dy_v, un%grid%dy_t, un%grid%area_v, un%grid%gphiv)
        END DO 
      END DO 
      DO j=2,jstop
        DO i=2,istop
          CALL bc_ssh_code(i, j, istp, ssha_t%data, ssha_t%grid%tmask)
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop
          CALL bc_solid_u_code(i, j, ua%data, ua%grid%tmask)
        END DO 
      END DO 
      DO j=1,jstop
        DO i=1,istop+1
          CALL bc_solid_v_code(i, j, va%data, va%grid%tmask)
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop
          CALL bc_flather_u_code(i, j, ua%data, hu%data, sshn_u%data, hu%grid%tmask)
        END DO 
      END DO 
      DO j=1,jstop
        DO i=1,istop+1
          CALL bc_flather_v_code(i, j, va%data, hv%data, sshn_v%data, hv%grid%tmask)
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop+1
          CALL field_copy_code(i, j, un%data, ua%data)
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop+1
          CALL field_copy_code(i, j, vn%data, va%data)
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop+1
          CALL field_copy_code(i, j, sshn_t%data, ssha_t%data)
        END DO 
      END DO 
      DO j=2,jstop
        DO i=2,istop-1
          CALL next_sshu_code(i, j, sshn_u%data, sshn_t%data, sshn_t%grid%tmask, sshn_t%grid%area_t, sshn_t%grid%area_u)
        END DO 
      END DO 
      DO j=2,jstop-1
        DO i=2,istop
          CALL next_sshv_code(i, j, sshn_v%data, sshn_t%data, sshn_t%grid%tmask, sshn_t%grid%area_t, sshn_t%grid%area_v)
        END DO 
      END DO 
    END SUBROUTINE invoke_0
  END MODULE psy_gocean2d