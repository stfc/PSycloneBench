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
          !$omp task out(ssha_t%data(i,j)) in(sshn_t%data(i,j), sshn_u%data(i,j), sshn_u%data(i-1,j), sshn_v%data(i,j), sshn_v%data(i,j-1), hu%data(i,j), hu%data(i-1,j), hv%data(i,j), hv%data(i,j-1), un%data(i,j), un%data(i-1,j), vn%data(i,j), vn%data(i,j-1))
          CALL continuity_code(i, j, ssha_t%data, sshn_t%data, sshn_u%data, sshn_v%data, hu%data, hv%data, un%data, vn%data, rdt, sshn_t%grid%area_t)
          !$omp end task
        END DO 
      END DO 
      DO j=2,jstop
        DO i=2,istop-1
          !$omp task inout(ua%data(i,j)), in(un%data(i,j), un%data(i+1,j), un%data(i-1,j), un%data(i,j-1), un%data(i,j+1), vn%data(i,j), vn%data(i,j-1), vn%data(i+1,j-1), vn%data(i+1,j), hu%data(i,j), hu%data(i,j-1), hu%data(i,j+1), hv%data(i,j), hv%data(i,j-1), hv%data(i+1,j-1), hv%data(i+1,j), ht%data(i,j), ht%data(i+1,j), ssha_u%data(i,j), sshn_t%data(i,j), sshn_t%data(i+1,j), sshn_u%data(i,j), sshn_u%data(i,j-1), sshn_u%data(i,j+1), sshn_v%data(i,j), sshn_v%data(i,j-1), sshn_v%data(i+1,j-1), sshn_v%data(i+1,j))
          CALL momentum_u_code(i, j, ua%data, un%data, vn%data, hu%data, hv%data, ht%data, ssha_u%data, sshn_t%data, sshn_u%data, sshn_v%data, un%grid%tmask, un%grid%dx_u, un%grid%dx_v, un%grid%dx_t, un%grid%dy_u, un%grid%dy_t, un%grid%area_u, un%grid%gphiu)
          !$omp end task
        END DO
      END DO 
      DO j=2,jstop-1
        DO i=2,istop
          !$omp task inout(va%data(i,j)) in(***TBD un, vn, hu, hv, ht, ssha_v, sshn_t, sshn_u, sshn_v)
          CALL momentum_v_code(i, j, va%data, un%data, vn%data, hu%data, hv%data, ht%data, ssha_v%data, sshn_t%data, sshn_u%data, sshn_v%data, un%grid%tmask, un%grid%dx_v, un%grid%dx_t, un%grid%dy_u, un%grid%dy_v, un%grid%dy_t, un%grid%area_v, un%grid%gphiv)
          !$omp end task
        END DO
      END DO 
      DO j=2,jstop
        DO i=2,istop
          !$omp task inout(ssha_t%data(i,j))
           CALL bc_ssh_code(i, j, istp, ssha_t%data, ssha_t%grid%tmask)
          !$omp end task
        END DO
      END DO 
      DO j=1,jstop+1
        DO i=1,istop
          !$omp task inout(ua%data(i,j))
          CALL bc_solid_u_code(i, j, ua%data, ua%grid%tmask)
          !$omp end task
        END DO 
      END DO 
      DO j=1,jstop
        DO i=1,istop+1
          !$omp task inout(va%data(i,j))
          CALL bc_solid_v_code(i, j, va%data, va%grid%tmask)
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop
          !$omp task inout(ua%data(i,j)) in(hu, sshn_u)???????? writer has a stencil for boundary
          CALL bc_flather_u_code(i, j, ua%data, hu%data, sshn_u%data, hu%grid%tmask)
          !$omp end task
        END DO 
      END DO 
      DO j=1,jstop
        DO i=1,istop+1
          !$omp task inout(va%data(i,j)) in(hv, sshn_v)???????? writer has a stencil for boundary
          CALL bc_flather_v_code(i, j, va%data, hv%data, sshn_v%data, hv%grid%tmask)
          !$omp end task
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop+1
          !$omp task out(un%data(i,j)) in(ua%data(i,j))
          CALL field_copy_code(i, j, un%data, ua%data)
          !$omp end task
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop+1
          !$omp task out(vn%data(i,j)) in(va%data(i,j))
          CALL field_copy_code(i, j, vn%data, va%data)
          !$omp end task
        END DO 
      END DO 
      DO j=1,jstop+1
        DO i=1,istop+1
          !$omp task out(sshn_t%data(i,j)) in(ssha_t%data(i,j))
          CALL field_copy_code(i, j, sshn_t%data, ssha_t%data)
          !$omp end task
        END DO 
      END DO 
      DO j=2,jstop
        DO i=2,istop-1
          !$omp task inout(sshn_u%data(i,j)) in(sshn_t%data(i,j), sshn_t%data(i+1,j))
          CALL next_sshu_code(i, j, sshn_u%data, sshn_t%data, sshn_t%grid%tmask, sshn_t%grid%area_t, sshn_t%grid%area_u)
          !$omp end task
        END DO 
      END DO 
      DO j=2,jstop-1
        DO i=2,istop
          !$omp task inout(sshn_v%data(i,j)) in(sshn_t%data(i,j), sshn_t%data(i,j+1))
          CALL next_sshv_code(i, j, sshn_v%data, sshn_t%data, sshn_t%grid%tmask, sshn_t%grid%area_t, sshn_t%grid%area_v)
          !$omp end task
        END DO 
      END DO
     
      ! synchronise at the end of the invoke to be safe
      !$omp taskwait

    END SUBROUTINE invoke_0
  END MODULE psy_gocean2d
