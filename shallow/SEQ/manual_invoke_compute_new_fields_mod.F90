MODULE manual_invoke_compute_new_fields_mod
  USE kind_params_mod
  IMPLICIT none
  PRIVATE

  PUBLIC manual_invoke_compute_new_fields

CONTAINS

  SUBROUTINE manual_invoke_compute_new_fields(unew, uold, &
                                              vnew, vold, &
                                              pnew, pold, &
                                              z, cufld, cvfld, hfld, tdt)
    USE compute_unew_mod, ONLY: compute_unew_code
    USE compute_vnew_mod, ONLY: compute_vnew_code
    USE compute_pnew_mod, ONLY: compute_pnew_code
    use topology_mod,     ONLY: cu_grid, cv_grid, ct_grid
    IMPLICIT none
    REAL(wp), INTENT(out), DIMENSION(:,:) :: unew, vnew, pnew
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: uold, vold, pold
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: z, cufld, cvfld, hfld
    REAL(wp), INTENT(in) :: tdt
    ! Locals
    integer :: i, j

    !CALL manual_invoke_compute_unew(unew, uold,  z, cv, h, tdt)
    DO J=cu_grid%jstart, cu_grid%jstop 
       DO I=cu_grid%istart, cu_grid%istop

          CALL compute_unew_code(i, j, unew, uold, &
                                 z, cvfld, hfld, tdt)

       END DO
    END DO
    !CALL manual_invoke_compute_vnew(vnew, vold,  z, cu, h, tdt)
    DO J=cv_grid%jstart, cv_grid%jstop
       DO I=cv_grid%istart, cv_grid%istop

          CALL compute_vnew_code(i, j, vnew, vold, &
                                 z, cufld, hfld, tdt)
       END DO
    END DO
    !CALL manual_invoke_compute_pnew(pnew, pold, cu, cv,    tdt)
    DO J=ct_grid%jstart, ct_grid%jstop
       DO I=ct_grid%istart, ct_grid%istop

          CALL compute_pnew_code(i, j, pnew, pold, &
                                 cufld, cvfld, tdt)
       END DO
    END DO

  END SUBROUTINE manual_invoke_compute_new_fields

END MODULE manual_invoke_compute_new_fields_mod
