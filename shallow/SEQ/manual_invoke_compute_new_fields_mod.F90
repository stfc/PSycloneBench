MODULE manual_invoke_compute_new_fields_mod
  USE kind_params_mod
  IMPLICIT none
  PRIVATE

  PUBLIC invoke_compute_new_fields

CONTAINS

  SUBROUTINE invoke_compute_new_fields(unew, uold, &
                                              vnew, vold, &
                                              pnew, pold, &
                                              z, cufld, cvfld, hfld, tdt)
    USE compute_unew_mod, ONLY: compute_unew_code
    USE compute_vnew_mod, ONLY: compute_vnew_code
    USE compute_pnew_mod, ONLY: compute_pnew_code
    use topology_mod,     ONLY: cu, cv, ct
    IMPLICIT none
    REAL(wp), INTENT(out), DIMENSION(:,:) :: unew, vnew, pnew
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: uold, vold, pold
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: z, cufld, cvfld, hfld
    REAL(wp), INTENT(in) :: tdt
    ! Locals
    integer :: i, j

    !CALL manual_invoke_compute_unew(unew, uold,  z, cv, h, tdt)
    DO J=cu%jstart, cu%jstop 
       DO I=cu%istart, cu%istop

          CALL compute_unew_code(i, j, unew, uold, &
                                 z, cvfld, hfld, tdt)

       END DO
    END DO
    !CALL manual_invoke_compute_vnew(vnew, vold,  z, cu, h, tdt)
    DO J=cv%jstart, cv%jstop
       DO I=cv%istart, cv%istop

          CALL compute_vnew_code(i, j, vnew, vold, &
                                 z, cufld, hfld, tdt)
       END DO
    END DO
    !CALL manual_invoke_compute_pnew(pnew, pold, cu, cv,    tdt)
    DO J=ct%jstart, ct%jstop
       DO I=ct%istart, ct%istop

          CALL compute_pnew_code(i, j, pnew, pold, &
                                 cufld, cvfld, tdt)
       END DO
    END DO

  END SUBROUTINE invoke_compute_new_fields

END MODULE manual_invoke_compute_new_fields_mod
