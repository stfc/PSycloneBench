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
    use compute_unew_mod, ONLY: compute_unew_code
    use compute_vnew_mod, ONLY: compute_vnew_code
    use compute_pnew_mod, ONLY: compute_pnew_code
    use field_mod
    implicit none
    type(r2d_field_type), intent(out) :: unew, vnew, pnew
    type(r2d_field_type), intent(in) :: uold, vold, pold
    type(r2d_field_type), intent(in) :: z, cufld, cvfld, hfld
    real(wp), intent(in) :: tdt
    ! Locals
    integer :: i, j

    !CALL manual_invoke_compute_unew(unew, uold,  z, cv, h, tdt)
    DO J=unew%internal%ystart, unew%internal%ystop 
       DO I=unew%internal%xstart, unew%internal%xstop

          CALL compute_unew_code(i, j, unew%data, uold%data, &
                                 z%data, cvfld%data, hfld%data, tdt)

       END DO
    END DO
    !CALL manual_invoke_compute_vnew(vnew, vold,  z, cu, h, tdt)
    DO J=vnew%internal%ystart, vnew%internal%ystop
       DO I=vnew%internal%xstart, vnew%internal%xstop

          CALL compute_vnew_code(i, j, vnew%data, vold%data, &
                                 z%data, cufld%data, hfld%data, tdt)
       END DO
    END DO
    !CALL manual_invoke_compute_pnew(pnew, pold, cu, cv,    tdt)
    DO J=pnew%internal%ystart, pnew%internal%ystop
       DO I=pnew%internal%xstart, pnew%internal%xstop

          CALL compute_pnew_code(i, j, pnew%data, pold%data, &
                                 cufld%data, cvfld%data, tdt)
       END DO
    END DO

  END SUBROUTINE manual_invoke_compute_new_fields

END MODULE manual_invoke_compute_new_fields_mod
