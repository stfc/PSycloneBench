MODULE manual_invoke_apply_bcs_mod
  USE kind_params
  IMPLICIT none
  PRIVATE

  PUBLIC manual_invoke_apply_bcs_uvtf
  PUBLIC manual_invoke_apply_bcs_uvt

CONTAINS

  !===================================================

  SUBROUTINE manual_invoke_apply_bcs_uvtf(ufield, vfield, tfield, ffield)
    USE apply_bcs_cf_mod, ONLY: manual_invoke_apply_bcs_cf
    use apply_bcs_ct_mod, only: manual_invoke_apply_bcs_ct
    use apply_bcs_cu_mod, only: manual_invoke_apply_bcs_cu
    use apply_bcs_cv_mod, only: manual_invoke_apply_bcs_cv
    IMPLICIT none
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: ufield, vfield, tfield, ffield

    call manual_invoke_apply_bcs_cu(ufield)
    call manual_invoke_apply_bcs_ct(tfield)
    call manual_invoke_apply_bcs_cv(vfield)
    call manual_invoke_apply_bcs_cf(ffield)

  END SUBROUTINE manual_invoke_apply_bcs_uvtf

  !===================================================

  subroutine manual_invoke_apply_bcs_uvt(ufield, vfield, tfield)
    use apply_bcs_ct_mod, only: manual_invoke_apply_bcs_ct
    use apply_bcs_cu_mod, only: manual_invoke_apply_bcs_cu
    use apply_bcs_cv_mod, only: manual_invoke_apply_bcs_cv
    implicit none
    real(wp), intent(inout), dimension(:,:) :: ufield, vfield, tfield

    call manual_invoke_apply_bcs_cu(ufield)
    call manual_invoke_apply_bcs_cv(vfield)
    call manual_invoke_apply_bcs_ct(tfield)

  end subroutine manual_invoke_apply_bcs_uvt

  !===================================================

END MODULE manual_invoke_apply_bcs_mod
