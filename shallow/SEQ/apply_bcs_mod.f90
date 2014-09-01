MODULE apply_bcs_mod
  USE kind_params_mod
  IMPLICIT none
  PRIVATE

  PUBLIC invoke_apply_bcs_uvtf
  PUBLIC invoke_apply_bcs_uvt

CONTAINS

  !===================================================

  SUBROUTINE invoke_apply_bcs_uvtf(ufield, vfield, tfield, ffield)
    USE apply_bcs_cf_mod, ONLY: invoke_apply_bcs_cf
    use apply_bcs_ct_mod, only: invoke_apply_bcs_ct
    use apply_bcs_cu_mod, only: invoke_apply_bcs_cu
    use apply_bcs_cv_mod, only: invoke_apply_bcs_cv
    IMPLICIT none
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: ufield, vfield, tfield, ffield

    call invoke_apply_bcs_cu(ufield)
    call invoke_apply_bcs_ct(tfield)
    call invoke_apply_bcs_cv(vfield)
    call invoke_apply_bcs_cf(ffield)

  END SUBROUTINE invoke_apply_bcs_uvtf

  !===================================================

  subroutine invoke_apply_bcs_uvt(ufield, vfield, tfield)
    use apply_bcs_ct_mod, only: invoke_apply_bcs_ct
    use apply_bcs_cu_mod, only: invoke_apply_bcs_cu
    use apply_bcs_cv_mod, only: invoke_apply_bcs_cv
    implicit none
    real(wp), intent(inout), dimension(:,:) :: ufield, vfield, tfield

    call invoke_apply_bcs_cu(ufield)
    call invoke_apply_bcs_cv(vfield)
    call invoke_apply_bcs_ct(tfield)

  end subroutine invoke_apply_bcs_uvt

  !===================================================

END MODULE apply_bcs_mod
