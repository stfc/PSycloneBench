MODULE manual_invoke_apply_bcs_mod
  USE kind_params
  IMPLICIT none
  PRIVATE

  PUBLIC manual_invoke_apply_bcs_uvtf
  PUBLIC manual_invoke_apply_bcs_uvt

CONTAINS

  !===================================================

  SUBROUTINE manual_invoke_apply_bcs_uvtf(ufield, vfield, tfield, ffield)
    USE apply_bcs_cf, ONLY: apply_bcs_cf_code
    USE apply_bcs_ct, ONLY: apply_bcs_ct_code
    USE apply_bcs_cu, ONLY: apply_bcs_cu_code
    USE apply_bcs_cv, ONLY: apply_bcs_cv_code
    IMPLICIT none
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: ufield, vfield, tfield, ffield
    ! Locals
    integer :: m, n, mp1, np1

    MP1 = SIZE(ufield, 1)
    NP1 = SIZE(ufield, 2)
    M = MP1 - 1
    N = NP1 - 1

    call apply_bcs_cu_code(n, mp1, np1, ufield)
    call apply_bcs_ct_code(n, mp1, np1, tfield)
    call apply_bcs_cv_code(m, mp1, np1, vfield)
    call apply_bcs_cf_code(mp1, np1, ffield)

  END SUBROUTINE manual_invoke_apply_bcs_uvtf

  !===================================================

  subroutine manual_invoke_apply_bcs_uvt(ufield, vfield, tfield)
    use apply_bcs_ct, only: apply_bcs_ct_code
    use apply_bcs_cu, only: apply_bcs_cu_code
    use apply_bcs_cv, only: apply_bcs_cv_code
    implicit none
    real(wp), intent(inout), dimension(:,:) :: ufield, vfield, tfield
    ! Locals
    integer :: m, n, mp1, np1

    MP1 = SIZE(ufield, 1)
    NP1 = SIZE(ufield, 2)
    M = MP1 - 1
    N = NP1 - 1

    call apply_bcs_cu_code(n, mp1, np1, ufield)
    call apply_bcs_cv_code(m, mp1, np1, vfield)
    call apply_bcs_ct_code(n, mp1, np1, tfield)

  end subroutine manual_invoke_apply_bcs_uvt

  !===================================================

END MODULE manual_invoke_apply_bcs_mod
