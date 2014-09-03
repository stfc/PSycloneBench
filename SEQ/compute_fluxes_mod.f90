module compute_fluxes_mod
  use kind_params_mod
  use topology_mod, only: cu, cv, ct, cf
  use compute_cu_mod, only: compute_cu_code
  use compute_cv_mod, only: compute_cv_code
  use compute_z_mod, only: compute_z_code
  use compute_h_mod, only: compute_h_code
  implicit none
  private

  public invoke_compute_fluxes

contains

  !===================================================

  subroutine invoke_compute_fluxes(cufld, cvfld, zfld, hfld, &
                                   pfld, ufld, vfld)
    implicit none
    real(wp), intent(inout), dimension(:,:) :: cufld, cvfld, zfld, hfld
    real(wp), intent(in),    dimension(:,:) :: pfld, ufld, vfld
    ! Locals
    integer :: i, j

    !CALL invoke_compute_cu(CU, P, U)
    ! Change to hard-code lower-loop bounds with variable
    ! trip count expressed in terms of M/N+/-1
    ! M/N obtained from topology  look-up
    do J=cu%jstart, cu%jstop
       do I=cu%istart, cu%istop

          call compute_cu_code(i, j, cufld, pfld, ufld)
       end do
    end do

    !CALL invoke_compute_cv(CV, P, V)
    do J=cv%jstart, cv%jstop !2, size(cv, 2)
       do I=cv%istart, cv%istop !1, size(cv, 1) - 1

          call compute_cv_code(i, j, cvfld, pfld, vfld)
       end do
    end do

    !CALL invoke_compute_z(z, P, U, V)
    do J=cf%jstart, cf%jstop !2, size(z, 2)
       do I=cf%istart, cf%istop !2, size(z, 1)

          call compute_z_code(i, j, zfld, pfld, ufld, vfld)
       end do
    end do

    !CALL invoke_compute_h(h, P, U, V)
    DO J=ct%jstart, ct%jstop !1, SIZE(h, 2) - 1
       DO I=ct%istart, ct%istop !1, SIZE(h, 1) - 1

          CALL compute_h_code(i, j, hfld, pfld, ufld, vfld)
       END DO
    END DO

  end subroutine invoke_compute_fluxes

  !===================================================

end module compute_fluxes_mod
