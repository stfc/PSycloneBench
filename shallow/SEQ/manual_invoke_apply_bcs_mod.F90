module manual_invoke_apply_bcs_mod
  use kind_params_mod
  use field_mod
  implicit none

  private

  public manual_invoke_apply_bcs

contains

  !===================================================

  subroutine manual_invoke_apply_bcs(field)
    implicit none
    type(r2d_field), intent(inout) :: field
    ! Locals
    integer :: ihalo

!DIR$ LOOP_INFO max_trips(2)
    do ihalo = 1, field%num_halos

      ! Copy from source to destination
      call copy_field(field,                    &
                      field%halo(ihalo)%source, &
                      field%halo(ihalo)%dest)
    end do

  end subroutine manual_invoke_apply_bcs

  !===================================================

end module manual_invoke_apply_bcs_mod
