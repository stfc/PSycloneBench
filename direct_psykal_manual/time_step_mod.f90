module time_step_mod
  implicit none

  private

  public invoke_time_step

contains

  subroutine invoke_time_step(istp, ssha, ssha_u, ssha_v, &
                              sshn, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)
    use field_mod
    use momentum_mod, only: invoke_momentum_u, invoke_momentum_v
    use continuity_mod, only: invoke_continuity
    use time_update_mod, only: invoke_next_sshu, invoke_next_sshv
    use boundary_conditions_mod
    implicit none
    integer,         intent(in) :: istp
    type(r2d_field), intent(inout) :: un, vn, sshn, sshn_u, sshn_v
    type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
    type(r2d_field), intent(in) :: hu, hv, ht

  CALL invoke_continuity(ssha, sshn, sshn_u, sshn_v, &
                         hu, hv, un, vn)

  CALL invoke_momentum_u(ua, un, vn, &
                         hu, hv, ht, &
                         ssha_u, sshn, sshn_u, sshn_v)

  CALL invoke_momentum_v(va, un, vn, &
                         hu, hv, ht, &
                         ssha_v, sshn, sshn_u, sshn_v)

  ! Apply open and solid boundary conditions
  CALL invoke_bc_ssh(istp, ssha)
  CALL invoke_bc_solid_u(ua)
  CALL invoke_bc_solid_v(va)
  CALL invoke_bc_flather_u(ua, hu, sshn_u)
  CALL invoke_bc_flather_v(va, hv, sshn_v)

  ! Time update of fields
  call copy_field(ua, un)
  call copy_field(va, vn)
  call copy_field(ssha, sshn)
  call invoke_next_sshu(sshn_u, sshn)
  call invoke_next_sshv(sshn_v, sshn)

  end subroutine invoke_time_step

end module time_step_mod
