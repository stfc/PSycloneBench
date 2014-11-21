module time_step_mod
  implicit none

  private

  public invoke_time_step

contains

  subroutine invoke_time_step(istp, ssha, ssha_u, ssha_v, &
                              sshn_t, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)
    use field_mod
    use model_mod, only: rdt
    use momentum_mod, only: momentum_u_code, momentum_v_code
    use continuity_mod, only: continuity_code
    use time_update_mod, only: next_sshu_code, next_sshv_code
    use boundary_conditions_mod
    implicit none
    integer,         intent(in)    :: istp
    type(r2d_field), intent(inout) :: un, vn, sshn_t, sshn_u, sshn_v
    type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
    type(r2d_field), intent(in)    :: hu, hv, ht
    ! Locals
    integer :: ji, jj

    do jj = ssha%internal%ystart, ssha%internal%ystop, 1
      do ji = ssha%internal%xstart, ssha%internal%xstop, 1

        call continuity_code(ji, jj,                             &
                             ssha%data, sshn_t%data,             &
                             sshn_u%data, sshn_v%data,           &
                             hu%data, hv%data, un%data, vn%data, &
                             rdt, ssha%grid%area_t)
      end do
    end do

    do jj = ua%internal%ystart, ua%internal%ystop, 1
      do ji = ua%internal%xstart, ua%internal%xstop, 1

        call momentum_u_code(ji, jj, &
                             ua%data, un%data, vn%data, &
                             hu%data, hv%data, ht%data, &
                             ssha_u%data, sshn_t%data,  &
                             sshn_u%data, sshn_v%data,  &
                             ua%grid%tmask,  &
                             ua%grid%dx_u,   &
                             ua%grid%dx_v,   &
                             ua%grid%dx_t,   &
                             ua%grid%dy_u,   &
                             ua%grid%dy_t,   &
                             ua%grid%area_u, &
                             ua%grid%gphiu)
      end do
    end do
 
    do jj = va%internal%ystart, va%internal%ystop, 1
      do ji = va%internal%xstart, va%internal%xstop, 1

        call momentum_v_code(ji, jj, &
                             va%data, un%data, vn%data, &
                             hu%data, hv%data, ht%data, &
                             ssha_v%data, sshn_t%data,      &
                             sshn_u%data, sshn_v%data,      &
                             va%grid%tmask, va%grid%dx_v,   &
                             va%grid%dx_t, &
                             va%grid%dy_u, va%grid%dy_v,    &
                             va%grid%dy_t,     &
                             va%grid%area_v, va%grid%gphiv)

      end do
    end do

    ! Apply open and solid boundary conditions

    DO jj = ssha%internal%ystart, ssha%internal%ystop
       DO ji = ssha%internal%xstart, ssha%internal%xstop
          call bc_ssh_code(ji, jj, &
                           istp, ssha%data, ssha%grid%tmask)
       END DO
    END DO


    do jj = ua%whole%ystart, ua%whole%ystop, 1
       do ji = ua%whole%xstart, ua%whole%xstop, 1
          call bc_solid_u_code(ji, jj, ua%data, ua%grid%tmask)
       end do
    end do

    do jj = va%whole%ystart, va%whole%ystop, 1
       do ji = va%whole%xstart, va%whole%xstop, 1
          call bc_solid_v_code(ji,jj,va%data,va%grid%tmask)
      end do
    end do

    DO jj = ua%whole%ystart, ua%whole%ystop, 1
       DO ji = ua%whole%xstart, ua%whole%xstop, 1
          call bc_flather_u_code(ji,jj, &
                                 ua%data, hu%data, sshn_u%data, &
                                 ua%grid%tmask)
       END DO
    END DO

    DO jj = va%whole%ystart, va%whole%ystop, 1
       DO ji = va%whole%xstart, va%whole%xstop, 1
          call bc_flather_v_code(ji,jj, &
                                 va%data, hv%data, sshn_v%data, &
                                 va%grid%tmask)
       END DO
    END DO

    ! Time update of fields
    call copy_field(ua, un)
    call copy_field(va, vn)
    call copy_field(ssha, sshn_t)

    do jj = sshn_u%internal%ystart, sshn_u%internal%ystop, 1
      do ji = sshn_u%internal%xstart, sshn_u%internal%xstop, 1

         call next_sshu_code(ji, jj, sshn_u%data, sshn_t%data, &
                            sshn_u%grid%tmask,                 &
                            sshn_u%grid%area_t, sshn_u%grid%area_u)
      end do
    end do

    do jj = sshn_v%internal%ystart, sshn_v%internal%ystop, 1
      do ji = sshn_v%internal%xstart, sshn_v%internal%xstop, 1

        call next_sshv_code(ji, jj,                   &
                            sshn_v%data, sshn_t%data, &
                            sshn_v%grid%tmask,        &
                            sshn_v%grid%area_t, sshn_v%grid%area_v)
      end do
    end do

  end subroutine invoke_time_step

end module time_step_mod
