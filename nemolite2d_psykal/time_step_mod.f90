module time_step_mod
  implicit none

  private

  public invoke_time_step

contains

  subroutine invoke_time_step(istp, ssha, ssha_u, ssha_v, &
                              sshn_t, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)
    use kind_params_mod
    use timing_mod
    use field_mod
    use grid_mod
    use model_mod,       only: rdt
    use physical_params_mod, only: g, omega, d2r
    use momentum_mod,    only: momentum_v_code
    use momentum_mod,    only: momentum_u_code
    use continuity_mod,  only: continuity_code
    use time_update_mod, only: next_sshu_code, next_sshv_code
    use boundary_conditions_mod
    implicit none
    integer,         intent(in)    :: istp
    type(r2d_field), intent(inout) :: un, vn, sshn_t, sshn_u, sshn_v
    type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
    type(r2d_field), intent(in)    :: hu, hv, ht
    ! Locals
    integer :: ji, jj

    !M  = ssha%grid%simulation_domain%xstop
    !N  = ssha%grid%simulation_domain%ystop

    !call timer_start('Continuity',idxt)

!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none), &
!$OMP private(ji,jj) shared(ssha, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, rdt)
    do jj = ssha%internal%ystart, ssha%internal%ystop, 1
      do ji = ssha%internal%xstart, ssha%internal%xstop, 1
!    do jj = 2, N, 1
!      do ji = 2, M, 1

        call continuity_code(ji, jj,                             &
                             ssha%data, sshn_t%data,             &
                             sshn_u%data, sshn_v%data,           &
                             hu%data, hv%data, un%data, vn%data, &
                             rdt, sshn_t%grid%area_t)
      end do
    end do
!$OMP END PARALLEL DO

    !call timer_stop(idxt)

    !call timer_start('Momentum',idxt)

!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none) private(ji,jj) &
!$OMP shared(ua, un, vn, hu, hv, ht, ssha_u, sshn_t, sshn_u, sshn_v)
    do jj = ua%internal%ystart, ua%internal%ystop, 1
      do ji = ua%internal%xstart, ua%internal%xstop, 1
!    do jj = 2, N, 1
!      do ji = 2, M-1, 1

        call momentum_u_code(ji, jj, &
                             ua%data, un%data, vn%data, &
                             hu%data, hv%data, ht%data, &
                             ssha_u%data, sshn_t%data,  &
                             sshn_u%data, sshn_v%data,  &
                             un%grid%tmask,  &
                             un%grid%dx_u,   &
                             un%grid%dx_v,   &
                             un%grid%dx_t,   &
                             un%grid%dy_u,   &
                             un%grid%dy_t,   &
                             un%grid%area_u, &
                             un%grid%gphiu)

      end do
    end do
!$OMP END PARALLEL DO

!    do jj = 2, N-1, 1
!      do ji = 2, M, 1
!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none) private(ji,jj) &
!$OMP shared(va, un, vn, hu, hv, ht, ssha_v, sshn_t, sshn_u, sshn_v)
    do jj = va%internal%ystart, va%internal%ystop, 1
      do ji = va%internal%xstart, va%internal%xstop, 1

        call momentum_v_code(ji, jj, &
                             va%data, un%data, vn%data, &
                             hu%data, hv%data, ht%data, &
                             ssha_v%data, sshn_t%data,  &
                             sshn_u%data, sshn_v%data,  &
                             vn%grid%tmask,    &
                             vn%grid%dx_v,     &
                             vn%grid%dx_t,     &
                             vn%grid%dy_u,     &
                             vn%grid%dy_v,     &
                             vn%grid%dy_t,     &
                             vn%grid%area_v,   &
                             vn%grid%gphiv)

      end do
    end do
!$OMP END PARALLEL DO

    !call timer_stop(idxt)

    ! Apply open and solid boundary conditions

    !call timer_start('BCs', idxt)

!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none) private(ji,jj) &
!$OMP shared(istp, ssha, sshn_t)
    DO jj = ssha%internal%ystart, ssha%internal%ystop 
       DO ji = ssha%internal%xstart, ssha%internal%xstop 
!    DO jj = 2, N
!       DO ji = 2, M
          call bc_ssh_code(ji, jj, &
                           istp, ssha%data, sshn_t%grid%tmask)
       END DO
    END DO
!$OMP END PARALLEL DO


!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none) private(ji,jj) &
!$OMP shared(ua)
    do jj = ua%whole%ystart, ua%whole%ystop, 1
       do ji = ua%whole%xstart, ua%whole%xstop, 1
!    do jj = 1, N+1, 1
!       do ji = 1, M, 1
          call bc_solid_u_code(ji, jj, &
                               ua%data, ua%grid%tmask)
       end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none) private(ji,jj) &
!$OMP shared(va)
    DO jj = va%whole%ystart, va%whole%ystop, 1 
       DO ji = va%whole%xstart, va%whole%xstop, 1
!    do jj = 1, N, 1
!       do ji = 1, M+1, 1
          call bc_solid_v_code(ji,jj, &
                               va%data, va%grid%tmask)
      end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none) private(ji,jj) &
!$OMP shared(ua, hu, sshn_u)
    DO jj = ua%whole%ystart, ua%whole%ystop, 1
       DO ji = ua%whole%xstart, ua%whole%xstop, 1
!    DO jj = 1, N+1, 1
!       DO ji = 1, M, 1
          call bc_flather_u_code(ji,jj, &
                                 ua%data, hu%data, sshn_u%data, &
                                 sshn_u%grid%tmask)
       END DO
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none) private(ji,jj) &
!$OMP shared(va, hv, sshn_v)
    DO jj = va%whole%ystart, va%whole%ystop, 1 
       DO ji = va%whole%xstart, va%whole%xstop, 1
!     DO jj = 1, N, 1
!       DO ji = 1, M+1, 1
          call bc_flather_v_code(ji,jj,                         &
                                 va%data, hv%data, sshn_v%data, &
                                 sshn_v%grid%tmask)
       END DO
    END DO
!$OMP END PARALLEL DO

    !call timer_stop(idxt)

    ! Time update of fields

    !call timer_start('Next', idxt)

    call copy_field(ua, un)
    call copy_field(va, vn)
    call copy_field(ssha, sshn_t)

!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none) private(ji,jj) &
!$OMP shared(sshn_u, sshn_t)
    DO jj = sshn_u%internal%ystart, sshn_u%internal%ystop 
       DO ji = sshn_u%internal%xstart, sshn_u%internal%xstop 
!    do jj = 2, N, 1
!      do ji = 2, M-1, 1

         call next_sshu_code(ji, jj, sshn_u%data, sshn_t%data, &
                             sshn_t%grid%tmask,                 &
                             sshn_t%grid%area_t, sshn_t%grid%area_u)
      end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(RUNTIME), default(none) private(ji,jj) &
!$OMP shared(sshn_v, sshn_t)
    DO jj = sshn_v%internal%ystart, sshn_v%internal%ystop 
       DO ji = sshn_v%internal%xstart, sshn_v%internal%xstop 
!    do jj = 2, N-1, 1
!      do ji = 2, M, 1

        call next_sshv_code(ji, jj,                   &
                            sshn_v%data, sshn_t%data, &
                            sshn_t%grid%tmask,        &
                            sshn_t%grid%area_t, sshn_t%grid%area_v)
      end do
    end do
!$OMP END PARALLEL DO

    !call timer_stop(idxt)

  end subroutine invoke_time_step

end module time_step_mod
