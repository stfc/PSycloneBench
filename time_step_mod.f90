module time_step_mod
  implicit none

  private

  public invoke_time_step

contains

  subroutine invoke_time_step(istp, ssha, ssha_u, ssha_v, &
                              sshn_t, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)
    use timing_mod
    use field_mod
    use grid_mod
    use model_mod,       only: rdt
    use momentum_mod,    only: momentum_u_code, momentum_v_code
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
    integer :: M, N, idxt

    M  = ssha%grid%simulation_domain%xstop
    N  = ssha%grid%simulation_domain%ystop

    ! In the general case we have to reason about whether or not the
    ! domain has PBCs and what sort of offset convention the kernels
    ! use. However, this is a middle layer specific to NEMOLite2D and
    ! therefore we know that we have no periodic BCs and are using a
    ! NE stagger
    !txstart = 2 ! grid%simulation_domain%xstart
    !tystart = 2 ! grid%simulation_domain%ystart

    !uxstart = 2 ! grid%simulation_domain%xstart
    !uxstop  = M - 1
    !uystart = 2 ! grid%simulation_domain%ystart
    !uystop  = N

    !vxstart = 2 ! grid%simulation_domain%xstart
    !vxstop  = M
    !vystart = 2 ! grid%simulation_domain%ystart
    !vystop  = N - 1

    !uwhole_xstart = 1 ! uxstart - NBOUNDARY
    !uwhole_xstop  = M ! uxstop  + NBOUNDARY
    !uwhole_ystart = 1 ! uystart - NBOUNDARY
    !uwhole_ystop  = N+1 ! uystop  + NBOUNDARY

    !vwhole_xstart = 1 ! vxstart - NBOUNDARY
    !vwhole_xstop  = M+1 ! vxstop  + NBOUNDARY
    !vwhole_ystart = 1 ! vystart - NBOUNDARY
    !vwhole_ystop  = N ! vystop  + NBOUNDARY

    call timer_start('Continuity',idxt)

!    do jj = ssha%internal%ystart, ssha%internal%ystop, 1
!      do ji = ssha%internal%xstart, ssha%internal%xstop, 1
    do jj = 2, N, 1
      do ji = 2, M, 1

        call continuity_code(ji, jj,                             &
                             ssha%data, sshn_t%data,             &
                             sshn_u%data, sshn_v%data,           &
                             hu%data, hv%data, un%data, vn%data, &
                             rdt, ssha%grid%area_t)
      end do
    end do

    call timer_stop(idxt)

    call timer_start('Momentum',idxt)

!    do jj = ua%internal%ystart, ua%internal%ystop, 1
!      do ji = ua%internal%xstart, ua%internal%xstop, 1
!dir$ safe_address
    do jj = 2, N, 1
! SIMD
      do ji = 2, M-1, 1

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
 
!dir$ safe_address
    do jj = 2, N-1, 1
! SIMD
      do ji = 2, M, 1

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

    call timer_stop(idxt)

    ! Apply open and solid boundary conditions

    call timer_start('BCs', idxt)

!    DO jj = ssha%internal%ystart, ssha%internal%ystop 
!       DO ji = ssha%internal%xstart, ssha%internal%xstop 
    DO jj = 2, N
! SIMD
       DO ji = 2, M
          call bc_ssh_code(ji, jj, &
                           istp, ssha%data, ssha%grid%tmask)
       END DO
    END DO


!    do jj = uwhole_ystart, uwhole_ystop, 1
!       do ji = uwhole_xstart, uwhole_xstop, 1
!dir$ safe_address
    do jj = 1, N+1, 1
       do ji = 1, M, 1
          call bc_solid_u_code(ji, jj, ua%data, ua%grid%tmask)
       end do
    end do

!    DO jj = va%whole%ystart, va%whole%ystop, 1 
!       DO ji = va%whole%xstart, va%whole%xstop, 1
!    do jj = vwhole_ystart, vwhole_ystop, 1
!       do ji = vwhole_xstart, vwhole_xstop, 1
!dir$ safe_address
    do jj = 1, N, 1
       do ji = 1, M+1, 1
          call bc_solid_v_code(ji,jj,va%data,va%grid%tmask)
      end do
    end do

!    DO jj = uwhole_ystart, uwhole_ystop, 1
!       DO ji = uwhole_xstart, uwhole_xstop, 1
!dir$ safe_address
    DO jj = 1, N+1, 1
       DO ji = 1, M, 1
          call bc_flather_u_code(ji,jj, &
                                 ua%data, hu%data, sshn_u%data, &
                                 ua%grid%tmask)
       END DO
    END DO

!    DO jj = va%whole%ystart, va%whole%ystop, 1 
!       DO ji = va%whole%xstart, va%whole%xstop, 1
!     DO jj = vwhole_ystart, vwhole_ystop, 1
!       DO ji = vwhole_xstart, vwhole_xstop, 1
!dir$ safe_address
     DO jj = 1, N, 1
       DO ji = 1, M+1, 1
          call bc_flather_v_code(ji,jj, &
                                 va%data, hv%data, sshn_v%data, &
                                 va%grid%tmask)
       END DO
    END DO

    call timer_stop(idxt)

    ! Time update of fields

    call timer_start('Next', idxt)

    call copy_field(ua, un)
    call copy_field(va, vn)
    call copy_field(ssha, sshn_t)

    do jj = 2, N, 1
      do ji = 2, M-1, 1

         call next_sshu_code(ji, jj, sshn_u%data, sshn_t%data, &
                            sshn_u%grid%tmask,                 &
                            sshn_u%grid%area_t, sshn_u%grid%area_u)
      end do
    end do

    do jj = 2, N-1, 1
      do ji = 2, M, 1

        call next_sshv_code(ji, jj,                   &
                            sshn_v%data, sshn_t%data, &
                            sshn_v%grid%tmask,        &
                            sshn_v%grid%area_t, sshn_v%grid%area_v)
      end do
    end do

    call timer_stop(idxt)

  end subroutine invoke_time_step

end module time_step_mod
