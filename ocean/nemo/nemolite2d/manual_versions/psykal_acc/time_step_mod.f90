module time_step_mod
  use kind_params_mod, only: wp
  implicit none

  private

  public invoke_time_step

contains

  subroutine invoke_time_step(istp, ssha, ssha_u, ssha_v, &
                              sshn_t, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)
    !use dl_timer
    use field_mod
!    use grid_mod
!    use model_mod,       only: rdt, cbfr, visc
!    use physical_params_mod, only: g, omega, d2r
!    use momentum_mod,    only: momentum_v_code
!    use momentum_mod,    only: momentum_u_code
!    use continuity_mod,  only: continuity_code
!    use time_update_mod, only: next_sshu_code, next_sshv_code
!    use boundary_conditions_mod
    implicit none
    integer,         intent(in)    :: istp
    type(r2d_field), intent(inout) :: un, vn, sshn_t, sshn_u, sshn_v
    type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
    type(r2d_field), intent(in)    :: hu, hv, ht

    call invoke_time_step1(istp, ssha, ssha_u, ssha_v, &
                              sshn_t, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)

    call invoke_time_step2(istp, ssha, ssha_u, ssha_v, &
                              sshn_t, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)

  end subroutine invoke_time_step

  subroutine invoke_time_step1(istp, ssha, ssha_u, ssha_v, &
                              sshn_t, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)
    use field_mod
    use grid_mod
    use model_mod,       only: rdt
    implicit none
    integer,         intent(in)    :: istp
    type(r2d_field), intent(inout) :: un, vn, sshn_t, sshn_u, sshn_v
    type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
    type(r2d_field), intent(in)    :: hu, hv, ht
    ! Locals
    integer :: ji, jj, jiu, jiv
    integer :: M, N
    integer :: cont_timer, mom_timer, bc_timer, next_timer

    M  = ssha%grid%simulation_domain%xstop
    N  = ssha%grid%simulation_domain%ystop

    !call timer_start(cont_timer, label='Continuity')

!    do jj = ssha%internal%ystart, ssha%internal%ystop, 1
!      do ji = ssha%internal%xstart, ssha%internal%xstop, 1
! Copy data to GPU. We use pcopyin so that if the data is already
    ! on the GPU then that copy is used.
    if(.not. ssha%data_on_device .or. &
       .not. sshn_t%data_on_device .or. &
       .not. sshn_u%data_on_device .or. &
       .not. sshn_v%data_on_device .or. &
       .not. hu%data_on_device .or. &
       .not. hv%data_on_device .or. &
       .not. un%data_on_device .or. &
       .not. vn%data_on_device)then
!$acc enter data if(.not. ssha%data_on_device) &
!$acc copyin(sshn_t, sshn_t%data, sshn_u, sshn_u%data, sshn_v, sshn_v%data, &
!$acc        ssha, ssha%data,  &
!$acc        hu, hu%data, hv, hv%data, ht, ht%data, &
!$acc        un, un%data, vn, vn%data, &
!$acc        sshn_t%grid, sshn_t%grid%tmask, sshn_t%grid%area_t)
       ssha%data_on_device = .true.
       sshn_t%data_on_device = .true.
       sshn_u%data_on_device = .true.
       sshn_v%data_on_device = .true.
       hu%data_on_device = .true.
       hv%data_on_device = .true.
       un%data_on_device = .true.
       vn%data_on_device = .true.
    end if
    
!$acc parallel default(present)
!$acc loop collapse(2)
    do jj = 2, N, 1
      do ji = 2, M, 1

        call continuity_code(ji, jj,                             &
                             ssha%data, sshn_t%data,             &
                             sshn_u%data, sshn_v%data,           &
                             hu%data, hv%data, un%data, vn%data, &
                             rdt, sshn_t%grid%area_t)
      end do
    end do
!$acc end parallel

    !call timer_stop(cont_timer)

  end subroutine invoke_time_step1

  subroutine invoke_time_step2(istp, ssha, ssha_u, ssha_v, &
                                sshn_t, sshn_u, sshn_v, &
                                hu, hv, ht, ua, va, un, vn)
    use field_mod
    use grid_mod
    use model_mod,       only: rdt, cbfr, visc
    use momentum_mod,    only: momentum_u_code
    implicit none
    integer,         intent(in)    :: istp
    type(r2d_field), intent(inout) :: un, vn, sshn_t, sshn_u, sshn_v
    type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
    type(r2d_field), intent(in)    :: hu, hv, ht
    ! Locals
    integer :: ji, jj, jiu, jiv
    integer :: M, N
    integer :: cont_timer, mom_timer, bc_timer, next_timer

    M  = ua%grid%simulation_domain%xstop
    N  = ua%grid%simulation_domain%ystop

    !call timer_start(cont_timer, label='Momentum-u')

!    do jj = ssha%internal%ystart, ssha%internal%ystop, 1
!      do ji = ssha%internal%xstart, ssha%internal%xstop, 1
! Copy data to GPU. We use pcopyin so that if the data is already
! on the GPU then that copy is used.
    if(.not. ua%data_on_device .or. &
       .not. un%data_on_device .or. &
       .not. vn%data_on_device .or. &
       .not. hu%data_on_device .or. &
       .not. hv%data_on_device .or. &
       .not. ht%data_on_device .or. &
       .not. ssha_u%data_on_device .or. &
       .not. sshn_t%data_on_device .or. &
       .not. sshn_u%data_on_device .or. &
       .not. sshn_v%data_on_device)then
!$acc enter data &
!$acc copyin(sshn_t, sshn_t%data, sshn_u, sshn_u%data, sshn_v, sshn_v%data, &
!$acc        ssha_u, ssha_u%data,  &
!$acc        hu, hu%data, hv, hv%data, ht, ht%data, &
!$acc        ua, ua%data, un, un%data, vn, vn%data, &
!$acc        sshn_t%grid, sshn_t%grid%tmask, sshn_t%grid%area_t, &
!$acc        sshn_t%grid%area_u, sshn_t%grid%dx_u, sshn_t%grid%dx_v,   &
!$acc        sshn_t%grid%dx_t, sshn_t%grid%dy_u, sshn_t%grid%dy_t,   &
!$acc        sshn_t%grid%gphiu, rdt, cbfr, visc)
       ua%data_on_device = .true.
       un%data_on_device = .true.
       vn%data_on_device = .true.
       hu%data_on_device = .true.
       hv%data_on_device = .true.
       ht%data_on_device = .true.
       ssha_u%data_on_device = .true.
       sshn_t%data_on_device = .true.
       sshn_u%data_on_device = .true.
       sshn_v%data_on_device = .true.
     end if
     
!$acc parallel default(present)
!$acc loop collapse(2)
    do jj = 2, N, 1
       do ji = 2, M-1, 1

        call momentum_u_code(ji, jj, &
                             ua%data, un%data, vn%data, &
                             hu%data, hv%data, ht%data, &
                             ssha_u%data, sshn_t%data,  &
                             sshn_u%data, sshn_v%data,  &
                             sshn_t%grid%tmask,  &
                             sshn_t%grid%dx_u,   &
                             sshn_t%grid%dx_v,   &
                             sshn_t%grid%dx_t,   &
                             sshn_t%grid%dy_u,   &
                             sshn_t%grid%dy_t,   &
                             sshn_t%grid%area_u, &
                             sshn_t%grid%gphiu)
       end do
    end do
!$acc end parallel

    !call timer_stop(cont_timer)

  end subroutine invoke_time_step2

  ! This routine has been 'module in-lined' to check that PGI compiler
  ! can cope with a !$acc routine so long as it's within the same
  ! module as the call site.
  subroutine continuity_code(ji, jj,                     &
                             ssha, sshn, sshn_u, sshn_v, &
                             hu, hv, un, vn, rdt, e12t)
    implicit none
!$acc routine seq
    integer,                  intent(in)  :: ji, jj
    real(wp),                 intent(in)  :: rdt
    real(wp), dimension(:,:), intent(in)  :: e12t
    real(wp), dimension(:,:), intent(out) :: ssha
    real(wp), dimension(:,:), intent(in)  :: sshn, sshn_u, sshn_v, &
                                             hu, hv, un, vn
    ! Locals
    real(wp) :: rtmp1, rtmp2, rtmp3, rtmp4

    rtmp1 = (sshn_u(ji  ,jj  ) + hu(ji  ,jj  )) * un(ji  ,jj  )
    rtmp2 = (sshn_u(ji-1,jj  ) + hu(ji-1,jj  )) * un(ji-1,jj  )
    rtmp3 = (sshn_v(ji  ,jj  ) + hv(ji  ,jj  )) * vn(ji  ,jj  )
    rtmp4 = (sshn_v(ji  ,jj-1) + hv(ji  ,jj-1)) * vn(ji  ,jj-1)

    ssha(ji,jj) = sshn(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * &
                    rdt / e12t(ji,jj)

  end subroutine continuity_code

end module time_step_mod
