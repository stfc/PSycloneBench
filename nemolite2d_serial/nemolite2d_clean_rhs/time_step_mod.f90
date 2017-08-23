module time_step_mod
  implicit none

  private

  public invoke_time_step

contains

  subroutine invoke_time_step(istp, ssha, ssha_u, ssha_v, &
                              sshn_t, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)
    use kind_params_mod
    use dl_timer
    use field_mod
    use grid_mod
    use model_mod,       only: rdt, cbfr, visc
    use physical_params_mod, only: g, omega, d2r
!    use momentum_mod,    only: momentum_v_code
!    use momentum_mod,    only: momentum_u_code
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
    ! Locals for momentum
    REAL(wp) :: u_e, u_w, v_n, v_s
    real(wp) :: v_nc, v_sc
    real(wp) :: depe, depw, deps, depn
    real(wp) :: hpg, adv, cor, vis
    real(wp) :: dudx_e, dudx_w, dudy_s, dudy_n
    real(wp) :: uu_e, uu_n, uu_s, uu_w
    real(wp) :: u_ec, u_wc, vv_e, vv_n, vv_s, vv_w
    real(wp) :: dvdx_e, dvdx_w, dvdy_n, dvdy_s

    ! end locals for momentum-u

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

    call timer_start(idxt, label='Continuity')

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

    call timer_start(idxt, label='Momentum')

!    do jj = ua%internal%ystart, ua%internal%ystop, 1
!      do ji = ua%internal%xstart, ua%internal%xstop, 1
!dir$ safe_address
    do jj = 2, N, 1
! SIMD
!dir$ vector always
      do ji = 2, M-1, 1

!!$        call momentum_u_code(ji, jj, &
!!$                             ua%data, un%data, vn%data, &
!!$                             hu%data, hv%data, ht%data, &
!!$                             ssha_u%data, sshn_t%data,  &
!!$                             sshn_u%data, sshn_v%data,  &
!!$                             un%grid%tmask,  &
!!$                             un%grid%dx_u,   &
!!$                             un%grid%dx_v,   &
!!$                             un%grid%dx_t,   &
!!$                             un%grid%dy_u,   &
!!$                             un%grid%dy_t,   &
!!$                             un%grid%area_u, &
!!$                             un%grid%gphiu)

    IF(un%grid%tmask(ji,jj) + un%grid%tmask(ji+1,jj) <= 0)  CYCLE   !jump over non-computational domain
    IF(un%grid%tmask(ji,jj) <= 0 .OR. un%grid%tmask(ji+1,jj) <= 0)  CYCLE !jump over boundary u

    u_e  = 0.5 * (un%data(ji,jj) + un%data(ji+1,jj)) * un%grid%dy_t(ji+1,jj)   !add length scale.
    depe = ht%data(ji+1,jj) + sshn_t%data(ji+1,jj)

    u_w  = 0.5 * (un%data(ji,jj) + un%data(ji-1,jj)) * un%grid%dy_t(ji,jj)     !add length scale
    depw = ht%data(ji,jj) + sshn_t%data(ji,jj)

    v_sc = 0.5_wp * (vn%data(ji,jj-1) + vn%data(ji+1,jj-1))
    v_s  = 0.5_wp * v_sc * (un%grid%dx_v(ji,jj-1) + un%grid%dx_v(ji+1,jj-1))
    deps = 0.5_wp * (hv%data(ji,jj-1) + sshn_v%data(ji,jj-1) + hv%data(ji+1,jj-1) + &
                     sshn_v%data(ji+1,jj-1))

    v_nc = 0.5_wp * (vn%data(ji,jj) + vn%data(ji+1,jj))
    v_n  = 0.5_wp * v_nc * (un%grid%dx_v(ji,jj) + un%grid%dx_v(ji+1,jj))
    depn = 0.5_wp * (hv%data(ji,jj) + sshn_v%data(ji,jj) + hv%data(ji+1,jj) + &
                     sshn_v%data(ji+1,jj))

    ! -advection (currently first order upwind)
    uu_w = (0.5_wp - SIGN(0.5_wp, u_w)) * un%data(ji,jj)              + & 
         & (0.5_wp + SIGN(0.5_wp, u_w)) * un%data(ji-1,jj) 
    uu_e = (0.5_wp + SIGN(0.5_wp, u_e)) * un%data(ji,jj)              + & 
         & (0.5_wp - SIGN(0.5_wp, u_e)) * un%data(ji+1,jj) 

    IF(un%grid%tmask(ji,jj-1) <=0 .OR. un%grid%tmask(ji+1,jj-1) <= 0) THEN   
       uu_s = (0.5_wp - SIGN(0.5_wp, v_s)) * un%data(ji,jj)   
    ELSE
       uu_s = (0.5_wp - SIGN(0.5_wp, v_s)) * un%data(ji,jj)              + & 
            & (0.5_wp + SIGN(0.5_wp, v_s)) * un%data(ji,jj-1) 
    END If

    IF(un%grid%tmask(ji,jj+1) <=0 .OR. un%grid%tmask(ji+1,jj+1) <= 0) THEN   
       uu_n = (0.5_wp + SIGN(0.5_wp, v_n)) * un%data(ji,jj)
    ELSE
       uu_n = (0.5_wp + SIGN(0.5_wp, v_n)) * un%data(ji,jj)              + & 
            & (0.5_wp - SIGN(0.5_wp, v_n)) * un%data(ji,jj+1)
    END IF

    adv = uu_w * u_w * depw - uu_e * u_e * depe + &
          uu_s * v_s * deps - uu_n * v_n * depn
    !end kernel u adv 

    ! -viscosity

    !kernel  u vis 
    dudx_e = (un%data(ji+1,jj) - un%data(ji,  jj)) / un%grid%dx_t(ji+1,jj) * &
             (ht%data(ji+1,jj) + sshn_t%data(ji+1,jj))
    dudx_w = (un%data(ji,  jj) - un%data(ji-1,jj)) / un%grid%dx_t(ji,  jj) * &
             (ht%data(ji,  jj) + sshn_t%data(ji,  jj))
    IF(un%grid%tmask(ji,jj-1) <=0 .OR. un%grid%tmask(ji+1,jj-1) <= 0) THEN   
       dudy_s = 0.0_wp !slip boundary
    ELSE
       dudy_s = (un%data(ji,jj) - un%data(ji,jj-1)) / (un%grid%dy_u(ji,jj) + un%grid%dy_u(ji,jj-1)) * &
            & (hu%data(ji,jj) + sshn_u%data(ji,jj) + hu%data(ji,jj-1) + sshn_u%data(ji,jj-1))
    END IF

    IF(un%grid%tmask(ji,jj+1) <= 0 .OR. un%grid%tmask(ji+1,jj+1) <= 0) THEN   
       dudy_n = 0.0_wp ! slip boundary
    ELSE
       dudy_n = (un%data(ji,jj+1) - un%data(ji,jj)) / (un%grid%dy_u(ji,jj) + un%grid%dy_u(ji,jj+1)) * &
            & (hu%data(ji,jj) + sshn_u%data(ji,jj) + hu%data(ji,jj+1) + sshn_u%data(ji,jj+1))
    END If

    vis = (dudx_e - dudx_w ) * un%grid%dy_u(ji,jj)  + &
         & (dudy_n - dudy_s ) * un%grid%dx_u(ji,jj) * 0.5_wp  
    vis = visc * vis   !visc will be an array visc(1:jpijglou) 
    !for variable viscosity, such as turbulent viscosity
    !End  kernel u vis 

    ! -Coriolis' force (can be implemented implicitly)
    !kernel cor 
    cor = 0.5_wp * (2._wp * omega * SIN(un%grid%gphiu(ji,jj) * d2r) * (v_sc + v_nc)) * &
         & un%grid%area_u(ji,jj) * (hu%data(ji,jj) + sshn_u%data(ji,jj))
    !end kernel cor 

    ! -pressure gradient
    !start kernel hpg 
    hpg = -g * (hu%data(ji,jj) + sshn_u%data(ji,jj)) * un%grid%dy_u(ji,jj) * &
           (sshn_t%data(ji+1,jj) - sshn_t%data(ji,jj))
    !end kernel hpg 
    ! -linear bottom friction (implemented implicitly.
    !kernel ua calculation 
    ua%data(ji,jj) = (un%data(ji,jj) * (hu%data(ji,jj) + sshn_u%data(ji,jj)) + rdt * &
                 (adv + vis + cor + hpg) / un%grid%area_u(ji,jj)) / &
                (hu%data(ji,jj) + ssha_u%data(ji,jj)) / (1.0_wp + cbfr * rdt) 

      end do
    end do
 
!dir$ safe_address
    do jj = 2, N-1, 1
! SIMD
!dir$ vector always
      do ji = 2, M, 1

!!$        call momentum_v_code(ji, jj, &
!!$                             va%data, un%data, vn%data, &
!!$                             hu%data, hv%data, ht%data, &
!!$                             ssha_v%data, sshn_t%data,  &
!!$                             sshn_u%data, sshn_v%data,  &
!!$                             vn%grid%tmask,    &
!!$                             vn%grid%dx_v,     &
!!$                             vn%grid%dx_t,     &
!!$                             vn%grid%dy_u,     &
!!$                             vn%grid%dy_v,     &
!!$                             vn%grid%dy_t,     &
!!$                             vn%grid%area_v,   &
!!$                             vn%grid%gphiv)

    IF(vn%grid%tmask(ji,jj) + vn%grid%tmask(ji+1,jj) <= 0)  cycle !jump over non-computatinal domain
    IF(vn%grid%tmask(ji,jj) <= 0 .OR. vn%grid%tmask(ji,jj+1) <= 0) cycle !jump over v boundary cells

    ! kernel v adv 
    v_n  = 0.5 * (vn%data(ji,jj) + vn%data(ji,jj+1)) * vn%grid%dx_t(ji,jj+1)  !add length scale.
    depn = ht%data(ji,jj+1) + sshn_t%data(ji,jj+1)

    v_s  = 0.5 * (vn%data(ji,jj) + vn%data(ji,jj-1)) * vn%grid%dx_t(ji,jj)    !add length scale
    deps = ht%data(ji,jj) + sshn_t%data(ji,jj)

    u_wc = 0.5_wp * (un%data(ji-1,jj) + un%data(ji-1,jj+1))
    u_w  = 0.5_wp * u_wc * (vn%grid%dy_u(ji-1,jj) + vn%grid%dy_u(ji-1,jj+1))
    depw = 0.50_wp * (hu%data(ji-1,jj) + sshn_u%data(ji-1,jj) + &
                      hu%data(ji-1,jj+1) + sshn_u%data(ji-1,jj+1))

    u_ec = 0.5_wp * (un%data(ji,jj) + un%data(ji,jj+1))
    u_e  = 0.5_wp * u_ec * (vn%grid%dy_u(ji,jj) + vn%grid%dy_u(ji,jj+1))
    depe = 0.50_wp * (hu%data(ji,jj) + sshn_u%data(ji,jj) + &
                      hu%data(ji,jj+1) + sshn_u%data(ji,jj+1))

    ! -advection (currently first order upwind)
    vv_s = (0.5_wp - SIGN(0.5_wp, v_s)) * vn%data(ji,jj)     + & 
         & (0.5_wp + SIGN(0.5_wp, v_s)) * vn%data(ji,jj-1) 
    vv_n = (0.5_wp + SIGN(0.5_wp, v_n)) * vn%data(ji,jj)     + & 
         & (0.5_wp - SIGN(0.5_wp, v_n)) * vn%data(ji,jj+1) 

    IF(vn%grid%tmask(ji-1,jj) <= 0 .OR. vn%grid%tmask(ji-1,jj+1) <= 0) THEN   
       vv_w = (0.5_wp - SIGN(0.5_wp, u_w)) * vn%data(ji,jj)  
    ELSE
       vv_w = (0.5_wp - SIGN(0.5_wp, u_w)) * vn%data(ji,jj)    + & 
            & (0.5_wp + SIGN(0.5_wp, u_w)) * vn%data(ji-1,jj) 
    END If

    IF(vn%grid%tmask(ji+1,jj) <= 0 .OR. vn%grid%tmask(ji+1,jj+1) <= 0) THEN
       vv_e = (0.5_wp + SIGN(0.5_wp, u_e)) * vn%data(ji,jj)
    ELSE
       vv_e = (0.5_wp + SIGN(0.5_wp, u_e)) * vn%data(ji,jj)  + & 
              (0.5_wp - SIGN(0.5_wp, u_e)) * vn%data(ji+1,jj)
    END IF

    adv = vv_w * u_w * depw - vv_e * u_e * depe + &
          vv_s * v_s * deps - vv_n * v_n * depn

    !end kernel v adv 

    ! -viscosity

    
    !kernel v dis 
    dvdy_n = (vn%data(ji,jj+1) - vn%data(ji,  jj)) / vn%grid%dy_t(ji,jj+1) * &
                          (ht%data(ji,jj+1) + sshn_t%data(ji,jj+1))
    dvdy_s = (vn%data(ji,  jj) - vn%data(ji,jj-1)) / vn%grid%dy_t(ji,  jj) * &
                          (ht%data(ji,  jj) + sshn_t%data(ji,  jj))

    IF(vn%grid%tmask(ji-1,jj) <= 0 .OR. vn%grid%tmask(ji-1,jj+1) <= 0) THEN
       dvdx_w = 0.0_wp !slip boundary
    ELSE
       dvdx_w = (vn%data(ji,jj) - vn%data(ji-1,jj)) / &
                (vn%grid%dx_v(ji,jj) + vn%grid%dx_v(ji-1,jj)) * &
                (hv%data(ji,jj) + sshn_v%data(ji,jj) + hv%data(ji-1,jj) + sshn_v%data(ji-1,jj))
    END IF

    IF(vn%grid%tmask(ji+1,jj) <= 0 .OR. vn%grid%tmask(ji+1,jj+1) <= 0) THEN
       dvdx_e = 0.0_wp ! slip boundary
    ELSE
       dvdx_e = (vn%data(ji+1,jj) - vn%data(ji,jj)) / (vn%grid%dx_v(ji,jj) + vn%grid%dx_v(ji+1,jj)) * &
                  (hv%data(ji,jj) + sshn_v%data(ji,jj) + hv%data(ji+1,jj) + sshn_v%data(ji+1,jj))
    END If

    vis = (dvdy_n - dvdy_s ) * vn%grid%dx_v(ji,jj)  + &
          (dvdx_e - dvdx_w ) * vn%grid%dy_v(ji,jj) * 0.5_wp  

    vis = visc * vis   !visc will be a array visc(1:jpijglou) 
    !for variable viscosity, such as turbulent viscosity
    !end kernel v dis 

    ! -Coriolis' force (can be implemented implicitly)
    !kernel v cor 
    cor = -0.5_wp*(2._wp * omega * SIN(vn%grid%gphiv(ji,jj) * d2r) * (u_ec + u_wc)) * &
               vn%grid%area_v(ji,jj) * (hv%data(ji,jj) + sshn_v%data(ji,jj))
    !end kernel v cor 

    ! -pressure gradient
    !kernel v hpg 
    hpg = -g * (hv%data(ji,jj) + sshn_v%data(ji,jj)) * vn%grid%dx_v(ji,jj) * &
           (sshn_t%data(ji,jj+1) - sshn_t%data(ji,jj))
    !kernel v hpg 

    ! -linear bottom friction (implemented implicitly.
    !kernel ua calculation 
    va%data(ji,jj) = (vn%data(ji,jj) * (hv%data(ji,jj) + sshn_v%data(ji,jj)) + &
                 rdt * (adv + vis + cor + hpg) / vn%grid%area_v(ji,jj) ) / &
                 ((hv%data(ji,jj) + ssha_v%data(ji,jj))) / (1.0_wp + cbfr * rdt) 

      end do
    end do

    call timer_stop(idxt)

    ! Apply open and solid boundary conditions

    call timer_start(idxt, label='BCs')

!    DO jj = ssha%internal%ystart, ssha%internal%ystop 
!       DO ji = ssha%internal%xstart, ssha%internal%xstop 
    DO jj = 2, N
! SIMD
       DO ji = 2, M
          call bc_ssh_code(ji, jj, &
                           istp, ssha%data, sshn_t%grid%tmask)
       END DO
    END DO


!    do jj = uwhole_ystart, uwhole_ystop, 1
!       do ji = uwhole_xstart, uwhole_xstop, 1
!dir$ safe_address
    do jj = 1, N+1, 1
       do ji = 1, M, 1
          call bc_solid_u_code(ji, jj, &
                               ua%data, va%grid%tmask)
       end do
    end do

!    DO jj = va%whole%ystart, va%whole%ystop, 1 
!       DO ji = va%whole%xstart, va%whole%xstop, 1
!    do jj = vwhole_ystart, vwhole_ystop, 1
!       do ji = vwhole_xstart, vwhole_xstop, 1
!dir$ safe_address
    do jj = 1, N, 1
       do ji = 1, M+1, 1
          call bc_solid_v_code(ji,jj, &
                               va%data, ua%grid%tmask)
      end do
    end do

!    DO jj = uwhole_ystart, uwhole_ystop, 1
!       DO ji = uwhole_xstart, uwhole_xstop, 1
!dir$ safe_address
    DO jj = 1, N+1, 1
       DO ji = 1, M, 1
          call bc_flather_u_code(ji,jj, &
                                 ua%data, hu%data, sshn_u%data, &
                                 sshn_u%grid%tmask)
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
                                 sshn_v%grid%tmask)
       END DO
    END DO

    call timer_stop(idxt)

    ! Time update of fields

    call timer_start(idxt, label='Next')

    call copy_field(ua, un)
    call copy_field(va, vn)
    call copy_field(ssha, sshn_t)

    do jj = 2, N, 1
      do ji = 2, M-1, 1

         call next_sshu_code(ji, jj, sshn_u%data, sshn_t%data, &
                            sshn_t%grid%tmask,                 &
                            sshn_t%grid%area_t, sshn_t%grid%area_u)
      end do
    end do

    do jj = 2, N-1, 1
      do ji = 2, M, 1

        call next_sshv_code(ji, jj,                   &
                            sshn_v%data, sshn_t%data, &
                            sshn_t%grid%tmask,        &
                            sshn_t%grid%area_t, sshn_t%grid%area_v)
      end do
    end do

    call timer_stop(idxt)

  end subroutine invoke_time_step

end module time_step_mod
