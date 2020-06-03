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
    ! use momentum_mod,    only: momentum_v_code
    ! use momentum_mod,    only: momentum_u_code
    ! use continuity_mod,  only: continuity_code
    ! use time_update_mod, only: next_sshu_code, next_sshv_code
    use boundary_conditions_mod
    implicit none
    real(go_wp),        intent(in)    :: istp
    type(r2d_field), intent(inout) :: un, vn, sshn_t, sshn_u, sshn_v
    type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
    type(r2d_field), intent(in)    :: hu, hv, ht
    ! Locals
    integer :: ji, jj, jiu, jiv
    integer :: M, N
    integer :: cont_timer, mom_timer, bc_timer, next_timer
    ! Locals for momentum
    REAL(go_wp) :: u_e, u_w, v_n, v_s
    real(go_wp) :: v_nc, v_sc
    real(go_wp) :: depe, depw, deps, depn
    real(go_wp) :: hpg, adv, cor, vis
    real(go_wp) :: dudx_e, dudx_w, dudy_s, dudy_n
    real(go_wp) :: uu_e, uu_n, uu_s, uu_w
    real(go_wp) :: u_ec, u_wc, vv_e, vv_n, vv_s, vv_w
    real(go_wp) :: dvdx_e, dvdx_w, dvdy_n, dvdy_s
    real(go_wp) :: rtmp1, rtmp2, rtmp3, rtmp4
    ! end locals for momentum
    ! Locals for BCs
    real(go_wp) :: amp_tide, omega_tide, rtime
    ! Locals for loop bounds
    integer :: xstart, xstop, ystart, ystop, uxstop, vystop
    integer :: whole_xstart, whole_xstop, whole_ystart, whole_ystop
    integer :: uwhole_xstop, vwhole_ystop

    ! In the general case we have to reason about whether or not the
    ! domain has PBCs and what sort of offset convention the kernels
    ! use. However, this is a middle layer specific to NEMOLite2D and
    ! therefore we know that we have no periodic BCs and are using a
    ! NE stagger
    xstart = ssha%grid%subdomain%internal%xstart
    ystart = ssha%grid%subdomain%internal%ystart
    xstop = ssha%grid%subdomain%internal%xstop
    ystop = ssha%grid%subdomain%internal%ystop

    uxstop = xstop - 1
    vystop = ystop - 1

    whole_xstart = xstart - NBOUNDARY
    whole_xstop  = xstop  + NBOUNDARY
    whole_ystart = ystart - NBOUNDARY
    whole_ystop  = ystop  + NBOUNDARY

    uwhole_xstop = uxstop + NBOUNDARY
    vwhole_ystop = vystop + NBOUNDARY

    call timer_start(cont_timer, label='Continuity')

    do jj = ystart, ystop
      do ji = xstart, xstop 

        ! call continuity_code(ji, jj,                             &
        !                      ssha%data, sshn_t%data,             &
        !                      sshn_u%data, sshn_v%data,           &
        !                      hu%data, hv%data, un%data, vn%data, &
        !                      sshn_t%grid%area_t)
         rtmp1 = (sshn_u%data(ji  ,jj ) + hu%data(ji  ,jj  ))*un%data(ji  ,jj)
         rtmp2 = (sshn_u%data(ji-1,jj ) + hu%data(ji-1,jj  ))*un%data(ji-1,jj)
         rtmp3 = (sshn_v%data(ji ,jj ) + hv%data(ji  ,jj  ))*vn%data(ji ,jj)
         rtmp4 = (sshn_v%data(ji ,jj-1) + hv%data(ji  ,jj-1))*vn%data(ji,jj-1)

         ssha%data(ji,jj) = sshn_t%data(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * &
                       rdt / sshn_t%grid%area_t(ji,jj)
      end do
    end do

    call timer_stop(cont_timer)

    call timer_start(mom_timer, label='Momentum')

    !dir$ safe_address
    do jj = ystart, ystop
      ! SIMD
      !dir$ vector always
      do ji = xstart, uxstop 

        ! call momentum_u_code(ji, jj, &
        !                      ua%data, un%data, vn%data, &
        !                      hu%data, hv%data, ht%data, &
        !                      ssha_u%data, sshn_t%data,  &
        !                      sshn_u%data, sshn_v%data,  &
        !                      un%grid%tmask,  &
        !                      un%grid%dx_u,   &
        !                      un%grid%dx_v,   &
        !                      un%grid%dx_t,   &
        !                      un%grid%dy_u,   &
        !                      un%grid%dy_t,   &
        !                      un%grid%area_u, &
        !                      un%grid%gphiu)
        
        ! Jump over non-computational domain
        IF(un%grid%tmask(ji,jj) + un%grid%tmask(ji+1,jj) <= 0)  CYCLE
        ! Jump over boundary u
        IF(un%grid%tmask(ji,jj) <= 0 .OR. un%grid%tmask(ji+1,jj) <= 0)  CYCLE

        ! Add length scale.
        u_e  = 0.5 * (un%data(ji,jj) + un%data(ji+1,jj)) * un%grid%dy_t(ji+1,jj)   
        depe = ht%data(ji+1,jj) + sshn_t%data(ji+1,jj)

        ! Add length scale
        u_w  = 0.5 * (un%data(ji,jj) + un%data(ji-1,jj)) * un%grid%dy_t(ji,jj)
        depw = ht%data(ji,jj) + sshn_t%data(ji,jj)

        v_sc = 0.5_go_wp * (vn%data(ji,jj-1) + vn%data(ji+1,jj-1))
        v_s  = 0.5_go_wp * v_sc * (un%grid%dx_v(ji,jj-1) + un%grid%dx_v(ji+1,jj-1))
        deps = 0.5_go_wp * (hv%data(ji,jj-1) + sshn_v%data(ji,jj-1) + hv%data(ji+1,jj-1) + &
                         sshn_v%data(ji+1,jj-1))

        v_nc = 0.5_go_wp * (vn%data(ji,jj) + vn%data(ji+1,jj))
        v_n  = 0.5_go_wp * v_nc * (un%grid%dx_v(ji,jj) + un%grid%dx_v(ji+1,jj))
        depn = 0.5_go_wp * (hv%data(ji,jj) + sshn_v%data(ji,jj) + hv%data(ji+1,jj) + &
                         sshn_v%data(ji+1,jj))

        ! -advection (currently first order upwind)
        uu_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * un%data(ji,jj)              + & 
             & (0.5_go_wp + SIGN(0.5_go_wp, u_w)) * un%data(ji-1,jj) 
        uu_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * un%data(ji,jj)              + & 
             & (0.5_go_wp - SIGN(0.5_go_wp, u_e)) * un%data(ji+1,jj) 

        IF(un%grid%tmask(ji,jj-1) <=0 .OR. un%grid%tmask(ji+1,jj-1) <= 0) THEN   
           uu_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * un%data(ji,jj)   
        ELSE
           uu_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * un%data(ji,jj)              + & 
                & (0.5_go_wp + SIGN(0.5_go_wp, v_s)) * un%data(ji,jj-1) 
        END If

        IF(un%grid%tmask(ji,jj+1) <=0 .OR. un%grid%tmask(ji+1,jj+1) <= 0) THEN   
           uu_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * un%data(ji,jj)
        ELSE
           uu_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * un%data(ji,jj)              + & 
                & (0.5_go_wp - SIGN(0.5_go_wp, v_n)) * un%data(ji,jj+1)
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
           dudy_s = 0.0_go_wp !slip boundary
        ELSE
           dudy_s = (un%data(ji,jj) - un%data(ji,jj-1)) / (un%grid%dy_u(ji,jj) + &
               un%grid%dy_u(ji,jj-1)) * (hu%data(ji,jj) + sshn_u%data(ji,jj) + &
               hu%data(ji,jj-1) + sshn_u%data(ji,jj-1))
        END IF

        IF(un%grid%tmask(ji,jj+1) <= 0 .OR. un%grid%tmask(ji+1,jj+1) <= 0) THEN   
           dudy_n = 0.0_go_wp ! slip boundary
        ELSE
           dudy_n = (un%data(ji,jj+1) - un%data(ji,jj)) / (un%grid%dy_u(ji,jj) + &
               un%grid%dy_u(ji,jj+1)) * (hu%data(ji,jj) + sshn_u%data(ji,jj) + &
               hu%data(ji,jj+1) + sshn_u%data(ji,jj+1))
        END If

        vis = (dudx_e - dudx_w ) * un%grid%dy_u(ji,jj)  + &
             & (dudy_n - dudy_s ) * un%grid%dx_u(ji,jj) * 0.5_go_wp  
        vis = visc * vis   !visc will be an array visc(1:jpijglou) 
        !for variable viscosity, such as turbulent viscosity
        !End  kernel u vis 

        ! -Coriolis' force (can be implemented implicitly)
        !kernel cor 
        cor = 0.5_go_wp * (2._go_wp * omega * SIN(un%grid%gphiu(ji,jj) * d2r) * (v_sc + v_nc)) * &
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
                    (hu%data(ji,jj) + ssha_u%data(ji,jj)) / (1.0_go_wp + cbfr * rdt) 

      end do
    end do
 
    !dir$ safe_address
    do jj = ystart, vystop
      ! SIMD
      !dir$ vector always
      do ji = xstart, xstop

        ! call momentum_v_code(ji, jj, &
        !                      va%data, un%data, vn%data, &
        !                      hu%data, hv%data, ht%data, &
        !                      ssha_v%data, sshn_t%data,  &
        !                      sshn_u%data, sshn_v%data,  &
        !                      vn%grid%tmask,    &
        !                      vn%grid%dx_v,     &
        !                      vn%grid%dx_t,     &
        !                      vn%grid%dy_u,     &
        !                      vn%grid%dy_v,     &
        !                      vn%grid%dy_t,     &
        !                      vn%grid%area_v,   &
        !                      vn%grid%gphiv)

        ! Jump over non-computatinal domain
        IF(vn%grid%tmask(ji,jj) + vn%grid%tmask(ji+1,jj) <= 0)  cycle
        ! Jump over v boundary cells
        IF(vn%grid%tmask(ji,jj) <= 0 .OR. vn%grid%tmask(ji,jj+1) <= 0) cycle

        ! kernel v adv
        ! add length scale.
        v_n  = 0.5 * (vn%data(ji,jj) + vn%data(ji,jj+1)) * vn%grid%dx_t(ji,jj+1)
        depn = ht%data(ji,jj+1) + sshn_t%data(ji,jj+1)

        ! Add length scale
        v_s  = 0.5 * (vn%data(ji,jj) + vn%data(ji,jj-1)) * vn%grid%dx_t(ji,jj)
        deps = ht%data(ji,jj) + sshn_t%data(ji,jj)

        u_wc = 0.5_go_wp * (un%data(ji-1,jj) + un%data(ji-1,jj+1))
        u_w  = 0.5_go_wp * u_wc * (vn%grid%dy_u(ji-1,jj) + vn%grid%dy_u(ji-1,jj+1))
        depw = 0.50_go_wp * (hu%data(ji-1,jj) + sshn_u%data(ji-1,jj) + &
                          hu%data(ji-1,jj+1) + sshn_u%data(ji-1,jj+1))

        u_ec = 0.5_go_wp * (un%data(ji,jj) + un%data(ji,jj+1))
        u_e  = 0.5_go_wp * u_ec * (vn%grid%dy_u(ji,jj) + vn%grid%dy_u(ji,jj+1))
        depe = 0.50_go_wp * (hu%data(ji,jj) + sshn_u%data(ji,jj) + &
                          hu%data(ji,jj+1) + sshn_u%data(ji,jj+1))

        ! -advection (currently first order upwind)
        vv_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * vn%data(ji,jj)     + & 
             & (0.5_go_wp + SIGN(0.5_go_wp, v_s)) * vn%data(ji,jj-1) 
        vv_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * vn%data(ji,jj)     + & 
             & (0.5_go_wp - SIGN(0.5_go_wp, v_n)) * vn%data(ji,jj+1) 

        IF(vn%grid%tmask(ji-1,jj) <= 0 .OR. vn%grid%tmask(ji-1,jj+1) <= 0) THEN   
           vv_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * vn%data(ji,jj)  
        ELSE
           vv_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * vn%data(ji,jj)    + & 
                & (0.5_go_wp + SIGN(0.5_go_wp, u_w)) * vn%data(ji-1,jj) 
        END If

        IF(vn%grid%tmask(ji+1,jj) <= 0 .OR. vn%grid%tmask(ji+1,jj+1) <= 0) THEN
           vv_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * vn%data(ji,jj)
        ELSE
           vv_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * vn%data(ji,jj)  + & 
                  (0.5_go_wp - SIGN(0.5_go_wp, u_e)) * vn%data(ji+1,jj)
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
           dvdx_w = 0.0_go_wp !slip boundary
        ELSE
           dvdx_w = (vn%data(ji,jj) - vn%data(ji-1,jj)) / &
                    (vn%grid%dx_v(ji,jj) + vn%grid%dx_v(ji-1,jj)) * &
                    (hv%data(ji,jj) + sshn_v%data(ji,jj) + hv%data(ji-1,jj) + &
                    sshn_v%data(ji-1,jj))
        END IF

        IF(vn%grid%tmask(ji+1,jj) <= 0 .OR. vn%grid%tmask(ji+1,jj+1) <= 0) THEN
           dvdx_e = 0.0_go_wp ! slip boundary
        ELSE
           dvdx_e = (vn%data(ji+1,jj) - vn%data(ji,jj)) / (vn%grid%dx_v(ji,jj) + &
               vn%grid%dx_v(ji+1,jj)) * (hv%data(ji,jj) + sshn_v%data(ji,jj) + &
               hv%data(ji+1,jj) + sshn_v%data(ji+1,jj))
        END If

        vis = (dvdy_n - dvdy_s ) * vn%grid%dx_v(ji,jj)  + &
              (dvdx_e - dvdx_w ) * vn%grid%dy_v(ji,jj) * 0.5_go_wp  

        vis = visc * vis   !visc will be a array visc(1:jpijglou) 
        !for variable viscosity, such as turbulent viscosity
        !end kernel v dis 

        ! -Coriolis' force (can be implemented implicitly)
        !kernel v cor 
        cor = -0.5_go_wp*(2._go_wp * omega * SIN(vn%grid%gphiv(ji,jj) * d2r) * &
            (u_ec + u_wc)) * vn%grid%area_v(ji,jj) * (hv%data(ji,jj) + sshn_v%data(ji,jj))
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
                     ((hv%data(ji,jj) + ssha_v%data(ji,jj))) / (1.0_go_wp + cbfr * rdt) 

      end do
    end do

    call timer_stop(mom_timer)

    ! Apply open and solid boundary conditions

    call timer_start(bc_timer, label='BCs')

    DO jj = ystart, ystop
      ! SIMD
      DO ji = xstart, xstop
        ! call bc_ssh_code(ji, jj, istp, ssha%data, sshn_t%grid%tmask)
        amp_tide   = 0.2_go_wp
        omega_tide = 2.0_go_wp * 3.14159_go_wp / (12.42_go_wp * 3600._go_wp)
        rtime = istp * rdt

        if(sshn_t%grid%tmask(ji,jj) <= 0) cycle

        IF     (sshn_t%grid%tmask(ji,jj-1) < 0) THEN
           ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
        ELSE IF(sshn_t%grid%tmask(ji,jj+1) < 0) THEN
           ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
        ELSE IF(sshn_t%grid%tmask(ji+1,jj) < 0) THEN
           ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
        ELSE IF(sshn_t%grid%tmask(ji-1,jj) < 0) THEN
           ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
        END IF

       END DO
    END DO


    !dir$ safe_address
    do jj = whole_ystart, whole_ystop
       do ji = whole_xstart, uwhole_xstop
         ! call bc_solid_u_code(ji, jj, ua%data, va%grid%tmask)

         if(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji+1,jj) == 0)then
            ua%data(ji,jj) = 0._go_wp
         end if
       end do
    end do

    !dir$ safe_address
    do jj = whole_ystart, vwhole_ystop
       do ji = whole_xstart, whole_xstop
         ! call bc_solid_v_code(ji,jj, va%data, ua%grid%tmask)
         if(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji,jj+1) == 0)then
           va%data(ji,jj) = 0._go_wp
          end if
      end do
    end do

    !dir$ safe_address
    do jj = whole_ystart, whole_ystop
       do ji = whole_xstart, uwhole_xstop
          ! call bc_flather_u_code(ji,jj, &
          !                        ua%data, hu%data, sshn_u%data, &
          !                        sshn_u%grid%tmask)
          ! Check whether this point lies within the domain
          if(sshn_t%grid%tmask(ji,jj) + sshn_t%grid%tmask(ji+1,jj) <= -1) cycle

          if(sshn_t%grid%tmask(ji,jj) < 0) then
             jiu = ji + 1
             ua%data(ji,jj) = ua%data(jiu,jj) + sqrt(g/hu%data(ji,jj))* &
                  (sshn_u%data(ji,jj) - sshn_u%data(jiu,jj))
          else if(sshn_t%grid%tmask(ji+1,jj )< 0) then
             jiu = ji - 1 
             ua%data(ji,jj) = ua%data(jiu,jj) + sqrt(g/hu%data(ji,jj)) * &
                  (sshn_u%data(ji,jj) - sshn_u%data(jiu,jj))
          end if
       END DO
    END DO

    !dir$ safe_address
    do jj = whole_ystart, vwhole_ystop
       do ji = whole_xstart, whole_xstop
          ! call bc_flather_v_code(ji,jj, &
          !                        va%data, hv%data, sshn_v%data, &
          !                        sshn_v%grid%tmask)
          IF(sshn_t%grid%tmask(ji,jj) + sshn_t%grid%tmask(ji,jj+1) <= -1) cycle
    
          IF(sshn_t%grid%tmask(ji,jj) < 0) THEN
             jiv = jj + 1
             va%data(ji,jj) = va%data(ji,jiv) + SQRT(g/hv%data(ji,jj)) * &
                  (sshn_v%data(ji,jj) - sshn_v%data(ji,jiv))
          ELSE IF(sshn_t%grid%tmask(ji,jj+1) < 0) THEN
             jiv = jj - 1 
             va%data(ji,jj) = va%data(ji,jiv) + SQRT(g/hv%data(ji,jj)) * &
                  (sshn_v%data(ji,jj) - sshn_v%data(ji,jiv))
          END IF
       END DO
    END DO

    call timer_stop(bc_timer)

    ! Time update of fields

    call timer_start(next_timer, label='Next')

    ! call copy_field(ua, un)
    ! call copy_field(va, vn)
    ! call copy_field(ssha, sshn_t)
    un%data = ua%data
    vn%data = va%data
    sshn_t%data = ssha%data

    !dir$ safe_address
    do jj = ystart, ystop
      !dir$ vector always
      do ji = xstart, uxstop
        ! call next_sshu_code(ji, jj, sshn_u%data, sshn_t%data, &
        !                     sshn_t%grid%tmask,                 &
        !                     sshn_t%grid%area_t, sshn_t%grid%area_u)

         if(sshn_t%grid%tmask(ji,jj) + &
            sshn_t%grid%tmask(ji+1,jj) <= 0) cycle !jump over non-computational domain

         IF(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji+1,jj) > 0) THEN
            rtmp1 = sshn_t%grid%area_t(ji,jj) * sshn_t%data(ji,jj) + &
                 sshn_t%grid%area_t(ji+1,jj) * sshn_t%data(ji+1,jj)
            sshn_u%data(ji,jj) = 0.5_go_wp * rtmp1 / sshn_t%grid%area_u(ji,jj) 
         ELSE IF(sshn_t%grid%tmask(ji,jj) <= 0) THEN
            sshn_u%data(ji,jj) = sshn_t%data(ji+1,jj)
         ELSE IF(sshn_t%grid%tmask(ji+1,jj) <= 0) THEN
            sshn_u%data(ji,jj) = sshn_t%data(ji,jj)
         END IF

      end do
    end do

    !dir$ safe_address
    do jj = ystart, vystop
      !dir$ vector always
      do ji = xstart, xstop
        ! call next_sshv_code(ji, jj,                   &
        !                     sshn_v%data, sshn_t%data, &
        !                     sshn_t%grid%tmask,        &
        !                     sshn_t%grid%area_t, sshn_t%grid%area_v)
 
         if(sshn_t%grid%tmask(ji,jj) + &
            sshn_t%grid%tmask(ji,jj+1) <= 0)  cycle !jump over non-computational domain
         if(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji,jj+1) > 0) then
            rtmp1 = sshn_t%grid%area_t(ji,jj)*sshn_t%data(ji,jj) + &
                 sshn_t%grid%area_t(ji,jj+1) * sshn_t%data(ji,jj+1)
            sshn_v%data(ji,jj) = 0.5_go_wp * rtmp1 / sshn_t%grid%area_v(ji,jj) 
         else if(sshn_t%grid%tmask(ji,jj) <= 0) then
            sshn_v%data(ji,jj) = sshn_t%data(ji,jj+1)
         else if(sshn_t%grid%tmask(ji,jj+1) <= 0) then
            sshn_v%data(ji,jj) = sshn_t%data(ji,jj)
         end if
      end do
    end do

    call timer_stop(next_timer)

  end subroutine invoke_time_step

end module time_step_mod
