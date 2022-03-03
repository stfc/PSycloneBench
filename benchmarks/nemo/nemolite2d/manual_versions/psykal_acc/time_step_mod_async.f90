module time_step_mod
  use kind_params_mod
  implicit none

  private

  public invoke_time_step

contains

  subroutine invoke_time_step(istp, ssha, ssha_u, ssha_v, &
                              sshn_t, sshn_u, sshn_v, &
                              hu, hv, ht, ua, va, un, vn)
    !> This routine is purely a wrapper that converts references to
    !! data arrays within derived types to straightforward data
    !! arrays in the called routine. This is to stop compilers
    !! from panicking at the sight of a '%' symbol!
    use field_mod
    use grid_mod
    implicit none
    integer,         intent(in)    :: istp
    type(r2d_field), intent(inout) :: un, vn, sshn_t, sshn_u, sshn_v
    type(r2d_field), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
    type(r2d_field), intent(inout) :: hu, hv, ht
    LOGICAL, save :: first_time=.true.
                     
    call invoke_time_step_arrays(istp,                                  &
                                 ua%grid%nx, ua%grid%ny,                &
                                 ua%grid%subdomain%internal%xstop,      &
                                 ua%grid%subdomain%internal%ystop,      &
                                 ua%grid%area_t,                        &
                                 ua%grid%area_u,                        &
                                 ua%grid%area_v,                        &
                                 ua%grid%dx_u,                          &
                                 ua%grid%dx_v,                          &
                                 ua%grid%dx_t,                          &
                                 ua%grid%dy_u,                          &
                                 ua%grid%dy_v,                          &
                                 ua%grid%dy_t,                          &
                                 ua%grid%gphiu,                         &
                                 ua%grid%gphiv,                         &
                                 ua%grid%tmask,                         &
                                 ssha%data, ssha_u%data, ssha_v%data,   &
                                 sshn_t%data, sshn_u%data, sshn_v%data, &
                                 hu%data, hv%data, ht%data,             &
                                 ua%data, va%data, un%data, vn%data, un%grid)

    if (first_time) then
        first_time = .false.
        ua%data_on_device     = .TRUE.
        va%data_on_device     = .TRUE.
        un%data_on_device     = .TRUE.
        vn%data_on_device     = .TRUE.
        ssha%data_on_device   = .TRUE.
        ssha_u%data_on_device = .TRUE.
        ssha_v%data_on_device = .TRUE.
        sshn_t%data_on_device = .TRUE.
        sshn_u%data_on_device = .TRUE.
        sshn_v%data_on_device = .TRUE.
        hu%data_on_device     = .TRUE.
        hv%data_on_device     = .TRUE.
        ht%data_on_device     = .TRUE.

        ! Specify device data retrieving methods
        ssha%read_from_device_f => read_openacc
        sshn_t%read_from_device_f => read_openacc
        sshn_u%read_from_device_f => read_openacc
        sshn_v%read_from_device_f => read_openacc
        un%read_from_device_f => read_openacc
        vn%read_from_device_f => read_openacc
        ua%read_from_device_f => read_openacc
        va%read_from_device_f => read_openacc
    endif

  end subroutine invoke_time_step

  !===============================================

  subroutine read_openacc(from, to, startx, starty, nx, ny, blocking)
    ! Function that specify how to retrieve the device data 'from' to a host
    ! location 'to', this function will be called by the infrastructure
    ! whenever the data is needed on the host.
    use iso_c_binding, only: c_ptr
    type(c_ptr), intent(in) :: from
    real(go_wp), dimension(:,:), intent(inout), target :: to
    integer, intent(in) :: startx, starty, nx, ny
    logical, intent(in) :: blocking

    ! Currently non-blocking reads are only requested by halo_exchanges which
    ! this manual version doesn't have, so it is safe to ignore. In the future
    ! we could use the 'async' clause if non-blocking synchronisations are
    ! needed.

    !$acc update host(to)
  end subroutine read_openacc

  subroutine invoke_time_step_arrays(istp, nx, ny, M, N,         &
                                     area_t, area_u, area_v,     &
                                     dx_u, dx_v, dx_t,           &
                                     dy_u, dy_v, dy_t,           &
                                     gphiu, gphiv,               &
                                     tmask,                      &
                                     ssha, ssha_u, ssha_v,       &
                                     sshn_t, sshn_u, sshn_v,     &
                                     hu, hv, ht, ua, va, un, vn, grid)
    use dl_timer
    use model_mod,           only: rdt, cbfr, visc
    use physical_params_mod, only: g, omega, d2r
    use grid_mod
    implicit none
    !> The current time step
    integer,         intent(in) :: istp
    !> The allocated dimensions of all fields
    integer,         intent(in) :: nx, ny
    !> The loop bounds for the simulated domain
    integer,         intent(in) :: M, N
    real(go_wp), dimension(nx,ny), intent(inout) :: un, vn, sshn_t, sshn_u, sshn_v
    real(go_wp), dimension(nx,ny), intent(inout) :: ua, va, ssha, ssha_u, ssha_v
    real(go_wp), dimension(nx,ny), intent(in)    :: hu, hv, ht
    real(go_wp), dimension(nx,ny), intent(in)    :: area_t, area_u, area_v
    real(go_wp), dimension(nx,ny), intent(in)    :: dx_t, dx_u, dx_v
    real(go_wp), dimension(nx,ny), intent(in)    :: dy_t, dy_u, dy_v
    real(go_wp), dimension(nx,ny), intent(in)    :: gphiu, gphiv
    integer,  dimension(nx,ny), intent(in)    :: tmask
    type(grid_type), intent(in) :: grid
    ! Locals
    integer :: it, ji, jj, jiu, jiv
    integer :: idxt
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

    integer, dimension(:,:), pointer :: tmaskptr

!    call timer_start('Continuity',idxt)
    tmaskptr => grid%get_tmask()


! Copy data to GPU. We use pcopyin so that if the data is already
! on the GPU then that copy is used.
!$acc enter data                             &
!$acc pcopyin(tmask, area_t, area_u, area_v, &
!$acc        dx_t, dx_u, dx_v,               &
!$acc        dy_t, dy_u, dy_v, gphiu, gphiv, &
!$acc        sshn_t, sshn_u, sshn_v,         &
!$acc        ssha, ssha_u, ssha_v,           &
!$acc        ht, hu, hv, ua, va, un, vn, rdt, tmaskptr)

!$acc parallel loop present(ssha, sshn_t, sshn_u, sshn_v, &
!$acc                  hu, hv, un, vn, area_t, rdt) async(1)
    do jj = 2, N, 1
      do ji = 2, M, 1

        !call continuity_code(ji, jj,                   &
        !                     ssha, sshn_t,             &
        !                     sshn_u, sshn_v,           &
        !                     hu, hv, un, vn,           &
        !                     rdt, area_t)
        rtmp1 = (sshn_u(ji  ,jj ) + hu(ji  ,jj  ))*un(ji  ,jj)
        rtmp2 = (sshn_u(ji-1,jj ) + hu(ji-1,jj  ))*un(ji-1,jj)
        rtmp3 = (sshn_v(ji ,jj )  + hv(ji  ,jj  ))*vn(ji ,jj)
        rtmp4 = (sshn_v(ji ,jj-1) + hv(ji  ,jj-1))*vn(ji,jj-1)

        ssha(ji,jj) = sshn_t(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * &
                      rdt / area_t(ji,jj)
      end do
    end do
!$acc end parallel

!    call timer_stop(idxt)

!    call timer_start('Momentum',idxt)

!$acc parallel loop present(ua, un, vn, hu, hv, ht,         &
!$acc                  ssha_u, sshn_t, sshn_u, sshn_v, &
!$acc                  tmask, area_u, gphiu,           &
!$acc                  dx_u, dx_v, dx_t,               &
!$acc                  dy_u, dy_t) async(2)
    do jj = 2, N, 1
      do ji = 2, M-1, 1
!!$        call momentum_u_code(ji, jj, &
!!$                             ua, un, vn, &
!!$                             hu, hv, ht, &
!!$                             ssha_u, sshn_t,  &
!!$                             sshn_u, sshn_v,  &
!!$                             un%grid%tmask,  &
!!$                             un%grid%dx_u,   &
!!$                             un%grid%dx_v,   &
!!$                             un%grid%dx_t,   &
!!$                             un%grid%dy_u,   &
!!$                             un%grid%dy_t,   &
!!$                             un%grid%area_u, &
!!$                             un%grid%gphiu)

    IF(tmask(ji,jj) + tmask(ji+1,jj) <= 0)  CYCLE   !jump over non-computational domain
    IF(tmask(ji,jj) <= 0 .OR. tmask(ji+1,jj) <= 0)  CYCLE !jump over boundary u

    u_e  = 0.5 * (un(ji,jj) + un(ji+1,jj)) * dy_t(ji+1,jj)   !add length scale.
    depe = ht(ji+1,jj) + sshn_t(ji+1,jj)

    u_w  = 0.5 * (un(ji,jj) + un(ji-1,jj)) * dy_t(ji,jj)     !add length scale
    depw = ht(ji,jj) + sshn_t(ji,jj)

    v_sc = 0.5_go_wp * (vn(ji,jj-1) + vn(ji+1,jj-1))
    v_s  = 0.5_go_wp * v_sc * (dx_v(ji,jj-1) + dx_v(ji+1,jj-1))
    deps = 0.5_go_wp * (hv(ji,jj-1) + sshn_v(ji,jj-1) + hv(ji+1,jj-1) + &
                     sshn_v(ji+1,jj-1))

    v_nc = 0.5_go_wp * (vn(ji,jj) + vn(ji+1,jj))
    v_n  = 0.5_go_wp * v_nc * (dx_v(ji,jj) + dx_v(ji+1,jj))
    depn = 0.5_go_wp * (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + &
                     sshn_v(ji+1,jj))

    ! -advection (currently first order upwind)
    uu_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * un(ji,jj)              + & 
         & (0.5_go_wp + SIGN(0.5_go_wp, u_w)) * un(ji-1,jj) 
    uu_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * un(ji,jj)              + & 
         & (0.5_go_wp - SIGN(0.5_go_wp, u_e)) * un(ji+1,jj) 

    IF(tmask(ji,jj-1) <=0 .OR. tmask(ji+1,jj-1) <= 0) THEN   
       uu_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * un(ji,jj)   
    ELSE
       uu_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * un(ji,jj)              + & 
            & (0.5_go_wp + SIGN(0.5_go_wp, v_s)) * un(ji,jj-1) 
    END If

    IF(tmask(ji,jj+1) <=0 .OR. tmask(ji+1,jj+1) <= 0) THEN   
       uu_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * un(ji,jj)
    ELSE
       uu_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * un(ji,jj)              + & 
            & (0.5_go_wp - SIGN(0.5_go_wp, v_n)) * un(ji,jj+1)
    END IF

    adv = uu_w * u_w * depw - uu_e * u_e * depe + &
          uu_s * v_s * deps - uu_n * v_n * depn
    !end kernel u adv 

    ! -viscosity

    !kernel  u vis 
    dudx_e = (un(ji+1,jj) - un(ji,  jj)) / dx_t(ji+1,jj) * &
             (ht(ji+1,jj) + sshn_t(ji+1,jj))
    dudx_w = (un(ji,  jj) - un(ji-1,jj)) / dx_t(ji,  jj) * &
             (ht(ji,  jj) + sshn_t(ji,  jj))
    IF(tmask(ji,jj-1) <=0 .OR. tmask(ji+1,jj-1) <= 0) THEN   
       dudy_s = 0.0_go_wp !slip boundary
    ELSE
       dudy_s = (un(ji,jj) - un(ji,jj-1)) / (dy_u(ji,jj) + dy_u(ji,jj-1)) * &
            & (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj-1) + sshn_u(ji,jj-1))
    END IF

    IF(tmask(ji,jj+1) <= 0 .OR. tmask(ji+1,jj+1) <= 0) THEN   
       dudy_n = 0.0_go_wp ! slip boundary
    ELSE
       dudy_n = (un(ji,jj+1) - un(ji,jj)) / (dy_u(ji,jj) + dy_u(ji,jj+1)) * &
            & (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj+1) + sshn_u(ji,jj+1))
    END If

    vis = (dudx_e - dudx_w ) * dy_u(ji,jj)  + &
         & (dudy_n - dudy_s ) * dx_u(ji,jj) * 0.5_go_wp  
    vis = visc * vis   !visc will be an array visc(1:jpijglou) 
    !for variable viscosity, such as turbulent viscosity
    !End  kernel u vis 

    ! -Coriolis' force (can be implemented implicitly)
    !kernel cor 
    cor = 0.5_go_wp * (2._go_wp * omega * SIN(gphiu(ji,jj) * d2r) * (v_sc + v_nc)) * &
         & area_u(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj))
    !end kernel cor 

    ! -pressure gradient
    !start kernel hpg 
    hpg = -g * (hu(ji,jj) + sshn_u(ji,jj)) * dy_u(ji,jj) * &
           (sshn_t(ji+1,jj) - sshn_t(ji,jj))
    !end kernel hpg 
    ! -linear bottom friction (implemented implicitly.
    !kernel ua calculation 
    ua(ji,jj) = (un(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj)) + rdt * &
                 (adv + vis + cor + hpg) / area_u(ji,jj)) / &
                (hu(ji,jj) + ssha_u(ji,jj)) / (1.0_go_wp + cbfr * rdt) 

      end do
    end do
!$acc end parallel


!$acc parallel loop present(tmask, va, ssha_v,area_v,gphiv,dy_v,hv,sshn_v,dx_v, &
!$acc                  dy_t,vn,sshn_u,hu,dy_u,un,ht,sshn_t,dx_t) async(3)
    do jj = 2, N-1, 1
      do ji = 2, M, 1

!!$        call momentum_v_code(ji, jj, &
!!$                             va, un, vn, &
!!$                             hu, hv, ht, &
!!$                             ssha_v, sshn_t,  &
!!$                             sshn_u, sshn_v,  &
!!$                             vn%grid%tmask,    &
!!$                             vn%grid%dx_v,     &
!!$                             vn%grid%dx_t,     &
!!$                             vn%grid%dy_u,     &
!!$                             vn%grid%dy_v,     &
!!$                             vn%grid%dy_t,     &
!!$                             vn%grid%area_v,   &
!!$                             vn%grid%gphiv)

    IF(tmask(ji,jj) + tmask(ji+1,jj) <= 0)  cycle !jump over non-computatinal domain
    IF(tmask(ji,jj) <= 0 .OR. tmask(ji,jj+1) <= 0) cycle !jump over v boundary cells

    ! kernel v adv 
    v_n  = 0.5 * (vn(ji,jj) + vn(ji,jj+1)) * dx_t(ji,jj+1)  !add length scale.
    depn = ht(ji,jj+1) + sshn_t(ji,jj+1)

    v_s  = 0.5 * (vn(ji,jj) + vn(ji,jj-1)) * dx_t(ji,jj)    !add length scale
    deps = ht(ji,jj) + sshn_t(ji,jj)

    u_wc = 0.5_go_wp * (un(ji-1,jj) + un(ji-1,jj+1))
    u_w  = 0.5_go_wp * u_wc * (dy_u(ji-1,jj) + dy_u(ji-1,jj+1))
    depw = 0.50_go_wp * (hu(ji-1,jj) + sshn_u(ji-1,jj) + &
                      hu(ji-1,jj+1) + sshn_u(ji-1,jj+1))

    u_ec = 0.5_go_wp * (un(ji,jj) + un(ji,jj+1))
    u_e  = 0.5_go_wp * u_ec * (dy_u(ji,jj) + dy_u(ji,jj+1))
    depe = 0.50_go_wp * (hu(ji,jj) + sshn_u(ji,jj) + &
                      hu(ji,jj+1) + sshn_u(ji,jj+1))

    ! -advection (currently first order upwind)
    vv_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * vn(ji,jj)     + & 
         & (0.5_go_wp + SIGN(0.5_go_wp, v_s)) * vn(ji,jj-1) 
    vv_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * vn(ji,jj)     + & 
         & (0.5_go_wp - SIGN(0.5_go_wp, v_n)) * vn(ji,jj+1) 

    IF(tmask(ji-1,jj) <= 0 .OR. tmask(ji-1,jj+1) <= 0) THEN   
       vv_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * vn(ji,jj)  
    ELSE
       vv_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * vn(ji,jj)    + & 
            & (0.5_go_wp + SIGN(0.5_go_wp, u_w)) * vn(ji-1,jj) 
    END If

    IF(tmask(ji+1,jj) <= 0 .OR. tmask(ji+1,jj+1) <= 0) THEN
       vv_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * vn(ji,jj)
    ELSE
       vv_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * vn(ji,jj)  + & 
              (0.5_go_wp - SIGN(0.5_go_wp, u_e)) * vn(ji+1,jj)
    END IF

    adv = vv_w * u_w * depw - vv_e * u_e * depe + &
          vv_s * v_s * deps - vv_n * v_n * depn

    !end kernel v adv 

    ! -viscosity

    
    !kernel v dis 
    dvdy_n = (vn(ji,jj+1) - vn(ji,  jj)) / dy_t(ji,jj+1) * &
                          (ht(ji,jj+1) + sshn_t(ji,jj+1))
    dvdy_s = (vn(ji,  jj) - vn(ji,jj-1)) / dy_t(ji,  jj) * &
                          (ht(ji,  jj) + sshn_t(ji,  jj))

    IF(tmask(ji-1,jj) <= 0 .OR. tmask(ji-1,jj+1) <= 0) THEN
       dvdx_w = 0.0_go_wp !slip boundary
    ELSE
       dvdx_w = (vn(ji,jj) - vn(ji-1,jj)) / &
                (dx_v(ji,jj) + dx_v(ji-1,jj)) * &
                (hv(ji,jj) + sshn_v(ji,jj) + hv(ji-1,jj) + sshn_v(ji-1,jj))
    END IF

    IF(tmask(ji+1,jj) <= 0 .OR. tmask(ji+1,jj+1) <= 0) THEN
       dvdx_e = 0.0_go_wp ! slip boundary
    ELSE
       dvdx_e = (vn(ji+1,jj) - vn(ji,jj)) / (dx_v(ji,jj) + dx_v(ji+1,jj)) * &
                  (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + sshn_v(ji+1,jj))
    END If

    vis = (dvdy_n - dvdy_s ) * dx_v(ji,jj)  + &
          (dvdx_e - dvdx_w ) * dy_v(ji,jj) * 0.5_go_wp  

    vis = visc * vis   !visc will be a array visc(1:jpijglou) 
    !for variable viscosity, such as turbulent viscosity
    !end kernel v dis 

    ! -Coriolis' force (can be implemented implicitly)
    !kernel v cor 
    cor = -0.5_go_wp*(2._go_wp * omega * SIN(gphiv(ji,jj) * d2r) * (u_ec + u_wc)) * &
               area_v(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj))
    !end kernel v cor 

    ! -pressure gradient
    !kernel v hpg 
    hpg = -g * (hv(ji,jj) + sshn_v(ji,jj)) * dx_v(ji,jj) * &
           (sshn_t(ji,jj+1) - sshn_t(ji,jj))
    !kernel v hpg 

    ! -linear bottom friction (implemented implicitly.
    !kernel ua calculation 
    va(ji,jj) = (vn(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj)) + &
                 rdt * (adv + vis + cor + hpg) / area_v(ji,jj) ) / &
                 ((hv(ji,jj) + ssha_v(ji,jj))) / (1.0_go_wp + cbfr * rdt) 

      end do
    end do
!$acc end parallel
!$acc acc wait

!    call timer_stop(idxt)

! We block here as, strictly speaking, a thread could enter the loop
! below and begin writing to ssha while another is still in the very
! first loop and is also writing to ssha.

    ! Apply open and solid boundary conditions

!    call timer_start('BCs', idxt)
!$acc parallel present(tmask, ssha)
!$acc loop collapse(2)
    DO jj = 2, N
! SIMD
       DO ji = 2, M
!          call bc_ssh_code(ji, jj, &
!                           istp, ssha, sshn_t%grid%tmask)

          amp_tide   = 0.2_go_wp
          omega_tide = 2.0_go_wp * 3.14159_go_wp / (12.42_go_wp * 3600._go_wp)
          rtime = real(istp, go_wp) * rdt

          if(tmask(ji,jj) <= 0) cycle

          IF     (tmask(ji,jj-1) < 0) THEN
             ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
          ELSE IF(tmask(ji,jj+1) < 0) THEN
             ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
          ELSE IF(tmask(ji+1,jj) < 0) THEN
             ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
          ELSE IF(tmask(ji-1,jj) < 0) THEN
             ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
          END IF

       END DO
    END DO
!$acc end parallel
! This loop only writes to ssha and subsequent loop does not use
! this field therefore we need not block.

!$acc parallel present(tmaskptr, ua)
!$acc loop collapse(2)
    do jj = 1, N+1, 1
       do ji = 1, M, 1
          !call bc_solid_u_code(ji, jj, &
          !                     ua, tmaskptr)

          !> \todo It's more compiler-friendly to separately compare these two
          !! integer masks with zero but that's a kernel-level optimisation.
          if(tmask(ji,jj) * tmask(ji+1,jj) == 0)then
             ua(ji,jj) = 0._go_wp
          end if

       end do
    end do
!$acc end parallel
! This loop only writes to ua and subsequent loop does not use this field
! or the preceeding ssha so no need to block.

!$acc parallel present(tmask, va)
!$acc loop collapse(2)
    do jj = 1, N, 1
       do ji = 1, M+1, 1
!          call bc_solid_v_code(ji,jj, &
!                               va, ua%grid%tmask)
          if(tmask(ji,jj) * tmask(ji,jj+1) == 0)then
             va(ji,jj) = 0._go_wp
          end if

      end do
    end do
!$acc end parallel
! We must block here as next loop reads and writes ua.


!dir$ safe_address
! Although there is an apparent loop-carried dependency in j in this
! loop, the tmask is always set-up such that each trip *will* be
! independent
!$acc parallel present(tmask, va, hv, sshn_v)
!$acc loop collapse(2) independent
     DO jj = 1, N, 1
       DO ji = 1, M+1, 1
!          call bc_flather_v_code(ji,jj, &
!                                 va, hv, sshn_v, &
!                                 sshn_v%grid%tmask)
          IF(tmask(ji,jj) + tmask(ji,jj+1) <= -1) cycle
    
          IF(tmask(ji,jj) < 0) THEN
             jiv = jj + 1
             va(ji,jj) = va(ji,jiv) + SQRT(g/hv(ji,jj)) * &
                  (sshn_v(ji,jj) - sshn_v(ji,jiv))
          ELSE IF(tmask(ji,jj+1) < 0) THEN
             jiv = jj - 1 
             va(ji,jj) = va(ji,jiv) + SQRT(g/hv(ji,jj)) * &
                  (sshn_v(ji,jj) - sshn_v(ji,jiv))
          END IF

       END DO
    END DO
!$acc end parallel

!$acc parallel present(tmask, ua, hu, sshn_u)
!$acc loop collapse(2) independent
    DO jj = 1, N+1, 1
       DO ji = 1, M, 1
!          call bc_flather_u_code(ji,jj, &
!                                 ua, hu, sshn_u, &
!                                 sshn_u%grid%tmask)
          ! Check whether this point lies within the domain
          if(tmask(ji,jj) + tmask(ji+1,jj) <= -1) cycle

          if(tmask(ji,jj) < 0) then
             ! Read from column to the right (East) of us
             jiu = ji + 1
             ua(ji,jj) = ua(jiu,jj) + sqrt(g/hu(ji,jj))* &
                  (sshn_u(ji,jj) - sshn_u(jiu,jj))
          else if(tmask(ji+1,jj )< 0) then
             ! Read from column to the left of us
             jiu = ji - 1 
             ua(ji,jj) = ua(jiu,jj) + sqrt(g/hu(ji,jj)) * &
                  (sshn_u(ji,jj) - sshn_u(jiu,jj))
          end if
       END DO
    END DO
!$acc end parallel
! This loop only writes to ua and following loop does not use that field
! so no need to block here.

!    call timer_stop(idxt)

! We must block here since following loop reads from ua and va.

    ! Time update of fields

!    call timer_start('Next', idxt)

!> \todo It would be more efficient to merge these copies into the 
!! subsequent loops that update ssh{u,v}
!    call copy_field(ua, un)
!    call copy_field(va, vn)
!    call copy_field(ssha, sshn_t)
!$acc parallel present(un,vn,ua,va,ssha,sshn_t)
!$acc loop collapse(2)
     do jj = 1, N+1, 1
       do ji = 1, M+1, 1
          un(ji,jj) = ua(ji,jj)
          vn(ji,jj) = va(ji,jj)
          sshn_t(ji,jj) = ssha(ji,jj)
       end do
    end do
!$acc end parallel
! We have to block here since sshn_t is used in the following loop.
! We could avoid this by altering the following two loop nests to read from
! ssha instead of sshn_t.

!$acc parallel present(tmask, sshn_u, sshn_t, area_u, area_t)
!$acc loop collapse(2)
    do jj = 2, N, 1
      do ji = 2, M-1, 1

!         call next_sshu_code(ji, jj, sshn_u, sshn_t, &
!                            sshn_t%grid%tmask,                 &
!                            sshn_t%grid%area_t, sshn_t%grid%area_u)

         if(tmask(ji,jj) + &
            tmask(ji+1,jj) <= 0) cycle !jump over non-computational domain

         IF(tmask(ji,jj) * tmask(ji+1,jj) > 0) THEN
            rtmp1 = area_t(ji,jj) * sshn_t(ji,jj) + &
                 area_t(ji+1,jj) * sshn_t(ji+1,jj)
            sshn_u(ji,jj) = 0.5_go_wp * rtmp1 / area_u(ji,jj) 
         ELSE IF(tmask(ji,jj) <= 0) THEN
            sshn_u(ji,jj) = sshn_t(ji+1,jj)
         ELSE IF(tmask(ji+1,jj) <= 0) THEN
            sshn_u(ji,jj) = sshn_t(ji,jj)
         END IF

      end do
    end do
!$acc end parallel

!$acc parallel present(tmask, area_t, area_v, sshn_t, sshn_v)
!$acc loop collapse(2)
    do jj = 2, N-1, 1
      do ji = 2, M, 1

!        call next_sshv_code(ji, jj,                   &
!                            sshn_v, sshn_t, &
!                            sshn_t%grid%tmask,        &
!                            sshn_t%grid%area_t, sshn_t%grid%area_v)
 
         if(tmask(ji,jj) + &
            tmask(ji,jj+1) <= 0)  cycle !jump over non-computational domain
         if(tmask(ji,jj) * tmask(ji,jj+1) > 0) then
            rtmp1 = area_t(ji,jj)*sshn_t(ji,jj) + &
                 area_t(ji,jj+1) * sshn_t(ji,jj+1)
            sshn_v(ji,jj) = 0.5_go_wp * rtmp1 / area_v(ji,jj) 
         else if(tmask(ji,jj) <= 0) then
            sshn_v(ji,jj) = sshn_t(ji,jj+1)
         else if(tmask(ji,jj+1) <= 0) then
            sshn_v(ji,jj) = sshn_t(ji,jj)
         end if
      end do
    end do
!$acc end parallel

!    call timer_stop(idxt)

  end subroutine invoke_time_step_arrays

end module time_step_mod
