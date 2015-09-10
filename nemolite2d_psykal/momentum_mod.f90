module momentum_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use grid_mod
  use field_mod
  implicit none

  private

  public invoke_momentum_u, invoke_momentum_v
  public momentum_u, momentum_v
  public momentum_u_code, momentum_v_code

  !=======================================

  type, extends(kernel_type) :: momentum_u
     type(arg), dimension(18) :: meta_args =  &
          (/ arg(READWRITE, CU, POINTWISE),  & ! ua
             arg(READ,      CU, POINTWISE),  & ! un
             arg(READ,      CV, POINTWISE),  & ! vn
             arg(READ,      CU, POINTWISE),  & ! hu
             arg(READ,      CV, POINTWISE),  & ! hv
             arg(READ,      CT, POINTWISE),  & ! ht
             arg(READ,      CU, POINTWISE),  & ! ssha_u
             arg(READ,      CT, POINTWISE),  & ! sshn_t
             arg(READ,      CU, POINTWISE),  & ! sshn_u
             arg(READ,      CV, POINTWISE),  & ! sshn_v
             arg(READ,      GRID_MASK_T),    &
             arg(READ,      GRID_DX_U),      &
             arg(READ,      GRID_DX_V),      &
             arg(READ,      GRID_DX_T),      &
             arg(READ,      GRID_DY_U),      &
             arg(READ,      GRID_DY_T),      &
             arg(READ,      GRID_AREA_U),    &
             arg(READ,      GRID_LAT_U)      &
           /)

     !> We update only those points within the internal region
     !! of the simulated domain.
     integer :: ITERATES_OVER = INTERNAL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => momentum_u_code
  end type momentum_u

  !=======================================

  type, extends(kernel_type) :: momentum_v
     type(arg), dimension(18) :: meta_args =  &
          (/ arg(READWRITE, CV, POINTWISE),  & ! va
             arg(READ,      CU, POINTWISE),  & ! un
             arg(READ,      CV, POINTWISE),  & ! vn
             arg(READ,      CU, POINTWISE),  & ! hu
             arg(READ,      CV, POINTWISE),  & ! hv
             arg(READ,      CT, POINTWISE),  & ! ht
             arg(READ,      CV, POINTWISE),  & ! ssha_v
             arg(READ,      CT, POINTWISE),  & ! sshn_t
             arg(READ,      CU, POINTWISE),  & ! sshn_u
             arg(READ,      CV, POINTWISE),  & ! sshn_v
             arg(READ,      GRID_MASK_T),    &
             arg(READ,      GRID_DX_V),      &
             arg(READ,      GRID_DX_T),      &
             arg(READ,      GRID_DY_U),      &
             arg(READ,      GRID_DY_V),      &
             arg(READ,      GRID_DY_T),      &
             arg(READ,      GRID_AREA_V),    &
             arg(READ,      GRID_LAT_V)      &
           /)

     !> We update only those points within the internal region
     !! of the simulated domain.
     integer :: ITERATES_OVER = INTERNAL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => momentum_v_code
  end type momentum_v

  !=======================================

contains

  !====================================================
 
  subroutine invoke_momentum_u(ua_fld, un_fld, vn_fld, &
                               hu_fld, hv_fld, ht_fld, ssha_u_fld, &
                               sshn_t_fld, sshn_u_fld, sshn_v_fld)
    implicit none
    type(r2d_field), intent(inout) :: ua_fld
    type(r2d_field), intent(in) :: un_fld, vn_fld
    type(r2d_field), intent(in) :: hu_fld, hv_fld, ht_fld
    type(r2d_field), intent(in) :: ssha_u_fld, sshn_t_fld, &
                                        sshn_u_fld, sshn_v_fld
    ! Locals
    integer :: ji, jj

    do jj = ua_fld%internal%ystart, ua_fld%internal%ystop, 1
      do ji = ua_fld%internal%xstart, ua_fld%internal%xstop, 1

        call momentum_u_code(ji, jj, &
                             ua_fld%data, un_fld%data, vn_fld%data, &
                             hu_fld%data, hv_fld%data, ht_fld%data, &
                             ssha_u_fld%data, sshn_t_fld%data,      &
                             sshn_u_fld%data, sshn_v_fld%data, &
                             ua_fld%grid%tmask,  &
                             ua_fld%grid%dx_u,   &
                             ua_fld%grid%dx_v,   &
                             ua_fld%grid%dx_t,   &
                             ua_fld%grid%dy_u,   &
                             ua_fld%grid%dy_t,   &
                             ua_fld%grid%area_u, &
                             ua_fld%grid%gphiu)
      end do
   end do

  end subroutine invoke_momentum_u

  !====================================================

  subroutine momentum_u_code(ji, jj, &
                             ua, un, vn, &
                             hu, hv, ht, ssha_u, &
                             sshn, sshn_u, sshn_v, &
                             tmask, e1u, e1v, e1t, e2u, e2t, e12u, gphiu)
    use physical_params_mod
    use model_mod, only: rdt, cbfr, visc
    implicit none
    integer, intent(in) :: ji, jj
    integer,  dimension(:,:), intent(in) :: tmask
    real(wp), dimension(:,:), intent(in) :: e1u, e1v, e1t, e12u, e2u, e2t, gphiu
    real(wp), dimension(:,:), intent(in) :: hu, hv, ht
    real(wp), dimension(:,:), intent(in) :: ssha_u, sshn, sshn_u, sshn_v
    real(wp), dimension(:,:), intent(in) :: un, vn
    real(wp), dimension(:,:), intent(out) :: ua
    ! Locals
    REAL(wp) :: u_e, u_w, v_n, v_s
    real(wp) :: v_nc, v_sc
    real(wp) :: depe, depw, deps, depn
    real(wp) :: hpg, adv, cor, vis
    real(wp) :: dudx_e, dudx_w
    real(wp) :: dudy_s, dudy_n
    real(wp) :: uu_e, uu_n, uu_s, uu_w

    IF(tmask(ji,jj) + tmask(ji+1,jj) <= 0)  RETURN   !jump over non-computational domain
    IF(tmask(ji,jj) <= 0 .OR. tmask(ji+1,jj) <= 0)  RETURN !jump over boundary u

    u_e  = 0.5 * (un(ji,jj) + un(ji+1,jj)) * e2t(ji+1,jj)   !add length scale.
    depe = ht(ji+1,jj) + sshn(ji+1,jj)

    u_w  = 0.5 * (un(ji,jj) + un(ji-1,jj)) * e2t(ji,jj)     !add length scale
    depw = ht(ji,jj) + sshn(ji,jj)

    v_sc = 0.5_wp * (vn(ji,jj-1) + vn(ji+1,jj-1))
    v_s  = 0.5_wp * v_sc * (e1v(ji,jj-1) + e1v(ji+1,jj-1))
    deps = 0.5_wp * (hv(ji,jj-1) + sshn_v(ji,jj-1) + hv(ji+1,jj-1) + &
                     sshn_v(ji+1,jj-1))

    v_nc = 0.5_wp * (vn(ji,jj) + vn(ji+1,jj))
    v_n  = 0.5_wp * v_nc * (e1v(ji,jj) + e1v(ji+1,jj))
    depn = 0.5_wp * (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + &
                     sshn_v(ji+1,jj))

    ! -advection (currently first order upwind)
    uu_w = (0.5_wp - SIGN(0.5_wp, u_w)) * un(ji,jj)              + & 
         & (0.5_wp + SIGN(0.5_wp, u_w)) * un(ji-1,jj) 
    uu_e = (0.5_wp + SIGN(0.5_wp, u_e)) * un(ji,jj)              + & 
         & (0.5_wp - SIGN(0.5_wp, u_e)) * un(ji+1,jj) 

    IF(tmask(ji,jj-1) <=0 .OR. tmask(ji+1,jj-1) <= 0) THEN   
       uu_s = (0.5_wp - SIGN(0.5_wp, v_s)) * un(ji,jj)   
    ELSE
       uu_s = (0.5_wp - SIGN(0.5_wp, v_s)) * un(ji,jj)              + & 
            & (0.5_wp + SIGN(0.5_wp, v_s)) * un(ji,jj-1) 
    END If

    IF(tmask(ji,jj+1) <=0 .OR. tmask(ji+1,jj+1) <= 0) THEN   
       uu_n = (0.5_wp + SIGN(0.5_wp, v_n)) * un(ji,jj)
    ELSE
       uu_n = (0.5_wp + SIGN(0.5_wp, v_n)) * un(ji,jj)              + & 
            & (0.5_wp - SIGN(0.5_wp, v_n)) * un(ji,jj+1)
    END IF

    adv = uu_w * u_w * depw - uu_e * u_e * depe + &
          uu_s * v_s * deps - uu_n * v_n * depn
    !end kernel u adv 

    ! -viscosity

    !kernel  u vis 
    dudx_e = (un(ji+1,jj) - un(ji,  jj)) / e1t(ji+1,jj) * &
             (ht(ji+1,jj) + sshn(ji+1,jj))
    dudx_w = (un(ji,  jj) - un(ji-1,jj)) / e1t(ji,  jj) * &
             (ht(ji,  jj) + sshn(ji,  jj))
    IF(tmask(ji,jj-1) <=0 .OR. tmask(ji+1,jj-1) <= 0) THEN   
       dudy_s = 0.0_wp !slip boundary
    ELSE
       dudy_s = (un(ji,jj) - un(ji,jj-1)) / (e2u(ji,jj) + e2u(ji,jj-1)) * &
            & (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj-1) + sshn_u(ji,jj-1))
    END IF

    IF(tmask(ji,jj+1) <= 0 .OR. tmask(ji+1,jj+1) <= 0) THEN   
       dudy_n = 0.0_wp ! slip boundary
    ELSE
       dudy_n = (un(ji,jj+1) - un(ji,jj)) / (e2u(ji,jj) + e2u(ji,jj+1)) * &
            & (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj+1) + sshn_u(ji,jj+1))
    END If

    vis = (dudx_e - dudx_w ) * e2u(ji,jj)  + &
         & (dudy_n - dudy_s ) * e1u(ji,jj) * 0.5_wp  
    vis = visc * vis   !visc will be an array visc(1:jpijglou) 
    !for variable viscosity, such as turbulent viscosity
    !End  kernel u vis 

    ! -Coriolis' force (can be implemented implicitly)
    !kernel cor 
    cor = 0.5_wp * (2._wp * omega * SIN(gphiu(ji,jj) * d2r) * (v_sc + v_nc)) * &
         & e12u(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj))
    !end kernel cor 

    ! -pressure gradient
    !start kernel hpg 
    hpg = -g * (hu(ji,jj) + sshn_u(ji,jj)) * e2u(ji,jj) * &
           (sshn(ji+1,jj) - sshn(ji,jj))
    !end kernel hpg 
    ! -linear bottom friction (implemented implicitly.
    !kernel ua calculation 
    ua(ji,jj) = (un(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj)) + rdt * &
                 (adv + vis + cor + hpg) / e12u(ji,jj)) / &
                (hu(ji,jj) + ssha_u(ji,jj)) / (1.0_wp + cbfr * rdt) 

  end subroutine momentum_u_code
 
  !====================================================
 
  subroutine invoke_momentum_v(va_fld, un_fld, vn_fld, &
                               hu_fld, hv_fld, ht_fld, ssha_v_fld, &
                               sshn_t_fld, sshn_u_fld, sshn_v_fld)
    implicit none
    type(r2d_field), intent(inout) :: va_fld
    type(r2d_field), intent(in) :: un_fld, vn_fld
    type(r2d_field), intent(in) :: hu_fld, hv_fld, ht_fld
    type(r2d_field), intent(in) :: ssha_v_fld, sshn_t_fld, &
                                   sshn_u_fld, sshn_v_fld
    ! Locals
    integer :: ji, jj

    do jj = va_fld%internal%ystart, va_fld%internal%ystop, 1
      do ji = va_fld%internal%xstart, va_fld%internal%xstop, 1

        call momentum_v_code(ji, jj, &
                             va_fld%data, un_fld%data, vn_fld%data, &
                             hu_fld%data, hv_fld%data, ht_fld%data, &
                             ssha_v_fld%data, sshn_t_fld%data,      &
                             sshn_u_fld%data, sshn_v_fld%data,      &
                             va_fld%grid%tmask, va_fld%grid%dx_v,   &
                             va_fld%grid%dx_t, &
                             va_fld%grid%dy_u, va_fld%grid%dy_v,    &
                             va_fld%grid%dy_t,     &
                             va_fld%grid%area_v, va_fld%grid%gphiv)

      end do
   end do

  end subroutine invoke_momentum_v

  !====================================================

  subroutine momentum_v_code(ji, jj, &
                             va, un, vn, &
                             hu, hv, ht, ssha_v, &
                             sshn, sshn_u, sshn_v, &
                             tmask, e1v, e1t, e2u, e2v, e2t, e12v, gphiv)

    use physical_params_mod
    use model_mod, only: rdt, cbfr, visc
    implicit none
    integer, intent(in) :: ji, jj
    integer,  dimension(:,:), intent(in) :: tmask
    real(wp), dimension(:,:), intent(in) :: e1v, e1t, e12v, e2u, e2v, e2t, gphiv
    real(wp), dimension(:,:), intent(in) :: hu, hv, ht
    real(wp), dimension(:,:), intent(in) :: ssha_v, sshn, sshn_u, sshn_v
    real(wp), dimension(:,:), intent(in) :: un, vn
    real(wp), dimension(:,:), intent(out) :: va
    ! Locals
    REAL(wp) :: u_e, u_w, v_n, v_s
    real(wp) :: u_ec, u_wc, vv_e, vv_n, vv_s, vv_w
    real(wp) :: depe, depw, deps, depn
    real(wp) :: hpg, adv, cor, vis
    real(wp) :: dvdx_e, dvdx_w, dvdy_n, dvdy_s

    IF(tmask(ji,jj) + tmask(ji+1,jj) <= 0)  return !jump over non-computatinal domain
    IF(tmask(ji,jj) <= 0 .OR. tmask(ji,jj+1) <= 0) return !jump over v boundary cells

    ! kernel v adv 
    v_n  = 0.5 * (vn(ji,jj) + vn(ji,jj+1)) * e1t(ji,jj+1)  !add length scale.
    depn = ht(ji,jj+1) + sshn(ji,jj+1)

    v_s  = 0.5 * (vn(ji,jj) + vn(ji,jj-1)) * e1t(ji,jj)    !add length scale
    deps = ht(ji,jj) + sshn(ji,jj)

    u_wc = 0.5_wp * (un(ji-1,jj) + un(ji-1,jj+1))
    u_w  = 0.5_wp * u_wc * (e2u(ji-1,jj) + e2u(ji-1,jj+1))
    depw = 0.50_wp * (hu(ji-1,jj) + sshn_u(ji-1,jj) + &
                      hu(ji-1,jj+1) + sshn_u(ji-1,jj+1))

    u_ec = 0.5_wp * (un(ji,jj) + un(ji,jj+1))
    u_e  = 0.5_wp * u_ec * (e2u(ji,jj) + e2u(ji,jj+1))
    depe = 0.50_wp * (hu(ji,jj) + sshn_u(ji,jj) + &
                      hu(ji,jj+1) + sshn_u(ji,jj+1))

    ! -advection (currently first order upwind)
    vv_s = (0.5_wp - SIGN(0.5_wp, v_s)) * vn(ji,jj)              + & 
         & (0.5_wp + SIGN(0.5_wp, v_s)) * vn(ji,jj-1) 
    vv_n = (0.5_wp + SIGN(0.5_wp, v_n)) * vn(ji,jj)              + & 
         & (0.5_wp - SIGN(0.5_wp, v_n)) * vn(ji,jj+1) 

    IF(tmask(ji-1,jj) <= 0 .OR. tmask(ji-1,jj+1) <= 0) THEN   
       vv_w = (0.5_wp - SIGN(0.5_wp, u_w)) * vn(ji,jj)  
    ELSE
       vv_w = (0.5_wp - SIGN(0.5_wp, u_w)) * vn(ji,jj)              + & 
            & (0.5_wp + SIGN(0.5_wp, u_w)) * vn(ji-1,jj) 
    END If

    IF(tmask(ji+1,jj) <= 0 .OR. tmask(ji+1,jj+1) <= 0) THEN
       vv_e = (0.5_wp + SIGN(0.5_wp, u_e)) * vn(ji,jj)
    ELSE
       vv_e = (0.5_wp + SIGN(0.5_wp, u_e)) * vn(ji,jj)              + & 
              (0.5_wp - SIGN(0.5_wp, u_e)) * vn(ji+1,jj)
    END IF

    adv = vv_w * u_w * depw - vv_e * u_e * depe + &
          vv_s * v_s * deps - vv_n * v_n * depn

    !end kernel v adv 

    ! -viscosity

    
    !kernel v dis 
    dvdy_n = (vn(ji,jj+1) - vn(ji,  jj)) / e2t(ji,jj+1) * &
                          (ht(ji,jj+1) + sshn(ji,jj+1))
    dvdy_s = (vn(ji,  jj) - vn(ji,jj-1)) / e2t(ji,  jj) * &
                          (ht(ji,  jj) + sshn(ji,  jj))

    IF(tmask(ji-1,jj) <= 0 .OR. tmask(ji-1,jj+1) <= 0) THEN
       dvdx_w = 0.0_wp !slip boundary
    ELSE
       dvdx_w = (vn(ji,jj) - vn(ji-1,jj)) / (e1v(ji,jj) + e1v(ji-1,jj)) * &
                (hv(ji,jj) + sshn_v(ji,jj) + hv(ji-1,jj) + sshn_v(ji-1,jj))
    END IF

    IF(tmask(ji+1,jj) <= 0 .OR. tmask(ji+1,jj+1) <= 0) THEN
       dvdx_e = 0.0_wp ! slip boundary
    ELSE
       dvdx_e = (vn(ji+1,jj) - vn(ji,jj)) / (e1v(ji,jj) + e1v(ji+1,jj)) * &
                  (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + sshn_v(ji+1,jj))
    END If

    vis = (dvdy_n - dvdy_s ) * e1v(ji,jj)  + &
          (dvdx_e - dvdx_w ) * e2v(ji,jj) * 0.5_wp  

    vis = visc * vis   !visc will be a array visc(1:jpijglou) 
    !for variable viscosity, such as turbulent viscosity
    !end kernel v dis 

    ! -Coriolis' force (can be implemented implicitly)
    !kernel v cor 
    cor = -0.5_wp*(2._wp * omega * SIN(gphiv(ji,jj) * d2r) * (u_ec + u_wc)) * &
               e12v(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj))
    !end kernel v cor 

    ! -pressure gradient
    !kernel v hpg 
    hpg = -g * (hv(ji,jj) + sshn_v(ji,jj)) * e1v(ji,jj) * &
           (sshn(ji,jj+1) - sshn(ji,jj))
    !kernel v hpg 

    ! -linear bottom friction (implemented implicitly.
    !kernel ua calculation 
    va(ji,jj) = (vn(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj)) + &
                 rdt * (adv + vis + cor + hpg) / e12v(ji,jj) ) / &
                 ((hv(ji,jj) + ssha_v(ji,jj))) / (1.0_wp + cbfr * rdt) 

  end subroutine momentum_v_code

  !====================================================

end module momentum_mod
