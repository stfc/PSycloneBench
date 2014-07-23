module momentum_mod
  use kind_params_mod
  use physical_params_mod
  implicit none

  private

  public momentum

contains

  subroutine momentum(grid, jpi, jpj, ua, va, un, vn, &
                      sshn, sshn_u, sshn_v, ssha_u, ssha_v, &
                      pt, hu, hv, ht)
    use grid_mod
    use model_mod, only: rdt, cbfr, visc
    implicit none
    type(grid_type), intent(in) :: grid
    integer, intent(in) :: jpi, jpj
    real(wp), dimension(:,:), intent(out) :: ua, va
    real(wp), dimension(:,:), intent(in)  :: un, vn
    real(wp), dimension(:,:), intent(in)  :: sshn, sshn_u, sshn_v, ssha_u, ssha_v
    real(wp), dimension(:,:), intent(in)  :: hu, hv, ht
    integer,  dimension(:,:), intent(in)  :: pt
    ! Locals
    REAL(wp) :: u_e, u_w
    REAL(wp) :: v_s, v_n
    REAL(wp) :: v_sc, v_nc, u_ec, u_wc
    REAL(wp) :: uu_e, uu_w, uu_s, uu_n
    REAL(wp) :: vv_e, vv_w, vv_s, vv_n
    REAL(wp) :: depe, depw, deps, depn
    REAL(wp) :: dudx_e, dudy_n, dvdx_e, dvdy_n
    REAL(wp) :: dudx_w, dudy_s, dvdx_w, dvdy_s

    REAL(wp) :: adv, vis, hpg, cor

    integer :: ji, jj

    ! u equation
    DO jj = 1, jpj
       DO ji = 1, jpi-1

          ! kernel u adv 

          IF(pt(ji,jj) + pt(ji+1,jj) <= 0)  CYCLE         !jump over non-computatinal domain
          IF(pt(ji,jj) <= 0 .OR. pt(ji+1,jj) <= 0)  CYCLE           !jump over boundary u

          u_e  = 0.5 * (un(ji,jj) + un(ji+1,jj)) * grid%e2t(ji+1,jj)      !add length scale.
          depe = ht(ji+1,jj) + sshn(ji+1,jj)

          u_w  = 0.5 * (un(ji,jj) + un(ji-1,jj)) * grid%e2t(ji,jj)        !add length scale
          depw = ht(ji,jj) + sshn(ji,jj)

          v_sc = 0.5_wp * (vn(ji,jj-1) + vn(ji+1,jj-1))
          v_s  = 0.5_wp * v_sc * (grid%e1v(ji,jj-1) + grid%e1v(ji+1,jj-1))
          deps = 0.5_wp * (hv(ji,jj-1) + sshn_v(ji,jj-1) + hv(ji+1,jj-1) + sshn_v(ji+1,jj-1))

          v_nc = 0.5_wp * (vn(ji,jj) + vn(ji+1,jj))
          v_n  = 0.5_wp * v_nc * (grid%e1v(ji,jj) + grid%e1v(ji+1,jj))
          depn = 0.5_wp * (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + sshn_v(ji+1,jj))

          ! -advection (currently first order upwind)
          uu_w = (0.5_wp - SIGN(0.5_wp, u_w)) * un(ji,jj)              + & 
               & (0.5_wp + SIGN(0.5_wp, u_w)) * un(ji-1,jj) 
          uu_e = (0.5_wp + SIGN(0.5_wp, u_e)) * un(ji,jj)              + & 
               & (0.5_wp - SIGN(0.5_wp, u_e)) * un(ji+1,jj) 

          IF(pt(ji,jj-1) <=0 .OR. pt(ji+1,jj-1) <= 0) THEN   
             uu_s = (0.5_wp - SIGN(0.5_wp, v_s)) * un(ji,jj)   
          ELSE
             uu_s = (0.5_wp - SIGN(0.5_wp, v_s)) * un(ji,jj)              + & 
                  & (0.5_wp + SIGN(0.5_wp, v_s)) * un(ji,jj-1) 
          END If

          IF(pt(ji,jj+1) <=0 .OR. pt(ji+1,jj+1) <= 0) THEN   
             uu_n = (0.5_wp + SIGN(0.5_wp, v_n)) * un(ji,jj)
          ELSE
             uu_n = (0.5_wp + SIGN(0.5_wp, v_n)) * un(ji,jj)              + & 
                  & (0.5_wp - SIGN(0.5_wp, v_n)) * un(ji,jj+1)
          END IF

          adv = uu_w * u_w * depw - uu_e * u_e * depe + uu_s * v_s * deps - uu_n * v_n * depn
          !end kernel u adv 

          ! -viscosity

          !kernel  u vis 
          dudx_e = (un(ji+1,jj) - un(ji,  jj)) / grid%e1t(ji+1,jj) * (ht(ji+1,jj) + sshn(ji+1,jj))
          dudx_w = (un(ji,  jj) - un(ji-1,jj)) / grid%e1t(ji,  jj) * (ht(ji,  jj) + sshn(ji,  jj))
          IF(pt(ji,jj-1) <=0 .OR. pt(ji+1,jj-1) <= 0) THEN   
             dudy_s = 0.0_wp !slip boundary
          ELSE
             dudy_s = (un(ji,jj) - un(ji,jj-1)) / (grid%e2u(ji,jj) + grid%e2u(ji,jj-1)) * &
                  & (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj-1) + sshn_u(ji,jj-1))
          END IF

          IF(pt(ji,jj+1) <= 0 .OR. pt(ji+1,jj+1) <= 0) THEN   
             dudy_n = 0.0_wp ! slip boundary
          ELSE
             dudy_n = (un(ji,jj+1) - un(ji,jj)) / (grid%e2u(ji,jj) + grid%e2u(ji,jj+1)) * &
                  & (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj+1) + sshn_u(ji,jj+1))
          END If

          vis = (dudx_e - dudx_w ) * grid%e2u(ji,jj)  + &
               & (dudy_n - dudy_s ) * grid%e1u(ji,jj) * 0.5_wp  
          vis = visc * vis   !visc will be an array visc(1:jpijglou) 
          !for variable viscosity, such as turbulent viscosity
          !End  kernel u vis 

          ! -Coriolis' force (can be implemented implicitly)
          !kernel cor 
          cor = 0.5_wp * (2._wp * omega * SIN(grid%gphiu(ji,jj) * d2r) * (v_sc + v_nc)) * &
               & grid%e12u(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj))
          !end kernel cor 

          ! -pressure gradient
          !start kernel hpg 
          hpg = -g * (hu(ji,jj) + sshn_u(ji,jj)) * grid%e2u(ji,jj) * (sshn(ji+1,jj) - sshn(ji,jj))
          !end kernel hpg 
          ! -linear bottom friction (implemented implicitly.
          !kernel ua calculation 
          ua(ji,jj) = (un(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj)) + rdt * (adv + vis + cor + hpg) / grid%e12u(ji,jj)) / &
               & (hu(ji,jj) + ssha_u(ji,jj)) / (1.0_wp + cbfr * rdt) 
          !end kernel ua 
       END DO
    END DO

    ! v equation
    DO jj = 1, jpj-1
       DO ji = 1, jpi

          IF(pt(ji,jj) + pt(ji+1,jj) <= 0)  CYCLE                    !jump over non-computatinal domain
          IF(pt(ji,jj) <= 0 .OR. pt(ji,jj+1) <= 0)  CYCLE                !jump over v boundary cells

          ! kernel v adv 
          v_n  = 0.5 * (vn(ji,jj) + vn(ji,jj+1)) * grid%e1t(ji,jj+1)  !add length scale.
          depn = ht(ji,jj+1) + sshn(ji,jj+1)

          v_s  = 0.5 * (vn(ji,jj) + vn(ji,jj-1)) * grid%e1t(ji,jj)    !add length scale
          deps = ht(ji,jj) + sshn(ji,jj)

          u_wc = 0.5_wp * (un(ji-1,jj) + un(ji-1,jj+1))
          u_w  = 0.5_wp * u_wc * (grid%e2u(ji-1,jj) + grid%e2u(ji-1,jj+1))
          depw = 0.50_wp * (hu(ji-1,jj) + sshn_u(ji-1,jj) + hu(ji-1,jj+1) + sshn_u(ji-1,jj+1))

          u_ec = 0.5_wp * (un(ji,jj) + un(ji,jj+1))
          u_e  = 0.5_wp * u_ec * (grid%e2u(ji,jj) + grid%e2u(ji,jj+1))
          depe = 0.50_wp * (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj+1) + sshn_u(ji,jj+1))

          ! -advection (currently first order upwind)
          vv_s = (0.5_wp - SIGN(0.5_wp, v_s)) * vn(ji,jj)              + & 
               & (0.5_wp + SIGN(0.5_wp, v_s)) * vn(ji,jj-1) 
          vv_n = (0.5_wp + SIGN(0.5_wp, v_n)) * vn(ji,jj)              + & 
               & (0.5_wp - SIGN(0.5_wp, v_n)) * vn(ji,jj+1) 

          IF(pt(ji-1,jj) <= 0 .OR. pt(ji-1,jj+1) <= 0) THEN   
             vv_w = (0.5_wp - SIGN(0.5_wp, u_w)) * vn(ji,jj)  
          ELSE
             vv_w = (0.5_wp - SIGN(0.5_wp, u_w)) * vn(ji,jj)              + & 
                  & (0.5_wp + SIGN(0.5_wp, u_w)) * vn(ji-1,jj) 
          END If

          IF(pt(ji+1,jj) <= 0 .OR. pt(ji+1,jj+1) <= 0) THEN
             vv_e = (0.5_wp + SIGN(0.5_wp, u_e)) * vn(ji,jj)
          ELSE
             vv_e = (0.5_wp + SIGN(0.5_wp, u_e)) * vn(ji,jj)              + & 
                  & (0.5_wp - SIGN(0.5_wp, u_e)) * vn(ji+1,jj)
          END IF

          adv = vv_w * u_w * depw - vv_e * u_e * depe + vv_s * v_s * deps - vv_n * v_n * depn

          !end kernel v adv 

          ! -viscosity


          !kernel v dis 
          dvdy_n = (vn(ji,jj+1) - vn(ji,  jj)) / grid%e2t(ji,jj+1) * (ht(ji,jj+1) + sshn(ji,jj+1))
          dvdy_s = (vn(ji,  jj) - vn(ji,jj-1)) / grid%e2t(ji,  jj) * (ht(ji,  jj) + sshn(ji,  jj))

          IF(pt(ji-1,jj) <= 0 .OR. pt(ji-1,jj+1) <= 0) THEN
             dvdx_w = 0.0_wp !slip boundary
          ELSE
             dvdx_w = (vn(ji,jj) - vn(ji-1,jj)) / (grid%e1v(ji,jj) + grid%e1v(ji-1,jj)) * &
                  & (hv(ji,jj) + sshn_v(ji,jj) + hv(ji-1,jj) + sshn_v(ji-1,jj))
          END IF

          IF(pt(ji+1,jj) <= 0 .OR. pt(ji+1,jj+1) <= 0) THEN
             dvdx_e = 0.0_wp ! slip boundary
          ELSE
             dvdx_e = (vn(ji+1,jj) - vn(ji,jj)) / (grid%e1v(ji,jj) + grid%e1v(ji+1,jj)) * &
                  & (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + sshn_v(ji+1,jj))
          END If

          vis = (dvdy_n - dvdy_s ) * grid%e1v(ji,jj)  + &
               & (dvdx_e - dvdx_w ) * grid%e2v(ji,jj) * 0.5_wp  

          vis = visc * vis   !visc will be a array visc(1:jpijglou) 
          !for variable viscosity, such as turbulent viscosity
          !end kernel v dis 

          ! -Coriolis' force (can be implemented implicitly)
          !kernel v cor 
          cor = -0.5_wp * (2._wp * omega * SIN(grid%gphiv(ji,jj) * d2r) * (u_ec + u_wc)) * &
               & grid%e12v(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj))
          !end kernel v cor 

          ! -pressure gradient
          !kernel v hpg 
          hpg = -g * (hv(ji,jj) + sshn_v(ji,jj)) * grid%e1v(ji,jj) * (sshn(ji,jj+1) - sshn(ji,jj))
          !kernel v hpg 

          ! -linear bottom friction (implemented implicitly.
          !kernel ua calculation 
          va(ji,jj) = (vn(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj)) + rdt * (adv + vis + cor + hpg) / grid%e12v(ji,jj) ) / &
               & ((hv(ji,jj) + ssha_v(ji,jj))) / (1.0_wp + cbfr * rdt) 
          !end kernel ua calculation 
       END DO
    END DO
  END SUBROUTINE momentum
 
end module momentum_mod
