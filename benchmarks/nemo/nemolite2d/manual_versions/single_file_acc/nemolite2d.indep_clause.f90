    module sizes_mod

         INTEGER,  PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)
         INTEGER,  PARAMETER :: dp = SELECTED_REAL_KIND(12, 307)
         INTEGER,  PARAMETER :: wp = dp         
                    
   end module sizes_mod          
                    
    module cuda_mod
         use cudafor
         use sizes_mod
         IMPLICIT NONE

     contains    

              
        ! Assumes padded allocations
        SUBROUTINE momentum_acc(un, vn, sshn, sshn_u, sshn_v, ssha_u, ssha_v, &
                                                   hu, hv, ht, e1v, e2v, e1u, e2u, e1t, e2t, &
                                                   gphiu, gphiv, e12u, e12v, pt, ua, va ,&
                                                   ssha, e12t, &
                                                   visc, omega, d2r, g, rdt, cbfr, jpi, jpj)
                                               
          REAL(wp), intent(in), device, dimension(-31:jpi,1:jpj) :: un, hu, gphiu, e1u, e2u, sshn_u, ssha_u, e12u
          REAL(wp), intent(in), device, dimension(1:jpi,0:jpj)   :: vn, hv, gphiv, e1v, e2v, sshn_v, ssha_v, e12v
          REAL(wp), intent(in), device, dimension(1:jpi,1:jpj)   :: sshn, ht, e1t, e2t, e12t
          REAL(wp), intent(out), device, dimension(1:jpi,1:jpj)  :: ssha
          
          
          INTEGER, intent(in), device, dimension(0:jpi+1+30,0:jpj+1) :: pt
          INTEGER,  value, intent(in) :: jpi, jpj
          REAL(wp), value, intent(in) :: visc, omega, d2r, g, rdt, cbfr
          
          REAL(wp), intent(out), device, dimension(-31:jpi,1:jpj) :: ua
          REAL(wp), intent(out), device, dimension(1:jpi,0:jpj) :: va
                                                   
          REAL(wp) :: rtmp1, rtmp2, rtmp3, rtmp4
                                       
          REAL(wp) :: u_e, u_w
          REAL(wp) :: v_s, v_n
          REAL(wp) :: v_sc, v_nc, u_ec, u_wc
          REAL(wp) :: uu_e, uu_w, uu_s, uu_n
          REAL(wp) :: vv_e, vv_w, vv_s, vv_n
          REAL(wp) :: depe, depw, deps, depn
          REAL(wp) :: dudx_e, dudy_n, dvdx_e, dvdy_n
          REAL(wp) :: dudx_w, dudy_s, dvdx_w, dvdy_s

          REAL(wp) :: adv, vis, hpg, cor
          INTEGER  :: ji, jj, js, jii, je                          !temporary loop index

          ! Reads: un, vn, sshn, sshn_u, sshn_v, ssha_u, ssha_v, hu, hv, ht, e1v, e2v, e1u, e2u, e1t, e2t, gphiu, gphiv, e12u, e12v (20)
          ! Reads (int): pt (1)
          ! Writes: ua, va (2)
          ! = 180 bytes
          ! For 1026x1026 = 181 MB.  

          ! u equation
          !$acc kernels present(un, vn, hu, ht, sshn_u, sshn, e1t, e2t, e2u)
          ! acc loop shortloop
          DO jj = 1, jpj 
          ! acc loop shortloop
          DO ji = 1, jpi
! kernel u adv 
             IF (ji <= jpi-1) THEN ! .AND. &
               ! (pt(ji,jj) + pt(ji+1,jj) > 0) .AND. &                 !jump over non-computatinal domain
               ! (pt(ji,jj) > 0 .AND. pt(ji+1,jj) > 0))  THEN                    !jump over boundary u
               
             ! IF(pt(ji,jj) + pt(ji+1,jj) <= 0)  CYCLE                    !jump over non-computatinal domain
             ! IF(pt(ji,jj) <= 0 .OR. pt(ji+1,jj) <= 0)  CYCLE                    !jump over boundary u

             u_e  = 0.5 * (un(ji,jj) + un(ji+1,jj)) * e2t(ji+1,jj)      !add length scale.
             depe = ht(ji+1,jj) + sshn(ji+1,jj)

             u_w  = 0.5 * (un(ji,jj) + un(ji-1,jj)) * e2t(ji,jj)        !add length scale
             depw = ht(ji,jj) + sshn(ji,jj)

             v_sc = 0.5_wp * (vn(ji,jj-1) + vn(ji+1,jj-1))
             v_s  = 0.5_wp * v_sc * (e1v(ji,jj-1) + e1v(ji+1,jj-1))
             deps = 0.5_wp * (hv(ji,jj-1) + sshn_v(ji,jj-1) + hv(ji+1,jj-1) + sshn_v(ji+1,jj-1))

             v_nc = 0.5_wp * (vn(ji,jj) + vn(ji+1,jj))
             v_n  = 0.5_wp * v_nc * (e1v(ji,jj) + e1v(ji+1,jj))
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
            dudx_e = (un(ji+1,jj) - un(ji,  jj)) / e1t(ji+1,jj) * (ht(ji+1,jj) + sshn(ji+1,jj))
            dudx_w = (un(ji,  jj) - un(ji-1,jj)) / e1t(ji,  jj) * (ht(ji,  jj) + sshn(ji,  jj))
            IF(pt(ji,jj-1) <=0 .OR. pt(ji+1,jj-1) <= 0) THEN   
              dudy_s = 0.0_wp !slip boundary
            ELSE
              dudy_s = (un(ji,jj) - un(ji,jj-1)) / (e2u(ji,jj) + e2u(ji,jj-1)) * &
                     & (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj-1) + sshn_u(ji,jj-1))
            END IF

            IF(pt(ji,jj+1) <= 0 .OR. pt(ji+1,jj+1) <= 0) THEN   
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
            hpg = -g * (hu(ji,jj) + sshn_u(ji,jj)) * e2u(ji,jj) * (sshn(ji+1,jj) - sshn(ji,jj))
!end kernel hpg 
            ! -linear bottom friction (implemented implicitly.
!kernel ua calculation 
            ua(ji,jj) = (un(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj)) + rdt * (adv + vis + cor + hpg) / e12u(ji,jj)) / &
                   & (hu(ji,jj) + ssha_u(ji,jj)) / (1.0_wp + cbfr * rdt) 
!end kernel ua 
            END IF
         END DO
       END DO

       ! v equation
       ! acc loop shortloop
       DO jj = 1, jpj 
          ! acc loop shortloop
          DO ji = 1, jpi
        
             IF (jj <= jpj-1) THEN ! .AND. &
               ! pt(ji,jj) + pt(ji+1,jj) >= 0 .AND. &                         !jump over non-computatinal domain
               ! (pt(ji,jj) > 0 .AND. pt(ji,jj+1) > 0))  THEN                !jump over v boundary cells
               
             ! IF(pt(ji,jj) + pt(ji+1,jj) <= 0)  CYCLE                    !jump over non-computatinal domain
             ! IF(pt(ji,jj) <= 0 .OR. pt(ji,jj+1) <= 0)  CYCLE                !jump over v boundary cells

! kernel v adv 
             v_n  = 0.5 * (vn(ji,jj) + vn(ji,jj+1)) * e1t(ji,jj+1)  !add length scale.
             depn = ht(ji,jj+1) + sshn(ji,jj+1)

             v_s  = 0.5 * (vn(ji,jj) + vn(ji,jj-1)) * e1t(ji,jj)    !add length scale
             deps = ht(ji,jj) + sshn(ji,jj)

             u_wc = 0.5_wp * (un(ji-1,jj) + un(ji-1,jj+1))
             u_w  = 0.5_wp * u_wc * (e2u(ji-1,jj) + e2u(ji-1,jj+1))
             depw = 0.50_wp * (hu(ji-1,jj) + sshn_u(ji-1,jj) + hu(ji-1,jj+1) + sshn_u(ji-1,jj+1))

             u_ec = 0.5_wp * (un(ji,jj) + un(ji,jj+1))
             u_e  = 0.5_wp * u_ec * (e2u(ji,jj) + e2u(ji,jj+1))
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
            dvdy_n = (vn(ji,jj+1) - vn(ji,  jj)) / e2t(ji,jj+1) * (ht(ji,jj+1) + sshn(ji,jj+1))
            dvdy_s = (vn(ji,  jj) - vn(ji,jj-1)) / e2t(ji,  jj) * (ht(ji,  jj) + sshn(ji,  jj))

            IF(pt(ji-1,jj) <= 0 .OR. pt(ji-1,jj+1) <= 0) THEN
              dvdx_w = 0.0_wp !slip boundary
            ELSE
              dvdx_w = (vn(ji,jj) - vn(ji-1,jj)) / (e1v(ji,jj) + e1v(ji-1,jj)) * &
                     & (hv(ji,jj) + sshn_v(ji,jj) + hv(ji-1,jj) + sshn_v(ji-1,jj))
            END IF

            IF(pt(ji+1,jj) <= 0 .OR. pt(ji+1,jj+1) <= 0) THEN
              dvdx_e = 0.0_wp ! slip boundary
            ELSE
              dvdx_e = (vn(ji+1,jj) - vn(ji,jj)) / (e1v(ji,jj) + e1v(ji+1,jj)) * &
                     & (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + sshn_v(ji+1,jj))
            END If

            vis = (dvdy_n - dvdy_s ) * e1v(ji,jj)  + &
                & (dvdx_e - dvdx_w ) * e2v(ji,jj) * 0.5_wp  

            vis = visc * vis   !visc will be a array visc(1:jpijglou) 
                               !for variable viscosity, such as turbulent viscosity
!end kernel v dis 

            ! -Coriolis' force (can be implemented implicitly)
!kernel v cor 
            cor = -0.5_wp * (2._wp * omega * SIN(gphiv(ji,jj) * d2r) * (u_ec + u_wc)) * &
                & e12v(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj))
!end kernel v cor 

            ! -pressure gradient
!kernel v hpg 
            hpg = -g * (hv(ji,jj) + sshn_v(ji,jj)) * e1v(ji,jj) * (sshn(ji,jj+1) - sshn(ji,jj))
!kernel v hpg 

            ! -linear bottom friction (implemented implicitly.
!kernel ua calculation 
            va(ji,jj) = (vn(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj)) + rdt * (adv + vis + cor + hpg) / e12v(ji,jj) ) / &
                   & ((hv(ji,jj) + ssha_v(ji,jj))) / (1.0_wp + cbfr * rdt) 
            END IF
!end kernel ua calculation 

              rtmp1 = (sshn_u(ji  ,jj  ) + hu(ji  ,jj  )) * un(ji  ,jj  )
              rtmp2 = (sshn_u(ji-1,jj  ) + hu(ji-1,jj  )) * un(ji-1,jj  )
              rtmp3 = (sshn_v(ji  ,jj  ) + hv(ji  ,jj  )) * vn(ji  ,jj  )
              rtmp4 = (sshn_v(ji  ,jj-1) + hv(ji  ,jj-1)) * vn(ji  ,jj-1)
              ssha(ji,jj) = sshn(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * rdt / e12t(ji,jj)

          ! END IF
          END DO
          END DO
          ! END DO
          !$acc end kernels
          !$acc wait
        END SUBROUTINE momentum_acc
          
        attributes(global) SUBROUTINE momentum_ker(un, vn, sshn, sshn_u, sshn_v, ssha_u, ssha_v, &
                                                   hu, hv, ht, e1v, e2v, e1u, e2u, e1t, e2t, &
                                                   gphiu, gphiv, e12u, e12v, pt, ua, va ,&
                                                   visc, omega, d2r, g, rdt, cbfr, jpi, jpj)
                                               
          REAL(wp), intent(in), device, dimension(-31:jpi,1:jpj) :: un, hu, gphiu, e1u, e2u, sshn_u, ssha_u, e12u
          REAL(wp), intent(in), device, dimension(1:jpi,0:jpj)    :: vn, hv, gphiv, e1v, e2v, sshn_v, ssha_v, e12v
          REAL(wp), intent(in), device, dimension(1:jpi,1:jpj)    :: sshn, ht, e1t, e2t
          
          
          INTEGER, intent(in), device, dimension(0:jpi+1+30,0:jpj+1) :: pt
          INTEGER,  value, intent(in) :: jpi, jpj
          REAL(wp), value, intent(in) :: visc, omega, d2r, g, rdt, cbfr
          
          REAL(wp), intent(out), device, dimension(-31:jpi,1:jpj) :: ua
          REAL(wp), intent(out), device, dimension(1:jpi,0:jpj) :: va
                          
                                                   
          REAL(wp) :: u_e, u_w
          REAL(wp) :: v_s, v_n
          REAL(wp) :: v_sc, v_nc, u_ec, u_wc
          REAL(wp) :: uu_e, uu_w, uu_s, uu_n
          REAL(wp) :: vv_e, vv_w, vv_s, vv_n
          REAL(wp) :: depe, depw, deps, depn
          REAL(wp) :: dudx_e, dudy_n, dvdx_e, dvdy_n
          REAL(wp) :: dudx_w, dudy_s, dvdx_w, dvdy_s

          REAL(wp) :: adv, vis, hpg, cor
          INTEGER  :: ji, jj                          !temporary loop index

          ji = blockDim%x * (blockIdx%x - 1) + threadIdx%x
          jj = blockDim%y * (blockIdx%y - 1) + threadIdx%y
          
          
          ! Reads: un, vn, sshn, sshn_u, sshn_v, ssha_u, ssha_v, hu, hv, ht, e1v, e2v, e1u, e2u, e1t, e2t, gphiu, gphiv, e12u, e12v (20)
          ! Reads (int): pt (1)
          ! Writes: ua, va (2)
          ! = 180 bytes
          ! For 1026x1026 = 181 MB, ie. 1ms if bandwidth bound.

          ! u equation
          ! !$acc kernels present(un, vn, hu, ht, sshn_u, sshn, e1t, e2t, e2u)
          ! !$acc loop collapse(2)
          ! DO jj = 1, jpj
          ! DO ji = 1, jpi
          IF (ji <= jpi .AND. jj <= jpj) THEN
! kernel u adv 
             IF (ji <= jpi-1 .AND. &
               (pt(ji,jj) + pt(ji+1,jj) > 0) .AND. &                 !jump over non-computatinal domain
               (pt(ji,jj) > 0 .AND. pt(ji+1,jj) > 0))  THEN                    !jump over boundary u
               
             ! IF(pt(ji,jj) + pt(ji+1,jj) <= 0)  CYCLE                    !jump over non-computatinal domain
             ! IF(pt(ji,jj) <= 0 .OR. pt(ji+1,jj) <= 0)  CYCLE                    !jump over boundary u

             u_e  = 0.5 * (un(ji,jj) + un(ji+1,jj)) * e2t(ji+1,jj)      !add length scale.
             depe = ht(ji+1,jj) + sshn(ji+1,jj)

             u_w  = 0.5 * (un(ji,jj) + un(ji-1,jj)) * e2t(ji,jj)        !add length scale
             depw = ht(ji,jj) + sshn(ji,jj)

             v_sc = 0.5_wp * (vn(ji,jj-1) + vn(ji+1,jj-1))
             v_s  = 0.5_wp * v_sc * (e1v(ji,jj-1) + e1v(ji+1,jj-1))
             deps = 0.5_wp * (hv(ji,jj-1) + sshn_v(ji,jj-1) + hv(ji+1,jj-1) + sshn_v(ji+1,jj-1))

             v_nc = 0.5_wp * (vn(ji,jj) + vn(ji+1,jj))
             v_n  = 0.5_wp * v_nc * (e1v(ji,jj) + e1v(ji+1,jj))
             depn = 0.5_wp * (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + sshn_v(ji+1,jj))

            ! -advection (currently first order upwind)
            uu_w = (0.5_wp - SIGN(0.5_wp, u_w)) * un(ji,jj)              + & 
                 & (0.5_wp + SIGN(0.5_wp, u_w)) * un(ji-1,jj) 
            uu_e = (0.5_wp + SIGN(0.5_wp, u_e)) * un(ji,jj)              + & 
                 & (0.5_wp - SIGN(0.5_wp, u_e)) * un(ji+1,jj) 

            ! if (u_w >= 0) THEN
               ! uu_w = un(ji,jj)  
            ! else
               ! uu_w = un(ji-1,jj)
            ! END IF
            ! if (u_e >= 0) THEN
               ! uu_e = un(ji,jj)  
            ! else 
               ! uu_e = un(ji+1,jj)
            ! END IF
                 
                 
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
            dudx_e = (un(ji+1,jj) - un(ji,  jj)) / e1t(ji+1,jj) * (ht(ji+1,jj) + sshn(ji+1,jj))
            dudx_w = (un(ji,  jj) - un(ji-1,jj)) / e1t(ji,  jj) * (ht(ji,  jj) + sshn(ji,  jj))
            IF(pt(ji,jj-1) <=0 .OR. pt(ji+1,jj-1) <= 0) THEN   
              dudy_s = 0.0_wp !slip boundary
            ELSE
              dudy_s = (un(ji,jj) - un(ji,jj-1)) / (e2u(ji,jj) + e2u(ji,jj-1)) * &
                     & (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj-1) + sshn_u(ji,jj-1))
            END IF

            IF(pt(ji,jj+1) <= 0 .OR. pt(ji+1,jj+1) <= 0) THEN   
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
            hpg = -g * (hu(ji,jj) + sshn_u(ji,jj)) * e2u(ji,jj) * (sshn(ji+1,jj) - sshn(ji,jj))
!end kernel hpg 
            ! -linear bottom friction (implemented implicitly.
!kernel ua calculation 
            ua(ji,jj) = (un(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj)) + rdt * (adv + vis + cor + hpg) / e12u(ji,jj)) / &
                   & (hu(ji,jj) + ssha_u(ji,jj)) / (1.0_wp + cbfr * rdt) 
!end kernel ua 
            END IF
          ! END DO
          ! END DO

          ! v equation
          ! DO jj = 1, jpj
          ! DO ji = 1, jpi

             IF (jj <= jpj-1 .AND. &
               pt(ji,jj) + pt(ji+1,jj) >= 0 .AND. &                         !jump over non-computatinal domain
               (pt(ji,jj) > 0 .AND. pt(ji,jj+1) > 0))  THEN                !jump over v boundary cells
               
             ! IF(pt(ji,jj) + pt(ji+1,jj) <= 0)  CYCLE                    !jump over non-computatinal domain
             ! IF(pt(ji,jj) <= 0 .OR. pt(ji,jj+1) <= 0)  CYCLE                !jump over v boundary cells

! kernel v adv 
             v_n  = 0.5 * (vn(ji,jj) + vn(ji,jj+1)) * e1t(ji,jj+1)  !add length scale.
             depn = ht(ji,jj+1) + sshn(ji,jj+1)

             v_s  = 0.5 * (vn(ji,jj) + vn(ji,jj-1)) * e1t(ji,jj)    !add length scale
             deps = ht(ji,jj) + sshn(ji,jj)

             u_wc = 0.5_wp * (un(ji-1,jj) + un(ji-1,jj+1))
             u_w  = 0.5_wp * u_wc * (e2u(ji-1,jj) + e2u(ji-1,jj+1))
             depw = 0.50_wp * (hu(ji-1,jj) + sshn_u(ji-1,jj) + hu(ji-1,jj+1) + sshn_u(ji-1,jj+1))

             u_ec = 0.5_wp * (un(ji,jj) + un(ji,jj+1))
             u_e  = 0.5_wp * u_ec * (e2u(ji,jj) + e2u(ji,jj+1))
             depe = 0.50_wp * (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj+1) + sshn_u(ji,jj+1))

            ! -advection (currently first order upwind)
            vv_s = (0.5_wp - SIGN(0.5_wp, v_s)) * vn(ji,jj)              + & 
                 & (0.5_wp + SIGN(0.5_wp, v_s)) * vn(ji,jj-1) 
            vv_n = (0.5_wp + SIGN(0.5_wp, v_n)) * vn(ji,jj)              + & 
                 & (0.5_wp - SIGN(0.5_wp, v_n)) * vn(ji,jj+1) 

            ! if (v_s >= 0) THEN
               ! vv_s = vn(ji,jj)  
            ! else 
               ! vv_s = vn(ji,jj-1)
            ! END IF
            ! if (v_n >= 0) THEN
               ! vv_n = vn(ji,jj)  
            ! else 
               ! vv_n = vn(ji,jj+1)
            ! END IF
                 
                 
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
            dvdy_n = (vn(ji,jj+1) - vn(ji,  jj)) / e2t(ji,jj+1) * (ht(ji,jj+1) + sshn(ji,jj+1))
            dvdy_s = (vn(ji,  jj) - vn(ji,jj-1)) / e2t(ji,  jj) * (ht(ji,  jj) + sshn(ji,  jj))

            IF(pt(ji-1,jj) <= 0 .OR. pt(ji-1,jj+1) <= 0) THEN
              dvdx_w = 0.0_wp !slip boundary
            ELSE
              dvdx_w = (vn(ji,jj) - vn(ji-1,jj)) / (e1v(ji,jj) + e1v(ji-1,jj)) * &
                     & (hv(ji,jj) + sshn_v(ji,jj) + hv(ji-1,jj) + sshn_v(ji-1,jj))
            END IF

            IF(pt(ji+1,jj) <= 0 .OR. pt(ji+1,jj+1) <= 0) THEN
              dvdx_e = 0.0_wp ! slip boundary
            ELSE
              dvdx_e = (vn(ji+1,jj) - vn(ji,jj)) / (e1v(ji,jj) + e1v(ji+1,jj)) * &
                     & (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + sshn_v(ji+1,jj))
            END If

            vis = (dvdy_n - dvdy_s ) * e1v(ji,jj)  + &
                & (dvdx_e - dvdx_w ) * e2v(ji,jj) * 0.5_wp  

            vis = visc * vis   !visc will be a array visc(1:jpijglou) 
                               !for variable viscosity, such as turbulent viscosity
!end kernel v dis 

            ! -Coriolis' force (can be implemented implicitly)
!kernel v cor 
            cor = -0.5_wp * (2._wp * omega * SIN(gphiv(ji,jj) * d2r) * (u_ec + u_wc)) * &
                & e12v(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj))
!end kernel v cor 

            ! -pressure gradient
!kernel v hpg 
            hpg = -g * (hv(ji,jj) + sshn_v(ji,jj)) * e1v(ji,jj) * (sshn(ji,jj+1) - sshn(ji,jj))
!kernel v hpg 

            ! -linear bottom friction (implemented implicitly.
!kernel ua calculation 
            va(ji,jj) = (vn(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj)) + rdt * (adv + vis + cor + hpg) / e12v(ji,jj) ) / &
                   & ((hv(ji,jj) + ssha_v(ji,jj))) / (1.0_wp + cbfr * rdt) 
            END IF
!end kernel ua calculation 
          ! END DO
          ! END DO
          END IF
          ! !$acc end kernels
          ! !$acc wait
        END SUBROUTINE momentum_ker
     
     end module cuda_mod


PROGRAM nemolite2d
         
         !!! A Horizontal 2D hydrodynamic ocean model which
         !!   1) using structured grid
         !!   2) using direct data addressig structures
         use timing_mod
         use field_mod
         use gocean_mod,      only: model_write_log
         use cuda_mod
         use cudafor
         use sizes_mod
         IMPLICIT NONE

         
         REAL(wp), PARAMETER :: pi    = 3.1415926535897932_wp  
         REAL(wp), PARAMETER :: g     = 9.80665_wp                          ! gravity constant
         REAL(wp), PARAMETER :: omega = 7.292116e-05_wp                     ! earth rotation speed (s^(-1))
         REAL(wp), PARAMETER :: d2r   = pi / 180._wp                        ! degree to radian
                                                                            
                                                                            
         INTEGER, ALLOCATABLE :: pt(:,:)                                    ! properties of t-cells 
                                                                            ! 1: water cell within computational domain
                                                                            ! 0: land cell
                                                                            !-1: water cell outside computational domain

         REAL(wp), ALLOCATABLE :: e1t(:,:), e2t(:,:), e1u(:,:), e2u(:,:) 
         REAL(wp), ALLOCATABLE :: e1f(:,:), e2f(:,:), e1v(:,:), e2v(:,:) 
         REAL(wp), ALLOCATABLE :: e12t(:,:), e12u(:,:), e12v(:,:)

         REAL(wp), ALLOCATABLE :: gphiu(:,:), gphiv(:,:), gphif(:,:)

         REAL(wp), ALLOCATABLE :: xt(:,:), yt(:,:)


         REAL(wp), ALLOCATABLE :: ht(:,:), hu(:,:), hv(:,:), hf(:,:) 

         REAL(wp), ALLOCATABLE :: sshb(:,:), sshb_u(:,:), sshb_v(:,:)
         REAL(wp), ALLOCATABLE :: sshn(:,:), sshn_u(:,:), sshn_v(:,:)
         REAL(wp), ALLOCATABLE :: ssha(:,:), ssha_u(:,:), ssha_v(:,:)

         REAL(wp), ALLOCATABLE :: un(:,:),  vn(:,:), ua(:,:),  va(:,:)

         INTEGER  :: jpiglo, jpjglo, jpi, jpj        !dimensions of grid
         INTEGER  :: jphgr_msh                       !type of grid
         INTEGER  :: nit000, nitend, irecord         !start-end and record time steps

         REAL(wp) :: dx, dy, dep_const               !regular grid size and constant depth
         REAL(wp) :: rdt                             !time step
                                                     
         REAL(wp) :: cbfr                            !bottom friction coefficient
         REAL(wp) :: visc                            !backgroud/constant viscosity 

         INTEGER  :: istp                            !time stepping index

         INTEGER  :: itmp1, itmp2                    !integer temporary vars
         integer :: idxt ! Index for main-loop timer

         !! read in model parameters
         CALL setup

         !! allocate memory and read in or setup model grid 
         CALL grid

         !! setup model initial condition
         CALL initialisation

         istp = 0
         CALL output

         
         call acc_init()
         
         !$acc data pcopy(pt) &
         !$acc      pcopy(e1t, e2t, e1u, e2u) &
         !$acc      pcopy(e1f, e2f, e1v, e2v) &
         !$acc      pcopy(e12t, e12u, e12v) &
         !$acc      pcopy(gphiu, gphiv, gphif) &
         !$acc      pcopy(xt, yt) &
         !$acc      pcopy(ht, hu, hv, hf ) &
         !$acc      pcopy(sshb, sshb_u, sshb_v) &
         !$acc      pcopy(sshn, sshn_u, sshn_v) &
         !$acc      pcopy(ssha, ssha_u, ssha_v) &
         !$acc      pcopy(un,  vn, ua,  va)
         
         call timer_start('Time-stepping', idxt, (nitend-nit000+1))

         !! time stepping 
         DO istp = nit000, nitend, 1
           !print*, 'istp == ', istp
           CALL step
         END DO
         
         call timer_stop(idxt)
         !$acc end data
         

         ! Compute and output some checksums for error checking
         call model_write_log("('ua checksum = ',E16.8)", &
                              field_checksum(ua(1:jpiglo,1:jpjglo)) )
         call model_write_log("('va checksum = ',E16.8)", &
                              field_checksum(va(1:jpiglo,1:jpjglo)))

         !! finalise the model run
         CALL finalisation

         ! WRITE(*,*) 'Simulation finished!!'
         ! WRITE(*,*) 
!-----------------------------------
CONTAINS         

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE setup
          INTEGER :: ios

          !! Read in model setup parameters 
          NAMELIST/namctl/ jpiglo, jpjglo, jphgr_msh, &
            &              dx    , dy    , dep_const, &
            &              nit000, nitend, irecord  , &
            &              rdt   , cbfr  , visc

          !! Default value

          jpiglo      =     258               !  number of columns of model grid
          jpjglo      =     258               !  number of rows of model grid
          jphgr_msh   =       1               !  type of grid (0: read in a data file; 1: setup with following parameters)
          dx          =   1000._wp            !  grid size in x direction (m)
          dy          =   1000._wp            !  grid size in y direction (m)
          dep_const   =    100._wp            !  constant depth (m)
          nit000      =       1               !  first time step
          nitend      =    1000               !  end time step
          irecord     =       1               !  intervals to save results
          rdt         =     10._wp            !  size of time step (second) 
          cbfr        =   0.001_wp            !  bottom friction coefficeint
          visc        =     100._wp            !  horizontal kinematic viscosity coefficient 
 
 
          OPEN(1, file='namelist', STATUS='OLD')
          REWIND(1)
          READ(1, NML=namctl, IOSTAT = ios, ERR = 901)
901       IF(ios /= 0) STOP "err found in reading namelist file"
          WRITE(*,NML=namctl)

          ! This correction is to allow the same format namelist file for this
          ! 'original' version of the code and the PSyKAl version. The original
          ! version implicitly adds on the boundaries to the specified domain
          ! whereas in the psykal version the user explicitly specifies the
          ! extent of the domain including boundaries. This is to ensure
          ! consistency with future developments where the T-mask defining the
          ! domain may be read from file.
          jpiglo = jpiglo - 2
          jpjglo = jpjglo - 2

          jpi = jpiglo
          jpj = jpjglo

          CLOSE(1)

        END SUBROUTINE setup

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE grid
           INTEGER  :: ji, jj                          !temporary loop index

        !  need a kernel for boundary /grid 
          !! Allocate working arrays
          !! Define (or read in) model grid

          !jphgr_msh = 0    ! read in this from a namelist file
          !jphgr_msh = 1    ! define manually 
          

          SELECT CASE( jphgr_msh)

          CASE(0) ! read in grid from a coordinate file

            ! to be added
            STOP "It is not ready to Read in grid from a file"
            CALL allocation
            ! add reading data part here
            ! the following variables/arrays are needed:
            ! jpi, jpj, 
            ! pt(0:jpi+1)
            ! e1t(jpi,    jpj), e2t(jpi,    jpj)
            ! e1u(0:jpi,  jpj), e2u(0:jpi,  jpj)
            ! e1v(jpi,  0:jpj), e2v(jpi,  0:jpj)
            ! e1f(0:jpi,0:jpj), e2f(0:jpi,0:jpj)
            ! gphiu(0:jpi,jpj), gphiv(jpi,0:jpj), gphif(0:jpi,0:jpj) 
            ! xt(jpi,jpj), yt(jpi,jpj)
            ! ht(jpi,jpj), hu(jpi,jpj), hv(jpi,jpj)

          CASE(1)

            !##### a manually defined grid

            ! -size of each grid cell
            ! -depth on each T points
            ! -grid dimension


            CALL allocation


            !Define Model solid/open Boundaries via the properties of t-cells
            DO jj = 0, jpj+1
              DO ji = 0, jpi+1
                pt(ji,jj) = 1                             ! all inner celles
              END DO
            END DO

            ! -define solid/open boundaries
            DO jj = 0, jpj+1
              pt(0,    jj) = 0                            ! west solid boundary
              pt(jpi+1,jj) = 0                            ! east solid boundary
            END Do

            DO ji = 0, jpi+1
              pt(ji,jpj+1) = 0                            ! north solid boundary
            END Do

            DO ji = 1, jpi
              pt(ji,0) = -1                               ! south open boundary
            END Do



            !horizontal grid cells

            ! -horizontal scales of grid cells
            ! -Latititude of grid points
            ! -horizontal coordinates of grid points


            e1t(1:jpi, 1:jpj)   = dx
            e2t(1:jpi, 1:jpj)   = dy
            e1u(0:jpi, 1:jpj)   = dx
            e2u(0:jpi, 1:jpj)   = dy
            e1v(1:jpi, 0:jpj)   = dx
            e2v(1:jpi, 0:jpj)   = dy
            e1f(0:jpi, 0:jpj)   = dx
            e2f(0:jpi, 0:jpj)   = dy


            ! -here is a f-plane testing case
            gphiu(0:jpi, 1:jpj) = 50._wp
            gphiv(1:jpi, 0:jpj) = 50._wp
            gphif(0:jpi, 0:jpj) = 50._wp


            xt(1,:) = 0.0_wp + 0.5_wp * e1t(1,:)
            yt(:,1) = 0.0_wp + 0.5_wp * e2t(:,1)

            DO ji = 2, jpi
              xt(ji,1:jpj) = xt(ji-1, 1:jpj) + dx
            END DO
            
            DO jj = 2, jpj
              yt(1:jpi,jj) = yt(1:jpi, jj-1) + dy
            END DO
            !Depth 

            ! -depth on grid points

            DO jj = 1, jpj
            DO ji = 1, jpi
              ht(ji,jj) = dep_const 
            END DO
            END DO

            DO jj = 1, jpj
            DO ji = 0, jpi
              hu(ji,jj) = dep_const 
            END DO
            END DO

            DO jj = 0, jpj
            DO ji = 1, jpi
              hv(ji,jj) = dep_const 
            END DO
            END DO
            
          CASE DEFAULT
            ! undefined jphgr_msh value
            ! add interrupt here
            STOP "Wrong grid defination type, check your setup !!!!"
          END SELECT

          ! calculate t,u,v cell area
          DO jj = 1, jpj
            DO ji = 1, jpi
              e12t(ji,jj) = e1t(ji,jj) * e2t(ji,jj)
            END DO
          END DO

          DO jj = 1, jpj
            DO ji = 0, jpi
              e12u(ji,jj) = e1u(ji,jj) * e2u(ji,jj)
            END DO
          END DO

          DO jj = 0, jpj
            DO ji = 1, jpi
              e12v(ji,jj) = e1v(ji,jj) * e2v(ji,jj)
            END DO
          END DO

        END SUBROUTINE grid

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE allocation
          !! Read in model setup parameters and allocate working arrays
          INTEGER :: ierr(11)
          REAL(wp) :: init_val

          IF (MODULO(jpi, 32) /= 0) THEN
             ALLOCATE(e1t(jpi,jpj), e2t(jpi,jpj), e1u(0:jpi,jpj), e2u(0:jpi,jpj), STAT=ierr(1))
             ALLOCATE(e1f(0:jpi,0:jpj), e2f(0:jpi,0:jpj), e1v(jpi,0:jpj), e2v(jpi,0:jpj), STAT=ierr(2)) 
             ALLOCATE(e12t(jpi,jpj), e12u(0:jpi,jpj), e12v(jpi,0:jpj), STAT=ierr(3))

             ALLOCATE(gphiu(0:jpi,jpj), gphiv(jpi,0:jpj), gphif(0:jpi,0:jpj), STAT=ierr(4))

             ALLOCATE(xt(jpi,jpj), yt(jpi,jpj), STAT=ierr(5))

             ALLOCATE(ht(jpi,jpj), hu(0:jpi,jpj), hv(jpi,0:jpj), hf(0:jpi,0:jpj), STAT=ierr(6))

             ALLOCATE(sshb(jpi,jpj), sshb_u(0:jpi,jpj), sshb_v(jpi,0:jpj), STAT=ierr(7))
             ALLOCATE(sshn(jpi,jpj), sshn_u(0:jpi,jpj), sshn_v(jpi,0:jpj), STAT=ierr(8))
             ALLOCATE(ssha(jpi,jpj), ssha_u(0:jpi,jpj), ssha_v(jpi,0:jpj), STAT=ierr(9))

             ALLOCATE(un(0:jpi,jpj), vn(jpi,0:jpj), ua(0:jpi,jpj), va(jpi,0:jpj), STAT=ierr(10))

             ALLOCATE(pt(0:jpi+1,0:jpj+1), STAT=ierr(11))

          ELSE 
             ! Pad the allocations. TODO: Make this more generic rather than lazy
             ALLOCATE(e1t(jpi,jpj), e2t(jpi,jpj), e1u(-31:jpi,jpj), e2u(-31:jpi,jpj), STAT=ierr(1))
             ALLOCATE(e1f(-31:jpi,0:jpj), e2f(-31:jpi,0:jpj), e1v(jpi,0:jpj), e2v(jpi,0:jpj), STAT=ierr(2)) 
             ALLOCATE(e12t(jpi,jpj), e12u(-31:jpi,jpj), e12v(jpi,0:jpj), STAT=ierr(3))

             ALLOCATE(gphiu(-31:jpi,jpj), gphiv(jpi,0:jpj), gphif(-31:jpi,0:jpj), STAT=ierr(4))

             ALLOCATE(xt(jpi,jpj), yt(jpi,jpj), STAT=ierr(5))

             ALLOCATE(ht(jpi,jpj), hu(-31:jpi,jpj), hv(jpi,0:jpj), hf(-31:jpi,0:jpj), STAT=ierr(6))

             ALLOCATE(sshb(jpi,jpj), sshb_u(-31:jpi,jpj), sshb_v(jpi,0:jpj), STAT=ierr(7))
             ALLOCATE(sshn(jpi,jpj), sshn_u(-31:jpi,jpj), sshn_v(jpi,0:jpj), STAT=ierr(8))
             ALLOCATE(ssha(jpi,jpj), ssha_u(-31:jpi,jpj), ssha_v(jpi,0:jpj), STAT=ierr(9))

             ALLOCATE(un(-31:jpi,jpj), vn(jpi,0:jpj), ua(-31:jpi,jpj), va(jpi,0:jpj), STAT=ierr(10))

             ALLOCATE(pt(0:jpi+1+30,0:jpj+1), STAT=ierr(11))
          END IF
          
          IF(ANY(ierr /= 0, 1)) STOP "in SUBROUTINE ALLOCATION: failed to allocate arrays"

          ! Initialise all solution arrays
          init_val = -999.0
          ht = init_val ; hu = init_val ; hv = init_val ; hf = init_val
          sshb = init_val ; sshb_u = init_val ; sshb_v = init_val
          sshn = init_val ; sshn_u = init_val ; sshn_v = init_val
          ssha = init_val ; ssha_u = init_val ; ssha_v = init_val
          un = init_val ; vn = init_val ; ua = init_val ; va = init_val
          

        END SUBROUTINE allocation


!+++++++++++++++++++++++++++++++++++


        SUBROUTINE initialisation
          INTEGER  :: ji, jj                          !temporary loop index
          REAL(wp) :: rtmp1, rtmp2, rtmp3, rtmp4      !real temporary variables 

          call timer_init()

          ! define (or read in) initil ssh and velocity fields
!         ! split this part into ssh, sshu, sshv, u, v kernels 
            DO ji=1,jpi
              DO jj =1, jpj
                sshn(ji,jj) = 0.0_wp
              END DO
            END DO

            DO ji=0,jpi
              DO jj =1, jpj
                itmp1 = min(ji+1,jpi)
                itmp2 = max(ji,1)
                rtmp1 = e12t(itmp1,jj) * sshn(itmp1,jj) + e12t(itmp2,jj) * sshn(itmp2,jj)
                sshn_u(ji,jj) = 0.5_wp * rtmp1 / e12u(ji,jj)
              END DO
            END DO

            DO ji=1,jpi
              DO jj =0, jpj
                itmp1 = min(jj+1,jpj)
                itmp2 = max(jj,1)
                rtmp1 = e12t(ji,itmp1) * sshn(ji,itmp1) + e12t(ji,itmp2) * sshn(ji,itmp2)
                sshn_v(ji,jj) = 0.5_wp * rtmp1 / e12v(ji,jj)
              END DO
            END DO

            
            DO jj =1, jpj
              DO ji=0,jpi
                 un(ji,jj) = 0._wp
              END DO
            END DO

            DO jj =0, jpj
              DO ji=1,jpi
                 vn(ji,jj) = 0._wp
              END DO
            END DO

! This call updates 'a' quantities which are not used before the code
! in step again updates them. Therefore I think it does nothing. ARP.
!            CALL bc(0._wp)


        END SUBROUTINE initialisation


!+++++++++++++++++++++++++++++++++++


        SUBROUTINE step
          REAL(wp) :: rtime

          rtime = REAL(istp, wp) * rdt

          CALL continuity
          CALL momentum
          CALL bc(rtime)  ! open and solid boundary condition
          CALL next
          IF(MOD(istp, irecord) == 0)  CALL output

        END SUBROUTINE step

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE continuity
          INTEGER  :: ji, jj                          !temporary loop index
          REAL(wp) :: rtmp1, rtmp2, rtmp3, rtmp4      !real temporary variables 
          return
!kernel continuity
          ! Going at 198 GB/s. Fusion?
          !$acc kernels
          ! acc loop collapse(2)
          DO jj = 1, jpj
            DO ji = 1, jpi
              rtmp1 = (sshn_u(ji  ,jj  ) + hu(ji  ,jj  )) * un(ji  ,jj  )
              rtmp2 = (sshn_u(ji-1,jj  ) + hu(ji-1,jj  )) * un(ji-1,jj  )
              rtmp3 = (sshn_v(ji  ,jj  ) + hv(ji  ,jj  )) * vn(ji  ,jj  )
              rtmp4 = (sshn_v(ji  ,jj-1) + hv(ji  ,jj-1)) * vn(ji  ,jj-1)
              ssha(ji,jj) = sshn(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * rdt / e12t(ji,jj)
            END DO
          END DO
          !$acc end kernels
!end kernel continuity
        END SUBROUTINE continuity

!+++++++++++++++++++++++++++++++++++
 
 
 
        SUBROUTINE momentum
          use cuda_mod
          
          TYPE(dim3) gridDim1
          TYPE(dim3) gridDim2
          TYPE(dim3) blockDim1
          TYPE(dim3) blockDim2
          
          integer istat

         INTERFACE
           SUBROUTINE momentum_c(un, vn, sshn, sshn_u, sshn_v, ssha_u, ssha_v, &
                                                   hu, hv, ht, e1v, e2v, e1u, e2u, e1t, e2t, &
                                                   gphiu, gphiv, e12u, e12v, pt, ua, va, &
                                                   ssha, e12t, &
                                                   visc, omega, d2r, g, rdt, cbfr, jpi, jpj) bind(C, name="momentum_c_")
             USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
             use sizes_mod
             IMPLICIT NONE
             REAL(wp), intent(in), device, dimension(:,:) :: un, hu, gphiu, e1u, e2u, sshn_u, ssha_u, e12u
             REAL(wp), intent(in), device, dimension(:,:)    :: vn, hv, gphiv, e1v, e2v, sshn_v, ssha_v, e12v
             REAL(wp), intent(in), device, dimension(:,:)    :: ssha, sshn, ht, e1t, e2t, e12t
             
             
             INTEGER(C_INT), intent(in), device, dimension(:,:) :: pt
             INTEGER(C_INT), value, intent(in) :: jpi, jpj
             REAL(wp), value, intent(in) :: visc, omega, d2r, g, rdt, cbfr
             
             REAL(wp), intent(out), device, dimension(:,:) :: ua
             REAL(wp), intent(out), device, dimension(:,:) :: va
           END SUBROUTINE momentum_c
         END INTERFACE      
          
          blockDim1%x = 128
          blockDim2%x = 128
          blockDim1%y = 1
          blockDim2%y = 1
          
          gridDim1%x = jpi / blockDim1%x
          gridDim2%x = jpi / blockDim2%x
          gridDim1%y = jpj / (blockDim1%y)
          gridDim2%y = jpj / (blockDim2%y)
          
          !$acc host_data use_device(un, vn, sshn, sshn_u, sshn_v, ssha_u, ssha_v, hu, hv, ht, e1v, e2v, e1u, e2u, e1t, e2t, gphiu, gphiv, e12u, e12v, pt, ua, va, ssha, e12t)
          !call momentum_ker <<< gridDim1, blockDim1 >>>(un, vn, sshn, sshn_u, sshn_v, ssha_u, ssha_v, &
          !                                             hu, hv, ht, e1v, e2v, e1u, e2u, e1t, e2t, &
          !                                             gphiu, gphiv, e12u, e12v, pt, ua, va, &
          !                                             visc, omega, d2r, g, rdt, cbfr, jpi, jpj)
          call momentum_acc(un, vn, sshn, sshn_u, sshn_v, ssha_u, ssha_v, &
                            hu, hv, ht, e1v, e2v, e1u, e2u, e1t, e2t, &
                            gphiu, gphiv, e12u, e12v, pt, ua, va, &
                            ssha, e12t, &
                            visc, omega, d2r, g, rdt, cbfr, jpi, jpj)
                                                      
          ! call momentum_c(un, vn, sshn, sshn_u, sshn_v, ssha_u, ssha_v, &
                                                      ! hu, hv, ht, e1v, e2v, e1u, e2u, e1t, e2t, &
                                                      ! gphiu, gphiv, e12u, e12v, pt, ua, va, &
                                                      ! ssha, e12t, &
                                                      ! visc, omega, d2r, g, rdt, cbfr, jpi, jpj)
          !$acc end host_data
          istat = cudaDeviceSynchronize()         
        END SUBROUTINE momentum

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE bc(rtime)
          REAL(wp), INTENT(IN) :: rtime
          REAL(wp) :: amp_tide, omega_tide
          INTEGER :: jiu, jiv
          INTEGER  :: ji, jj                          !temporary loop index

          !open boundary condition of clamped ssh

!kernel ssh clamped obc
            amp_tide   = 0.2_wp
            omega_tide = 2.0_wp * 3.14159_wp / (12.42_wp * 3600._wp)
            !$acc kernels
            ! acc loop collapse(2)
            DO jj = 1, jpj  
              DO ji = 1, jpi 
                IF(pt(ji,jj) <= 0) CYCLE
                IF     (pt(ji,jj-1) < 0) THEN
                  ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
                ELSE IF(pt(ji,jj+1) < 0) THEN
                  ssha(ji,jj) = amp_tide * sin(omega_tide * rtime) 
                ELSE IF(pt(ji+1,jj) < 0) THEN
                  ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
                ELSE IF(pt(ji-1,jj) < 0) THEN
                  ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
                END IF
              END DO
            END DO
!end kernel ssh clamped obc


! kernel"solid boundary conditions for u-velocity" 
            ! !$acc loop collapse(2)
            ! DO jj = 1, jpj
              ! DO ji = 0, jpi
                ! IF(pt(ji,jj) * pt(ji+1,jj) == 0) ua(ji,jj) = 0._wp
              ! END DO
            ! END DO
!end kernel "solid boundary conditions for u-velocity" 

!kernel "solid boundary conditions for v-velocity" 
            ! !$acc loop collapse(2)
            ! DO jj = 0, jpj
              ! DO ji = 1, jpi
                ! IF(pt(ji,jj) * pt(ji,jj+1) == 0) va(ji,jj) = 0._wp
              ! END DO
            ! END DO
!end kernel "solid boundary conditions for v-velocity" 

          !                                            Du                 Dssh
          !start of "Flather open boundary condition [---- = sqrt(g/H) * ------]" Kernel
          !                                            Dn                 Dn
          ! ua and va in du/dn should be the specified tidal forcing


! kernel Flather u 
            DO jj = 1, jpj
              !$acc loop independent
              DO ji = 0, jpi  
                IF(pt(ji,jj) * pt(ji+1,jj) == 0) ua(ji,jj) = 0._wp
                IF(pt(ji,jj) + pt(ji+1,jj) <= -1) CYCLE                         ! not in the domain
                IF(pt(ji,jj) < 0) THEN
                  jiu = ji + 1
                  ua(ji,jj) = ua(jiu,jj) + SQRT(g/hu(ji,jj)) * (sshn_u(ji,jj) - sshn_u(jiu,jj))
                ELSE IF(pt(ji+1,jj )< 0) THEN
                  jiu = ji - 1 
                  ua(ji,jj) = ua(jiu,jj) + SQRT(g/hu(ji,jj)) * (sshn_u(ji,jj) - sshn_u(jiu,jj))
                END IF
              END DO
            END DO
!end kernel flather u .

!kernel Flather v 
            !$acc loop independent
            DO jj = 0, jpj
              DO ji = 1, jpi
                IF(pt(ji,jj) * pt(ji,jj+1) == 0) va(ji,jj) = 0._wp
                IF(pt(ji,jj) + pt(ji,jj+1) <= -1) CYCLE                         ! not in the domain
                IF(pt(ji,jj) < 0) THEN
                  jiv = jj + 1
                  va(ji,jj) = va(ji,jiv) + SQRT(g/hv(ji,jj)) * (sshn_v(ji,jj) - sshn_v(ji,jiv))
                ELSE IF(pt(ji,jj+1) < 0) THEN
                  jiv = jj - 1 
                  va(ji,jj) = va(ji,jiv) + SQRT(g/hv(ji,jj)) * (sshn_v(ji,jj) - sshn_v(ji,jiv))
                END IF
              END DO
            END DO
            !$acc end kernels
!end kernel flather v .

        END SUBROUTINE bc

!+++++++++++++++++++++++++++++++++++


        SUBROUTINE next
          ! update the now-velocity and ssh
          INTEGER  :: ji, jj                          !temporary loop index
          REAL(wp) :: rtmp1, rtmp2, rtmp3, rtmp4      !real temporary variables 
  
! kernel  un updating
          !$acc kernels present(e12t, sshn)
          ! acc loop collapse(2)
          DO jj = 1, jpj
            DO ji = 0, jpi
              un(ji,jj)   = ua(ji,jj)
            END DO
          END DO
! end of kernel sshn_u updating.

! kernel vn updating
          ! acc loop collapse(2)
          DO jj = 0, jpj
            DO ji = 1, jpi
              vn(ji,jj)   = va(ji,jj)
            END DO
          END DO
! end kernel vn updating.

! kernel sshn updating
          ! acc loop collapse(2)
          DO jj = 1, jpj
            DO ji = 1, jpi
              sshn(ji,jj) = ssha(ji,jj)
            END DO
          END DO
! end kernel sshn_u updating.

! kernel sshn_u updating
          ! acc loop collapse(2)
          DO jj = 1, jpj
            DO ji = 0, jpi
              IF(pt(ji,jj) + pt(ji+1,jj) <= 0)  CYCLE                              !jump over non-computational domain
              IF(pt(ji,jj) * pt(ji+1,jj) > 0) THEN
                rtmp1 = e12t(ji,jj) * sshn(ji,jj) + e12t(ji+1,jj) * sshn(ji+1,jj)
                sshn_u(ji,jj) = 0.5_wp * rtmp1 / e12u(ji,jj) 
              ELSE IF(pt(ji,jj) <= 0) THEN
                sshn_u(ji,jj) = sshn(ji+1,jj)
              ELSE IF(pt(ji+1,jj) <= 0) THEN
                sshn_u(ji,jj) = sshn(ji,jj)
              END IF
            END DO
          END DO
! end kernel sshn_u updating.

! kernel: sshn_v updating
          ! acc loop collapse(2)
          DO jj = 0, jpj
            DO ji = 1, jpi
              IF(pt(ji,jj) + pt(ji,jj+1) <= 0)  CYCLE                              !jump over non-computational domain
              IF(pt(ji,jj) * pt(ji,jj+1) > 0) THEN
                rtmp1 = e12t(ji,jj) * sshn(ji,jj) + e12t(ji,jj+1) * sshn(ji,jj+1)
                sshn_v(ji,jj) = 0.5_wp * rtmp1 / e12v(ji,jj) 
              ELSE IF(pt(ji,jj) <= 0) THEN
                sshn_v(ji,jj) = sshn(ji,jj+1)
              ELSE IF(pt(ji,jj+1) <= 0) THEN
                sshn_v(ji,jj) = sshn(ji,jj)
              END If
            END DO
          END DO
          !$acc end kernels
! end kernel sshn_v updating.
            
        END SUBROUTINE next

!+++++++++++++++++++++++++++++++++++


        SUBROUTINE output
          INTEGER  :: ji, jj                          !temporary loop index
          REAL(wp) :: rtmp1, rtmp2, rtmp3, rtmp4      !real temporary variables 

          ! output model results
          CHARACTER(len=5) :: fname
          
          return
          
          WRITE(fname, '(I5.5)') istp
          OPEN(1, file='go2d_'//fname//'.dat', STATUS='UNKNOWN')
          !OPEN(1, file='go2d_'//fname//'.csv', STATUS='UNKNOWN')
          REWIND(1)

          DO jj = 1, jpj
            DO ji = 1, jpi
              rtmp1 = 0.5_wp * (un(ji-1,jj) + un(ji,jj))
              rtmp2 = 0.5_wp * (vn(ji,jj-1) + vn(ji,jj))

              ! write "x-coord, y-coord, depth, ssh, u-velocity, v-velocity" to ASCII files

              !WRITE(1,'(2f20.3, 2f15.4, 2e18.3)')  &            
              !WRITE(1,'(f20.3,'','',f20.3,'','',f15.4,'','',f15.4,'','',f18.3,'','',f18.3)') &
              write(1,'(6e16.7)') &
                   & xt(ji,jj), yt(ji,jj), ht(ji,jj), sshn(ji,jj),rtmp1, rtmp2 
            END DO
            WRITE(1,*)
          END DO
          
          CLOSE(1)

        END SUBROUTINE output

!+++++++++++++++++++++++++++++++++++


        SUBROUTINE finalisation
          INTEGER :: ierr(11)

          DEALLOCATE(e1t, e2t, e1u, e2u, STAT=ierr(1))
          DEALLOCATE(e1f, e2f, e1v, e2v, STAT=ierr(2)) 
          DEALLOCATE(e12t, e12u, e12v, STAT=ierr(3))

          DEALLOCATE(gphiu, gphiv, gphif, STAT=ierr(4))

          DEALLOCATE(xt, yt, STAT=ierr(5))

          DEALLOCATE(ht, hu, hv, hf, STAT=ierr(6))

          DEALLOCATE(sshb, sshb_u, sshb_v, STAT=ierr(7))
          DEALLOCATE(sshn, sshn_u, sshn_v, STAT=ierr(8))
          DEALLOCATE(ssha, ssha_u, ssha_v, STAT=ierr(9))

          DEALLOCATE(un, vn, ua, va, STAT=ierr(10))

          DEALLOCATE(pt, STAT=ierr(11))

          !close opened files
          !send out some ending information
          IF(ANY(ierr /= 0, 1)) STOP "in SUBROUTINE finalisation: failed to deallocate arrays"

         call timer_report()

        END SUBROUTINE finalisation

!+++++++++++++++++++++++++++++++++++

    END PROGRAM nemolite2d
