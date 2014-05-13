    PROGRAM gocean2d
         !!! A Horizontal 2D hydrodynamic ocean model which
         !!   1) using structured grid, but
         !!   2) using in-direct data structures

         IMPLICIT NONE

         INTEGER,  PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)
         INTEGER,  PARAMETER :: dp = SELECTED_REAL_KIND(12, 307)
         INTEGER,  PARAMETER :: wp = dp

         REAL(wp), PARAMETER :: pi    = 3.1415926535897932_wp  
         REAL(wp), PARAMETER :: g     = 9.80665_wp              ! gravity constant
         REAL(wp), PARAMETER :: omega = 7.292116e-05_wp         ! earth rotation speed (s^(-1))
         REAL(wp), PARAMETER :: d2r   = pi / 180._wp            ! degree to radian

         INTEGER,  ALLOCATABLE :: tt_w(:), tt_e(:), tt_n(:), tt_s(:) 
         INTEGER,  ALLOCATABLE :: tu_w(:), tu_e(:), tv_n(:), tv_s(:) 
         INTEGER,  ALLOCATABLE :: ut_w(:), ut_e(:), vt_n(:), vt_s(:) 

         REAL(wp), ALLOCATABLE :: e1t(:), e2t(:), e1u(:), e2u(:) 
         REAL(wp), ALLOCATABLE :: e1f(:), e2f(:), e1v(:), e2v(:) 
         REAL(wp), ALLOCATABLE :: e1e2t(:), e1e2u(:), e1e2v(:)

         REAL(wp), ALLOCATABLE :: gphiu(:), gphiv(:), gphif(:)

         REAL(wp), ALLOCATABLE :: xt(:), yt(:), xf(:), yf(:), ff(:)
         REAL(wp), ALLOCATABLE :: xu(:), yu(:), xv(:), yv(:) 

         REAL(wp), ALLOCATABLE :: ht(:), hu(:), hv(:), hf(:) 

         REAL(wp), ALLOCATABLE :: sshb(:), sshb_u(:), sshb_v(:)
         REAL(wp), ALLOCATABLE :: sshn(:), sshn_u(:), sshn_v(:)
         REAL(wp), ALLOCATABLE :: ssha(:), ssha_u(:), ssha_v(:)

         REAL(wp), ALLOCATABLE :: un(:),  vn(:), ua(:),  va(:)

         INTEGER  :: jpiglo, jpjglo, jpijglot, jpijglof, jpijglou, jpijglov !dimensions of grid
         INTEGER  :: jphgr_msh                                              !type of grid
         INTEGER  :: nit000, nitend, irecord                                !start-end and record time steps

         REAL(wp) :: dx, dy, dep_const                                      !regular grid size and constant depth
         REAL(wp) :: rdt                                                    !time step
                                                                    
         REAL(wp) :: cbfr                                                   !bottom friction coefficient
         REAL(wp) :: visc                                                   !backgroud/constant viscosity 

         INTEGER  :: istp                                                   !time stepping index
         INTEGER  :: ji, jj, jk                                             !temporary loop index

         INTEGER  :: itmp1, itmp2, itmp3, itmp4, itmp5, itmp6               !integer temporary variables
         REAL(wp) :: rtmp1, rtmp2, rtmp3, rtmp4, rtmp5, rtmp6               !real    temporary variables 


         !! read in model parameters
         CALL setup

         !! allocate memory and read in or setup model grid 
         CALL grid

         !! setup model initial condition
         CALL initialisation

         !! time stepping 
         DO istp = nit000, nitend, 1
           CALL step
         END DO

         !! finalise the model run
         CALL finalisation

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

          OPEN(1, file='namelist', STATUS='OLD')
          REWIND(1)
          READ(1, namctl, IOSTAT = ios, ERR = 901)
901       IF(ios /= 0) STOP "err found in reading namelist file"

          CLOSE(1)

        END SUBROUTINE setup

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE grid

          !! Allocate working arrays
          !! Define (or read in) model grid

          !jphgr_msh = 0    ! read in this from a namelist file
          !jphgr_msh = 1    ! define manually 
          

          SELECT CASE( jphgr_msh)

          CASE(0) ! read in grid from a coordinate file

            ! to be added
            STOP "It is not ready to Read in grid from a file"
            CALL allocation

          CASE(1)

            !##### a manually defined grid

            ! -size of each grid cell
            ! -depth on each T points
            ! -grid dimension

            jpijglot = jpiglo * jpjglo
            jpijglof = (jpiglo + 1) * (jpjglo + 1)
            jpijglou = (jpiglo + 1) * jpjglo 
            jpijglov = jpiglo + (jpjglo + 1)

            CALL allocation

            ! -grid topo info

            tt_w(1:jpijglot) = 0        !West  T-cell Neighbour of T-cell
            tt_e(1:jpijglot) = 0        !East  T-cell Neighbour of T-cell
            tt_n(1:jpijglot) = 0        !North T-cell Neighbour of T-cell
            tt_s(1:jpijglot) = 0        !South T-cell Neighbour of T-cell

            tu_w(1:jpijglot) = 0        !West  U-cell Neighbour of T-cell
            tu_e(1:jpijglot) = 0        !East  U-cell Neighbour of T-cell
            tv_n(1:jpijglot) = 0        !North V-cell Neighbour of T-cell
            tv_s(1:jpijglot) = 0        !South V-cell Neighbour of T-cell

            ut_w(1:jpijglou) = 0        !West  T-cell Neighbour of U-cell
            ut_e(1:jpijglou) = 0        !East  T-cell Neighbour of U-cell
            vt_n(1:jpijglov) = 0        !North T-cell Neighbour of V-cell
            vt_s(1:jpijglov) = 0        !South T-cell Neighbour of V-cell

            DO ji = 1, jpijglot
              !West&East T-cell Neighbours of T-cell
              itmp1 = MOD(ji, jpiglo)
              IF(itmp1 /= 1) THEN
                tt_w(ji) = ji - 1
              END IF

              IF(itmp1 /= 0) THEN
                tt_e(ji) = ji + 1
              END IF

              !North&South T-cell Neighbours of T-cell
              IF(ji > jpiglo) THEN
                tt_s(ji) = ji - jpiglo
              END IF

              IF(ji <= jpijglot - jpiglo) THEN
                tt_n(ji) = ji + jpiglo
              END IF

              !West&East U-cell Neighbours of T-cell
              itmp1 = (ji - 1) / jpiglo + 1
              itmp2 = ji + itmp1 
              tu_e(ji) = itmp2
              tu_w(ji) = itmp2 - 1

              ut_w(itmp2)     = ji
              ut_e(itmp2 - 1) = ji


              !North&South V-cell Neighbours of T-cell
              tv_s(ji) = ji 
              tv_n(ji) = ji + jpiglo

              vt_n(ji)          = ji
              vt_s(ji + jpiglo) = ji
            END DO

            !Boundary cell information
            
            ! -solid Boundary cells:   tt_w/e(:) = 0
            !                          tt_n/s(:) = 0
            ! -open Boundary cells:    tt_w/e(:) = -1
            !                          tt_n/s(:) = -1

            ! -manually define a south open boundary (tt_s/n  -1, tt_e/w -1)
            tt_s(1:jpiglo) = -1

            !horizontal grid cells

            ! -horizontal scales of grid cells
            ! -Latititude of grid points
            ! -horizontal coordinates of grid points

            e1t(1:jpijglot)   = dx
            e2t(1:jpijglot)   = dy
            e1u(1:jpijglou)   = dx
            e2u(1:jpijglou)   = dy
            e1v(1:jpijglov)   = dx
            e2v(1:jpijglov)   = dy
            e1f(1:jpijglof)   = dx
            e2f(1:jpijglof)   = dy
            e1e2t(1:jpijglot) = e1t(1:jpijglot) * e2t(1:jpijglot) 
            e1e2u(1:jpijglou) = e1u(1:jpijglou) * e2u(1:jpijglou) 
            e1e2v(1:jpijglov) = e1v(1:jpijglov) * e2v(1:jpijglov) 


            ! -here is a f-plane testing case
            gphiu(1:jpijglou) = 50._wp
            gphiv(1:jpijglov) = 50._wp
            gphif(1:jpijglof) = 50._wp

            DO ji = 1, jpijglof
              xf(ji) = MOD(ji-1, jpiglo+1) * dx
              yf(ji) = ((ji-1) / (jpiglo + 1)  + 1.0_wp) * dy
            END DO

            DO ji = 1, jpijglot
              xt(ji) = (MOD(ji-1,jpiglo) + 0.5_wp) * dx
              yt(ji) = ((ji-1) / jpiglo  + 0.5_wp) * dy
            END DO
            
            DO ji = 1, jpijglou
              xu(ji) = (MOD(ji-1, jpiglo+1) + 0.0_wp) * dx
              yu(ji) = ((ji-1) / (jpiglo+1) + 0.5_wp) * dy
            END DO

            DO ji = 1, jpijglov
              xv(ji) = (MOD(ji-1,jpiglo) + 0.5_wp) * dx
              yv(ji) = ((ji-1) / jpiglo  + 0.0_wp) * dy
            END DO


            !Depth 

            ! -depth on grid points

            DO ji = 1, jpijglof
              hf(ji) = dep_const 
            END DO

            DO ji = 1, jpijglot
              ht(ji) = dep_const
            END DO
            
            DO ji = 1, jpijglou
              hu(ji) = dep_const
            END DO

            DO ji = 1, jpijglov
              hv(ji) = dep_const
            END DO

          CASE DEFAULT
            ! undefined jphgr_msh value
            ! add interrupt here
            STOP "Wrong grid defination type, check your setup !!!!"
          END SELECT

        END SUBROUTINE grid

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE allocation
          !! Read in model setup parameters and allocate working arrays
          INTEGER :: ierr(14)
          
          ALLOCATE(tt_w(jpijglot), tt_e(jpijglot), tt_n(jpijglot), tt_s(jpijglot), STAT=ierr(1)) 
          ALLOCATE(tu_w(jpijglot), tu_e(jpijglot), tv_n(jpijglot), tv_s(jpijglot), STAT=ierr(2))
          ALLOCATE(ut_w(jpijglou), ut_e(jpijglou), vt_n(jpijglov), vt_s(jpijglov), STAT=ierr(3)) 

          ALLOCATE(e1t(jpijglot), e2t(jpijglot), e1u(jpijglou), e2u(jpijglou), STAT=ierr(4))
          ALLOCATE(e1f(jpijglof), e2f(jpijglof), e1v(jpijglov), e2v(jpijglov), STAT=ierr(5)) 
          ALLOCATE(e1e2t(jpijglot), e1e2u(jpijglou), e1e2v(jpijglov), STAT=ierr(6))

          ALLOCATE(gphiu(jpijglou), gphiv(jpijglov), gphif(jpijglot), STAT=ierr(7))

          ALLOCATE(xt(jpijglot), yt(jpijglot), xf(jpijglot), yf(jpijglot), ff(jpijglot), STAT=ierr(8))
          ALLOCATE(xu(jpijglou), yu(jpijglou), xv(jpijglov), yv(jpijglov), STAT=ierr(9))

          ALLOCATE(ht(jpijglot), hu(jpijglou), hv(jpijglov), hf(jpijglot), STAT=ierr(10))

          ALLOCATE(sshb(jpijglot), sshb_u(jpijglou), sshb_v(jpijglov), STAT=ierr(11))
          ALLOCATE(sshn(jpijglot), sshn_u(jpijglou), sshn_v(jpijglov), STAT=ierr(12))
          ALLOCATE(ssha(jpijglot), ssha_u(jpijglou), ssha_v(jpijglov), STAT=ierr(13))

          ALLOCATE(un(jpijglou), vn(jpijglov), ua(jpijglou), va(jpijglov), STAT=ierr(14))

          IF(ANY(ierr /= 0, 1)) STOP "in SUBROUTINE ALLOCATION: failed to allocate arrays"

        END SUBROUTINE allocation


!+++++++++++++++++++++++++++++++++++


        SUBROUTINE initialisation
          ! define (or read in) initil ssh and velocity fields

            DO ji = 1, jpijglot
              sshn(ji) = 0._wp
            END DO

            DO ji = 1, jpijglou
              itmp1 = ut_e(ji)
              itmp2 = ut_w(ji)
              IF(itmp1 * itmp2 > 0) THEN
                rtmp1 = e1e2t(itmp1) * sshn(itmp1) + e1e2t(itmp2) * sshn(itmp2)
                sshn_u(ji) = 0.5_wp * rtmp1 / e1e2u(ji) 
              END IF

              IF(itmp1 <= 0) THEN
                sshn_u(ji) = sshn(itmp2)
              END IF

              IF(itmp2 <=0) THEN
                sshn_u(ji) = sshn(itmp1)
              END IF
            END DO

            DO ji = 1, jpijglov
              itmp1 = vt_n(ji)
              itmp2 = vt_s(ji)
              IF(itmp1 * itmp2 > 0) THEN
                rtmp1 = e1e2t(itmp1) * sshn(itmp1) + e1e2t(itmp2) * sshn(itmp2)
                sshn_v(ji) = 0.5_wp * rtmp1 / e1e2v(ji) 
              END IF

              IF(itmp1 <= 0) THEN
                sshn_v(ji) = sshn(itmp2)
              END IF

              IF(itmp2 <= 0) THEN
                sshn_v(ji) = sshn(itmp1)
              END If
            END DO
            
            DO ji = 1, jpijglou
              un(ji) = 0._wp
            END DO

            DO ji = 1, jpijglov
              vn(ji) = 0._wp
            END DO

            CALL bc(0._wp)


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
          INTEGER :: jie, jiw, jin, jis
          DO ji = 1, jpijglot
             jie = tu_e(ji)
             jiw = tu_w(ji)
             jin = tv_n(ji)
             jis = tv_s(ji)
             rtmp1 = (sshn_u(jie) + hu(jie)) * un(jie)
             rtmp2 = (sshn_u(jiw) + hu(jiw)) * un(jiw)
             rtmp3 = (sshn_v(jin) + hv(jin)) * vn(jin)
             rtmp4 = (sshn_v(jis) + hv(jis)) * vn(jis)
             ssha(ji) = sshn(ji) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * rdt / e1e2t(ji)
          END DO
        END SUBROUTINE continuity

!+++++++++++++++++++++++++++++++++++
 
        SUBROUTINE momentum
          INTEGER  :: jie, jiw, jin, jis
          INTEGER  :: jue, juw, jvn, jvs
          INTEGER  :: jve, jvw
          REAL(wp) :: zadv, zhpg, zcor, zdiff 
          REAL(wp) :: u_e, u_w, u_s, u_n
          REAL(wp) :: v_e, v_w, v_s, v_n
          REAL(wp) :: uu_e, uu_w, uu_s, uu_n
          REAL(wp) :: vv_e, vv_w, vv_s, vv_n
          REAL(wp) :: depe, depw, deps, depn
          REAL(wp) :: dudx_e, dudy_n, dvdx_e, dvdy_n
          REAL(wp) :: dudx_w, dudy_s, dvdx_w, dvdy_s

          REAL(wp) :: adv, vis, hpg, cor, bfr

          ! u equation
          DO ji = 1, jpijglou
             jie = ut_e(ji)
             jiw = ut_w(ji)
             IF( jie <= 0 .OR. jiw <= 0) CYCLE

             jue = tu_e(jie)
             juw = ji
             u_e  = 0.5 * (un(jue) + un(juw)) * e2t(jie)  !add length scale.
             depe = ht(jie) + sshn(jie)

             jue = ji
             juw = tu_w(jiw)
             u_w = 0.5 * (un(jue) + un(juw)) * e2t(jiw)   !add length scale
             depw = ht(jiw) + sshn(jiw)

             IF(ji <= jpiglo) THEN
               v_s  = 0._wp
               deps = 0._wp                    !can be any value as v_s=0
             ELSE
               jvs = tv_s(jiw)
               v_s  = 0.25_wp * (vn(jvs) + vn(jvs+1)) * (e1v(jvs) + e1v(jvs+1))
               deps = 0.5_wp * (hv(jvs) + sshn_v(jvs) + hv(jvs+1) + sshn_v(jvs+1))
             END IF

             IF(ji > jpijglou - jpiglo -1) THEN
               v_n  = 0._wp
               depn = 0._wp
             ELSE
               jvn = tv_n(jiw)
               v_n  = 0.25_wp * (vn(jvn) + vn(jvn+1)) * (e1v(jvn) + e1v(jvn+1))
               depn = 0.5_wp * (hv(jvn) + sshn_v(jvn) + hv(jvn+1) + sshn_v(jvn+1))
             END If

            ! -advection (currently first order upwind)
            uu_w = (0.5_wp - SIGN(0.5_wp, u_w)) * un(ji)              + & 
                 & (0.5_wp + SIGN(0.5_wp, u_w)) * un(ji - 1) 
            uu_e = (0.5_wp + SIGN(0.5_wp, u_e)) * un(ji)              + & 
                 & (0.5_wp - SIGN(0.5_wp, u_e)) * un(ji + 1) 
            uu_s = (0.5_wp - SIGN(0.5_wp, v_s)) * un(ji)              + & 
                 & (0.5_wp + SIGN(0.5_wp, v_s)) * un(ji - jpiglo - 1) 
            uu_n = (0.5_wp + SIGN(0.5_wp, v_n)) * un(ji)              + & 
                 & (0.5_wp - SIGN(0.5_wp, v_n)) * un(ji + jpiglo + 1) 

            adv = uu_w * u_w * depw - uu_e * u_e * depe + uu_s * v_s * deps - uu_n * v_n * depn

            ! -viscosity
            itmp1 = ji + jpiglo + 1
            itmp2 = ji - jpiglo - 1

            dudx_e = (un(ji+1) - un(ji))   / e1t(jie)
            dudx_w = (un(ji)   - un(ji-1)) / e1t(jiw)
            dudy_n = 2.0_wp * (un(itmp1) - un(ji))    / (e2u(ji) + e2u(itmp1))
            dudy_s = 2.0_wp * (un(ji)    - un(itmp2)) / (e2u(ji) + e2u(itmp2))

            vis = (dudx_e * (ht(jie) + sshn(jie)) - dudx_w * (ht(jiw) + sshn(jiw))) * e2u(ji)  + &
                & (dudy_n * (hu(ji) + hu(itmp1) + sshn_u(ji) + sshn_u(itmp1)) -              &
                &  dudy_s * (hu(ji) + hu(itmp2) + sshn_u(ji) + sshn_u(itmp2))) * e1u(ji) * 0.5_wp  
            vis = visc * vis   !visc will be a array visc(1:jpijglou) 
                               !for variable viscosity, such as turbulent viscosity

            ! -Coriolis' force (can be implemented implicitly)
            cor = 0.5_wp * (2._wp * omega * SIN(gphiu(ji) * d2r) * (v_s + v_n)) * &
                & e1u(ji) * e2u(ji) * (hu(ji) + sshn_u(ji))

            ! -pressure gradient
            hpg = -g * (hu(ji) + sshn_u(ji)) * e2u(ji) * (sshn(jie) - sshn(jiw))

            ! -linear bottom friction (implemented implicitly.

            ua(ji) = (un(ji) * (hu(ji) + sshn_u(ji)) + rdt * (adv + vis + cor + hpg) / e1e2u(ji)) / &
                   & (hu(ji) + ssha_u(ji)) / (1.0_wp + cbfr) 

          END DO

          ! v equation
          DO ji = 1, jpijglov
             jis = vt_s(ji)
             jin = vt_n(ji)
             IF( jis <= 0 .OR. jin <= 0) CYCLE

             jvn = tv_n(jin)
             jvs = ji
             v_n  = 0.5 * (vn(jvn) + vn(jvs)) * e1t(jin)  !add length scale.
             depn = ht(jin) + sshn(jin)

             jvn = ji
             jvs = tv_s(jis)
             v_s = 0.5 * (vn(jvn) + vn(jvs)) * e1t(jis)   !add length scale
             deps = ht(jis) + sshn(jis)

             IF(MOD(ji,jpiglo) == 1) THEN
               u_w  = 0._wp
               depw = 0._wp                    !can be any value as v_s=0
             ELSE
               jvw = tu_w(jis)
               u_w = 0.25_wp * (un(jvw) + un(jvw+jpiglo+1)) * &
                   &           (e2u(jvw) + e2u(jvw+jpiglo+1))
               depw = 0.5_wp * (hu(jvw) + sshn_u(jvw) + hu(jvw+jpiglo+1) + sshn_u(jvw+jpiglo+1))
             END IF

             IF(MOD(ji,jpiglo) == 0) THEN
               u_e  = 0._wp
               depe = 0._wp                    !can be any value as v_s=0
             ELSE
               jve = tu_e(jis)
               u_e = 0.25_wp * (un(jve) + un(jve+jpiglo+1)) * &
                   &           (e2t(jve) + e2t(jve+jpiglo+1))
               depe = 0.5_wp * (hu(jve) + sshn_u(jve) + hu(jve+jpiglo+1) + sshn_u(ji+jpiglo+1))
             END IF

            ! -advection (currently first order upwind)
            vv_s = (0.5_wp - SIGN(0.5_wp, v_s)) * vn(ji)              + & 
                 & (0.5_wp + SIGN(0.5_wp, v_s)) * vn(ji - jpiglo) 
            vv_n = (0.5_wp + SIGN(0.5_wp, v_n)) * vn(ji)              + & 
                 & (0.5_wp - SIGN(0.5_wp, v_n)) * vn(ji + jpiglo) 
            vv_w = (0.5_wp - SIGN(0.5_wp, u_w)) * vn(ji)              + & 
                 & (0.5_wp + SIGN(0.5_wp, u_w)) * vn(ji - 1) 
            vv_e = (0.5_wp + SIGN(0.5_wp, u_e)) * vn(ji)              + & 
                 & (0.5_wp - SIGN(0.5_wp, u_e)) * vn(ji + 1) 

            adv = vv_w * u_w * depw - vv_e * u_e * depe + vv_s * v_s * deps - vv_n * v_n * depn

            ! -viscosity
            itmp1 = ji + jpiglo 
            itmp2 = ji - jpiglo

            dvdy_n = (vn(itmp1) - vn(ji))    / e2t(jin)
            dvdy_s = (vn(ji)    - vn(itmp2)) / e2t(jis)
            dvdx_e = 2.0_wp * (vn(ji + 1) - vn(ji))     / (e1v(ji) + e1v(ji + 1))
            dvdx_w = 2.0_wp * (vn(ji)     - vn(ji - 1)) / (e1v(ji) + e1v(ji - 1))

            vis = (dvdy_n * (ht(jin) + sshn(jin)) - dvdy_s * (ht(jis) + sshn(jis))) * e1v(ji)  + &
                & (dvdx_e * (hv(ji) + hv(ji + 1) + sshn_v(ji) + sshn_v(ji + 1)) -              &
                &  dvdx_w * (hv(ji) + hv(ji - 1) + sshn_v(ji) + sshn_v(ji - 1))) * e1u(ji) * 0.5_wp  

            vis = visc * vis   !visc will be a array visc(1:jpijglou) 
                               !for variable viscosity, such as turbulent viscosity


            ! -Coriolis' force (can be implemented implicitly)
            cor = -0.25_wp * (2._wp * omega * SIN(gphiv(ji) * d2r) * (u_e + u_w)) * &
                & e1v(ji) * e2v(ji) * (hv(ji) + sshn_v(ji))

            ! -pressure gradient
            hpg = -g * (hv(ji) + sshn_v(ji)) * e1u(ji) * (sshn(jin) - sshn(jis))

            ! -linear bottom friction (implemented implicitly.

            va(ji) = (vn(ji) * (hv(ji) + sshn_v(ji)) + rdt * (adv + vis + cor + hpg) / e1e2v(ji) ) / &
                   & ((hv(ji) + ssha_v(ji))) / (1.0_wp + cbfr) 

          END DO
        END SUBROUTINE momentum

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE bc(rtime)
          REAL(wp), INTENT(IN) :: rtime
          REAL(wp) :: amp_tide, omega_tide
          INTEGER :: jiu, jiv

          !open boundary condition of clamped ssh
            amp_tide   = 1.0_wp
            omega_tide = 2.0_wp * 3.14159_wp / (12.42_wp * 3600._wp)
            DO ji = 1, jpijglot  ! here can be improved
              IF  (tt_s(ji) < 0) THEN
                ssha(ji)    = amp_tide * sin(omega_tide * rtime)
              ELSE IF(tt_n(ji) < 0) THEN
                ssha(ji) = amp_tide * sin(omega_tide * rtime)
              ELSE IF(tt_e(ji) < 0) THEN
                ssha(ji) = amp_tide * sin(omega_tide * rtime)
              ELSE IF(tt_w(ji) < 0) THEN
                ssha(ji) = amp_tide * sin(omega_tide * rtime)
              END IF
            END DO
            

          !solid boundary conditions for velocity
            DO ji = 1, jpijglou
              IF(ut_e(ji) == 0 .OR. ut_w(ji) == 0) ua(ji) = 0._wp
            END DO

            DO ji = 1, jpijglov
              IF(vt_s(ji) == 0 .OR. vt_n(ji) == 0) va(ji) = 0._wp
            END DO
          !                                  Du                 Dssh
          !Flather open boundary condition [---- = sqrt(g/H) * ------]
          !                                  Dn                 Dn
            DO ji = 1, jpijglou  
              IF  (vt_s(ji) < 0) THEN
                jiv = ji + jpiglo
                va(ji) = va(jiv) + SQRT(g/hv(ji)) * (ssha(ji) - ssha(jiv))
              ELSE IF(vt_n(ji) < 0) THEN
                jiv = ji - jpiglo 
                va(ji) = va(jiv) + SQRT(g/hv(ji)) * (ssha(ji) - ssha(jiv))
              END IF
            END DO

            DO ji = 1, jpijglov  
              IF  (ut_w(ji) < 0) THEN
                jiu = ji + 1
                ua(ji) = ua(jiu) + SQRT(g/hu(ji)) * (ssha(ji) - ssha(jiu))
              ELSE IF(ut_e(ji) < 0) THEN
                jiu = ji - 1
                ua(ji) = ua(jiu) + SQRT(g/hu(ji)) * (ssha(ji) - ssha(jiu))
              END IF
            END DO

        END SUBROUTINE bc

!+++++++++++++++++++++++++++++++++++


        SUBROUTINE next
          ! update the now-velocity and ssh
          sshn(1:jpijglot) = ssha(1:jpijglot)
          un(1:jpijglou)   = ua(1:jpijglou)
          vn(1:jpijglov)   = va(1:jpijglov)

          DO ji = 1, jpijglou
            itmp1 = ut_e(ji)
            itmp2 = ut_w(ji)
            IF(itmp1 * itmp2 > 0) THEN
              rtmp1 = e1e2t(itmp1) * sshn(itmp1) + e1e2t(itmp2) * sshn(itmp2)
              sshn_u(ji) = 0.5_wp * rtmp1 / e1e2u(ji) 
            END IF

            IF(itmp1 <= 0) THEN
              sshn_u(ji) = sshn(itmp2)
            END IF

            IF(itmp2 <=0) THEN
              sshn_u(ji) = sshn(itmp1)
            END IF
          END DO

          DO ji = 1, jpijglov
            itmp1 = vt_n(ji)
            itmp2 = vt_s(ji)
            IF(itmp1 * itmp2 > 0) THEN
              rtmp1 = e1e2t(itmp1) * sshn(itmp1) + e1e2t(itmp2) * sshn(itmp2)
              sshn_v(ji) = 0.5_wp * rtmp1 / e1e2u(ji) 
            END IF

            IF(itmp1 <= 0) THEN
              sshn_v(ji) = sshn(itmp2)
            END IF

            IF(itmp2 <= 0) THEN
              sshn_v(ji) = sshn(itmp1)
            END If
          END DO
            
        END SUBROUTINE next

!+++++++++++++++++++++++++++++++++++


        SUBROUTINE output

          ! output model results
          CHARACTER(len=5) :: fname
          WRITE(fname, '(I4.4)') istp
          OPEN(1, file='go2d_'//fname//'.dat', STATUS='UNKNOWN')
          REWIND(1)

          DO ji = 1, jpijglot

            rtmp1 = 0.5_wp * (un(tu_e(ji)) + un(tu_w(ji)))
            rtmp2 = 0.5_wp * (vn(tv_n(ji)) + vn(tv_s(ji)))

            WRITE(1,'(2f10.3, 2f9.2, 2e8.3)') xt(ji), yt(ji), ht(ji), sshn(ji),rtmp1, rtmp2 

          END DO
          
          CLOSE(1)

        END SUBROUTINE output

!+++++++++++++++++++++++++++++++++++


        SUBROUTINE finalisation
          DEALLOCATE(tt_w, tt_e, tt_n, tt_s) 
          DEALLOCATE(tu_w, tu_e, tv_n, tv_s)
          DEALLOCATE(ut_w, ut_e, vt_n, vt_s) 

          DEALLOCATE(e1t, e2t, e1u, e2u)
          DEALLOCATE(e1f, e2f, e1v, e2v) 
          DEALLOCATE(e1e2t,  e1e2u, e1e2v)

          DEALLOCATE(gphiu,  gphiv, gphif)

          DEALLOCATE(xt,  yt, xf, yf, ff)
          DEALLOCATE(xu,  yu, xv, yv)

          DEALLOCATE(ht,  hu, hv, hf)

          DEALLOCATE(sshb,  sshb_u, sshb_v)
          DEALLOCATE(sshn,  sshn_u, sshn_v)
          DEALLOCATE(ssha,  ssha_u, ssha_v)

          DEALLOCATE(un,  vn, ua,  va)
          !close opened files
          !send out some ending information
          WRITE(*,*) 'Simulation finished!!'

        END SUBROUTINE finalisation

!+++++++++++++++++++++++++++++++++++

    END PROGRAM gocean2D

