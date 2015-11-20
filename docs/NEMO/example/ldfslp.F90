MODULE ldfslp
   !!======================================================================
   !!                       ***  MODULE  ldfslp  ***
   !! Ocean physics: slopes of neutral surfaces
   !!======================================================================
   !! History :  OPA  ! 1994-12  (G. Madec, M. Imbard)  Original code
   !!            8.0  ! 1997-06  (G. Madec)  optimization, lbc
   !!            8.1  ! 1999-10  (A. Jouzeau)  NEW profile in the mixed layer
   !!   NEMO     1.0  ! 2002-10  (G. Madec)  Free form, F90
   !!             -   ! 2005-10  (A. Beckmann)  correction for s-coordinates
   !!            3.3  ! 2010-10  (G. Nurser, C. Harris, G. Madec)  add Griffies operator
   !!             -   ! 2010-11  (F. Dupond, G. Madec)  bug correction in slopes just below the ML
   !!            3.7  ! 2013-12  (F. Lemarie, G. Madec)  add limiter on triad slopes
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_slp       : calculates the slopes of neutral surface   (Madec operator)
   !!   ldf_slp_triad : calculates the triads of isoneutral slopes (Griffies operator)
   !!   ldf_slp_mxl   : calculates the slopes at the base of the mixed layer (Madec operator)
   !!   ldf_slp_init  : initialization of the slopes computation
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE phycst         ! physical constants
   USE zdfmxl         ! mixed layer depth
   USE eosbn2         ! equation of states
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distribued memory computing library
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
   USE wrk_nemo       ! work arrays
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_slp         ! routine called by step.F90
   PUBLIC   ldf_slp_triad   ! routine called by step.F90
   PUBLIC   ldf_slp_init    ! routine called by nemogcm.F90

   LOGICAL , PUBLIC ::   l_ldfslp = .FALSE.     !: slopes flag

   LOGICAL , PUBLIC ::   ln_traldf_iso   = .TRUE.       !: iso-neutral direction
   LOGICAL , PUBLIC ::   ln_traldf_triad = .FALSE.      !: griffies triad scheme

   LOGICAL , PUBLIC ::   ln_triad_iso    = .FALSE.      !: pure horizontal mixing in ML
   LOGICAL , PUBLIC ::   ln_botmix_triad = .FALSE.      !: mixing on bottom
   REAL(wp), PUBLIC ::   rn_sw_triad     = 1._wp        !: =1 switching triads ; =0 all four triads used 
   REAL(wp), PUBLIC ::   rn_slpmax       = 0.01_wp      !: slope limit

   LOGICAL , PUBLIC ::   l_grad_zps = .FALSE.           !: special treatment for Horz Tgradients w partial steps (triad operator)
   
   !                                                     !! Classic operator (Madec)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   uslp, wslpi          !: i_slope at U- and W-points
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   vslp, wslpj          !: j-slope at V- and W-points
   !                                                     !! triad operator (Griffies)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   wslp2                !: wslp**2 from Griffies quarter cells
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::   triadi_g, triadj_g   !: skew flux  slopes relative to geopotentials
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::   triadi  , triadj     !: isoneutral slopes relative to model-coordinate
   !                                                     !! both operators
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   ah_wslp2             !: ah * slope^2 at w-point
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)     ::   akz                  !: stabilizing vertical diffusivity
   
   !                                                     !! Madec operator
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   omlmask           ! mask of the surface mixed layer at T-pt
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   uslpml, wslpiml   ! i_slope at U- and W-points just below the mixed layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   vslpml, wslpjml   ! j_slope at V- and W-points just below the mixed layer

   REAL(wp) ::   repsln = 1.e-25_wp       ! tiny value used as minium of di(rho), dj(rho) and dk(rho)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2014)
   !! $Id: ldfslp.F90 5836 2015-10-26 14:49:40Z cetlod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_slp( kt, prd, pn2 )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE ldf_slp  ***
      !!
      !! ** Purpose :   Compute the slopes of neutral surface (slope of isopycnal
      !!              surfaces referenced locally) (ln_traldf_iso=T).
      !!
      !! ** Method  :   The slope in the i-direction is computed at U- and
      !!      W-points (uslp, wslpi) and the slope in the j-direction is
      !!      computed at V- and W-points (vslp, wslpj).
      !!      They are bounded by 1/100 over the whole ocean, and within the
      !!      surface layer they are bounded by the distance to the surface
      !!      ( slope<= depth/l  where l is the length scale of horizontal
      !!      diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
      !!      of 10cm/s)
      !!        A horizontal shapiro filter is applied to the slopes
      !!        ln_sco=T, s-coordinate, add to the previously computed slopes
      !!      the slope of the model level surface.
      !!        macro-tasked on horizontal slab (jk-loop)  (2, jpk-1)
      !!      [slopes already set to zero at level 1, and to zero or the ocean
      !!      bottom slope (ln_sco=T) at level jpk in inildf]
      !!
      !! ** Action : - uslp, wslpi, and vslp, wslpj, the i- and  j-slopes
      !!               of now neutral surfaces at u-, w- and v- w-points, resp.
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in)                   ::   kt    ! ocean time-step index
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   prd   ! in situ density
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   pn2   ! Brunt-Vaisala frequency (locally ref.)
      !!
      INTEGER  ::   ji , jj , jk    ! dummy loop indices
      INTEGER  ::   ii0, ii1, iku   ! temporary integer
      INTEGER  ::   ij0, ij1, ikv   ! temporary integer
      REAL(wp) ::   zeps, zm1_g, zm1_2g, z1_16, zcofw, z1_slpmax ! local scalars
      REAL(wp) ::   zci, zfi, zau, zbu, zai, zbi   !   -      -
      REAL(wp) ::   zcj, zfj, zav, zbv, zaj, zbj   !   -      -
      REAL(wp) ::   zck, zfk,      zbw             !   -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zwz, zww
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zdzr
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zgru, zgrv
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ldf_slp')
      !
      CALL wrk_alloc( jpi,jpj,jpk, zwz, zww, zdzr, zgru, zgrv )

      zeps   =  1.e-20_wp        !==   Local constant initialization   ==!
      z1_16  =  1.0_wp / 16._wp
      zm1_g  = -1.0_wp / grav
      zm1_2g = -0.5_wp / grav
      z1_slpmax = 1._wp / rn_slpmax
      !
      zww(:,:,:) = 0._wp
      zwz(:,:,:) = 0._wp
      !
      DO jk = 1, jpk             !==   i- & j-gradient of density   ==!
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zgru(ji,jj,jk) = umask(ji,jj,jk) * ( prd(ji+1,jj  ,jk) - prd(ji,jj,jk) )
               zgrv(ji,jj,jk) = vmask(ji,jj,jk) * ( prd(ji  ,jj+1,jk) - prd(ji,jj,jk) )
            END DO
         END DO
      END DO
      IF( ln_zps ) THEN                           ! partial steps correction at the bottom ocean level
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zgru(ji,jj,mbku(ji,jj)) = gru(ji,jj)
               zgrv(ji,jj,mbkv(ji,jj)) = grv(ji,jj)
            END DO
         END DO
      ENDIF
      !
      zdzr(:,:,1) = 0._wp        !==   Local vertical density gradient at T-point   == !   (evaluated from N^2)
      DO jk = 2, jpkm1
         !                                ! zdzr = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point
         !                                !   trick: tmask(ik  )  = 0   =>   all pn2   = 0   =>   zdzr = 0
         !                                !    else  tmask(ik+1)  = 0   =>   pn2(ik+1) = 0   =>   zdzr divides by 1
         !                                !          umask(ik+1) /= 0   =>   all pn2  /= 0   =>   zdzr divides by 2
         !                                ! NB: 1/(tmask+1) = (1-.5*tmask)  substitute a / by a *  ==> faster
         zdzr(:,:,jk) = zm1_g * ( prd(:,:,jk) + 1._wp )              &
            &                 * ( pn2(:,:,jk) + pn2(:,:,jk+1) ) * ( 1._wp - 0.5_wp * tmask(:,:,jk+1) )
      END DO
      !
      !                          !==   Slopes just below the mixed layer   ==!
      CALL ldf_slp_mxl( prd, pn2, zgru, zgrv, zdzr )        ! output: uslpml, vslpml, wslpiml, wslpjml


      ! I.  slopes at u and v point      | uslp = d/di( prd ) / d/dz( prd )
      ! ===========================      | vslp = d/dj( prd ) / d/dz( prd )
      !
      DO jk = 2, jpkm1                            !* Slopes at u and v points
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               !                                      ! horizontal and vertical density gradient at u- and v-points
               zau = zgru(ji,jj,jk) * r1_e1u(ji,jj)
               zav = zgrv(ji,jj,jk) * r1_e2v(ji,jj)
               zbu = 0.5_wp * ( zdzr(ji,jj,jk) + zdzr(ji+1,jj  ,jk) )
               zbv = 0.5_wp * ( zdzr(ji,jj,jk) + zdzr(ji  ,jj+1,jk) )
               !                                      ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
               !                                      ! + kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
               zbu = MIN(  zbu, - z1_slpmax * ABS( zau ) , -7.e+3_wp/fse3u(ji,jj,jk)* ABS( zau )  )
               zbv = MIN(  zbv, - z1_slpmax * ABS( zav ) , -7.e+3_wp/fse3v(ji,jj,jk)* ABS( zav )  )
               !                                      ! uslp and vslp output in zwz and zww, resp.
               zfi = MAX( omlmask(ji,jj,jk), omlmask(ji+1,jj,jk) )
               zfj = MAX( omlmask(ji,jj,jk), omlmask(ji,jj+1,jk) )
               zwz(ji,jj,jk) = ( ( 1. - zfi) * zau / ( zbu - zeps )                                              &
                  &                   + zfi  * uslpml(ji,jj)                                                     &
                  &                          * 0.5_wp * ( fsdept(ji+1,jj,jk)+fsdept(ji,jj,jk)-fse3u(ji,jj,1) )   &
                  &                          / MAX( hmlpt(ji,jj), hmlpt(ji+1,jj), 5._wp ) ) * umask(ji,jj,jk)
               zww(ji,jj,jk) = ( ( 1. - zfj) * zav / ( zbv - zeps )                                              &
                  &                   + zfj  * vslpml(ji,jj)                                                     &
                  &                          * 0.5_wp * ( fsdept(ji,jj+1,jk)+fsdept(ji,jj,jk)-fse3v(ji,jj,1) )   &
                  &                          / MAX( hmlpt(ji,jj), hmlpt(ji,jj+1), 5. ) ) * vmask(ji,jj,jk)
!!gm  modif to suppress omlmask.... (as in Griffies case)
!               !                                         ! jk must be >= ML level for zf=1. otherwise  zf=0.
!               zfi = REAL( 1 - 1/(1 + jk / MAX( nmln(ji+1,jj), nmln(ji,jj) ) ), wp )
!               zfj = REAL( 1 - 1/(1 + jk / MAX( nmln(ji,jj+1), nmln(ji,jj) ) ), wp )
!               zci = 0.5 * ( fsdept(ji+1,jj,jk)+fsdept(ji,jj,jk) ) / MAX( hmlpt(ji,jj), hmlpt(ji+1,jj), 10. ) )
!               zcj = 0.5 * ( fsdept(ji,jj+1,jk)+fsdept(ji,jj,jk) ) / MAX( hmlpt(ji,jj), hmlpt(ji,jj+1), 10. ) )
!               zwz(ji,jj,jk) = ( zfi * zai / ( zbi - zeps ) + ( 1._wp - zfi ) * wslpiml(ji,jj) * zci ) * tmask(ji,jj,jk)
!               zww(ji,jj,jk) = ( zfj * zaj / ( zbj - zeps ) + ( 1._wp - zfj ) * wslpjml(ji,jj) * zcj ) * tmask(ji,jj,jk)
!!gm end modif
            END DO
         END DO
      END DO
      CALL lbc_lnk( zwz, 'U', -1. )   ;   CALL lbc_lnk( zww, 'V', -1. )      ! lateral boundary conditions
      !
      !                                            !* horizontal Shapiro filter
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1, MAX(1, jpj-3)                        ! rows jj=2 and =jpjm1 only
            DO ji = 2, jpim1
               uslp(ji,jj,jk) = z1_16 * (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)      &
                  &                       +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)      &
                  &                       + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)      &
                  &                       +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )    &
                  &                       + 4.*  zwz(ji  ,jj  ,jk)                       )
               vslp(ji,jj,jk) = z1_16 * (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)      &
                  &                       +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)      &
                  &                       + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)      &
                  &                       +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )    &
                  &                       + 4.*  zww(ji,jj    ,jk)                       )
            END DO
         END DO
         DO jj = 3, jpj-2                               ! other rows
            DO ji = fs_2, fs_jpim1   ! vector opt.
               uslp(ji,jj,jk) = z1_16 * (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)      &
                  &                       +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)      &
                  &                       + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)      &
                  &                       +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )    &
                  &                       + 4.*  zwz(ji  ,jj  ,jk)                       )
               vslp(ji,jj,jk) = z1_16 * (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)      &
                  &                       +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)      &
                  &                       + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)      &
                  &                       +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )    &
                  &                       + 4.*  zww(ji,jj    ,jk)                       )
            END DO
         END DO
         !                                        !* decrease along coastal boundaries
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               uslp(ji,jj,jk) = uslp(ji,jj,jk) * ( umask(ji,jj+1,jk) + umask(ji,jj-1,jk  ) ) * 0.5_wp   &
                  &                            * ( umask(ji,jj  ,jk) + umask(ji,jj  ,jk+1) ) * 0.5_wp
               vslp(ji,jj,jk) = vslp(ji,jj,jk) * ( vmask(ji+1,jj,jk) + vmask(ji-1,jj,jk  ) ) * 0.5_wp   &
                  &                            * ( vmask(ji  ,jj,jk) + vmask(ji  ,jj,jk+1) ) * 0.5_wp
            END DO
         END DO
      END DO


      ! II.  slopes at w point           | wslpi = mij( d/di( prd ) / d/dz( prd )
      ! ===========================      | wslpj = mij( d/dj( prd ) / d/dz( prd )
      !
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               !                                  !* Local vertical density gradient evaluated from N^2
               zbw = zm1_2g * pn2 (ji,jj,jk) * ( prd (ji,jj,jk) + prd (ji,jj,jk-1) + 2. )
               !                                  !* Slopes at w point
               !                                        ! i- & j-gradient of density at w-points
               zci = MAX(  umask(ji-1,jj,jk  ) + umask(ji,jj,jk  )           &
                  &      + umask(ji-1,jj,jk-1) + umask(ji,jj,jk-1) , zeps  ) * e1t(ji,jj)
               zcj = MAX(  vmask(ji,jj-1,jk  ) + vmask(ji,jj,jk-1)           &
                  &      + vmask(ji,jj-1,jk-1) + vmask(ji,jj,jk  ) , zeps  ) * e2t(ji,jj)
               zai =    (  zgru (ji-1,jj,jk  ) + zgru (ji,jj,jk-1)           &
                  &      + zgru (ji-1,jj,jk-1) + zgru (ji,jj,jk  )   ) / zci * tmask (ji,jj,jk)
               zaj =    (  zgrv (ji,jj-1,jk  ) + zgrv (ji,jj,jk-1)           &
                  &      + zgrv (ji,jj-1,jk-1) + zgrv (ji,jj,jk  )   ) / zcj * tmask (ji,jj,jk)
               !                                        ! bound the slopes: abs(zw.)<= 1/100 and zb..<0.
               !                                        ! + kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
               zbi = MIN( zbw ,- 100._wp* ABS( zai ) , -7.e+3_wp/fse3w(ji,jj,jk)* ABS( zai )  )
               zbj = MIN( zbw , -100._wp* ABS( zaj ) , -7.e+3_wp/fse3w(ji,jj,jk)* ABS( zaj )  )
               !                                        ! wslpi and wslpj with ML flattening (output in zwz and zww, resp.)
               zfk = MAX( omlmask(ji,jj,jk), omlmask(ji,jj,jk-1) )   ! zfk=1 in the ML otherwise zfk=0
               zck = fsdepw(ji,jj,jk) / MAX( hmlp(ji,jj), 10._wp )
               zwz(ji,jj,jk) = (  zai / ( zbi - zeps ) * ( 1._wp - zfk ) + zck * wslpiml(ji,jj) * zfk  ) * tmask(ji,jj,jk)
               zww(ji,jj,jk) = (  zaj / ( zbj - zeps ) * ( 1._wp - zfk ) + zck * wslpjml(ji,jj) * zfk  ) * tmask(ji,jj,jk)

!!gm  modif to suppress omlmask....  (as in Griffies operator)
!               !                                         ! jk must be >= ML level for zfk=1. otherwise  zfk=0.
!               zfk = REAL( 1 - 1/(1 + jk / nmln(ji+1,jj)), wp )
!               zck = fsdepw(ji,jj,jk)    / MAX( hmlp(ji,jj), 10. )
!               zwz(ji,jj,jk) = ( zfk * zai / ( zbi - zeps ) + ( 1._wp - zfk ) * wslpiml(ji,jj) * zck ) * tmask(ji,jj,jk)
!               zww(ji,jj,jk) = ( zfk * zaj / ( zbj - zeps ) + ( 1._wp - zfk ) * wslpjml(ji,jj) * zck ) * tmask(ji,jj,jk)
!!gm end modif
            END DO
         END DO
      END DO
      CALL lbc_lnk( zwz, 'T', -1. )   ;    CALL lbc_lnk( zww, 'T', -1. )      ! lateral boundary conditions
      !
      !                                           !* horizontal Shapiro filter
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1, MAX(1, jpj-3)                        ! rows jj=2 and =jpjm1 only
            DO ji = 2, jpim1
               zcofw = wmask(ji,jj,jk) * z1_16
               wslpi(ji,jj,jk) = (         zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)     &
                    &               +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)     &
                    &               + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)     &
                    &               +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )   &
                    &               + 4.*  zwz(ji  ,jj  ,jk)                         ) * zcofw

               wslpj(ji,jj,jk) = (         zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)     &
                    &               +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)     &
                    &               + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)     &
                    &               +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )   &
                    &               + 4.*  zww(ji  ,jj  ,jk)                         ) * zcofw
            END DO
         END DO
         DO jj = 3, jpj-2                               ! other rows
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcofw = wmask(ji,jj,jk) * z1_16
               wslpi(ji,jj,jk) = (         zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)     &
                    &               +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)     &
                    &               + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)     &
                    &               +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )   &
                    &               + 4.*  zwz(ji  ,jj  ,jk)                         ) * zcofw

               wslpj(ji,jj,jk) = (         zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)     &
                    &               +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)     &
                    &               + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)     &
                    &               +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )   &
                    &               + 4.*  zww(ji  ,jj  ,jk)                         ) * zcofw
            END DO
         END DO
         !                                        !* decrease in vicinity of topography
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zck =   ( umask(ji,jj,jk) + umask(ji-1,jj,jk) )   &
                  &  * ( vmask(ji,jj,jk) + vmask(ji,jj-1,jk) ) * 0.25
               wslpi(ji,jj,jk) = wslpi(ji,jj,jk) * zck
               wslpj(ji,jj,jk) = wslpj(ji,jj,jk) * zck
            END DO
         END DO
      END DO

      ! IV. Lateral boundary conditions
      ! ===============================
      CALL lbc_lnk( uslp , 'U', -1. )      ;      CALL lbc_lnk( vslp , 'V', -1. )
      CALL lbc_lnk( wslpi, 'W', -1. )      ;      CALL lbc_lnk( wslpj, 'W', -1. )


      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=uslp , clinfo1=' slp  - u : ', tab3d_2=vslp,  clinfo2=' v : ', kdim=jpk)
         CALL prt_ctl(tab3d_1=wslpi, clinfo1=' slp  - wi: ', tab3d_2=wslpj, clinfo2=' wj: ', kdim=jpk)
      ENDIF
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zwz, zww, zdzr, zgru, zgrv )
      !
      IF( nn_timing == 1 )  CALL timing_stop('ldf_slp')
      !
   END SUBROUTINE ldf_slp


   SUBROUTINE ldf_slp_triad ( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE ldf_slp_triad  ***
      !!
      !! ** Purpose :   Compute the squared slopes of neutral surfaces (slope
      !!      of iso-pycnal surfaces referenced locally) (ln_traldf_triad=T)
      !!      at W-points using the Griffies quarter-cells.
      !!
      !! ** Method  :   calculates alpha and beta at T-points
      !!
      !! ** Action : - triadi_g, triadj_g   T-pts i- and j-slope triads relative to geopot. (used for eiv)
      !!             - triadi , triadj    T-pts i- and j-slope triads relative to model-coordinate
      !!             - wslp2              squared slope of neutral surfaces at w-points.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt             ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jl, ip, jp, kp  ! dummy loop indices
      INTEGER  ::   iku, ikv                    ! local integer
      REAL(wp) ::   zfacti, zfactj              ! local scalars
      REAL(wp) ::   znot_thru_surface           ! local scalars
      REAL(wp) ::   zdit, zdis, zdkt, zbu, zbti, zisw
      REAL(wp) ::   zdjt, zdjs, zdks, zbv, zbtj, zjsw
      REAL(wp) ::   zdxrho_raw, zti_coord, zti_raw, zti_lim, zti_g_raw, zti_g_lim
      REAL(wp) ::   zdyrho_raw, ztj_coord, ztj_raw, ztj_lim, ztj_g_raw, ztj_g_lim
      REAL(wp) ::   zdzrho_raw
      REAL(wp) ::   zbeta0, ze3_e1, ze3_e2
      REAL(wp), POINTER, DIMENSION(:,:)     ::   z1_mlbw
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zalbet
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zdxrho , zdyrho, zdzrho     ! Horizontal and vertical density gradients
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zti_mlb, ztj_mlb            ! for Griffies operator only
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ldf_slp_triad')
      !
      CALL wrk_alloc( jpi,jpj, z1_mlbw )
      CALL wrk_alloc( jpi,jpj,jpk, zalbet )
      CALL wrk_alloc( jpi,jpj,jpk,2, zdxrho , zdyrho, zdzrho,              klstart = 0  )
      CALL wrk_alloc( jpi,jpj,  2,2, zti_mlb, ztj_mlb,        kkstart = 0, klstart = 0  )
      !
      !--------------------------------!
      !  Some preliminary calculation  !
      !--------------------------------!
      !
      DO jl = 0, 1                            !==  unmasked before density i- j-, k-gradients  ==!
         !
         ip = jl   ;   jp = jl                ! guaranteed nonzero gradients ( absolute value larger than repsln)
         DO jk = 1, jpkm1                     ! done each pair of triad
            DO jj = 1, jpjm1                  ! NB: not masked ==>  a minimum value is set
               DO ji = 1, fs_jpim1            ! vector opt.
                  zdit = ( tsb(ji+1,jj,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )    ! i-gradient of T & S at u-point
                  zdis = ( tsb(ji+1,jj,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  zdjt = ( tsb(ji,jj+1,jk,jp_tem) - tsb(ji,jj,jk,jp_tem) )    ! j-gradient of T & S at v-point
                  zdjs = ( tsb(ji,jj+1,jk,jp_sal) - tsb(ji,jj,jk,jp_sal) )
                  zdxrho_raw = ( - rab_b(ji+ip,jj   ,jk,jp_tem) * zdit + rab_b(ji+ip,jj   ,jk,jp_sal) * zdis ) * r1_e1u(ji,jj)
                  zdyrho_raw = ( - rab_b(ji   ,jj+jp,jk,jp_tem) * zdjt + rab_b(ji   ,jj+jp,jk,jp_sal) * zdjs ) * r1_e2v(ji,jj)
                  zdxrho(ji+ip,jj   ,jk,1-ip) = SIGN(  MAX( repsln, ABS( zdxrho_raw ) ), zdxrho_raw  )   ! keep the sign
                  zdyrho(ji   ,jj+jp,jk,1-jp) = SIGN(  MAX( repsln, ABS( zdyrho_raw ) ), zdyrho_raw  )
               END DO
            END DO
         END DO
         !
         IF( ln_zps .AND. l_grad_zps ) THEN     ! partial steps: correction of i- & j-grad on bottom
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  iku  = mbku(ji,jj)          ;   ikv  = mbkv(ji,jj)             ! last ocean level (u- & v-points)
                  zdit = gtsu(ji,jj,jp_tem)   ;   zdjt = gtsv(ji,jj,jp_tem)      ! i- & j-gradient of Temperature
                  zdis = gtsu(ji,jj,jp_sal)   ;   zdjs = gtsv(ji,jj,jp_sal)      ! i- & j-gradient of Salinity
                  zdxrho_raw = ( - rab_b(ji+ip,jj   ,iku,jp_tem) * zdit + rab_b(ji+ip,jj   ,iku,jp_sal) * zdis ) * r1_e1u(ji,jj)
                  zdyrho_raw = ( - rab_b(ji   ,jj+jp,ikv,jp_tem) * zdjt + rab_b(ji   ,jj+jp,ikv,jp_sal) * zdjs ) * r1_e2v(ji,jj)
                  zdxrho(ji+ip,jj   ,iku,1-ip) = SIGN( MAX( repsln, ABS( zdxrho_raw ) ), zdxrho_raw )   ! keep the sign
                  zdyrho(ji   ,jj+jp,ikv,1-jp) = SIGN( MAX( repsln, ABS( zdyrho_raw ) ), zdyrho_raw )
               END DO
            END DO
         ENDIF
         !
      END DO

      DO kp = 0, 1                            !==  unmasked before density i- j-, k-gradients  ==!
         DO jk = 1, jpkm1                     ! done each pair of triad
            DO jj = 1, jpj                    ! NB: not masked ==>  a minimum value is set
               DO ji = 1, jpi                 ! vector opt.
                  IF( jk+kp > 1 ) THEN        ! k-gradient of T & S a jk+kp
                     zdkt = ( tsb(ji,jj,jk+kp-1,jp_tem) - tsb(ji,jj,jk+kp,jp_tem) )
                     zdks = ( tsb(ji,jj,jk+kp-1,jp_sal) - tsb(ji,jj,jk+kp,jp_sal) )
                  ELSE
                     zdkt = 0._wp                                             ! 1st level gradient set to zero
                     zdks = 0._wp
                  ENDIF
                  zdzrho_raw = ( - zalbet(ji,jj,jk) * zdkt + zbeta0*zdks ) / fse3w(ji,jj,jk+kp)
                  zdzrho(ji,jj,jk,kp) = - MIN( - repsln , zdzrho_raw )    ! force zdzrho >= repsln
                 END DO
            END DO
         END DO
      END DO
      !
      DO jj = 1, jpj                          !==  Reciprocal depth of the w-point below ML base  ==!
         DO ji = 1, jpi
            jk = MIN( nmln(ji,jj), mbkt(ji,jj) ) + 1     ! MIN in case ML depth is the ocean depth
            z1_mlbw(ji,jj) = 1._wp / fsdepw(ji,jj,jk)
         END DO
      END DO
      !
      !                                       !==  intialisations to zero  ==!
      !
      wslp2  (:,:,:)     = 0._wp              ! wslp2 will be cumulated 3D field set to zero
      triadi_g(:,:,1,:,:) = 0._wp   ;   triadi_g(:,:,jpk,:,:) = 0._wp   ! set surface and bottom slope to zero
      triadj_g(:,:,1,:,:) = 0._wp   ;   triadj_g(:,:,jpk,:,:) = 0._wp
      !!gm _iso set to zero missing
      triadi  (:,:,1,:,:) = 0._wp   ;   triadj  (:,:,jpk,:,:) = 0._wp   ! set surface and bottom slope to zero
      triadj  (:,:,1,:,:) = 0._wp   ;   triadj  (:,:,jpk,:,:) = 0._wp

      !-------------------------------------!
      !  Triads just below the Mixed Layer  !
      !-------------------------------------!
      !
      DO jl = 0, 1                            ! calculate slope of the 4 triads immediately ONE level below mixed-layer base
         DO kp = 0, 1                         ! with only the slope-max limit   and   MASKED
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  ip = jl   ;   jp = jl
                  !
                  jk = nmln(ji+ip,jj) + 1
                  IF( jk > mbkt(ji+ip,jj) ) THEN   ! ML reaches bottom
                     zti_mlb(ji+ip,jj   ,1-ip,kp) = 0.0_wp
                  ELSE                             
                     ! Add s-coordinate slope at t-points (do this by *subtracting* gradient of depth)
                     zti_g_raw = (  zdxrho(ji+ip,jj,jk-kp,1-ip) / zdzrho(ji+ip,jj,jk-kp,kp)      &
                        &          - ( fsdept(ji+1,jj,jk-kp) - fsdept(ji,jj,jk-kp) ) * r1_e1u(ji,jj)  ) * umask(ji,jj,jk)
                     ze3_e1    =  fse3w(ji+ip,jj,jk-kp) * r1_e1u(ji,jj) 
                     zti_mlb(ji+ip,jj   ,1-ip,kp) = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1  , ABS( zti_g_raw ) ), zti_g_raw )
                  ENDIF
                  !
                  jk = nmln(ji,jj+jp) + 1
                  IF( jk .GT. mbkt(ji,jj+jp) ) THEN  !ML reaches bottom
                     ztj_mlb(ji   ,jj+jp,1-jp,kp) = 0.0_wp
                  ELSE
                     ztj_g_raw = (  zdyrho(ji,jj+jp,jk-kp,1-jp) / zdzrho(ji,jj+jp,jk-kp,kp)      &
                        &      - ( fsdept(ji,jj+1,jk-kp) - fsdept(ji,jj,jk-kp) ) / e2v(ji,jj)  ) * vmask(ji,jj,jk)
                     ze3_e2    =  fse3w(ji,jj+jp,jk-kp) / e2v(ji,jj)
                     ztj_mlb(ji   ,jj+jp,1-jp,kp) = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2  , ABS( ztj_g_raw ) ), ztj_g_raw )
                  ENDIF
               END DO
            END DO
         END DO
      END DO

      !-------------------------------------!
      !  Triads with surface limits         !
      !-------------------------------------!
      !
      DO kp = 0, 1                            ! k-index of triads
         DO jl = 0, 1
            ip = jl   ;   jp = jl             ! i- and j-indices of triads (i-k and j-k planes)
            DO jk = 1, jpkm1
               ! Must mask contribution to slope from dz/dx at constant s for triads jk=1,kp=0 that poke up though ocean surface
               znot_thru_surface = REAL( 1-1/(jk+kp), wp )  !jk+kp=1,=0.; otherwise=1.0
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1         ! vector opt.
                     !
                     ! Calculate slope relative to geopotentials used for GM skew fluxes
                     ! Add s-coordinate slope at t-points (do this by *subtracting* gradient of depth)
                     ! Limit by slope *relative to geopotentials* by rn_slpmax, and mask by psi-point
                     ! masked by umask taken at the level of dz(rho)
                     !
                     ! raw slopes: unmasked unbounded slopes (relative to geopotential (zti_g) and model surface (zti)
                     !
                     zti_raw   = zdxrho(ji+ip,jj   ,jk,1-ip) / zdzrho(ji+ip,jj   ,jk,kp)                   ! unmasked
                     ztj_raw   = zdyrho(ji   ,jj+jp,jk,1-jp) / zdzrho(ji   ,jj+jp,jk,kp)
                     !
                     ! Must mask contribution to slope for triad jk=1,kp=0 that poke up though ocean surface
                     zti_coord = znot_thru_surface * ( fsdept(ji+1,jj  ,jk) - fsdept(ji,jj,jk) ) * r1_e1u(ji,jj)
                     ztj_coord = znot_thru_surface * ( fsdept(ji  ,jj+1,jk) - fsdept(ji,jj,jk) ) * r1_e2v(ji,jj)     ! unmasked
                     zti_g_raw = zti_raw - zti_coord      ! ref to geopot surfaces
                     ztj_g_raw = ztj_raw - ztj_coord
                     ! additional limit required in bilaplacian case
                     ze3_e1    = fse3w(ji+ip,jj   ,jk+kp) * r1_e1u(ji,jj)
                     ze3_e2    = fse3w(ji   ,jj+jp,jk+kp) * r1_e2v(ji,jj)
                     ! NB: hard coded factor 5 (can be a namelist parameter...)
                     zti_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e1, ABS( zti_g_raw ) ), zti_g_raw )
                     ztj_g_lim = SIGN( MIN( rn_slpmax, 5.0_wp * ze3_e2, ABS( ztj_g_raw ) ), ztj_g_raw )
                     !
                     ! Below  ML use limited zti_g as is & mask
                     ! Inside ML replace by linearly reducing sx_mlb towards surface & mask
                     !
                     zfacti = REAL( 1 - 1/(1 + (jk+kp-1)/nmln(ji+ip,jj)), wp )  ! k index of uppermost point(s) of triad is jk+kp-1
                     zfactj = REAL( 1 - 1/(1 + (jk+kp-1)/nmln(ji,jj+jp)), wp )  ! must be .ge. nmln(ji,jj) for zfact=1
                     !                                                          !                   otherwise  zfact=0
                     zti_g_lim =          ( zfacti   * zti_g_lim                       &
                        &      + ( 1._wp - zfacti ) * zti_mlb(ji+ip,jj,1-ip,kp)   &
                        &                           * fsdepw(ji+ip,jj,jk+kp) * z1_mlbw(ji+ip,jj) ) * umask(ji,jj,jk+kp)
                     ztj_g_lim =          ( zfactj   * ztj_g_lim                       &
                        &      + ( 1._wp - zfactj ) * ztj_mlb(ji,jj+jp,1-jp,kp)   &
                        &                           * fsdepw(ji,jj+jp,jk+kp) * z1_mlbw(ji,jj+jp) ) * vmask(ji,jj,jk+kp)
                     !
                     triadi_g(ji+ip,jj   ,jk,1-ip,kp) = zti_g_lim
                     triadj_g(ji   ,jj+jp,jk,1-jp,kp) = ztj_g_lim
                     !
                     ! Get coefficients of isoneutral diffusion tensor
                     ! 1. Utilise gradients *relative* to s-coordinate, so add t-point slopes (*subtract* depth gradients)
                     ! 2. We require that isoneutral diffusion  gives no vertical buoyancy flux
                     !     i.e. 33 term = (real slope* 31, 13 terms)
                     ! To do this, retain limited sx**2  in vertical flux, but divide by real slope for 13/31 terms
                     ! Equivalent to tapering A_iso = sx_limited**2/(real slope)**2
                     !
                     zti_lim  = ( zti_g_lim + zti_coord ) * umask(ji,jj,jk+kp)    ! remove coordinate slope => relative to coordinate surfaces
                     ztj_lim  = ( ztj_g_lim + ztj_coord ) * vmask(ji,jj,jk+kp)
                     !
                     IF( ln_triad_iso ) THEN
                        zti_raw = zti_lim*zti_lim / zti_raw
                        ztj_raw = ztj_lim*ztj_lim / ztj_raw
                        zti_raw = SIGN( MIN( ABS(zti_lim), ABS( zti_raw ) ), zti_raw )
                        ztj_raw = SIGN( MIN( ABS(ztj_lim), ABS( ztj_raw ) ), ztj_raw )
                        zti_lim = zfacti * zti_lim + ( 1._wp - zfacti ) * zti_raw
                        ztj_lim = zfactj * ztj_lim + ( 1._wp - zfactj ) * ztj_raw
                     ENDIF
                     !                                      ! switching triad scheme 
                     zisw = (rn_sw_triad - 1._wp ) + rn_sw_triad    &
                        &            * 2._wp * ABS( 0.5_wp - kp - ( 0.5_wp - ip ) * SIGN( 1._wp , zdxrho(ji+ip,jj,jk,1-ip) )  )
                     zjsw = (rn_sw_triad - 1._wp ) + rn_sw_triad    &
                        &            * 2._wp * ABS( 0.5_wp - kp - ( 0.5_wp - jp ) * SIGN( 1._wp , zdyrho(ji,jj+jp,jk,1-jp) )  )
                     !
                     triadi(ji+ip,jj   ,jk,1-ip,kp) = zti_lim * zisw
                     triadj(ji   ,jj+jp,jk,1-jp,kp) = ztj_lim * zjsw
                     !
                     zbu  = e1e2u(ji   ,jj   ) * fse3u(ji   ,jj   ,jk   )
                     zbv  = e1e2v(ji   ,jj   ) * fse3v(ji   ,jj   ,jk   )
                     zbti = e1e2t(ji+ip,jj   ) * fse3w(ji+ip,jj   ,jk+kp)
                     zbtj = e1e2t(ji   ,jj+jp) * fse3w(ji   ,jj+jp,jk+kp)
                     !
                     wslp2(ji+ip,jj,jk+kp) = wslp2(ji+ip,jj,jk+kp) + 0.25_wp * zbu / zbti * zti_g_lim*zti_g_lim      ! masked
                     wslp2(ji,jj+jp,jk+kp) = wslp2(ji,jj+jp,jk+kp) + 0.25_wp * zbv / zbtj * ztj_g_lim*ztj_g_lim
                  END DO
               END DO
            END DO
         END DO
      END DO
      !
      wslp2(:,:,1) = 0._wp                ! force the surface wslp to zero

      CALL lbc_lnk( wslp2, 'W', 1. )      ! lateral boundary confition on wslp2 only   ==>>> gm : necessary ? to be checked
      !
      CALL wrk_dealloc( jpi,jpj, z1_mlbw )
      CALL wrk_dealloc( jpi,jpj,jpk, zalbet )
      CALL wrk_dealloc( jpi,jpj,jpk,2, zdxrho , zdyrho, zdzrho,              klstart = 0  )
      CALL wrk_dealloc( jpi,jpj,  2,2, zti_mlb, ztj_mlb,        kkstart = 0, klstart = 0  )
      !
      IF( nn_timing == 1 )  CALL timing_stop('ldf_slp_triad')
      !
   END SUBROUTINE ldf_slp_triad


   SUBROUTINE ldf_slp_mxl( prd, pn2, p_gru, p_grv, p_dzr )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_slp_mxl  ***
      !!
      !! ** Purpose :   Compute the slopes of iso-neutral surface just below
      !!              the mixed layer.
      !!
      !! ** Method  :   The slope in the i-direction is computed at u- & w-points
      !!              (uslpml, wslpiml) and the slope in the j-direction is computed
      !!              at v- and w-points (vslpml, wslpjml) with the same bounds as
      !!              in ldf_slp.
      !!
      !! ** Action  :   uslpml, wslpiml :  i- &  j-slopes of neutral surfaces
      !!                vslpml, wslpjml    just below the mixed layer
      !!                omlmask         :  mixed layer mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   prd            ! in situ density
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   pn2            ! Brunt-Vaisala frequency (locally ref.)
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   p_gru, p_grv   ! i- & j-gradient of density (u- & v-pts)
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   p_dzr          ! z-gradient of density      (T-point)
      !!
      INTEGER  ::   ji , jj , jk                   ! dummy loop indices
      INTEGER  ::   iku, ikv, ik, ikm1             ! local integers
      REAL(wp) ::   zeps, zm1_g, zm1_2g, z1_slpmax ! local scalars
      REAL(wp) ::   zci, zfi, zau, zbu, zai, zbi   !   -      -
      REAL(wp) ::   zcj, zfj, zav, zbv, zaj, zbj   !   -      -
      REAL(wp) ::   zck, zfk,      zbw             !   -      -
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ldf_slp_mxl')
      !
      zeps   =  1.e-20_wp        !==   Local constant initialization   ==!
      zm1_g  = -1.0_wp / grav
      zm1_2g = -0.5_wp / grav
      z1_slpmax = 1._wp / rn_slpmax
      !
      uslpml (1,:) = 0._wp      ;      uslpml (jpi,:) = 0._wp
      vslpml (1,:) = 0._wp      ;      vslpml (jpi,:) = 0._wp
      wslpiml(1,:) = 0._wp      ;      wslpiml(jpi,:) = 0._wp
      wslpjml(1,:) = 0._wp      ;      wslpjml(jpi,:) = 0._wp
      !
      !                                            !==   surface mixed layer mask   !
      DO jk = 1, jpk                               ! =1 inside the mixed layer, =0 otherwise
         DO jj = 1, jpj
            DO ji = 1, jpi
               ik = nmln(ji,jj) - 1
               IF( jk <= ik ) THEN   ;   omlmask(ji,jj,jk) = 1._wp
               ELSE                  ;   omlmask(ji,jj,jk) = 0._wp
               ENDIF
            END DO
         END DO
      END DO


      ! Slopes of isopycnal surfaces just before bottom of mixed layer
      ! --------------------------------------------------------------
      ! The slope are computed as in the 3D case.
      ! A key point here is the definition of the mixed layer at u- and v-points.
      ! It is assumed to be the maximum of the two neighbouring T-point mixed layer depth.
      ! Otherwise, a n2 value inside the mixed layer can be involved in the computation
      ! of the slope, resulting in a too steep diagnosed slope and thus a spurious eddy
      ! induce velocity field near the base of the mixed layer.
      !-----------------------------------------------------------------------
      !
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            !                        !==   Slope at u- & v-points just below the Mixed Layer   ==!
            !
            !                        !- vertical density gradient for u- and v-slopes (from dzr at T-point)
            iku = MIN(  MAX( 1, nmln(ji,jj) , nmln(ji+1,jj) ) , jpkm1  )   ! ML (MAX of T-pts, bound by jpkm1)
            ikv = MIN(  MAX( 1, nmln(ji,jj) , nmln(ji,jj+1) ) , jpkm1  )   !
            zbu = 0.5_wp * ( p_dzr(ji,jj,iku) + p_dzr(ji+1,jj  ,iku) )
            zbv = 0.5_wp * ( p_dzr(ji,jj,ikv) + p_dzr(ji  ,jj+1,ikv) )
            !                        !- horizontal density gradient at u- & v-points
            zau = p_gru(ji,jj,iku) * r1_e1u(ji,jj)
            zav = p_grv(ji,jj,ikv) * r1_e2v(ji,jj)
            !                        !- bound the slopes: abs(zw.)<= 1/100 and zb..<0
            !                           kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
            zbu = MIN(  zbu , - z1_slpmax * ABS( zau ) , -7.e+3_wp/fse3u(ji,jj,iku)* ABS( zau )  )
            zbv = MIN(  zbv , - z1_slpmax * ABS( zav ) , -7.e+3_wp/fse3v(ji,jj,ikv)* ABS( zav )  )
            !                        !- Slope at u- & v-points (uslpml, vslpml)
            uslpml(ji,jj) = zau / ( zbu - zeps ) * umask(ji,jj,iku)
            vslpml(ji,jj) = zav / ( zbv - zeps ) * vmask(ji,jj,ikv)
            !
            !                        !==   i- & j-slopes at w-points just below the Mixed Layer   ==!
            !
            ik   = MIN( nmln(ji,jj) + 1, jpk )
            ikm1 = MAX( 1, ik-1 )
            !                        !- vertical density gradient for w-slope (from N^2)
            zbw = zm1_2g * pn2 (ji,jj,ik) * ( prd (ji,jj,ik) + prd (ji,jj,ikm1) + 2. )
            !                        !- horizontal density i- & j-gradient at w-points
            zci = MAX(   umask(ji-1,jj,ik  ) + umask(ji,jj,ik  )           &
               &       + umask(ji-1,jj,ikm1) + umask(ji,jj,ikm1) , zeps  ) * e1t(ji,jj)
            zcj = MAX(   vmask(ji,jj-1,ik  ) + vmask(ji,jj,ik  )           &
               &       + vmask(ji,jj-1,ikm1) + vmask(ji,jj,ikm1) , zeps  ) * e2t(ji,jj)
            zai =    (   p_gru(ji-1,jj,ik  ) + p_gru(ji,jj,ik)           &
               &       + p_gru(ji-1,jj,ikm1) + p_gru(ji,jj,ikm1  )  ) / zci  * tmask(ji,jj,ik)
            zaj =    (   p_grv(ji,jj-1,ik  ) + p_grv(ji,jj,ik  )           &
               &       + p_grv(ji,jj-1,ikm1) + p_grv(ji,jj,ikm1)  ) / zcj  * tmask(ji,jj,ik)
            !                        !- bound the slopes: abs(zw.)<= 1/100 and zb..<0.
            !                           kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
            zbi = MIN(  zbw , -100._wp* ABS( zai ) , -7.e+3_wp/fse3w(ji,jj,ik)* ABS( zai )  )
            zbj = MIN(  zbw , -100._wp* ABS( zaj ) , -7.e+3_wp/fse3w(ji,jj,ik)* ABS( zaj )  )
            !                        !- i- & j-slope at w-points (wslpiml, wslpjml)
            wslpiml(ji,jj) = zai / ( zbi - zeps ) * tmask (ji,jj,ik)
            wslpjml(ji,jj) = zaj / ( zbj - zeps ) * tmask (ji,jj,ik)
         END DO
      END DO
      !!gm this lbc_lnk should be useless....
      CALL lbc_lnk( uslpml , 'U', -1. )   ;   CALL lbc_lnk( vslpml , 'V', -1. )   ! lateral boundary cond. (sign change)
      CALL lbc_lnk( wslpiml, 'W', -1. )   ;   CALL lbc_lnk( wslpjml, 'W', -1. )   ! lateral boundary conditions
      !
      IF( nn_timing == 1 )  CALL timing_stop('ldf_slp_mxl')
      !
   END SUBROUTINE ldf_slp_mxl


   SUBROUTINE ldf_slp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_slp_init  ***
      !!
      !! ** Purpose :   Initialization for the isopycnal slopes computation
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      INTEGER ::   ierr         ! local integer
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ldf_slp_init')
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_slp_init : direction of lateral mixing'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      ALLOCATE( ah_wslp2(jpi,jpj,jpk) , akz(jpi,jpj,jpk) , STAT=ierr )
      IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'ldf_slp_init : unable to allocate ah_slp2 or akz' )
      !
      IF( ln_traldf_triad ) THEN        ! Griffies operator : triad of slopes
         IF(lwp) WRITE(numout,*) '              Griffies (triad) operator initialisation'
         ALLOCATE( triadi_g(jpi,jpj,jpk,0:1,0:1) , triadj_g(jpi,jpj,jpk,0:1,0:1) ,     &
            &      triadi  (jpi,jpj,jpk,0:1,0:1) , triadj  (jpi,jpj,jpk,0:1,0:1) ,     &
            &      wslp2   (jpi,jpj,jpk)                                         , STAT=ierr )
         IF( ierr > 0      )   CALL ctl_stop( 'STOP', 'ldf_slp_init : unable to allocate Griffies operator slope' )
         IF( ln_dynldf_iso )   CALL ctl_stop( 'ldf_slp_init: Griffies operator on momentum not supported' )
         !
      ELSE                             ! Madec operator : slopes at u-, v-, and w-points
         IF(lwp) WRITE(numout,*) '              Madec operator initialisation'
         ALLOCATE( omlmask(jpi,jpj,jpk) ,                                                                        &
            &      uslp(jpi,jpj,jpk) , uslpml(jpi,jpj) , wslpi(jpi,jpj,jpk) , wslpiml(jpi,jpj) ,     &
            &      vslp(jpi,jpj,jpk) , vslpml(jpi,jpj) , wslpj(jpi,jpj,jpk) , wslpjml(jpi,jpj) , STAT=ierr )
         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'ldf_slp_init : unable to allocate Madec operator slope ' )

         ! Direction of lateral diffusion (tracers and/or momentum)
         ! ------------------------------
         uslp (:,:,:) = 0._wp   ;   uslpml (:,:) = 0._wp      ! set the slope to zero (even in s-coordinates)
         vslp (:,:,:) = 0._wp   ;   vslpml (:,:) = 0._wp
         wslpi(:,:,:) = 0._wp   ;   wslpiml(:,:) = 0._wp
         wslpj(:,:,:) = 0._wp   ;   wslpjml(:,:) = 0._wp

         !!gm I no longer understand this.....
!!gm         IF( (ln_traldf_hor .OR. ln_dynldf_hor) .AND. .NOT. (lk_vvl .AND. ln_rstart) ) THEN
!            IF(lwp)   WRITE(numout,*) '          Horizontal mixing in s-coordinate: slope = slope of s-surfaces'
!
!            ! geopotential diffusion in s-coordinates on tracers and/or momentum
!            ! The slopes of s-surfaces are computed once (no call to ldfslp in step)
!            ! The slopes for momentum diffusion are i- or j- averaged of those on tracers
!
!            ! set the slope of diffusion to the slope of s-surfaces
!            !      ( c a u t i o n : minus sign as fsdep has positive value )
!            DO jk = 1, jpk
!               DO jj = 2, jpjm1
!                  DO ji = fs_2, fs_jpim1   ! vector opt.
!                     uslp (ji,jj,jk) = - ( fsdept(ji+1,jj,jk) - fsdept(ji ,jj ,jk) ) * r1_e1u(ji,jj) * umask(ji,jj,jk)
!                     vslp (ji,jj,jk) = - ( fsdept(ji,jj+1,jk) - fsdept(ji ,jj ,jk) ) * r1_e2v(ji,jj) * vmask(ji,jj,jk)
!                     wslpi(ji,jj,jk) = - ( fsdepw(ji+1,jj,jk) - fsdepw(ji-1,jj,jk) ) * r1_e1t(ji,jj) * wmask(ji,jj,jk) * 0.5
!                     wslpj(ji,jj,jk) = - ( fsdepw(ji,jj+1,jk) - fsdepw(ji,jj-1,jk) ) * r1_e2t(ji,jj) * wmask(ji,jj,jk) * 0.5
!                  END DO
!               END DO
!            END DO
!            CALL lbc_lnk( uslp , 'U', -1. )   ;   CALL lbc_lnk( vslp , 'V', -1. )      ! Lateral boundary conditions
!            CALL lbc_lnk( wslpi, 'W', -1. )   ;   CALL lbc_lnk( wslpj, 'W', -1. )
!!gm         ENDIF
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('ldf_slp_init')
      !
   END SUBROUTINE ldf_slp_init

   !!======================================================================
END MODULE ldfslp
