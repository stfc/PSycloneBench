!> Module containing a subroutine that implements the computation performed by the
!  tracer-advection benchmark. This is done to ease the conversion to SIR (which
!  doesn't support the various bits of general-purpose code needed to set-up and
!  shutdown the benchmark program).
module tra_adv_compute_mod

  implicit none

contains

  subroutine tra_adv_compute_01_fort(zind, tsn, ztfreez, rnfmsk, rnfmsk_z, upsmsk, tmask, jpi, jpj, jpk, iter)   
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout):: zind
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in)   :: tsn, tmask
      REAL*8, ALLOCATABLE, DIMENSION(:,:), intent(in) :: ztfreez, rnfmsk, upsmsk
      REAL*8, ALLOCATABLE, DIMENSION(:), intent(in) :: rnfmsk_z

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter

      REAL*8 :: zice
      INTEGER :: ji, jj, jk

      CHARACTER(len=64) :: fieldname     

      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zice = 0.d0
               IF( tsn(ji,jj,jk) <= ztfreez(ji,jj) + 0.1d0 ) THEN   ;   zice = 1.d0
               ELSE                                                 ;   zice = 0.d0
               ENDIF

               zind(ji,jj,jk) = MAX (   &
                  rnfmsk(ji,jj) * rnfmsk_z(jk),      & 
                  upsmsk(ji,jj)               ,      &
                  zice                               &
                  &                  ) * tmask(ji,jj,jk)
               zind(ji,jj,jk) = 1 - zind(ji,jj,jk)           
            END DO
         END DO
      END DO

    end subroutine tra_adv_compute_01_fort

    subroutine tra_adv_compute_02_fort(zwx, zwy, umask, vmask, mydomain, jpi, jpj, jpk, iter)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout) :: zwx, zwy
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: umask, vmask, mydomain

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter

      INTEGER :: ji, jj, jk     

      zwx(:,:,jpk) = 0.e0   ;   zwy(:,:,jpk) = 0.e0

      DO jk = 1, jpk-1
         DO jj = 1, jpj-1
            DO ji = 1, jpi-1
               zwx(ji,jj,jk) = umask(ji,jj,jk) * ( mydomain(ji+1,jj,jk) - mydomain(ji,jj,jk) )
               zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( mydomain(ji,jj+1,jk) - mydomain(ji,jj,jk) )
            END DO
         END DO
      END DO     

   end subroutine tra_adv_compute_02_fort

   subroutine tra_adv_compute_03_fort(zslpx, zslpy, zwx, zwy, jpi, jpj, jpk, iter)   
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout) :: zslpx, zslpy
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: zwx, zwy

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter

      INTEGER :: ji, jj, jk

      zslpx(:,:,jpk) = 0.e0   ;   zslpy(:,:,jpk) = 0.e0
      
      DO jk = 1, jpk-1
         DO jj = 2, jpj
            DO ji = 2, jpi 
               zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji-1,jj  ,jk) )   &
                  &            * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji-1,jj  ,jk) ) )
               zslpy(ji,jj,jk) =                    ( zwy(ji,jj,jk) + zwy(ji  ,jj-1,jk) )   &
                  &            * ( 0.25d0 + SIGN( 0.25d0, zwy(ji,jj,jk) * zwy(ji  ,jj-1,jk) ) )
            END DO
         END DO
      END DO

   end subroutine tra_adv_compute_03_fort

   subroutine tra_adv_compute_04_fort(zslpx, zslpy, zwx, zwy, jpi, jpj, jpk, iter)   
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout) :: zslpx, zslpy
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: zwx, zwy

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter

      INTEGER :: ji, jj, jk

      DO jk = 1, jpk-1    
         DO jj = 2, jpj
            DO ji = 2, jpi
               zslpx(ji,jj,jk) = SIGN( 1.d0, zslpx(ji,jj,jk) ) * MIN(    ABS( zslpx(ji  ,jj,jk) ),   &
                  &                                                2.d0*ABS( zwx  (ji-1,jj,jk) ),   &
                  &                                                2.d0*ABS( zwx  (ji  ,jj,jk) ) )
               zslpy(ji,jj,jk) = SIGN( 1.d0, zslpy(ji,jj,jk) ) * MIN(    ABS( zslpy(ji,jj  ,jk) ),   &
                  &                                                2.d0*ABS( zwy  (ji,jj-1,jk) ),   &
                  &                                                2.d0*ABS( zwy  (ji,jj  ,jk) ) )
            END DO
         END DO
      END DO

    end subroutine tra_adv_compute_04_fort

    subroutine tra_adv_compute_05_fort(zwx, zwy, mydomain, zind, zslpx, zslpy, pun, pvn, jpi, jpj, jpk, iter)   
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout) :: zwx, zwy
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: mydomain, zind, zslpx, zslpy
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: pun, pvn

      REAL*8 :: z0u, zalpha, zu, zzwx, zzwy, z0v, zv

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter            

      INTEGER :: ji, jj, jk

      DO jk = 1, jpk-1
         DO jj = 2, jpj-1
            DO ji = 2, jpi-1
               z0u = SIGN( 0.5d0, pun(ji,jj,jk) )
               zalpha = 0.5d0 - z0u
               zu  = z0u - 0.5d0 * pun(ji,jj,jk) * 1.0
               
               zzwx = mydomain(ji+1,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji+1,jj,jk))
               zzwy = mydomain(ji  ,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji  ,jj,jk))
               
               zwx(ji,jj,jk) = pun(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
               
               z0v = SIGN( 0.5d0, pvn(ji,jj,jk) )
               zalpha = 0.5d0 - z0v
               zv  = z0v - 0.5d0 * pvn(ji,jj,jk) * 1.0
               
               zzwx = mydomain(ji,jj+1,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj+1,jk))
               zzwy = mydomain(ji,jj  ,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj  ,jk))
               
               zwy(ji,jj,jk) = pvn(ji,jj,jk) * ( zalpha * zzwx + (1.d0-zalpha) * zzwy )
            END DO
         END DO
      END DO

    end subroutine tra_adv_compute_05_fort

    subroutine tra_adv_compute_06_fort(mydomain, zwx, zwy, jpi, jpj, jpk, iter)   

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout) :: mydomain
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: zwx, zwy

      REAL*8 :: zbtr, ztra

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter      

      INTEGER :: ji, jj, jk
         
      DO jk = 1, jpk-1
         DO jj = 2, jpj-1     
            DO ji = 2, jpi-1
               zbtr = 1.
               ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                  &               + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
               mydomain(ji,jj,jk) = mydomain(ji,jj,jk) + ztra
            END DO
         END DO
      END DO

   end subroutine tra_adv_compute_06_fort

   subroutine tra_adv_compute_07_fort(zwx, tmask, mydomain, jpi, jpj, jpk, iter)   

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout) :: zwx
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: tmask, mydomain      

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter     

      INTEGER :: ji, jj, jk
    
      zwx (:,:, 1 ) = 0.e0    ;    zwx (:,:,jpk) = 0.e0

      DO jk = 2, jpk-1   
         zwx(:,:,jk) = tmask(:,:,jk) * ( mydomain(:,:,jk-1) - mydomain(:,:,jk) )
      END DO

    end subroutine tra_adv_compute_07_fort

    subroutine tra_adv_compute_08_fort(zslpx, zwx, jpi, jpj, jpk, iter)   

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout) :: zslpx
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: zwx

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter     

      INTEGER :: ji, jj, jk

      zslpx(:,:,1) = 0.e0

      DO jk = 2, jpk-1    
         DO jj = 1, jpj
            DO ji = 1, jpi
               zslpx(ji,jj,jk) =                    ( zwx(ji,jj,jk) + zwx(ji,jj,jk+1) )   &
                  &            * ( 0.25d0 + SIGN( 0.25d0, zwx(ji,jj,jk) * zwx(ji,jj,jk+1) ) )
            END DO
         END DO
      END DO

      DO jk = 2, jpk-1     
         DO jj = 1, jpj
            DO ji = 1, jpi
               zslpx(ji,jj,jk) = SIGN( 1.d0, zslpx(ji,jj,jk) ) * MIN( ABS( zslpx(ji,jj,jk  ) ), &
                  &                                               2.d0*ABS( zwx  (ji,jj,jk+1) ),   &
                  &                                               2.d0*ABS( zwx  (ji,jj,jk  ) )  )
            END DO
         END DO
      END DO

   end subroutine tra_adv_compute_08_fort

   subroutine tra_adv_compute_09_fort(zwx, pwn, mydomain, zslpx, zind, jpi, jpj, jpk, iter)   
    
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout) :: zwx
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: pwn, mydomain, zslpx, zind

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter     

      INTEGER :: ji, jj, jk

      REAL*8 :: z0w, zalpha, zw, zzwx, zzwy

      zwx(:,:, 1 ) = pwn(:,:,1) * mydomain(:,:,1)

      !DO jk = 1, jpk-1
      !   DO jj = 2, jpj-1     
      !      DO ji = 2, jpi-1
      !         z0w = SIGN( 0.5d0, pwn(ji,jj,jk+1) )
      !         zalpha = 0.5d0 + z0w
      !         zw  = z0w - 0.5d0 * pwn(ji,jj,jk+1) * zdt * zbtr
      !         
      !         zzwx = mydomain(ji,jj,jk+1) + zind(ji,jj,jk) * (zw * zslpx(ji,jj,jk+1))
      !         zzwy = mydomain(ji,jj,jk  ) + zind(ji,jj,jk) * (zw * zslpx(ji,jj,jk  ))
      !         
      !         zwx(ji,jj,jk+1) = pwn(ji,jj,jk+1) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
      !      END DO
      !   END DO
      !END DO
      DO jk = 2, jpk
         DO jj = 2, jpj-1     
            DO ji = 2, jpi-1
               z0w = SIGN( 0.5d0, pwn(ji,jj,jk) )
               zalpha = 0.5d0 + z0w
               zw  = z0w - 0.5d0 * pwn(ji,jj,jk) * 1.0 * 1.0
               
               zzwx = mydomain(ji,jj,jk  ) + zind(ji,jj,jk-1) * (zw * zslpx(ji,jj,jk  ))
               zzwy = mydomain(ji,jj,jk-1) + zind(ji,jj,jk-1) * (zw * zslpx(ji,jj,jk-1))
               
               zwx(ji,jj,jk) = pwn(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
            END DO
         END DO
      END DO

    end subroutine tra_adv_compute_09_fort

    subroutine tra_adv_compute_10_fort(mydomain, zwx, jpi, jpj, jpk, iter)   

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout) :: mydomain
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in) :: zwx

      INTEGER, INTENT(IN) :: jpi, jpj, jpk
      INTEGER, INTENT(IN) :: iter     

      REAL*8 :: ztra
      
      INTEGER :: ji, jj, jk

      DO jk = 1, jpk-1
         DO jj = 2, jpj-1     
            DO ji = 2, jpi-1
               ztra = -1.0 * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
               mydomain(ji,jj,jk) = ztra
            END DO
         END DO
      END DO

   end subroutine tra_adv_compute_10_fort

end module tra_adv_compute_mod
