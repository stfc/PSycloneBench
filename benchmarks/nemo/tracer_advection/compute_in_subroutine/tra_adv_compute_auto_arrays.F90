!> Module containing a subroutine that implements the computation performed by the
!  tracer-advection benchmark. This is done to ease the conversion to SIR (which
!  doesn't support the various bits of general-purpose code needed to set-up and
!  shutdown the benchmark program).
module tra_adv_compute_mod

  implicit none

contains

  subroutine tra_adv_compute(ztfreez, pun, pvn, pwn, umask, vmask,           &
                             tmask, rnfmsk, rnfmsk_z, upsmsk, mydomain, tsn, &
                             jpi, jpj, jpk)
    implicit none

    REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(in)   :: pun, pvn, pwn, umask, vmask, &
                                                           tmask, tsn
    REAL*8, ALLOCATABLE, DIMENSION(:,:), intent(in) :: ztfreez, rnfmsk, upsmsk
    REAL*8, ALLOCATABLE, DIMENSION(:), intent(in) :: rnfmsk_z
    REAL*8, ALLOCATABLE, DIMENSION(:,:,:), intent(inout):: mydomain
    INTEGER, INTENT(IN) :: jpi, jpj, jpk
    
    ! local variables
    REAL*8, DIMENSION(jpi,jpj,jpk) :: zslpx, zslpy, zwx, zwy, zind
    REAL*8                         :: zu, z0u, zzwx, zv, z0v, zzwy, ztra, zbtr, zdt, &
         zalpha, zice, zw, z0w
    INTEGER                        :: ji, jj, jk

    DO jk = 1, jpk
       DO jj = 1, jpj
          DO ji = 1, jpi
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

    zwx(:,:,jpk) = 0.e0   ;   zwy(:,:,jpk) = 0.e0

    DO jk = 1, jpk-1
       DO jj = 1, jpj-1
          DO ji = 1, jpi-1
             zwx(ji,jj,jk) = umask(ji,jj,jk) * ( mydomain(ji+1,jj,jk) - mydomain(ji,jj,jk) )
             zwy(ji,jj,jk) = vmask(ji,jj,jk) * ( mydomain(ji,jj+1,jk) - mydomain(ji,jj,jk) )
          END DO
       END DO
    END DO

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

    DO jk = 1, jpk-1
       zdt  = 1
       DO jj = 2, jpj-1
          DO ji = 2, jpi-1
             z0u = SIGN( 0.5d0, pun(ji,jj,jk) )
             zalpha = 0.5d0 - z0u
             zu  = z0u - 0.5d0 * pun(ji,jj,jk) * zdt
             
             zzwx = mydomain(ji+1,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji+1,jj,jk))
             zzwy = mydomain(ji  ,jj,jk) + zind(ji,jj,jk) * (zu * zslpx(ji  ,jj,jk))
             
             zwx(ji,jj,jk) = pun(ji,jj,jk) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
             
             z0v = SIGN( 0.5d0, pvn(ji,jj,jk) )
             zalpha = 0.5d0 - z0v
             zv  = z0v - 0.5d0 * pvn(ji,jj,jk) * zdt
             
             zzwx = mydomain(ji,jj+1,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj+1,jk))
             zzwy = mydomain(ji,jj  ,jk) + zind(ji,jj,jk) * (zv * zslpy(ji,jj  ,jk))
             
             zwy(ji,jj,jk) = pvn(ji,jj,jk) * ( zalpha * zzwx + (1.d0-zalpha) * zzwy )
          END DO
       END DO
    END DO
    
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
    
    zwx (:,:, 1 ) = 0.e0    ;    zwx (:,:,jpk) = 0.e0
    
    DO jk = 2, jpk-1   
       zwx(:,:,jk) = tmask(:,:,jk) * ( mydomain(:,:,jk-1) - mydomain(:,:,jk) )
    END DO

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
    
    zwx(:,:, 1 ) = pwn(:,:,1) * mydomain(:,:,1)

    zdt  = 1
    zbtr = 1.
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1     
          DO ji = 2, jpi-1
             z0w = SIGN( 0.5d0, pwn(ji,jj,jk+1) )
             zalpha = 0.5d0 + z0w
             zw  = z0w - 0.5d0 * pwn(ji,jj,jk+1) * zdt * zbtr
             
             zzwx = mydomain(ji,jj,jk+1) + zind(ji,jj,jk) * (zw * zslpx(ji,jj,jk+1))
             zzwy = mydomain(ji,jj,jk  ) + zind(ji,jj,jk) * (zw * zslpx(ji,jj,jk  ))
             
             zwx(ji,jj,jk+1) = pwn(ji,jj,jk+1) * ( zalpha * zzwx + (1.-zalpha) * zzwy )
          END DO
       END DO
    END DO

    zbtr = 1.
    DO jk = 1, jpk-1
       DO jj = 2, jpj-1     
          DO ji = 2, jpi-1
             ztra = -zbtr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
             mydomain(ji,jj,jk) = ztra
          END DO
       END DO
    END DO

  end subroutine tra_adv_compute

end module tra_adv_compute_mod
