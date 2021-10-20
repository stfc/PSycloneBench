   !!=====================================================================================
   !! ***  traadv kernel extracted from the NEMO software (http://www.nemo-ocean.eu ) ***
   !! ***          governed by the CeCILL licence (http://www.cecill.info)            ***
   !!                                                   
   !! ***                             IS-ENES2 - CMCC/STFC                            ***
   !!=====================================================================================
PROGRAM tra_adv
   USE dl_timer, only: timer_init, timer_register, timer_start, timer_stop, &
       timer_report
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tsn 
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: pun, pvn, pwn
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mydomain, zslpx, zslpy, zwx, zwy, umask, vmask, tmask, zind
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: ztfreez, rnfmsk, upsmsk
   REAL*8, ALLOCATABLE, SAVE, DIMENSION(:)       :: rnfmsk_z
   REAL*8                                        :: zice, zu, z0u, zzwx, zv, z0v, zzwy, ztra, zbtr, zdt, zalpha
   REAL*8                                        :: r, checksum
   REAL*8                                        :: zw, z0w
   INTEGER                                       :: jpi, jpj, jpk, ji, jj, jk, jt
   INTEGER*8                                     :: itn_count
   CHARACTER(len=10)                             :: env
   !> Timer indexes, one for initialisation, one for the 'time-stepping'
   INTEGER :: init_timer, step_timer

   CALL get_environment_variable("JPI", env)
   READ ( env, '(i10)' ) jpi
   CALL get_environment_variable("JPJ", env)
   READ ( env, '(i10)' ) jpj
   CALL get_environment_variable("JPK", env)
   READ ( env, '(i10)' ) jpk
   CALL get_environment_variable("IT", env)
   READ ( env, '(i10)' ) itn_count

   ! Set-up our timers

   CALL timer_init()
   CALL timer_register(init_timer, label='Initialisation')
   CALL timer_register(step_timer, label='Time-stepping', num_repeats=itn_count)

   ! Initialisation

   call timer_start(init_timer)

   ALLOCATE( mydomain (jpi,jpj,jpk))
   ALLOCATE( zwx (jpi,jpj,jpk))
   ALLOCATE( zwy (jpi,jpj,jpk))
   ALLOCATE( zslpx (jpi,jpj,jpk))
   ALLOCATE( zslpy (jpi,jpj,jpk))
   ALLOCATE( pun (jpi,jpj,jpk))
   ALLOCATE( pvn (jpi,jpj,jpk))
   ALLOCATE( pwn (jpi,jpj,jpk))
   ALLOCATE( umask (jpi,jpj,jpk))
   ALLOCATE( vmask (jpi,jpj,jpk))
   ALLOCATE( tmask (jpi,jpj,jpk))
   ALLOCATE( zind (jpi,jpj,jpk))
   ALLOCATE( ztfreez (jpi,jpj))
   ALLOCATE( rnfmsk (jpi,jpj))
   ALLOCATE( upsmsk (jpi,jpj))
   ALLOCATE( rnfmsk_z (jpk))
   ALLOCATE( tsn(jpi,jpj,jpk))


! arrays initialization

   r = jpi*jpj*jpk

   ! the following three lines can be uncommented to randomize arrays initialization
   !call random_seed()
   !call random_number(r)
   !r = r*jpi*jpj*jpk

   DO jk = 1, jpk
      DO jj = 1, jpj
          DO ji = 1, jpi
              umask(ji,jj,jk) = ji*jj*jk/r
              mydomain(ji,jj,jk) =ji*jj*jk/r
              pun(ji,jj,jk) =ji*jj*jk/r
              pvn(ji,jj,jk) =ji*jj*jk/r
              pwn(ji,jj,jk) =ji*jj*jk/r
              vmask(ji,jj,jk)= ji*jj*jk/r
              tsn(ji,jj,jk)= ji*jj*jk/r
              tmask(ji,jj,jk)= ji*jj*jk/r
          END DO
      END DO
   END DO

   r = jpi*jpj
   DO jj=1, jpj
      DO ji=1, jpi
         ztfreez(ji,jj) = ji*jj/r
         upsmsk(ji,jj) = ji*jj/r
         rnfmsk(ji,jj) = ji*jj/r
      END DO
   END DO

   DO jk=1, jpk
      rnfmsk_z(jk)=jk/jpk
   END DO

   call timer_stop(init_timer)

!***********************
!* Start of the symphony
!***********************
   call timer_start(step_timer)

   DO jt = 1, itn_count
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
               ztra = - zbtr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
               mydomain(ji,jj,jk) = ztra
            END DO
         END DO
      END DO
  END DO

  call timer_stop(step_timer)

  OPEN(unit = 24, file = 'output.dat', form='formatted')

  checksum = 0.0d0
  DO jk = 1, jpk-1
     DO jj = 2, jpj-1
        DO ji = 2, jpi-1
           checksum = checksum + mydomain(ji,jj,jk)
           write(24,*) mydomain(ji,jj,jk)
        END DO
     END DO
  END DO

  write(*, "('Checksum for domain ', 2(I4, ' x'), I4, ' (',I4,' iterations) = ',E23.16)") &
       jpi, jpj, jpk, itn_count, checksum

  close(24)

  deallocate( mydomain )
  deallocate( zwx )
  deallocate( zwy )
  deallocate( zslpx )
  deallocate( zslpy )
  deallocate( pun )
  deallocate( pvn )
  deallocate( pwn )
  deallocate( umask)
  deallocate( vmask)
  deallocate( tmask)
  deallocate( zind )
  deallocate( ztfreez )
  deallocate( rnfmsk)
  deallocate( upsmsk)
  deallocate( rnfmsk_z)
  deallocate( tsn)

  call timer_report()

end program tra_adv
