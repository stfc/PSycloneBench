program tracer_advection
  USE dl_timer, only: timer_init, timer_register, timer_start, timer_stop, timer_report
  USE tra_adv_compute_mod, only: tra_adv_compute
  implicit none
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tsn 
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: pun, pvn, pwn
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mydomain, umask, vmask, tmask, zind
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: ztfreez, rnfmsk, upsmsk
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:)       :: rnfmsk_z
  REAL*8                                        :: r
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

  ALLOCATE( mydomain (jpi,jpj,jpk), &
            pun (jpi,jpj,jpk), &
            pvn (jpi,jpj,jpk), &
            pwn (jpi,jpj,jpk), &
            umask (jpi,jpj,jpk), &
            vmask (jpi,jpj,jpk), &
            tmask (jpi,jpj,jpk), &
            zind (jpi,jpj,jpk), &
            ztfreez (jpi,jpj), &
            rnfmsk (jpi,jpj), &
            upsmsk (jpi,jpj), &
            rnfmsk_z (jpk), &
            tsn(jpi,jpj,jpk))

  ! Array initialization

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

  call timer_start(step_timer)

  do jt = 1, itn_count
     call tra_adv_compute(ztfreez, pun, pvn, pwn, umask, vmask, tmask, rnfmsk, rnfmsk_z, upsmsk, mydomain, tsn, jpi, jpj, jpk)
  end do

  call timer_stop(step_timer)

  OPEN(unit = 24, file = 'output.dat', form='formatted')
  
  DO jk = 1, jpk-1
     DO jj = 2, jpj-1
        DO ji = 2, jpi-1
           write(24,*) mydomain(ji,jj,jk)
        END DO
     END DO
  END DO

  CLOSE(24)

  DEALLOCATE( mydomain )
  DEALLOCATE( pun )
  DEALLOCATE( pvn )
  DEALLOCATE( pwn )
  DEALLOCATE( umask)
  DEALLOCATE( vmask)
  DEALLOCATE( tmask)
  DEALLOCATE( zind )
  DEALLOCATE( ztfreez )
  DEALLOCATE( rnfmsk)
  DEALLOCATE( upsmsk)
  DEALLOCATE( rnfmsk_z)
  DEALLOCATE( tsn)

  CALL timer_report()

end program tracer_advection
