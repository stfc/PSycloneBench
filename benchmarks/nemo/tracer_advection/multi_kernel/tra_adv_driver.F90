program tracer_advection
  USE dl_timer, only: timer_init, timer_register, timer_start, timer_stop, timer_report
  use tra_adv_compute_mod, only: tra_adv_compute_01_fort, tra_adv_compute_02_fort, tra_adv_compute_03_fort, & 
                                 tra_adv_compute_04_fort, tra_adv_compute_05_fort, tra_adv_compute_06_fort, &
                                 tra_adv_compute_07_fort, tra_adv_compute_08_fort, tra_adv_compute_09_fort, &
                                 tra_adv_compute_10_fort
  implicit none
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tsn 
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: pun, pvn, pwn
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mydomain, umask, vmask, tmask, zind
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: zslpx, zslpy, zwx, zwy
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:,:)     :: ztfreez, rnfmsk, upsmsk
  REAL*8, ALLOCATABLE, SAVE, DIMENSION(:)       :: rnfmsk_z
  REAL*8                                        :: r, checksum
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
            tsn(jpi,jpj,jpk), &
            zslpx(jpi,jpj,jpk), &
            zslpy(jpi,jpj,jpk), &
            zwx(jpi,jpj,jpk), &
            zwy(jpi,jpj,jpk))
  
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
      call tra_adv_compute_01_fort(zind, tsn, ztfreez, rnfmsk, rnfmsk_z, upsmsk, tmask, jpi, jpj, jpk, jt)
      call tra_adv_compute_02_fort(zwx, zwy, umask, vmask, mydomain, jpi, jpj, jpk, jt)
      call tra_adv_compute_03_fort(zslpx, zslpy, zwx, zwy, jpi, jpj, jpk, jt)
      call tra_adv_compute_04_fort(zslpx, zslpy, zwx, zwy, jpi, jpj, jpk, jt) 
      call tra_adv_compute_05_fort(zwx, zwy, mydomain, zind, zslpx, zslpy, pun, pvn, jpi, jpj, jpk, jt)
      call tra_adv_compute_06_fort(mydomain, zwx, zwy, jpi, jpj, jpk, jt)
      call tra_adv_compute_07_fort(zwx, tmask, mydomain, jpi, jpj, jpk, jt)
      call tra_adv_compute_08_fort(zslpx, zwx, jpi, jpj, jpk, jt)   
      call tra_adv_compute_09_fort(zwx, pwn, mydomain, zslpx, zind, jpi, jpj, jpk, jt)  
      call tra_adv_compute_10_fort(mydomain, zwx, jpi, jpj, jpk, jt)
  end do

  call timer_stop(step_timer)

  ! Output final field and compute checksum

  open(unit = 24, file = 'output.dat', form='formatted')

  checksum = 0.0d0
  do jk = 1, jpk-1
     do jj = 2, jpj-1
        do ji = 2, jpi-1
           checksum = checksum + mydomain(ji,jj,jk)
           write(24,*) mydomain(ji,jj,jk)
        end do
     end do
  end do

  write(*, "('Checksum for domain ', 2(I4, ' x'), I4, ' (',I4,' iterations) = ',E23.16)") &
       jpi, jpj, jpk, itn_count, checksum

  close(24)

  deallocate( mydomain )
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
  deallocate( zslpx)
  deallocate( zslpy)
  deallocate( zwx)
  deallocate( zwy)

  call timer_report()

end program tracer_advection
