program gocean2d
  use grid_mod
  use field_mod
  use model_mod
  !!! A Horizontal 2D hydrodynamic ocean model which
  !!   1) using structured grid
  !!   2) using direct data addressig structures

  implicit none
                                                                            
                                                                            
  INTEGER, ALLOCATABLE :: pt(:,:)   ! properties of t-cells 
                                    ! 1: water cell within computational domain
                                    ! 0: land cell
                                    !-1: water cell outside computational domain
  
  REAL(wp), ALLOCATABLE :: ht(:,:), hu(:,:), hv(:,:), hf(:,:) 
  
  REAL(wp), ALLOCATABLE :: sshb(:,:), sshb_u(:,:), sshb_v(:,:)
  REAL(wp), ALLOCATABLE :: sshn(:,:), sshn_u(:,:), sshn_v(:,:)
  REAL(wp), ALLOCATABLE :: ssha(:,:), ssha_u(:,:), ssha_v(:,:)
  
  REAL(wp), ALLOCATABLE :: un(:,:),  vn(:,:)
  REAL(wp), ALLOCATABLE :: ua(:,:),  va(:,:)
         
  !> The grid on which our fields are defined
  type(grid_type), target :: model_grid
  type(r2d_field_type) :: blah

  ! time stepping index
  integer  :: istp   

  ! Create the model grid
  model_grid = grid_type(ARAKAWA_C)

  !! read in model parameters and read in or setup model grid 
  CALL model_init(model_grid)

  !! setup model initial condition
  CALL initialisation

  !! time stepping 
  DO istp = nit000, nitend, 1
     print*, 'istp == ', istp
     CALL step
  END DO

  !! finalise the model run
  CALL finalisation
  
  WRITE(*,*) 'Simulation finished!!'
end PROGRAM gocean2d

!+++++++++++++++++++++++++++++++++++

SUBROUTINE grid
  !  need a kernel for boundary /grid 
  !! Allocate working arrays
  !! Define (or read in) model grid

  !jphgr_msh = 0    ! read in this from a namelist file
  !jphgr_msh = 1    ! define manually 
          

  SELECT CASE( jphgr_msh)

  CASE(0) ! read in grid from a coordinate file

     ! to be added
     STOP "Reading the grid definition from file is not yet supported!"
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


     xt(1,1) = 0.0_wp + 0.5_wp * e1t(1,1)
     yt(1,1) = 0.0_wp + 0.5_wp * e2t(1,1)

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
     STOP "Wrong grid definition type, check your setup !!!!"
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
          INTEGER :: ierr(6)
          

          ALLOCATE(ht(jpi,jpj), hu(0:jpi,jpj), hv(jpi,0:jpj), hf(0:jpi,0:jpj), STAT=ierr(1))

          ALLOCATE(sshb(jpi,jpj), sshb_u(0:jpi,jpj), sshb_v(jpi,0:jpj), STAT=ierr(2))
          ALLOCATE(sshn(jpi,jpj), sshn_u(0:jpi,jpj), sshn_v(jpi,0:jpj), STAT=ierr(3))
          ALLOCATE(ssha(jpi,jpj), ssha_u(0:jpi,jpj), ssha_v(jpi,0:jpj), STAT=ierr(4))

          ALLOCATE(un(0:jpi,jpj), vn(jpi,0:jpj), ua(0:jpi,jpj), va(jpi,0:jpj), STAT=ierr(5))

          ALLOCATE(pt(0:jpi+1,0:jpj+1), STAT=ierr(6))

          IF(ANY(ierr /= 0, 1)) STOP "in SUBROUTINE ALLOCATION: failed to allocate arrays"

        END SUBROUTINE allocation


!+++++++++++++++++++++++++++++++++++


        SUBROUTINE initialisation
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

            CALL bc(0._wp)


        END SUBROUTINE initialisation


!+++++++++++++++++++++++++++++++++++


        SUBROUTINE step
          use momentum_mod, only: momentum
          implicit none
          REAL(wp) :: rtime

          rtime = REAL(istp, wp) * rdt

          CALL continuity
          CALL momentum
          CALL bc(rtime)  ! open and solid boundary condition
          CALL next
          IF(MOD(istp, irecord) == 0)  CALL output

        END SUBROUTINE step

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE continuity(ssha, sshn, sshn_u, sshn_v, &
                              hu, hv, un, vn, &
                              e12t, rdt)
          implicit none
          real(wp), dimension(:,:), intent(out) :: ssha
          real(wp), dimension(:,:), intent(in)  :: sshn, sshn_u
          real(wp), dimension(:,:), intent(in)  :: hu, hv
          real(wp), dimension(:,:), intent(in)  :: un, vn
          real(wp), dimension(:,:), intent(in)  :: e12t
          real(wp),                 intent(in)  :: rdt
          ! Locals
          integer :: ji, jj
          real(wp) :: rtmp1, rtmp2, rmp3, rtmp4

!kernel continuity
          DO jj = 1, jpj
            DO ji = 1, jpi
              rtmp1 = (sshn_u(ji  ,jj  ) + hu(ji  ,jj  )) * un(ji  ,jj  )
              rtmp2 = (sshn_u(ji-1,jj  ) + hu(ji-1,jj  )) * un(ji-1,jj  )
              rtmp3 = (sshn_v(ji  ,jj  ) + hv(ji  ,jj  )) * vn(ji  ,jj  )
              rtmp4 = (sshn_v(ji  ,jj-1) + hv(ji  ,jj-1)) * vn(ji  ,jj-1)
              ssha(ji,jj) = sshn(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * rdt / e12t(ji,jj)
            END DO
          END DO
!end kernel continuity
        END SUBROUTINE continuity

!+++++++++++++++++++++++++++++++++++

        SUBROUTINE bc(rtime)
          REAL(wp), INTENT(IN) :: rtime
          REAL(wp) :: amp_tide, omega_tide
          INTEGER :: jiu, jiv

          !open boundary condition of clamped ssh

!kernel ssh clamped obc
            amp_tide   = 0.2_wp
            omega_tide = 2.0_wp * 3.14159_wp / (12.42_wp * 3600._wp)
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
            DO jj = 1, jpj
              DO ji = 0, jpi
                IF(pt(ji,jj) * pt(ji+1,jj) == 0) ua(ji,jj) = 0._wp
              END DO
            END DO
!end kernel "solid boundary conditions for u-velocity" 

!kernel "solid boundary conditions for v-velocity" 
            DO jj = 0, jpj
              DO ji = 1, jpi
                IF(pt(ji,jj) * pt(ji,jj+1) == 0) va(ji,jj) = 0._wp
              END DO
            END DO
!end kernel "solid boundary conditions for v-velocity" 

          !                                            Du                 Dssh
          !start of "Flather open boundary condition [---- = sqrt(g/H) * ------]" Kernel
          !                                            Dn                 Dn
          ! ua and va in du/dn should be the specified tidal forcing


! kernel Flather u 
            DO jj = 1, jpj
              DO ji = 0, jpi  
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
            DO jj = 0, jpj
              DO ji = 1, jpi
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
!end kernel flather v .

        END SUBROUTINE bc

!+++++++++++++++++++++++++++++++++++


        SUBROUTINE next
          ! update the now-velocity and ssh

  
! kernel  un updating
          DO jj = 1, jpj
            DO ji = 0, jpi
              un(ji,jj)   = ua(ji,jj)
            END DO
          END DO
! end of kernel sshn_u updating.

! kernel vn updating
          DO jj = 0, jpj
            DO ji = 1, jpi
              vn(ji,jj)   = va(ji,jj)
            END DO
          END DO
! end kernel vn updating.

! kernel sshn updating
          DO jj = 1, jpj
            DO ji = 1, jpi
              sshn(ji,jj) = ssha(ji,jj)
            END DO
          END DO
! end kernel sshn_u updating.

! kernel sshn_u updating
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
! end kernel sshn_v updating.
            
        END SUBROUTINE next

!+++++++++++++++++++++++++++++++++++


        SUBROUTINE output

          ! output model results
          CHARACTER(len=5) :: fname
          WRITE(fname, '(I5.5)') istp
          !OPEN(1, file='go2d_'//fname//'.dat', STATUS='UNKNOWN')
          OPEN(1, file='go2d_'//fname//'.csv', STATUS='UNKNOWN')
          REWIND(1)

          DO jj = 1, jpj
            DO ji = 1, jpi
              rtmp1 = 0.5_wp * (un(ji-1,jj) + un(ji,jj))
              rtmp2 = 0.5_wp * (vn(ji,jj-1) + vn(ji,jj))

              ! write "x-coord, y-coord, depth, ssh, u-velocity, v-velocity" to ASCII files

              !WRITE(1,'(2f20.3, 2f15.4, 2e18.3)')  &            
              WRITE(1,'(f20.3,'','',f20.3,'','',f15.4,'','',f15.4,'','',f18.3,'','',f18.3)') &
                   & xt(ji,jj), yt(ji,jj), ht(ji,jj), sshn(ji,jj),rtmp1, rtmp2 
            END DO
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

        END SUBROUTINE finalisation

!+++++++++++++++++++++++++++++++++++

    END PROGRAM gocean2D

