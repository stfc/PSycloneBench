module boundary_conditions_mod
  implicit none

contains
  
  !================================================

  subroutine boundary_conditions_init(grid)
    use kind_params_mod
    use grid_mod
    use model_mod, only: jpi, jpj, jphgr_msh, pt, ht, hu, hv, &
                         dep_const
    implicit none
    type(grid_type), intent(in) :: grid
    ! Locals
    integer :: ji, jj
    !  need a kernel for boundary /grid 
    !! Allocate working arrays
    !! Define (or read in) model grid

    !jphgr_msh = 0    ! read in this from a namelist file
    !jphgr_msh = 1    ! define manually 
          
    SELECT CASE( jphgr_msh)

    CASE(0) ! read in grid from a coordinate file

     ! to be added
     STOP "Reading the grid definition from file is not yet supported!"
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

end subroutine boundary_conditions_init

  !===========================================================

  subroutine bc(rtime)
    use kind_params_mod
    use physical_params_mod
    use model_mod, only: jpi, jpj, ua, va, hu, hv, pt
    use model_mod, only: ssha, sshn_u, sshn_v
    implicit none
    real(wp), intent(in) :: rtime
    real(wp) :: amp_tide, omega_tide
    integer :: jiu, jiv
    integer :: ji, jj

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

end module boundary_conditions_mod
