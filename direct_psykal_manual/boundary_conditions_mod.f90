module boundary_conditions_mod
  implicit none

contains
  
  !================================================
  
  subroutine boundary_conditions_init(grid)
    use kind_params_mod
    use grid_mod
    use model_mod, only: jphgr_msh, ht, hu, hv, &
                         dep_const
    implicit none
    type(grid_type), intent(inout) :: grid
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

       ! Define Model solid/open Boundaries via the properties of t-cells

       grid%tmask(:,:) = 1                 ! all inner cells

       ! -define solid/open boundaries
       DO jj = 1, grid%ny
          grid%tmask(      1, jj) = 0              ! west solid boundary
          grid%tmask(grid%nx, jj) = 0              ! east solid boundary
       END Do

       DO ji = 1, grid%nx
          grid%tmask(ji,grid%ny) = 0               ! north solid boundary
       END Do

       DO ji = 1, grid%nx
          grid%tmask(ji,1) = -1                    ! south open boundary
       END Do

       ! The actual part of this domain that is simulated. The outer-most 
       ! rows and columns of T points are not in the domain but define the
       ! hard or open boundaries.
       grid%simulation_domain%xstart = 2
       grid%simulation_domain%xstop  = grid%nx - 1
       grid%simulation_domain%ystart = 2
       grid%simulation_domain%ystop  = grid%ny - 1
       grid%simulation_domain%nx = grid%nx - 2
       grid%simulation_domain%ny = grid%ny - 2

       !Depth 

       ! -depth on grid points

       DO jj = 1, grid%ny
          DO ji = 1, grid%nx
             ht(ji,jj) = dep_const 
          END DO
       END DO

       DO jj = 1, grid%ny
          DO ji = 1, grid%nx
             hu(ji,jj) = dep_const 
          END DO
       END DO

       DO jj = 1, grid%ny
          DO ji = 1, grid%nx
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

  subroutine bc(rtime, grid, sshn_u, sshn_v, ssha, ua, va, hu, hv)
    use kind_params_mod
    use physical_params_mod
    use grid_mod
    implicit none
    real(wp),        intent(in) :: rtime
    type(grid_type), intent(in) :: grid
    real(wp), dimension(:,:), intent(in)    :: sshn_u, sshn_v
    real(wp), dimension(:,:), intent(inout) :: ssha
    real(wp), dimension(:,:), intent(in)    :: hu, hv
    real(wp), dimension(:,:), intent(inout) :: ua, va
    ! Locals
    real(wp) :: amp_tide, omega_tide
    integer :: jiu, jiv
    integer :: ji, jj

    !open boundary condition of clamped ssh

    !kernel ssh clamped obc
    amp_tide   = 0.2_wp
    omega_tide = 2.0_wp * 3.14159_wp / (12.42_wp * 3600._wp)
    DO jj = 1, grid%ny
       DO ji = 1, grid%nx
          IF(grid%tmask(ji,jj) <= 0) CYCLE
          IF     (grid%tmask(ji,jj-1) < 0) THEN
             ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
          ELSE IF(grid%tmask(ji,jj+1) < 0) THEN
             ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
          ELSE IF(grid%tmask(ji+1,jj) < 0) THEN
             ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
          ELSE IF(grid%tmask(ji-1,jj) < 0) THEN
             ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
          END IF
       END DO
    END DO
    !end kernel ssh clamped obc


    ! kernel"solid boundary conditions for u-velocity" 
    DO jj = 1, grid%ny
       DO ji = 1, grid%nx-1
          IF(grid%tmask(ji,jj) * grid%tmask(ji+1,jj) == 0) ua(ji,jj) = 0._wp
       END DO
    END DO
    !end kernel "solid boundary conditions for u-velocity" 

    !kernel "solid boundary conditions for v-velocity" 
    DO jj = 1, grid%ny-1
       DO ji = 1, grid%nx
          IF(grid%tmask(ji,jj) * grid%tmask(ji,jj+1) == 0) va(ji,jj) = 0._wp
       END DO
    END DO
    !end kernel "solid boundary conditions for v-velocity" 

    !                                            Du                 Dssh
    !start of "Flather open boundary condition [---- = sqrt(g/H) * ------]" Kernel
    !                                            Dn                 Dn
    ! ua and va in du/dn should be the specified tidal forcing


    ! kernel Flather u 
    DO jj = 1, grid%ny
       DO ji = 1, grid%nx-1  

          IF(grid%tmask(ji,jj) + grid%tmask(ji+1,jj) <= -1) CYCLE                         ! not in the domain

          IF(grid%tmask(ji,jj) < 0) THEN
             jiu = ji + 1
             ua(ji,jj) = ua(jiu,jj) + SQRT(g/hu(ji,jj)) * (sshn_u(ji,jj) - sshn_u(jiu,jj))
          ELSE IF(grid%tmask(ji+1,jj )< 0) THEN
             jiu = ji - 1 
             ua(ji,jj) = ua(jiu,jj) + SQRT(g/hu(ji,jj)) * (sshn_u(ji,jj) - sshn_u(jiu,jj))
          END IF
       END DO
    END DO
    !end kernel flather u .

    !kernel Flather v 
    DO jj = 1, grid%ny-1
       DO ji = 1, grid%nx

          ! Check whether this point is inside the simulated domain
          IF(grid%tmask(ji,jj) + grid%tmask(ji,jj+1) <= -1) CYCLE

          IF(grid%tmask(ji,jj) < 0) THEN
             jiv = jj + 1
             va(ji,jj) = va(ji,jiv) + SQRT(g/hv(ji,jj)) * (sshn_v(ji,jj) - sshn_v(ji,jiv))
          ELSE IF(grid%tmask(ji,jj+1) < 0) THEN
             jiv = jj - 1 
             va(ji,jj) = va(ji,jiv) + SQRT(g/hv(ji,jj)) * (sshn_v(ji,jj) - sshn_v(ji,jiv))
          END IF
       END DO
    END DO
    !end kernel flather v .

  END SUBROUTINE bc

end module boundary_conditions_mod
