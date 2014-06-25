program shallow

!> @mainpage BENCHMARK WEATHER PREDICTION PROGRAM FOR COMPARING THE
!! PREFORMANCE OF CURRENT SUPERCOMPUTERS. THE MODEL IS
!! BASED OF THE PAPER - THE DYNAMICS OF FINITE-DIFFERENCE
!! MODELS OF THE SHALLOW-WATER EQUATIONS, BY ROBERT SADOURNY
!! J. ATM. SCIENCES, VOL 32, NO 4, APRIL 1975.
!!     
!! CODE BY PAUL N. SWARZTRAUBER, NATIONAL CENTER FOR
!! ATMOSPHERIC RESEARCH, BOULDER, CO,  OCTOBER 1984.
!! Modified by Juliana Rew, NCAR, January 2006
!
!     In this version, shallow4.f, initial and calculated values
!     of U, V, and P are written to a netCDF file
!     for later use in visualizing the results. The netCDF data
!     management library is freely available from
!     http://www.unidata.ucar.edu/software/netcdf
!     This code is still serial but has been brought up to modern
!     Fortran constructs and uses portable intrinsic Fortran 90 timing routines
!     This can be compiled on the IBM SP using:
!     xlf90 -qmaxmem=-1 -g -o shallow4 -qfixed=132 -qsclk=micro \
!     -I/usr/local/include shallow4.f -L/usr/local/lib32/r4i4 -l netcdf
!     where the -L and -I point to local installation of netCDF
!     
!     Changes from shallow4.f (Annette Osprey, January 2010):
!     - Converted to free-form fortran 90.  
!     - Some tidying up of old commented-out code.   
!     - Explicit type declarations.
!     - Variables n, m, itmax and mprint read in from namelist. 
!     - Dynamic array allocation.
!     - Only write to netcdf at mprint timesteps.
!     - Don't write wrap-around points to NetCDF file.
!     - Use 8-byte reals. 
!
!     This version heavily modified as part of the GOcean-2D project
!     with the mantra "all computation must occur in a kernel."
!     Andrew Porter, April 2014

  use shallow_io_mod
  use timing_mod
  use model_mod
  use grid_mod
  use field_mod
  use initial_conditions_mod
  use time_smooth_mod,  ONLY: manual_invoke_time_smooth
  use manual_invoke_apply_bcs_mod, ONLY: manual_invoke_apply_bcs
  use compute_cu_mod, ONLY: manual_invoke_compute_cu
  use compute_cv_mod, ONLY: manual_invoke_compute_cv
  use compute_z_mod,  ONLY: manual_invoke_compute_z
  use compute_h_mod,  ONLY: manual_invoke_compute_h
  use manual_invoke_compute_new_fields_mod, ONLY: manual_invoke_compute_new_fields
  implicit none

  type(grid_type), target :: model_grid
  !> Pressure at {current,previous,next} time step
  type(r2d_field_type) :: p_fld, pold_fld, pnew_fld
  !> Velocity in x direction at {current,previous,next} time step
  type(r2d_field_type) :: u_fld, uold_fld, unew_fld
  !> Velocity in x direction at {current,previous,next} time step
  type(r2d_field_type) :: v_fld, vold_fld, vnew_fld
  !> Mass flux in x and y directions
  type(r2d_field_type) :: cu_fld, cv_fld
  !> Potential vorticity
  type(r2d_field_type) :: z_fld
  !> Surface height
  type(r2d_field_type) :: h_fld
  !> Stream function
  type(r2d_field_type) :: psi_fld
  REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: psi  

  !> Checksum used for each array
  REAL(KIND=8) :: csum

  !> Loop counter for time-stepping loop
  INTEGER :: ncycle
   
  !> Integer tags for timers
  INTEGER :: idxt0, idxt1

  !  ** Initialisations of model parameters (dt etc) ** 
  CALL model_init()
 
  ! Create the model grid
  model_grid = grid_type(ARAKAWA_C)

  ! Create fields on this grid
  p_fld    = r2d_field_type(model_grid, &
                            T_POINTS,   &
                            BC_PERIODIC)
  pold_fld = r2d_field_type(model_grid, &
                            T_POINTS,   &
                            BC_PERIODIC)
  pnew_fld = r2d_field_type(model_grid, &
                            T_POINTS,   &
                            BC_PERIODIC)

  u_fld    = r2d_field_type(model_grid, &
                            U_POINTS,   &
                            BC_PERIODIC)
  uold_fld = r2d_field_type(model_grid, &
                            U_POINTS,   &
                            BC_PERIODIC)
  unew_fld = r2d_field_type(model_grid, &
                            U_POINTS,   &
                            BC_PERIODIC)

  v_fld    = r2d_field_type(model_grid, &
                            V_POINTS,   &
                            BC_PERIODIC)
  vold_fld = r2d_field_type(model_grid, &
                            V_POINTS,   &
                            BC_PERIODIC)
  vnew_fld = r2d_field_type(model_grid, &
                            V_POINTS,   &
                            BC_PERIODIC)

  cu_fld = r2d_field_type(model_grid, &
                          U_POINTS,   &
                          BC_PERIODIC)

  cv_fld = r2d_field_type(model_grid, &
                          V_POINTS,   &
                          BC_PERIODIC)

  z_fld = r2d_field_type(model_grid, &
                         F_POINTS,   &
                         BC_PERIODIC)

  h_fld = r2d_field_type(model_grid, &
                         T_POINTS,   &
                         BC_PERIODIC)

  psi_fld = r2d_field_type(model_grid,   &
                           ALL_POINTS,   &
                           BC_NONE)

  ! NOTE BELOW THAT TWO DELTA T (TDT) IS SET TO DT ON THE FIRST
  ! CYCLE AFTER WHICH IT IS RESET TO DT+DT.
  ! dt and tdt are examples of fields that are actually a 
  ! single parameter.
  CALL copy_field(dt, tdt)

  !     INITIAL VALUES OF THE STREAM FUNCTION AND P

  CALL init_initial_condition_params()
  CALL invoke_init_stream_fn_kernel(psi_fld)
  CALL init_pressure(p_fld)

  !     INITIALIZE VELOCITIES
 
  !> \todo Remove need to pass m and n to init_velocity_u()
  CALL init_velocity_u(u_fld, psi_fld, m, n)
  CALL init_velocity_v(v_fld, psi_fld)

  !     PERIODIC CONTINUATION
  CALL manual_invoke_apply_bcs(u_fld)
  CALL manual_invoke_apply_bcs(v_fld)

  ! Initialise fields that will hold data at previous time step
  CALL copy_field(u_fld, uold_fld)
  CALL copy_field(v_fld, vold_fld)
  CALL copy_field(p_fld, pold_fld)
     
  ! Write intial values of p, u, and v into a netCDF file   
  CALL model_write(0, p_fld, u_fld, v_fld)

  !     Start timer
  CALL timer_start('Time-stepping',idxt0)

  !  ** Start of time loop ** 
  DO ncycle=1,itmax
    
    ! COMPUTE CAPITAL U, CAPITAL V, Z AND H

    CALL timer_start('Compute c{u,v},z,h', idxt1)

    CALL manual_invoke_compute_cu(CU_fld, p_fld, u_fld)
    CALL manual_invoke_compute_cv(CV_fld, p_fld, v_fld)
    CALL manual_invoke_compute_z(z_fld, p_fld, u_fld, v_fld)
    CALL manual_invoke_compute_h(h_fld, p_fld, u_fld, v_fld)

    CALL timer_stop(idxt1)

    ! PERIODIC CONTINUATION

    CALL timer_start('PBCs-1',idxt1)
    ! This call could be generated automatically by PSyclone
    CALL manual_invoke_apply_bcs(CU_fld)
    CALL manual_invoke_apply_bcs(CV_fld)
    CALL manual_invoke_apply_bcs(H_fld)
    CALL manual_invoke_apply_bcs(Z_fld)
    CALL timer_stop(idxt1)

    ! COMPUTE NEW VALUES U,V AND P

    CALL timer_start('Compute new fields', idxt1)
    CALL manual_invoke_compute_new_fields(unew_fld, uold_fld, &
                                          vnew_fld, vold_fld, &
                                          pnew_fld, pold_fld, &
                                          z_fld, cu_fld, cv_fld, &
                                          h_fld, tdt%data)
    CALL timer_stop(idxt1)

    ! PERIODIC CONTINUATION
    CALL timer_start('PBCs-2',idxt1)
    ! This call could be generated by PSyclone
    CALL manual_invoke_apply_bcs(UNEW_fld)
    CALL manual_invoke_apply_bcs(VNEW_fld)
    CALL manual_invoke_apply_bcs(PNEW_fld)
    CALL timer_stop(idxt1)

    ! Time is in seconds but we never actually need it
    !CALL increment(time, dt)

    CALL model_write(ncycle, p_fld, u_fld, v_fld)

    ! TIME SMOOTHING AND UPDATE FOR NEXT CYCLE
    IF(NCYCLE .GT. 1) then

      CALL timer_start('Time smoothing',idxt1)

      CALL manual_invoke_time_smooth(u_fld, UNEW_fld, UOLD_fld)
      CALL manual_invoke_time_smooth(v_fld, VNEW_fld, VOLD_fld)
      CALL manual_invoke_time_smooth(p_fld, PNEW_fld, POLD_fld)

      CALL timer_stop(idxt1)

    ELSE ! ncycle == 1

      ! Make TDT actually = 2*DT
      CALL increment(tdt, tdt)

    ENDIF ! ncycle > 1

    CALL timer_start('Field copy',idxt1)

    CALL copy_field(UNEW_fld, U_fld)
    CALL copy_field(VNEW_fld, V_fld)
    CALL copy_field(PNEW_fld, p_fld)

    CALL timer_stop(idxt1)

  END DO

  !  ** End of time loop ** 

  CALL timer_stop(idxt0)

  CALL compute_checksum(pnew_fld%data, csum)
  CALL model_write_log("('P CHECKSUM after ',I6,' steps = ',E15.7)", &
                       itmax, csum)

  CALL compute_checksum(unew_fld%data, csum)
  CALL model_write_log("('U CHECKSUM after ',I6,' steps = ',E15.7)", &
                       itmax, csum)

  CALL compute_checksum(vnew_fld%data, csum)
  CALL model_write_log("('V CHECKSUM after ',I6,' steps = ',E15.7)", &
                       itmax, csum)

  CALL model_finalise()

CONTAINS

  !===================================================

  SUBROUTINE compute_checksum(field, val)
    IMPLICIT none
    REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: field
    REAL(KIND=8), INTENT(out) :: val

    val = SUM(field)

  END SUBROUTINE compute_checksum

  !===================================================

END PROGRAM shallow
