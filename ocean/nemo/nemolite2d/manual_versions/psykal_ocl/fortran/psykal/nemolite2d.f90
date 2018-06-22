!
! uses module clfortran.mod
!
program nemolite2d
  use iso_c_binding
  use clfortran
  use opencl_utils_mod
  use kernel_args_mod
  use grid_mod
  use field_mod
  use initialisation_mod, only: initialisation
  use model_mod
  use gocean2d_io_mod, only: model_write
  use gocean_mod,      only: model_write_log
  use psykalite_mod, only: invoke_kernels
  implicit none

  integer irec, i, iallocerr

  character(len=CL_UTIL_STR_LEN) :: version_str, filename

  integer(c_int32_t) :: ierr, arg_idx
  integer(c_size_t) :: size_in_bytes, tmask_size_in_bytes
  integer(c_int32_t), target :: status, istp

  !> The grid on which our fields are defined
  type(grid_type), target :: model_grid
  !> Current ('now') sea-surface height at different grid points
  type(r2d_field), target :: sshn_u_fld, sshn_v_fld, sshn_t_fld
  !> 'After' sea-surface height at different grid points
  type(r2d_field), target :: ssha_u_fld, ssha_v_fld, ssha_t_fld
  !> Distance from sea-bed to mean sea level at the different grid points.
  !! This is not time varying.
  type(r2d_field), target :: ht_fld, hu_fld, hv_fld
  !> Current ('now') velocity components
  type(r2d_field), target :: un_fld, vn_fld
  !> 'After' velocity components
  type(r2d_field), target :: ua_fld, va_fld

  enum, bind(c)
     enumerator :: K_CONTINUITY = 1 ! Start from 1 rather than 0
     enumerator K_MOM_U
     enumerator K_MOM_V
     enumerator K_BC_SSH
     enumerator K_BC_SOLID_U
     enumerator K_BC_SOLID_V
     enumerator K_BC_FLATHER_U
     enumerator K_BC_FLATHER_V
     enumerator K_NEXT_SSH_U
     enumerator K_NEXT_SSH_V
     ! Add any new kernels before this line
     enumerator K_NUM_KERNELS_PLUS_ONE
  end enum
  integer, PARAMETER :: K_NUM_KERNELS = K_NUM_KERNELS_PLUS_ONE - 1

  !> The name of each kernel. This must be the name of the OpenCL
  !! routine.
  character(len=20), dimension(K_NUM_KERNELS) :: kernel_names = &
                 [character(len=20) :: "continuity_code", &
                                       "momentum_u_code", &
                                       "momentum_v_code", &
                                       "bc_ssh_code",     &
                                       "bc_solid_u_code", &
                                       "bc_solid_v_code", &
                                       "bc_flather_u_code", &
                                       "bc_flather_v_code", &
                                       "next_sshu_code", &
                                       "next_sshv_code"]    

  ! C_LOC determines the C address of an object
  ! The variable type must be either a pointer or a target
  integer, target :: num_buffers
  integer(c_int32_t), target :: build_log
  integer(c_intptr_t), target :: ssh_events(2)
  integer(c_intptr_t), target :: device
  integer(c_intptr_t), target :: context
  !> Array of kernel objects used in the program
  integer(c_intptr_t), target :: kernels(K_NUM_KERNELS)
  integer(c_intptr_t), target :: tmask_device
  ! Scratch space for logging messages
  character(len=160) :: log_str

  ! Create the model grid. We use a NE offset (i.e. the U, V and F
  ! points immediately to the North and East of a T point all have the
  ! same i,j index).  This is the same offset scheme as used by NEMO.
  model_grid = grid_type(ARAKAWA_C, &
  !  BC_PERIODIC, BC_NON_PERIODIC ??
                         (/BC_EXTERNAL,BC_EXTERNAL,BC_NONE/), &
                         OFFSET_NE)

  !! read in model parameters and configure the model grid 
  CALL model_init(model_grid)

  ! Create fields on this grid

  ! Sea-surface height now (current time step)
  sshn_u_fld = r2d_field(model_grid, U_POINTS)
  sshn_v_fld = r2d_field(model_grid, V_POINTS)
  sshn_t_fld = r2d_field(model_grid, T_POINTS)

  ! Sea-surface height 'after' (next time step)
  ssha_u_fld = r2d_field(model_grid, U_POINTS)
  ssha_v_fld = r2d_field(model_grid, V_POINTS)
  ssha_t_fld = r2d_field(model_grid, T_POINTS)

  ! Distance from sea-bed to mean sea level
  hu_fld = r2d_field(model_grid, U_POINTS)
  hv_fld = r2d_field(model_grid, V_POINTS)
  ht_fld = r2d_field(model_grid, T_POINTS)

  ! Velocity components now (current time step)
  un_fld = r2d_field(model_grid, U_POINTS)
  vn_fld = r2d_field(model_grid, V_POINTS)

  ! Velocity components 'after' (next time step)
  ua_fld = r2d_field(model_grid, U_POINTS)
  va_fld = r2d_field(model_grid, V_POINTS)

  !! setup model initial conditions
  call initialisation(ht_fld, hu_fld, hv_fld, &
                      sshn_u_fld, sshn_v_fld, sshn_t_fld, &
                      un_fld, vn_fld)

  call model_write(model_grid, 0, ht_fld, sshn_t_fld, un_fld, vn_fld)

  write(log_str, "('Simulation domain = (',I4,':',I4,',',I4,':',I4,')')") &
                       model_grid%simulation_domain%xstart, &
                       model_grid%simulation_domain%xstop,  &
                       model_grid%simulation_domain%ystart, &
                       model_grid%simulation_domain%ystop
  call model_write_log("((A))", TRIM(log_str))

  ! Create the OpenCL kernels needed by this model
  filename = "../kernels/nemolite2d_kernels.aocx"
  call gocean_add_kernels(K_NUM_KERNELS, kernel_names, filename)

  do istp = nit000, nitend
     
     call invoke_kernels()

  end do ! End of time-stepping loop
  
  ! Dump the final fields to disk
  call model_write(model_grid, nitend, ht_fld, sshn_t_fld, un_fld, vn_fld)

  ! Clean-up.
  !> \todo call this from gocean_finalise() once that's on master
  call gocean_release()

 end program nemolite2d
