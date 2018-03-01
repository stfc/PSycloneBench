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
!  use boundary_conditions_mod
  use gocean2d_io_mod, only: model_write
  use gocean_mod,      only: model_write_log
  implicit none

  integer irec, i, iallocerr
  ! The number of OpenCL command queues we will use
  integer, parameter :: NUM_CMD_QUEUES = 2

  character(len=CL_UTIL_STR_LEN) :: version_str, filename

  integer(c_int32_t) :: ierr, arg_idx
  integer(c_size_t) :: size_in_bytes, tmask_size_in_bytes
  integer(c_int32_t), target :: status

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
  integer(c_size_t),target :: globalsize(2), localsize(2)
  integer(c_int32_t), target :: build_log
  integer(c_intptr_t), allocatable, target :: write_events(:)
  integer(c_intptr_t), target :: device
  integer(c_intptr_t), target :: context, prog
  !> Array of kernel objects used in the program
  integer(c_intptr_t), target :: kernels(K_NUM_KERNELS)
  !> Array of command queues - used to achieve concurrent execution
  integer(c_intptr_t), target :: cmd_queues(NUM_CMD_QUEUES)
  ! Pointers to device memory
  integer(c_intptr_t), target :: ssha_device, ssha_u_device, ssha_v_device
  integer(c_intptr_t), target :: sshn_device, sshn_u_device, sshn_v_device
  integer(c_intptr_t), target :: un_device, vn_device
  integer(c_intptr_t), target :: ua_device, va_device
  integer(c_intptr_t), target :: hu_device, hv_device, ht_device
  integer(c_intptr_t), target :: e1u_device, e1v_device, e1t_device
  integer(c_intptr_t), target :: e2u_device, e2v_device, e2t_device
  integer(c_intptr_t), target :: e12u_device, e12v_device, e12t_device
  integer(c_intptr_t), target :: gphiu_device, gphiv_device
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

  
  ! Initialise the OpenCL device
  call init_device(device, version_str, context)

  do i=1, NUM_CMD_QUEUES
     cmd_queues(i) = clCreateCommandQueue(context, device, &
                                          CL_QUEUE_PROFILING_ENABLE, ierr)
     call check_status('clCreateCommandQueue', ierr)
  end do

  filename = "../kernels/nemolite2d_kernels.aocx"
  prog = get_program(context, device, version_str, filename)

  do i=1, K_NUM_KERNELS
     kernels(i) = get_kernel(prog, kernel_names(i))
  end do

  ! Release the program now that we've created the kernels
  ierr = clReleaseProgram(prog)
  call check_status('clReleaseProgram', ierr)

  size_in_bytes = int(model_grid%nx*model_grid%ny, 8)*8_8
  ! tmask is integer
  tmask_size_in_bytes = int(model_grid%nx*model_grid%ny, 8)*4_8
  ! Size of an element, typically: 
  ! 4_8 for integer or real 
  ! 8_8 for double precision or complex 
  ! 16_8 for double complex
  ! second 8 indicates kind=8 (same is true for 8 in int() argument list)

  ! allocate device memory
  num_buffers = 0
  ssha_device = clCreateBuffer(context, CL_MEM_READ_WRITE, &
                               size_in_bytes, C_NULL_PTR, ierr)
  call check_status('clCreateBuffer', ierr)
  num_buffers = num_buffers + 1
  ssha_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, &
                               size_in_bytes, C_NULL_PTR, ierr)
  call check_status('clCreateBuffer', ierr)
  num_buffers = num_buffers + 1
  ssha_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, &
                               size_in_bytes, C_NULL_PTR, ierr)
  call check_status('clCreateBuffer', ierr)
  num_buffers = num_buffers + 1

  sshn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			       C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  sshn_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
				 C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  sshn_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
				 C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1

  hu_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  hv_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  ht_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1

  ! Velocity fields
  ua_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  va_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  un_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  vn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  gphiu_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
	 		        C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  gphiv_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
	 		        C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1

  ! Mesh scale factors
  e1u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  e1v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  e1t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  e2u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  e2v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  e2t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  e12u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  e12v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  e12t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1

  tmask_device = clCreateBuffer(context, CL_MEM_READ_ONLY, &
                                tmask_size_in_bytes,       &
                                C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1

  ! copy data to device memory
  ! 0_8 is a zero of kind=8
  ! Create an array to store the event associated with each write
  !   to the device
  allocate(write_events(num_buffers), stat=iallocerr)
  if(iallocerr /= 0) stop 'Allocate write_events'
  write (*,*) "Allocated space for ", num_buffers, " events"

  num_buffers = 1;
  ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_device, CL_TRUE, 0_8, &
 			      size_in_bytes, C_LOC(ssha_t_fld%data), 0,   &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_u_device, CL_TRUE, 0_8, &
			      size_in_bytes, C_LOC(ssha_u_fld%data), 0,   &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_v_device, CL_TRUE, 0_8, &
			      size_in_bytes, C_LOC(ssha_v_fld%data), 0,   &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_device, CL_TRUE, 0_8, &
			      size_in_bytes, C_LOC(sshn_t_fld%data), 0,   &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u_device, CL_TRUE,0_8, &
			      size_in_bytes, C_LOC(sshn_u_fld%data), 0,  &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v_device, CL_TRUE,0_8, &
			      size_in_bytes, C_LOC(sshn_v_fld%data), 0,  &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), hu_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(hu_fld%data), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), hv_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(hv_fld%data), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), ht_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(ht_fld%data), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), ua_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(ua_fld%data), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), va_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(va_fld%data), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), un_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(un_fld%data), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), vn_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(vn_fld%data), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)

  ! Mesh properties
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), tmask_device, CL_TRUE, 0_8,      &
			      tmask_size_in_bytes, C_LOC(model_grid%tmask), 0,&
                              C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e1u_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%dx_u), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e1v_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%dx_v), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e1t_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%dx_t), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e2u_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%dy_u), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e2v_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%dy_v), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e2t_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%dy_t), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e12u_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%area_u), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e12v_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%area_v), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e12t_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%area_t), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), gphiu_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%gphiu), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), gphiv_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(model_grid%gphiv), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)

  ! Wait for the writes
  write(*,*) "Sent ", num_buffers, " buffers to device"
  ierr = clWaitForEvents(num_buffers, C_LOC(write_events));
  call check_status("clWaitForEvents", ierr)

  ! Set-up the kernel arguments
  call set_continuity_args(kernels(K_CONTINUITY), &
                           model_grid%nx,         &
                           ssha_device,           &
                           sshn_device,           &
                           sshn_u_device,         &
                           sshn_v_device,         &
                           hu_device,             &
                           hv_device,             &
                           un_device,             &
                           vn_device,             &
                           rdt,                   &
                           e12t_device)

  call set_momu_args(kernels(K_MOM_U), &
                     model_grid%nx,    &
		     ua_device, &
		     un_device, &
		     vn_device, &
		     hu_device, &
		     hv_device, &
		     ht_device, &
		     ssha_u_device, &
		     sshn_device,   &
		     sshn_u_device, &
		     sshn_v_device, &
		     tmask_device, &
		     e1u_device, e1v_device, &
		     e1t_device, e2u_device, &
		     e2t_device, e12u_device,&
		     gphiu_device,           &
		     rdt, cbfr, visc)

  call set_momv_args(kernels(K_MOM_V), &
                     model_grid%nx,    &
		     va_device, &
		     un_device, &
		     vn_device, &
		     hu_device, &
		     hv_device, &
		     ht_device, &
		     ssha_v_device, &
		     sshn_device,   &
		     sshn_u_device, &
		     sshn_v_device, &
		     tmask_device,  &
		     e1v_device, e1t_device, e2u_device, &
		     e2v_device, e2t_device, e12v_device,&
		     gphiv_device, &
		     rdt, cbfr, visc)

  call set_bc_ssh_args(kernels(K_BC_SSH), &
                       model_grid%nx,     &
                       ssha_device,       &
                       tmask_device,      &
                       rdt)
  call set_bc_solid_u_args(kernels(K_BC_SOLID_U), &
                           model_grid%nx,         &
                           ua_device, tmask_device)
  call set_bc_solid_v_args(kernels(K_BC_SOLID_V), &
                           model_grid%nx,         &
                           va_device, tmask_device)
  call set_bc_flather_u_args(kernels(K_BC_FLATHER_U), &
                             model_grid%nx,           &
                             ua_device, hu_device, sshn_u_device, tmask_device)
  call set_bc_flather_v_args(kernels(K_BC_FLATHER_V), &
                             model_grid%nx,           &
                             va_device, hv_device, sshn_v_device, tmask_device)
  call set_next_sshu_args(kernels(K_NEXT_SSH_U), &
                          model_grid%nx,         &
                          sshn_u_device, sshn_device, &
                          tmask_device, e12t_device, e12u_device)
  call set_next_sshv_args(kernels(K_NEXT_SSH_V), &
                          model_grid%nx,         &
                          sshn_v_device, sshn_device, &
                          tmask_device, e12t_device, e12v_device)

  globalsize(1) = model_grid%nx
  globalsize(2) = model_grid%ny

  ! Execute the Continuity kernel
  ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(K_CONTINUITY), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR)
  call check_status('clEnqueueNDRangeKernel', ierr)

  ! Execute the u-Momentum kernel
  ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernels(K_MOM_U), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR)
  call check_status('clEnqueueNDRangeKernel', ierr)

  ! Execute the v-Momentum kernel
  ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernels(K_MOM_V), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR, C_NULL_PTR)
  call check_status('clEnqueueNDRangeKernel', ierr)
  
  ierr = clEnqueueNDRangeKernel(cmd_queues(2), kernels(K_BC_SSH), &
       2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
       0, C_NULL_PTR, C_NULL_PTR)
  call check_status("clEnqueueNDRangeKernel(bc-ssh)", ierr)

  ! Apply boundary conditions
  ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(K_BC_SOLID_U), &
       2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
       0, C_NULL_PTR, C_NULL_PTR)
  call check_status("clEnqueueNDRangeKernel(bc-solid-u)", ierr)
  
  ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(K_BC_SOLID_V), &
       2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
       0, C_NULL_PTR, C_NULL_PTR)
  call check_status("clEnqueueNDRangeKernel(bc-solid-v)", ierr)
  
  ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(K_BC_FLATHER_U), &
       2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
       0, C_NULL_PTR, C_NULL_PTR);
  call check_status("clEnqueueNDRangeKernel(bc-flather-u)", ierr)
  
  ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(K_BC_FLATHER_V), &
       2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
       0, C_NULL_PTR, C_NULL_PTR);
  call check_status("clEnqueueNDRangeKernel(bc-flather-v)", ierr)

  ! Copy 'after' fields to 'now' fields !
  ierr = clEnqueueCopyBuffer(cmd_queues(1), ua_device, un_device, 0_8, 0_8, &
       size_in_bytes,0,C_NULL_PTR,C_NULL_PTR)
  call check_status("clEnqueueCopyBuffer", ierr)
  ierr = clEnqueueCopyBuffer(cmd_queues(1), va_device, vn_device, 0_8, 0_8, &
       size_in_bytes,0,C_NULL_PTR,C_NULL_PTR)
  call check_status("clEnqueueCopyBuffer", ierr)
  ierr = clEnqueueCopyBuffer(cmd_queues(1), ssha_device, sshn_device, 0_8, 0_8, &
       size_in_bytes,0,C_NULL_PTR,C_NULL_PTR)
  call check_status("clEnqueueCopyBuffer", ierr)

  ! Update of sshu and sshv fields
  ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(K_NEXT_SSH_U), &
       2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
       0, C_NULL_PTR, C_NULL_PTR)
  call check_status("clEnqueueNDRangeKernel(next-sshu)", ierr)
                                 
  ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(K_NEXT_SSH_V), &
       2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
       0, C_NULL_PTR, C_NULL_PTR)
  call check_status("clEnqueueNDRangeKernel(next-sshv)", ierr)

  do i=1, NUM_CMD_QUEUES
     ierr=clFinish(cmd_queues(i))
     call check_status('clFinish', ierr)
  end do
  
  ! read the resulting vector from device memory
  ierr = clEnqueueReadBuffer(cmd_queues(1), ssha_device, &
                             CL_TRUE, 0_8, size_in_bytes, &
                             C_LOC(ssha_t_fld%data), &
                             0, C_NULL_PTR, C_LOC(write_events(1)))
  if (ierr.ne.0) stop 'clEnqueueReadBuffer'

  ierr = clWaitForEvents(1, C_LOC(write_events))
  call check_status('clWaitForEvents', ierr)

  do i=1, K_NUM_KERNELS
     ierr = clReleaseKernel(kernels(i))
     call check_status('clReleaseKernel', ierr)
  end do
  do i=1, NUM_CMD_QUEUES
     ierr=clReleaseCommandQueue(cmd_queues(i))
     if (ierr.ne.0) stop 'clReleaseCommandQueue'
  end do
  ierr=clReleaseContext(context)
  if (ierr.ne.0) stop 'clReleaseContext'

 end program nemolite2d
