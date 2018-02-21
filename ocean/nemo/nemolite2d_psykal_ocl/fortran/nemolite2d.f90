!
! uses module clfortran.mod
!
program nemolite2d
  use iso_c_binding
  use clfortran
  use opencl_utils_mod
  use kernel_args_mod, only: set_continuity_args
  implicit none

  integer irec, i, iallocerr
  integer, parameter :: wp = kind(1.0d0)
  ! The number of OpenCL command queues we will use
  integer, parameter :: NUM_CMD_QUEUES = 2

  character(len=CL_UTIL_STR_LEN) :: version_str, filename

  integer(c_int32_t) :: ierr, arg_idx
  integer(c_size_t) :: size_in_bytes
  integer(c_int32_t), target :: status

  real(kind=wp), allocatable, dimension(:,:), target :: un, vn, sshn
  real(kind=wp), allocatable, dimension(:,:), target :: sshn_u, sshn_v, ssha
  real(kind=wp), allocatable, dimension(:,:), target :: hu, hv, e12t

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

  ! C_LOC determines the C address of an object
  ! The variable type must be either a pointer or a target
  integer, target :: nx, ny, num_buffers
  real(kind=wp), target :: rdt
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
  integer(c_intptr_t), target :: ssha_device
  integer(c_intptr_t), target :: sshn_device, sshn_u_device, sshn_v_device
  integer(c_intptr_t), target :: un_device, vn_device
  integer(c_intptr_t), target :: hu_device, hv_device
  integer(c_intptr_t), target :: e12t_device

  nx = 128
  ny = 128
  allocate(un(nx, ny), vn(nx, ny), sshn_u(nx, ny),   &
           sshn_v(nx, ny), sshn(nx,ny), ssha(nx, ny), hu(nx, ny), &
           hv(nx, ny), e12t(nx, ny), stat=iallocerr)
  if(iallocerr /= 0) stop 'Allocate of fields failed'
  un = 1.0d0
  vn = 1.0d0
  sshn_u = 1.0d0
  sshn_v = 1.0d0
  ssha = 1.0d0
  hu = 1.0d0
  hv = 1.0d0
  
  ! Initialise the OpenCL device
  call init_device(device, version_str, context)

  do i=1, NUM_CMD_QUEUES
     cmd_queues(i) = clCreateCommandQueue(context, device, &
                                          CL_QUEUE_PROFILING_ENABLE, ierr)
     call check_status('clCreateCommandQueue', ierr)
  end do

  filename = "../kernels/nemolite2d_kernels.aocx"
  prog = get_program(context, device, version_str, filename)

  kernels(K_CONTINUITY) = get_kernel(prog, 'continuity_code')

  ierr=clReleaseProgram(prog)
  if (ierr.ne.0) stop 'clReleaseProgram'

  size_in_bytes = int(nx*ny,8)*8_8
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

  ! Velocity fields
  un_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  vn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  call check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1

  e12t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
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
 			      size_in_bytes, C_LOC(ssha), 0,            &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_device, CL_TRUE, 0_8, &
			      size_in_bytes, C_LOC(sshn), 0,            &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u_device, CL_TRUE,0_8, &
			      size_in_bytes, C_LOC(sshn_u), 0,           &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v_device, CL_TRUE,0_8, &
			      size_in_bytes, C_LOC(sshn_v), 0,           &
			      C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), hu_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(hu), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), hv_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(hv), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), un_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(un), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), vn_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(vn), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queues(1), e12t_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(e12t), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  call check_status("clEnqueueWriteBuffer", ierr)

  ! Wait for the writes
  write(*,*) "Sent ", num_buffers, " buffers to device"
  ierr = clWaitForEvents(num_buffers, C_LOC(write_events));
  call check_status("clWaitForEvents", ierr)

  ! Set-up the kernel arguments
  call set_continuity_args(kernels(K_CONTINUITY), &
                           nx,                    &
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
  globalsize(1) = nx
  globalsize(2) = ny

  ! execute the kernel
  ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernels(K_CONTINUITY), &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR,C_NULL_PTR)
  call check_status('clEnqueueNDRangeKernel', ierr)

  ierr=clFinish(cmd_queues(1))
  call check_status('clFinish', ierr)

  ! read the resulting vector from device memory
  ierr = clEnqueueReadBuffer(cmd_queues(1), ssha_device, &
                             CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha), &
                             0, C_NULL_PTR, C_LOC(write_events(1)))
  if (ierr.ne.0) stop 'clEnqueueReadBuffer'

  ierr = clWaitForEvents(1, C_LOC(write_events))
  call check_status('clWaitForEvents', ierr)

  ierr=clReleaseKernel(kernels(K_CONTINUITY))
  if (ierr.ne.0) stop 'clReleaseKernel'
  do i=1, NUM_CMD_QUEUES
     ierr=clReleaseCommandQueue(cmd_queues(i))
     if (ierr.ne.0) stop 'clReleaseCommandQueue'
  end do
  ierr=clReleaseContext(context)
  if (ierr.ne.0) stop 'clReleaseContext'

 end program nemolite2d
