! gfortran -fcheck=all -ffixed-form -fbacktrace -L/usr/lib64/nvidia -lOpenCL -o sum sum.f
! srun --gres=gpu ./sum
!
! uses module clfortran.mod
!

program sum
  use clfortran
  use ISO_C_BINDING
  implicit none

  integer irec,i,iallocerr,iplatform,idevice
  integer, parameter :: iunit=10
  integer, parameter :: wp = kind(1.0d0)

  integer(c_int32_t) :: ierr,num_devices,num_platforms, arg_idx
  integer(c_size_t) :: iret,size_in_bytes,zero_size= 0
  integer(c_int32_t), target :: status
  integer(c_size_t), target :: binary_size
  character, dimension(1) :: char
  character(len=1024) :: options,kernel_name

  real(kind=wp), allocatable, dimension(:,:), target :: un, vn, sshn
  real(kind=wp), allocatable, dimension(:,:), target :: sshn_u, sshn_v, ssha
  real(kind=wp), allocatable, dimension(:,:), target :: hu, hv, e12t

! C_LOC determines the C address of an object
! The variable type must be either a pointer or a target
  integer, target :: isize, nx, ny, num_buffers
  real(kind=wp), target :: rdt
  integer(c_size_t),target :: globalsize(2), localsize(2)
  integer(c_int32_t), target :: device_cu, build_log
  integer(c_intptr_t), allocatable, target :: &
       platform_ids(:), device_ids(:), write_events(:)
  integer(c_intptr_t), target :: &
       ctx_props(3),context,cmd_queue,prog,kernel,cl_vec1,cl_vec2
  integer(c_intptr_t), target :: ssha_device
  integer(c_intptr_t), target :: sshn_device, sshn_u_device, sshn_v_device
  integer(c_intptr_t), target :: un_device, vn_device
  integer(c_intptr_t), target :: hu_device, hv_device
  integer(c_intptr_t), target :: e12t_device
  character(len=1,kind=c_char), allocatable, target :: &
       source(:),device_name(:)
  character(len=1,kind=c_char), target :: &
       source2(1:1024),retinfo(1:1024), &
       c_options(1:1024),c_kernel_name(1:1024)
  real, allocatable, target :: vec1(:), vec2(:)
  type(c_ptr), target :: psource

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

  ierr = clGetPlatformIDs(0, C_NULL_PTR, num_platforms)
  ierr = check_status('clGetPlatformIDs', ierr)
  if (num_platforms < 1)then
     write (*,*) "Failed to get any OpenCL platform IDs"
     stop
  end if
  print '(a,i2)','Num Platforms: ',num_platforms

  allocate(platform_ids(num_platforms), stat=iallocerr)
  if (iallocerr.ne.0) stop 'memory allocation error'

  ! whenever "&" appears in C subroutine (address-of) call,
  ! then C_LOC has to be used in Fortran
  ierr = clGetPlatformIDs(num_platforms, C_LOC(platform_ids), &
                          num_platforms)
  if (ierr.ne.CL_SUCCESS) stop 'clGetPlatformIDs'

! Get device IDs only for platform 1
  iplatform=1

  ierr=clGetDeviceIDs(platform_ids(iplatform),CL_DEVICE_TYPE_ALL, &
          0,C_NULL_PTR,num_devices)
  if ((ierr.ne.CL_SUCCESS).or.(num_devices.lt.1))then
     stop 'clGetDeviceIDs'
  end if
  print '(a,i2)','Num Devices: ',num_devices
  allocate(device_ids(num_devices),stat=iallocerr)
  if (iallocerr.ne.0) stop 'memory allocation error'

  ! whenever "&" appears in C subroutine (address-off) call,
  ! then C_LOC has to be used in Fortran
  ierr = clGetDeviceIDs(platform_ids(iplatform),CL_DEVICE_TYPE_ALL, &
                        num_devices,C_LOC(device_ids),num_devices)
  if (ierr.ne.CL_SUCCESS) stop 'clGetDeviceIDs'

! Get device info only for device 1
  idevice=1

  ierr=clGetDeviceInfo(device_ids(idevice), &
          CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(device_cu), &
          C_LOC(device_cu),iret)
  if (ierr.ne.CL_SUCCESS) stop 'clGetDeviceInfo'
  ierr=clGetDeviceInfo(device_ids(idevice), &
       CL_DEVICE_NAME,zero_size,C_NULL_PTR,iret)
  if (ierr.ne.CL_SUCCESS) stop 'clGetDeviceInfo'
  allocate(device_name(iret),stat=iallocerr)
  if (iallocerr.ne.0) stop 'allocate'
  ierr=clGetDeviceInfo(device_ids(idevice), &
       CL_DEVICE_NAME,sizeof(device_name),C_LOC(device_name),iret)
  if (ierr.ne.CL_SUCCESS) stop 'clGetDeviceInfo'
  write (*,'(a,i2,a,i3,a)',advance='no') &
          ' Device (#',idevice,', Compute Units: ',device_cu,') - '
  print *,device_name(1:iret)
  deallocate(device_name)

  print '(a,i2,a)', 'Creating context for: ', num_devices,' devices'
  print '(a,i2)', 'for platform: ',iplatform
  ctx_props(1)=CL_CONTEXT_PLATFORM
  ctx_props(2)=platform_ids(iplatform)
  ctx_props(3)=0
  context=clCreateContext(C_LOC(ctx_props),num_devices, &
       C_LOC(device_ids),C_NULL_FUNPTR,C_NULL_PTR,ierr)
  if (ierr.ne.CL_SUCCESS) stop 'clCreateContext'

  cmd_queue=clCreateCommandQueue(context,device_ids(idevice), &
       CL_QUEUE_PROFILING_ENABLE,ierr)
  if (ierr.ne.CL_SUCCESS) stop 'clCreateCommandQueue'

! read kernel from disk
  open(iunit,file='../nemolite2d_kernels.aocx',access='direct', &
       status='old',action='read',iostat=ierr,recl=1)
  if (ierr.ne.0) stop 'Cannot open file ../nemolite2d_kernels.aocx'
  irec=1
  do
     read(iunit,rec=irec,iostat=ierr) char
     if (ierr.ne.0) exit
     irec=irec+1
  end do

  if (irec.eq.0) stop 'nothing read'
  allocate(source(irec+1),stat=iallocerr)
  if (iallocerr.ne.0) stop 'allocate'
  do i=1,irec
     read(iunit,rec=i,iostat=ierr) source(i:i)
  enddo
  close(iunit)

  print '(a,i7)','size of source code in bytes: ',irec

  psource=C_LOC(source) ! pointer to source code
  binary_size = irec
  prog = clCreateProgramWithBinary(context, 1, C_LOC(device_ids(idevice)), &
                                   C_LOC(binary_size), C_LOC(psource), &
                                   C_NULL_PTR, ierr)

  ierr = check_status('clCreateProgramWithSource', ierr)

  ! check if program has uploaded successfully to CL device
  !ierr=clGetProgramInfo(prog,CL_PROGRAM_SOURCE, &
  !     sizeof(source2),C_LOC(source2),iret)
  !if (ierr.ne.CL_SUCCESS) stop 'clGetProgramInfo'
  !print '(a)','**** code retrieved from device start ****'
  !print '(1024a)',source2(1:min(iret,1024))
  !print '(a)','**** code retrieved from device end ****'

  options = "" !'-cl-opt-disable' ! compiler options
  irec = len(trim(options))
  do i=1, irec
     c_options(i)=options(i:i)
  enddo
  c_options(irec+1) = C_NULL_CHAR
  ierr=clBuildProgram(prog, 0, C_NULL_PTR, C_LOC(c_options), &
                      C_NULL_FUNPTR,C_NULL_PTR)
  if (ierr.ne.CL_SUCCESS) then
     print *,'clBuildProgram',ierr
     ierr=clGetProgramBuildInfo(prog,device_ids(idevice), &
          CL_PROGRAM_BUILD_LOG,sizeof(retinfo),C_LOC(retinfo),iret)
     if (ierr.ne.0) stop 'clGetProgramBuildInfo'
     print '(a)','build log start'
     print '(1024a)',retinfo(1:min(iret,1024))
     print '(a)','build log end'
     stop
  endif

  kernel_name = 'continuity_code'
  irec = len(trim(kernel_name))
  do i=1, irec
     c_kernel_name(i)=kernel_name(i:i)
  enddo
  c_kernel_name(irec+1)=C_NULL_CHAR
  kernel=clCreateKernel(prog, C_LOC(c_kernel_name), ierr)
  if (ierr.ne.0) stop 'clCreateKernel'

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
  ierr = check_status('clCreateBuffer', ierr)
  num_buffers = num_buffers + 1
  sshn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			       C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  sshn_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
				 C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  sshn_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
				 C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  hu_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  hv_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1

  ! Velocity fields
  un_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1
  vn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1

  e12t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, size_in_bytes, &
                               C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  num_buffers = num_buffers + 1

  ! copy data to device memory
  ! 0_8 is a zero of kind=8
  ! Create an array to store the event associated with each write
  !   to the device
  allocate(write_events(num_buffers), stat=iallocerr)
  if(iallocerr /= 0) stop 'Allocate write_events'
  write (*,*) "Allocated space for ", num_buffers, " events"

  num_buffers = 1;
  ierr = clEnqueueWriteBuffer(cmd_queue, ssha_device, CL_TRUE, 0_8, &
			     size_in_bytes, C_LOC(ssha), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  ierr = check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queue, sshn_device, CL_TRUE, 0_8, &
			     size_in_bytes, C_LOC(sshn), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  ierr = check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queue, sshn_u_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(sshn_u), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  ierr = check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queue, sshn_v_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(sshn_v), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  ierr = check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queue, hu_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(hu), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  ierr = check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queue, hv_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(hv), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  ierr = check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queue, un_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(un), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  ierr = check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queue, vn_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(vn), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  ierr = check_status("clEnqueueWriteBuffer", ierr)
  num_buffers = num_buffers + 1
  ierr = clEnqueueWriteBuffer(cmd_queue, e12t_device, CL_TRUE,0_8, &
			     size_in_bytes, C_LOC(e12t), 0, &
			     C_NULL_PTR, C_LOC(write_events(num_buffers)))
  ierr = check_status("clEnqueueWriteBuffer", ierr)

  ! Wait for the writes
  write(*,*) "Sent ", num_buffers, " buffers to device"
  ierr = clWaitForEvents(num_buffers, C_LOC(write_events));
  ierr = check_status("clWaitForEvents", ierr)

  ! set the kernel arguments
  arg_idx = 0
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(nx), C_LOC(nx))
  ierr = check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(ssha_device), &
		       C_LOC(ssha_device))
  ierr = check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(sshn_device), &
		       C_LOC(sshn_device))
  ierr = check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(sshn_u_device), &
		       C_LOC(sshn_u_device))
  ierr = check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(sshn_v_device), &
		       C_LOC(sshn_v_device))
  ierr =  check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(hu_device), &
		       C_LOC(hu_device))
  ierr = check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(hv_device), &
		       C_LOC(hv_device))
  ierr = check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(un_device), &
		       C_LOC(un_device))
  ierr = check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(vn_device), &
		       C_LOC(vn_device))
  ierr = check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(rdt), C_LOC(rdt))
  ierr = check_status("clSetKernelArg", ierr)
  arg_idx = arg_idx + 1
  ierr = clSetKernelArg(kernel, arg_idx, sizeof(e12t_device), &
                        C_LOC(e12t_device))
  ierr = check_status("clSetKernelArg", ierr)

  ! get the local size for the kernel
  !ierr=clGetKernelWorkGroupInfo(kernel,device_ids(idevice), &
  !        CL_KERNEL_WORK_GROUP_SIZE,sizeof(localsize), &
  !        C_LOC(localsize),iret)
  !if (ierr.ne.0) stop 'clGetKernelWorkGroupInfo'
  !globalsize=int(isize,8)
  !if (mod(globalsize,localsize).ne.0) globalsize= &
  !     globalsize+localsize-mod(globalsize,localsize) 

  globalsize(1) = nx
  globalsize(2) = ny

  ! execute the kernel
  ierr = clEnqueueNDRangeKernel(cmd_queue, kernel, &
          2, C_NULL_PTR, C_LOC(globalsize), C_NULL_PTR, &
          0, C_NULL_PTR,C_NULL_PTR)
  if (ierr.ne.0) stop 'clEnqueueNDRangeKernel'
  ierr=clFinish(cmd_queue)
  if (ierr.ne.0) stop 'clFinish'

  ! read the resulting vector from device memory
  ierr = clEnqueueReadBuffer(cmd_queue, ssha_device, &
                             CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha), &
                             0, C_NULL_PTR, C_NULL_PTR)
  if (ierr.ne.0) stop 'clEnqueueReadBuffer'

  ierr=clReleaseKernel(kernel)
  if (ierr.ne.0) stop 'clReleaseKernel'
  ierr=clReleaseCommandQueue(cmd_queue)
  if (ierr.ne.0) stop 'clReleaseCommandQueue'
  ierr=clReleaseContext(context)
  if (ierr.ne.0) stop 'clReleaseContext'

contains

function check_status(text, ierr)
  use clfortran
  implicit none
  integer :: check_status
  character(len=*), intent(in) :: text
  integer, intent(in) :: ierr

  logical, parameter :: verbose = .TRUE.

  if(ierr /= CL_SUCCESS)then
    write(*,'("Hit error: ",(A),": ",(A))') text, OCL_GetErrorString(ierr)
    stop
  end if
  if(verbose)then
    write(*,'("Called ",(A)," OK")') text 
  end if
  check_status = 0
end function check_status
  
function OCL_GetErrorString(error)
  use clfortran
  implicit none
  character(len=64) :: OCL_GetErrorString
  integer, intent(in) :: error
  select case(error)

    case (CL_SUCCESS)
        OCL_GetErrorString = "CL_SUCCESS"
    case (CL_DEVICE_NOT_FOUND)
        OCL_GetErrorString = "CL_DEVICE_NOT_FOUND"
    case (CL_DEVICE_NOT_AVAILABLE)
        OCL_GetErrorString = "CL_DEVICE_NOT_AVAILABLE"
    case (CL_COMPILER_NOT_AVAILABLE)
        OCL_GetErrorString = "CL_COMPILER_NOT_AVAILABLE"
    case (CL_MEM_OBJECT_ALLOCATION_FAILURE)
        OCL_GetErrorString = "CL_MEM_OBJECT_ALLOCATION_FAILURE"
    case (CL_OUT_OF_RESOURCES)
        OCL_GetErrorString = "CL_OUT_OF_RESOURCES"
    case (CL_OUT_OF_HOST_MEMORY)
        OCL_GetErrorString = "CL_OUT_OF_HOST_MEMORY"
    case (CL_PROFILING_INFO_NOT_AVAILABLE)
        OCL_GetErrorString = "CL_PROFILING_INFO_NOT_AVAILABLE"
    case (CL_MEM_COPY_OVERLAP)
        OCL_GetErrorString = "CL_MEM_COPY_OVERLAP"
    case (CL_IMAGE_FORMAT_MISMATCH)
        OCL_GetErrorString = "CL_IMAGE_FORMAT_MISMATCH"
    case (CL_IMAGE_FORMAT_NOT_SUPPORTED)
        OCL_GetErrorString = "CL_IMAGE_FORMAT_NOT_SUPPORTED"
    case (CL_BUILD_PROGRAM_FAILURE)
        OCL_GetErrorString = "CL_BUILD_PROGRAM_FAILURE"
    case (CL_MAP_FAILURE)
        OCL_GetErrorString = "CL_MAP_FAILURE"
    case (CL_INVALID_VALUE)
        OCL_GetErrorString = "CL_INVALID_VALUE"
    case (CL_INVALID_DEVICE_TYPE)
        OCL_GetErrorString = "CL_INVALID_DEVICE_TYPE"
    case (CL_INVALID_PLATFORM)
        OCL_GetErrorString = "CL_INVALID_PLATFORM"
    case (CL_INVALID_DEVICE)
        OCL_GetErrorString = "CL_INVALID_DEVICE"
    case (CL_INVALID_CONTEXT)
        OCL_GetErrorString = "CL_INVALID_CONTEXT"
    case (CL_INVALID_QUEUE_PROPERTIES)
        OCL_GetErrorString = "CL_INVALID_QUEUE_PROPERTIES"
    case (CL_INVALID_COMMAND_QUEUE)
        OCL_GetErrorString = "CL_INVALID_COMMAND_QUEUE"
    case (CL_INVALID_HOST_PTR)
        OCL_GetErrorString = "CL_INVALID_HOST_PTR"
    case (CL_INVALID_MEM_OBJECT)
        OCL_GetErrorString = "CL_INVALID_MEM_OBJECT"
    case (CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
        OCL_GetErrorString = "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"
    case (CL_INVALID_IMAGE_SIZE)
        OCL_GetErrorString = "CL_INVALID_IMAGE_SIZE"
    case (CL_INVALID_SAMPLER)
        OCL_GetErrorString = "CL_INVALID_SAMPLER"
    case (CL_INVALID_BINARY)
        OCL_GetErrorString = "CL_INVALID_BINARY"
    case (CL_INVALID_BUILD_OPTIONS)
        OCL_GetErrorString = "CL_INVALID_BUILD_OPTIONS"
    case (CL_INVALID_PROGRAM)
        OCL_GetErrorString = "CL_INVALID_PROGRAM"
    case (CL_INVALID_PROGRAM_EXECUTABLE)
        OCL_GetErrorString = "CL_INVALID_PROGRAM_EXECUTABLE"
    case (CL_INVALID_KERNEL_NAME)
        OCL_GetErrorString = "CL_INVALID_KERNEL_NAME"
    case (CL_INVALID_KERNEL_DEFINITION)
        OCL_GetErrorString = "CL_INVALID_KERNEL_DEFINITION"
    case (CL_INVALID_KERNEL)
        OCL_GetErrorString = "CL_INVALID_KERNEL"
    case (CL_INVALID_ARG_INDEX)
        OCL_GetErrorString = "CL_INVALID_ARG_INDEX"
    case (CL_INVALID_ARG_VALUE)
        OCL_GetErrorString = "CL_INVALID_ARG_VALUE"
    case (CL_INVALID_ARG_SIZE)
        OCL_GetErrorString = "CL_INVALID_ARG_SIZE"
    case (CL_INVALID_KERNEL_ARGS)
        OCL_GetErrorString = "CL_INVALID_KERNEL_ARGS"
    case (CL_INVALID_WORK_DIMENSION)
        OCL_GetErrorString = "CL_INVALID_WORK_DIMENSION"
    case (CL_INVALID_WORK_GROUP_SIZE)
        OCL_GetErrorString = "CL_INVALID_WORK_GROUP_SIZE"
    case (CL_INVALID_WORK_ITEM_SIZE)
        OCL_GetErrorString = "CL_INVALID_WORK_ITEM_SIZE"
    case (CL_INVALID_GLOBAL_OFFSET)
        OCL_GetErrorString = "CL_INVALID_GLOBAL_OFFSET"
    case (CL_INVALID_EVENT_WAIT_LIST)
        OCL_GetErrorString = "CL_INVALID_EVENT_WAIT_LIST"
    case (CL_INVALID_EVENT)
        OCL_GetErrorString = "CL_INVALID_EVENT"
    case (CL_INVALID_OPERATION)
        OCL_GetErrorString = "CL_INVALID_OPERATION"
    case (CL_INVALID_GL_OBJECT)
        OCL_GetErrorString = "CL_INVALID_GL_OBJECT"
    case (CL_INVALID_BUFFER_SIZE)
        OCL_GetErrorString = "CL_INVALID_BUFFER_SIZE"
    case (CL_INVALID_MIP_LEVEL)
        OCL_GetErrorString = "CL_INVALID_MIP_LEVEL"
    case(CL_INVALID_GLOBAL_WORK_SIZE)
        OCL_GetErrorString = "CL_INVALID_GLOBAL_WORK_SIZE"
    case default
        OCL_GetErrorString = "unknown error code"
     end select
   end function OCL_GetErrorString

end program sum






! Argument lists need to be handled with care. Neither the C nor the Fortran compiler checks for mismatched argument types, or even a mismatch in the number of arguments. Some bizarre run-time errors therefore arise. Keep in mind that Fortran passes arguments by reference whereas C passes arguments by value. When you put a variable name into a function call from Fortran, the corresponding C function receives a pointer to that variable. Similarly, when calling a Fortran subroutine from C, you must explicitly pass addresses rather than values in the argument list.

! When passing arrays, remember that C arrays start with subscript zero. Fortran stores multidimensional arrays in column-major order (first index varies fastest) whereas C stores them in row-major order (last index varies fastest).

! Passing character strings is a special problem. C knows it has come to the end of string when it hits a null character, but Fortran uses the declared length of a string. The C-Fortran interface provides an extra argument for each character string in a C argument list to receive the declared length of a string when called from Fortran. Consider the following Fortran fragment:



