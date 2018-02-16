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

  real(kind=wp), allocatable, dimension(:,:) :: un, vn, sshn_t, sshn_u, sshn_v
  real(kind=wp), allocatable, dimension(:,:) :: ua, va, ssha, ssha_u, ssha_v
  real(kind=wp), allocatable, dimension(:,:) :: hu, hv, ht

! C_LOC determines the C address of an object
! The variable type must be either a pointer or a target
  integer, target :: isize, nx, ny
  real(kind=wp), target :: rdt
  integer(c_size_t),target :: globalsize,localsize
  integer(c_int32_t), target :: device_cu,build_log
  integer(c_intptr_t), allocatable, target :: &
       platform_ids(:),device_ids(:)
  integer(c_intptr_t), target :: &
       ctx_props(3),context,cmd_queue,prog,kernel,cl_vec1,cl_vec2
  integer(c_intptr_t), target :: ssha_device, ssha_u_device, ssha_v_device
  integer(c_intptr_t), target :: sshn_device, sshn_u_device, sshn_v_device
  integer(c_intptr_t), target :: un_device, vn_device
  integer(c_intptr_t), target :: ht_device, hu_device, hv_device
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
  allocate(un(nx, ny), vn(nx, ny), sshn_t(nx, ny), sshn_u(nx, ny), &
           sshn_v(nx, ny), ua(nx, ny), va(nx, ny), ssha(nx, ny), &
           ssha_u(nx, ny), ssha_v(nx, ny), hu(nx, ny), hv(nx, ny), ht(nx, ny))
  un = 1.0d0
  vn = 1.0d0
  sshn_t = 1.0d0
  sshn_u = 1.0d0
  sshn_v = 1.0d0
  ua = 1.0d0
  va = 1.0d0
  ssha = 1.0d0
  ssha_u = 1.0d0
  ssha_v = 1.0d0
  hu = 1.0d0
  hv = 1.0d0
  ht = 1.0d0

  ierr=clGetPlatformIDs(0,C_NULL_PTR,num_platforms)
  if ((ierr.ne.CL_SUCCESS).or.(num_platforms.lt.1))then
     print *,'clGetPlatformIDs',ierr
  end if
  print '(a,i2)','Num Platforms: ',num_platforms
  allocate(platform_ids(num_platforms),stat=iallocerr)
  if (iallocerr.ne.0) stop 'memory allocation error'
! whenever "&" appears in C subroutine (address-of) call,
! then C_LOC has to be used in Fortran
  ierr=clGetPlatformIDs(num_platforms,C_LOC(platform_ids),&
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
  ierr=clGetDeviceIDs(platform_ids(iplatform),CL_DEVICE_TYPE_ALL, &
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

  print '(a,i4)','size of source code in bytes: ',irec
! in C, strings end with c_null_char
  source(irec+1)=C_NULL_CHAR

  psource=C_LOC(source) ! pointer to source code
  binary_size = irec
  prog = clCreateProgramWithBinary(context, 1, C_LOC(device_ids(idevice)), &
                                   C_LOC(binary_size), psource, C_LOC(status), &
                                   ierr)
  !prog = clCreateProgramWithSource(context, 1, C_LOC(psource), &
  !                                 C_NULL_PTR, ierr)
  if (ierr.ne.CL_SUCCESS) stop 'clCreateProgramWithSource'

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
  ssha_device = clCreateBuffer(context, CL_MEM_READ_WRITE, &
                               size_in_bytes, C_NULL_PTR, ierr)
  ierr = check_status('clCreateBuffer', ierr)
  ssha_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
				 C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  ssha_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
				 C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  sshn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			       C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  sshn_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
				 C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  sshn_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
				 C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  hu_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  hv_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  ht_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  
  ! Velocity fields
  un_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)
  vn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, size_in_bytes, &
			     C_NULL_PTR, ierr)
  ierr = check_status("clCreateBuffer", ierr)


! copy data to device memory
  ierr=clEnqueueWriteBuffer(cmd_queue,cl_vec1, &
       CL_TRUE,0_8,size_in_bytes,C_LOC(vec1), &
       0,C_NULL_PTR,C_NULL_PTR)
! 0_8 is a zero of kind=8
  if (ierr.ne.0) stop 'clEnqueueWriteBuffer'
  ierr=clEnqueueWriteBuffer(cmd_queue,cl_vec2, &
          CL_TRUE,0_8,size_in_bytes,C_LOC(vec2), &
          0,C_NULL_PTR,C_NULL_PTR)
  if (ierr.ne.0) stop 'clEnqueueWriteBuffer'


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


  ierr=clSetKernelArg(kernel,0,sizeof(isize),C_LOC(isize))
  if (ierr.ne.0) stop 'clSetKernelArg'
  ierr=clSetKernelArg(kernel,1,sizeof(cl_vec1),C_LOC(cl_vec1))
  if (ierr.ne.0) stop 'clSetKernelArg'
  ierr=clSetKernelArg(kernel,2,sizeof(cl_vec2),C_LOC(cl_vec2))
  if (ierr.ne.0) stop 'clSetKernelArg'

! get the local size for the kernel
  ierr=clGetKernelWorkGroupInfo(kernel,device_ids(idevice), &
          CL_KERNEL_WORK_GROUP_SIZE,sizeof(localsize), &
          C_LOC(localsize),iret)
  if (ierr.ne.0) stop 'clGetKernelWorkGroupInfo'
  globalsize=int(isize,8)
  if (mod(globalsize,localsize).ne.0) globalsize= &
       globalsize+localsize-mod(globalsize,localsize) 


! execute the kernel
  ierr=clEnqueueNDRangeKernel(cmd_queue,kernel, &
          1,C_NULL_PTR,C_LOC(globalsize),C_LOC(localsize), &
          0,C_NULL_PTR,C_NULL_PTR)
  if (ierr.ne.0) stop 'clEnqueueNDRangeKernel'
  ierr=clFinish(cmd_queue)
  if (ierr.ne.0) stop 'clFinish'

  print '(a)','sent to device:'
  do i=1,10
     print '(2e12.5)',vec1(i),vec2(i)
  enddo

! read the resulting vector from device memory
  ierr=clEnqueueReadBuffer(cmd_queue,cl_vec2, &
          CL_TRUE,0_8,size_in_bytes,C_LOC(vec2), &
          0,C_NULL_PTR,C_NULL_PTR)
  if (ierr.ne.0) stop 'clEnqueueReadBuffer'

  print '(a)','retrieved from device:'
  do i=1,10
     print '(12x,e12.5)',vec2(i)
  enddo

  ierr=clReleaseKernel(kernel)
  if (ierr.ne.0) stop 'clReleaseKernel'
  ierr=clReleaseCommandQueue(cmd_queue)
  if (ierr.ne.0) stop 'clReleaseCommandQueue'
  ierr=clReleaseContext(context)
  if (ierr.ne.0) stop 'clReleaseContext'

contains

function check_status(text, ierr, verbose)
  use clfortran
  implicit none
  integer :: check_status
  character(len=*), intent(in) :: text
  integer, intent(in) :: ierr
  logical, optional, intent(in) :: verbose
  if(ierr /= CL_SUCCESS)then
    write(*,'("Hit error: ",(A),": ",(A))') text, "OCL error string goes here" !OCL_GetErrorString(ierr)
    stop
  end if
  if(verbose)then
    write(*,'("Called ",(A)," OK\n")') text 
  end if
  check_status = 0
end function check_status
  
end program sum






! Argument lists need to be handled with care. Neither the C nor the Fortran compiler checks for mismatched argument types, or even a mismatch in the number of arguments. Some bizarre run-time errors therefore arise. Keep in mind that Fortran passes arguments by reference whereas C passes arguments by value. When you put a variable name into a function call from Fortran, the corresponding C function receives a pointer to that variable. Similarly, when calling a Fortran subroutine from C, you must explicitly pass addresses rather than values in the argument list.

! When passing arrays, remember that C arrays start with subscript zero. Fortran stores multidimensional arrays in column-major order (first index varies fastest) whereas C stores them in row-major order (last index varies fastest).

! Passing character strings is a special problem. C knows it has come to the end of string when it hits a null character, but Fortran uses the declared length of a string. The C-Fortran interface provides an extra argument for each character string in a C argument list to receive the declared length of a string when called from Fortran. Consider the following Fortran fragment:



