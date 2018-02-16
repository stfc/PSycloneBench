module opencl_utils_mod
  implicit none

  integer, parameter :: CL_UTIL_STR_LEN = 64

contains
  
  subroutine init_device(device, version_str, context)
    use clfortran
    use iso_c_binding
    integer(c_intptr_t), intent(inout) :: device, context 
    character(len=CL_UTIL_STR_LEN), intent(inout) :: version_str

    integer :: iplatform, idevice, iallocerr, irec
    integer(c_intptr_t), target :: ctx_props(3)
    integer(c_int32_t), target :: device_cu
    integer(c_size_t) :: iret, zero_size = 0
    integer(c_int32_t) :: ierr, num_devices, num_platforms
    integer(c_intptr_t), allocatable, target :: &
       platform_ids(:), device_ids(:)
    character(len=1,kind=c_char), allocatable, target :: device_name(:)

    ierr = clGetPlatformIDs(0, C_NULL_PTR, num_platforms)
    call check_status('clGetPlatformIDs', ierr)
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
    call check_status('clGetPlatformIDs', ierr)

    ! Get device IDs only for platform 1
    iplatform=1

    ierr=clGetDeviceIDs(platform_ids(iplatform), CL_DEVICE_TYPE_ALL, &
                        0, C_NULL_PTR, num_devices)
    call check_status('clGetDeviceIDs', ierr)
    if (num_devices < 1)then
       stop 'Failed to find any OpenCL devices'
    end if
    print '(a,i2)','Num Devices: ',num_devices

    allocate(device_ids(num_devices), stat=iallocerr)
    if (iallocerr.ne.0) stop 'memory allocation error'

    ! whenever "&" appears in C subroutine (address-off) call,
    ! then C_LOC has to be used in Fortran
    ierr = clGetDeviceIDs(platform_ids(iplatform), CL_DEVICE_TYPE_ALL, &
                          num_devices, C_LOC(device_ids), num_devices)
    call check_status('clGetDeviceIDs', ierr)

    ! Get device info only for device 1
    idevice=1
    device = device_ids(idevice)

    ierr=clGetDeviceInfo(device_ids(idevice), &
                         CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(device_cu), &
                         C_LOC(device_cu), iret)
    call check_status('clGetDeviceInfo', ierr)
    ierr=clGetDeviceInfo(device_ids(idevice), &
                         CL_DEVICE_NAME, zero_size, C_NULL_PTR,iret)
    call check_status('clGetDeviceInfo', ierr)

    allocate(device_name(iret), stat=iallocerr)
    if (iallocerr.ne.0) stop 'allocate'

    ierr=clGetDeviceInfo(device_ids(idevice), CL_DEVICE_NAME, &
                         sizeof(device_name), C_LOC(device_name), iret)
    if (ierr.ne.CL_SUCCESS) stop 'clGetDeviceInfo'

    write (*,'(a,i2,a,i3,a)',advance='no') &
         ' Device (#',idevice,', Compute Units: ',device_cu,') - '
    print *,device_name(1:iret)
    deallocate(device_name)

    print '(a,i2,a)', 'Creating context for: ', num_devices,' devices'
    print '(a,i2)', 'for platform: ',iplatform
    ctx_props(1) = CL_CONTEXT_PLATFORM
    ctx_props(2) = platform_ids(iplatform)
    ctx_props(3) = 0
    context = clCreateContext(C_LOC(ctx_props), num_devices, &
                              C_LOC(device_ids), C_NULL_FUNPTR, C_NULL_PTR, &
                              ierr)
    call check_status('clCreateContext', ierr)

  end subroutine init_device

  subroutine check_status(text, ierr)
    use clfortran
    implicit none
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

  end subroutine check_status
  
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

end module opencl_utils_mod
