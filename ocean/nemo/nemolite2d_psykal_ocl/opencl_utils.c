#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include "opencl_utils.h"

/** Maximum number of OpenCL devices we will query */
#define MAX_DEVICES 4

#define MAX_SOURCE_SIZE (0x100000)
#define VERBOSE 1

/** Query the available OpenCL devices and choose one */
void init_device(cl_device_id *device,
		 char *version_str,
		 cl_context *context)
{
  /** The version of OpenCL supported by the selected device */
  int cl_version;
  cl_device_id device_ids[MAX_DEVICES];
  cl_platform_id platform_ids[MAX_DEVICES];
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  cl_int ret;
  int idev;

  /* Get Platform and Device Info */
  ret = clGetPlatformIDs(MAX_DEVICES, platform_ids, &ret_num_platforms);
  check_status("clGetPlatformIDs", ret);
  fprintf(stdout, "Have %d platforms.\n", ret_num_platforms);
  char result_str[128];
  cl_device_fp_config fp_config;
  cl_device_type device_type;
  size_t result_len;
  for(idev=0;idev<ret_num_platforms;idev++){
    ret = clGetPlatformInfo(platform_ids[idev],
			    CL_PLATFORM_NAME,
			    (size_t)128, (void *)result_str, &result_len);
    fprintf(stdout, "Platform %d (id=%ld) is: %s\n",
	    idev, (long)(platform_ids[idev]), result_str);
  }

  ret = clGetDeviceIDs(platform_ids[0], CL_DEVICE_TYPE_DEFAULT, MAX_DEVICES,
		       device_ids, &ret_num_devices);
  check_status("clGetDeviceIDs", ret);
  fprintf(stdout, "Have %d devices\n", ret_num_devices);
  for (idev=0; idev<ret_num_devices; idev++){
    ret = clGetDeviceInfo(device_ids[idev], CL_DEVICE_NAME,
			  (size_t)128, result_str, &result_len);
    ret = clGetDeviceInfo(device_ids[idev], CL_DEVICE_TYPE,
			  (size_t)(sizeof(cl_device_type)), &device_type,
				   &result_len);
    ret = clGetDeviceInfo(device_ids[idev], CL_DEVICE_VERSION,
			  (size_t)128, version_str, &result_len);
#ifdef CL_DEVICE_DOUBLE_FP_CONFIG
    ret = clGetDeviceInfo(device_ids[idev], CL_DEVICE_DOUBLE_FP_CONFIG,
    			  (size_t)(sizeof(cl_device_fp_config)),
    			  &fp_config, &result_len);
#else
    /* The Intel/Altera OpenCL SDK is only version 1.0 and that doesn't
     have the CL_DEVICE_DOUBLE_FP_CONFIG property */
    fp_config = 0;
#endif
    size_t wg_size;
    ret = clGetDeviceInfo(device_ids[idev], CL_DEVICE_MAX_WORK_GROUP_SIZE,
			  (size_t)(sizeof(size_t)),
			  &wg_size, &result_len);
    cl_uint ndims;
    ret = clGetDeviceInfo(device_ids[idev], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
			  sizeof(cl_uint), &ndims, &result_len);
    size_t *max_sizes;
    max_sizes = (size_t*)malloc(ndims*sizeof(size_t));
    ret = clGetDeviceInfo(device_ids[idev], CL_DEVICE_MAX_WORK_ITEM_SIZES,
			  ndims*sizeof(size_t), max_sizes, &result_len);
    cl_uint max_units;
    ret = clGetDeviceInfo(device_ids[idev], CL_DEVICE_MAX_COMPUTE_UNITS,
			  sizeof(cl_uint), &max_units, &result_len);

    fprintf(stdout, "Device %d is: %s, type=%d, version=%s\n",
	    idev, result_str, (int)(device_type), version_str);
    if((int)fp_config == 0){
      fprintf(stdout, "             double precision NOT supported\n");
    }
    else{
      fprintf(stdout, "             double precision supported\n");
    }
    fprintf(stdout, "             max work group size = %ld\n", (long)wg_size);
    fprintf(stdout, "             max size of each work item dimension: "
	    "%ld %ld\n", max_sizes[0], max_sizes[1]);
    fprintf(stdout, "             max compute units = %ld\n", (long)max_units);
    free(max_sizes);
  }
  /* Choose device 0 */
  idev = 0;
  *device = device_ids[idev];

  /* Check what version of OpenCL is supported */
  if(strstr(version_str, "OpenCL 1.0")){
    cl_version = 100;
  }
  else if(strstr(version_str, "OpenCL 1.2")){
    cl_version = 120;
  }
  else if(strstr(version_str, "OpenCL 2.0")){
    cl_version = 200;
  }
  else{
    fprintf(stderr, "Unsupported OpenCL version: %s\n", version_str);
    exit(1);
  }

  /* Create OpenCL context for just 1 device */
  cl_context_properties cl_props[3];
  /* The list of properties for this context is zero terminated */
  cl_props[0] = CL_CONTEXT_PLATFORM;
  cl_props[1] = (cl_context_properties)(platform_ids[0]);
  cl_props[2] = 0;

  *context = clCreateContext(cl_props, 1, device, NULL, NULL, &ret);
  check_status("clCreateContext", ret);

}


const char* OCL_GetErrorString(cl_int error)
{
    switch (error)
    {
    case CL_SUCCESS:
        return "CL_SUCCESS";
    case CL_DEVICE_NOT_FOUND:
        return "CL_DEVICE_NOT_FOUND";
    case CL_DEVICE_NOT_AVAILABLE:
        return "CL_DEVICE_NOT_AVAILABLE";
    case CL_COMPILER_NOT_AVAILABLE:
        return "CL_COMPILER_NOT_AVAILABLE";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:
        return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
    case CL_OUT_OF_RESOURCES:
        return "CL_OUT_OF_RESOURCES";
    case CL_OUT_OF_HOST_MEMORY:
        return "CL_OUT_OF_HOST_MEMORY";
    case CL_PROFILING_INFO_NOT_AVAILABLE:
        return "CL_PROFILING_INFO_NOT_AVAILABLE";
    case CL_MEM_COPY_OVERLAP:
        return "CL_MEM_COPY_OVERLAP";
    case CL_IMAGE_FORMAT_MISMATCH:
        return "CL_IMAGE_FORMAT_MISMATCH";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:
        return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
    case CL_BUILD_PROGRAM_FAILURE:
        return "CL_BUILD_PROGRAM_FAILURE";
    case CL_MAP_FAILURE:
        return "CL_MAP_FAILURE";
    case CL_INVALID_VALUE:
        return "CL_INVALID_VALUE";
    case CL_INVALID_DEVICE_TYPE:
        return "CL_INVALID_DEVICE_TYPE";
    case CL_INVALID_PLATFORM:
        return "CL_INVALID_PLATFORM";
    case CL_INVALID_DEVICE:
        return "CL_INVALID_DEVICE";
    case CL_INVALID_CONTEXT:
        return "CL_INVALID_CONTEXT";
    case CL_INVALID_QUEUE_PROPERTIES:
        return "CL_INVALID_QUEUE_PROPERTIES";
    case CL_INVALID_COMMAND_QUEUE:
        return "CL_INVALID_COMMAND_QUEUE";
    case CL_INVALID_HOST_PTR:
        return "CL_INVALID_HOST_PTR";
    case CL_INVALID_MEM_OBJECT:
        return "CL_INVALID_MEM_OBJECT";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
        return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
    case CL_INVALID_IMAGE_SIZE:
        return "CL_INVALID_IMAGE_SIZE";
    case CL_INVALID_SAMPLER:
        return "CL_INVALID_SAMPLER";
    case CL_INVALID_BINARY:
        return "CL_INVALID_BINARY";
    case CL_INVALID_BUILD_OPTIONS:
        return "CL_INVALID_BUILD_OPTIONS";
    case CL_INVALID_PROGRAM:
        return "CL_INVALID_PROGRAM";
    case CL_INVALID_PROGRAM_EXECUTABLE:
        return "CL_INVALID_PROGRAM_EXECUTABLE";
    case CL_INVALID_KERNEL_NAME:
        return "CL_INVALID_KERNEL_NAME";
    case CL_INVALID_KERNEL_DEFINITION:
        return "CL_INVALID_KERNEL_DEFINITION";
    case CL_INVALID_KERNEL:
        return "CL_INVALID_KERNEL";
    case CL_INVALID_ARG_INDEX:
        return "CL_INVALID_ARG_INDEX";
    case CL_INVALID_ARG_VALUE:
        return "CL_INVALID_ARG_VALUE";
    case CL_INVALID_ARG_SIZE:
        return "CL_INVALID_ARG_SIZE";
    case CL_INVALID_KERNEL_ARGS:
        return "CL_INVALID_KERNEL_ARGS";
    case CL_INVALID_WORK_DIMENSION:
        return "CL_INVALID_WORK_DIMENSION";
    case CL_INVALID_WORK_GROUP_SIZE:
        return "CL_INVALID_WORK_GROUP_SIZE";
    case CL_INVALID_WORK_ITEM_SIZE:
        return "CL_INVALID_WORK_ITEM_SIZE";
    case CL_INVALID_GLOBAL_OFFSET:
        return "CL_INVALID_GLOBAL_OFFSET";
    case CL_INVALID_EVENT_WAIT_LIST:
        return "CL_INVALID_EVENT_WAIT_LIST";
    case CL_INVALID_EVENT:
        return "CL_INVALID_EVENT";
    case CL_INVALID_OPERATION:
        return "CL_INVALID_OPERATION";
    case CL_INVALID_GL_OBJECT:
        return "CL_INVALID_GL_OBJECT";
    case CL_INVALID_BUFFER_SIZE:
        return "CL_INVALID_BUFFER_SIZE";
    case CL_INVALID_MIP_LEVEL:
        return "CL_INVALID_MIP_LEVEL";
    case CL_INVALID_GLOBAL_WORK_SIZE:
        return "CL_INVALID_GLOBAL_WORK_SIZE";
        // unknown
    default:
        return "unknown error code";
    }
}

void check_status(const char *text, cl_int err){
  if(err != CL_SUCCESS){
    fprintf(stderr, "Hit error: %s: %s\n", text, OCL_GetErrorString(err));
    exit(1);
  }
  if(VERBOSE){
    fprintf(stdout, "Called %s OK\n", text); 
  }
}

/** Creates an OpenCL kernel for the supplied context and device. If the
 device is an FPGA then the kernel must be pre-compiled. */
cl_kernel get_kernel(cl_context *context, cl_device_id *device,
		     const char *version_str,
		     const char *filename,
		     const char *kernel_name){
  /* Holds return value of calls to OpenCL API */
  cl_int ret;
  cl_kernel kernel = NULL;
  cl_program program;
  
  if( strstr(version_str, "FPGA SDK") ){
    program = get_binary_kernel(context, device, filename);
  }
  else{
    program = get_source_kernel(context, device, filename);
  }

  /* Create OpenCL Kernel */
  kernel = clCreateKernel(program, kernel_name, &ret);
  check_status("clCreateKernel", ret);

  return kernel;
}

/** Creates an OpenCL kernel by compiling it from the supplied source */
cl_program get_source_kernel(cl_context *context,
			     cl_device_id *device,
			     const char *filename){
  FILE *fp;
  char *source_str;
  size_t source_size;
  /* Holds return value of calls to OpenCL API */
  cl_int ret;

  /* Load the source code containing the kernel*/
  fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel source: %s.\n", filename);
    exit(1);
  }
  source_str = (char*)malloc(MAX_SOURCE_SIZE);
  source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose( fp );

  /* Create Kernel Program from the source */
  cl_program program = clCreateProgramWithSource(*context, 1,
						 (const char **)&source_str,
						 (const size_t *)&source_size,
						 &ret);
  check_status("clCreateProgramWithSource", ret);

  /* Build Kernel Program */
  char *build_options = "-cl-mad-enable -cl-fast-relaxed-math";
  ret = clBuildProgram(program, 1, device, build_options, NULL, NULL);
  if(ret == CL_BUILD_PROGRAM_FAILURE){
    char *build_log;
    size_t log_size;
    clGetProgramBuildInfo(program, *device, CL_PROGRAM_BUILD_LOG,
			  0, NULL, &log_size);
    build_log = (char *)malloc(log_size+1);
    clGetProgramBuildInfo(program, *device, CL_PROGRAM_BUILD_LOG,
			  log_size, build_log, NULL);
    build_log[log_size] = '\0';
    fprintf(stderr, "%s\n", build_log);
    free(build_log);
    exit(1);
  }
  check_status("clBuildProgram", ret);

  /* Clean up */
  free(source_str);

  return program;
}
  
cl_program get_binary_kernel(cl_context *context,
			     cl_device_id *device,
			     const char *filename){
  FILE *fp;
  const int num_binaries = 1;
  unsigned char *binary_buffers[num_binaries];
  size_t binary_sizes[num_binaries];
  cl_int binary_status[num_binaries];
  cl_int ret_codes[num_binaries];
  /* Modified filename of the kernel binary (as opposed to source) */
  char bname[STD_STRING_LEN];
  char *ptr;
  /* Holds return value of calls to OpenCL API */
  cl_int ret;
  /* Modify the name of the kernel to point to a pre-compiled .aocx file */
  strcpy(bname, filename);
  
  if(ptr = strstr(bname, ".c")){
    /* We've been given the name of the kernel source file. We change the
       suffix to get the name of the compiled version. */
    sprintf(ptr, ".aocx");
  }
  else if(!strstr(bname, ".aocx")){
    fprintf(stderr,
	    "ERROR: get_binary_kernel: supplied filename (%s) is not a c "
	    "source (.c) file or a compiled (.aocx) file\n", bname);
    exit(1);
  }

  /* Open and read the file containing the pre-compiled kernel */  
  fp = fopen(bname, "rb");
  if (!fp) {
    fprintf(stderr,
	    "ERROR: get_binary_kernel: Failed to load pre-compiled kernel file: "
	    "%s.\n", filename);
    exit(1);
  }
  fseek(fp, 0, SEEK_END);
  binary_sizes[0] = ftell(fp);
  binary_buffers[0] = (unsigned char*)malloc(
				  sizeof(unsigned char)*binary_sizes[0]);
  rewind(fp);
  fread(binary_buffers[0], binary_sizes[0], 1, fp);
  fclose(fp);
  fprintf(stdout, "Read %d bytes for binary %s\n", binary_sizes[0], bname);

  /* Create the program object from the loaded binary */
  cl_program program = clCreateProgramWithBinary(*context, 1, device,
						 binary_sizes, binary_buffers,
						 binary_status, &ret);
  check_status("clCreateProgramWithBinary", ret);

  /* Clean up */
  for(int ibuf=0; ibuf<num_binaries; ibuf++){
    free(binary_buffers[ibuf]);
  }

  return program;
}

/** Returns the duration of the supplied OpenCL event in nanoseconds.
 Requires OpenCL profiling to have been enabled on the queue that performed
 the operation to which the event corresponds. */
cl_ulong duration_ns(cl_event event){
  cl_ulong start_time, end_time;
  cl_int ret;

  ret = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START,
				sizeof(cl_ulong), (void*)&start_time, NULL);
  check_status("clGetEventProfilingInfo", ret);

  ret = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END,
				sizeof(cl_ulong), (void*)&end_time, NULL);
  check_status("clGetEventProfilingInfo", ret);

  return (end_time - start_time);
}
