#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Header for the C version of the kernel */
#include "continuity.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define MAX_SOURCE_SIZE (0x100000)
/** Maximum number of OpenCL devices we will query */
#define MAX_DEVICES 4
#define VERBOSE 1

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

void check_status(char *text, cl_int err){
  if(err != CL_SUCCESS){
    fprintf(stderr, "Hit error: %s: %s\n", text, OCL_GetErrorString(err));
    exit(1);
  }
  if(VERBOSE){
    fprintf(stdout, "Called %s OK\n", text); 
  }
}

/** Compute a checksum for a double precision array of nx*ny values */
double checksum(double *array, int nx, int ny, int xstart, int ystart){
  int i, j, jidx;
  double sum = 0.0;
  for(j=ystart; j<ny; j++){
    jidx = j*nx;
    for(i=xstart; i<nx; i++){
      sum += array[i+jidx];
    }
  }
  return sum;
}

int main(){
  /** The version of OpenCL supported by the selected device */
  int cl_version;
  cl_device_id device_ids[MAX_DEVICES];
  cl_device_id *device;
  int idev;
  cl_context context = NULL;
  cl_command_queue command_queue = NULL;
  cl_mem ssha_device = NULL;
  cl_mem sshn_device = NULL;
  cl_mem sshn_u_device = NULL;
  cl_mem sshn_v_device = NULL;
  cl_mem hu_device = NULL;
  cl_mem hv_device = NULL;
  cl_mem un_device = NULL;
  cl_mem vn_device = NULL;
  cl_mem e12t_device = NULL;
  cl_program program = NULL;
  cl_kernel kernel = NULL;
  cl_platform_id platform_ids[MAX_DEVICES];
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  cl_int ret;
  cl_event event;
  cl_command_queue_properties queue_properties;

  cl_int nx = 128;
  cl_int ny = 128;
  int ji, jj, idx;
  int buff_size;
  cl_double *ssha, *sshn, *sshn_u, *sshn_v;
  cl_double *hu, *hv, *un, *vn;
  cl_double rdt;
  cl_double *e12t;
  cl_double dep_const = 2.0;
  cl_double dx = 0.5, dy = 0.5;
  /** For computing checksums for validation */
  double sum;
  rdt = 0.5;

  /*------------------------------------------------------------*/
  /* OpenCL initialisation */
  FILE *fp;
  char fileName[] = "./continuity_kern.c";
  char *source_str;
  size_t source_size;

  /* Load the source code containing the kernel*/
  fp = fopen(fileName, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel %s.\n", fileName);
    exit(1);
  }
  source_str = (char*)malloc(MAX_SOURCE_SIZE);
  source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose( fp );
  
  /* Get Platform and Device Info */
  ret = clGetPlatformIDs(MAX_DEVICES, platform_ids, &ret_num_platforms);
  check_status("clGetPlatformIDs", ret);
  fprintf(stdout, "Have %d platforms.\n", ret_num_platforms);
  char result_str[128], version_str[128];
  cl_device_fp_config fp_config;
  cl_device_type device_type;
  size_t result_len;
  for(idev=0;idev<ret_num_platforms;idev++){
    ret = clGetPlatformInfo(platform_ids[idev],
			    CL_PLATFORM_NAME,
			    (size_t)128, (void *)result_str, &result_len);
    fprintf(stdout, "Platform %d (id=%d) is: %s\n",
	    idev, platform_ids[idev], result_str);
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
			  (size_t)128, &version_str, &result_len);
    ret = clGetDeviceInfo(device_ids[idev], CL_DEVICE_DOUBLE_FP_CONFIG,
			  (size_t)(sizeof(cl_device_fp_config)),
			  &fp_config, &result_len);
    fprintf(stdout, "Device %d is: %s, type=%d, version=%s\n",
	    idev, result_str, (int)(device_type), version_str);
    if((int)fp_config == 0){
      fprintf(stdout, "             double precision NOT supported\n");
    }
    else{
      fprintf(stdout, "             double precision supported\n");
    }
  }
  /* Choose device 0 */
  idev = 0;
  device = &(device_ids[idev]);

  /* Check what version of OpenCL is supported */
  if(strstr(version_str, "OpenCL 1.2")){
    cl_version = 12;
  }
  else if(strstr(version_str, "OpenCL 2.0")){
    cl_version = 20;
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
  context = clCreateContext(cl_props, 1, device, NULL, NULL, &ret);
  check_status("clCreateContext", ret);
  
  /* Create Command Queue with properties set to NULL */
  if(cl_version == 12){
    /* NVIDIA only support OpenCL 1.2 so we get a seg. fault if we attempt
     to call the ...WithProperties version of this routine */
    command_queue = clCreateCommandQueue(context,
					 *device,
					 0, &ret);
    check_status("clCreateCommandQueue", ret);
  }
  else if(cl_version == 20){
    command_queue = clCreateCommandQueueWithProperties(context,
						       *device,
						       NULL, &ret);
    check_status("clCreateCommandQueueWithProperties", ret);
  }

  /* Create Kernel Program from the source */
  program = clCreateProgramWithSource(context, 1, (const char **)&source_str,
				      (const size_t *)&source_size, &ret);
  check_status("clCreateProgramWithSource", ret);

  /* Build Kernel Program */
  ret = clBuildProgram(program, 1, device, NULL, NULL, NULL);
  if(ret == CL_BUILD_PROGRAM_FAILURE){
    char *build_log;
    size_t log_size;
    clGetProgramBuildInfo(program, *device, CL_PROGRAM_BUILD_LOG,
			  0, NULL, &log_size);
    build_log = malloc(log_size+1);
    clGetProgramBuildInfo(program, *device, CL_PROGRAM_BUILD_LOG,
			  log_size, build_log, NULL);
    build_log[log_size] = '\0';
    fprintf(stderr, "%s\n", build_log);
    exit(1);
  }
  check_status("clBuildProgram", ret);

  /* Create OpenCL Kernel */
  kernel = clCreateKernel(program, "continuity_code", &ret);
  check_status("clCreateKernel", ret);

  /* Create Device Memory Buffers */
  buff_size = nx*ny*sizeof(cl_double);
  ssha_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);
  sshn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);
  sshn_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);
  sshn_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);
  hu_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);
  hv_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);
  un_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);
  vn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);
  e12t_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);

  fprintf(stdout, "Created device buffers OK\n");

  /* Set OpenCL Kernel Parameters */
  ret = clSetKernelArg(kernel, 0, sizeof(cl_int), (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&ssha_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&un_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&vn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 9, sizeof(cl_double), (void *)&rdt);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&e12t_device);
  check_status("clSetKernelArg", ret);


  /*------------------------------------------------------------*/
  /* Field initialisation on host */
  
  ssha = malloc(buff_size);
  sshn = malloc(buff_size);
  sshn_u = malloc(buff_size);
  sshn_v = malloc(buff_size);
  hu = malloc(buff_size);
  hv = malloc(buff_size);
  un = malloc(buff_size);
  vn = malloc(buff_size);
  e12t = malloc(buff_size);

  for(jj=0;jj<ny;jj++){
    for(ji=0;ji<nx;ji++){
      idx = jj*nx + ji;
      hu[idx] = dep_const;
      hv[idx] = dep_const;
      un[idx] = ji;
      vn[idx] = jj;
      sshn_u[idx] = 0.0;
      sshn_v[idx] = 2.0;
      sshn[idx] = 1.0;
      e12t[idx] = dx*dy;
      ssha[idx] = 0.0;
    }
  }

  /* Copy data to device synchronously */
  ret = clEnqueueWriteBuffer(command_queue, ssha_device, 1, 0,
			     (size_t)buff_size, (void *)ssha, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, sshn_device, 1, 0,
			     (size_t)buff_size, (void *)sshn, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, sshn_u_device, 1, 0,
			     (size_t)buff_size, (void *)sshn_u, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, sshn_v_device, 1, 0,
			     (size_t)buff_size, (void *)sshn_v, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, hu_device, 1, 0,
			     (size_t)buff_size, (void *)hu, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, hv_device, 1, 0,
			     (size_t)buff_size, (void *)hv, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, un_device, 1, 0,
			     (size_t)buff_size, (void *)un, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, vn_device, 1, 0,
			     (size_t)buff_size, (void *)vn, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e12t_device, 1, 0,
			     (size_t)buff_size, (void *)e12t, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);

  /*------------------------------------------------------------*/
  /* Run the kernel */
  size_t global_size[2] = {nx, ny}; 
  clEnqueueNDRangeKernel(command_queue, kernel, 2, 0,
			 global_size, NULL, 0,0,0);
  
  /* Run the kernel on the CPU */
  for(jj=1;jj<ny;jj++){
    for(ji=1;ji<nx;ji++){
      continuity_code(ji, jj, nx,                     
		      ssha, sshn, sshn_u, sshn_v, hu, hv,
		      un, vn, rdt, e12t);
    }
  }
  printf("ssha[1,1] [1,2] = %e, %e\n", ssha[nx+1], ssha[nx+2]);

  sum = checksum(ssha, nx, ny, 1, 1);
  fprintf(stdout, "Checksum on CPU = %e\n", sum);

  /* Copy data back from device, synchronously */
  clEnqueueReadBuffer(command_queue, ssha_device, 1, 0,
		      (size_t)buff_size, (void *)ssha, 0, NULL, NULL);
  check_status("clEnqueueReadBuffer", ret);

  /* Compute and output a checksum */
  sum = checksum(ssha, nx, ny, 1, 1);
  printf("ssha[1,1] [1,2] = %e, %e\n", ssha[nx+1], ssha[nx+2]);
  fprintf(stdout, "Checksum from OpenCL = %e\n", sum);

  /* Clean up */
  ret = clFlush(command_queue);
  ret = clFinish(command_queue);
  ret = clReleaseKernel(kernel);
  ret = clReleaseProgram(program);
  ret = clReleaseMemObject(ssha_device);
  ret = clReleaseMemObject(sshn_device);
  ret = clReleaseMemObject(sshn_u_device);
  ret = clReleaseMemObject(sshn_v_device);
  ret = clReleaseMemObject(hu_device);
  ret = clReleaseMemObject(hv_device);
  ret = clReleaseMemObject(un_device);
  ret = clReleaseMemObject(vn_device);
  ret = clReleaseMemObject(e12t_device);

  ret = clReleaseCommandQueue(command_queue);
  ret = clReleaseContext(context);

  free(source_str);

}
