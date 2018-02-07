#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nemolite2d_utils.h"
#include "opencl_utils.h"
#include "timing.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

/** Maximum number of OpenCL devices we will query */
#define MAX_DEVICES 4

/** The number of concurrent kernels in our design */
#define NUM_KERNELS 2


/** Top-level driver program. Queries the hardware to find OpenCL devices,
    creates OpenCL kernels and runs them. Also runs the same kernels on
    the CPU and compares the results. */
int main(){
  /** The version of OpenCL supported by the selected device */
  int cl_version;
  cl_device_id device_ids[MAX_DEVICES];
  cl_device_id *device;
  int idev;
  cl_context context = NULL;
  /* In order to run multiple kernels concurrently, we must have
     multiple command queues (since the Intel OpenCL SDK only supports
     in-order queues) */
  cl_command_queue command_queue[NUM_KERNELS];
  cl_command_queue_properties queue_properties = 0;
  cl_bool profiling_enabled = 0;

  /* Buffers on the device */
  cl_mem ssha_device = NULL;
  cl_mem sshn_device = NULL;

  cl_program program = NULL;

  /* Our OpenCL kernel objects */
  cl_kernel cont_kernel = NULL;
  cl_kernel next_sshu_kernel  = NULL;

  cl_platform_id platform_ids[MAX_DEVICES];
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  cl_int ret;

  /** Default problem size. May be overridden by setting
      NEMOLITE2D_N{X,Y} environment variables. */
  char *env_string;
  /** Name of file containing single, compiled image for multiple
      kernels */
  char *image_file = NULL;
  cl_int nx = 127;
  cl_int ny = 127;
    int xstart, xstop, ystart, ystop;
  /** Our time-step index (passed into BCs kernel) */
  cl_int istep;
  /** Number of time-steps to do. May be overridden by setting
      NEMOLITE2D_NSTEPS environment variable. */
  cl_int nsteps = 2000;
  int ji, jj, idx;
  int buff_size;
  /** Sea-surface height */
  cl_double *ssha, *ssha_u, *ssha_v, *sshn, *sshn_u, *sshn_v;
  cl_double *hu, *hv, *ht, *un, *vn, *ua, *va;
  cl_double *gphiu, *gphiv;
  cl_double *e1u, *e1v, *e1t, *e2u, *e2v, *e2t, *e12t, *e12u, *e12v;
  /** T-point mask */
  cl_int *tmask;
  cl_double dep_const = 100.0;
  /** Horizontal grid resolution */
  cl_double dx = 1000.0, dy = 1000.0;
  /** For computing checksums for validation */
  double cpu_sum[3], ocl_sum[3];
  /** Time step (s) */
  cl_double rdt = 20.0;
  /** horiz. kinematic viscosity coeff. */
  cl_double visc = 0.1;
  /** Coefficient of bottom friction */
  cl_double cbfr = 0.00015;

  TimerInit();

  /*------------------------------------------------------------*/
  /* Run-time configuration of the benchmark */
  if(getenv("NEMOLITE2D_PROFILING")){
    profiling_enabled = CL_TRUE;
    /* We create the OpenCL command queue with profiling enabled */
    queue_properties = CL_QUEUE_PROFILING_ENABLE;
  }
  if( (env_string = getenv("NEMOLITE2D_NSTEPS")) ){
    if(sscanf(env_string, "%d", &nsteps) != 1){
      fprintf(stderr,
	      "Error parsing NEMOLITE2D_NSTEPS environment variable (%s)\n",
	      env_string);
      exit(1);
    }
    if(nsteps < 1){
      fprintf(stderr, "Error, NEMOLITE2D_NSTEPS must be >= 1 but got %d\n",
	      nsteps);
      exit(1);
    }
  }
  if( (env_string = getenv("NEMOLITE2D_NX")) ){
    if(sscanf(env_string, "%d", &nx) != 1){
      fprintf(stderr, "Error parsing NEMOLITE2D_NX environment variable (%s)\n",
	      env_string);
      exit(1);
    }
  }
  if( (env_string = getenv("NEMOLITE2D_NY")) ){
    if(sscanf(env_string, "%d", &ny) != 1){
      fprintf(stderr, "Error parsing NEMOLITE2D_NY environment variable (%s)\n",
	      env_string);
      exit(1);
    }
  }
  /* Check to see whether we should get our kernels from a single image file */
  if( (env_string = getenv("NEMOLITE2D_SINGLE_IMAGE")) ){
    if(sscanf(env_string, "%ms", &image_file) != 1){
      fprintf(stderr,
	      "Error parsing NEMOLITE2D_SINGLE_IMAGE environment "
	      "variable (%s)\n",
	      env_string);
      exit(1);
    }
  }
  /* Extend domain by one in each dimension to allow for staggering */
  nx += 1;
  ny += 1;

  /*------------------------------------------------------------*/
  /* OpenCL initialisation */

  TimerStart("OCL Init");
  
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
			  (size_t)128, &version_str, &result_len);
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
    fprintf(stdout, "             max size of each work item dimension: %ld %ld\n", max_sizes[0], max_sizes[1]);
    fprintf(stdout, "             max compute units = %ld\n", (long)max_units);
    free(max_sizes);
  }
  /* Choose device 0 */
  idev = 0;
  device = &(device_ids[idev]);

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
  context = clCreateContext(cl_props, 1, device, NULL, NULL, &ret);
  check_status("clCreateContext", ret);
  
  /* Create Command Queue with properties set to NULL */
  for(ji=0; ji<NUM_KERNELS; ji++){
    /* The Intel/Altera OpenCL SDK is only version 1.0 */
    /* NVIDIA only support OpenCL 1.2 so we get a seg. fault if we attempt
       to call the ...WithProperties version of this routine */
    command_queue[ji] = clCreateCommandQueue(context,
					     *device,
					     queue_properties, &ret);
    check_status("clCreateCommandQueue", ret);
  }

  /* Create OpenCL Kernels and associated event objects (latter used
   to obtain detailed timing information). */
  cl_event cont_evt;
  if(image_file){
    cont_kernel = get_kernel(&context, device, version_str,
			     image_file,
			     "channel_write");
  }
  else{
    fprintf(stderr, "Please set NEMOLITE2D_SINGLE_IMAGE to point to the "
	    ".aocx containing the compiled kernels\n");
    exit(1);
  }
  cl_event next_sshu_evt;
  if(image_file){
    next_sshu_kernel = get_kernel(&context, device, version_str,
				  image_file,
				  "channel_read");
  }
  else{
    fprintf(stderr, "Please set NEMOLITE2D_SINGLE_IMAGE to point to the "
	    ".aocx containing the compiled kernels\n");
    exit(1);
  }

  /* Create Device Memory Buffers */
  int num_buffers = 0;
  buff_size = nx*ny*sizeof(cl_double);
  ssha_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			       NULL, &ret);
  num_buffers++;
  sshn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			       NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);

  fprintf(stdout, "Created %d device buffers OK\n", num_buffers);

  int arg_idx = 0;
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_int), (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_int), (void *)&ny);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ssha_device);
  check_status("clSetKernelArg", ret);
  arg_idx = 0;
  ret = clSetKernelArg(next_sshu_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshu_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&ny);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_device);
  check_status("clSetKernelArg", ret);
  TimerStop();

  /*------------------------------------------------------------*/
  /* Field initialisation on host */
  
  init_fields(nx, ny, dx, dy, dep_const,
	      &xstart, &xstop, &ystart, &ystop,
	      &ssha,  &ssha_u, &ssha_v,
	      &sshn,  &sshn_u, &sshn_v,
	      &hu,  &hv,  &ht,
	      &un,  &vn,  &ua, &va,
	      &gphiu, &gphiv,
	      &e1u,  &e1v,  &e1t,
	      &e2u,  &e2v,  &e2t,
	      &e12t, &e12u, &e12v,
	      &tmask);
  
  write_ifield("tmask.dat", nx, ny, 0, 0, tmask);
  
  /*------------------------------------------------------------*/
  /* Copy data to device synchronously */
  TimerStart("Write buffers to device");

  /* Create an array to store the event associated with each write
     to the device */
  cl_event *write_events = (cl_event*)malloc(num_buffers*sizeof(cl_event));
  int buf_idx = 0;
  ret = clEnqueueWriteBuffer(command_queue[0], ssha_device, 1, 0,
			     (size_t)buff_size, (void *)ssha, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], sshn_device, 1, 0,
			     (size_t)buff_size, (void *)sshn, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);

  ret = clWaitForEvents(num_buffers, write_events);
  check_status("clWaitForEvents", ret);

  TimerStop();
  
  /*------------------------------------------------------------*/
  /* Run the kernels */

  TimerStart("Time-stepping, OpenCL");
  
  for(istep=1; istep<=nsteps; istep++){

    ret = clEnqueueTask(command_queue[0], cont_kernel, 0,
			NULL, &cont_evt);
    check_status("clEnqueueTask(Continuity)", ret);

    /* Update of sshu field */
    ret = clEnqueueTask(command_queue[1], next_sshu_kernel, 0,
			NULL, &next_sshu_evt);
    check_status("clEnqueueTask(next-sshu)", ret);
  }

  ret = clWaitForEvents(1, &cont_evt);
  check_status("clWaitForEvents", ret);
  fprintf(stdout, "First kernel done!\n");
  /* Block on the execution of the last kernel */
  ret = clWaitForEvents(1, &next_sshu_evt);
  check_status("clWaitForEvents", ret);

  TimerStop();
    
  /* Copy data back from device, synchronously. */
  cl_event read_events[3];
  int nread = 0;
  ret = clEnqueueReadBuffer(command_queue[0], ssha_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)ssha, 0, NULL,
			    &(read_events[0]));
  nread++;
  check_status("clEnqueueReadBuffer", ret);
  clWaitForEvents(nread, read_events);
  check_status("clWaitForEvents", ret);

  /* Dump final fields computed on OpenCL device */
  write_field("ssha_ocl.dat", nx, ny, 0, 0, ssha);
  write_field("ua_ocl.dat", nx, ny, 0, 0, ua);
  write_field("va_ocl.dat", nx, ny, 0, 0, va);

  /* Extract profiling info from the OpenCL runtime. Note that this will
     be the time taken during the most recent execution of each kernel. */
  if(profiling_enabled){
    fprintf(stdout, "Time spent in Continuity kern = %e us\n",
  	    duration_ns(cont_evt)*0.001);
    fprintf(stdout, "Time spent in next-sshu kern = %e us\n",
            duration_ns(next_sshu_evt)*0.001);
  }

  /* Generate report on host timings */
  TimerReport();

  /* Clean up */
  for(ji=0; ji<NUM_KERNELS; ji++){
    ret = clFlush(command_queue[ji]);
    ret = clFinish(command_queue[ji]);
  }
  ret = clReleaseKernel(cont_kernel);
  ret = clReleaseKernel(next_sshu_kernel);

  ret = clReleaseProgram(program);
  ret = clReleaseMemObject(ssha_device);
  ret = clReleaseMemObject(sshn_device);

  for(ji=0; ji<NUM_KERNELS; ji++){
    ret = clReleaseCommandQueue(command_queue[ji]);
  }
  ret = clReleaseContext(context);

  if(image_file){
    free(image_file);
  }
}
