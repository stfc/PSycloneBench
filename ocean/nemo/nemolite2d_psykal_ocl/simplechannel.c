#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nemolite2d_utils.h"
#include "opencl_utils.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

//#include "AOCLUtils/aocl_utils.h"

//using namespace aocl_utils;

/** The number of concurrent kernels in our design */
#define NUM_KERNELS 2

/** Top-level driver program. Queries the hardware to find OpenCL devices,
    creates OpenCL kernels and runs them. Also runs the same kernels on
    the CPU and compares the results. */
int main(){
  /** Pointer to the OpenCL device (hardware) that we will use */
  cl_device_id device = NULL;
  cl_context context = NULL;
  /* String holding info on chosen device */
  char version_str[128];

  /* In order to run multiple kernels concurrently, we must have
     multiple command queues (since the Intel OpenCL SDK only supports
     in-order queues) */
  cl_command_queue write_queue;
  cl_command_queue read_queue;
  cl_command_queue_properties queue_properties = 0;
  cl_bool profiling_enabled = 0;

  /* Buffers on the device */
  cl_mem ssha_device = NULL;
  cl_mem sshn_device = NULL;

  cl_program program = NULL;

  cl_int ret;

  /* Our OpenCL kernel objects */
  cl_kernel read_kernel = NULL;
  cl_kernel write_kernel = NULL;

  /* Ptr used when getting env variables */
  char *env_string;
  /** Name of file containing single, compiled image for multiple
      kernels */
  char *image_file = NULL;

  /** Default problem size. May be overridden by setting
      NEMOLITE2D_N{X,Y} environment variables. */
  cl_int nx = 127;
  cl_int ny = 127;
  cl_int n_points, n_points_out;
  int ji, jj, idx;
  int buff_size;
  /** Sea-surface height */
  cl_int *ssha, *sshn;
  /** For computing checksums for validation */
  double cpu_sum[3], ocl_sum[3];

  /*------------------------------------------------------------*/
  /* Run-time configuration of the benchmark */
  if(getenv("NEMOLITE2D_PROFILING")){
    profiling_enabled = CL_TRUE;
    /* We create the OpenCL command queue with profiling enabled */
    queue_properties = CL_QUEUE_PROFILING_ENABLE;
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

  init_device(&device, version_str, &context);

  /* Create Command Queue with properties set to NULL */
  /* The Intel/Altera OpenCL SDK is only version 1.0 */
  /* NVIDIA only support OpenCL 1.2 so we get a seg. fault if we attempt
     to call the ...WithProperties version of this routine */
  write_queue = clCreateCommandQueue(context, device,
				     CL_QUEUE_PROFILING_ENABLE, &ret);
  check_status("clCreateCommandQueue", ret);
  read_queue = clCreateCommandQueue(context, device,
				    CL_QUEUE_PROFILING_ENABLE, &ret);
  check_status("clCreateCommandQueue", ret);
  
  /* Create the OpenCL program object */
  if(image_file){
    program = get_program(context, &device, version_str, image_file);
  }
  else{
    fprintf(stderr, "Please set NEMOLITE2D_SINGLE_IMAGE to point to the "
	    ".aocx file containing the compiled kernels\n");
    exit(1);
  }

  /* Create the kernels associated with the program object */
  write_kernel = clCreateKernel(program, "channel_write", &ret);
  check_status("clCreateCommandQueue", ret);
  read_kernel = clCreateKernel(program, "channel_read", &ret);
  check_status("clCreateCommandQueue", ret);

  buff_size = nx*ny*sizeof(cl_int);

  ssha = (cl_int*)malloc(buff_size);
  sshn = (cl_int*)malloc(buff_size);

  /* Create Device Memory Buffers */
  int num_buffers = 0;
  ssha_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			       NULL, &ret);
  num_buffers++;
  sshn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			       NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);

  fprintf(stdout, "Created %d device buffers OK\n", num_buffers);
  
  /*------------------------------------------------------------*/
  /* Set kernel arguments */
  
  ret = clSetKernelArg(write_kernel, 0, sizeof(cl_int), (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(write_kernel, 1, sizeof(cl_int), (void *)&ny);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(write_kernel, 2, sizeof(cl_mem),
		       (void *)&ssha_device);
  check_status("clSetKernelArg", ret);

  ret = clSetKernelArg(read_kernel, 0, sizeof(cl_int), (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(read_kernel, 1, sizeof(cl_int), (void *)&ny);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(read_kernel, 2, sizeof(cl_mem), (void *)&sshn_device);
  check_status("clSetKernelArg", ret);

  /*------------------------------------------------------------*/
  /* Copy data to device synchronously */

  /* Set a recognisable first value in the input buffer */
  ssha[0] = 99;
  
  /* Create an array to store the event associated with each write
     to the device */
  cl_event *write_events = (cl_event*)malloc(num_buffers*sizeof(cl_event));
  int buf_idx = 0;
  ret = clEnqueueWriteBuffer(write_queue, ssha_device, 1, 0,
			     (size_t)buff_size, (void *)ssha, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(write_queue, sshn_device, 1, 0,
			     (size_t)buff_size, (void *)sshn, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);

  ret = clWaitForEvents(num_buffers, write_events);
  check_status("clWaitForEvents", ret);

  /*------------------------------------------------------------*/
  /* Run the kernels */

  /* Write to channel */
  ret = clEnqueueTask(write_queue, write_kernel, 0, NULL, NULL);
  check_status("clEnqueueTask(write_kernel)", ret);

  /* Read from channel */
  ret = clEnqueueTask(read_queue, read_kernel, 0, NULL, NULL);
  check_status("clEnqueueTask(read_kernel)", ret);
  
  fprintf(stdout, "Waiting for kernels to finish...\n");
  
  ret = clFinish(write_queue);
  check_status("clFinish", ret);
  ret = clFinish(read_queue);
  check_status("clFinish", ret);

  fprintf(stdout, "Kernels done!\n");

  /* Copy data back from device, synchronously. */
  cl_event read_events[3];
  ret = clEnqueueReadBuffer(read_queue, sshn_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)sshn, 0, NULL,
			    &(read_events[0]));
  check_status("clEnqueueReadBuffer", ret);
  clWaitForEvents(1, read_events);
  check_status("clWaitForEvents", ret);

  fprintf(stdout, "First value in read buffer = %d\n", sshn[0]);
  fprintf(stdout, "Second value in read buffer = %d\n", sshn[1]);
  
  /* Clean up */
  ret = clFlush(read_queue);
  ret = clFinish(read_queue);

  ret = clFlush(write_queue);
  ret = clFinish(write_queue);
  
  ret = clReleaseKernel(write_kernel);
  ret = clReleaseKernel(read_kernel);

  ret = clReleaseProgram(program);
  ret = clReleaseMemObject(ssha_device);
  ret = clReleaseMemObject(sshn_device);

  ret = clReleaseCommandQueue(write_queue);
  ret = clReleaseCommandQueue(read_queue);

  ret = clReleaseContext(context);

  if(image_file){
    free(image_file);
  }
}
