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

/** The number of concurrent kernels in our design */
#define NUM_KERNELS 2

/** Top-level driver program. Queries the hardware to find OpenCL devices,
    creates OpenCL kernels and runs them. Also runs the same kernels on
    the CPU and compares the results. */
int main(){
  /** Pointer to the OpenCL device (hardware) that we will use */
  cl_device_id device;
  cl_context context = NULL;
  /* String holding info on chosen device */
  char version_str[128];

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

  cl_int ret;

  /* Our OpenCL kernel objects */
  cl_kernel read_kernel = NULL;
  cl_kernel write_kernel  = NULL;

  /* Ptr used when getting env variables */
  char *env_string;
  /** Name of file containing single, compiled image for multiple
      kernels */
  char *image_file = NULL;

  /** Default problem size. May be overridden by setting
      NEMOLITE2D_N{X,Y} environment variables. */
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

  init_device(&device, version_str, &context);
  
  /* Create Command Queue with properties set to NULL */
  for(ji=0; ji<NUM_KERNELS; ji++){
    /* The Intel/Altera OpenCL SDK is only version 1.0 */
    /* NVIDIA only support OpenCL 1.2 so we get a seg. fault if we attempt
       to call the ...WithProperties version of this routine */
    command_queue[ji] = clCreateCommandQueue(context,
					     device,
					     queue_properties, &ret);
    check_status("clCreateCommandQueue", ret);
  }

  /* Create OpenCL Kernels and associated event objects (latter used
   to obtain detailed timing information). */
  cl_event cont_evt;
  if(image_file){
    write_kernel = get_kernel(&context, &device, version_str,
			     image_file,
			     "channel_write");
  }
  else{
    fprintf(stderr, "Please set NEMOLITE2D_SINGLE_IMAGE to point to the "
	    ".aocx file containing the compiled kernels\n");
    exit(1);
  }
  cl_event next_sshu_evt;
  if(image_file){
    read_kernel = get_kernel(&context, &device, version_str,
			     image_file,
			     "channel_read");
  }
  else{
    fprintf(stderr, "Please set NEMOLITE2D_SINGLE_IMAGE to point to the "
	    ".aocx file containing the compiled kernels\n");
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

  /* Set a recognisable first value in the input buffer */
  ssha[0] = 99.0;
  
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

  /*------------------------------------------------------------*/
  /* Run the kernels */

  for(istep=1; istep<=nsteps; istep++){

    /* Write to channel */
    ret = clEnqueueTask(command_queue[1], write_kernel, 0,
			NULL, NULL);
    check_status("clEnqueueTask(write_kernel)", ret);

    /* Read from channel */
    ret = clEnqueueTask(command_queue[0], read_kernel, 0,
			NULL, NULL);
    check_status("clEnqueueTask(read_kernel)", ret);
  }
  
  fprintf(stdout, "Waiting for kernels to finish...\n");
  
  for(int i=0; i<NUM_KERNELS; ++i) {
    ret = clFinish(command_queue[i]);
    check_status("clFinish", ret);
  }

  fprintf(stdout, "Kernels done!\n");

  /* Copy data back from device, synchronously. */
  cl_event read_events[3];
  int nread = 0;
  ret = clEnqueueReadBuffer(command_queue[0], sshn_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)sshn, 0, NULL,
			    &(read_events[0]));
  nread++;
  check_status("clEnqueueReadBuffer", ret);
  clWaitForEvents(nread, read_events);
  check_status("clWaitForEvents", ret);

  fprintf(stdout, "First value in read buffer = %lf\n", sshn[0]);
  
  /* Dump final fields computed on OpenCL device */
  write_field("sshn_ocl.dat", nx, ny, 0, 0, sshn);

  /* Extract profiling info from the OpenCL runtime. Note that this will
     be the time taken during the most recent execution of each kernel. */
  if(profiling_enabled){
    fprintf(stdout, "Time spent in Continuity kern = %e us\n",
  	    duration_ns(cont_evt)*0.001);
    fprintf(stdout, "Time spent in next-sshu kern = %e us\n",
            duration_ns(next_sshu_evt)*0.001);
  }

  /* Clean up */
  for(ji=0; ji<NUM_KERNELS; ji++){
    ret = clFlush(command_queue[ji]);
    ret = clFinish(command_queue[ji]);
  }
  ret = clReleaseKernel(write_kernel);
  ret = clReleaseKernel(read_kernel);

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
