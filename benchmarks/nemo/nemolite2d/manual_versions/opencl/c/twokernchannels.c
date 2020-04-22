#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Headers for the C versions of the kernels */
#include "continuity.h"
#include "time_update.h"

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
  cl_device_id device = NULL;
  cl_context context = NULL;
  /* String holding info on chosen device */
  char version_str[STD_STRING_LEN];

  /* In order to run multiple kernels concurrently, we must have
     multiple command queues (since the Intel OpenCL SDK only supports
     in-order queues) */
  cl_command_queue command_queue[NUM_KERNELS];
  cl_command_queue_properties queue_properties = 0;
  cl_bool profiling_enabled = 0;

  /* Buffers on the device */
  cl_mem ssha_device = NULL;
  cl_mem ssha_u_device = NULL;
  cl_mem ssha_v_device = NULL;
  cl_mem sshn_device = NULL;
  cl_mem sshn_u_device = NULL;
  cl_mem sshn_v_device = NULL;
  cl_mem ht_device = NULL;
  cl_mem hu_device = NULL;
  cl_mem hv_device = NULL;
  cl_mem un_device = NULL;
  cl_mem vn_device = NULL;
  cl_mem ua_device = NULL;
  cl_mem va_device = NULL;
  cl_mem e1u_device = NULL;
  cl_mem e1v_device = NULL;
  cl_mem e1t_device = NULL;
  cl_mem e2u_device = NULL;
  cl_mem e2v_device = NULL;
  cl_mem e2t_device = NULL;
  cl_mem e12u_device = NULL;
  cl_mem e12v_device = NULL;
  cl_mem e12t_device = NULL;
  cl_mem gphiu_device = NULL;
  cl_mem gphiv_device = NULL;
  cl_mem tmask_device = NULL;

  cl_program program = NULL;

  /* Our OpenCL kernel objects */
  cl_kernel cont_kernel = NULL;
  cl_kernel next_sshu_kernel  = NULL;

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
  
  init_device(&device, version_str, &context);
  
  /* Create Command Queue with properties set to NULL */
  for(ji=0; ji<NUM_KERNELS; ji++){
    /* The Intel/Altera OpenCL SDK is only version 1.0 */
    /* NVIDIA only support OpenCL 1.2 so we get a seg. fault if we attempt
       to call the ...WithProperties version of this routine */
    command_queue[ji] = clCreateCommandQueue(context, device,
					     queue_properties, &ret);
    check_status("clCreateCommandQueue", ret);
  }

  /* Create OpenCL Kernels and associated event objects (latter used
   to obtain detailed timing information). */
  cl_event cont_evt;
  cl_event next_sshu_evt;
  if(image_file){
    /* We expect all kernels to have been compiled into a single
       image */
    program = get_program(context, &device, version_str, image_file);
  }
  else{
    fprintf(stderr, "Please set NEMOLITE2D_SINGLE_IMAGE to point to the "
	    ".aocx file containing the compiled kernels\n");
    exit(1);
  }

  /* Create the kernels associated with the program object */
  cont_kernel = clCreateKernel(program, "continuity_code", &ret);
  check_status("clCreateKernel", ret);
  next_sshu_kernel = clCreateKernel(program, "next_sshu_code", &ret);
  check_status("clCreateKernel", ret);

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
  sshn_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
				 NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  sshn_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
				 NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  hu_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  hv_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  
  /* Velocity fields */
  un_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  vn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);

  /* Mesh scale factors */
  e12u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                          NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  e12t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                          NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  tmask_device = clCreateBuffer(context, CL_MEM_READ_ONLY,
				(size_t)(nx*ny*sizeof(cl_int)), NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);

  fprintf(stdout, "Created %d device buffers OK\n", num_buffers);

  /* Set OpenCL Kernel Parameters for Continuity */
  set_args_continuity(cont_kernel, &nx,
		      &ssha_device, &sshn_device,
		      &sshn_u_device, &sshn_v_device,
		      &hu_device, &hv_device,
		      &un_device, &vn_device,
		      &rdt, &e12t_device); 
  
  /* Set OpenCL Kernel Parameters for next_sshu kernel */
  set_args_next_sshu(next_sshu_kernel,
		     &nx, &sshn_u_device, &sshn_device, &tmask_device,
		     &e12t_device, &e12u_device);
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
  ret = clEnqueueWriteBuffer(command_queue[0], sshn_u_device, 1, 0,
			     (size_t)buff_size, (void *)sshn_u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], sshn_v_device, 1, 0,
			     (size_t)buff_size, (void *)sshn_v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], hu_device, 1, 0,
			     (size_t)buff_size, (void *)hu, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], hv_device, 1, 0,
			     (size_t)buff_size, (void *)hv, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], un_device, 1, 0,
			     (size_t)buff_size, (void *)un, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], vn_device, 1, 0,
			     (size_t)buff_size, (void *)vn, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], e12u_device, 1, 0,
			     (size_t)buff_size, (void *)e12u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], e12t_device, 1, 0,
			     (size_t)buff_size, (void *)e12t, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], tmask_device, 1, 0,
			     (size_t)(nx*ny*sizeof(cl_int)), (void *)tmask, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);

  ret = clWaitForEvents(num_buffers, write_events);
  check_status("clWaitForEvents", ret);

  TimerStop();
  
  /*------------------------------------------------------------*/
  /* Run the kernels */

  // Thread block size 
  size_t global_size[2] = {(size_t)nx, (size_t)ny};
  size_t local_size[2] = {64, 1};

  TimerStart("Time-stepping, OpenCL");
  
  for(istep=1; istep<=nsteps; istep++){

    ret = clEnqueueNDRangeKernel(command_queue[0], cont_kernel, 2, 0,
    				 global_size, NULL, 0, NULL, &cont_evt);
    check_status("clEnqueueNDRangeKernel(Continuity)", ret);

    /* Update of sshu field */
    ret = clEnqueueNDRangeKernel(command_queue[1], next_sshu_kernel, 2, 0,
				 global_size, NULL, 0, 0, &next_sshu_evt);
    check_status("clEnqueueNDRangeKernel(next-sshu)", ret);
  }

  /* Block on the execution of the last kernel */
  ret = clWaitForEvents(1, &next_sshu_evt);
  check_status("clWaitForEvents", ret);

  TimerStop();
  
  /* Run the kernels on the CPU */

  TimerStart("Time-stepping, C");
  
  for(istep=1; istep<=nsteps; istep++){
    
    for(jj=ystart; jj<ystop; jj++){
      for(ji=xstart; ji<xstop; ji++){
    	continuity_code(ji, jj, nx,
    			ssha, sshn, sshn_u, sshn_v, hu, hv,
    			un, vn, rdt, e12t);
      }
    }

    for(jj=ystart; jj<ystop; jj++){
      for(ji=xstart; ji<xstop-1; ji++){
	next_sshu_code(ji, jj, nx,
		       sshn_u, sshn, tmask, e12t, e12u);
      }
    }
  }

  TimerStop();

  write_field("ssha_cpu.dat", nx, ny, 0, 0, ssha);
  write_field("ua_cpu.dat", nx, ny, 0, 0, ua);
  write_field("va_cpu.dat", nx, ny, 0, 0, va);
  
  /* Dump final fields computed on CPU */
  cpu_sum[0] = checksum(ssha, nx, xstop, ystop, xstart, ystart);
  cpu_sum[1] = checksum(ua, nx, xstop-1, ystop, xstart, ystart);
  cpu_sum[2] = checksum(va, nx, xstop, ystop-1, xstart, ystart);

  /* Copy data back from device, synchronously. */
  cl_event read_events[3];
  int nread = 0;
  ret = clEnqueueReadBuffer(command_queue[0], sshn_u_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)sshn_u, 0, NULL,
			    &(read_events[0]));
  nread++;
  check_status("clEnqueueReadBuffer", ret);
  clWaitForEvents(nread, read_events);
  check_status("clWaitForEvents", ret);

  /* Compute and output a checksum */
  ocl_sum[0] = checksum(ssha, nx, xstop, ystop, xstart, ystart);

  /* Dump final fields computed on OpenCL device */
  write_field("sshn_u_ocl.dat", nx, ny, 0, 0, sshn_u);
  write_field("ua_ocl.dat", nx, ny, 0, 0, ua);
  write_field("va_ocl.dat", nx, ny, 0, 0, va);

  fprintf(stdout, "After %d time-steps for %d x %d domain:\n",
	  nsteps, nx-1, ny-1);
  fprintf(stdout, "ssha checksums (CPU, OpenCL) = %e, %e. Diff = %e\n",
          cpu_sum[0], ocl_sum[0], cpu_sum[0]-ocl_sum[0]);
  ocl_sum[1] = checksum(ua, nx, xstop-1, ystop, xstart, ystart);
  fprintf(stdout, "ua checksums (CPU, OpenCL) = %e, %e. Diff = %e\n",
          cpu_sum[1], ocl_sum[1], cpu_sum[1]-ocl_sum[1]);
  ocl_sum[2] = checksum(va, nx, xstop, ystop-1, xstart, ystart);
  fprintf(stdout, "va checksums (CPU, OpenCL) = %e, %e. Diff = %e\n",
          cpu_sum[2], ocl_sum[2], cpu_sum[2]-ocl_sum[2]);

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
  ret = clReleaseMemObject(sshn_u_device);
  ret = clReleaseMemObject(sshn_v_device);
  ret = clReleaseMemObject(hu_device);
  ret = clReleaseMemObject(hv_device);
  ret = clReleaseMemObject(un_device);
  ret = clReleaseMemObject(vn_device);
  ret = clReleaseMemObject(e12t_device);
  ret = clReleaseMemObject(e12u_device);
  ret = clReleaseMemObject(tmask_device);

  for(ji=0; ji<NUM_KERNELS; ji++){
    ret = clReleaseCommandQueue(command_queue[ji]);
  }
  ret = clReleaseContext(context);

  if(image_file){
    free(image_file);
  }
}
