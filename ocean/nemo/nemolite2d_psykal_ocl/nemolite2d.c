#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Headers for the C versions of the kernels */
#include "continuity.h"
#include "momentum.h"
#include "boundary_conditions.h"
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

/* Number of OpenCL command queues we will use */
#define NUM_QUEUES 2

enum KERNELS {
  K_CONTINUITY,
  K_MOM_U,
  K_MOM_V,
  K_BC_SSH,
  K_BC_SOLID_U,
  K_BC_SOLID_V,
  K_BC_FLATHER_U,
  K_BC_FLATHER_V,
  K_NEXT_SSH_U,
  K_NEXT_SSH_V,
  K_NUM_KERNELS
};

/* Names of all of our kernels */
static const char* kernel_names[K_NUM_KERNELS] =
{
  "continuity_code",
  "momentum_u_code",
  "momentum_v_code",
  "bc_ssh_code",
  "bc_solid_u_code",
  "bc_solid_v_code",
  "bc_flather_u_code",
  "bc_flather_v_code",
  "next_sshu_code",
  "next_sshv_code"
};

/* The source files containing each kernel (only used if not compiling
   them to a single image) */
static const char* kernel_files[K_NUM_KERNELS] =
{
  "./continuity_kern.c",
  "./momentum_u_kern.c",
  "./momentum_v_kern.c",
  "./boundary_conditions_kern.c",
  "./boundary_conditions_kern.c",
  "./boundary_conditions_kern.c",
  "./boundary_conditions_kern.c",
  "./boundary_conditions_kern.c",
  "./time_update_kern.c",
  "./time_update_kern.c"   
};

/** Top-level driver program. Queries the hardware to find OpenCL devices,
    creates OpenCL kernels and runs them. Also runs the same kernels on
    the CPU and compares the results. */
int main(){
  cl_device_id device;
  cl_context context = NULL;
  /* String holding info on chosen device */
  char version_str[STD_STRING_LEN];

  /* In order to run multiple kernels concurrently, we must have
     multiple command queues (since the Intel OpenCL SDK only supports
     in-order queues) */
  cl_command_queue command_queue[NUM_QUEUES];
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
  cl_kernel clkernel[K_NUM_KERNELS];
  cl_event clkernevt[K_NUM_KERNELS];

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
  int ji, jj, ikern;
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
    /* the %ms below instructs sscanf to allocate memory as required */
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
  
  /* Create Command Queues */
  /* The Intel/Altera OpenCL SDK is only version 1.0 */
  /* NVIDIA only support OpenCL 1.2 so we get a seg. fault if we attempt
     to call the ...WithProperties version of this routine */
  for(ji=0; ji<NUM_QUEUES; ji++){
    /* The Intel/Altera OpenCL SDK is only version 1.0 */
    /* NVIDIA only support OpenCL 1.2 so we get a seg. fault if we attempt
       to call the ...WithProperties version of this routine */
    command_queue[ji] = clCreateCommandQueue(context, device,
					     queue_properties, &ret);
    check_status("clCreateCommandQueue", ret);
  }

  // Create our Program object (contains all of the individual kernels)
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
    
  /* Create OpenCL Kernels and associated event objects (latter used
   to obtain detailed timing information). */
  for(ikern=0; ikern<K_NUM_KERNELS; ikern++){

    if(!image_file){
      program = get_program(context, &device, version_str,
			    kernel_files[ikern]);
    }
    fprintf(stdout, "Creating kernel %s...\n", kernel_names[ikern]);
    clkernel[ikern] = clCreateKernel(program, kernel_names[ikern], &ret);
    check_status("clCreateKernel", ret);
  }

  /* Create Device Memory Buffers */
  int num_buffers = 0;
  buff_size = nx*ny*sizeof(cl_double);
  ssha_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			       NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  ssha_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
				 NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  ssha_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
				 NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
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
  check_status("clCreateBuffer", ret);
  ht_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  
  /* Velocity fields */
  un_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  vn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
#ifndef SINGLE_KERNEL
  ua_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  va_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
#endif

  /* Mesh scale factors */
  e1t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			      NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  e1u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			      NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  e1v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			      NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  e2u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			      NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  e2v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			      NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  e2t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			      NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  e12u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                          NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);

  e12v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                          NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);

  tmask_device = clCreateBuffer(context, CL_MEM_READ_ONLY,
				(size_t)(nx*ny*sizeof(cl_int)), NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  e12t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                          NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);

  /* Coriolis parameters */
  gphiu_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
				NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  gphiv_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
				NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);

  fprintf(stdout, "Created %d device buffers OK\n", num_buffers);

  /* Set OpenCL Kernel Parameters for Continuity */
  set_args_continuity(clkernel[K_CONTINUITY],
		      &nx,
		      &ssha_device, &sshn_device,
		      &sshn_u_device, &sshn_v_device,
		      &hu_device, &hv_device,
		      &un_device, &vn_device,
		      &rdt, &e12t_device); 
  /* Set OpenCL Kernel Parameters for Momentum-u */
  set_args_momu(clkernel[K_MOM_U],
		&nx,
		&ua_device, &un_device, &vn_device,
		&hu_device, &hv_device, &ht_device,
		&ssha_u_device, &sshn_device,
		&sshn_u_device, &sshn_v_device,
		&tmask_device,
		&e1u_device, &e1v_device,
		&e1t_device, &e2u_device,
		&e2t_device, &e12u_device,
		&gphiu_device,
		&rdt, &cbfr, &visc);
  /* Set OpenCL Kernel Parameters for Momentum-v */
  set_args_momv(clkernel[K_MOM_V],
		&nx,
		&va_device, &un_device, &vn_device,
		&hu_device, &hv_device, &ht_device,
		&ssha_v_device, &sshn_device,
		&sshn_u_device, &sshn_v_device,
		&tmask_device,
		&e1v_device, &e1t_device,
		&e2u_device, &e2v_device,
		&e2t_device, &e12v_device,
		&gphiv_device,
		&rdt, &cbfr, &visc);

  /* Set kernel arguments for bc_ssh */
  int arg_idx = 0;
  ret = clSetKernelArg(clkernel[K_BC_SSH], arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  arg_idx++; // Skip time-step argument here - do in time-stepping loop
  /* ret = clSetKernelArg(clkernel[K_BC_SSH], arg_idx++, sizeof(cl_int), */
  /* 		       (void *)&istep); */
  /* check_status("clSetKernelArg", ret); */
  ret = clSetKernelArg(clkernel[K_BC_SSH], arg_idx++, sizeof(cl_mem),
		       (void *)&ssha_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_SSH], arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_SSH], arg_idx++, sizeof(cl_double),
		       (void *)&rdt);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc-ssh kernel\n", arg_idx);

  /* Set OpenCL Kernel Parameters for bc_solid_u kernel */
  arg_idx = 0;
  ret = clSetKernelArg(clkernel[K_BC_SOLID_U], arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_SOLID_U], arg_idx++, sizeof(cl_mem),
		       (void *)&ua_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_SOLID_U], arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc_solid_u kernel\n", arg_idx);

  /* Set OpenCL Kernel Parameters for bc_solid_v kernel */
  arg_idx = 0;
  ret = clSetKernelArg(clkernel[K_BC_SOLID_V], arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_SOLID_V], arg_idx++, sizeof(cl_mem),
		       (void *)&va_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_SOLID_V], arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc_solid_v kernel\n", arg_idx);

  /* Set OpenCL Kernel Parameters for bc_flather_u kernel */
  arg_idx = 0;
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_U], arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_U], arg_idx++, sizeof(cl_mem),
		       (void *)&ua_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_U], arg_idx++, sizeof(cl_mem),
		       (void *)&hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_U], arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_U], arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc_flather_u kernel\n", arg_idx);

  /* Set OpenCL Kernel Parameters for bc_flather_v kernel */
  arg_idx = 0;
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_V], arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_V], arg_idx++, sizeof(cl_mem),
		       (void *)&va_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_V], arg_idx++, sizeof(cl_mem),
		       (void *)&hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_V], arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(clkernel[K_BC_FLATHER_V], arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc_flather_v kernel\n", arg_idx);
  
  /* Set OpenCL Kernel Parameters for next_sshu kernel */
  set_args_next_sshu(clkernel[K_NEXT_SSH_U],
		     &nx, &sshn_u_device, &sshn_device, &tmask_device,
		     &e12t_device, &e12u_device);
  
  /* Set OpenCL Kernel Parameters for next_sshv kernel */
  set_args_next_sshv(clkernel[K_NEXT_SSH_V],
		     &nx, &sshn_v_device, &sshn_device, &tmask_device,
		     &e12t_device, &e12v_device);
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
#ifndef SINGLE_KERNEL
  ret = clEnqueueWriteBuffer(command_queue[0], ssha_u_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], ssha_v_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
#endif
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
#ifndef SINGLE_KERNEL
  ret = clEnqueueWriteBuffer(command_queue[0], ht_device, 1, 0,
			     (size_t)buff_size, (void *)ht, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
#endif
  ret = clEnqueueWriteBuffer(command_queue[0], un_device, 1, 0,
			     (size_t)buff_size, (void *)un, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], vn_device, 1, 0,
			     (size_t)buff_size, (void *)vn, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
#ifndef SINGLE_KERNEL
  ret = clEnqueueWriteBuffer(command_queue[0], ua_device, 1, 0,
			     (size_t)buff_size, (void *)ua, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], va_device, 1, 0,
			     (size_t)buff_size, (void *)va, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], e1u_device, 1, 0,
			     (size_t)buff_size, (void *)e1u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], e1v_device, 1, 0,
			     (size_t)buff_size, (void *)e1v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], e1t_device, 1, 0,
			     (size_t)buff_size, (void *)e1t, 0,
			     NULL, &(write_events[buf_idx++]));
  ret = clEnqueueWriteBuffer(command_queue[0], e2u_device, 1, 0,
			     (size_t)buff_size, (void *)e2u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], e2v_device, 1, 0,
			     (size_t)buff_size, (void *)e2v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], e2t_device, 1, 0,
			     (size_t)buff_size, (void *)e2t, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], e12u_device, 1, 0,
			     (size_t)buff_size, (void *)e12u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], e12v_device, 1, 0,
			     (size_t)buff_size, (void *)e12v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
#endif
  ret = clEnqueueWriteBuffer(command_queue[0], e12t_device, 1, 0,
			     (size_t)buff_size, (void *)e12t, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
#ifndef SINGLE_KERNEL
  ret = clEnqueueWriteBuffer(command_queue[0], gphiu_device, 1, 0,
			     (size_t)buff_size, (void *)gphiu, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], gphiv_device, 1, 0,
			     (size_t)buff_size, (void *)gphiv, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue[0], tmask_device, 1, 0,
			     (size_t)(nx*ny*sizeof(cl_int)), (void *)tmask, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
#endif
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

#ifndef SINGLE_KERNEL
    /* Update the time-step index argument to the bc-ssh kernel */
    ret = clSetKernelArg(clkernel[K_BC_SSH], 1, sizeof(cl_int),
			 (void *)&istep);
    check_status("clSetKernelArg", ret);
#endif

    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_CONTINUITY], 2, 0,
    				 global_size, NULL, 0, NULL,
				 &(clkernevt[K_CONTINUITY]));
    check_status("clEnqueueNDRangeKernel(Continuity)", ret);

#ifndef SINGLE_KERNEL
    /* Experimentation on a K20c showed that this work-group size was
       optimal for that device */
    local_size[0] = 64;
    local_size[1] = 1;
    if (local_size[0] > nx) local_size[0] = nx;
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_MOM_U], 2, 0,
				 global_size, local_size, 0, NULL,
				 &(clkernevt[K_MOM_U]));
    check_status("clEnqueueNDRangeKernel(Mom-u)", ret);
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_MOM_V], 2, 0,
				 global_size, local_size, 0, NULL,
				 &(clkernevt[K_MOM_V]));
    check_status("clEnqueueNDRangeKernel(Mom-v)", ret);
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_SSH], 2, NULL,
    				 global_size, NULL, 0,0,
				 &(clkernevt[K_BC_SSH]));
    check_status("clEnqueueNDRangeKernel(bc-ssh)", ret);

    /* Apply boundary conditions */
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_SOLID_U], 2, 0,
				 global_size, NULL, 0, 0,
				 &(clkernevt[K_BC_SOLID_U]));
    check_status("clEnqueueNDRangeKernel(bc-solid-u)", ret);
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_SOLID_V], 2, 0,
    				 global_size, NULL, 0, 0,
				 &(clkernevt[K_BC_SOLID_V]));
    check_status("clEnqueueNDRangeKernel(bc-solid-v)", ret);
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_FLATHER_U], 2, 0,
				 global_size, NULL, 0, 0,
				 &(clkernevt[K_BC_FLATHER_U]));
    check_status("clEnqueueNDRangeKernel(bc-flather-u)", ret);
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_FLATHER_V], 2, 0,
    				 global_size, NULL, 0, 0,
				 &(clkernevt[K_BC_FLATHER_V]));
    check_status("clEnqueueNDRangeKernel(bc-flather-v)", ret);

    /* Copy 'after' fields to 'now' fields */
    ret = clEnqueueCopyBuffer(command_queue[0], ua_device, un_device, 0, 0,
			      buff_size,0, NULL, NULL);
    check_status("clEnqueueCopyBuffer", ret);
    ret = clEnqueueCopyBuffer(command_queue[0], va_device, vn_device, 0, 0,
			      buff_size,0, NULL, NULL);
    check_status("clEnqueueCopyBuffer", ret);
    ret = clEnqueueCopyBuffer(command_queue[0], ssha_device, sshn_device, 0, 0,
			      buff_size,0, NULL, NULL);
    check_status("clEnqueueCopyBuffer", ret);

    /* Update of sshu and sshv fields */
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_NEXT_SSH_U], 2, 0,
				 global_size, NULL, 0, 0,
				 &(clkernevt[K_NEXT_SSH_U]));
    check_status("clEnqueueNDRangeKernel(next-sshu)", ret);
  
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_NEXT_SSH_U], 2, 0,
				 global_size, NULL, 0, 0,
				 &(clkernevt[K_NEXT_SSH_V]));
    check_status("clEnqueueNDRangeKernel(next-sshv)", ret);
#endif
  }

  /* Block on the execution of the last kernel */
  ret = clWaitForEvents(1, &(clkernevt[K_NEXT_SSH_V]));
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
#ifndef SINGLE_KERNEL
    for(jj=ystart; jj<ystop; jj++){
      for(ji=xstart; ji<xstop-1; ji++){
	momentum_u_code(ji, jj, nx,                     
			ua, un, vn, hu, hv, ht,
			ssha_u, sshn, sshn_u, sshn_v, tmask,
			e1u, e1v, e1t, e2u, e2t, e12u, gphiu,
			rdt, cbfr, visc);
      }
    }
    for(jj=ystart; jj<ystop-1; jj++){
      for(ji=xstart; ji<xstop; ji++){
	momentum_v_code(ji, jj, nx,                     
			va, un, vn, hu, hv, ht,
			ssha_v, sshn, sshn_u, sshn_v, tmask,
			e1v, e1t, e2u, e2v, e2t, e12v, gphiv,
			rdt, cbfr, visc);
      }
    }
    for(jj=ystart; jj<ystop; jj++){
      for(ji=xstart; ji<xstop; ji++){
    	bc_ssh_code(ji, jj, nx,
    		    istep, ssha, tmask, rdt);
      }
    }
    /* Upper loop limit for jj should encompass whole domain here (i.e. be
       greater than any of the limits for the previous loops) */
    for(jj=ystart-1; jj<ystop+1; jj++){
      for(ji=xstart-1; ji<xstop; ji++){
	bc_solid_u_code(ji, jj, nx,
			ua, tmask);
      }
    }
    /* Upper loop limit for ji should encompass whole domain here (i.e. be
       greater than any of the limits for the previous loops) */
    for(jj=ystart-1; jj<ystop; jj++){
      for(ji=xstart-1; ji<xstop+1; ji++){
    	bc_solid_v_code(ji, jj, nx,
    			va, tmask);
      }
    }
    /* Upper loop limit for jj should encompass whole domain here (i.e. be
       greater than any of the limits for the previous loops) */
    for(jj=ystart-1; jj<ystop+1; jj++){
      for(ji=xstart-1; ji<xstop; ji++){
	bc_flather_u_code(ji, jj, nx,
			  ua, hu, sshn_u, tmask);
      }
    }
    /* Upper loop limit for ji should encompass whole domain here (i.e. be
       greater than any of the limits for the previous loops) */
    for(jj=ystart-1; jj<ystop; jj++){
      for(ji=xstart-1; ji<xstop+1; ji++){
    	bc_flather_v_code(ji, jj, nx,
    			  va, hv, sshn_v, tmask);
      }
    }
    /* Copy 'after' fields to 'now' fields */
    memcpy(un, ua, buff_size);
    memcpy(vn, va, buff_size);
    memcpy(sshn, ssha, buff_size);

    for(jj=ystart; jj<ystop; jj++){
      for(ji=xstart; ji<xstop-1; ji++){
	next_sshu_code(ji, jj, nx,
		       sshn_u, sshn, tmask, e12t, e12u);
      }
    }
    for(jj=ystart; jj<ystop-1; jj++){
      for(ji=xstart;ji<xstop; ji++){
	next_sshv_code(ji, jj, nx,
		       sshn_v, sshn, tmask, e12t, e12v);
      }
    }
#endif
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
  ret = clEnqueueReadBuffer(command_queue[0], ssha_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)ssha, 0, NULL,
			    &(read_events[0]));
  nread++;
  check_status("clEnqueueReadBuffer", ret);
#ifndef SINGLE_KERNEL
  ret = clEnqueueReadBuffer(command_queue[0], ua_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)ua, 0, NULL,
			    &(read_events[1]));
  nread++;
  check_status("clEnqueueReadBuffer", ret);
  ret = clEnqueueReadBuffer(command_queue[0], va_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)va, 0, NULL,
			    &(read_events[2]));
  nread++;
  check_status("clEnqueueReadBuffer", ret);
#endif
  clWaitForEvents(nread, read_events);
  check_status("clWaitForEvents", ret);

  /* Compute and output a checksum */
  ocl_sum[0] = checksum(ssha, nx, xstop, ystop, xstart, ystart);

  /* Dump final fields computed on OpenCL device */
  write_field("ssha_ocl.dat", nx, ny, 0, 0, ssha);
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
  	    duration_ns(clkernevt[K_CONTINUITY])*0.001);
    fprintf(stdout, "Time spent in Momentum-u kern = %e us\n",
            duration_ns(clkernevt[K_MOM_U])*0.001);
    fprintf(stdout, "Time spent in Momentum-v kern = %e us\n",
            duration_ns(clkernevt[K_MOM_V])*0.001);
  }

  /* Generate report on host timings */
  TimerReport();

  /* Clean up */
  for(ji=0; ji<NUM_QUEUES; ji++){
    ret = clFlush(command_queue[ji]);
    ret = clFinish(command_queue[ji]);
  }
  for(ikern=0; ikern<K_NUM_KERNELS; ikern++){
    ret = clReleaseKernel(clkernel[ikern]);
  }

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
  ret = clReleaseMemObject(ssha_u_device);
  ret = clReleaseMemObject(ssha_v_device);
  ret = clReleaseMemObject(ht_device);
  ret = clReleaseMemObject(ua_device);
  ret = clReleaseMemObject(va_device);
  ret = clReleaseMemObject(e1u_device);
  ret = clReleaseMemObject(e1v_device);
  ret = clReleaseMemObject(e1t_device);
  ret = clReleaseMemObject(e2u_device);
  ret = clReleaseMemObject(e2v_device);
  ret = clReleaseMemObject(e2t_device);
  ret = clReleaseMemObject(e12u_device);
  ret = clReleaseMemObject(e12v_device);
  ret = clReleaseMemObject(gphiu_device);
  ret = clReleaseMemObject(gphiv_device);
  ret = clReleaseMemObject(tmask_device);
  for(ji=0; ji<NUM_QUEUES; ji++){
    ret = clReleaseCommandQueue(command_queue[ji]);
  }
  ret = clReleaseContext(context);

  if(image_file){
    free(image_file);
  }
}
