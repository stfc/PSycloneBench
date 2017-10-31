#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Headers for the C versions of the kernels */
#include "continuity.h"
#include "momentum.h"
#include "boundary_conditions.h"
#include "time_update.h"

#include "opencl_utils.h"
#include "timing.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

/** Maximum number of OpenCL devices we will query */
#define MAX_DEVICES 4

/** Compute a checksum for a double precision array of unknown no. of rows
    but with each row containing width entries */
double checksum(double *array, int width,
		int nx, int ny, int xstart, int ystart){
  int i, j, jidx;
  double sum = 0.0;
  for(j=ystart; j<ny; j++){
    jidx = j*width;
    for(i=xstart; i<nx; i++){
      sum += array[i+jidx];
    }
  }
  return sum;
}

/** Write the supplied integer field data to the specified file. Data
    formatted for use with gnuplot's splot command. */
void write_ifield(const char *filename, int nx, int ny,
		 int xstart, int ystart, int *field){
  int ji, jj, idx;
  FILE *fp = fopen(filename, "w");
  if(!fp){
    fprintf(stderr, "write_ifield: failed to open file %s\n", filename);
    return;
  }

  idx = 0;
  for(jj=ystart; jj<ny; jj++){
    for(ji=xstart; ji<nx; ji++){
      fprintf(fp, "%d %d\n", ji, field[idx++]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

/** Write the supplied double-precision field data to the specified
    file. Data formatted for use with gnuplot's splot command. */
void write_field(const char *filename, int nx, int ny,
		 int xstart, int ystart, double *field){
  int ji, jj, idx;
  FILE *fp = fopen(filename, "w");
  if(!fp){
    fprintf(stderr, "write_field: failed to open file %s\n", filename);
    return;
  }

  idx = 0;
  for(jj=ystart; jj<ny; jj++){
    for(ji=xstart; ji<nx; ji++){
      fprintf(fp, "%d %e\n", ji, field[idx++]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
}

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
  cl_command_queue command_queue = NULL;
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
  cl_kernel momu_kernel = NULL;
  cl_kernel momv_kernel = NULL;
  cl_kernel bc_ssh_kernel = NULL;
  cl_kernel bc_solid_u_kernel  = NULL;
  cl_kernel bc_solid_v_kernel  = NULL;
  cl_kernel bc_flather_u_kernel  = NULL;
  cl_kernel bc_flather_v_kernel  = NULL;
  cl_kernel next_sshu_kernel  = NULL;
  cl_kernel next_sshv_kernel  = NULL;

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
  if(cl_version < 200){
    /* The Intel/Altera OpenCL SDK is only version 1.0 */
    /* NVIDIA only support OpenCL 1.2 so we get a seg. fault if we attempt
       to call the ...WithProperties version of this routine */
    command_queue = clCreateCommandQueue(context,
					 *device,
					 queue_properties, &ret);
    check_status("clCreateCommandQueue", ret);
  }
  else{
    //command_queue = clCreateCommandQueueWithProperties(context,
    //						       *device,
    //					       &queue_properties, &ret);
    //check_status("clCreateCommandQueueWithProperties", ret);
    fprintf(stderr, "Implement use of clCreateCommandQueueWithProperties!\n");
    exit(1);
  }
    
  /* Create OpenCL Kernels and associated event objects (latter used
   to obtain detailed timing information). */
  cl_event cont_evt;
  cont_kernel = get_kernel(&context, device, version_str,
			   "./continuity_kern.c", "continuity_code");
  cl_event momu_evt;
  momu_kernel = get_kernel(&context, device, version_str,
			   "./momentum_kern.c", "momentum_u_code");
  cl_event momv_evt;
  momv_kernel = get_kernel(&context, device, version_str,
			   "./momentum_kern.c", "momentum_v_code");
  cl_event bcssh_evt;
  bc_ssh_kernel = get_kernel(&context, device, version_str,
			     "./boundary_conditions_kern.c", "bc_ssh_code");
  cl_event solidu_evt;
  bc_solid_u_kernel = get_kernel(&context, device, version_str,
				 "./boundary_conditions_kern.c",
				 "bc_solid_u_code");
  cl_event solidv_evt;
  bc_solid_v_kernel = get_kernel(&context, device, version_str,
				 "./boundary_conditions_kern.c",
				 "bc_solid_v_code");
  cl_event flatheru_evt;
  bc_flather_u_kernel = get_kernel(&context, device, version_str,
				   "./boundary_conditions_kern.c",
				   "bc_flather_u_code");
  cl_event flatherv_evt;
  bc_flather_v_kernel = get_kernel(&context, device, version_str,
				   "./boundary_conditions_kern.c",
				   "bc_flather_v_code");
  cl_event next_sshu_evt;
  next_sshu_kernel = get_kernel(&context, device, version_str,
				"./time_update_kern.c",
				"next_sshu_code");
  cl_event next_sshv_evt;
  next_sshv_kernel = get_kernel(&context, device, version_str,
				"./time_update_kern.c",
				"next_sshv_code");

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
  ua_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  va_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  num_buffers++;
  check_status("clCreateBuffer", ret);
  
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
  e12t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
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
  int arg_idx = 0;
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_int), (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ssha_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&un_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&vn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_double),
		       (void *)&rdt);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e12t_device);
  check_status("clSetKernelArg", ret);
  fprintf(stdout, "Set %d arguments for Continuity kernel\n", arg_idx);

  /* Set OpenCL Kernel Parameters for Momentum-u */
  set_args_momu(momu_kernel,
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
  set_args_momv(momv_kernel,
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
  arg_idx = 0;
  ret = clSetKernelArg(bc_ssh_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  arg_idx++; // Skip time-step argument here - do in time-stepping loop
  /* ret = clSetKernelArg(bc_ssh_kernel, arg_idx++, sizeof(cl_int), */
  /* 		       (void *)&istep); */
  /* check_status("clSetKernelArg", ret); */
  ret = clSetKernelArg(bc_ssh_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ssha_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_ssh_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_ssh_kernel, arg_idx++, sizeof(cl_double),
		       (void *)&rdt);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc-ssh kernel\n", arg_idx);

  /* Set OpenCL Kernel Parameters for bc_solid_u kernel */
  arg_idx = 0;
  ret = clSetKernelArg(bc_solid_u_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_solid_u_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ua_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_solid_u_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc_solid_u kernel\n", arg_idx);

  /* Set OpenCL Kernel Parameters for bc_solid_v kernel */
  arg_idx = 0;
  ret = clSetKernelArg(bc_solid_v_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_solid_v_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&va_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_solid_v_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc_solid_v kernel\n", arg_idx);

  /* Set OpenCL Kernel Parameters for bc_flather_u kernel */
  arg_idx = 0;
  ret = clSetKernelArg(bc_flather_u_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_flather_u_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ua_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_flather_u_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_flather_u_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_flather_u_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc_flather_u kernel\n", arg_idx);

  /* Set OpenCL Kernel Parameters for bc_flather_v kernel */
  arg_idx = 0;
  ret = clSetKernelArg(bc_flather_v_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_flather_v_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&va_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_flather_v_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_flather_v_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_flather_v_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for bc_flather_v kernel\n", arg_idx);
  
  /* Set OpenCL Kernel Parameters for next_sshu kernel */
  set_args_next_sshv(next_sshu_kernel,
		     &nx, &sshn_u_device, &sshn_device, &tmask_device,
		     &e12t_device, &e12u_device);
  
  /* Set OpenCL Kernel Parameters for next_sshv kernel */
  set_args_next_sshv(next_sshv_kernel,
		     &nx, &sshn_v_device, &sshn_device, &tmask_device,
		     &e12t_device, &e12v_device);

  TimerStop();

  /*------------------------------------------------------------*/
  /* Field initialisation on host */
  
  ssha = (cl_double*)malloc(buff_size);
  ssha_u = (cl_double*)malloc(buff_size);
  ssha_v = (cl_double*)malloc(buff_size);
  sshn = (cl_double*)malloc(buff_size);
  sshn_u = (cl_double*)malloc(buff_size);
  sshn_v = (cl_double*)malloc(buff_size);
  hu = (cl_double*)malloc(buff_size);
  hv = (cl_double*)malloc(buff_size);
  ht = (cl_double*)malloc(buff_size);
  un = (cl_double*)malloc(buff_size);
  vn = (cl_double*)malloc(buff_size);
  ua = (cl_double*)malloc(buff_size);
  va = (cl_double*)malloc(buff_size);
  e1u = (cl_double*)malloc(buff_size);
  e1v = (cl_double*)malloc(buff_size);
  e1t = (cl_double*)malloc(buff_size);
  e2u = (cl_double*)malloc(buff_size);
  e2v = (cl_double*)malloc(buff_size);
  e2t = (cl_double*)malloc(buff_size);
  e12u = (cl_double*)malloc(buff_size);
  e12v = (cl_double*)malloc(buff_size);
  e12t = (cl_double*)malloc(buff_size);
  tmask = (cl_int*)malloc(nx*ny*sizeof(cl_int));

  gphiu = (cl_double*)malloc(buff_size);
  gphiv = (cl_double*)malloc(buff_size);

  int xstart = 1;
  int xstop = nx - 1;
  int ystart = 1;
  int ystop = ny - 1;

  for(jj=0;jj<ny;jj++){
    for(ji=0;ji<nx;ji++){
      idx = jj*nx + ji;
      hu[idx] = dep_const;
      hv[idx] = dep_const;
      ht[idx] = dep_const;
      un[idx] = 0.0;
      vn[idx] = 0.0;
      sshn_u[idx] = 0.0;
      sshn_v[idx] = 0.0;
      sshn[idx] = 0.0;
      ssha[idx] = 0.0;
      // Grid properties
      e1u[idx] = dx;
      e1v[idx] = dx;
      e1t[idx] = dx;
      e2u[idx] = dy;
      e2v[idx] = dy;
      e2t[idx] = dy;
      e12u[idx] = dx*dy;
      e12v[idx] = e12u[idx];
      e12t[idx] = e12u[idx];
      // f-plane test case (constant Coriolis parameter)
      gphiu[idx] = 50.0;
      gphiv[idx] = 50.0;
      // All inner cells
      tmask[idx] = 1;
    }
  }
  for(jj=0;jj<ny;jj++){
    idx = jj*nx;
    // West solid boundary
    for(ji=0; ji<xstart; ji++){
      tmask[idx+ji] = 0;
    }
    // East solid boundary
    for(ji=xstop; ji<nx; ji++){
      tmask[idx+ji] = 0;
    }
  }
  // Southern open boundary
  for(jj=0; jj<ystart; jj++){
    idx = jj*nx;
    for(ji=0;ji<nx;ji++){
      tmask[idx + ji] = -1;
    }
  }
  // North solid boundary
  for(jj=ystop; jj<ny; jj++){
    idx = jj*nx;
    for(ji=0;ji<nx;ji++){
      tmask[idx + ji] = 0;
    }
  }
  
  write_ifield("tmask.dat", nx, ny, 0, 0, tmask);
  
  /*------------------------------------------------------------*/
  /* Copy data to device synchronously */
  TimerStart("Write buffers to device");

  /* Create an array to store the event associated with each write
     to the device */
  cl_event *write_events = (cl_event*)malloc(num_buffers*sizeof(cl_event));
  int buf_idx = 0;
  ret = clEnqueueWriteBuffer(command_queue, ssha_device, 1, 0,
			     (size_t)buff_size, (void *)ssha, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, ssha_u_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, ssha_v_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, sshn_device, 1, 0,
			     (size_t)buff_size, (void *)sshn, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, sshn_u_device, 1, 0,
			     (size_t)buff_size, (void *)sshn_u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, sshn_v_device, 1, 0,
			     (size_t)buff_size, (void *)sshn_v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, hu_device, 1, 0,
			     (size_t)buff_size, (void *)hu, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, hv_device, 1, 0,
			     (size_t)buff_size, (void *)hv, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, ht_device, 1, 0,
			     (size_t)buff_size, (void *)ht, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, un_device, 1, 0,
			     (size_t)buff_size, (void *)un, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, vn_device, 1, 0,
			     (size_t)buff_size, (void *)vn, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, ua_device, 1, 0,
			     (size_t)buff_size, (void *)ua, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, va_device, 1, 0,
			     (size_t)buff_size, (void *)va, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e1u_device, 1, 0,
			     (size_t)buff_size, (void *)e1u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e1v_device, 1, 0,
			     (size_t)buff_size, (void *)e1v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e1t_device, 1, 0,
			     (size_t)buff_size, (void *)e1t, 0,
			     NULL, &(write_events[buf_idx++]));
  ret = clEnqueueWriteBuffer(command_queue, e2u_device, 1, 0,
			     (size_t)buff_size, (void *)e2u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e2v_device, 1, 0,
			     (size_t)buff_size, (void *)e2v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e2t_device, 1, 0,
			     (size_t)buff_size, (void *)e2t, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e12u_device, 1, 0,
			     (size_t)buff_size, (void *)e12u, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e12v_device, 1, 0,
			     (size_t)buff_size, (void *)e12v, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e12t_device, 1, 0,
			     (size_t)buff_size, (void *)e12t, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, gphiu_device, 1, 0,
			     (size_t)buff_size, (void *)gphiu, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, gphiv_device, 1, 0,
			     (size_t)buff_size, (void *)gphiv, 0,
			     NULL, &(write_events[buf_idx++]));
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, tmask_device, 1, 0,
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

    /* Update the time-step index argument to the bc-ssh kernel */
    ret = clSetKernelArg(bc_ssh_kernel, 1, sizeof(cl_int),
			 (void *)&istep);
    check_status("clSetKernelArg", ret);

    ret = clEnqueueNDRangeKernel(command_queue, cont_kernel, 2, 0,
    				 global_size, NULL, 0, NULL, &cont_evt);
    check_status("clEnqueueNDRangeKernel(Continuity)", ret);

    /* Experimentation on a K20c showed that this work-group size was
       optimal for that device */
    local_size[0] = 64;
    local_size[1] = 1;
    if (local_size[0] > nx) local_size[0] = nx;
    ret = clEnqueueNDRangeKernel(command_queue, momu_kernel, 2, 0,
				 global_size, local_size, 0, NULL, &momu_evt);
    check_status("clEnqueueNDRangeKernel(Mom-u)", ret);
    ret = clEnqueueNDRangeKernel(command_queue, momv_kernel, 2, 0,
				 global_size, local_size, 0, NULL, &momv_evt);
    check_status("clEnqueueNDRangeKernel(Mom-v)", ret);
    ret = clEnqueueNDRangeKernel(command_queue, bc_ssh_kernel, 2, NULL,
    				 global_size, NULL, 0,0, &bcssh_evt);
    check_status("clEnqueueNDRangeKernel(bc-ssh)", ret);

    /* Apply boundary conditions */
    ret = clEnqueueNDRangeKernel(command_queue, bc_solid_u_kernel, 2, 0,
				 global_size, NULL, 0,0, &solidu_evt);
    check_status("clEnqueueNDRangeKernel(bc-solid-u)", ret);
    ret = clEnqueueNDRangeKernel(command_queue, bc_solid_v_kernel, 2, 0,
    				 global_size, NULL, 0,0, &solidv_evt);
    check_status("clEnqueueNDRangeKernel(bc-solid-v)", ret);
    ret = clEnqueueNDRangeKernel(command_queue, bc_flather_u_kernel, 2, 0,
				 global_size, NULL, 0,0, &flatheru_evt);
    check_status("clEnqueueNDRangeKernel(bc-flather-u)", ret);
    ret = clEnqueueNDRangeKernel(command_queue, bc_flather_v_kernel, 2, 0,
    				 global_size, NULL, 0,0, &flatherv_evt);
    check_status("clEnqueueNDRangeKernel(bc-flather-v)", ret);

    /* Copy 'after' fields to 'now' fields */
    ret = clEnqueueCopyBuffer(command_queue, ua_device, un_device, 0, 0,
			      buff_size,0, NULL, NULL);
    check_status("clEnqueueCopyBuffer", ret);
    ret = clEnqueueCopyBuffer(command_queue, va_device, vn_device, 0, 0,
			      buff_size,0, NULL, NULL);
    check_status("clEnqueueCopyBuffer", ret);
    ret = clEnqueueCopyBuffer(command_queue, ssha_device, sshn_device, 0, 0,
			      buff_size,0, NULL, NULL);
    check_status("clEnqueueCopyBuffer", ret);

    /* Update of sshu and sshv fields */
    ret = clEnqueueNDRangeKernel(command_queue, next_sshu_kernel, 2, 0,
				 global_size, NULL, 0, 0, &next_sshu_evt);
    check_status("clEnqueueNDRangeKernel(next-sshu)", ret);
  
    ret = clEnqueueNDRangeKernel(command_queue, next_sshv_kernel, 2, 0,
				 global_size, NULL, 0, 0, &next_sshv_evt);
    check_status("clEnqueueNDRangeKernel(next-sshv)", ret);
  }

  /* Block on the execution of the last kernel */
  ret = clWaitForEvents(1, &next_sshv_evt);
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
  clEnqueueReadBuffer(command_queue, ssha_device, CL_TRUE, 0,
		      (size_t)buff_size, (void *)ssha, 0, NULL,
		      &(read_events[0]));
  check_status("clEnqueueReadBuffer", ret);
  clEnqueueReadBuffer(command_queue, ua_device, CL_TRUE, 0,
		      (size_t)buff_size, (void *)ua, 0, NULL,
		      &(read_events[1]));
  check_status("clEnqueueReadBuffer", ret);
  clEnqueueReadBuffer(command_queue, va_device, CL_TRUE, 0,
		      (size_t)buff_size, (void *)va, 0, NULL,
		      &(read_events[2]));
  check_status("clEnqueueReadBuffer", ret);

  clWaitForEvents(3, read_events);
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
  	    duration_ns(cont_evt)*0.001);
    fprintf(stdout, "Time spent in Momentum-u kern = %e us\n",
            duration_ns(momu_evt)*0.001);
    fprintf(stdout, "Time spent in Momentum-v kern = %e us\n",
            duration_ns(momv_evt)*0.001);
  }

  /* Generate report on host timings */
  TimerReport();

  /* Clean up */
  ret = clFlush(command_queue);
  ret = clFinish(command_queue);
  ret = clReleaseKernel(cont_kernel);
  ret = clReleaseKernel(momu_kernel);
  ret = clReleaseKernel(momv_kernel);
  ret = clReleaseKernel(bc_ssh_kernel);
  ret = clReleaseKernel(bc_solid_u_kernel);
  ret = clReleaseKernel(bc_solid_v_kernel);
  ret = clReleaseKernel(bc_flather_u_kernel);
  ret = clReleaseKernel(bc_flather_v_kernel);
  ret = clReleaseKernel(next_sshu_kernel);
  ret = clReleaseKernel(next_sshv_kernel);

  ret = clReleaseProgram(program);
  ret = clReleaseMemObject(ssha_u_device);
  ret = clReleaseMemObject(ssha_v_device);
  ret = clReleaseMemObject(ssha_device);
  ret = clReleaseMemObject(sshn_device);
  ret = clReleaseMemObject(sshn_u_device);
  ret = clReleaseMemObject(sshn_v_device);
  ret = clReleaseMemObject(ht_device);
  ret = clReleaseMemObject(hu_device);
  ret = clReleaseMemObject(hv_device);
  ret = clReleaseMemObject(ua_device);
  ret = clReleaseMemObject(va_device);
  ret = clReleaseMemObject(un_device);
  ret = clReleaseMemObject(vn_device);
  ret = clReleaseMemObject(e1u_device);
  ret = clReleaseMemObject(e1v_device);
  ret = clReleaseMemObject(e1t_device);
  ret = clReleaseMemObject(e2u_device);
  ret = clReleaseMemObject(e2v_device);
  ret = clReleaseMemObject(e2t_device);
  ret = clReleaseMemObject(e12u_device);
  ret = clReleaseMemObject(e12v_device);
  ret = clReleaseMemObject(e12t_device);
  ret = clReleaseMemObject(gphiu_device);
  ret = clReleaseMemObject(gphiv_device);
  ret = clReleaseMemObject(tmask_device);

  ret = clReleaseCommandQueue(command_queue);
  ret = clReleaseContext(context);

  if(image_file){
    free(image_file);
  }
}
