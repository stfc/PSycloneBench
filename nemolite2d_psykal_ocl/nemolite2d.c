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

void write_field(char *filename, int nx, int ny,
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

int main(){
  /** The version of OpenCL supported by the selected device */
  int cl_version;
  cl_device_id device_ids[MAX_DEVICES];
  cl_device_id *device;
  int idev;
  cl_context context = NULL;
  cl_command_queue command_queue = NULL;

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
  cl_event event;
  //cl_command_queue_properties queue_properties;
  /** Problem size */
  cl_int nx = 128;
  cl_int ny = 128;
  /* Extend domain by one in each dimension to allow for staggering */
  nx += 1;
  ny += 1;
  /** Our time-step index (passed into BCs kernel) */
  cl_int istep;
  /** Number of time-steps to do */
  cl_int nsteps = 1;
  int ji, jj, idx;
  int buff_size;
  /** Sea-surface height */
  cl_double *ssha, *ssha_u, *ssha_v, *sshn, *sshn_u, *sshn_v;
  cl_double *hu, *hv, *ht, *un, *vn, *ua, *va;
  cl_double *gphiu, *gphiv;
  cl_double *e1u, *e1v, *e1t, *e2u, *e2v, *e2t, *e12t, *e12u, *e12v;
  /** T-point mask */
  cl_int *tmask;
  cl_double dep_const = 2.0;
  /** Horizontal grid resolution */
  cl_double dx = 0.5, dy = 0.5;
  /** For computing checksums for validation */
  double sum;
  /** Reciprocal of the time step */
  cl_double rdt = 0.5;
  /** horiz. kinematic viscosity coeff. */
  cl_double visc = 100.0;
  /** Coefficient of bottom friction */
  cl_double cbfr = 0.001;

  TimerInit();

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

  /* Create OpenCL Kernels */
  cont_kernel = build_kernel(&context, device,
			     "./continuity_kern.c", "continuity_code");
  momu_kernel = build_kernel(&context, device,
			     "./momentum_kern.c", "momentum_u_code");
  momv_kernel = build_kernel(&context, device,
			     "./momentum_kern.c", "momentum_v_code");
  bc_ssh_kernel = build_kernel(&context, device,
  			       "./boundary_conditions_kern.c", "bc_ssh_code");
  bc_ssh_kernel = build_kernel(&context, device,
  			       "./boundary_conditions_kern.c", "bc_ssh_code");
  bc_solid_u_kernel = build_kernel(&context, device,
				   "./boundary_conditions_kern.c",
				   "bc_solid_u_code");
  bc_solid_v_kernel = build_kernel(&context, device,
				   "./boundary_conditions_kern.c",
				   "bc_solid_v_code");
  bc_flather_u_kernel = build_kernel(&context, device,
				     "./boundary_conditions_kern.c",
				     "bc_flather_u_code");
  bc_flather_v_kernel = build_kernel(&context, device,
				     "./boundary_conditions_kern.c",
				     "bc_flather_v_code");
  next_sshu_kernel = build_kernel(&context, device,
				  "./time_update_kern.c",
				  "next_sshu_code");
  next_sshv_kernel = build_kernel(&context, device,
				  "./time_update_kern.c",
				  "next_sshv_code");

  /* Create Device Memory Buffers */
  buff_size = nx*ny*sizeof(cl_double);
  ssha_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			       NULL, &ret);
  check_status("clCreateBuffer", ret);
  ssha_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
				 NULL, &ret);
  check_status("clCreateBuffer", ret);
  ssha_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
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
  ht_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  check_status("clCreateBuffer", ret);

  /* Velocity fields */
  un_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  check_status("clCreateBuffer", ret);
  vn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  check_status("clCreateBuffer", ret);
  ua_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  check_status("clCreateBuffer", ret);
  va_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			     NULL, &ret);
  check_status("clCreateBuffer", ret);
  
  /* Mesh scale factors */
  e1t_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			      NULL, &ret);
  check_status("clCreateBuffer", ret);
  e1u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			      NULL, &ret);
  check_status("clCreateBuffer", ret);
  e1v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			      NULL, &ret);
  check_status("clCreateBuffer", ret);
  e2u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			      NULL, &ret);
  check_status("clCreateBuffer", ret);
  e2v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			      NULL, &ret);
  check_status("clCreateBuffer", ret);
  e2t_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			      NULL, &ret);
  check_status("clCreateBuffer", ret);
  e12t_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);
  e12u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);

  e12v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                          NULL, &ret);
  check_status("clCreateBuffer", ret);

  tmask_device = clCreateBuffer(context, CL_MEM_READ_WRITE,
				nx*ny*sizeof(cl_int), NULL, &ret);
  check_status("clCreateBuffer", ret);

  /* Coriolis parameters */
  gphiu_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
				NULL, &ret);
  check_status("clCreateBuffer", ret);
  gphiv_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
				NULL, &ret);
  check_status("clCreateBuffer", ret);
  fprintf(stdout, "Created device buffers OK\n");

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
  arg_idx = 0;
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_int), (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ua_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&un_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&vn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ht_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ssha_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e1u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e1v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e1t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e2u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e2t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e12u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&gphiu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_double),
		       (void *)&rdt);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_double),
		       (void *)&cbfr);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momu_kernel, arg_idx++, sizeof(cl_double),
		       (void *)&visc);
  check_status("clSetKernelArg", ret);
  fprintf(stdout, "Set %d arguments for Momentum-u kernel\n", arg_idx);
  
  /* Set OpenCL Kernel Parameters for Momentum-v */
  arg_idx = 0;
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_int), (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&va_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&un_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&vn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ht_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&ssha_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e1v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e1t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e2u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e2v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e2t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e12v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&gphiv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_double),
		       (void *)&rdt);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_double),
		       (void *)&cbfr);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(momv_kernel, arg_idx++, sizeof(cl_double),
		       (void *)&visc);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for Momentum-v kernel\n", arg_idx);

  /* Set kernel arguments for bc_ssh */
  arg_idx = 0;
  ret = clSetKernelArg(bc_ssh_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(bc_ssh_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&istep);
  check_status("clSetKernelArg", ret);
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
  arg_idx = 0;
  ret = clSetKernelArg(next_sshu_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e12t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshu_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e12u_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for next_sshu kernel\n", arg_idx);
  
  /* Set OpenCL Kernel Parameters for next_sshv kernel */
  arg_idx = 0;
  ret = clSetKernelArg(next_sshv_kernel, arg_idx++, sizeof(cl_int),
		       (void *)&nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e12t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(next_sshv_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)&e12v_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for next_sshv kernel\n", arg_idx);

  TimerStop();

  /*------------------------------------------------------------*/
  /* Field initialisation on host */
  
  ssha = malloc(buff_size);
  ssha_u = malloc(buff_size);
  ssha_v = malloc(buff_size);
  sshn = malloc(buff_size);
  sshn_u = malloc(buff_size);
  sshn_v = malloc(buff_size);
  hu = malloc(buff_size);
  hv = malloc(buff_size);
  ht = malloc(buff_size);
  un = malloc(buff_size);
  vn = malloc(buff_size);
  ua = malloc(buff_size);
  va = malloc(buff_size);
  e1u = malloc(buff_size);
  e1v = malloc(buff_size);
  e1t = malloc(buff_size);
  e2u = malloc(buff_size);
  e2v = malloc(buff_size);
  e2t = malloc(buff_size);
  e12u = malloc(buff_size);
  e12v = malloc(buff_size);
  e12t = malloc(buff_size);
  tmask = malloc(nx*ny*sizeof(cl_int));

  gphiu = malloc(buff_size);
  gphiv = malloc(buff_size);

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
      un[idx] = ji;
      vn[idx] = jj;
      sshn_u[idx] = 0.0;
      sshn_v[idx] = 2.0;
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
    tmask[idx] = 0;
    // East solid boundary
    for(ji=xstop-1; ji<nx; ji++){
      tmask[idx+ji] = 0;
    }
  }
  // Southern open boundary
  for(ji=0;ji<nx;ji++){
    tmask[ji] = -1;
  }
  
  /*------------------------------------------------------------*/
  /* Copy data to device synchronously */
  TimerStart("Write buffers to device");
  
  ret = clEnqueueWriteBuffer(command_queue, ssha_device, 1, 0,
			     (size_t)buff_size, (void *)ssha, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, ssha_u_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_u, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, ssha_v_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_v, 0,
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
  ret = clEnqueueWriteBuffer(command_queue, ht_device, 1, 0,
			     (size_t)buff_size, (void *)ht, 0,
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
  ret = clEnqueueWriteBuffer(command_queue, ua_device, 1, 0,
			     (size_t)buff_size, (void *)ua, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, va_device, 1, 0,
			     (size_t)buff_size, (void *)va, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e1u_device, 1, 0,
			     (size_t)buff_size, (void *)e1u, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e1v_device, 1, 0,
			     (size_t)buff_size, (void *)e1v, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e1t_device, 1, 0,
			     (size_t)buff_size, (void *)e1t, 0,
			     NULL, &event);
  ret = clEnqueueWriteBuffer(command_queue, e2u_device, 1, 0,
			     (size_t)buff_size, (void *)e2u, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e2v_device, 1, 0,
			     (size_t)buff_size, (void *)e2v, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e2t_device, 1, 0,
			     (size_t)buff_size, (void *)e2t, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e12u_device, 1, 0,
			     (size_t)buff_size, (void *)e12t, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e12v_device, 1, 0,
			     (size_t)buff_size, (void *)e12t, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, e12t_device, 1, 0,
			     (size_t)buff_size, (void *)e12t, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, gphiu_device, 1, 0,
			     (size_t)buff_size, (void *)gphiu, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, gphiv_device, 1, 0,
			     (size_t)buff_size, (void *)gphiv, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);
  ret = clEnqueueWriteBuffer(command_queue, tmask_device, 1, 0,
			     (size_t)(nx*ny*sizeof(cl_int)), (void *)tmask, 0,
			     NULL, &event);
  check_status("clEnqueueWriteBuffer", ret);

  TimerStop();
  
  /*------------------------------------------------------------*/
  /* Run the kernels */
  size_t global_size[2] = {xstop, ystop};

  TimerStart("Time-stepping, OpenCL");
  
  for(istep=1; istep<=nsteps; istep++){
    
    ret = clEnqueueNDRangeKernel(command_queue, cont_kernel, 2, 0,
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);
    ret = clEnqueueNDRangeKernel(command_queue, momu_kernel, 2, 0,
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);
    ret = clEnqueueNDRangeKernel(command_queue, momv_kernel, 2, 0,
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);
    ret = clEnqueueNDRangeKernel(command_queue, bc_ssh_kernel, 2, 0,
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);
    ret = clEnqueueNDRangeKernel(command_queue, bc_solid_u_kernel, 2, 0,
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);
    ret = clEnqueueNDRangeKernel(command_queue, bc_solid_v_kernel, 2, 0,
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);
    ret = clEnqueueNDRangeKernel(command_queue, bc_flather_u_kernel, 2, 0,
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);
    ret = clEnqueueNDRangeKernel(command_queue, bc_flather_v_kernel, 2, 0,
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);

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
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);
  
    ret = clEnqueueNDRangeKernel(command_queue, next_sshv_kernel, 2, 0,
				 global_size, NULL, 0,0,0);
    check_status("clEnqueueNDRangeKernel", ret);
  }

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

  printf("ssha[1,1] [1,2] = %e, %e\n", ssha[nx*ystart+xstart],
	 ssha[nx*ystart+xstart+1]);

  write_field("ssha_cpu.dat", nx, ny, 1, 1, ssha);
  
  sum = checksum(ssha, nx, xstop, ystop, xstart, ystart);
  fprintf(stdout, "ssha checksum on CPU = %e\n", sum);
  sum = checksum(ua, nx, xstop-1, ystop, xstart, ystart);
  fprintf(stdout, "ua checksum on CPU = %e\n", sum);
  sum = checksum(va, nx, xstop, ystop-1, xstart, ystart);
  fprintf(stdout, "va checksum on CPU = %e\n", sum);

  /* Copy data back from device, synchronously */
  clEnqueueReadBuffer(command_queue, ssha_device, 1, 0,
		      (size_t)buff_size, (void *)ssha, 0, NULL, NULL);
  check_status("clEnqueueReadBuffer", ret);
  clEnqueueReadBuffer(command_queue, ua_device, 1, 0,
		      (size_t)buff_size, (void *)ua, 0, NULL, NULL);
  check_status("clEnqueueReadBuffer", ret);
  clEnqueueReadBuffer(command_queue, va_device, 1, 0,
		      (size_t)buff_size, (void *)va, 0, NULL, NULL);
  check_status("clEnqueueReadBuffer", ret);

  /* Compute and output a checksum */
  sum = checksum(ssha, nx, xstop, ystop, xstart, ystart);
  printf("ssha[1,1] [1,2] = %e, %e\n", ssha[nx*ystart+xstart],
	 ssha[nx*ystart+xstart+1]);

  write_field("ssha_ocl.dat", nx, ny, 1, 1, ssha);

  fprintf(stdout, "ssha checksum from OpenCL = %e\n", sum);
  sum = checksum(ua, nx, xstop-1, ystop, xstart, ystart);
  fprintf(stdout, "ua checksum from OpenCL = %e\n", sum);
  sum = checksum(va, nx, xstop, ystop-1, xstart, ystart);
  fprintf(stdout, "va checksum from OpenCL = %e\n", sum);

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

}
