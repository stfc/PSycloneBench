#ifndef __OPENCL_VERSION__  // If its not an OpenCL Kernel
#include <stdio.h>

#ifdef OPENCL_HOST // If it is OpenCL infrastructure

void set_args_next_sshv(cl_kernel kern,
			cl_int *nx,
			cl_int *xstop,
			cl_mem* sshn_v_device,
			cl_mem* sshn_device,
			cl_mem* tmask_device,
			cl_mem* e12t_device,
			cl_mem* e12v_device){
  cl_int ret;
  int arg_idx = 0;

  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_int),
		       (void *)nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_int),
		       (void *)xstop);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)e12t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)e12v_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for next_sshv kernel\n", arg_idx);

}

void set_args_next_sshu(cl_kernel kern,
			cl_int *nx,
			cl_int *xstop,
			cl_mem *sshn_u_device,
			cl_mem *sshn_device,
			cl_mem *tmask_device,
			cl_mem *e12t_device,
			cl_mem *e12u_device){
  cl_int ret;
  int arg_idx = 0;
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_int), (void *)nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_int), (void *)xstop);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)e12t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)e12u_device);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for next_sshu kernel\n", arg_idx);
}

#endif
#endif


#ifdef __OPENCL_VERSION__
__kernel void next_sshu_code(int width, int xstop,
			     __global double* restrict sshn_u,
			     __global double* restrict sshn,
			     __global int* restrict tmask,
			     __global double* restrict e12t,
			     __global double* restrict e12u){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if (ji > xstop) return;
#else
#if defined(KOKKOS_INLINE_FUNCTION)
KOKKOS_INLINE_FUNCTION
#else
inline
#endif
void next_sshu_code(int ji, int jj, int width,
		      double* sshn_u,
		      const double* sshn,
		      const int* tmask,
		      const double* e12t,
		      const double* e12u){
#endif
  double rtmp1;
  int idx = jj*width + ji;
  int idxip1 = idx + 1;

  if(tmask[idx] + tmask[idxip1] <= 0)return; // jump over non-computational domain

  if(tmask[idx] * tmask[idxip1] > 0){
    rtmp1 = e12t[idx] * sshn[idx] + e12t[idxip1] * sshn[idxip1];
    sshn_u[idx] = 0.5 * rtmp1 / e12u[idx] ;
  }
  else if(tmask[idx] <= 0){
    sshn_u[idx] = sshn[idxip1];
  }
  else if(tmask[idxip1] <= 0){
      sshn_u[idx] = sshn[idx];
  }

}

#ifdef __OPENCL_VERSION__
__kernel void next_sshv_code(int width, int xstop,
			     __global double* restrict sshn_v,
			     __global double* restrict sshn,
			     __global int* restrict tmask,
			     __global double* restrict e12t,
			     __global double* restrict e12v){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if (ji > xstop) return;
#else
#if defined(KOKKOS_INLINE_FUNCTION)
KOKKOS_INLINE_FUNCTION
#else
inline
#endif
void next_sshv_code(int ji, int jj, int width,
		    double* sshn_v,
		    double* sshn,
		    int* tmask,
		    double* e12t,
		    double* e12v){
#endif
  double rtmp1;
  int idx = jj*width + ji;
  int idxjp1 = idx + width;
    
  if((tmask[idx] + tmask[idxjp1]) <= 0)return; //jump over non-computational domain
  if((tmask[idx] * tmask[idxjp1]) > 0){
    rtmp1 = e12t[idx] * sshn[idx] + e12t[idxjp1] * sshn[idxjp1];
    sshn_v[idx] = 0.5 * rtmp1 / e12v[idx] ;
  }
  else if(tmask[idx] <= 0){
    sshn_v[idx] = sshn[idxjp1];
  }
  else if(tmask[idxjp1] <= 0){
    sshn_v[idx] = sshn[idx];
  }
  
}
