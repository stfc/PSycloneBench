#ifndef __OPENCL_VERSION__
#include <stdio.h>
#else
#include "opencl_utils.h"
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

void set_args_next_sshv(cl_kernel kern,
			cl_int *nx,
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
			cl_mem *sshn_u_device,
			cl_mem *sshn_device,
			cl_mem *tmask_device,
			cl_mem *e12t_device,
			cl_mem *e12u_device){
  cl_int ret;
  int arg_idx = 0;
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_int), (void *)nx);
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

/*
  type, extends(kernel_type) :: next_sshu
     type(arg), dimension(5) :: meta_args =  &
          (/ arg(READWRITE, CU, POINTWISE),  &
             arg(READ,      CU, POINTWISE),  &
             arg(READ,      GRID_MASK_T),    &
             arg(READ,      GRID_AREA_T),    &
             arg(READ,      GRID_AREA_U)     &
           /)

     !> We update only those points within the internal region
     !! of the simulated domain.
     integer :: ITERATES_OVER = INTERNAL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => next_sshu_code
  end type next_sshu
*/

#ifdef __OPENCL_VERSION__
__kernel void next_sshu_code(int width,
			     __global double* restrict sshn_u,
			     __global double* restrict sshn,
			     __global int* restrict tmask,
			     __global double* restrict e12t,
			     __global double* restrict e12u){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
inline void next_sshu_code(int ji, int jj, int width,
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

/*
  type, extends(kernel_type) :: next_sshv
     type(arg), dimension(5) :: meta_args =  &
          (/ arg(READWRITE, CV, POINTWISE),  &
             arg(READ,      CV, POINTWISE),  &
             arg(READ,      GRID_MASK_T),    &
             arg(READ,      GRID_AREA_T),    &
             arg(READ,      GRID_AREA_V)     &
           /)

     !> We update only those points within the internal region
     !! of the simulated domain.
     integer :: ITERATES_OVER = INTERNAL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => next_sshv_code
  end type next_sshv
*/
  
#ifdef __OPENCL_VERSION__
__kernel void next_sshv_code(int width,
			     __global double* restrict sshn_v,
			     __global double* restrict sshn,
			     __global int* restrict tmask,
			     __global double* restrict e12t,
			     __global double* restrict e12v){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
inline void next_sshv_code(int ji, int jj, int width,
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
