#ifndef __OPENCL_VERSION__
#include <stdio.h>
#include "opencl_utils.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

void set_args_continuity(cl_kernel cont_kernel,
			 cl_int *nx,
			 cl_mem *ssha_device,
			 cl_mem *sshn_device,
			 cl_mem *sshn_u_device,
			 cl_mem *sshn_v_device,
			 cl_mem *hu_device,
			 cl_mem *hv_device,
			 cl_mem *un_device,
			 cl_mem *vn_device,
			 cl_double *rdt,
			 cl_mem *e12t_device){
  cl_int ret;
  int arg_idx = 0;
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_int), (void *)nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)ssha_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)un_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)vn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_double),
		       (void *)rdt);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)e12t_device);
  check_status("clSetKernelArg", ret);
  
  fprintf(stdout, "Set %d arguments for Continuity kernel\n", arg_idx);
}
#endif
/*
  type, extends(kernel_type) :: continuity
     type(arg), dimension(10) :: meta_args =    &
          (/ arg(WRITE, CT, POINTWISE),        & ! ssha
             arg(READ,  CT, POINTWISE),        & ! sshn
             arg(READ,  CU, POINTWISE),        & ! sshn_u
             arg(READ,  CV, POINTWISE),        & ! sshn_v
             arg(READ,  CU, POINTWISE),        & ! hu
             arg(READ,  CV, POINTWISE),        & ! hv
             arg(READ,  CU, POINTWISE),        & ! un
             arg(READ,  CV, POINTWISE),        & ! vn
             arg(READ,  TIME_STEP),            &
             arg(READ,  GRID_AREA_T)           &
           /)
     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => continuity_code
  end type continuity
*/

#ifdef __OPENCL_VERSION__
/** Interface to OpenCL version of kernel */
__kernel void continuity_code(int width,                     
			      __global double* restrict ssha,
			      __global double* restrict sshn,
			      __global double* restrict sshn_u,
			      __global double* restrict sshn_v,
			      __global double* restrict hu,
			      __global double* restrict hv,
			      __global double* restrict un,
			      __global double* restrict vn,
			      double rdt,
			      __global double* restrict e12t){
    int ji = get_global_id(0);
    int jj = get_global_id(1);
#else

/** Interface to standard C version of kernel */
void continuity_code(int ji, int jj,
		     int width,                     
		     double *ssha,
		     double *sshn,
		     double *sshn_u,
		     double *sshn_v,
		     double* hu,
		     double *hv,
		     double *un,
		     double *vn,
		     double rdt,
		     double *e12t){
#endif
    /* Locals */
    double rtmp1, rtmp2, rtmp3, rtmp4;
    int idxim1, idxjm1;
    int idx = jj*width + ji;

    if(jj == 0)return;

    idxim1 = idx - 1;
    idxjm1 = idx - width;

    rtmp1 = (sshn_u[idx] + hu[idx]) * un[idx];
    rtmp2 = (sshn_u[idxim1] + hu[idxim1]) * un[idxim1];
    rtmp3 = (sshn_v[idx] + hv[idx]) * vn[idx];
    rtmp4 = (sshn_v[idxjm1] + hv[idxjm1]) * vn[idxjm1];

    ssha[idx] = sshn[idx] + (rtmp2 - rtmp1 + rtmp4 - rtmp3) *
      rdt / e12t[idx];
    // Following line is for testing only
    //ssha[idx] = (double)idx;
  }
