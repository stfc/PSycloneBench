#ifndef __OPENCL_VERSION__
// This header isn't available/required in OpenCL
#include <math.h>
#include <stdio.h>
#else
#include "opencl_utils.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

/** Set the arguments for the OpenCL kernel */
void set_args_momu(cl_kernel kern,
		   cl_int *nx,
		   cl_mem *ua_device,
		   cl_mem *un_device,
		   cl_mem *vn_device,
		   cl_mem *hu_device,
		   cl_mem *hv_device,
		   cl_mem *ht_device,
		   cl_mem *ssha_u_device,
		   cl_mem *sshn_device,
		   cl_mem *sshn_u_device,
		   cl_mem *sshn_v_device,
		   cl_mem *tmask_device,
		   cl_mem *e1u_device, cl_mem *e1v_device,
		   cl_mem *e1t_device, cl_mem *e2u_device,
		   cl_mem *e2t_device, cl_mem *e12u_device,
		   cl_mem *gphiu_device,
		   cl_double *rdt, cl_double *cbfr, cl_double *visc){
  cl_int ret;
  cl_int arg_idx = 0;
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_int), (void *)nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)ua_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)un_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)vn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)ht_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)ssha_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)e1u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)e1v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)e1t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)e2u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)e2t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)e12u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem),
		       (void *)gphiu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_double),
		       (void *)rdt);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_double),
		       (void *)cbfr);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_double),
		       (void *)visc);
  check_status("clSetKernelArg", ret);
  fprintf(stdout, "Set %d arguments for Momentum-u kernel\n", arg_idx);
}

#endif

/*
  type, extends(kernel_type) :: momentum_u
     type(arg), dimension(18) :: meta_args =  &
          (/ arg(READWRITE, CU, POINTWISE),  & ! ua
             arg(READ,      CU, POINTWISE),  & ! un
             arg(READ,      CV, POINTWISE),  & ! vn
             arg(READ,      CU, POINTWISE),  & ! hu
             arg(READ,      CV, POINTWISE),  & ! hv
             arg(READ,      CT, POINTWISE),  & ! ht
             arg(READ,      CU, POINTWISE),  & ! ssha_u
             arg(READ,      CT, POINTWISE),  & ! sshn_t
             arg(READ,      CU, POINTWISE),  & ! sshn_u
             arg(READ,      CV, POINTWISE),  & ! sshn_v
             arg(READ,      GRID_MASK_T),    &
             arg(READ,      GRID_DX_U),      &
             arg(READ,      GRID_DX_V),      &
             arg(READ,      GRID_DX_T),      &
             arg(READ,      GRID_DY_U),      &
             arg(READ,      GRID_DY_T),      &
             arg(READ,      GRID_AREA_U),    &
             arg(READ,      GRID_LAT_U)      &
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
    procedure, nopass :: code => momentum_u_code
  end type momentum_u
*/


#ifdef __OPENCL_VERSION__
/** Interface to OpenCL version of kernel */
__kernel void momentum_u_code(int width,
			      __global double* restrict ua,
			      const __global double* restrict un,
			      const __global double* restrict vn,
			      const __global double* restrict hu,
			      const __global double* restrict hv,
			      const __global double* restrict ht,
			      const __global double* restrict ssha_u,
			      const __global double* restrict sshn,
			      const __global double* restrict sshn_u,
			      const __global double* restrict sshn_v,
			      const __global int* restrict tmask,
			      const __global double* restrict e1u,
			      const __global double* restrict e1v,
			      const __global double* restrict e1t,
			      const __global double* restrict e2u,
			      const __global double* restrict e2t,
			      const __global double* restrict e12u,
			      const __global double* restrict gphiu,
			      double rdt, double cbfr, double visc){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
/** Interface to standard C version of kernel */
inline void momentum_u_code(int ji, int jj, int width,
		     double *ua, double *un, double *vn,
		     double *hu, double *hv, double *ht, double *ssha_u,
		     double *sshn, double *sshn_u, double *sshn_v,
		     int *tmask,
		     double *e1u, double *e1v, double *e1t,
		     double *e2u, double *e2t, double *e12u, double *gphiu,
		     double rdt, double cbfr, double visc, double omega, double d2r, double g){
#endif
  double u_e, u_w, v_n, v_s;
  double v_nc, v_sc;
  double depe, depw, deps, depn;
  double hpg, adv, cor, vis;
  double dudx_e, dudx_w;
  double dudy_s, dudy_n;
  double uu_e, uu_n, uu_s, uu_w;
  
  int idxim1, idxjm1, idxip1, idxjp1, idxip1jm1;
  int idx = jj*width + ji;
  
#ifdef __OPENCL_VERSION__
  //  int nrow = (int)get_global_size(1);
  //if(ji==0 || ji > (width-2))return;
  //if(jj==0 || jj > (nrow-2))return;
#endif

#ifdef SIMPLE_MOMENTUM
  ua[idx] = (un[idx] * (hu[idx] + sshn_u[idx]) + rdt * 
	     (1.0) / e12u[idx]) / 
    (hu[idx] + ssha_u[idx]) / (1.0 + cbfr * rdt) ;
#else
  idxim1 = idx - 1;
  idxip1 = idx + 1;
  idxjm1 = idx - width;
  idxjp1 = idx + width;
  idxip1jm1 = idx - width + 1;

  if(tmask[idx] + tmask[idxip1] <= 0) return; // jump over non-computational domain
  if(tmask[idx] <= 0 || tmask[idxip1] <= 0) return; // jump over boundary u

  u_e  = 0.5 * (un[idx] + un[idxip1]) * e2t[idxip1];   // add length scale.
  depe = ht[idxip1] + sshn[idxip1];

  u_w  = 0.5 * (un[idx] + un[idxim1]) * e2t[idx];   // add length scale
  depw = ht[idx] + sshn[idx];

  v_sc = 0.5 * (vn[idxjm1] + vn[idxip1jm1]);
  v_s  = 0.5 * v_sc * (e1v[idxjm1] + e1v[idxip1jm1]);
  deps = 0.5 * (hv[idxjm1] + sshn_v[idxjm1] + hv[idxip1jm1] +
		sshn_v[idxip1jm1]);

  v_nc = 0.5 * (vn[idx] + vn[idxip1]);
  v_n  = 0.5 * v_nc * (e1v[idx] + e1v[idxip1]);
  depn = 0.5 * (hv[idx] + sshn_v[idx] + hv[idxip1] +
		sshn_v[idxip1]);

  // -advection (currently first order upwind)
  uu_w = (0.5 - copysign(0.5, u_w)) * un[idx] +
    (0.5 + copysign(0.5, u_w)) * un[idxim1] ;
  uu_e = (0.5 + copysign(0.5, u_e)) * un[idx] +
    (0.5 - copysign(0.5, u_e)) * un[idxip1] ;

  if(tmask[idxjm1] <=0 || tmask[idxip1jm1] <= 0){
    uu_s = (0.5 - copysign(0.5, v_s)) * un[idx];
  }
  else{
    uu_s = (0.5 - copysign(0.5, v_s)) * un[idx] +
      (0.5 + copysign(0.5, v_s)) * un[idxjm1] ;
  }

  if(tmask[idxjp1] <=0 || tmask[idxjp1+1] <= 0){
    uu_n = (0.5 + copysign(0.5, v_n)) * un[idx];
  }
  else{
    uu_n = (0.5 + copysign(0.5, v_n)) * un[idx] +
      (0.5 - copysign(0.5, v_n)) * un[idxjp1];
  }

  adv = uu_w * u_w * depw - uu_e * u_e * depe +
    uu_s * v_s * deps - uu_n * v_n * depn;

  // -viscosity

  dudx_e = (un[idxip1] - un[idx]) / e1t[idxip1]*(ht[idxip1] + sshn[idxip1]);
  dudx_w = (un[idx] - un[idxim1]) / e1t[idx]*(ht[idx] + sshn[idx]);
  if(tmask[idxjm1] <=0 || tmask[idxip1jm1] <= 0){
    dudy_s = 0.0; // slip boundary
  }
  else{
    dudy_s = (un[idx] - un[idxjm1]) / (e2u[idx] + e2u[idxjm1]) * 
      (hu[idx] + sshn_u[idx] + hu[idxjm1] + sshn_u[idxjm1]);
  }

  if(tmask[idxjp1] <= 0 || tmask[idxjp1+1] <= 0){
    dudy_n = 0.0; // slip boundary
  }
  else{
    dudy_n = (un[idxjp1] - un[idx]) / (e2u[idx] + e2u[idxjp1]) * 
      (hu[idx] + sshn_u[idx] + hu[idxjp1] + sshn_u[idxjp1]);
  }

  vis = (dudx_e - dudx_w ) * e2u[idx]  + 
    (dudy_n - dudy_s ) * e1u[idx] * 0.5;
  vis = visc * vis;   // visc will be an array visc(1:jpijglou) 
  // for variable viscosity, such as turbulent viscosity

  // -Coriolis' force (can be implemented implicitly)
  cor = 0.5 * (2. * omega * sin(gphiu[idx] * d2r) * (v_sc + v_nc)) * 
    e12u[idx] * (hu[idx] + sshn_u[idx]);

  // -pressure gradient
  hpg = -g * (hu[idx] + sshn_u[idx]) * e2u[idx] * 
    (sshn[idxip1] - sshn[idx]);

  // -linear bottom friction (implemented implicitly.
  ua[idx] = (un[idx] * (hu[idx] + sshn_u[idx]) + rdt * 
                 (adv + vis + cor + hpg) / e12u[idx]) / 
    (hu[idx] + ssha_u[idx]) / (1.0 + cbfr * rdt) ;
#endif
}
