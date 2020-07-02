#ifndef __OPENCL_VERSION__  // If its not an OpenCL Kernel
#include <stdio.h>
#include <math.h>

#ifdef OPENCL_HOST // If it is OpenCL infrastructure

void set_args_momv(cl_kernel kern,
           cl_int *nx,
           cl_int *xstop,
           cl_mem *va_device,
           cl_mem *un_device,
           cl_mem *vn_device,
           cl_mem *hu_device,
           cl_mem *hv_device,
           cl_mem *ht_device,
           cl_mem *ssha_v_device,
           cl_mem *sshn_device,
           cl_mem *sshn_u_device,
           cl_mem *sshn_v_device,
           cl_mem *tmask_device,
           cl_mem *e1v_device, cl_mem *e1t_device,
           cl_mem *e2u_device, cl_mem *e2v_device,
           cl_mem *e2t_device, cl_mem *e12v_device,
           cl_mem *gphiv_device,
           cl_double *rdt, cl_double *cbfr, cl_double *visc,
           cl_double *omega, cl_double *d2r, cl_double *g){
  cl_int ret;
  cl_int arg_idx = 0;

  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_int), (void *)nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_int), (void *)xstop);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)va_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)un_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)vn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)ht_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)ssha_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)tmask_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)e1v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)e1t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)e2u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)e2v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)e2t_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)e12v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_mem), (void *)gphiv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_double), (void *)rdt);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_double), (void *)cbfr);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_double), (void *)visc);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_double), (void *)omega);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_double), (void *)d2r);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(kern, arg_idx++, sizeof(cl_double), (void *)g);
  check_status("clSetKernelArg", ret);

  fprintf(stdout, "Set %d arguments for Momentum-v kernel\n", arg_idx);

}

#endif
#endif


#ifdef __OPENCL_VERSION__
/** Interface to OpenCL version of kernel */
__kernel void momentum_v_code(int width, int xstop,
                  __global double* restrict va,
                  const __global double* restrict un,
                  const __global double* restrict vn, 
                  const __global double* restrict hu,
                  const __global double* restrict hv,
                  const __global double* restrict ht,
                  const __global double* restrict ssha_v, 
                  const __global double* restrict sshn,
                  const __global double* restrict sshn_u,
                  const __global double* restrict sshn_v, 
                  const __global int* restrict tmask,
                  const __global double* restrict e1v,
                  const __global double* restrict e1t,
                  const __global double* restrict e2u,
                  const __global double* restrict e2v,
                  const __global double* restrict e2t,
                  const __global double* restrict e12v,
                  const __global double* restrict gphiv,
                  double rdt, double cbfr, double visc,
                  double omega, double d2r, double g){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if(ji > xstop)return;
#else
/** Interface to standard C version of kernel */
#if defined(KOKKOS_INLINE_FUNCTION)
KOKKOS_INLINE_FUNCTION
#else
inline
#endif
void momentum_v_code(int ji, int jj, int width,
             double *va, double *un, double *vn, 
             double *hu, double *hv, double *ht, double *ssha_v, 
             double *sshn, double *sshn_u, double *sshn_v, 
             int *tmask, double *e1v, double *e1t, double *e2u,
             double *e2v, double *e2t, double *e12v, double *gphiv,
             double rdt, double cbfr, double visc,
             double omega, double d2r, double g){
#endif

  double u_e, u_w, v_n, v_s;
  double u_ec, u_wc, vv_e, vv_n, vv_s, vv_w;
  double depe, depw, deps, depn;
  double hpg, adv, cor, vis;
  double dvdx_e, dvdx_w, dvdy_n, dvdy_s;
  
  int idxim1, idxjm1, idxip1, idxjp1, idxim1jp1;
  int idx = jj*width + ji;
  
#ifdef SIMPLE_MOMENTUM
  va[idx] = (vn[idx] * (hv[idx] + sshn_v[idx]) + 
           rdt * (1.0) / e12v[idx] ) / 
    ((hv[idx] + ssha_v[idx])) / (1.0 + cbfr * rdt) ;
#else
  idxim1 = idx - 1;
  idxip1 = idx + 1;
  idxjm1 = idx - width;
  idxjp1 = idx + width;
  idxim1jp1 = idx + width - 1;

  if(tmask[idx] + tmask[idxip1] <= 0)  return; // jump over non-computatinal domain
  if(tmask[idx] <= 0 || tmask[idxjp1] <= 0) return; // jump over v boundary cells

  // kernel v adv 
  v_n  = 0.5 * (vn[idx] + vn[idxjp1]) * e1t[idxjp1]; // add length scale.
  depn = ht[idxjp1] + sshn[idxjp1];

  v_s  = 0.5 * (vn[idx] + vn[idxjm1]) * e1t[idx]; // add length scale
  deps = ht[idx] + sshn[idx];

  u_wc = 0.5 * (un[idxim1] + un[idxim1jp1]);
  u_w  = 0.5 * u_wc * (e2u[idxim1] + e2u[idxim1jp1]);
  depw = 0.50 * (hu[idxim1] + sshn_u[idxim1] + 
         hu[idxim1jp1] + sshn_u[idxim1jp1]);

  u_ec = 0.5 * (un[idx] + un[idxjp1]);
  u_e  = 0.5 * u_ec * (e2u[idx] + e2u[idxjp1]);
  depe = 0.50 * (hu[idx] + sshn_u[idx] + 
         hu[idxjp1] + sshn_u[idxjp1]);

  // -advection (currently first order upwind)
  vv_s = (0.5 - copysign(0.5, v_s)) * vn[idx] +  
    (0.5 + copysign(0.5, v_s)) * vn[idxjm1] ;
  vv_n = (0.5 + copysign(0.5, v_n)) * vn[idx] +  
    (0.5 - copysign(0.5, v_n)) * vn[idxjp1] ;

  if(tmask[idxim1] <= 0 || tmask[idxim1jp1] <= 0){
    vv_w = (0.5 - copysign(0.5, u_w)) * vn[idx];
  }
  else{
    vv_w = (0.5 - copysign(0.5, u_w)) * vn[idx] +  
      (0.5 + copysign(0.5, u_w)) * vn[idxim1] ;
  }

  if(tmask[idxip1] <= 0 || tmask[idxjp1+1] <= 0){
    vv_e = (0.5 + copysign(0.5, u_e)) * vn[idx];
  }
  else{
    vv_e = (0.5 + copysign(0.5, u_e)) * vn[idx] +
      (0.5 - copysign(0.5, u_e)) * vn[idxip1];
  }

  adv = vv_w * u_w * depw - vv_e * u_e * depe + 
    vv_s * v_s * deps - vv_n * v_n * depn;

  // -viscosity
    
  dvdy_n = (vn[idxjp1] - vn[idx]) / e2t[idxjp1] * 
    (ht[idxjp1] + sshn[idxjp1]);
  dvdy_s = (vn[idx] - vn[idxjm1]) / e2t[idx] * 
    (ht[idx] + sshn[idx]);

  if(tmask[idxim1] <= 0 || tmask[idxim1jp1] <= 0){
    dvdx_w = 0.0; // slip boundary
  }
  else{
    dvdx_w = (vn[idx] - vn[idxim1]) / (e1v[idx] + e1v[idxim1]) * 
      (hv[idx] + sshn_v[idx] + hv[idxim1] + sshn_v[idxim1]);
  }

  if(tmask[idxip1] <= 0 || tmask[idxjp1+1] <= 0){
    dvdx_e = 0.0; // slip boundary
  }
  else{
    dvdx_e = (vn[idxip1] - vn[idx]) / (e1v[idx] + e1v[idxip1]) * 
      (hv[idx] + sshn_v[idx] + hv[idxip1] + sshn_v[idxip1]);
  }

  vis = (dvdy_n - dvdy_s ) * e1v[idx]  + 
    (dvdx_e - dvdx_w ) * e2v[idx] * 0.5  ;

  vis = visc * vis;   // visc will be a array visc(1:jpijglou) 
  // for variable viscosity, such as turbulent viscosity

  // -Coriolis' force (can be implemented implicitly)
  cor = -0.5*(2. * omega * sin(gphiv[idx] * d2r) * (u_ec + u_wc)) * 
    e12v[idx] * (hv[idx] + sshn_v[idx]);

  // -pressure gradient
  hpg = -g * (hv[idx] + sshn_v[idx]) * e1v[idx] * 
    (sshn[idxjp1] - sshn[idx]);

  // -linear bottom friction (implemented implicitly.
  va[idx] = (vn[idx] * (hv[idx] + sshn_v[idx]) + 
           rdt * (adv + vis + cor + hpg) / e12v[idx] ) / 
    ((hv[idx] + ssha_v[idx])) / (1.0 + cbfr * rdt) ;
#endif
}
