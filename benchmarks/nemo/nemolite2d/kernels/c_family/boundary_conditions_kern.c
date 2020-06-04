#ifndef __OPENCL_VERSION__  // If its not an OpenCL Kernel
#include <stdio.h>
#include <math.h>

#ifdef OPENCL_HOST // If it is OpenCL infrastructure

#endif
#endif

#ifdef __OPENCL_VERSION__
__kernel void bc_ssh_code(int width,
              int istep,
              __global double* restrict ssha,
              __global int* restrict tmask,
              double rdt){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  int nrow = (int)get_global_size(1);
  if(ji==0 || ji > (width-2))return;
  if(jj==0 || jj > (nrow-2))return;
#else
inline void bc_ssh_code(int ji, int jj, int width,
         int istep, double *ssha, int *tmask, double rdt){
#endif
  int idx = jj*width + ji;

  double amp_tide, omega_tide, rtime;

  amp_tide   = 0.2;
  omega_tide = 2.0 * 3.14159 / (12.42 * 3600.0);
  rtime = istep * rdt;
  
  if(tmask[idx] <= 0) return;
  
  if(tmask[idx-width] < 0){
    ssha[idx] = amp_tide * sin(omega_tide * rtime);
  }
  else if(tmask[idx+width] < 0){
    ssha[idx] = amp_tide * sin(omega_tide * rtime);
  }
  else if(tmask[idx+1] < 0){
    ssha[idx] = amp_tide * sin(omega_tide * rtime);
  }
  else if(tmask[idx-1] < 0){
    ssha[idx] = amp_tide * sin(omega_tide * rtime);
  }
}
  

    /** Kernel to apply solid boundary conditions for u-velocity */
#ifdef __OPENCL_VERSION__
__kernel void bc_solid_u_code(int width,
                  __global double* restrict ua,
                  __global int* restrict tmask){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if(ji > (width-2))return;
#else
inline void bc_solid_u_code(int ji, int jj, int width, double *ua, int *tmask){
#endif
  int idx = jj*width + ji;

  if(tmask[idx] * tmask[idx+1] == 0){
    ua[idx] = 0.0;
  }

}
  
  /** Kernel to apply solid boundary conditions for v-velocity */
#ifdef __OPENCL_VERSION__
__kernel void bc_solid_v_code(int width,
                  __global double* restrict va,
                  __global int* restrict tmask){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  int nrow = (int)get_global_size(1);
  if(jj > (nrow-2))return;
#else
inline void bc_solid_v_code(int ji, int jj, int width, double *va, int *tmask){
#endif
  int idx = jj*width + ji;

  if(tmask[idx] * tmask[idx+width] == 0){
    va[idx] = 0.0;
  }

}
  
/** Kernel to apply Flather condition to U */
#ifdef __OPENCL_VERSION__
__kernel void bc_flather_u_code(int width,
                __global double* restrict ua,
                __global double* restrict hu,
                __global double* restrict sshn_u,
                __global int* restrict tmask,
                double g){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if(ji > (width-2))return;
#else
inline void bc_flather_u_code(int ji, int jj, int width,
               double *ua, double *hu, double *sshn_u, int *tmask, double g){
#endif
  int idx = jj*width + ji;

  /*                                  Du                 Dssh
    Flather open boundary condition [---- = sqrt(g/H) * ------]
                                      Dn                 Dn
    ua and va in du/dn should be the specified tidal forcing */

  // Check whether this point lies within the domain
  if(tmask[idx] + tmask[idx+1] <= -1) return;

  if(tmask[idx] < 0){
    ua[idx] = ua[idx+1] +
      sqrt(g/hu[idx]) * (sshn_u[idx] - sshn_u[idx+1]);
  }
  else if(tmask[idx+1]< 0){
    ua[idx] = ua[idx-1] + sqrt(g/hu[idx]) *
     (sshn_u[idx] - sshn_u[idx-1]);
  }
  
}

  /** Kernel to apply Flather boundary condition to v component
      of velocity */
#ifdef __OPENCL_VERSION__
__kernel void bc_flather_v_code(int width,
                __global double* restrict va,
                __global double* restrict hv, 
                __global double* restrict sshn_v, 
                __global int* restrict tmask,
                double g){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  int nrow = (int)get_global_size(1);
  if(jj > (nrow-2))return;
#else
inline void bc_flather_v_code(int ji, int jj, int width,
               double *va, double *hv, double *sshn_v, int *tmask, double g){
#endif
  int idx = jj*width + ji;

  /* Check whether this point is inside the simulated domain
     \todo I could set-up a V-mask using exactly the same code structure
     as below. Could then apply the BC and multiply by V-mask and thus
     remove conditionals => get vectorisation.*/
  if(tmask[idx] + tmask[idx+width] <= -1) return;
    
  if(tmask[idx] < 0){
    va[idx] = va[idx+width] + sqrt(g/hv[idx]) *
     (sshn_v[idx] - sshn_v[idx+width]);
  }
  else if(tmask[idx+width] < 0){
    va[idx] = va[idx-width] + sqrt(g/hv[idx]) *
      (sshn_v[idx] - sshn_v[idx-width]);
  }

}
