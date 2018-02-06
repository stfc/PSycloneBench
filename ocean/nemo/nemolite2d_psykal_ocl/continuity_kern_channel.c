
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

#ifdef INTELFPGA_CL
#pragma OPENCL EXTENSION cl_intel_channels : enable
channel double ssh_channel __attribute__((depth(1000)));
#endif

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
  double rtmp1, rtmp2, rtmp3, rtmp4, ssha_idx;
    int idxim1, idxjm1;
    int idx = jj*width + ji;

    //ARPDBG to avoid stalling on a channel if(jj == 0)return;
    /*
    idxim1 = idx - 1;
    idxjm1 = idx - width;

    rtmp1 = (sshn_u[idx] + hu[idx]) * un[idx];
    rtmp2 = (sshn_u[idxim1] + hu[idxim1]) * un[idxim1];
    rtmp3 = (sshn_v[idx] + hv[idx]) * vn[idx];
    rtmp4 = (sshn_v[idxjm1] + hv[idxjm1]) * vn[idxjm1];

    ssha_idx = sshn[idx] + (rtmp2 - rtmp1 + rtmp4 - rtmp3) *
      rdt / e12t[idx];
    */

    // Following line is for testing only
    ssha_idx = (double)idx;

#ifdef __OPENCL_VERSION__
    // Write the value we've just computed to the channel
    // This call is specific to the Intel OpenCL SDK
    write_channel_intel(ssh_channel, ssha_idx);
    //mem_fence(CLK_CHANNEL_MEM_FENCE);
#endif

  }
