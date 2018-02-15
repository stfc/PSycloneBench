
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

#pragma OPENCL EXTENSION cl_intel_channels : enable

__kernel void next_sshu_code(int width,
			     __global double* restrict sshn_u,
			     __global double* restrict sshn,
			     __global int* restrict tmask,
			     __global double* restrict e12t,
			     __global double* restrict e12u){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
  void next_sshu_code(const int ji, const int jj, const int width,
		      double* sshn_u,
		      const double* sshn,
		      const int* tmask,
		      const double* e12t,
		      const double* e12u){
#endif
  double rtmp1;
  int idx = jj*width + ji;
  int idxip1 = idx + 1;
  double ssh = 0.0;
#ifdef INTELFPGA_CL
  // the ssh_channel is declared in continuity_kern_channel.c
  ssh = read_channel_intel(ssh_channel);
#endif
  sshn_u[idx] = 10.0 * ssh;
  /*
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
  */
}
