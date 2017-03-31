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


// void invoke_continuity(*ssha, *sshn, *sshn_u, *sshn_v, *hu, *hv, *un, *vn){
//   //use model_mod, only: rdt
//   //type(r2d_field),     intent(in) :: hu, hv, un, vn
//   /* Locals */
//   int ji, jj;

//   for(jj = ssha%internal%ystart, ssha%internal%ystop){
//     for(ji = ssha%internal%xstart, ssha%internal%xstop){

//       continuity_code(ji, jj,
// 		      ssha->data, sshn->data,
// 		      sshn_u->data, sshn_v->data,
// 		      hu->data, hv->data, un->data, vn->data,
// 		      rdt, ssha->grid.area_t);
//     }
//   }

// }

__kernel void continuity_code(int width,                     
			    __global double *ssha,
			    __global double *sshn,
			    __global double *sshn_u,
			    __global double *sshn_v,
			    __global double* hu,
			    __global double *hv,
			    __global double *un,
			    __global double *vn,
			    double rdt,
			    __global double *e12t){
    /* Locals */
    double rtmp1, rtmp2, rtmp3, rtmp4;
    int idx, idxim1, idxjm1;
    int ji = 1, jj = 1;

    idx = jj*width + ji;
    idxim1 = idx - 1;
    idxjm1 = idx - width;
    rtmp1 = (sshn_u[idx] + hu[idx]) * un[idx];
    rtmp2 = (sshn_u[idxim1] + hu[idxim1]) * un[idxim1];
    rtmp3 = (sshn_v[idx] + hv[idx]) * vn[idx];
    rtmp4 = (sshn_v[idxjm1] + hv[idxjm1]) * vn[idxjm1];

    ssha[idx] = sshn[idx] + (rtmp2 - rtmp1 + rtmp4 - rtmp3) *
      rdt / e12t[idx];

  }
