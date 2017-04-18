#ifndef __OPENCL_VERSION__
// This header isn't available in OpenCL
#include <math.h>
#endif

#define PI 3.1415926535897932
// Acceleration due to gravity (ms^-2)
#define G 9.80665
// Earth rotation speed (s^(-1))
#define OMEGA 7.292116e-05
// Degree to radian
#define d2r PI/180.0

// Although C99 has copysign(), this isn't necessarily available in
// OpenCL
//double sign(double a, double b){
//  if(b < 0.0)return -1.0*a;
//  return a;
//}

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


/*
  type, extends(kernel_type) :: momentum_v
     type(arg), dimension(18) :: meta_args =  &
          (/ arg(READWRITE, CV, POINTWISE),  & ! va
             arg(READ,      CU, POINTWISE),  & ! un
             arg(READ,      CV, POINTWISE),  & ! vn
             arg(READ,      CU, POINTWISE),  & ! hu
             arg(READ,      CV, POINTWISE),  & ! hv
             arg(READ,      CT, POINTWISE),  & ! ht
             arg(READ,      CV, POINTWISE),  & ! ssha_v
             arg(READ,      CT, POINTWISE),  & ! sshn_t
             arg(READ,      CU, POINTWISE),  & ! sshn_u
             arg(READ,      CV, POINTWISE),  & ! sshn_v
             arg(READ,      GRID_MASK_T),    &
             arg(READ,      GRID_DX_V),      &
             arg(READ,      GRID_DX_T),      &
             arg(READ,      GRID_DY_U),      &
             arg(READ,      GRID_DY_V),      &
             arg(READ,      GRID_DY_T),      &
             arg(READ,      GRID_AREA_V),    &
             arg(READ,      GRID_LAT_V)      &
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
    procedure, nopass :: code => momentum_v_code
  end type momentum_v
*/

/*
  subroutine invoke_momentum_u(ua_fld, un_fld, vn_fld, &
                               hu_fld, hv_fld, ht_fld, ssha_u_fld, &
                               sshn_t_fld, sshn_u_fld, sshn_v_fld)
    implicit none
    type(r2d_field), intent(inout) :: ua_fld
    type(r2d_field), intent(in) :: un_fld, vn_fld
    type(r2d_field), intent(in) :: hu_fld, hv_fld, ht_fld
    type(r2d_field), intent(in) :: ssha_u_fld, sshn_t_fld, &
                                        sshn_u_fld, sshn_v_fld
    ! Locals
    integer :: ji, jj

    do jj = ua_fld%internal%ystart, ua_fld%internal%ystop, 1
      do ji = ua_fld%internal%xstart, ua_fld%internal%xstop, 1

        call momentum_u_code(ji, jj, &
                             ua_fld%data, un_fld%data, vn_fld%data, &
                             hu_fld%data, hv_fld%data, ht_fld%data, &
                             ssha_u_fld%data, sshn_t_fld%data,      &
                             sshn_u_fld%data, sshn_v_fld%data, &
                             ua_fld%grid%tmask,  &
                             ua_fld%grid%dx_u,   &
                             ua_fld%grid%dx_v,   &
                             ua_fld%grid%dx_t,   &
                             ua_fld%grid%dy_u,   &
                             ua_fld%grid%dy_t,   &
                             ua_fld%grid%area_u, &
                             ua_fld%grid%gphiu)
      end do
   end do

  end subroutine invoke_momentum_u

*/

#ifdef __OPENCL_VERSION__
/** Interface to OpenCL version of kernel */
__kernel void momentum_u_code(int width,
			      __global double *ua,
			      __global double *un,
			      __global double *vn,
			      __global double *hu,
			      __global double *hv,
			      __global double *ht,
			      __global double *ssha_u,
			      __global double *sshn,
			      __global double *sshn_u,
			      __global double *sshn_v,
			      __global int *tmask,
			      __global double *e1u,
			      __global double *e1v,
			      __global double *e1t,
			      __global double *e2u,
			      __global double *e2t,
			      __global double *e12u,
			      __global double *gphiu,
			      double rdt, double cbfr, double visc){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
/** Interface to standard C version of kernel */
void momentum_u_code(int ji, int jj, int width,
		     double *ua, double *un, double *vn,
		     double *hu, double *hv, double *ht, double *ssha_u,
		     double *sshn, double *sshn_u, double *sshn_v,
		     int *tmask,
		     double *e1u, double *e1v, double *e1t,
		     double *e2u, double *e2t, double *e12u, double *gphiu,
		     double rdt, double cbfr, double visc){
#endif
  //    use physical_params_mod
  //  use model_mod, only: rdt, cbfr, visc
  double u_e, u_w, v_n, v_s;
  double v_nc, v_sc;
  double depe, depw, deps, depn;
  double hpg, adv, cor, vis;
  double dudx_e, dudx_w;
  double dudy_s, dudy_n;
  double uu_e, uu_n, uu_s, uu_w;
  
  int idxim1, idxjm1, idxip1, idxjp1, idxip1jm1;
  int idx = jj*width + ji;

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
  cor = 0.5 * (2. * OMEGA * sin(gphiu[idx] * d2r) * (v_sc + v_nc)) * 
    e12u[idx] * (hu[idx] + sshn_u[idx]);

  // -pressure gradient
  hpg = -G * (hu[idx] + sshn_u[idx]) * e2u[idx] * 
    (sshn[idxip1] - sshn[idx]);

  // -linear bottom friction (implemented implicitly.
  ua[idx] = (un[idx] * (hu[idx] + sshn_u[idx]) + rdt * 
                 (adv + vis + cor + hpg) / e12u[idx]) / 
    (hu[idx] + ssha_u[idx]) / (1.0 + cbfr * rdt) ;

}
 
/*
  subroutine invoke_momentum_v(va_fld, un_fld, vn_fld, 
                               hu_fld, hv_fld, ht_fld, ssha_v_fld, 
                               sshn_t_fld, sshn_u_fld, sshn_v_fld)
    implicit none
    type(r2d_field), intent(inout) :: va_fld
    type(r2d_field), intent(in) :: un_fld, vn_fld
    type(r2d_field), intent(in) :: hu_fld, hv_fld, ht_fld
    type(r2d_field), intent(in) :: ssha_v_fld, sshn_t_fld, 
                                   sshn_u_fld, sshn_v_fld
    ! Locals
    integer :: ji, jj

    do jj = va_fld%internal%ystart, va_fld%internal%ystop, 1
      do ji = va_fld%internal%xstart, va_fld%internal%xstop, 1

        call momentum_v_code(ji, jj, 
                             va_fld%data, un_fld%data, vn_fld%data, 
                             hu_fld%data, hv_fld%data, ht_fld%data, 
                             ssha_v_fld%data, sshn_t_fld%data,      
                             sshn_u_fld%data, sshn_v_fld%data,      
                             va_fld%grid%tmask, va_fld%grid%dx_v,   
                             va_fld%grid%dx_t, 
                             va_fld%grid%dy_u, va_fld%grid%dy_v,    
                             va_fld%grid%dy_t,     
                             va_fld%grid%area_v, va_fld%grid%gphiv)

      end do
   end do

  end subroutine invoke_momentum_v
*/

#ifdef __OPENCL_VERSION__
/** Interface to OpenCL version of kernel */
__kernel void momentum_v_code(int width,
			      __global double *va,
			      __global double *un,
			      __global double *vn, 
			      __global double *hu,
			      __global double *hv,
			      __global double *ht,
			      __global double *ssha_v, 
			      __global double *sshn,
			      __global double *sshn_u,
			      __global double *sshn_v, 
			      __global int *tmask,
			      __global double *e1v,
			      __global double *e1t,
			      __global double *e2u,
			      __global double *e2v,
			      __global double *e2t,
			      __global double *e12v,
			      __global double *gphiv,
			      double rdt, double cbfr, double visc){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
/** Interface to standard C version of kernel */
void momentum_v_code(int ji, int jj, int width,
		     double *va, double *un, double *vn, 
		     double *hu, double *hv, double *ht, double *ssha_v, 
		     double *sshn, double *sshn_u, double *sshn_v, 
		     int *tmask, double *e1v, double *e1t, double *e2u,
		     double *e2v, double *e2t, double *e12v, double *gphiv,
		     double rdt, double cbfr, double visc){
#endif

  //use physical_params_mod
  //use model_mod, only: rdt, cbfr, visc
  double u_e, u_w, v_n, v_s;
  double u_ec, u_wc, vv_e, vv_n, vv_s, vv_w;
  double depe, depw, deps, depn;
  double hpg, adv, cor, vis;
  double dvdx_e, dvdx_w, dvdy_n, dvdy_s;
  
  int idxim1, idxjm1, idxip1, idxjp1, idxip1jm1, idxim1jp1;
  int idx = jj*width + ji;

  idxim1 = idx - 1;
  idxip1 = idx + 1;
  idxjm1 = idx - width;
  idxjp1 = idx + width;
  idxip1jm1 = idx - width + 1;
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
  cor = -0.5*(2. * OMEGA * sin(gphiv[idx] * d2r) * (u_ec + u_wc)) * 
    e12v[idx] * (hv[idx] + sshn_v[idx]);

  // -pressure gradient
  hpg = -G * (hv[idx] + sshn_v[idx]) * e1v[idx] * 
    (sshn[idxjp1] - sshn[idx]);

  // -linear bottom friction (implemented implicitly.
  va[idx] = (vn[idx] * (hv[idx] + sshn_v[idx]) + 
	       rdt * (adv + vis + cor + hpg) / e12v[idx] ) / 
    ((hv[idx] + ssha_v[idx])) / (1.0 + cbfr * rdt) ;

}
