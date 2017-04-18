#include <stdio.h>

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
__kernel void momentum_u_code(__global double *ua,
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
			      __global double *gphiu){
#else
/** Interface to standard C version of kernel */
void momentum_u_code(int ji, int jj,
		     double *ua, double *un, double *vn,
		     double *hu, double *hv, double *ht, double *ssha_u,
		     double *sshn, double *sshn_u, double *sshn_v,
		     int *tmask,
		     double *e1u, double *e1v, double *e1t,
		     double *e2u, double *e2t, double *e12u, double *gphiu){
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

  if(tmask(ji,jj) + tmask(ji+1,jj) <= 0) return; // jump over non-computational domain
  if(tmask(ji,jj) <= 0 || tmask(ji+1,jj) <= 0) return; // jump over boundary u

  u_e  = 0.5 * (un(ji,jj) + un(ji+1,jj)) * e2t(ji+1,jj);   // add length scale.
  depe = ht(ji+1,jj) + sshn(ji+1,jj);

  u_w  = 0.5 * (un(ji,jj) + un(ji-1,jj)) * e2t(ji,jj);   // add length scale
  depw = ht(ji,jj) + sshn(ji,jj);

  v_sc = 0.5 * (vn(ji,jj-1) + vn(ji+1,jj-1));
  v_s  = 0.5 * v_sc * (e1v(ji,jj-1) + e1v(ji+1,jj-1));
  deps = 0.5 * (hv(ji,jj-1) + sshn_v(ji,jj-1) + hv(ji+1,jj-1) +
		sshn_v(ji+1,jj-1));

  v_nc = 0.5 * (vn(ji,jj) + vn(ji+1,jj));
  v_n  = 0.5 * v_nc * (e1v(ji,jj) + e1v(ji+1,jj));
  depn = 0.5 * (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) +
		sshn_v(ji+1,jj));

  // -advection (currently first order upwind)
  uu_w = (0.5 - SIGN(0.5, u_w)) * un(ji,jj) +
    (0.5 + SIGN(0.5, u_w)) * un(ji-1,jj) ;
  uu_e = (0.5 + SIGN(0.5, u_e)) * un(ji,jj) +
    (0.5 - SIGN(0.5, u_e)) * un(ji+1,jj) ;

  if(tmask(ji,jj-1) <=0 || tmask(ji+1,jj-1) <= 0){
    uu_s = (0.5_wp - SIGN(0.5, v_s)) * un(ji,jj);
  }
  else{
    uu_s = (0.5 - SIGN(0.5, v_s)) * un(ji,jj) +
      (0.5 + SIGN(0.5, v_s)) * un(ji,jj-1) ;
  }

  if(tmask(ji,jj+1) <=0 || tmask(ji+1,jj+1) <= 0){
    uu_n = (0.5 + SIGN(0.5, v_n)) * un(ji,jj);
  }
  else{
    uu_n = (0.5 + SIGN(0.5, v_n)) * un(ji,jj) +
      (0.5 - SIGN(0.5, v_n)) * un(ji,jj+1);
  }

  adv = uu_w * u_w * depw - uu_e * u_e * depe +
    uu_s * v_s * deps - uu_n * v_n * depn;

  // -viscosity

  dudx_e = (un(ji+1,jj) - un(ji,  jj)) / e1t(ji+1,jj) * 
    (ht(ji+1,jj) + sshn(ji+1,jj));
  dudx_w = (un(ji,  jj) - un(ji-1,jj)) / e1t(ji,  jj) * 
    (ht(ji,  jj) + sshn(ji,  jj));
  if(tmask(ji,jj-1) <=0 || tmask(ji+1,jj-1) <= 0){
    dudy_s = 0.0; // slip boundary
  }
  else{
    dudy_s = (un(ji,jj) - un(ji,jj-1)) / (e2u(ji,jj) + e2u(ji,jj-1)) * 
      (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj-1) + sshn_u(ji,jj-1));
  }

  if(tmask(ji,jj+1) <= 0 || tmask(ji+1,jj+1) <= 0){
    dudy_n = 0.0; // slip boundary
  }
  else{
    dudy_n = (un(ji,jj+1) - un(ji,jj)) / (e2u(ji,jj) + e2u(ji,jj+1)) * 
      (hu(ji,jj) + sshn_u(ji,jj) + hu(ji,jj+1) + sshn_u(ji,jj+1));
  }

  vis = (dudx_e - dudx_w ) * e2u(ji,jj)  + 
    (dudy_n - dudy_s ) * e1u(ji,jj) * 0.5;
  vis = visc * vis;   // visc will be an array visc(1:jpijglou) 
  // for variable viscosity, such as turbulent viscosity

  // -Coriolis' force (can be implemented implicitly)
  cor = 0.5 * (2. * omega * SIN(gphiu(ji,jj) * d2r) * (v_sc + v_nc)) * 
    e12u(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj));

  // -pressure gradient
  hpg = -g * (hu(ji,jj) + sshn_u(ji,jj)) * e2u(ji,jj) * 
    (sshn(ji+1,jj) - sshn(ji,jj));

  // -linear bottom friction (implemented implicitly.
  ua(ji,jj) = (un(ji,jj) * (hu(ji,jj) + sshn_u(ji,jj)) + rdt * 
                 (adv + vis + cor + hpg) / e12u(ji,jj)) / 
    (hu(ji,jj) + ssha_u(ji,jj)) / (1.0 + cbfr * rdt) ;

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
__kernel void momentum_v_code(__global double *va,
			      __global double *un,
			      __global double *vn, 
			      __global double *hu,
			      __global double *hv,
			      __global double *ht,
			      __global double *ssha_v, 
			      __global double *sshn,
			      __global double *sshn_u,
			      __global double *sshn_v, 
			      __global double *tmask,
			      __global double *e1v,
			      __global double *e1t,
			      __global double *e2u,
			      __global double *e2v,
			      __global double *e2t,
			      __global double *e12v,
			      __global double *gphiv){
#else
/** Interface to standard C version of kernel */
void momentum_v_code(int ji, int jj, 
		     double *va, double *un, double *vn, 
		     double *hu, double *hv, double *ht, double *ssha_v, 
		     double *sshn, double *sshn_u, double *sshn_v, 
		     int *tmask, double *e1v, double *e1t, double *e2u,
		     double *e2v, double *e2t, double *e12v, double *gphiv){
#endif

  //use physical_params_mod
  //use model_mod, only: rdt, cbfr, visc
  double u_e, u_w, v_n, v_s;
  double u_ec, u_wc, vv_e, vv_n, vv_s, vv_w;
  double depe, depw, deps, depn;
  double hpg, adv, cor, vis;
  double dvdx_e, dvdx_w, dvdy_n, dvdy_s;

  if(tmask(ji,jj) + tmask(ji+1,jj) <= 0)  return; // jump over non-computatinal domain
  if(tmask(ji,jj) <= 0 || tmask(ji,jj+1) <= 0) return; // jump over v boundary cells

  // kernel v adv 
  v_n  = 0.5 * (vn(ji,jj) + vn(ji,jj+1)) * e1t(ji,jj+1); // add length scale.
  depn = ht(ji,jj+1) + sshn(ji,jj+1);

  v_s  = 0.5 * (vn(ji,jj) + vn(ji,jj-1)) * e1t(ji,jj); // add length scale
  deps = ht(ji,jj) + sshn(ji,jj);

  u_wc = 0.5 * (un(ji-1,jj) + un(ji-1,jj+1));
  u_w  = 0.5 * u_wc * (e2u(ji-1,jj) + e2u(ji-1,jj+1));
  depw = 0.50 * (hu(ji-1,jj) + sshn_u(ji-1,jj) + 
		 hu(ji-1,jj+1) + sshn_u(ji-1,jj+1));

  u_ec = 0.5 * (un(ji,jj) + un(ji,jj+1));
  u_e  = 0.5 * u_ec * (e2u(ji,jj) + e2u(ji,jj+1));
  depe = 0.50 * (hu(ji,jj) + sshn_u(ji,jj) + 
		 hu(ji,jj+1) + sshn_u(ji,jj+1));

  // -advection (currently first order upwind)
  vv_s = (0.5 - SIGN(0.5, v_s)) * vn(ji,jj)              +  
    (0.5 + SIGN(0.5, v_s)) * vn(ji,jj-1) ;
  vv_n = (0.5 + SIGN(0.5, v_n)) * vn(ji,jj)              +  
    (0.5 - SIGN(0.5, v_n)) * vn(ji,jj+1) ;

  if(tmask(ji-1,jj) <= 0 || tmask(ji-1,jj+1) <= 0){
    vv_w = (0.5 - SIGN(0.5, u_w)) * vn(ji,jj);
  }
  else{
    vv_w = (0.5 - SIGN(0.5, u_w)) * vn(ji,jj) +  
      (0.5 + SIGN(0.5, u_w)) * vn(ji-1,jj) ;
  }

  if(tmask(ji+1,jj) <= 0 || tmask(ji+1,jj+1) <= 0){
    vv_e = (0.5 + SIGN(0.5, u_e)) * vn(ji,jj);
  }
  else{
    vv_e = (0.5 + SIGN(0.5, u_e)) * vn(ji,jj) +
      (0.5 - SIGN(0.5, u_e)) * vn(ji+1,jj);
  }

  adv = vv_w * u_w * depw - vv_e * u_e * depe + 
    vv_s * v_s * deps - vv_n * v_n * depn;

  // -viscosity

    
  dvdy_n = (vn(ji,jj+1) - vn(ji,  jj)) / e2t(ji,jj+1) * 
    (ht(ji,jj+1) + sshn(ji,jj+1));
  dvdy_s = (vn(ji,  jj) - vn(ji,jj-1)) / e2t(ji,  jj) * 
    (ht(ji,  jj) + sshn(ji,  jj));

  if(tmask(ji-1,jj) <= 0 || tmask(ji-1,jj+1) <= 0){
    dvdx_w = 0.0; // slip boundary
  }
  else{
    dvdx_w = (vn(ji,jj) - vn(ji-1,jj)) / (e1v(ji,jj) + e1v(ji-1,jj)) * 
      (hv(ji,jj) + sshn_v(ji,jj) + hv(ji-1,jj) + sshn_v(ji-1,jj));
  }

  if(tmask(ji+1,jj) <= 0 || tmask(ji+1,jj+1) <= 0){
    dvdx_e = 0.0; // slip boundary
  }
  else{
    dvdx_e = (vn(ji+1,jj) - vn(ji,jj)) / (e1v(ji,jj) + e1v(ji+1,jj)) * 
      (hv(ji,jj) + sshn_v(ji,jj) + hv(ji+1,jj) + sshn_v(ji+1,jj));
  }

  vis = (dvdy_n - dvdy_s ) * e1v(ji,jj)  + 
    (dvdx_e - dvdx_w ) * e2v(ji,jj) * 0.5  ;

  vis = visc * vis;   // visc will be a array visc(1:jpijglou) 
  // for variable viscosity, such as turbulent viscosity

  // -Coriolis' force (can be implemented implicitly)
  cor = -0.5*(2. * omega * SIN(gphiv(ji,jj) * d2r) * (u_ec + u_wc)) * 
    e12v(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj));

  // -pressure gradient
  hpg = -g * (hv(ji,jj) + sshn_v(ji,jj)) * e1v(ji,jj) * 
    (sshn(ji,jj+1) - sshn(ji,jj));

  // -linear bottom friction (implemented implicitly.
  va(ji,jj) = (vn(ji,jj) * (hv(ji,jj) + sshn_v(ji,jj)) + 
	       rdt * (adv + vis + cor + hpg) / e12v(ji,jj) ) / 
    ((hv(ji,jj) + ssha_v(ji,jj))) / (1.0 + cbfr * rdt) ;

}
