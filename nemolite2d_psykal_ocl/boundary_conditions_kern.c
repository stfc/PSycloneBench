#ifndef __OPENCL_VERSION__
// This header isn't available/required in OpenCL
#include <math.h>
#endif
#include "physical_params.h"
/*

  type, extends(kernel_type) :: bc_ssh
     type(arg), dimension(3) :: meta_args =        &
          (/ arg(READ,      I_SCALAR, POINTWISE),  &
             arg(READWRITE, CT,       POINTWISE),  &
             arg(READ,      GRID_MASK_T)           &
           /)

     !> Although this is a boundary-conditions kernel, it only
     !! acts on the internal points of the domain
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
    procedure, nopass :: code => bc_ssh_code
  end type bc_ssh

*/

/*

  type, extends(kernel_type) :: bc_solid_u
     type(arg), dimension(2) :: meta_args =  &
          (/ arg(WRITE, CU, POINTWISE),      &
             arg(READ,      GRID_MASK_T)     &
           /)

     !> This is a boundary-conditions kernel and therefore
     !! acts on all points of the domain rather than just
     !! those that are internal
     integer :: ITERATES_OVER = ALL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => bc_solid_u_code
  end type bc_solid_u
*/

/*
  type, extends(kernel_type) :: bc_solid_v
     type(arg), dimension(2) :: meta_args =  &
          (/ arg(WRITE, CV, POINTWISE),      &
             arg(READ,      GRID_MASK_T)     &
           /)

     !> This is a boundary-conditions kernel and therefore
     !! acts on all points of the domain rather than just
     !! those that are internal
     integer :: ITERATES_OVER = ALL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => bc_solid_v_code
  end type bc_solid_v
*/

/*
  type, extends(kernel_type) :: bc_flather_u
     type(arg), dimension(4) :: meta_args =  &
          (/ arg(READWRITE, CU, POINTWISE),  & ! ua
             arg(READ,      CU, POINTWISE),  & ! hu
             arg(READ,      CU, POINTWISE),  & ! sshn_u
             arg(READ,      GRID_MASK_T)     &
           /)

     !> This is a boundary-conditions kernel and therefore
     !! acts on all points of the domain rather than just
     !! those that are internal
     integer :: ITERATES_OVER = ALL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => bc_flather_u_code
  end type bc_flather_u
*/

/*
  type, extends(kernel_type) :: bc_flather_v
     type(arg), dimension(4) :: meta_args =  &
          (/ arg(READWRITE, CV, POINTWISE),  & ! va
             arg(READ,      CV, POINTWISE),  & ! hv
             arg(READ,      CV, POINTWISE),  & ! sshn_v
             arg(READ,      GRID_MASK_T)     &
           /)

     !> This is a boundary-conditions kernel and therefore
     !! acts on all points of the domain rather than just
     !! those that are internal
     integer :: ITERATES_OVER = ALL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => bc_flather_v_code
  end type bc_flather_v
*/

/*
  subroutine invoke_bc_ssh(istep, ssha)
    implicit none
    integer,            intent(in)    :: istep
    type(r2d_field),    intent(inout) :: ssha
    ! Locals
    integer  :: ji, jj

    DO jj = ssha%internal%ystart, ssha%internal%ystop
       DO ji = ssha%internal%xstart, ssha%internal%xstop
          call bc_ssh_code(ji, jj, &
                           istep, ssha%data, ssha%grid%tmask)
       END DO
    END DO

  end subroutine invoke_bc_ssh
*/

#ifdef __OPENCL_VERSION__
__kernel void bc_ssh_code(int width,
			  int istep,
			  __global double *ssha,
			  __global int *tmask,
			  double rdt){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
void bc_ssh_code(int ji, int jj, int width,
		 int istep, double *ssha, int *tmask, double rdt){
#endif
  int idx = jj*width + ji;

  double amp_tide, omega_tide, rtime;

  amp_tide   = 0.2;
  omega_tide = 2.0 * 3.14159 / (12.42 * 3600.0);
  rtime = (double)istep * rdt;

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
  
  /*
  !> Manual version of code to invoke kernel that applies solid 
  !! boundary conditions for u-velocity
  subroutine invoke_bc_solid_u(ua)
    implicit none
    type(r2d_field), intent(inout) :: ua
    ! Locals
    integer  :: ji, jj

! Original loop was:
!            DO jj = 1, jpj
!              DO ji = 0, jpi
! In original code, tmask is declared with one more row and column than
! any other field. ji==jpi IS last column of u field.
! 1/ How do I determine the full range of array indices to loop over for ua?
! 2/ If I do that, is tmask(ji+1,jj) going to stay within bounds?
    do jj = ua%whole%ystart, ua%whole%ystop, 1
       do ji = ua%whole%xstart, ua%whole%xstop, 1
          call bc_solid_u_code(ji, jj, ua%data, ua%grid%tmask)
       end do
    end do

  end subroutine invoke_bc_solid_u
*/

    /** Kernel to apply solid boundary conditions for u-velocity */
#ifdef __OPENCL_VERSION__
__kernel void bc_solid_u_code(int width,
			      __global double *ua,
			      __global int *tmask){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
  void bc_solid_u_code(int ji, int jj, int width, double *ua, int *tmask){
#endif
  int idx = jj*width + ji;

  if(tmask[idx] * tmask[idx+1] == 0){
    ua[idx] = 0.0;
  }

}
  
  /*
  !> Manual version of code to invoke the kernel
  !! that applies the solid-bc to a field on V pts.
  subroutine invoke_bc_solid_v(va)
    implicit none
    type(r2d_field), intent(inout) :: va
    ! Locals
    integer  :: ji, jj

    do jj = va%whole%ystart, va%whole%ystop, 1
       do ji = va%whole%xstart, va%whole%xstop, 1
          call bc_solid_v_code(ji,jj,va%data,va%grid%tmask)
      end do
    end do

  end subroutine invoke_bc_solid_v
*/
  
  /** Kernel to apply solid boundary conditions for v-velocity */
#ifdef __OPENCL_VERSION__
__kernel void bc_solid_v_code(int width,
			      __global double *va, __global int *tmask){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
void bc_solid_v_code(int ji, int jj, int width, double *va, int *tmask){
#endif
  int idx = jj*width + ji;

  if(tmask[idx] * tmask[idx+width] == 0){
    va[idx] = 0.0;
  }

}
  
  /*
  !>                                  Du                 Dssh
  !!Flather open boundary condition [---- = sqrt(g/H) * ------]
  !!                                  Dn                 Dn
  !! ua and va in du/dn should be the specified tidal forcing
  subroutine invoke_bc_flather_u(ua, hu, sshn_u)
    implicit none
    type(r2d_field), intent(inout) :: ua
    type(r2d_field), intent(in) :: hu, sshn_u
    ! Locals
    integer  :: ji, jj

    ! Original loop was:
    !            DO jj = 1, jpj
    !              DO ji = 0, jpi  
    DO jj = ua%whole%ystart, ua%whole%ystop, 1
       DO ji = ua%whole%xstart, ua%whole%xstop, 1
          call bc_flather_u_code(ji,jj, &
                                 ua%data, hu%data, sshn_u%data, &
                                 ua%grid%tmask)
       END DO
    END DO
  
  end subroutine invoke_bc_flather_u
*/  

/** Kernel to apply Flather condition to U */
#ifdef __OPENCL_VERSION__
__kernel void bc_flather_u_code(int width,
				__global double *ua,
				__global double *hu,
				__global double *sshn_u,
				__global int *tmask){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
void bc_flather_u_code(int ji, int jj, int width,
		       double *ua, double *hu, double *sshn_u, int *tmask){
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
      sqrt(G/hu[idx]) * (sshn_u[idx] - sshn_u[idx+1]);
  }
  else if(tmask[idx+1]< 0){
    ua[idx] = ua[idx-1] + sqrt(G/hu[idx]) *
	 (sshn_u[idx] - sshn_u[idx-1]);
  }
  
}

  /*
  !> Manual version of code to invoke the kernel for applying the
  !! Flather boundary condition to the v component of velocity.
  subroutine invoke_bc_flather_v(va, hv, sshn_v)
    implicit none
    type(r2d_field), intent(inout) :: va
    type(r2d_field), intent(in)    :: hv, sshn_v
    ! Locals
    integer  :: ji, jj

    !kernel Flather v 

    DO jj = va%whole%ystart, va%whole%ystop, 1
       DO ji = va%whole%xstart, va%whole%xstop, 1
          call bc_flather_v_code(ji,jj, &
                                 va%data, hv%data, sshn_v%data, &
                                 va%grid%tmask)
       END DO
    END DO

  end subroutine invoke_bc_flather_v
*/

  /** Kernel to apply Flather boundary condition to v component
      of velocity */
#ifdef __OPENCL_VERSION__
__kernel void bc_flather_v_code(int width,
				__global double *va,
				__global double *hv, 
				__global double *sshn_v, 
				__global int *tmask){
  int ji = get_global_id(0);
  int jj = get_global_id(1);
#else
void bc_flather_v_code(int ji, int jj, int width,
		       double *va, double *hv, double *sshn_v, int *tmask){
#endif
  int idx = jj*width + ji;

  /* Check whether this point is inside the simulated domain
     \todo I could set-up a V-mask using exactly the same code structure
     as below. Could then apply the BC and multiply by V-mask and thus
     remove conditionals => get vectorisation.*/
  if(tmask[idx] + tmask[idx+width] <= -1) return;
    
  if(tmask[idx] < 0){
    va[idx] = va[idx+width] + sqrt(G/hv[idx]) *
	 (sshn_v[idx] - sshn_v[idx+width]);
  }
  else if(tmask[idx+width] < 0){
    va[idx] = va[idx-width] + sqrt(G/hv[idx]) *
      (sshn_v[idx] - sshn_v[idx-width]);
  }

}

/*
  void setup_vmask_code(ji, jj, vmask, tmask){
    integer, intent(in) :: ji, jj
    integer,  dimension(:,:), intent(inout) :: vmask
    integer,  dimension(:,:), intent(in) :: tmask
    integer :: jiv

       vmask[idx] = 0;

if(tmask[idx] + tmask(ji,jj+1) <= -1) return;
    
    IF(tmask[idx] < 0) THEN
  jiv = jj + 1;
    ELSE IF(tmask(ji,jj+1) < 0) THEN
    jiv = jj - 1 ;
    END IF

    vmask(ji,jiv) = 1;

	  }
*/
