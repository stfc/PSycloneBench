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

/*
  type, extends(kernel_type) :: next_sshv
     type(arg), dimension(5) :: meta_args =  &
          (/ arg(READWRITE, CV, POINTWISE),  &
             arg(READ,      CV, POINTWISE),  &
             arg(READ,      GRID_MASK_T),    &
             arg(READ,      GRID_AREA_T),    &
             arg(READ,      GRID_AREA_V)     &
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
    procedure, nopass :: code => next_sshv_code
  end type next_sshv
*/

/*
  subroutine invoke_next_sshu(sshn_u, sshn)
    implicit none
    type(r2d_field), intent(inout) :: sshn_u
    type(r2d_field), intent(in)    :: sshn
    ! Locals
    integer :: ji, jj

    do jj = sshn_u%internal%ystart, sshn_u%internal%ystop, 1
      do ji = sshn_u%internal%xstart, sshn_u%internal%xstop, 1

        call next_sshu_code(ji, jj, sshn_u%data, sshn%data, &
                            sshn_u%grid%tmask, &
                            sshn_u%grid%area_t, sshn_u%grid%area_u)
      end do
    end do

  end subroutine invoke_next_sshu
*/

void next_sshu_code(int ji, int jj, int width,
		    double *sshn_u, double *sshn,
		    int *tmask, double *e12t, double *e12u){
  double rtmp1;
  int idx = jj*width + ji;
  int idxip1 = idx + 1;

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

}
    
  /*
  subroutine invoke_next_sshv(sshn_v, sshn)
    implicit none
    type(r2d_field), intent(inout) :: sshn_v
    type(r2d_field), intent(in)    :: sshn
    ! Locals
    integer :: ji, jj

    do jj = sshn_v%internal%ystart, sshn_v%internal%ystop, 1
      do ji = sshn_v%internal%xstart, sshn_v%internal%xstop, 1

        call next_sshv_code(ji, jj,  &
                            sshn_v%data, sshn%data, &
                            sshn_v%grid%tmask,      &
                            sshn_v%grid%area_t, sshn_v%grid%area_v)
      end do
    end do

  end subroutine invoke_next_sshv
  */
  
void next_sshv_code(int ji, int jj, int width,
		    double *sshn_v, double *sshn, int *tmask,
		    double *e12t, double *e12v){
  double rtmp1;
  int idx = jj*width + ji;
  int idxjp1 = idx + width;
    
  if((tmask[idx] + tmask[idxjp1]) <= 0)return; //jump over non-computational domain
  if((tmask[idx] * tmask[idxjp1]) > 0){
    rtmp1 = e12t[idx] * sshn[idx] + e12t[idxjp1] * sshn[idxjp1];
    sshn_v[idx] = 0.5 * rtmp1 / e12v[idx] ;
  }
  else if(tmask[idx] <= 0){
    sshn_v[idx] = sshn[idxjp1];
  }
  else if(tmask[idxjp1] <= 0){
    sshn_v[idx] = sshn[idx];
  }
  
}
