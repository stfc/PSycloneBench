module time_update_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use grid_mod
  use field_mod
  implicit none

  !=======================================

  type, extends(kernel_type) :: next_sshu
     type(arg), dimension(5) :: meta_args =            &
          (/ arg(READWRITE, CU, POINTWISE),            &
             arg(READ,      CU, STENCIL(000,011,000),  &
             arg(READ,      GRID_MASK_T),              &
             arg(READ,      GRID_AREA_T),              &
             arg(READ,      GRID_AREA_U)               &
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

  !=======================================

  type, extends(kernel_type) :: next_sshv
     type(arg), dimension(5) :: meta_args =  &
          (/ arg(READWRITE, CV, POINTWISE),  &
             arg(READ,      CV, STENCIL(010,010,000)),  &
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

contains

  !================================================

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

  !================================================

  subroutine next_sshu_code(ji,jj, sshn_u, sshn, &
                            tmask,e12t,e12u)
    implicit none
    integer,                  intent(in)    :: ji, jj
    integer,  dimension(:,:), intent(in)    :: tmask
    real(wp), dimension(:,:), intent(in)    :: e12t, e12u
    real(wp), dimension(:,:), intent(inout) :: sshn_u
    real(wp), dimension(:,:), intent(in)    :: sshn
    ! Locals
    real(wp) :: rtmp1

    if(tmask(ji,jj) + tmask(ji+1,jj) <= 0)  return   !jump over non-computational domain

    IF(tmask(ji,jj) * tmask(ji+1,jj) > 0) THEN
      rtmp1 = e12t(ji,jj) * sshn(ji,jj) + e12t(ji+1,jj) * sshn(ji+1,jj)
      sshn_u(ji,jj) = 0.5_wp * rtmp1 / e12u(ji,jj) 
    ELSE IF(tmask(ji,jj) <= 0) THEN
      sshn_u(ji,jj) = sshn(ji+1,jj)
    ELSE IF(tmask(ji+1,jj) <= 0) THEN
      sshn_u(ji,jj) = sshn(ji,jj)
    END IF

  end subroutine next_sshu_code
    
  !================================================

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

  !================================================

  subroutine next_sshv_code(ji, jj, &
                            sshn_v, sshn, tmask, e12t, e12v)
    implicit none
    integer,                  intent(in)    :: ji, jj
    integer,  dimension(:,:), intent(in)    :: tmask
    real(wp), dimension(:,:), intent(in)    :: e12t, e12v
    real(wp), dimension(:,:), intent(inout) :: sshn_v
    real(wp), dimension(:,:), intent(in)    :: sshn
    ! Locals
    real(wp) :: rtmp1

    if(tmask(ji,jj) + tmask(ji,jj+1) <= 0)  return !jump over non-computational domain
    if(tmask(ji,jj) * tmask(ji,jj+1) > 0) then
      rtmp1 = e12t(ji,jj) * sshn(ji,jj) + e12t(ji,jj+1) * sshn(ji,jj+1)
      sshn_v(ji,jj) = 0.5_wp * rtmp1 / e12v(ji,jj) 
    else if(tmask(ji,jj) <= 0) then
      sshn_v(ji,jj) = sshn(ji,jj+1)
    else if(tmask(ji,jj+1) <= 0) then
      sshn_v(ji,jj) = sshn(ji,jj)
    end if

  end subroutine next_sshv_code

  !================================================

end module time_update_mod
