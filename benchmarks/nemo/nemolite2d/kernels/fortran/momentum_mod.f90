module momentum_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use grid_mod
  use field_mod
  implicit none

  private

  public invoke_momentum_u, invoke_momentum_v
  public momentum_u, momentum_v
  public momentum_u_code, momentum_v_code

  !=======================================

  type, extends(kernel_type) :: momentum_u
     type(go_arg), dimension(20) :: meta_args =  &
          (/ go_arg(GO_READWRITE, GO_CU, GO_POINTWISE),  & ! ua
             go_arg(GO_READ,      GO_CU, GO_STENCIL(010,111,010)), & ! un
             go_arg(GO_READ,      GO_CV, GO_STENCIL(000,011,011)), & ! vn
             go_arg(GO_READ,      GO_CU, GO_STENCIL(010,010,010)), & ! hu
             go_arg(GO_READ,      GO_CV, GO_STENCIL(000,011,011)), & ! hv
             go_arg(GO_READ,      GO_CT, GO_STENCIL(000,011,000)), & ! ht
             go_arg(GO_READ,      GO_CU, GO_POINTWISE),  & ! ssha_u
             go_arg(GO_READ,      GO_CT, GO_STENCIL(000,011,000)),  & ! sshn_t
             go_arg(GO_READ,      GO_CU, GO_STENCIL(010,010,010)),  & ! sshn_u
             go_arg(GO_READ,      GO_CV, GO_STENCIL(000,011,011)),  & ! sshn_v
              go_arg(GO_READ,  GO_I_SCALAR, GO_POINTWISE),  &
              go_arg(GO_READ,  GO_I_SCALAR, GO_POINTWISE),  &
             go_arg(GO_READ,      GO_GRID_MASK_T),    &
             go_arg(GO_READ,      GO_GRID_DX_U),      &
             go_arg(GO_READ,      GO_GRID_DX_V),      &
             go_arg(GO_READ,      GO_GRID_DX_T),      &
             go_arg(GO_READ,      GO_GRID_DY_U),      &
             go_arg(GO_READ,      GO_GRID_DY_T),      &
             go_arg(GO_READ,      GO_GRID_AREA_U),    &
             go_arg(GO_READ,      GO_GRID_LAT_U)      &
           /)

     !> We update only those points within the internal region
     !! of the simulated domain.
     integer :: ITERATES_OVER = GO_INTERNAL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = GO_OFFSET_NE

  contains
    procedure, nopass :: code => momentum_u_code
  end type momentum_u

  !=======================================

  type, extends(kernel_type) :: momentum_v
     type(go_arg), dimension(20) :: meta_args =  &
          (/ go_arg(GO_READWRITE, GO_CV, GO_POINTWISE),  & ! va
             go_arg(GO_READ,      GO_CU, GO_STENCIL(110,110,000)), & ! un
             go_arg(GO_READ,      GO_CV, GO_STENCIL(010,111,010)), & ! vn
             go_arg(GO_READ,      GO_CU, GO_STENCIL(110,110,000)), & ! hu
             go_arg(GO_READ,      GO_CV, GO_STENCIL(000,111,000)), & ! hv
             go_arg(GO_READ,      GO_CT, GO_STENCIL(010,010,000)), & ! ht
             go_arg(GO_READ,      GO_CV, GO_POINTWISE),  & ! ssha_v
             go_arg(GO_READ,      GO_CT, GO_STENCIL(010,010,000)),  & ! sshn_t
             go_arg(GO_READ,      GO_CU, GO_STENCIL(110,110,010)),  & ! sshn_u
             go_arg(GO_READ,      GO_CV, GO_STENCIL(000,111,000)),  & ! sshn_v
              go_arg(GO_READ,  GO_I_SCALAR, GO_POINTWISE),  &
              go_arg(GO_READ,  GO_I_SCALAR, GO_POINTWISE),  &
             go_arg(GO_READ,      GO_GRID_MASK_T),    &
             go_arg(GO_READ,      GO_GRID_DX_V),      &
             go_arg(GO_READ,      GO_GRID_DX_T),      &
             go_arg(GO_READ,      GO_GRID_DY_U),      &
             go_arg(GO_READ,      GO_GRID_DY_V),      &
             go_arg(GO_READ,      GO_GRID_DY_T),      &
             go_arg(GO_READ,      GO_GRID_AREA_V),    &
             go_arg(GO_READ,      GO_GRID_LAT_V)      &
           /)

     !> We update only those points within the internal region
     !! of the simulated domain.
     integer :: ITERATES_OVER = GO_INTERNAL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = GO_OFFSET_NE

  contains
    procedure, nopass :: code => momentum_v_code
  end type momentum_v

  !=======================================

contains

  !====================================================
 
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
                             32, 2, &
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

  !====================================================

  subroutine momentum_u_code(ji, jj, &
                             ua, un, vn, &
                             hu, hv, ht, ssha_u, &
                             sshn, sshn_u, sshn_v, ct1, ct2, &
                             tmask, e1u, e1v, e1t, e2u, e2t, e12u, gphiu)
    use physical_params_mod, only: omega, d2r, g
    use model_mod, only: rdt, cbfr, visc
    implicit none
    integer, intent(in) :: ji, jj, ct1, ct2
    integer,  dimension(:,:), intent(in) :: tmask
    real(go_wp), dimension(:,:), intent(in) :: e1u, e1v, e1t, e12u, e2u, e2t, gphiu
    real(go_wp), dimension(:,:), intent(in) :: hu, hv, ht
    real(go_wp), dimension(:,:), intent(in) :: ssha_u, sshn, sshn_u, sshn_v
    real(go_wp), dimension(:,:), intent(in) :: un, vn
    real(go_wp), dimension(:,:), intent(inout) :: ua
    ! Locals
    REAL(go_wp) :: u_e, u_w, v_n, v_s
    real(go_wp) :: v_nc, v_sc
    real(go_wp) :: depe, depw, deps, depn
    real(go_wp) :: hpg, adv, cor, vis
    real(go_wp) :: dudx_e, dudx_w
    real(go_wp) :: dudy_s, dudy_n
    real(go_wp) :: uu_e, uu_n, uu_s, uu_w

    integer :: i

    do i = 1, ct1
        ua(ji,jj) = ua(ji,jj) / ct2
    enddo


  end subroutine momentum_u_code
 
  !====================================================
 
  subroutine invoke_momentum_v(va_fld, un_fld, vn_fld, &
                               hu_fld, hv_fld, ht_fld, ssha_v_fld, &
                               sshn_t_fld, sshn_u_fld, sshn_v_fld)
    implicit none
    type(r2d_field), intent(inout) :: va_fld
    type(r2d_field), intent(in) :: un_fld, vn_fld
    type(r2d_field), intent(in) :: hu_fld, hv_fld, ht_fld
    type(r2d_field), intent(in) :: ssha_v_fld, sshn_t_fld, &
                                   sshn_u_fld, sshn_v_fld
    ! Locals
    integer :: ji, jj

    do jj = va_fld%internal%ystart, va_fld%internal%ystop, 1
      do ji = va_fld%internal%xstart, va_fld%internal%xstop, 1

        call momentum_v_code(ji, jj, &
                             va_fld%data, un_fld%data, vn_fld%data, &
                             hu_fld%data, hv_fld%data, ht_fld%data, &
                             ssha_v_fld%data, sshn_t_fld%data,      &
                             sshn_u_fld%data, sshn_v_fld%data,      &
                             32, 2, &
                             va_fld%grid%tmask, va_fld%grid%dx_v,   &
                             va_fld%grid%dx_t, &
                             va_fld%grid%dy_u, va_fld%grid%dy_v,    &
                             va_fld%grid%dy_t,     &
                             va_fld%grid%area_v, va_fld%grid%gphiv)

      end do
   end do

  end subroutine invoke_momentum_v

  !====================================================

  subroutine kernel(ji, jj, array, ct1, ct2)
    implicit none
    integer, intent(in) :: ji, jj, ct1, ct2
    real(go_wp), dimension(:,:), intent(inout) :: array
    integer :: i

    do i = 1, ct1
        va(ji,jj) = va(ji,jj) / ct2
    enddo
  end subroutine kernel

  !====================================================

end module momentum_mod
