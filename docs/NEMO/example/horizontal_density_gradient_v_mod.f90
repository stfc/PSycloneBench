module horizontal_density_gradient_v_mod

  implicit none

  private
  public horizontal_density_gradient_v, horizontal_density_gradient_v_code

  ! metadata
  type, extends(kernel_type) :: horizontal_density_gradient_v
     type(arg), dimension(3) :: meta_args = &
          (/ arg(WRITE, CV),                & ! zgrv
             arg(READ,  CT, STENCIL(N,1)),  & ! prd
             arg(READ,  GRID_MASK_V)        & ! vmask
          /)
     integer :: iterates_over = 2D_LOCAL_POINTS
     integer :: index_offset  = OFFSET_NE
   contains
     procedure, nopass :: code => horizontal_density_gradient_v_code
  end type horizontal_density_gradient_v
	
contains

  subroutine horizontal_density_gradient_v_code(ji,jj,jk,zgrv,prd,vmask)

    integer,  intent(in)                    :: ji,jj,jk
    real(wp), intent(in),  dimension(:,:,:) :: prd, vmask
    real(wp), intent(out), dimension(:,:,:) :: zgrv

    zgrv(ji,jj,jk) = vmask(ji,jj,jk) * ( prd(ji,jj+1,jk) - prd(ji,jj,jk) )

  end subroutine horizontal_density_gradient_v_code

end module horizontal_density_gradient_v_mod
