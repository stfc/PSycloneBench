module horizontal_density_gradient_u_mod

  implicit none

  private
  public horizontal_density_gradient_u, horizontal_density_gradient_u_code

  ! metadata
  type, extends(kernel_type) :: horizontal_density_gradient_u
     type(arg), dimension(3) :: meta_args = &
          (/ arg(WRITE, CU),                & ! zgru
             arg(READ,  CT, STENCIL(E,1)),  & ! prd
             arg(READ,  GRID_MASK_U)        & ! umask
          /)
     integer :: iterates_over = 3D_LOCAL_POINTS
     integer :: index_offset  = OFFSET_NE
   contains
     procedure, nopass :: code => horizontal_density_gradient_u_code
  end type horizontal_density_gradient_u
	
contains

  subroutine horizontal_density_gradient_u_code(ji,jj,jk,zgru,prd,umask)

    integer,  intent(in)                    :: ji,jj,jk
    real(wp), intent(in),  dimension(:,:,:) :: prd, umask
    real(wp), intent(out), dimension(:,:,:) :: zgru

    zgru(ji,jj,jk) = umask(ji,jj,jk) * ( prd(ji+1,jj,jk) - prd(ji,jj,jk) )

  end subroutine horizontal_density_gradient_u_code

end module horizontal_density_gradient_u_mod
