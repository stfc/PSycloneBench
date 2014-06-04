!> Module to contain all mesh-specific data
module mesh_mod
  use kind_params_mod, only: wp
  implicit none

  private

  !> Extents of the grid. Note that the actual arrays used in
  !! the computation have extents one greater than this in
  !! each dimension. Where the different mesh points
  !! sit on this grid is defined in the topology module.
  integer :: nx, ny

  !> Grid spacings in x and y
  real(wp) :: dx, dy
  !> 4.0/dx and 4.0/dy
  real(wp) :: fsdx, fsdy

  public dx, dy, fsdx, fsdy
  public mesh_init

CONTAINS

  subroutine mesh_init(m, n, ldx, ldy)
    implicit none
    integer,  intent(in) :: m, n
    real(wp), intent(in) :: ldx, ldy

    nx = m
    ny = n

    dx = ldx
    dy = ldy

    fsdx = 4./DX
    fsdy = 4./DY

  end subroutine mesh_init

end module mesh_mod
