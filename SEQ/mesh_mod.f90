!> Module to contain all mesh-specific data
MODULE mesh_mod
  IMPLICIT none

  PRIVATE

  !> Extents of the grid. Note that the actual grid used in
  !! the computation has extents one greater than this in
  !! each dimension. Where the different mesh points
  !! sit on this grid is defined in the topology module.
  INTEGER :: nx, ny

  !> Grid spacings in x and y
  REAL(KIND=8) :: dx, dy
  !> 4.0/dx and 4.0/dy
  REAL(KIND=8) :: fsdx, fsdy

  PUBLIC dx, dy, fsdx, fsdy
  PUBLIC mesh_init

CONTAINS

  subroutine mesh_init(m, n, ldx, ldy)
    implicit none
    integer,      intent(in) :: m, n
    real(kind=8), intent(in) :: ldx, ldy

    nx = m
    ny = n

    dx = ldx
    dy = ldy

    fsdx = 4./DX
    fsdy = 4./DY

  END SUBROUTINE mesh_init

END MODULE mesh_mod
