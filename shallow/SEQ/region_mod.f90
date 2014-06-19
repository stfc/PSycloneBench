module region_mod
  implicit none

  !> Specify a region on the simulation grid
  type :: region_type
     integer :: xstart, xstop
     integer :: ystart, ystop
  end type region_type

end module region_mod
