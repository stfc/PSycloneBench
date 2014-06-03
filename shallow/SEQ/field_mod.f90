module field_mod
  use kind_params_mod
  use topology_mod, only: region
  implicit none

  ! A field is defined at certain points on the grid, e.g. at U
  ! or V points.
  ! Does that mean a field extends a grid type?
  ! Or that a field carries with it a description of the grid
  ! on which it lives? If the latter then we have to query
  ! the field to see what type of grid it's on.

  TYPE :: grid_type
     !> Total number of grid points
     INTEGER :: npts

  END type grid_type

  TYPE, EXTENDS(grid_type) :: rquad_type
     !> Extent of grid in x
     INTEGER :: nx
     !> Extent of grid in y
     INTEGER :: ny
     !> Grid spacing in x
     REAL(wp) :: dx
     !> Grid spacing in y
     REAL(wp) :: dy

  END TYPE rquad_type

  TYPE, EXTENDS(rquad_type) :: cgrid_type

     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: data

  END type cgrid_type

  TYPE, EXTENDS(cgrid_type) :: tpt_type

  END TYPE tpt_type

  TYPE :: field_type
     TYPE(grid_type) :: grid
  END type field_type

  TYPE, EXTENDS(field_type) :: scalar_field_type
     REAL(wp) :: data
  END TYPE scalar_field_type

  TYPE, EXTENDS(field_type) :: r2d_field_type
     REAL(wp), DIMENSION(:,:), ALLOCATABLE :: data
  END type r2d_field_type

  INTERFACE set
     MODULE PROCEDURE set_scalar_field
  END INTERFACE set

  INTERFACE copy_field
     MODULE PROCEDURE copy_scalar_field, copy_2dfield, copy_2dfield_patch
  END INTERFACE copy_field

  INTERFACE increment
     MODULE PROCEDURE increment_scalar_field
  END INTERFACE increment

CONTAINS

  !===================================================

!  function field_create(mesh, field_type)
!    implicit none
!    type(field_type) :: field_create
!    type(mesh_type)  :: mesh
!    
!  end function field_create

  !===================================================

  SUBROUTINE copy_scalar_field(field_in, field_out)
    IMPLICIT none
    type(scalar_field_type), INTENT(in) :: field_in
    type(scalar_field_type), INTENT(out) :: field_out

    field_out = field_in

  END SUBROUTINE copy_scalar_field

  !===================================================

  SUBROUTINE copy_2dfield(field_in, field_out)
    IMPLICIT none
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: field_in
    REAL(wp), INTENT(out), DIMENSION(:,:) :: field_out
        
    field_out(:,:) = field_in(:,:)
        
  END SUBROUTINE copy_2dfield

  !===================================================

  SUBROUTINE copy_2dfield_patch(field, src, dest)
    IMPLICIT none
    real(wp),     intent(inout), dimension(:,:) :: field
    type(region), intent(in) :: src, dest

    field(dest%istart:dest%istop,dest%jstart:dest%jstop) = &
         field(src%istart:src%istop,src%jstart:src%jstop)
        
  END SUBROUTINE copy_2dfield_patch

  !===================================================

  SUBROUTINE increment_scalar_field(field, incr)
    IMPLICIT none
    TYPE(scalar_field_type), INTENT(inout) :: field
    TYPE(scalar_field_type), INTENT(in)    :: incr

    field%data = field%data + incr%data

  END SUBROUTINE increment_scalar_field

  !===================================================

  SUBROUTINE set_scalar_field(fld, val)
    IMPLICIT none
    TYPE(scalar_field_type), INTENT(out) :: fld
    REAL(wp), INTENT(in) :: val
    fld%data = val
  END SUBROUTINE set_scalar_field

  !===================================================

!  SUBROUTINE set_2dreal(fld, val)
!    IMPLICIT none
!    TYPE(reg_field_type), INTENT(out) :: fld
!    REAL(wp), INTENT(in) :: val
!
!    fld%data = val
END MODULE field_mod
