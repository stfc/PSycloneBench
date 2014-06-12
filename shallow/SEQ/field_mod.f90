module field_mod
  use kind_params_mod
  use topology_mod
  implicit none

  ! A field is defined at certain points on the grid, e.g. at U
  ! or V points.
  ! Does that mean a field extends a grid type?
  ! Or that a field carries with it a description of the grid
  ! on which it lives? If the latter then we have to query
  ! the field to see what type of grid it's on.

  !> Labels for the different types of 'external' data source
  integer, parameter :: INTERNAL= 0
  integer, parameter :: REMOTE  = 1
  integer, parameter :: COUPLER = 2
  integer, parameter :: FILE    = 3

  ! Apart from MPI/thread halos, the application of external
  ! boundary conditions is algorithm-driven.
  ! e.g. the use of periodic boundary conditions is a key
  ! physical characteristic of the model. Therefore, I'm not
  ! sure that we need to capture the vagaries of where data
  ! might come from.

  type, extends(topology_type) :: rquad_type
     !> Total number of grid points
     integer :: npts
     !> Extent of grid in x
     INTEGER :: nx
     !> Extent of grid in y
     INTEGER :: ny
     !> Grid spacing in x
     REAL(wp) :: dx
     !> Grid spacing in y
     REAL(wp) :: dy
  end type rquad_type

  type :: field_type
      !> The grid on which this field is defined
     type(rquad_type) :: grid
     !> The internal region of this field
     type(region) :: IntRegion
  end type field_type

  type, extends(field_type) :: scalar_field_type
     real(wp) :: data
  end TYPE scalar_field_type

  type, extends(field_type) :: r2d_field_type
     !> Array holding the actual field values
     real(wp), dimension(:,:), allocatable :: data
  end type r2d_field_type

  !> A field defined on the U points of the Arakawa C Grid
  type, extends(r2d_field_type) :: cu_field_type

  end type cu_field_type

  type, extends(r2d_field_type) :: cv_field_type

  end type cv_field_type

  type, extends(r2d_field_type) :: ct_field_type

  end type ct_field_type

  type, extends(r2d_field_type) :: cf_field_type

  end type cf_field_type

  interface set
     module procedure set_scalar_field
  end interface set

  interface copy_field
     module procedure copy_scalar_field, copy_2dfield, copy_2dfield_patch
  end interface copy_field

  interface increment
     module procedure increment_scalar_field
  end interface increment

contains

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
    class(field_type), INTENT(out) :: fld
    REAL(wp), INTENT(in) :: val

    select type(fld)
    type is (scalar_field_type)
       fld%data = val
    type is (r2d_field_type)
       fld%data = val
    class default
    end select

  END SUBROUTINE set_scalar_field

  !===================================================

end module field_mod
