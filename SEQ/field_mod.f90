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

  !> Base type to capture the information needed about
  !! an external region. This is a region of the field that
  !! is not computed but receives data from elsewhere.
  type :: ExtRegionType
     !> What type of source this external region has
     integer      :: SrcType
     !> Where in the field the data is to be put
     type(region) :: dest     
  end type ExtRegionType

  !> Internal source of data - i.e. data is to be copied
  !! from another part of this same field.
  type, extends(ExtRegionType) :: InternalSourceType
     !> Location of source data within the field
     type(region) :: src
  end type InternalSourceType

  !> Remote source of data - envisioned as another MPI
  !! process.
  type, extends(ExtRegionType) :: RemoteSourceType
     !> MPI rank of process that owns data
     integer      :: iproc
     !> Region of data required from remote process
     type(region) :: src
  end type RemoteSourceType

  !> External region to be populated with data read from file
  type, extends(ExtRegionType) :: FileSourceType
     !> Unit number of file
     integer            :: unit_no
     !> Name of the field in the file to read
     character(len=128) :: field_name
  end type FileSourceType

  type :: grid_type
     !> Total number of grid points
     integer :: npts
  end type grid_type

  type, extends(grid_type) :: rquad_type
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
     !> The number of external regions that this field has
     integer :: NumExtRegions
     !> Array of external regions
     class(ExtRegionType), allocatable, dimension(:) :: ExtRegion
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

!  SUBROUTINE set_2dreal(fld, val)
!    IMPLICIT none
!    TYPE(reg_field_type), INTENT(out) :: fld
!    REAL(wp), INTENT(in) :: val
!
!    fld%data = val

  !===================================================

  subroutine update_external_regions(fld)
    implicit none
    class(field_type), intent(inout) :: fld
    ! Locals
    integer :: iregion

    do iregion=1, fld%NumExtRegions, 1
       write(*,*) 'External region ',iregion,' has type: ',&
                  fld%ExtRegion(iregion)%SrcType

       select case (fld%ExtRegion(iregion)%SrcType)
       case (INTERNAL)
       case (REMOTE)
       case (COUPLER)
       case (FILE)
       case default
       end select
    end do

  end subroutine update_external_regions

end module field_mod
