module field_mod
  use kind_params_mod
  use region_mod
  use halo_mod
  use grid_mod
  implicit none

  integer, parameter :: U_POINTS = 0
  integer, parameter :: V_POINTS = 1
  integer, parameter :: T_POINTS = 2
  integer, parameter :: F_POINTS = 3

  type :: field_type
     !> Which mesh points is the field defined on
     integer :: defined_on
     !> The grid on which this field is defined
     type(grid_type), pointer :: grid
     !> The internal region of this field
     type(region_type) :: internal
     !> The number of halo regions that this field has.
     !! Halo region values are not computed but copied
     !! from elsewhere.
     integer :: num_halos
     !> Array of objects describing the halos belonging to this field.
     type(halo_type), dimension(:), allocatable :: halo
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

  interface field_type
     module procedure field_constructor
  end interface field_type

contains

  !===================================================

  function field_constructor(grid_ptr, mesh_pts) result(self)
    implicit none
    type(grid_type), intent(in), pointer :: grid_ptr
    integer,         intent(in)          :: mesh_pts
    type(field_type), pointer :: self

    select case(mesh_pts)

    case(U_POINTS)
       self = cu_field_constructor(grid_ptr)
    case(V_POINTS)
       self%defined_on = mesh_pts
    case(T_POINTS)
       self%defined_on = mesh_pts
    case(F_POINTS)
       self%defined_on = mesh_pts
    
    case default
       stop "field_constructor: ERROR: invalid specifier for type of mesh points"
    end select

  end function field_constructor

  !===================================================

  function cu_field_constructor(grid_ptr) result(self)
    implicit none
    type(grid_type), intent(in), pointer :: grid_ptr
    type(field_type), target :: self
    ! Locals
    integer :: M, N

    self%defined_on = U_POINTS

    M = grid_ptr%nx - 1
    N = grid_ptr%ny - 1

    ! When updating a quantity on U points we write to:
    ! (using 'x' to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  x  x  x   j=1
    self%internal%xstart = 2
    self%internal%xstop  = M+1
    self%internal%ystart = 1
    self%internal%ystop  = N

    ! When applying periodic (wrap-around) boundary conditions (PBCs)
    ! we must fill the regions marked with an 'o' above.
    ! This looks like (using x to indicate a location that is written
    ! first and y a location that is written second):
    !
    !  i=1   i=M
    ! _ y  y  y  y 
    ! /|x  o  o  o   j=N 
    ! | x  o  o  o
    ! \ x  o  o  o
    !  \x_ o  o  o   j=1
    !   |\______/ 

    ! In array notation this looks like:
    !
    ! (1    , 1:N) = (M+1  , 1:N)
    ! (1:M+1, n+1) = (1:M+1, 1)

    self%num_halos = 2
    allocate(self%halo(self%num_halos))

    self%halo(1)%dest%xstart = 1   ; self%halo(1)%dest%xstop = 1
    self%halo(1)%dest%ystart = 1   ; self%halo(1)%dest%ystop = N
    self%halo(1)%source%xstart  = M+1 ; self%halo(1)%source%xstop  = M+1
    self%halo(1)%source%ystart  = 1   ; self%halo(1)%source%ystop  = N

    self%halo(2)%dest%xstart = 1   ; self%halo(2)%dest%xstop = M+1
    self%halo(2)%dest%ystart = N+1 ; self%halo(2)%dest%ystop = N+1
    self%halo(2)%source%xstart = 1 ; self%halo(2)%source%xstop  = M+1
    self%halo(2)%source%ystart = 1 ; self%halo(2)%source%ystop  = 1

  end function cu_field_constructor

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
        
  end subroutine copy_2dfield

  !===================================================

  subroutine copy_2dfield_patch(field, src, dest)
    implicit none
    real(wp),          intent(inout), dimension(:,:) :: field
    type(region_type), intent(in)                    :: src, dest

    field(dest%xstart:dest%xstop,dest%ystart:dest%ystop) = &
     field(src%xstart:src%xstop ,src%ystart:src%ystop)
        
  end subroutine copy_2dfield_patch

  !===================================================

  subroutine increment_scalar_field(field, incr)
    implicit none
    type(scalar_field_type), INTENT(inout) :: field
    type(scalar_field_type), INTENT(in)    :: incr

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
