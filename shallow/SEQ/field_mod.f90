module field_mod
  use kind_params_mod
  use region_mod
  use halo_mod
  use grid_mod
  implicit none

  private

  ! Enumeration of grid-point types on the Arakawa C grid. A
  ! field lives on one of these types.
  integer, public, parameter :: U_POINTS = 0
  integer, public, parameter :: V_POINTS = 1
  integer, public, parameter :: T_POINTS = 2
  integer, public, parameter :: F_POINTS = 3

  ! Enumeration of boundary-condition types
  !> Field has periodic boundary conditions (in both x and y)
  integer, public, parameter :: BC_PERIODIC = 0
  !> Field has external boundary conditions. This is a placeholder
  !! really as this is a complex area.
  integer, public, parameter :: BC_EXTERNAL = 1

  type, public :: field_type
     !> Which mesh points the field is defined upon
     integer :: defined_on
     !> The grid on which this field is defined
     type(grid_type), pointer :: grid
     !> The type of boundary conditions applied to this field
     integer :: boundary_conditions
     !> The internal region of this field
     type(region_type) :: internal
     !> The number of halo regions that this field has.
     !! Halo region values are not computed but copied
     !! from elsewhere.
     integer :: num_halos
     !> Array of objects describing the halos belonging to this field.
     type(halo_type), dimension(:), allocatable :: halo
  end type field_type

  type, public, extends(field_type) :: scalar_field_type
     real(wp) :: data
  end TYPE scalar_field_type

  type, public, extends(field_type) :: r2d_field_type
     !> Array holding the actual field values
     real(wp), dimension(:,:), allocatable :: data
  end type r2d_field_type

  interface set
     module procedure set_scalar_field
  end interface set

  interface copy_field
     module procedure copy_scalar_field, copy_2dfield_array, &
                      copy_2dfield_array_patch, copy_2dfield_patch
  end interface copy_field

  interface increment
     module procedure increment_scalar_field
  end interface increment

  interface field_type
     module procedure field_constructor
  end interface field_type

  ! User-defined constructor for r2d_field_type objects
  interface r2d_field_type
     module procedure r2d_field_constructor
  end interface r2d_field_type

contains

  !===================================================

  function field_constructor(grid_ptr,    &
                             grid_points, &
                             boundary_conditions) result(self)
    implicit none
    ! Arguments
    !> Pointer to the grid on which this field lives
    type(grid_type), intent(in), pointer :: grid_ptr
    !> Which grid-point type the field is defined on
    integer,         intent(in)          :: grid_points
    !> The boundary conditions that this field is subject to
    integer,         intent(in)          :: boundary_conditions
    ! Local declarations
    type(field_type) :: self

  end function field_constructor

  !===================================================

  function r2d_field_constructor(grid_ptr,    &
                                 grid_points, &
                                 boundary_conditions) result(self)
    implicit none
    ! Arguments
    !> Pointer to the grid on which this field lives
    type(grid_type), intent(in), pointer :: grid_ptr
    !> Which grid-point type the field is defined on
    integer,         intent(in)          :: grid_points
    !> The boundary conditions that this field is subject to
    integer,         intent(in)          :: boundary_conditions
    ! Local declarations
    type(field_type) :: self

    self%boundary_conditions = boundary_conditions
    ! Set this field's grid pointer to point to the grid pointed to
    ! by the supplied grid_ptr argument
    self%grid => grid_ptr

    select case(grid_points)

    case(U_POINTS)
       call cu_field_init(self)
    case(V_POINTS)
       call cv_field_init(self)
    case(T_POINTS)
       call ct_field_init(self)
    case(F_POINTS)
       call cf_field_init(self)
    case default
       stop "r2d_field_constructor: ERROR: invalid specifier for type of mesh points"
    end select


  end function r2d_field_constructor

  !===================================================

  subroutine cu_field_init(fld)
    implicit none
    type(field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    fld%defined_on = U_POINTS

    M = fld%grid%nx - 1
    N = fld%grid%ny - 1

    ! When updating a quantity on U points we write to:
    ! (using 'x' to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  x  x  x   j=1
    fld%internal%xstart = 2
    fld%internal%xstop  = M+1
    fld%internal%ystart = 1
    fld%internal%ystop  = N

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

    fld%num_halos = 2
    allocate(fld%halo(fld%num_halos))

    fld%halo(1)%dest%xstart   = 1   ; fld%halo(1)%dest%xstop = 1
    fld%halo(1)%dest%ystart   = 1   ; fld%halo(1)%dest%ystop = N
    fld%halo(1)%source%xstart = M+1 ; fld%halo(1)%source%xstop = M+1
    fld%halo(1)%source%ystart = 1   ; fld%halo(1)%source%ystop = N

    fld%halo(2)%dest%xstart   = 1   ; fld%halo(2)%dest%xstop = M+1
    fld%halo(2)%dest%ystart   = N+1 ; fld%halo(2)%dest%ystop = N+1
    fld%halo(2)%source%xstart = 1 ; fld%halo(2)%source%xstop = M+1
    fld%halo(2)%source%ystart = 1 ; fld%halo(2)%source%ystop = 1

  end subroutine cu_field_init

  !===================================================

  subroutine cv_field_init(fld)
    implicit none
    type(field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    fld%defined_on = V_POINTS

    M = fld%grid%nx - 1
    N = fld%grid%ny - 1

    ! When updating a quantity on V points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  x  x  x  o 
    !  x  x  x  o   j=N
    !  x  x  x  o
    !  o  o  o  o   j=1
    fld%internal%xstart = 1
    fld%internal%xstop  = M
    fld%internal%ystart = 2
    fld%internal%ystop  = N+1

    ! When applying periodic (wrap-around) boundary conditions (PBCs)
    ! we must fill the regions marked with an 'o' above.
    ! This looks like (using x to indicate a location that is written
    ! first and y a location that is written second):
    !
    !  i=1   i=M
    !  -o  o  o  y 
    ! / o  o  o  y   j=N 
    ! | o  o  o  y
    ! \ o  o  o  y
    !  \x  x  x _y   j=1
    !    \______/|

    ! In array notation this looks like:
    ! First row = last row
    ! field(1:M    ,1:1  ) = field(1:M,NP1:NP1)
    ! Last col = first col
    ! field(MP1:MP1,1:NP1) = field(1:1,  1:NP1)

    fld%num_halos = 2
    ALLOCATE(fld%halo(fld%num_halos))

    fld%halo(1)%dest%xstart = 1   ; fld%halo(1)%dest%xstop = M
    fld%halo(1)%dest%ystart = 1   ; fld%halo(1)%dest%ystop = 1
    fld%halo(1)%source%xstart  = 1   ; fld%halo(1)%source%xstop  = M
    fld%halo(1)%source%ystart  = N+1 ; fld%halo(1)%source%ystop  = N+1

    fld%halo(2)%dest%xstart = M+1 ; fld%halo(2)%dest%xstop = M+1
    fld%halo(2)%dest%ystart = 1   ; fld%halo(2)%dest%ystop = N+1
    fld%halo(2)%source%xstart  = 1   ; fld%halo(2)%source%xstop  = 1
    fld%halo(2)%source%ystart  = 1   ; fld%halo(2)%source%ystop  = N+1

  end subroutine cv_field_init

  !===================================================

  subroutine ct_field_init(fld)
    implicit none
    type(field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    fld%defined_on = T_POINTS

    M = fld%grid%nx - 1
    N = fld%grid%ny - 1

    ! When updating a quantity on T points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  x  x  x  o   j=N
    !  x  x  x  o
    !  x  x  x  o   j=1
    fld%internal%xstart = 1
    fld%internal%xstop  = M
    fld%internal%ystart = 1
    fld%internal%ystop  = N

    ! When applying periodic (wrap-around) boundary conditions
    ! (PBCs) we must fill the regions marked with an 'o' above.
    ! This looks like (using x to indicate a location that is 
    ! written first and y a location that is written second):
    !
    !  i=1   i=M
    ! _ y  y  y  y 
    ! /|o  o  o  x   j=N 
    ! | o  o  o  x
    ! \ o  o  o  x
    !  \o  o  o _x   j=1
    !    \______/|

    ! In array notation this looks like:
    ! Last col = first col
    ! field(MP1:MP1,  1:N  ) = field(1:1  ,1:N)
    ! Last row = first row
    ! field(1  :MP1,NP1:NP1) = field(1:MP1,1:1)

    fld%num_halos = 2
    ALLOCATE( fld%halo(fld%num_halos) )

    fld%halo(1)%dest%xstart = M+1 ; fld%halo(1)%dest%xstop = M+1
    fld%halo(1)%dest%ystart = 1   ; fld%halo(1)%dest%ystop = N
    fld%halo(1)%source%xstart  = 1   ; fld%halo(1)%source%xstop  = 1
    fld%halo(1)%source%ystart  = 1   ; fld%halo(1)%source%ystop  = N

    fld%halo(2)%dest%xstart = 1   ; fld%halo(2)%dest%xstop = M+1
    fld%halo(2)%dest%ystart = N+1 ; fld%halo(2)%dest%ystop = N+1
    fld%halo(2)%source%xstart  = 1   ; fld%halo(2)%source%xstop  = M+1
    fld%halo(2)%source%ystart  = 1   ; fld%halo(2)%source%ystop  = 1

  end subroutine ct_field_init

  !===================================================

  subroutine cf_field_init(fld)
    implicit none
    type(field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    fld%defined_on = F_POINTS

    M = fld%grid%nx - 1
    N = fld%grid%ny - 1

    ! When updating a quantity on F points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  x  x  x 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  o  o  o   j=1
    fld%internal%xstart = 2
    fld%internal%xstop  = M+1
    fld%internal%ystart = 2
    fld%internal%ystop  = N+1

    ! When applying periodic (wrap-around) boundary conditions
    ! (PBCs) we must fill the regions marked with an 'o' above.
    ! This looks like (using x to indicate a location that is 
    ! written first and y a location that is written second):
    !
    !  i=1   i=M
    ! .-x  o  o  o 
    ! | x  o  o  o   j=N 
    ! | x  o  o  o
    ! | x  o  o  o
    ! ->y_ y  y  y   j=1
    !   |\______/ 

    ! In array notation this looks like:
    ! First col = last col
    ! field(1:1  ,2:NP1) = field(MP1:MP1,  2:NP1)
    ! First row = last row
    ! field(1:MP1,1:1  ) = field(  1:MP1,NP1:NP1)

    fld%num_halos = 2
    ALLOCATE( fld%halo(fld%num_halos) )

    fld%halo(1)%dest%xstart = 1   ; fld%halo(1)%dest%xstop = 1
    fld%halo(1)%dest%ystart = 2   ; fld%halo(1)%dest%ystop = N+1
    fld%halo(1)%source%xstart  = M+1 ; fld%halo(1)%source%xstop  = M+1
    fld%halo(1)%source%ystart  = 2   ; fld%halo(1)%source%ystop  = N+1

    fld%halo(2)%dest%xstart = 1   ; fld%halo(2)%dest%xstop = M+1
    fld%halo(2)%dest%ystart = 1   ; fld%halo(2)%dest%ystop = 1
    fld%halo(2)%source%xstart  = 1   ; fld%halo(2)%source%xstop  = M+1
    fld%halo(2)%source%ystart  = N+1 ; fld%halo(2)%source%ystop  = N+1

  end subroutine cf_field_init

  !===================================================

  SUBROUTINE copy_scalar_field(field_in, field_out)
    IMPLICIT none
    type(scalar_field_type), INTENT(in) :: field_in
    type(scalar_field_type), INTENT(out) :: field_out

    field_out = field_in

  END SUBROUTINE copy_scalar_field

  !===================================================

  SUBROUTINE copy_2dfield_array(field_in, field_out)
    IMPLICIT none
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: field_in
    REAL(wp), INTENT(out), DIMENSION(:,:) :: field_out
        
    field_out(:,:) = field_in(:,:)
        
  end subroutine copy_2dfield_array

  !===================================================

  !> Copy from one patch in an array to another patch
  !! Will be superceded by interface using field object
  !! instead of array data.
  subroutine copy_2dfield_array_patch(field, src, dest)
    implicit none
    real(wp),          intent(inout), dimension(:,:) :: field
    type(region_type), intent(in)                    :: src, dest

    field(dest%xstart:dest%xstop,dest%ystart:dest%ystop) = &
     field(src%xstart:src%xstop ,src%ystart:src%ystop)
        
  end subroutine copy_2dfield_array_patch

  !===================================================

  !> Copy from one patch to another in a field
  subroutine copy_2dfield_patch(field, src, dest)
    implicit none
    type(r2d_field_type), intent(inout) :: field
    type(region_type),    intent(in)    :: src, dest

    field%data(dest%xstart:dest%xstop,dest%ystart:dest%ystop) = &
     field%data(src%xstart:src%xstop ,src%ystart:src%ystop)
        
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
