module field_mod
  use kind_params_mod
  use region_mod
  use halo_mod
  use grid_mod
  implicit none

  private

  ! Enumeration of grid-point types on the Arakawa C grid. A
  ! field lives on one of these types.
  integer, public, parameter :: U_POINTS   = 0
  integer, public, parameter :: V_POINTS   = 1
  integer, public, parameter :: T_POINTS   = 2
  integer, public, parameter :: F_POINTS   = 3
  !> A field that lives on all grid-points of the grid
  integer, public, parameter :: ALL_POINTS = 4

  ! Enumeration of boundary-condition types
  !> Field has periodic boundary conditions (in both x and y)
  integer, public, parameter :: BC_PERIODIC = 0
  !> Field has external boundary conditions. This is a placeholder
  !! really as this is a complex area.
  integer, public, parameter :: BC_EXTERNAL = 1
  !> Field has no boundary conditions
  integer, public, parameter :: BC_NONE = 2

  type, public :: field_type
     !> Which mesh points the field is defined upon
     integer :: defined_on
     !> The grid on which this field is defined
     type(grid_type), pointer :: grid
     !> The type of boundary conditions applied to this field
     !! in the x, y and z dimensions. Note that at this stage
     !! this is really only required for Periodic Boundary
     !! conditions.
     integer, dimension(3) :: boundary_conditions
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

  !> Interface for the copy_field operation. Overloaded to take
  !! a scalar, an array or an r2d_field_type.
  !! \todo Remove support for raw arrays from this interface.
  interface copy_field
     module procedure copy_scalar_field,                            &
                      copy_2dfield_array, copy_2dfield_array_patch, &
                      copy_2dfield, copy_2dfield_patch
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

  public increment
  public copy_field
  public set
  public field_checksum

! Grid points on an Arakawa C grid with NE staggering are arranged like so:
!
! v(1,ny)-----f(1,ny)-- -v(i-1,ny)--f(i-1,ny)--v(i,ny)----f(i,ny)-  -v(nx,ny)---f(nx,ny)  
! |           |          |          |          |          |          |          |        
! |           |          |          |          |          |          |          |        
! T[1,ny]-----u(1,ny)-- -T(i-1,ny)--u(i-1,ny)--T(i,ny)----u(i,ny)-  -T(nx,ny)---u(nx,ny)  
! |           |          |          |          |          |          |          |        
! |           |          |          |          |          |          |          |        
! v(1,j)------f(1,j)--- -v(i-1,j)---f(i-1,j)---v(i,j)-----f(i,j)--  -v(nx,j)----f(nx,j)   
! |           |          |          |          |          |          |          |        
! |           |          |          |          |          |          |          |        
! T[1,j]------u(1,j)--- -T(i-1,j)---u(i-1,j)---T(i,j)-----u(i,j)--  -T(nx,j)----u(nx,j)   
! |           |          |          |          |          |          |          |        
! |           |          |          |          |          |          |          |        
! v(1,j-1)----f(1,j-1)- -v(i-1,j-1)-f(i-1,j-1)-v(i,j-1)---f(i,j-1)- -v(nx,j-1)--f(nx,j-1) 
! |           |          |          |          |          |          |          |        
! |           |          |          |          |          |          |          |        
! T[1,j-1]----u(1,j-1)- -T(i-1,j-1)-u(i-1,j-1)-T(i,j-1)---u(i,j-1)- -T(nx,j-1)--u(nx,j-1) 
! |           |          |          |          |          |          |          |        
! |           |          |          |          |          |          |          |        
! v(1,1)------f(1,1)--- -v(i-1,1)---f(i-1,1)---v(i,1)-----f(i,1)--  -v(nx,1)----f(nx,1)   
! |           |          |          |          |          |          |          |        
! |           |          |          |          |          |          |          |        
! T[1,1]      u(1,1)     T(i-1,1)---u(i-1,1)---T(i,1)-----u(i,1)--  -T(nx,1)----u(nx,1)   

contains

  !===================================================

  function field_constructor(grid_ptr,    &
                             grid_points, &
                             boundary_conditions) result(self)
    implicit none
    ! Arguments
    !> Pointer to the grid on which this field lives
    type(grid_type),       intent(in), pointer :: grid_ptr
    !> Which grid-point type the field is defined on
    integer,               intent(in)          :: grid_points
    !> The boundary conditions that this field is subject to
    integer, dimension(3), intent(in)          :: boundary_conditions
    ! Local declarations
    type(field_type) :: self

    stop 'field_constructor: ERROR: I should not have been called!'

  end function field_constructor

  !===================================================

  function r2d_field_constructor(grid,    &
                                 grid_points, &
                                 boundary_conditions) result(self)
    implicit none
    ! Arguments
    !> Pointer to the grid on which this field lives
    type(grid_type), intent(in), target  :: grid
    !> Which grid-point type the field is defined on
    integer,         intent(in)          :: grid_points
    !> The boundary conditions that this field is subject to in the
    !! x and y dimensions
    integer, dimension(2), intent(in)    :: boundary_conditions
    ! Local declarations
    type(r2d_field_type) :: self
    integer :: ierr

    self%boundary_conditions(1:2) = boundary_conditions(1:2)
    ! This is the constructor for a 2D field so we obviously
    ! don't have any BC's for the 3rd dimension.
    self%boundary_conditions(3) = BC_NONE

    ! Set this field's grid pointer to point to the grid pointed to
    ! by the supplied grid_ptr argument
    self%grid => grid

    allocate(self%data(self%grid%nx,self%grid%ny), Stat=ierr)
    if(ierr /= 0)then
       stop 'r2d_field_constructor: ERROR: failed to allocate field'
    end if

    select case(grid_points)

    case(U_POINTS)
       call cu_field_init(self)
    case(V_POINTS)
       call cv_field_init(self)
    case(T_POINTS)
       call ct_field_init(self)
    case(F_POINTS)
       call cf_field_init(self)
    case(ALL_POINTS)
       call field_init(self)
    case default
       stop 'r2d_field_constructor: ERROR: invalid specifier for type of mesh points'
    end select


    ! Compute and store dimensions of internal region of field
    self%internal%nx = self%internal%xstop - self%internal%xstart + 1
    self%internal%ny = self%internal%ystop - self%internal%ystart + 1

  end function r2d_field_constructor

  !===================================================

  subroutine field_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    fld%defined_on = ALL_POINTS

    M = fld%grid%nx
    N = fld%grid%ny

    ! An 'all points' field is defined upon every point in the grid
    fld%internal%xstart = 1
    fld%internal%xstop  = M
    fld%internal%ystart = 1
    fld%internal%ystop  = N

    ! We have no halo regions
    fld%num_halos = 0

  end subroutine field_init

  !===================================================

  subroutine cu_field_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld

    fld%defined_on = U_POINTS

    select case(fld%grid%stagger)

    case(STAGGER_SW)
       call cu_sw_init(fld)

    case(STAGGER_NE)
       call cu_ne_init(fld)

    case default
       stop 'cu_field_init: ERROR - unsupported stagger!'

    end select

  end subroutine cu_field_init

  !===================================================

  subroutine cu_sw_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    M = fld%grid%nx
    N = fld%grid%ny

    ! Set up a field defined on U points when the grid has
    ! a South-West staggering:

    !   vi-1j+1--fij+1---vij+1---fi+1j+1
    !   |        |       |       |
    !   |        |       |       |
    !   Ti-1j----uij-----Tij-----ui+1j
    !   |        |       |       |
    !   |        |       |       |
    !   vi-1j----fij-----vij-----fi+1j
    !   |        |       |       |
    !   |        |       |       |
    !   Ti-1j-1--uij-1---Tij-1---ui+1j-1

    ! When updating a quantity on U points with this staggering
    ! we write to (using 'x' to indicate a location that is written,
    !                    'b' a boundary point and
    !                    'o' a point that is external to the domain):
    !
    ! i=1         i=M
    !  o  b  b  b  b   j=N 
    !  o  b  x  x  b
    !  o  b  x  x  b
    !  o  b  x  x  b
    !  o  b  b  b  b   j=1
    !
    if(fld%boundary_conditions(1) == BC_PERIODIC)then
       ! When implementing periodic boundary conditions, all
       ! mesh point types have the same extents as the grid of
       ! T points. We then have a halo of width 1 on either side
       ! of the domain.
       fld%internal%xstart = 2
       fld%internal%xstop  = M-1
    else
       fld%internal%xstart = 3
       fld%internal%xstop  = M-1
    end if

    fld%internal%ystart = 2
    fld%internal%ystop  = N-1

    ! When applying periodic (wrap-around) boundary conditions (PBCs)
    ! we must fill the regions marked with 'b' above.
    ! This looks like (using x to indicate a location that is written
    ! first and y a location that is written second):
    !
    !  i=2      i=M
    ! _ y  y  y  y   j=N  
    ! /|x  o  o  o
    ! | x  o  o  o
    ! \ x  o  o  o
    !  \x_ o  o  o   j=1
    !   |\______/ 

    ! In array notation this looks like:
    !
    ! (2    , 1:N-1) = (M  , 1:N-1)
    ! (2:M, N) = (2:M, 1)

    fld%num_halos = 2
    allocate(fld%halo(fld%num_halos))

    fld%halo(1)%dest%xstart   = 1   ; fld%halo(1)%dest%xstop = 1
    fld%halo(1)%dest%ystart   = 1   ; fld%halo(1)%dest%ystop = N
    fld%halo(1)%source%xstart = M-1 ; fld%halo(1)%source%xstop = M-1
    fld%halo(1)%source%ystart = 1   ; fld%halo(1)%source%ystop = N

    fld%halo(2)%dest%xstart   = 1 ; fld%halo(2)%dest%xstop = M
    fld%halo(2)%dest%ystart   = N ; fld%halo(2)%dest%ystop = N
    fld%halo(2)%source%xstart = 1 ; fld%halo(2)%source%xstop = M
    fld%halo(2)%source%ystart = 2 ; fld%halo(2)%source%ystop = 2

  end subroutine cu_sw_init

  !===================================================

  subroutine cu_ne_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    M = fld%grid%nx
    N = fld%grid%ny

    ! Set up a field defined on U points when the grid has
    ! a North-East staggering:

    ! It is the T points that define the whole domain and we are
    ! simulating a region within this domain. As a minimum, we will
    ! require one external T point around the whole perimeter of the
    ! simulated domain in order to specify boundary conditions. The U
    ! pts on the boundary will then lie between the last external and
    ! first internal T points. 
    ! With a (N)E stagger this means:

    ! ji indexing: 
    ! Lowermost i index of the u points will be the same as the T's.
    ! i.e. if we start at 1 then T(1,:) are external and u(1,:) are 
    ! boundary points too.
    ! However, the U points with ji==nx will lie outside the model
    ! domain. U points with ji==nx-1 will be the Eastern-most *boundary*
    ! points.
    ! jj indexing:
    ! Lowermost j index of the U points - U pts with jj the same as external T points 
    ! will also be external to domain and therefore unused. U points with jj one greater 
    ! than lowest ext. T pts will be *boundary* points. U pts with jj==ny will be 
    ! boundary points.

    ! When updating a quantity on U points with this staggering
    ! we write to (using 'x' to indicate a location that is written and 'b' a boundary 
    ! point):
    !
    ! i= 1          nx-1  nx
    !    b   b   b   b    o   ny
    !    b   x   x   b    o 
    !    b   x   x   b    o 
    !    b   x   x   b    o   
    !    b   b   b   b    o   1
    !                         j

    ! i.e. fld(2:M,2:N+1) = ...

    fld%internal%xstart = 2
    fld%internal%xstop  = M-2
    fld%internal%ystart = 2
    fld%internal%ystop  = N-1

!> \todo Is this concept of halo definitions useful?
    fld%num_halos = 0
!    allocate(fld%halo(fld%num_halos))

!    fld%halo(1)%dest%xstart   = 1   ; fld%halo(1)%dest%xstop = 1
!    fld%halo(1)%dest%ystart   = 1   ; fld%halo(1)%dest%ystop = N
!    fld%halo(1)%source%xstart = M+1 ; fld%halo(1)%source%xstop = M+1
!    fld%halo(1)%source%ystart = 1   ; fld%halo(1)%source%ystop = N

!    fld%halo(2)%dest%xstart   = 1   ; fld%halo(2)%dest%xstop = M+1
!    fld%halo(2)%dest%ystart   = N+1 ; fld%halo(2)%dest%ystop = N+1
!    fld%halo(2)%source%xstart = 1 ; fld%halo(2)%source%xstop = M+1
!    fld%halo(2)%source%ystart = 1 ; fld%halo(2)%source%ystop = 1

  end subroutine cu_ne_init

  !===================================================

  subroutine cv_field_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld

    fld%defined_on = V_POINTS

    select case(fld%grid%stagger)

    case(STAGGER_SW)
       call cv_sw_init(fld)

    case(STAGGER_NE)
       call cv_ne_init(fld)

    case default
       stop 'cv_field_init: ERROR - unsupported stagger!'

    end select

  end subroutine cv_field_init

  !===================================================

  subroutine cv_sw_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    M = fld%grid%nx
    N = fld%grid%ny

    ! When updating a quantity on V points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1      i=M
    !  b  b  b  b   j=N 
    !  b  x  x  b
    !  b  x  x  b
    !  b  b  b  b
    !  o  o  o  o   j=1
    fld%internal%xstart = 2
    fld%internal%xstop  = M-1
    if(fld%boundary_conditions(2) == BC_PERIODIC)then
       ! When implementing periodic boundary conditions, all
       ! mesh point types have the same extents as the grid of
       ! T points. We then have a halo of width 1 on either side
       ! of the domain.
       fld%internal%ystart = 2
       fld%internal%ystop  = N-1
    else
       fld%internal%ystart = 3
       fld%internal%ystop  = N-1
    endif

    ! When applying periodic (wrap-around) boundary conditions (PBCs)
    ! we must fill the regions marked with 'b' above.
    ! This looks like (using x to indicate a location that is written
    ! first and y a location that is written second):
    !
    !  i=1     i=M
    !  -o  o  o  y   j=N  
    ! / o  o  o  y
    ! | o  o  o  y
    ! \ o  o  o  y
    !  \x  x  x _y   j=2
    !    \______/|

    ! In array notation this looks like:
    ! First row = last internal row
    ! field(1:M    ,1:1  ) = field(1:M,N-1:N-1)
    ! Last col = first internal col
    ! field(M:M,1:N) = field(2:2,  1:N)

    fld%num_halos = 2
    ALLOCATE(fld%halo(fld%num_halos))

    fld%halo(1)%dest%xstart = 1   ; fld%halo(1)%dest%xstop = M
    fld%halo(1)%dest%ystart = 1   ; fld%halo(1)%dest%ystop = 1
    fld%halo(1)%source%xstart  = 1   ; fld%halo(1)%source%xstop  = M
    fld%halo(1)%source%ystart  = N-1 ; fld%halo(1)%source%ystop  = N-1

    fld%halo(2)%dest%xstart = M ; fld%halo(2)%dest%xstop = M
    fld%halo(2)%dest%ystart = 1   ; fld%halo(2)%dest%ystop = N
    fld%halo(2)%source%xstart  = 2   ; fld%halo(2)%source%xstop  = 2
    fld%halo(2)%source%ystart  = 1   ; fld%halo(2)%source%ystop  = N

  end subroutine cv_sw_init

  !===================================================

  subroutine cv_ne_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    M = fld%grid%nx
    N = fld%grid%ny

    ! ji indexing:
    ! Lowermost ji index of the V points will be the same as the T's.
    ! If the domain starts at 1 then T(1,:) are external and v(1,:)
    ! are boundary points.
    ! Uppermost ji index is nx. T(nx,:) are external and v(nx,:)
    ! are boundary points.

    ! jj indexing:
    ! If domain starts at 1 then T(:,1) are external and V(:,1) are
    ! boundary points.
    ! Uppermost jj index is ny. T(ny,:) are external and so are V(ny,:)
    ! (see diagram at start of module). It is V(ny-1,:) that are the 
    ! boundary points.

    ! When updating a quantity on V points with this staggering
    ! we write to (using 'x' to indicate a location that is written):
    !
    ! i=1       Nx
    !  o  o  o  o   Ny
    !  b  b  b  b   Ny-1
    !  b  x  x  b
    !  b  x  x  b
    !  b  b  b  b   j=1
    !

    fld%internal%xstart = 2
    fld%internal%xstop  = M-1
    fld%internal%ystart = 2
    fld%internal%ystop  = N-2

  end subroutine cv_ne_init

  !===================================================

  subroutine ct_field_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld

    fld%defined_on = T_POINTS

    select case(fld%grid%stagger)

    case(STAGGER_SW)
       call ct_sw_init(fld)

    case(STAGGER_NE)
       call ct_ne_init(fld)

    case default
       stop 'ct_field_init: ERROR - unsupported stagger!'

    end select

  end subroutine ct_field_init

  !===================================================

  subroutine ct_sw_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    M = fld%grid%nx
    N = fld%grid%ny

    ! When updating a quantity on T points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1      i=M
    !  b  b  b  b   j=N 
    !  b  x  x  b
    !  b  x  x  b
    !  b  b  b  b   j=1

    fld%internal%xstart = 2
    fld%internal%xstop  = M-1
    fld%internal%ystart = 2
    fld%internal%ystop  = N-1

    ! When applying periodic (wrap-around) boundary conditions
    ! (PBCs) we must fill the regions marked with 'b' above.
    ! This looks like (using x to indicate a location that is 
    ! written first and y a location that is written second):
    !
    !  i=1      i=M
    ! _ y  y  y  y   j=N  
    ! /|o  o  o  x
    ! | o  o  o  x
    ! \ o  o  o  x
    !  \o  o  o _x   j=1
    !    \______/|

    ! In array notation this looks like:
    ! Last col = first col
    ! field(M:M,  1:N-1  ) = field(1:1  ,1:N-1)
    ! Last row = first row
    ! field(1:M,N:N) = field(1:M,1:1)

    fld%num_halos = 2
    ALLOCATE( fld%halo(fld%num_halos) )

    fld%halo(1)%dest%xstart = M ; fld%halo(1)%dest%xstop = M
    fld%halo(1)%dest%ystart = 1   ; fld%halo(1)%dest%ystop = N-1
    fld%halo(1)%source%xstart  = 2   ; fld%halo(1)%source%xstop  = 2
    fld%halo(1)%source%ystart  = 1   ; fld%halo(1)%source%ystop  = N-1

    fld%halo(2)%dest%xstart = 1   ; fld%halo(2)%dest%xstop = M
    fld%halo(2)%dest%ystart = N ; fld%halo(2)%dest%ystop = N
    fld%halo(2)%source%xstart  = 1   ; fld%halo(2)%source%xstop  = M
    fld%halo(2)%source%ystart  = 2   ; fld%halo(2)%source%ystop  = 2

  end subroutine ct_sw_init

  !===================================================

  subroutine ct_ne_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    M = fld%grid%nx
    N = fld%grid%ny

    ! The mesh of T points defines the simulation domain. 
    ! Currently we assume a shell of thickness one around the actual
    ! simulation domain - this is the minimum required to specify
    ! boundary conditions.

    ! When updating a quantity on T points with a NE staggering
    ! we write to (using x to indicate a location that is written):
    !
    ! i=1          Nx
    !  b  b  b  b  b Ny
    !  b  x  x  x  b
    !  b  x  x  x  b 
    !  b  x  x  x  b
    !  b  b  b  b  b j=1
    fld%internal%nx = M-2
    fld%internal%ny = N-2

    fld%internal%xstart = 2
    fld%internal%xstop  = M-1
    fld%internal%ystart = 2
    fld%internal%ystop  = N-1

  end subroutine ct_ne_init

  !===================================================

  subroutine cf_field_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld

    fld%defined_on = F_POINTS

    select case(fld%grid%stagger)

    case(STAGGER_SW)
       call cf_sw_init(fld)

    case(STAGGER_NE)
       call cf_ne_init(fld)

    case default
       stop 'cf_field_init: ERROR - unsupported stagger!'

    end select

  end subroutine cf_field_init

  !===================================================

  subroutine cf_sw_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    M = fld%grid%nx
    N = fld%grid%ny

    ! When updating a quantity on F points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1         i=M
    !  o  b  b  b  b   j=N 
    !  o  b  x  x  b
    !  o  b  x  x  b
    !  o  b  b  b  b
    !  o  o  o  o  o   j=1
    if(fld%boundary_conditions(1) == BC_PERIODIC)then
       fld%internal%xstart = 2
       fld%internal%xstop  = M-1
    else
       fld%internal%xstart = 3
       fld%internal%xstop  = M-1
    end if

    if(fld%boundary_conditions(2) == BC_PERIODIC)then
       fld%internal%ystart = 2
       fld%internal%ystop  = N-1
    else
       fld%internal%ystart = 3
       fld%internal%ystop  = N-1
    end if

    ! When applying periodic (wrap-around) boundary conditions
    ! (PBCs) we must fill the regions marked with 'b' above.
    ! This looks like (using x to indicate a location that is 
    ! written first and y a location that is written second):
    !
    !  i=2      i=M
    ! .-x  o  o  o   j=N  
    ! | x  o  o  o
    ! | x  o  o  o
    ! | x  o  o  o
    ! ->y_ y  y  y   j=2
    !   |\______/ 

    ! In array notation this looks like:
    ! First col = last col
    ! field(2:2, 2:N) = field(M:M, 2:N)
    ! First row = last row
    ! field(2:M, 2:2) = field(2:M, N:N)

    fld%num_halos = 2
    ALLOCATE( fld%halo(fld%num_halos) )

    fld%halo(1)%dest%xstart = 1   ; fld%halo(1)%dest%xstop = 1
    fld%halo(1)%dest%ystart = 1   ; fld%halo(1)%dest%ystop = N
    fld%halo(1)%source%xstart  = M-1 ; fld%halo(1)%source%xstop  = M-1
    fld%halo(1)%source%ystart  = 1   ; fld%halo(1)%source%ystop  = N

    fld%halo(2)%dest%xstart = 1   ; fld%halo(2)%dest%xstop = M
    fld%halo(2)%dest%ystart = 1   ; fld%halo(2)%dest%ystop = 1
    fld%halo(2)%source%xstart  = 1  ; fld%halo(2)%source%xstop  = M
    fld%halo(2)%source%ystart  = N-1 ; fld%halo(2)%source%ystop  = N-1

  end subroutine cf_sw_init

  !===================================================

  subroutine cf_ne_init(fld)
    implicit none
    type(r2d_field_type), intent(inout) :: fld
    ! Locals
    integer :: M, N

    M = fld%grid%nx
    N = fld%grid%ny

    ! When updating a quantity on F points we write to:
    ! (using x to indicate a location that is written
    !        b a boundary point - defined by ext. b.c.
    !        o a point that is external to the domain):
    !
    ! i=1       Nx
    !  o  o  o  o   Ny
    !  b  b  b  o   
    !  b  x  b  o   
    !  b  x  b  o
    !  b  b  b  o   j=1

    fld%internal%xstart = 2
    fld%internal%xstop  = M-2
    fld%internal%ystart = 2
    fld%internal%ystop  = N-2

  end subroutine cf_ne_init

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

  SUBROUTINE copy_2dfield(field_in, field_out)
    IMPLICIT none
    type(r2d_field_type), intent(in)    :: field_in
    type(r2d_field_type), intent(inout) :: field_out
        
    field_out%data(:,:) = field_in%data(:,:)
        
  end subroutine copy_2dfield

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

  function field_checksum(field) result(val)
    implicit none
    type(r2d_field_type), intent(in) :: field
    real(wp) :: val

    val = SUM(field%data(field%internal%xstart:field%internal%xstop, &
                         field%internal%ystart:field%internal%ystop))

  end function field_checksum

  !===================================================

end module field_mod
