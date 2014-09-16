module field_mod
  use kind_params_mod
  use region_mod
  use halo_mod
  use grid_mod
  use gocean_mod, only: gocean_stop
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

  type, public :: field_type
     !> Which mesh points the field is defined upon
     integer :: defined_on
     !> The grid on which this field is defined
     type(grid_type), pointer :: grid
     !> The internal region of this field
     type(region_type) :: internal
     !> The whole region covered by this field - includes 
     !! boundary points
     type(region_type) :: whole
     !> The number of halo regions that this field has.
     !! Halo region values are not computed but copied
     !! from elsewhere.
     integer :: num_halos
     !> Array of objects describing the halos belonging to this field.
     type(halo_type), dimension(:), allocatable :: halo
  end type field_type

  type, public, extends(field_type) :: scalar_field
     real(wp) :: data
  end TYPE scalar_field

  type, public, extends(field_type) :: r2d_field
     !> Array holding the actual field values
     real(wp), dimension(:,:), allocatable :: data
  end type r2d_field

  !interface set_field
  !   module procedure set_scalar_field
  !end interface set_field

  !> Interface for the copy_field operation. Overloaded to take
  !! a scalar, an array or an r2d_field type.
  !! \todo Remove support for raw arrays from this interface.
  interface copy_field
     module procedure copy_scalar_field,                            &
                      copy_2dfield_array, copy_2dfield_array_patch, &
                      copy_2dfield, copy_2dfield_patch
  end interface copy_field

  interface increment_field
     module procedure increment_scalar_field, increment_scalar_field_r8
  end interface increment_field

!  interface field_type
!     module procedure field_constructor
!  end interface field_type

  ! User-defined constructor for r2d_field type objects
  interface r2d_field
     module procedure r2d_field_constructor
  end interface r2d_field

  public increment_field
  public copy_field
  public set_field
  public field_checksum

! Grid points on an Arakawa C grid with NE staggering are arranged like so:
!
!v(1,ny)----f(1,ny)---v(i-1,ny)--f(i-1,ny)--v(i,ny)----f(i,ny)--v(nx,ny)---f(nx,ny)  
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!T[1,ny]----u(1,ny)---T(i-1,ny)--u(i-1,ny)--T(i,ny)----u(i,ny)--T(nx,ny)---u(nx,ny)  
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!v(1,j)-----f(1,j)----v(i-1,j)---f(i-1,j)---v(i,j)-----f(i,j)---v(nx,j)----f(nx,j)   
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!T[1,j]-----u(1,j)----T(i-1,j)---u(i-1,j)---T(i,j)-----u(i,j)---T(nx,j)----u(nx,j)   
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!v(1,j-1)---f(1,j-1)--v(i-1,j-1)-f(i-1,j-1)-v(i,j-1)---f(i,j-1)-v(nx,j-1)--f(nx,j-1) 
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!T[1,j-1]---u(1,j-1)--T(i-1,j-1)-u(i-1,j-1)-T(i,j-1)---u(i,j-1)-T(nx,j-1)--u(nx,j-1) 
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!v(1,1)-----f(1,1)----v(i-1,1)---f(i-1,1)---v(i,1)-----f(i,1)---v(nx,1)----f(nx,1)   
!|          |         |          |          |          |        |          |        
!|          |         |          |          |          |        |          |        
!T[1,1]     u(1,1)    T(i-1,1)---u(i-1,1)---T(i,1)-----u(i,1)---T(nx,1)----u(nx,1)   

  !> The no. of cols/rows used to define boundary data in the absence
  !! of periodic BCs.
  integer, parameter :: NBOUNDARY = 1

contains

  !===================================================

!  function field_constructor(grid_ptr,    &
!                             grid_points) result(self)
!    implicit none
!    ! Arguments
!    !> Pointer to the grid on which this field lives
!    type(grid_type),       intent(in), pointer :: grid_ptr
!    !> Which grid-point type the field is defined on
!    integer,               intent(in)          :: grid_points
!    ! Local declarations
!    type(field_type) :: self
!
!    stop 'field_constructor: ERROR: I should not have been called!'
!
!  end function field_constructor

  !===================================================

  function r2d_field_constructor(grid,    &
                                 grid_points) result(self)
    implicit none
    ! Arguments
    !> Pointer to the grid on which this field lives
    type(grid_type), intent(in), target  :: grid
    !> Which grid-point type the field is defined on
    integer,         intent(in)          :: grid_points
    ! Local declarations
    type(r2d_field) :: self
    integer :: ierr

    ! Set this field's grid pointer to point to the grid pointed to
    ! by the supplied grid_ptr argument
    self%grid => grid

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
       stop 'r2d_field_constructor: ERROR: invalid specifier for '//&
            'type of mesh points'
    end select

    ! Compute and store dimensions of internal region of field
    self%internal%nx = self%internal%xstop - self%internal%xstart + 1
    self%internal%ny = self%internal%ystop - self%internal%ystart + 1

    ! In addition to the 'internal region' of the field, we may have
    ! external points that define B.C.'s or that act as halos. Here
    ! we store the full extent of the field, inclusive of such
    ! points.
    !> \todo Replace the use of NBOUNDARY here with info. computed
    !! from the T-point mask.
    if(self%grid%boundary_conditions(1) /= BC_PERIODIC)then
       self%whole%xstart = self%internal%xstart - NBOUNDARY
       self%whole%xstop  = self%internal%xstop  + NBOUNDARY
    else
       self%whole%xstart = self%internal%xstart - NBOUNDARY
       self%whole%xstop  = self%internal%xstop  + NBOUNDARY
    end if
    if(self%grid%boundary_conditions(2) /= BC_PERIODIC)then
       self%whole%ystart = self%internal%ystart - NBOUNDARY
       self%whole%ystop  = self%internal%ystop  + NBOUNDARY
    else
       self%whole%ystart = self%internal%ystart - NBOUNDARY
       self%whole%ystop  = self%internal%ystop  + NBOUNDARY
    end if

    write(*,*) 'allocating field with bounds: (', &
               1,':', self%internal%xstop+1, ',', &
               1,':',self%internal%ystop+1,')'
    !allocate(self%data(self%internal%xstart-1:self%internal%xstop+1, &
    !                   self%internal%ystart-1:self%internal%ystop+1),&
    !                   Stat=ierr)
    ! Allocating with a lower bound != 1 causes problems whenever
    ! array passed as assumed-shape dummy argument because lower
    ! bounds default to 1 in called unit.
    ! However, all loops will be in the generated, middle layer and
    ! the generator knows the array bounds. This may give us the
    ! ability to solve this problem.
    allocate(self%data(1:self%internal%xstop+1, &
                       1:self%internal%ystop+1),&
                       Stat=ierr)
    if(ierr /= 0)then
       call gocean_stop('r2d_field_constructor: ERROR: failed to allocate field')
    end if

  end function r2d_field_constructor

  !===================================================

  subroutine field_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld
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
    type(r2d_field), intent(inout) :: fld

    fld%defined_on = U_POINTS

    select case(fld%grid%stagger)

    case(STAGGER_SW)
       call cu_sw_init(fld)

    case(STAGGER_NE)
       call cu_ne_init(fld)

    case default
       call gocean_stop('cu_field_init: ERROR - unsupported stagger!')

    end select

  end subroutine cu_field_init

  !===================================================

  subroutine cu_sw_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld

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

    !
    if(fld%grid%boundary_conditions(1) == BC_PERIODIC)then
       ! When implementing periodic boundary conditions, all mesh
       ! point types have the same extents as the grid of T points. We
       ! then have a halo of width grid_mod::HALO_WIDTH_X on either
       ! side of the domain.
       ! When updating a quantity on U points we write to:
       ! (using 'x' to indicate a location that is written and 'a' an
       ! additional point that shallow dispenses with):
       !
       ! i=istart i=istop
       !  o  o  o  o  a
       !  o  x  x  x  a  j=ystop
       !  o  x  x  x  a
       !  o  x  x  x  a  j=ystart
       !  a  a  a  a  a

       fld%internal%xstart = fld%grid%simulation_domain%xstart
       fld%internal%xstop  = fld%grid%simulation_domain%xstop
    else
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
       fld%internal%xstart = fld%grid%simulation_domain%xstart + 1
       fld%internal%xstop  = fld%grid%simulation_domain%xstop
    end if

    fld%internal%ystart = fld%grid%simulation_domain%ystart
    fld%internal%ystop  = fld%grid%simulation_domain%ystop

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

    call init_periodic_bc_halos(fld)

  end subroutine cu_sw_init

  !===================================================

  subroutine cu_ne_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld

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
    ! Lowermost j index of the U points - U pts with jj the same as
    ! external T points will also be external to domain and therefore
    ! unused. U points with jj one greater than lowest ext. T pts will
    ! be *boundary* points. U pts with jj==ny will be boundary points.

    ! When updating a quantity on U points with this staggering
    ! we write to (using 'x' to indicate a location that is written and 
    ! 'b' a boundary point):
    !
    ! i= 1          nx-1  nx
    !    b   b   b   b    o   ny
    !    b   x   x   b    o 
    !    b   x   x   b    o 
    !    b   x   x   b    o   
    !    b   b   b   b    o   1
    !                         j

    ! i.e. fld(2:M,2:N+1) = ...

    if(fld%grid%boundary_conditions(1) /= BC_PERIODIC)then
      ! If we do not have periodic boundary conditions then we do
      ! not need to allow for boundary points here - they are
      ! already contained within the region defined by T mask.
      ! The T mask has been used to determine the grid%simulation_domain
      ! which describes the area on the grid that is actually being
      ! modelled (as opposed to having values supplied from B.C.'s etc.)
      fld%internal%xstart = fld%grid%simulation_domain%xstart
      fld%internal%xstop  = fld%grid%simulation_domain%xstop - 1
    else
      write(*,*) 'ERROR: cu_ne_init: implement periodic boundary conditions!'
      stop
    end if
    if(fld%grid%boundary_conditions(2) /= BC_PERIODIC)then
      fld%internal%ystart = fld%grid%simulation_domain%ystart
      fld%internal%ystop  = fld%grid%simulation_domain%ystop
    else
      write(*,*) 'ERROR: cu_ne_init: implement periodic boundary conditions!'
      stop
    end if

!> \todo Is this concept of halo definitions useful?
    fld%num_halos = 0

  end subroutine cu_ne_init

  !===================================================

  subroutine cv_field_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld

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
    type(r2d_field), intent(inout) :: fld

    if(fld%grid%boundary_conditions(2) == BC_PERIODIC)then
       ! When implementing periodic boundary conditions, all
       ! mesh point types have the same extents as the grid of
       ! T points. We then have a halo of width 1 on either side
       ! of the domain.
       fld%internal%xstart = fld%grid%simulation_domain%xstart
       fld%internal%xstop  = fld%grid%simulation_domain%xstop

       fld%internal%ystart = fld%grid%simulation_domain%ystart
       fld%internal%ystop  = fld%grid%simulation_domain%ystop
    else
       ! When updating a quantity on V points we write to:
       ! (using x to indicate a location that is written):
       !
       ! i=1      i=M
       !  b  b  b  b   j=N 
       !  b  x  x  b
       !  b  x  x  b
       !  b  b  b  b
       !  o  o  o  o   j=1
       ! We're not offset from the T points in the x dimension so
       ! we have the same x bounds.
       fld%internal%xstart = fld%grid%simulation_domain%xstart
       fld%internal%xstop  = fld%grid%simulation_domain%xstop

       fld%internal%ystart = fld%grid%simulation_domain%ystart + 1
       fld%internal%ystop  = fld%grid%simulation_domain%ystop
       call gocean_stop('cv_sw_init: IMPLEMENT non-periodic BCs!')
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

    call init_periodic_bc_halos(fld)

  end subroutine cv_sw_init

  !===================================================

  subroutine cv_ne_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld

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

    if(fld%grid%boundary_conditions(1) /= BC_PERIODIC)then
      ! If we do not have periodic boundary conditions then we do
      ! not need to allow for boundary points here - they are
      ! already contained within the region.
      fld%internal%xstart = fld%grid%simulation_domain%xstart
      fld%internal%xstop  = fld%grid%simulation_domain%xstop
    else
      call gocean_stop('ERROR: cv_ne_init: implement periodic BCs!')
    end if

    if(fld%grid%boundary_conditions(2) /= BC_PERIODIC)then
      fld%internal%ystart = fld%grid%simulation_domain%ystart
      fld%internal%ystop  = fld%grid%simulation_domain%ystop - 1
    else
      call gocean_stop('ERROR: cv_ne_init: implement periodic BCs!')
    end if

  end subroutine cv_ne_init

  !===================================================

  subroutine ct_field_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld

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
    type(r2d_field), intent(inout) :: fld

    ! When updating a quantity on T points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1      i=M
    !  b  b  b  b   j=N 
    !  b  x  x  b
    !  b  x  x  b
    !  b  b  b  b   j=1

    fld%internal%xstart = fld%grid%simulation_domain%xstart
    fld%internal%xstop  = fld%grid%simulation_domain%xstop
    fld%internal%ystart = fld%grid%simulation_domain%ystart
    fld%internal%ystop  = fld%grid%simulation_domain%ystop

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

    call init_periodic_bc_halos(fld)

  end subroutine ct_sw_init

  !===================================================

  subroutine ct_ne_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld

    ! When updating a quantity on T points with a NE staggering
    ! we write to (using x to indicate a location that is written):
    !
    ! i=1          Nx
    !  b  b  b  b  b Ny
    !  b  x  x  x  b
    !  b  x  x  x  b 
    !  b  x  x  x  b
    !  b  b  b  b  b j=1

    if(fld%grid%boundary_conditions(1) /= BC_PERIODIC)then
      ! If we do not have periodic boundary conditions then we do
      ! not need to allow for boundary points here - they are
      ! already contained within the region.
      ! Start and stop are just the same as those calculated from the T mask
      ! earlier because this is a field on T points.
      fld%internal%xstart = fld%grid%simulation_domain%xstart
      fld%internal%xstop  = fld%grid%simulation_domain%xstop
    else
      call gocean_stop('ERROR: ct_ne_init: implement periodic BCs!')
    end if

    if(fld%grid%boundary_conditions(2) /= BC_PERIODIC)then
      ! Start and stop are just the same as those calculated from the T mask
      ! earlier because this is a field on T points.
      fld%internal%ystart = fld%grid%simulation_domain%ystart
      fld%internal%ystop  = fld%grid%simulation_domain%ystop
    else
      call gocean_stop('ERROR: ct_ne_init: implement periodic BCs!')
    end if

  end subroutine ct_ne_init

  !===================================================

  subroutine cf_field_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld

    fld%defined_on = F_POINTS

    select case(fld%grid%stagger)

    case(STAGGER_SW)
       call cf_sw_init(fld)

    case(STAGGER_NE)
       call cf_ne_init(fld)

    case default
       call gocean_stop('cf_field_init: ERROR - unsupported stagger!')

    end select

  end subroutine cf_field_init

  !===================================================

  subroutine cf_sw_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld

    ! When updating a quantity on F points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1         i=M
    !  o  b  b  b  b   j=N 
    !  o  b  x  x  b
    !  o  b  x  x  b
    !  o  b  b  b  b
    !  o  o  o  o  o   j=1
    if(fld%grid%boundary_conditions(1) == BC_PERIODIC)then
       fld%internal%xstart = fld%grid%simulation_domain%xstart + 1
       fld%internal%xstop  = fld%internal%xstart + fld%grid%simulation_domain%nx - 1
    else
       fld%internal%xstart = fld%grid%simulation_domain%xstart + 1
       fld%internal%xstop  = fld%grid%simulation_domain%xstop
       ! I think these are correct but we stop because I've not properly
       ! gone through the coding.
       call gocean_stop('cf_sw_init: CHECK non-periodic BCs!')
    end if

    if(fld%grid%boundary_conditions(2) == BC_PERIODIC)then
       fld%internal%ystart = fld%grid%simulation_domain%ystart + 1
       fld%internal%ystop  = fld%internal%ystart + fld%grid%simulation_domain%ny - 1
    else
       fld%internal%ystart = fld%grid%simulation_domain%ystart + 1
       fld%internal%ystop  = fld%grid%simulation_domain%ystop
       ! I think these are correct but we stop because I've not properly
       ! gone through the coding.
       call gocean_stop('cf_sw_init: CHECK non-periodic BCs!')
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

    call init_periodic_bc_halos(fld)

  end subroutine cf_sw_init

  !===================================================

  subroutine cf_ne_init(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld

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

    if(fld%grid%boundary_conditions(1) /= BC_PERIODIC)then
      ! If we do not have periodic boundary conditions then we do
      ! not need to allow for boundary points here - they are
      ! already contained within the region.
      fld%internal%xstart = fld%grid%simulation_domain%xstart
      fld%internal%xstop  = fld%grid%simulation_domain%xstop - 1
    else
      call gocean_stop('ERROR: cf_ne_init: implement periodic BCs!')
      stop
    end if

    if(fld%grid%boundary_conditions(2) /= BC_PERIODIC)then
      fld%internal%ystart = fld%grid%simulation_domain%ystart
      fld%internal%ystop  = fld%grid%simulation_domain%ystop - 1
    else
      call gocean_stop('ERROR: cf_ne_init: implement periodic BCs!')
    end if

  end subroutine cf_ne_init

  !===================================================

  SUBROUTINE copy_scalar_field(field_in, field_out)
    IMPLICIT none
    type(scalar_field), INTENT(in) :: field_in
    type(scalar_field), INTENT(out) :: field_out

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
    type(r2d_field), intent(in)    :: field_in
    type(r2d_field), intent(inout) :: field_out
        
    field_out%data(:,:) = field_in%data(:,:)
        
  end subroutine copy_2dfield

  !===================================================

  !> Copy from one patch to another in a field
  subroutine copy_2dfield_patch(field, src, dest)
    implicit none
    type(r2d_field), intent(inout) :: field
    type(region_type),    intent(in)    :: src, dest

    field%data(dest%xstart:dest%xstop,dest%ystart:dest%ystop) = &
     field%data(src%xstart:src%xstop ,src%ystart:src%ystop)
        
  end subroutine copy_2dfield_patch

  !===================================================

  subroutine increment_scalar_field(field, incr)
    implicit none
    type(scalar_field), intent(inout) :: field
    type(scalar_field), intent(in)    :: incr

    field%data = field%data + incr%data

  END SUBROUTINE increment_scalar_field

  !===================================================

  subroutine increment_scalar_field_r8(field, incr)
    implicit none
    type(scalar_field), intent(inout) :: field
    real(wp), intent(in)    :: incr

    field%data = field%data + incr

  END SUBROUTINE increment_scalar_field_r8

  !===================================================

  SUBROUTINE set_field(fld, val)
    implicit none
    class(field_type), INTENT(out) :: fld
    real(wp), INTENT(in) :: val

    select type(fld)
    type is (scalar_field)
       fld%data = val
    type is (r2d_field)
       fld%data = val
    class default
    end select

  END SUBROUTINE set_field

  !===================================================

  function field_checksum(field) result(val)
    implicit none
    type(r2d_field), intent(in) :: field
    real(wp) :: val

    val = SUM( ABS(field%data(field%internal%xstart:field%internal%xstop, &
                              field%internal%ystart:field%internal%ystop)) )

  end function field_checksum

  !===================================================

  subroutine init_periodic_bc_halos(fld)
    implicit none
    type(r2d_field), intent(inout) :: fld
    ! Locals
    integer :: ihalo

    ! Check whether we have PBCs in x AND y dimensions
    fld%num_halos = 0
    if( fld%grid%boundary_conditions(1) == BC_PERIODIC )then
       fld%num_halos = fld%num_halos + 2
    end if
    if( fld%grid%boundary_conditions(2) == BC_PERIODIC )then
       fld%num_halos = fld%num_halos + 2
    end if

    allocate( fld%halo(fld%num_halos) )

    ihalo = 0
    if( fld%grid%boundary_conditions(1) == BC_PERIODIC )then
       ! E-most column set to W-most internal column
       ihalo = ihalo + 1
       fld%halo(ihalo)%dest%xstart = fld%internal%xstop + 1
       fld%halo(ihalo)%dest%xstop  = fld%internal%xstop + 1
       fld%halo(ihalo)%dest%ystart = fld%internal%ystart   
       fld%halo(ihalo)%dest%ystop  = fld%internal%ystop

       fld%halo(ihalo)%source%xstart = fld%internal%xstart   
       fld%halo(ihalo)%source%xstop  = fld%internal%xstart
       fld%halo(ihalo)%source%ystart = fld%internal%ystart  
       fld%halo(ihalo)%source%ystop  = fld%internal%ystop

       ! W-most column set to E-most internal column
       ihalo = ihalo + 1
       fld%halo(ihalo)%dest%xstart = fld%internal%xstart-1   
       fld%halo(ihalo)%dest%xstop  = fld%internal%xstart-1
       fld%halo(ihalo)%dest%ystart = fld%internal%ystart   
       fld%halo(ihalo)%dest%ystop  = fld%internal%ystop

       fld%halo(ihalo)%source%xstart = fld%internal%xstop 
       fld%halo(ihalo)%source%xstop  = fld%internal%xstop
       fld%halo(ihalo)%source%ystart = fld%internal%ystart
       fld%halo(ihalo)%source%ystop  = fld%internal%ystop
    end if

    if( fld%grid%boundary_conditions(2) == BC_PERIODIC )then
       ! N-most row set to S-most internal row
       ihalo = ihalo + 1
       fld%halo(ihalo)%dest%xstart = fld%internal%xstart - 1   
       fld%halo(ihalo)%dest%xstop  = fld%internal%xstop  + 1
       fld%halo(ihalo)%dest%ystart = fld%internal%ystop+1   
       fld%halo(ihalo)%dest%ystop  = fld%internal%ystop+1

       fld%halo(ihalo)%source%xstart = fld%internal%xstart - 1
       fld%halo(ihalo)%source%xstop  = fld%internal%xstop  + 1
       fld%halo(ihalo)%source%ystart = fld%internal%ystart 
       fld%halo(ihalo)%source%ystop  = fld%internal%ystart

       ! S-most row set to N-most internal row
       ihalo = ihalo + 1
       fld%halo(ihalo)%dest%xstart = fld%internal%xstart - 1   
       fld%halo(ihalo)%dest%xstop  = fld%internal%xstop  + 1
       fld%halo(ihalo)%dest%ystart = fld%internal%ystart - 1   
       fld%halo(ihalo)%dest%ystop  = fld%internal%ystart - 1

       fld%halo(ihalo)%source%xstart = fld%internal%xstart - 1 
       fld%halo(ihalo)%source%xstop  = fld%internal%xstop  + 1
       fld%halo(ihalo)%source%ystart = fld%internal%ystop 
       fld%halo(ihalo)%source%ystop  = fld%internal%ystop
    end if

  end subroutine init_periodic_bc_halos

end module field_mod
