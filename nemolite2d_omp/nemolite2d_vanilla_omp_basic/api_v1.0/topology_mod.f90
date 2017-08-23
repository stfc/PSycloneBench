!> Module to encapsulate the topology of the staggered Arakawa C grid
!! used in the finite difference scheme.
!! Arrays are allocated with extents (M+1,N+1) but our fields are really
!! only MxN.
module topology_mod
  implicit none

  private

  !> Specifies the region of the grid occupied by the
  !! grid-point type and the associated halos
  type, extends(region) :: topology_type
     !> No. of halo regions associated with this field
     integer :: nhalos
     !> Array of halo region definitions
     type(halo_type), allocatable, dimension(:) :: halo
  end type topology_type

  type(topology_type) :: cu_grid, cv_grid, ct_grid, cf_grid

  public region
  public topology_type
  public cu_grid, cv_grid, ct_grid, cf_grid
  public topology_init

contains

  !> Routine to set-up information on where the different
  !! mesh point types (U,V,T,F) sit on the computational mesh.
  subroutine topology_init(nx, ny)
    implicit none
    !> Extent in x and y of the computational mesh (arrays).
    !! Fields sit on a sub-array within these extents.
    integer, intent(in) :: nx, ny
    ! Locals
    integer :: M, N

    M = nx - 1
    N = ny - 1

    ! When updating a quantity on U points we write to:
    ! (using 'x' to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  x  x  x   j=1
    cu_grid%istart = 2
    cu_grid%istop  = M+1
    cu_grid%jstart = 1
    cu_grid%jstop  = N

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

    cu_grid%nhalos = 2
    ALLOCATE(cu_grid%halo(cu_grid%nhalos))

    cu_grid%halo(1)%dest%istart = 1   ; cu_grid%halo(1)%dest%istop = 1
    cu_grid%halo(1)%dest%jstart = 1   ; cu_grid%halo(1)%dest%jstop = N
    cu_grid%halo(1)%src%istart  = M+1 ; cu_grid%halo(1)%src%istop  = M+1
    cu_grid%halo(1)%src%jstart  = 1   ; cu_grid%halo(1)%src%jstop  = N

    cu_grid%halo(2)%dest%istart = 1   ; cu_grid%halo(2)%dest%istop = M+1
    cu_grid%halo(2)%dest%jstart = N+1 ; cu_grid%halo(2)%dest%jstop = N+1
    cu_grid%halo(2)%src%istart  = 1   ; cu_grid%halo(2)%src%istop  = M+1
    cu_grid%halo(2)%src%jstart  = 1   ; cu_grid%halo(2)%src%jstop  = 1


    ! When updating a quantity on V points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  x  x  x  o 
    !  x  x  x  o   j=N
    !  x  x  x  o
    !  o  o  o  o   j=1
    cv_grid%istart = 1
    cv_grid%istop  = nx-1
    cv_grid%jstart = 2
    cv_grid%jstop  = ny

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

    cv_grid%nhalos = 2
    ALLOCATE(cv_grid%halo(cv_grid%nhalos))

    cv_grid%halo(1)%dest%istart = 1   ; cv_grid%halo(1)%dest%istop = M
    cv_grid%halo(1)%dest%jstart = 1   ; cv_grid%halo(1)%dest%jstop = 1
    cv_grid%halo(1)%src%istart  = 1   ; cv_grid%halo(1)%src%istop  = M
    cv_grid%halo(1)%src%jstart  = N+1 ; cv_grid%halo(1)%src%jstop  = N+1

    cv_grid%halo(2)%dest%istart = M+1 ; cv_grid%halo(2)%dest%istop = M+1
    cv_grid%halo(2)%dest%jstart = 1   ; cv_grid%halo(2)%dest%jstop = N+1
    cv_grid%halo(2)%src%istart  = 1   ; cv_grid%halo(2)%src%istop  = 1
    cv_grid%halo(2)%src%jstart  = 1   ; cv_grid%halo(2)%src%jstop  = N+1

    ! When updating a quantity on T points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  x  x  x  o   j=N
    !  x  x  x  o
    !  x  x  x  o   j=1
    ct_grid%istart = 1
    ct_grid%istop  = nx-1
    ct_grid%jstart = 1
    ct_grid%jstop  = ny-1

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

    ct_grid%nhalos = 2
    ALLOCATE( ct_grid%halo(ct_grid%nhalos) )

    ct_grid%halo(1)%dest%istart = M+1 ; ct_grid%halo(1)%dest%istop = M+1
    ct_grid%halo(1)%dest%jstart = 1   ; ct_grid%halo(1)%dest%jstop = N
    ct_grid%halo(1)%src%istart  = 1   ; ct_grid%halo(1)%src%istop  = 1
    ct_grid%halo(1)%src%jstart  = 1   ; ct_grid%halo(1)%src%jstop  = N

    ct_grid%halo(2)%dest%istart = 1   ; ct_grid%halo(2)%dest%istop = M+1
    ct_grid%halo(2)%dest%jstart = N+1 ; ct_grid%halo(2)%dest%jstop = N+1
    ct_grid%halo(2)%src%istart  = 1   ; ct_grid%halo(2)%src%istop  = M+1
    ct_grid%halo(2)%src%jstart  = 1   ; ct_grid%halo(2)%src%jstop  = 1

    ! When updating a quantity on F points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  x  x  x 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  o  o  o   j=1
    cf_grid%istart = 2
    cf_grid%istop  = nx
    cf_grid%jstart = 2
    cf_grid%jstop  = ny

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

    cf_grid%nhalos = 2
    ALLOCATE( cf_grid%halo(cf_grid%nhalos) )

    cf_grid%halo(1)%dest%istart = 1   ; cf_grid%halo(1)%dest%istop = 1
    cf_grid%halo(1)%dest%jstart = 2   ; cf_grid%halo(1)%dest%jstop = N+1
    cf_grid%halo(1)%src%istart  = M+1 ; cf_grid%halo(1)%src%istop  = M+1
    cf_grid%halo(1)%src%jstart  = 2   ; cf_grid%halo(1)%src%jstop  = N+1

    cf_grid%halo(2)%dest%istart = 1   ; cf_grid%halo(2)%dest%istop = M+1
    cf_grid%halo(2)%dest%jstart = 1   ; cf_grid%halo(2)%dest%jstop = 1
    cf_grid%halo(2)%src%istart  = 1   ; cf_grid%halo(2)%src%istop  = M+1
    cf_grid%halo(2)%src%jstart  = N+1 ; cf_grid%halo(2)%src%jstop  = N+1

  end subroutine topology_init

end module topology_mod
