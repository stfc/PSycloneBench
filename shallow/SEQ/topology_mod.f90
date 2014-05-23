!> Module to encapsulate the topology of the staggered Arakawa C grid
!! used in the finite difference scheme.
!! Arrays are allocated with extents (M+1,N+1) but our fields are really
!! only MxN.
module topology_mod
  implicit none

  private

  !> Specify a region on the simulation grid
  type :: region
     integer :: istart, istop
     integer :: jstart, jstop
  end type region

  !> Specify the source and destination regions
  !! of a halo
  type :: halo_type
     type(region) :: src
     type(region) :: dest
  end type halo_type

  !> Specifies the region of the grid occupied by the
  !! grid-point type and the associated halos
  type, extends(region) :: topology_type
     !> No. of halo regions associated with this field
     integer :: nhalos
     !> Array of halo region definitions
     type(halo_type), allocatable, dimension(:) :: halo
  end type topology_type

  type(topology_type) :: cu, cv, ct, cf

  public cu, cv, ct, cf
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
    cu%istart = 2
    cu%istop  = M+1
    cu%jstart = 1
    cu%jstop  = N

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

    cu%nhalos = 2
    ALLOCATE(cu%halo(cu%nhalos))

    cu%halo(1)%dest%istart = 1   ; cu%halo(1)%dest%istop = 1
    cu%halo(1)%dest%jstart = 1   ; cu%halo(1)%dest%jstop = N
    cu%halo(1)%src%istart  = M+1 ; cu%halo(1)%src%istop  = M+1
    cu%halo(1)%src%jstart  = 1   ; cu%halo(1)%src%jstop  = N

    cu%halo(2)%dest%istart = 1   ; cu%halo(2)%dest%istop = M+1
    cu%halo(2)%dest%jstart = N+1 ; cu%halo(2)%dest%jstop = N+1
    cu%halo(2)%src%istart  = 1   ; cu%halo(2)%src%istop  = M+1
    cu%halo(2)%src%jstart  = 1   ; cu%halo(2)%src%jstop  = 1


    ! When updating a quantity on V points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  x  x  x  o 
    !  x  x  x  o   j=N
    !  x  x  x  o
    !  o  o  o  o   j=1
    cv%istart = 1
    cv%istop  = nx-1
    cv%jstart = 2
    cv%jstop  = ny

    ! When updating a quantity on T points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  x  x  x  o   j=N
    !  x  x  x  o
    !  x  x  x  o   j=1
    ct%istart = 1
    ct%istop  = nx-1
    ct%jstart = 1
    ct%jstop  = ny-1

    ! When updating a quantity on F points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  x  x  x 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  o  o  o   j=1
    cf%istart = 2
    cf%istop  = nx
    cf%jstart = 2
    cf%jstop  = ny

  end subroutine topology_init

end module topology_mod
