!> Module to encapsulate the topology of the staggered C grid
!! used in the finite difference scheme.
module topology_mod
  IMPLICIT none

  private

  type :: topology_type

     integer :: istart, istop
     integer :: jstart, jstop
  end type topology_type

  type(topology_type) :: cu, cv, ct, cf

  public cu, cv, ct, cf
  public topology_init

contains

  !> Routine to set-up information on where the different
  !! mesh point types (U,V,T,F) sit on the computational mesh.
  subroutine topology_init(nx, ny)
    implicit none
    ! Extent in x and y of the computational mesh. 
    integer, intent(in) :: nx, ny

    ! When updating a quantity on U points we write to:
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  o  x  x  x   j=N
    !  o  x  x  x
    !  o  x  x  x   j=1
    cu%istart = 2
    cu%istop  = nx
    cu%jstart = 1
    cu%jstop  = ny-1

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
