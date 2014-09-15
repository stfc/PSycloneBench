module compute_pnew_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use field_mod
  implicit none

  private

  public manual_invoke_compute_pnew
  public compute_pnew_type, compute_pnew_code

  TYPE, EXTENDS(kernel_type) :: compute_pnew_type
     TYPE(arg), DIMENSION(5) :: meta_args =    &
          (/ arg(WRITE, CT, POINTWISE),        & ! pnew
             arg(READ,  CT, POINTWISE),        & ! pold
             arg(READ,  CU, POINTWISE),        & ! cu
             arg(READ,  CV, POINTWISE),        & ! cv
             arg(READ,  R,  POINTWISE)         & ! tdt
           /)
     !> This kernel operates on fields that live on an
     !! orthogonal, regular grid.
     integer :: GRID_TYPE = ORTHOGONAL_REGULAR

     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     INTEGER :: ITERATES_OVER = DOFS

     !> This kernel is written assuming that the internal
     !! regions for each grid-point type begin at the
     !! following locations on the grid:
     integer :: tstart(2) = (/ 1, 1 /)
     integer :: ustart(2) = (/ 2, 1 /)
     integer :: vstart(2) = (/ 1, 2 /)
     integer :: fstart(2) = (/ 2, 2 /)
  CONTAINS
    procedure, nopass :: code => compute_pnew_code
  END TYPE compute_pnew_type

contains

  !===================================================

  subroutine manual_invoke_compute_pnew(pnew, pold, cu, cv, tdt)
    implicit none
    type(r2d_field_type), intent(inout) :: pnew
    type(r2d_field_type), intent(in)    :: pold, cu, cv
    real(wp), intent(in) :: tdt
    ! Locals
    integer :: I, J
    real(wp) :: dx, dy
    integer, dimension(2) :: kern_ushift, kern_vshift
    integer, dimension(2) :: ushift, vshift
    type(compute_pnew_type) :: this_kernel

    ! Note that we do not loop over the full extent of the field.
    ! Fields are allocated with extents (M+1,N+1).
    ! Presumably the extra row and column are needed for periodic BCs.
    ! We are updating a quantity on CT.
    ! This loop writes to pnew(1:M,1:N) so this looks like
    ! (using x to indicate a location that is written):
    !
    ! i=1   i=M
    !  o  o  o  o 
    !  x  x  x  o   j=N
    !  x  x  x  o
    !  x  x  x  o   j=1

    ! Original code looked like:
    !
    !     DO J=1,N
    !        DO I=1,M
    !           PNEW(I,J) = POLD(I,J)-TDTSDX*(CU(I+1,J)-CU(I,J))   & 
    !               -TDTSDY*(CV(I,J+1)-CV(I,J))
    !        END DO
    !     END DO

    ! pnew(i,j) depends upon:
    !   pold(i,j)                                : CT
    !   cu(i,j), cu(i+1,j)                       : CU
    !    => lateral CU neighbours of the CT pt being updated 
    !   cv(i,j), cv(i,j+1)                       : CT
    !    => vertical CV neighbours of the CT pt being updated

    !   x-------vij+1---fi+1j+1
    !   |       |       |
    !   |       |       |
    !   uij-----Tij-----ui+1j
    !   |       |       |
    !   |       |       |
    !   fij-----vij-----fi+1j
    !   |       |       |
    !   |       |       |
    !   uij-1- -Tij-1---ui+1j-1
    !
    dx = pnew%grid%dx
    dy = pnew%grid%dy

    ! The shifts from T to xx pts assumed by this kernel
    kern_ushift(1) = this_kernel%ustart(1) - this_kernel%tstart(1)
    kern_ushift(2) = this_kernel%ustart(2) - this_kernel%tstart(2)
    kern_vshift(1) = this_kernel%vstart(1) - this_kernel%tstart(1)
    kern_vshift(2) = this_kernel%vstart(2) - this_kernel%tstart(2)

    ! The shifts we must apply to account for the fact that our
    ! fields may not have the same relative positioning as
    ! assumed by the kernel
    ushift(1) = cu%internal%xstart - pnew%internal%xstart - kern_ushift(1)
    ushift(2) = cu%internal%ystart - pnew%internal%ystart - kern_ushift(2)
    vshift(1) = cv%internal%xstart - pnew%internal%xstart - kern_vshift(1)
    vshift(2) = cv%internal%ystart - pnew%internal%ystart - kern_vshift(2)

    DO J=pnew%internal%ystart, pnew%internal%ystop, 1
       DO I=pnew%internal%xstart, pnew%internal%xstop, 1

          CALL compute_pnew_code(i, j, dx, dy, &
                                 pnew%data, pold%data, &
                                 cu%data, cv%data, tdt)
       END DO
    END DO

  END SUBROUTINE manual_invoke_compute_pnew

  !===================================================

  subroutine compute_pnew_code(i, j, dx, dy, &
                               pnew, pold, cu, cv, tdt)
    implicit none
    integer,  intent(in) :: I, J
    real(wp), intent(in) :: dx, dy
    real(wp), intent(out), dimension(:,:) :: pnew
    real(wp), intent(in),  dimension(:,:) :: pold, cu, cv
    real(wp), intent(in) :: tdt
    ! Locals
    real(wp) :: tdtsdx, tdtsdy

    !> These quantities are computed here because tdt is not
    !! constant. (It is == dt for first time step, 2xdt for
    !! all remaining time steps.)
    tdtsdx = tdt/dx
    tdtsdy = tdt/dy

    PNEW(I,J) = POLD(I,J)-TDTSDX*(CU(I,J)-CU(I-1,J))   & 
                         -TDTSDY*(CV(I,J)-CV(I,J-1))

  END SUBROUTINE compute_pnew_code

END MODULE compute_pnew_mod
