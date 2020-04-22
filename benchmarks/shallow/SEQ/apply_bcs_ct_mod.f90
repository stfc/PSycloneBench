!> \brief Apply boundary conditions for field on CT
!! \detail Applies cyclic boundary conditions for a 
!! field defined on CT.
module apply_bcs_ct_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  implicit none

  private

  public invoke_apply_bcs_ct
  public apply_bcs_ct_type, apply_bcs_ct_code

  type, extends(kernel_type) :: apply_bcs_ct_type
     type(arg), dimension(1) :: meta_args =    &
          (/ arg(READWRITE, CT, POINTWISE)     & ! field
           /)
     !> This kernel writes only to internal points of the
     !! simulation domain.
     integer :: ITERATES_OVER = INTERNAL_PTS
  contains
    procedure, nopass :: code => apply_bcs_ct_code
  end type apply_bcs_ct_type

contains

  !===================================================

  !> Manual implementation of the code needed to invoke
  !! apply_bcs_ct_code().
  subroutine invoke_apply_bcs_ct(fld)
    use field_mod
    implicit none
    type(r2d_field), intent(inout) :: fld
    ! Locals
    integer :: ihalo

    ! Note that we do not loop over the full extent of the field.
    ! Fields are allocated with extents (M+1,N+1).
    ! Presumably the extra row and column are needed for periodic BCs.
    ! We are updating the last row and the last column of a quantity on CT.
    ! This looks like (using x to indicate a location that is written first
    ! and y a location that is written second):
    !
    !  i=1   i=M
    ! _ y  y  y  y 
    ! /|o  o  o  x   j=N 
    ! | o  o  o  x
    ! \ o  o  o  x
    !  \o  o  o _x   j=1
    !    \______/|

    !   vi-1j+1--fij+1---vij+1---fi+1j+1-vi+1j+1
    !   |        |       |       |        |    
    !   |        |       |       |        |    
    !   Ti-1j----uij-----Tij-----ui+1j --Ti+1j--
    !   |        |       |       |        |    
    !   |        |       |       |        |    
    !   vi-1j----fij-----vij-----fi+1j --vi+1j--
    !   |        |       |       |        |    
    !   |        |       |       |        |    
    !   Ti-1j-1--uij-1---Tij-1---ui+1j-1-Ti+1j-1
    !

!DIR$ LOOP_INFO max_trips(2)
    do ihalo = 1, fld%num_halos

      ! Copy from source to destination
      call copy_field(fld, fld%halo(ihalo)%source, fld%halo(ihalo)%dest)
  
    end do

  end subroutine invoke_apply_bcs_ct

  !===================================================

  !> Apply cyclic boundary conditions to field on CT
  subroutine apply_bcs_ct_code(n, mp1, np1, field)
    implicit none
    integer,  intent(in) :: n, mp1, np1
    real(wp), intent(inout), dimension(:,:) :: field

    ! Last col = first col
    field(MP1,1:N) = field(1,  1:N)
    ! Last row = first row
    field(1:MP1,NP1) = field(1:MP1,1)

  end subroutine apply_bcs_ct_code

end module apply_bcs_ct_mod
