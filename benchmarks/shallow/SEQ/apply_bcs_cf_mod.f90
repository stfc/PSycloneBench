!> \brief Apply boundary conditions for field on CF
!! \detail Applies cyclic boundary conditions for a 
!! field defined on CF.
module apply_bcs_cf_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  implicit none

  private

  public invoke_apply_bcs_cf
  public apply_bcs_cf_type, apply_bcs_cf_code

  type, extends(kernel_type) :: apply_bcs_cf_type
     type(arg), dimension(1) :: meta_args =    &
          (/ arg(READWRITE, CF, POINTWISE)     & ! field
           /)
     !> This kernel writes only to internal points of the
     !! simulation domain.
     integer :: ITERATES_OVER = INTERNAL_PTS
  contains
    procedure, nopass :: code => apply_bcs_cf_code
  end type apply_bcs_cf_type

contains

  !===================================================

  !> Manual implementation of the code needed to invoke
  !! apply_bcs_cf_code().
  subroutine invoke_apply_bcs_cf(field)
    use field_mod
    implicit none
    type(r2d_field), intent(inout) :: field
    ! Locals
    integer :: ihalo

    ! Note that we do not loop over the full extent of the array.
    ! Arrays are allocated with extents (M+1,N+1) but our field is really
    ! only MxN. The CF field sits at (2:M+1,2:N+1) in the array - this offset
    ! is due to the staggering of the fields on the Arakawa C grid.
    ! The extra row and column are used for periodic BCs.
    ! We are updating the first column and last row of a quantity on CF.
    ! This looks like (using x to indicate a location that is written first
    ! and y a location that is written second):
    !
    !  i=1   i=M
    ! .-x  o  o  o 
    ! | x  o  o  o   j=N 
    ! | x  o  o  o
    ! | x  o  o  o
    ! ->y_ y  y  y   j=1
    !   |\______/ 

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
    do ihalo = 1, field%num_halos

      ! Copy from source to destination
      call copy_field(field, &
                      field%halo(ihalo)%source, &
                      field%halo(ihalo)%dest)

    end do

  end subroutine invoke_apply_bcs_cf

  !===================================================

  !> Apply cyclic boundary conditions to field on CF
  subroutine apply_bcs_cf_code(mp1, np1, field)
    implicit none
    integer,  intent(in) :: mp1, np1
    real(wp), intent(inout), dimension(:,:) :: field

    ! First col = last col
    field(1,    2:NP1) = field(MP1,  2:NP1)
    ! First row = last row
    field(1:MP1,1)     = field(1:MP1,NP1)

  end subroutine apply_bcs_cf_code

end module apply_bcs_cf_mod
