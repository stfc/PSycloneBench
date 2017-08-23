module ssh_update_mod

  type, extends(kernel_type) :: sshn
     type(arg), dimension(3) :: meta_args =    &
          (/ arg(WRITE, CV, POINTWISE),        & ! cv
             arg(READ,  CT, POINTWISE),        & ! p
             arg(READ,  CV, POINTWISE)         & ! v
           /)
     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS
  contains
    procedure, nopass :: code => sshn_code
  end type sshn

contains

  subroutine sshn_code

  end subroutine sshn_code

end module ssh_update_mod
