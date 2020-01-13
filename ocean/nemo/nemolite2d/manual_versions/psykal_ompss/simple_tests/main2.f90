subroutine algorithm
  ! example algorithm code
  use field_mod, only : field
  use psy_layer, only : psy
  implicit none
  type(field) :: a

  call psy(a)

end subroutine algorithm
