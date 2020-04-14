module psy_layer
contains
subroutine psy(a)
  use field_mod, only : field
  implicit none
  type(field) :: a

  call psy_target(a%data)

end subroutine psy
end module psy_layer
