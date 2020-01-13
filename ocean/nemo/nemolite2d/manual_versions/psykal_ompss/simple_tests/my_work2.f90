module kernel_work2
  contains
subroutine work2(i,j,a)
  implicit none
  real :: a(2,2)
  integer :: i,j
  print *, "world", i,j,a(i,j)
end subroutine work2
end module kernel_work2
