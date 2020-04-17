module kernel_work1
  contains
subroutine work1(i,j,a)
  implicit none
  real :: a(2,2)
  integer :: i,j
  print *, "hello", i,j
  a(i,j)=i+j
end subroutine work1
end module
