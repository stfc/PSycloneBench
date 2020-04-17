program example3
  ! array dependencies. All of the first loop must complete before
  ! starting the second loop
  implicit none
  integer :: i
  real :: a(4)
  a = 0.0
  do i=1,4
     !$omp task out(a)
     call my_work1(i,a)
     !$omp end task
  end do

  do i=1,4
     !$omp task in(a)
     call my_work2(i,a)
     !$omp end task
  end do
end program example3

subroutine my_work1(i,a)
  implicit none
  real :: a(4)
  integer :: i
  print *, "hello", i
  a(i)=i
end subroutine my_work1

subroutine my_work2(i,a)
  implicit none
  real :: a(i)
  integer :: i
  print *, "world", i,a(i)
end subroutine my_work2
