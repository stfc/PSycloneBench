program example1
  ! hello world independent tasks
  implicit none
  integer :: i
  do i=1,4
     !$omp task
     print *, "hello world", i
     !$omp end task
  end do
end program example1
