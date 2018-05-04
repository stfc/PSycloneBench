program example2
  ! no explicit dependencies so add a barrier
  implicit none
  integer :: i
  do i=1,4
     !$omp task
     call my_work1(i)
     !$omp end task
  end do
  !$omp taskwait
  do i=1,4
     !$omp task
     call my_work2(i)
     !$omp end task
  end do
end program example2

subroutine my_work1(i)
  integer :: i
  print *, "hello", i
end subroutine my_work1

subroutine my_work2(i)
  integer :: i
  print *, "world", i
end subroutine my_work2
