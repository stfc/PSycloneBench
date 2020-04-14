subroutine psy_target(a)
  ! full field dependencies with gocean data types and doubly nested
  ! loops
  use kernel_work1, only : work1
  use kernel_work2, only : work2
  implicit none
  integer :: i, j
  real :: a(2,2)

  a = 0.0

  do j=1,2
     do i=1,2
        !$omp task out(a(i,j))
        call work1(i,j,a)
        !$omp end task
     end do
  end do

  do j=1,2
     do i=1,2
        !$omp task in(a(i,j))
        call work2(i,j,a)
        !$omp end task
     end do
  end do

  !$omp taskwait

end subroutine psy_target
