module psy_layer
contains
subroutine psy(a)
  ! full field dependencies with gocean data types and doubly nested
  ! loops
  use field_mod, only : field
  use kernel_work1, only : work1
  use kernel_work2, only : work2
  implicit none
  integer :: i, j
  type(field) :: a

  a%data = 0.0

  do j=1,2
     do i=1,2
        !$omp task out(a%data(i,j))
        call work1(i,j,a%data)
        !$omp end task
     end do
  end do

  do j=1,2
     do i=1,2
        !$omp task in(a%data(i,j))
        call work2(i,j,a%data)
        !$omp end task
     end do
  end do

  !$omp taskwait

end subroutine psy
end module psy_layer
