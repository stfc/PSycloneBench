program example8
  ! index field dependencies with gocean data types and doubly nested
  ! loops and PSy code split into two layers
  implicit none
  integer :: i, j
  type field
     real :: data(2,2)
  end type field
  type(field), target :: a

  real, dimension(:,:), pointer :: array_a

  array_a => a%data

  call control(array_a)

end program example8

subroutine control(array_a)

  real :: array_a(2,2)

  a = 0.0
    
  do j=1,2
     do i=1,2
        !$omp task out(array_a(i,j))
        call my_work1(i,j,array_a)
        !$omp end task
     end do
  end do

  do j=1,2
     do i=1,2
        !$omp task in(array_a(i,j))
        call my_work2(i,j,array_a)
        !$omp end task
     end do
  end do

  ! we get a seg fault without this taskwait at the end
  !$omp taskwait

end subroutine control

subroutine my_work1(i,j,a)
  implicit none
  real :: a(2,2)
  integer :: i,j
  print *, "hello", i,j
  a(i,j)=i+j
end subroutine my_work1

subroutine my_work2(i,j,a)
  implicit none
  real :: a(2,2)
  integer :: i,j
  print *, "world", i,j,a(i,j)
end subroutine my_work2
