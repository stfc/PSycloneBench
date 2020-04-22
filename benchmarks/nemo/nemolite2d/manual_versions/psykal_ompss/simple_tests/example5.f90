program example5
  ! full field dependencies with gocean data types and doubly nested
  ! loops
  implicit none
  integer :: i, j
  type field
     real :: data(2,2)
  end type field
  type(field) :: a
  a%data = 0.0
  do j=1,2
     do i=1,2
        !$omp task out(a)
        call my_work1(i,j,a%data)
        !$omp end task
     end do
  end do

  do j=1,2
     do i=1,2
        !$omp task in(a)
        call my_work2(i,j,a%data)
        !$omp end task
     end do
  end do
end program example5

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
