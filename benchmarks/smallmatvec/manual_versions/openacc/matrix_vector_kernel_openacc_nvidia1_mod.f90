module matrix_vector_kernel_mod

use constants_mod,           only : r_def, i_def

implicit none

contains

subroutine matrix_vector_kernel_code(cell,              &
                                     nlayers,           &
                                     lhs, x,            &
                                     ncell_3d,          &
                                     matrix,            &
                                     ndf1, undf1, map1, &
                                     ndf2, undf2, map2)

!$acc routine vector
  implicit none

  ! Arguments
  integer(kind=i_def),                  intent(in) :: cell, nlayers, ncell_3d
  integer(kind=i_def),                  intent(in) :: undf1, ndf1
  integer(kind=i_def),                  intent(in) :: undf2, ndf2
  integer(kind=i_def), dimension(ndf1), intent(in) :: map1
  integer(kind=i_def), dimension(ndf2), intent(in) :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  ! Internal variables
  integer(kind=i_def)               :: df, k, ik
  real(kind=r_def), dimension(ndf2) :: x_e
  real(kind=r_def) :: lhs_e
  integer :: i, j

  ik = (cell-1)*nlayers + 1

  !$acc loop vector independent
  do k = 0, nlayers-1
    do i = 0, ndf1, 1
      lhs_e=0.0_r_def
      do j = 1, ndf2, 1
        lhs_e=lhs_e + matrix(i,j,ik + k) * x(map2(j)+k)
      enddo
      !$acc atomic update
      lhs(map1(i)+k) = lhs(map1(i)+k) + lhs_e
    end do
  end do

end subroutine matrix_vector_kernel_code

end module matrix_vector_kernel_mod
