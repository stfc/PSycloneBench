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
  integer(kind=i_def)               :: df1, df2, k, ik

    ik = (cell-1)*nlayers

    do df1 = 1, ndf1
        do df2 = 1, ndf2
            !$OMP SIMD
            do k = 1, nlayers
                lhs(map1(df1)+k-1) = lhs(map1(df1)+k-1) + matrix(ik+k,df1,df2) * x(map2(df2)+k-1)
            end do
        end do
    end do

end subroutine matrix_vector_kernel_code

end module matrix_vector_kernel_mod
