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

  use libxsmm, only : libxsmm_dmmfunction, &
                      libxsmm_dispatch, libxsmm_prefetch_none, &
                      libxsmm_dmmcall, libxsmm_gemm_flag_none

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
  real(kind=r_def), dimension(ndf1) :: lhs_e

  ! libxmm matrix multiplication kernel handle
  type(libxsmm_dmmfunction) :: xmm

  ! Prepare kernel outside loop to reduce overhead
  ! Available prefetch strategies:
  ! libxsmm_prefetch_none
  ! libxsmm_prefetch_auto
  ! libxsmm_prefetch_al1
  ! libxsmm_prefetch_al2
  call libxsmm_dispatch(xmm, ndf1, 1, ndf2, ndf1, ndf2, ndf1, 1.0_r_def, &
                        0.0_r_def, libxsmm_gemm_flag_none, libxsmm_prefetch_none)

  do k = 0, nlayers-1
    do df = 1, ndf2
      x_e(df) = x(map2(df)+k)
    end do
    ik = (cell-1)*nlayers + k + 1
    ! lhs_e = matmul(matrix(:,:,ik),x_e)
    call libxsmm_dmmcall(xmm, matrix(:,:,ik), x_e, lhs_e)
    do df = 1,ndf1
       lhs(map1(df)+k) = lhs(map1(df)+k) + lhs_e(df)
    end do
  end do

end subroutine matrix_vector_kernel_code

end module matrix_vector_kernel_mod
