!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
! LICENCE.original is available from the Met Office Science Repository Service:
! https://code.metoffice.gov.uk/trac/lfric/browser/LFRic/trunk/LICENCE.original
! -----------------------------------------------------------------------------
! BSD 3-Clause License
!
! Modifications copyright (c) 2020, Science and Technology Facilities Council
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! * Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!
! * Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
!
! * Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! -----------------------------------------------------------------------------
! Modified by W. Hayek, NIWA, New Zealand
! Modified by R. W. Ford, STFC Daresbury Lab

! A copy of the LFRic matrix vector kernel with the kernel metadata
! (and associated use statements) removed, comments describing the
! interface removed and the code restructured by replacing the matmul
! intrinsic with an equivalent libxsmm call and by adding a
! just-in-time libxsmm dispatch call to the start of the kernel.

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
