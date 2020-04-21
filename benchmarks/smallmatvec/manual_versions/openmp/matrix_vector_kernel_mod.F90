!-----------------------------------------------------------------------------
!! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
!! For further details please refer to the file LICENCE.MetOffice which you
!! should have received as part of this distribution.
!!-----------------------------------------------------------------------------

#ifndef INNERREPS
#define INNERREPS 1
#endif

module matrix_vector_kernel_mod
use argument_mod,            only : arg_type,                               &
                                    GH_FIELD, GH_OPERATOR, GH_READ, GH_INC, &
                                    ANY_SPACE_1, ANY_SPACE_2,               &
                                    CELLS 
use constants_mod,           only : r_def
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: matrix_vector_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,    GH_INC,  ANY_SPACE_1),                    &  
       arg_type(GH_FIELD,    GH_READ, ANY_SPACE_2),                    &
       arg_type(GH_OPERATOR, GH_READ, ANY_SPACE_1, ANY_SPACE_2)        &
       /)
  integer :: iterates_over = CELLS
contains
!procedure, nopass ::matrix_vector_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface matrix_vector_kernel_type
   module procedure matrix_vector_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!public matrix_vector_code
contains

  type(matrix_vector_kernel_type) function matrix_vector_kernel_constructor() result(self)
  return
end function matrix_vector_kernel_constructor

!> @brief Computes lhs = matrix*x
!> @param[in] cell Horizontal cell index
!! @param[in] nlayers Number of layers
!! @param[in] ncell_3d Total number of cells
!! @param[in] ndf1 Number of degrees of freedom per cell for the output field
!! @param[in] undf1 Unique number of degrees of freedom  for the output field
!! @param[in] map1 Dofmap for the cell at the base of the column for the output field
!! @param[in] map2 Dofmap for the cell at the base of the column for the input field
!! @param[in] ndf2 Number of degrees of freedom per cell for the input field
!! @param[in] undf2 Unique number of degrees of freedom for the input field 
!! @param[in] x Input data
!> @param[inout] lhs Output lhs (A*x)
!! @param[in] matrix Local matrix assembly form of the operator A 
! Original implementation as found in LFRic, it uses the matmul intrinsic
subroutine matrix_vector_code_original(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, k, ik , df2
  real(kind=r_def), dimension(ndf2) :: x_e
  real(kind=r_def), dimension(ndf1) :: lhs_e
 
  do k = 0, nlayers-1
    do df = 1, ndf2  
       x_e(df) = x(map2(df)+k)
    end do

    do df = 1, ndf1
       lhs_e(df) = lhs(map1(df)+k)
    end do

    ik = (cell-1)*nlayers + k + 1

    lhs_e = matmul(matrix(:,:,ik),x_e)
    
    do df = 1,ndf1
       lhs(map1(df)+k) = lhs(map1(df)+k) + lhs_e(df) 
    end do
 end do
 
end subroutine matrix_vector_code_original

! This implementation fuses the gather-computation-scatter loops. The matrix
! multiplication operation is not done in a single step but interleaving the
! row-column operations of all matmuls that use the same matrix. This is done
! to minimize memory operations.
subroutine matrix_vector_code_fuseinterleaved(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, k, ik , df2

    ik = (cell-1)*nlayers

    do k = 1, nlayers
        do df = 1, ndf1
            do df2 = 1, ndf2
                lhs(map1(df)+k-1) = lhs(map1(df)+k-1) + matrix(df,df2,ik+k) * x(map2(df2)+k-1)
            end do
        end do
    end do

end subroutine matrix_vector_code_fuseinterleaved

! This implementations takes the fuseinterleaved and moves the k loop (layers)
! to the inner location. It also adds a SIMD pragma for this loop.
subroutine matrix_vector_code_kinner(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, k, ik , df2

    ik = (cell-1)*nlayers

    do df = 1, ndf1
        do df2 = 1, ndf2
            !$OMP SIMD
            do k = 1, nlayers
                lhs(map1(df)+k-1) = lhs(map1(df)+k-1) + matrix(df,df2,ik+k) * x(map2(df2)+k-1)
            end do
        end do
    end do

end subroutine matrix_vector_code_kinner

! This implementation uses the same loop order from the kinner implementation
! but expects a data-layout that matches the traversal order (with layers
! data-points in the contiguous dimension).
subroutine matrix_vector_code_nlayersf(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(nlayers,ndf2,ndf1,ncell_3d/nlayers), intent(in)    :: matrix
  real(kind=r_def) :: temp

  !Internal variables
  integer                           :: df, k, df2, m1,m2

    do df = 1, ndf1
        m1 = map1(df)
        do df2 = 1, ndf2
            m2 = map2(df2)
            !$OMP SIMD
            do k = 1, nlayers
                lhs(m1+k-1) = lhs(m1+k-1) + matrix(k,df2,df,cell) * x(m2+k-1)
            end do
        end do
    end do

end subroutine matrix_vector_code_nlayersf

! Similar to nlayersf but 
subroutine matrix_vector_code_nlayersf2(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(nlayers,ndf2,ndf1,ncell_3d/nlayers), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, k, df2, m1,m2
  real(kind=r_def), dimension(nlayers, ndf2)   :: tempx

    do df2 = 1, ndf2
        m2 = map2(df2)
        do k = 1, nlayers
            tempx(k,df2) = x(m2+k-1)
        enddo
    enddo

    do df = 1, ndf1
        m1 = map1(df)
        do df2 = 1, ndf2
            !$OMP SIMD
            do k = 1, nlayers
                lhs(m1+k-1) = lhs(m1+k-1) + matrix(k,df2,df,cell) * tempx(k,df2)
            end do
        end do
    end do

end subroutine matrix_vector_code_nlayersf2


subroutine matrix_vector_code_nlayersf_atomics(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(nlayers,ndf2,ndf1,ncell_3d/nlayers), intent(in)    :: matrix
  real(kind=r_def) :: temp

  !Internal variables
  integer                           :: df, k, df2, m1,m2

    do df = 1, ndf1
        m1 = map1(df)
        do df2 = 1, ndf2
            m2 = map2(df2)
            !$OMP SIMD
            do k = 1, nlayers
                !$OMP ATOMIC
                lhs(m1+k-1) = lhs(m1+k-1) + matrix(k,df2,df,cell) * x(m2+k-1)
            end do
        end do
    end do

end subroutine matrix_vector_code_nlayersf_atomics


! Similar to nlayersf but repeats the inner loop INNERREPS times in order to
! simulate the performance of multiple fused matrix multiplications if those
! could be fused up to the inner loop.
subroutine matrix_vector_code_nlayersf_moreops(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(nlayers,ndf2,ndf1,ncell_3d/nlayers), intent(in)    :: matrix
  real(kind=r_def) :: temp

  !Internal variables
  integer :: df, k, df2, m1, m2, rep

    do df = 1, ndf1
        m1 = map1(df)
        do df2 = 1, ndf2
            m2 = map2(df2)
            !$OMP SIMD
            do k = 1, nlayers
                do rep = 1, INNERREPS
                    lhs(m1+k-1) = lhs(m1+k-1) + matrix(k,df2,df,cell) * x(m2+k-1)
                end do
            enddo
        end do
    end do

end subroutine matrix_vector_code_nlayersf_moreops


! Similar to nlayersf but it tries to take advantage that ndf1==8 and there
! is a contiguity btw map1(1:4) and map(4:8).
! Only works when ndf1 is exactly 8
subroutine matrix_vector_code_nlayersf_split(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(nlayers,ndf2,ndf1,ncell_3d/nlayers), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, k, df2, m1_1, m1_2, m2

    do df = 1, ndf1/2
        m1_1 = map1(df)
        m1_2 = map1(df+4)
        do df2 = 1, ndf2
            m2 = map2(df2)
            !$OMP SIMD
            do k = 1, nlayers
                lhs(m1_1+k-1) = lhs(m1_1+k-1) + matrix(k,df2,df,cell) * x(m2+k-1)
                lhs(m1_2+k-1) = lhs(m1_2+k-1) + matrix(k,df2,df+4,cell) * x(m2+k-1)
            end do
        end do
    end do

end subroutine matrix_vector_code_nlayersf_split


! Like the original but manually blocking the outer loop with blocks of size
! 8 (vlen variable) to match the vector length. This this blocks are iterated
! as an inner loop and marked with the pragma SIMD.
subroutine matrix_vector_code_vlen(cell,        &
                              nlayers,     &
                              lhs, x,      &
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)

  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  !Internal variables
  integer, parameter :: vlen = 8
  integer                           :: df, k, ik, i,j,l
  real(kind=r_def), dimension(vlen , ndf2) :: x_e
  real(kind=r_def), dimension(vlen , ndf1) :: lhs_e

  ik = (cell-1)*nlayers
  do k = 0, nlayers-1, vlen
    do df = 1, ndf2
        !$OMP SIMD
        do l = 1,vlen
            x_e(l,df) = x(map2(df)+k+l-1)
        end do
    end do

    lhs_e(:,:) = 0.0_r_def

    do j= 1,ndf1
        do i = 1,ndf2
            !$OMP SIMD
            do l = 1,vlen
                lhs_e(l,j) = lhs_e(l,j) + matrix(j,i,ik+k+l) * x_e(l,i)
            end do
        end do
    end do

    do df = 1, ndf1
        !$OMP SIMD
        do l = 1,vlen
            lhs(map1(df)+k+l-1) = lhs(map1(df)+k+l-1) + lhs_e(l,df)
        end do
    end do
  end do

end subroutine matrix_vector_code_vlen

! Like the original implementation but replacing the matmul intrinsic with
! the specific number of row-column multiplication and addition operations
! needed.
subroutine matrix_vector_code_specialized_save(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, df2, k, ik
  real(kind=r_def), dimension(ndf2, nlayers) :: x_e
  real(kind=r_def), dimension(ndf1, nlayers) :: lhs_e
 
    do k = 0, nlayers-1
        do df = 1, ndf2  
            x_e(df,k+1) = x(map2(df)+k)
        end do
    enddo

    ik = (cell-1)*nlayers

    if (ndf2 == 6) then
        do k = 1, nlayers
            do df = 1, ndf1
                lhs_e(df,k) = matrix(df,1,ik+k) * x_e(1,k) &
                            + matrix(df,2,ik+k) * x_e(2,k) &
                            + matrix(df,3,ik+k) * x_e(3,k) &
                            + matrix(df,4,ik+k) * x_e(4,k) &
                            + matrix(df,5,ik+k) * x_e(5,k) &
                            + matrix(df,6,ik+k) * x_e(6,k)
            end do
        end do
    end if


    do k = 0, nlayers-1
        do df = 1,ndf1
            lhs(map1(df)+k) = lhs(map1(df)+k) + lhs_e(df,k+1) 
        end do
    enddo
 
end subroutine matrix_vector_code_specialized_save


! Like the original implementation but replacing the matmul intrinsic with
! the specific number of row-column multiplication and addition operations
! needed. Also fuses the gather-compute-scatter loops.
subroutine matrix_vector_code_specialized(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, df2, k, ik
 
    ik = (cell-1)*nlayers

    if (ndf2 == 6) then
        do k = 1, nlayers
            do df = 1, ndf1
                lhs(map1(df)+k-1) = lhs(map1(df)+k-1) &
                    + matrix(df,1,ik+k) * x(map2(1)+k-1) &
                    + matrix(df,2,ik+k) * x(map2(2)+k-1) &
                    + matrix(df,3,ik+k) * x(map2(3)+k-1) &
                    + matrix(df,4,ik+k) * x(map2(4)+k-1) &
                    + matrix(df,5,ik+k) * x(map2(5)+k-1) &
                    + matrix(df,6,ik+k) * x(map2(6)+k-1)
            end do
        end do
    end if

end subroutine matrix_vector_code_specialized

! Implementation using BLAS dgemv, to use uncomment library call and build
! with the appropriate compiler flags in the Makefile.
subroutine matrix_vector_code_dgemv(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, k, ik , df2
  real(kind=r_def), dimension(ndf2) :: x_e
  real(kind=r_def), dimension(ndf1) :: lhs_e
 
  do k = 0, nlayers-1
    do df = 1, ndf2  
       x_e(df) = x(map2(df)+k)
    end do

    do df = 1, ndf1
       lhs_e(df) = lhs(map1(df)+k)
    end do

    ik = (cell-1)*nlayers + k + 1

    !call dgemv('N', ndf1, ndf2, 1.0_r_def, matrix(:,:,ik), ndf1, x_e, 1, 1.0_r_def, lhs_e, 1)

    do df = 1,ndf1
       !lhs(map1(df)+k) = lhs(map1(df)+k) + lhs_e(df) 
       lhs(map1(df)+k) = lhs_e(df) 
    end do
 end do
 
end subroutine matrix_vector_code_dgemv

! Implementation using BLAS dgemm, to use uncomment library call and build
! with the appropriate compiler flags in the Makefile.
subroutine matrix_vector_code_dgemm(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, k, ik , df2
  real(kind=r_def), dimension(ndf2,1) :: x_e
  real(kind=r_def), dimension(ndf1,1) :: lhs_e
  real(kind=r_def), parameter :: alpha = 1.0_r_def
  real(kind=r_def), parameter :: beta = 1.0_r_def
 

  do k = 0, nlayers-1
    do df = 1, ndf2  
        x_e(df,1) = x(map2(df)+k)
    end do

    do df = 1, ndf1
        lhs_e(df,1) = lhs(map1(df)+k)
    end do

    ik = (cell-1)*nlayers + k + 1

    ! dgemv (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
    ! call dgemv('N', ndf1, ndf2, 1.0_r_def, matrix(:,:,ik), ndf1, x_e, 1, 1.0_r_def, lhs_e, 1)

    ! dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    !call dgemm('N', 'N', ndf1, 1, ndf2, alpha, matrix(:,:,ik), ndf1, x_e, ndf2, beta, lhs_e, ndf1)


    do df = 1, ndf1
        ! lhs(map1(df)+k) = lhs(map1(df)+k) + lhs_e(df,1)
        lhs(map1(df)+k) = lhs_e(df,1)
    end do
 end do
 
end subroutine matrix_vector_code_dgemm

! Implementation using libxsmm, to use uncomment library call and build
! with the appropriate compiler flags in the Makefile.
subroutine matrix_vector_code_libxsmm(cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 
  !Arguments
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(ndf1,ndf2,ncell_3d), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, k, ik , df2
  real(kind=r_def), dimension(ndf2,1) :: x_e
  real(kind=r_def), dimension(ndf1,1) :: lhs_e
  real(kind=r_def), parameter :: alpha = 1.0_r_def
  real(kind=r_def), parameter :: beta = 1.0_r_def
 

  do k = 0, nlayers-1
    do df = 1, ndf2  
        x_e(df,1) = x(map2(df)+k)
    end do

    do df = 1, ndf1
        lhs_e(df,1) = lhs(map1(df)+k)
    end do

    ik = (cell-1)*nlayers + k + 1

    ! dgemv (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
    ! call dgemv('N', ndf1, ndf2, 1.0_r_def, matrix(:,:,ik), ndf1, x_e, 1, 1.0_r_def, lhs_e, 1)

    ! dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    !call libxsmm_dgemm('N', 'N', ndf1, 1, ndf2, alpha, matrix(:,:,ik), ndf1, x_e, ndf2, beta, lhs_e, ndf1)

    do df = 1, ndf1
        ! lhs(map1(df)+k) = lhs(map1(df)+k) + lhs_e(df,1)
        lhs(map1(df)+k) = lhs_e(df,1)
    end do
 end do
 
end subroutine matrix_vector_code_libxsmm

end module matrix_vector_kernel_mod
