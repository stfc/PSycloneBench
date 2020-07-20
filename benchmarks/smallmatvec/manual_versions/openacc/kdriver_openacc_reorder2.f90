!-----------------------------------------------------------------------------
! Copyright (c) 2019,  Met Office, on behalf of HMSO and Queen's Printer
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
! Modified by S. Siso and R. W. Ford, STFC Daresbury Lab

! This is a modified version of the driver for the PSyKE-generated
! matrix-vector benchmark.
!
! The PSyKE version reads an input file to create its initial
! conditions and is therefore restricted to the data sizes provided in
! the input file. This version generates the initial conditions in the
! driver and can therefore be used to evaluate different numbers of
! cells in the horizontal and different numbers of levels in the
! vertical.
!
! To facilitate this, this driver expects two command line arguments,
! the first gives the number of cell columns in one of the two
! dimensions and the second gives the number of levels in a
! column. The total number of cell columns is derived by squaring the
! first argument.
!
! For example, "kdriver 30 30" would result in 900 columns of cells
! and 30 levels per column
!
program kdriver

  use constants_mod, only : i_def, r_def, str_max_filename
  use matrix_vector_kernel_mod
  use omp_lib

  implicit none

  ! the arguments and array sizes needed to read/write the arrays.
  integer(kind=i_def) :: ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, ndf_any_space_2_x, undf_any_space_2_x
  integer(kind=i_def) :: ncell, nlayers, ncell_3d, ncolour
  integer(kind=i_def), allocatable :: ncp_colour(:), cmap(:,:)
  integer(kind=i_def), allocatable,target :: map_any_space_1_theta_adv_term(:,:)
  integer(kind=i_def), allocatable,target :: map_any_space_2_x(:,:)
  real(kind=r_def), allocatable :: theta_adv_term_data(:), x_data(:), ans_data(:)
  real(kind=r_def), allocatable :: ptheta_2_local_stencil(:,:,:)
  real(kind=r_def), allocatable :: nlayers_first(:,:,:,:)
  integer(kind=i_def) :: df, df2, lp, nthreads
  real(kind=r_def) :: start, end

  integer(kind=i_def), pointer :: map1(:)=>null(), map2(:)=>null()
  integer(kind=i_def) :: colour, cell, cmap_tmp, map1_tmp, map2_tmp
  integer(kind=i_def) :: c1_index, c2_index, c3_index, c4_index

  integer(kind=i_def) :: nsize, el, i1, i2
  character(len=32) :: arg

  integer :: i, j, k, ik
  ! cmdline argument parsing
  if (command_argument_count() /= 2) then
     write(*,*) "Error parsing <columns> <layers> arguments"
     call exit(-1)
  end if
  call get_command_argument(1, arg)
  read(arg,*) nsize
  call get_command_argument(2, arg)
  read(arg,*) nlayers
  write(*,*) "Generating grid with size ", nsize, "and ", nlayers," layers"

  ! Generate scalars
  ndf_any_space_2_x = 6
  ndf_any_space_1_theta_adv_term = 8
  ncolour = 4

  ncell = nsize * nsize
  undf_any_space_1_theta_adv_term = (nsize+1)*(nsize+1)*(nlayers+1)
  undf_any_space_2_x = (2*(nsize+1)*(nsize)*nlayers) + (nsize*nsize*(nlayers+1))
  ncell_3d = nlayers * ncell

  ! Allocate arrays
  allocate( map_any_space_1_theta_adv_term(ndf_any_space_1_theta_adv_term,ncell))
  allocate( map_any_space_2_x(ndf_any_space_2_x,ncell) )
  allocate( theta_adv_term_data(undf_any_space_1_theta_adv_term) )
  allocate( x_data(undf_any_space_2_x) )
  allocate( ptheta_2_local_stencil( ndf_any_space_1_theta_adv_term, ndf_any_space_2_x, ncell_3d) )
  allocate( nlayers_first(nlayers, ndf_any_space_2_x, ndf_any_space_1_theta_adv_term, ncell) )
  allocate( ans_data(undf_any_space_1_theta_adv_term) )
  allocate( ncp_colour(ncolour) )

  ! Populate arrays
  do i2 = 0, nsize - 1
     do i1 = 0, nsize - 1
        ! Populate  map_any_space_1_theta_adv_term (8 vertices)
        do el = 1, 2
           map_any_space_1_theta_adv_term( 1 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = (i2 * (nsize+1) + i1) * (nlayers + 1) + el
           map_any_space_1_theta_adv_term( 2 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = (i2 * (nsize+1) + i1+1) * (nlayers + 1) + el
           map_any_space_1_theta_adv_term( 3 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = ((i2+1) * (nsize+1) + i1+1) * (nlayers + 1) + el
           map_any_space_1_theta_adv_term( 4 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = ((i2+1) * (nsize+1) + i1) * (nlayers + 1) + el
        enddo

        !Populate  map_any_space_2_x (6 surfaces)
        map_any_space_2_x(1, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + 1
        map_any_space_2_x(2, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + nlayers + 1
        if (i1.eq.nsize-1) then
           map_any_space_2_x(3, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + (i1+1) * (3*nlayers + 1) + 1
        else
           map_any_space_2_x(3, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + (i1+1) * (3*nlayers + 1) + nlayers + 1
        endif
        if (i2.eq.nsize-1) then
           map_any_space_2_x(4, (i2 * nsize) + i1 + 1 ) = (i2+1) * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * nlayers + 1
        else
           map_any_space_2_x(4, (i2 * nsize) + i1 + 1 ) = (i2+1) * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + 1
        endif
        map_any_space_2_x(5, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + 2*nlayers + 1
        map_any_space_2_x(6, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + 2*nlayers + 2
     enddo
  enddo

  ! Number of cells per color = num times per row * num times per column
  ncp_colour(1) = ((nsize/2)+mod(nsize,2)) * ((nsize/2)+mod(nsize,2))
  ncp_colour(2) = ( nsize/2 )              * ((nsize/2)+mod(nsize,2))
  ncp_colour(3) = ((nsize/2)+mod(nsize,2)) * (nsize/2)
  ncp_colour(4) = ( nsize/2 )              * (nsize/2)
  
  allocate(cmap(ncolour,maxval(ncp_colour))) !color map (colour, cell number)
  
  c1_index = 1
  c2_index = 1
  c3_index = 1
  c4_index = 1
  do i2 = 0, nsize - 1
     do i1 = 0, nsize - 1
        if (mod(i1,2).eq.0) then
           if (mod(i2,2).eq.0) then
              cmap(1,c1_index) = (i2 * nsize) + i1 + 1
              c1_index = c1_index + 1
           else
              cmap(2,c2_index) = (i2 * nsize) + i1 + 1
              c2_index = c2_index + 1
           endif
        else
           if (mod(i2,2).eq.0) then
              cmap(3,c3_index) = (i2 * nsize) + i1 + 1
              c3_index = c3_index + 1
           else
              cmap(4,c4_index) = (i2 * nsize) + i1 + 1
              c4_index = c4_index + 1
           endif
        endif
     enddo
  enddo
  
  ! Random data
  call random_number(theta_adv_term_data)
  call random_number(x_data)
  call random_number(ptheta_2_local_stencil)
  theta_adv_term_data = theta_adv_term_data * 100000
  x_data = x_data * 100000
  ptheta_2_local_stencil = ptheta_2_local_stencil * 100000
  ans_data =1

  ! Starting computation phase
  write(*,*) "Starting mv with:"
  write(*,*) "Colors=", ncolour, " -> ", ncp_colour
  write(*,*) "Layers=", nlayers
  write(*,*) "Horitzontal Cells=", ncell_3d/nlayers, ncell
  write(*,*) "3D Cells=", ncell_3d
  write(*,*) "ndf1=", ndf_any_space_1_theta_adv_term, ", undf1=", undf_any_space_1_theta_adv_term
  write(*,*) "ndf2=", ndf_any_space_2_x, ", undf2=", undf_any_space_2_x

  start = omp_get_wtime()

  do lp = 1, 10000
     do colour = 1, ncolour
          !$acc parallel loop independent
          do cell = 1, ncp_colour(colour)

             ! MANUAL INLINE
             cmap_tmp = cmap(colour,cell)
             ik = (cmap_tmp-1)*nlayers

             do i = 1, ndf_any_space_1_theta_adv_term
                map1_tmp = map_any_space_1_theta_adv_term(i,cmap_tmp)
                do j = 1, ndf_any_space_2_x
                   !map2_tmp = map_any_space_2_x(j,cmap_tmp)
                   !$acc loop vector independent
                   do k = 1, nlayers
                      !theta_adv_term_data(map1_tmp+k-1) = theta_adv_term_data(map1_tmp+k-1) + ptheta_2_local_stencil(ik+k,i,j) * x_data(map2_tmp+k-1)
                      theta_adv_term_data(map1_tmp+k-1) = theta_adv_term_data(map1_tmp+k-1) + ptheta_2_local_stencil(ik+k,i,j) * x_data(map_any_space_2_x(j,cmap_tmp)+k-1)
                      !theta_adv_term_data(map_any_space_1_theta_adv_term(i,cmap_tmp)+k-1) = theta_adv_term_data(map_any_space_1_theta_adv_term(i,cmap_tmp)+k-1) + ptheta_2_local_stencil(ik+k,i,j) * x_data(map_any_space_2_x(j,cmap_tmp)+k-1)
                   end do
                end do
             end do
          end do
          !$acc end parallel loop
     end do
  end do

  end = omp_get_wtime()
  write(*,*) "T=",nthreads, "Loop time=",end - start, " s"

  ! print out reduction sum
  write(*,*) "Reduction value:", sum(theta_adv_term_data)

  deallocate( map_any_space_1_theta_adv_term, map_any_space_2_x )
  deallocate( theta_adv_term_data, x_data, ptheta_2_local_stencil )

end program kdriver
