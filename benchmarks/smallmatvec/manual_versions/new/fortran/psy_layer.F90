module psy_layer
    implicit none
    integer, parameter :: r_def = 8      !< Default real kind for application.
contains

subroutine run_psy_layer( &
        ! benchmark parameters
        traverse, niters, ncell, nlayers, ncell_3d, &
        ! lhs
        lhs, map_lhs, ndf_lhs, undf_lhs, &
        ! matrix
        matrix, matrix_kinner, &
        ! x
        x, map_x, ndf_x, undf_x, &
        ! colour map
        ncolour, ncp_colour, cmap)

    character(len=*), intent(in) :: traverse
    integer, intent(in) :: niters, nlayers, ncell, ncell_3d
    integer, intent(in) :: ndf_lhs, undf_lhs, ndf_x, undf_x
    real(kind=r_def), dimension(undf_lhs),               intent(inout)   :: lhs
    real(kind=r_def), dimension(undf_x),                 intent(in)      :: x
    real(kind=r_def), dimension(ndf_lhs,ndf_x,ncell_3d), intent(in)      :: matrix
    real(kind=r_def), dimension(nlayers,ndf_lhs,ndf_x,ncell), intent(in) :: matrix_kinner
    integer, intent(in), allocatable, target :: map_lhs(:,:)
    integer, intent(in), allocatable, target :: map_x(:,:)
    integer, intent(in) :: ncolour
    integer, intent(in), dimension(:) :: ncp_colour
    integer, intent(in), dimension(:,:) :: cmap

    integer :: iter, cell, colour, ccell
#ifdef TARGET_GPU
   !$omp target enter data map(to: ncp_colour, cmap, x, matrix, map_x, lhs, map_lhs, matrix_kinner)
#endif 
    if (traverse.eq."linear") then  
        write(*,*) "Lineal traversing Version"
        do iter = 1, niters
#ifdef TARGET_GPU
            !$omp target loop
#else 
            !$omp parallel do default(shared), private(cell)
#endif
            do cell = 1, ncell
                call matrix_vector_code_optimised_atomic( &
                        cell, nlayers, &
                        lhs, x, ncell_3d, matrix, &
                        ndf_lhs, undf_lhs, map_lhs(:,cell), &
                        ndf_x, undf_x, map_x(:,cell) )
            enddo
#ifdef TARGET_GPU
            !$omp end target loop
#else
            !$omp end parallel do
#endif
        enddo
    elseif (traverse.eq."colouring") then
        write(*,*) "Starting computation with colouring"
        do iter = 1, niters
            do colour = 1, ncolour
#ifdef TARGET_GPU
            !$omp target loop
#else 
            !$omp parallel do default(shared), private(ccell, cell)
#endif
                do ccell = 1, ncp_colour(colour)
                    cell = cmap(colour, ccell)
                    call matrix_vector_code_optimised( &
                        cell, nlayers, &
                        lhs, x, ncell_3d, matrix, &
                        ndf_lhs, undf_lhs, map_lhs(:,cell), &
                        ndf_x, undf_x, map_x(:,cell) )
                enddo
#ifdef TARGET_GPU
           !$omp end target loop
#else
           !$omp end parallel do
#endif
            enddo
        enddo

   elseif (traverse.eq."linear-kinner") then
        write(*,*) "Starting computation with linear and kinner"
        do iter = 1, niters
#ifdef TARGET_GPU
        !$omp target loop
#else
        !$omp parallel do default(shared), private(cell)
#endif
            do cell = 1, ncell
                call matrix_vector_code_kinner_atomics( &
                        cell, nlayers, &
                        lhs, x, ncell_3d, matrix_kinner, &
                        ndf_lhs, undf_lhs, map_lhs(:,cell), &
                        ndf_x, undf_x, map_x(:,cell) )
            enddo
#ifdef TARGET_GPU
           !$omp end target loop
#else
           !$omp end parallel do
#endif
        enddo
 
    elseif (traverse.eq."colouring-kinner") then
        write(*,*) "Starting computation with colouring and kinner"
        do iter = 1, niters
            do colour = 1, ncolour

#ifdef TARGET_GPU
        !$omp target loop
#else
        !$omp parallel do default(shared), private(cell, ccell)
#endif
                do ccell = 1, ncp_colour(colour)
                    cell = cmap(colour, ccell)
                    call matrix_vector_code_kinner_atomics( &
                        cell, nlayers, &
                        lhs, x, ncell_3d, matrix_kinner, &
                        ndf_lhs, undf_lhs, map_lhs(:,cell), &
                        ndf_x, undf_x, map_x(:,cell) )
                enddo
#ifdef TARGET_GPU
           !$omp end target loop
#else
           !$omp end parallel do
#endif
            enddo
        enddo

    else
        write(*,*) "Not implemented:", traverse
    endif
#ifdef TARGET_GPU
        !$omp target update from (lhs)
#endif

end subroutine run_psy_layer

!> @brief Computes lhs = matrix*x
!! @param[in] cell Horizontal cell index
!! @param[in] nlayers Number of layers
!! @param[inout] lhs Output lhs (A*x)
!! @param[in] x Input data
!! @param[in] ncell_3d Total number of cells
!! @param[in] matrix Local matrix assembly form of the operator A
!! @param[in] ndf1 Number of degrees of freedom per cell for the output field
!! @param[in] undf1 Unique number of degrees of freedom  for the output field
!! @param[in] map1 Dofmap for the cell at the base of the column for the
!! @param[in] ndf2 Number of degrees of freedom per cell for the input field
!! @param[in] undf2 Unique number of degrees of freedom for the input field
!! @param[in] map2 Dofmap for the cell at the base of the column for the input
subroutine matrix_vector_code_original( &
                              cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 

! !$omp declare target
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
  integer                           :: df, k, ik !, df2
  real(kind=r_def), dimension(ndf2) :: x_e
  real(kind=r_def), dimension(ndf1) :: lhs_e
  
  do k = 0, nlayers-1
    do df = 1, ndf2  
       x_e(df) = x(map2(df)+k)
    end do

    ik = (cell-1)*nlayers + k + 1

    lhs_e = matmul(matrix(:,:,ik),x_e)
    
    do df = 1,ndf1
       lhs(map1(df)+k) = lhs(map1(df)+k) + lhs_e(df) 
    end do
  end do

  ! Optimised implementation (remove temporals by interleaving, bring k inside)
  ! ik = (cell-1)*nlayers
  ! do df = 1, ndf1
  !     do df2 = 1, ndf2
  !         !$OMP SIMD
  !         do k = 1, nlayers
  !             lhs(map1(df)+k-1) = lhs(map1(df)+k-1) + matrix(df,df2,ik+k) * x(map2(df2)+k-1)
  !         end do
  !     end do
  ! end do
end subroutine matrix_vector_code_original

subroutine matrix_vector_code_optimised( &
                              cell,        &
                              nlayers,     &
                              lhs, x,      & 
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
 

! !$omp declare target
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
  integer                           :: df, k, ik, df2
  real(kind=r_def), dimension(ndf2) :: x_e
  real(kind=r_def), dimension(ndf1) :: lhs_e

   !Optimised implementation (remove temporals by interleaving, bring k inside)
   ik = (cell-1)*nlayers
   do df = 1, ndf1
       do df2 = 1, ndf2
           !$OMP SIMD
           do k = 1, nlayers
               lhs(map1(df)+k-1) = lhs(map1(df)+k-1) + matrix(df,df2,ik+k) * x(map2(df2)+k-1)
           end do
       end do
   end do
end subroutine matrix_vector_code_optimised

subroutine matrix_vector_code_optimised_atomic( &
                              cell,        &
                              nlayers,     &
                              lhs, x,      &
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)


! !$omp declare target
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
  integer                           :: df, k, ik, df2
  real(kind=r_def), dimension(ndf2) :: x_e
  real(kind=r_def), dimension(ndf1) :: lhs_e

   !Optimised implementation (remove temporals by interleaving, bring k inside)
   ik = (cell-1)*nlayers
   do df = 1, ndf1
       do df2 = 1, ndf2
    !       !$OMP SIMD
           do k = 1, nlayers
               !$OMP ATOMIC
               lhs(map1(df)+k-1) = lhs(map1(df)+k-1) + matrix(df,df2,ik+k) * x(map2(df2)+k-1)
           end do
       end do
   end do
end subroutine matrix_vector_code_optimised_atomic

subroutine matrix_vector_code_original_atomic( &
                              cell,        &
                              nlayers,     &
                              lhs, x,      &
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)

  !!$omp declare target
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
  integer                           :: df, k, ik
  real(kind=r_def), dimension(ndf2) :: x_e
  real(kind=r_def), dimension(ndf1) :: lhs_e
  
  do k = 0, nlayers-1
    do df = 1, ndf2
       x_e(df) = x(map2(df)+k)
    end do

    ik = (cell-1)*nlayers + k + 1

    lhs_e = matmul(matrix(:,:,ik),x_e)

    do df = 1,ndf1
       !$omp atomic
       lhs(map1(df)+k) = lhs(map1(df)+k) + lhs_e(df)
    end do
 end do
end subroutine matrix_vector_code_original_atomic

subroutine matrix_vector_code_kinner_atomics(cell,        &
                              nlayers,     &
                              lhs, x,      &
                              ncell_3d,    &
                              matrix,      &
                              ndf1, undf1, map1, &
                              ndf2, undf2, map2)
  integer,                   intent(in)    :: cell, nlayers, ncell_3d
  integer,                   intent(in)    :: undf1, ndf1
  integer,                   intent(in)    :: undf2, ndf2
  integer, dimension(ndf1),  intent(in)    :: map1
  integer, dimension(ndf2),  intent(in)    :: map2
  real(kind=r_def), dimension(undf2),              intent(in)    :: x
  real(kind=r_def), dimension(undf1),              intent(inout) :: lhs
  real(kind=r_def), dimension(nlayers,ndf1,ndf2,ncell_3d/nlayers), intent(in)    :: matrix

  !Internal variables
  integer                           :: df, k, df2, m1,m2

  do df2 = 1, ndf2
      do df = 1, ndf1
          !$OMP LOOP BIND(PARALLEL) 
          do k = 1, nlayers
              !$OMP ATOMIC
              lhs(map1(df)+k-1) = lhs(map1(df)+k-1) + matrix(k,df,df2,cell) * x(map2(df2)+k-1)
          end do
      end do
  end do

end subroutine matrix_vector_code_kinner_atomics

end module psy_layer
