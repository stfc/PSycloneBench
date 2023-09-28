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
        matrix, &
        ! x
        x, map_x, ndf_x, undf_x, &
        ! colour map
        ncolour, ncp_colour, cmap)

    character(len=*), intent(in) :: traverse
    integer, intent(in) :: niters, nlayers, ncell, ncell_3d
    integer, intent(in) :: ndf_lhs, undf_lhs, ndf_x, undf_x
    real(kind=r_def), dimension(undf_lhs),               intent(inout) :: lhs
    real(kind=r_def), dimension(undf_x),                 intent(in)    :: x
    real(kind=r_def), dimension(ndf_lhs,ndf_x,ncell_3d), intent(in)    :: matrix
    integer, intent(in), allocatable, target :: map_lhs(:,:)
    integer, intent(in), allocatable, target :: map_x(:,:)
    integer, intent(in) :: ncolour
    integer, intent(in), dimension(:) :: ncp_colour
    integer, intent(in), dimension(:,:) :: cmap

    integer :: iter, cell, colour, ccell
        
    if (traverse.eq."linear") then  
        write(*,*) "Lineal traversing Version"
        do iter = 1, niters
            ! openmp CPU
            !$omp parallel do default(shared), private(cell)
            do cell = 1, ncell
                call matrix_vector_code_original( &
                        cell, nlayers, &
                        lhs, x, ncell_3d, matrix, &
                        ndf_lhs, undf_lhs, map_lhs(:,cell), &
                        ndf_x, undf_x, map_x(:,cell) )
            enddo
            !$omp end parallel do
        enddo
    elseif (traverse.eq."colouring") then
        write(*,*) "Starting computation with colouring"
        do iter = 1, niters
            do colour = 1, ncolour
                !$omp parallel do, private(ccell, cell)
                do ccell = 1, ncp_colour(colour)
                    cell = cmap(colour, ccell)
                    call matrix_vector_code_original( &
                        cell, nlayers, &
                        lhs, x, ncell_3d, matrix, &
                        ndf_lhs, undf_lhs, map_lhs(:,cell), &
                        ndf_x, undf_x, map_x(:,cell) )
                enddo
                !$omp end parallel do
            enddo
        enddo
    else
        write(*,*) "Not implemented:", traverse
    endif
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

end module psy_layer
