module psy_layer
    implicit none
    integer, parameter :: r_def = 8      !< Default real kind for application.

    ! Fortran to C wrapper interface using iso_c_bindings.
    interface
        subroutine wrapper_c_psy_layer( &
            ! benchmark parameters
            traverse, niters, ncell, nlayers, ncell_3d, &
            ! lhs
            lhs, map_lhs, ndf_lhs, undf_lhs, &
            ! matrix
            matrix, &
            ! x
            x, map_x, ndf_x, undf_x &
        ) bind (C, name="c_psy_layer")
            use iso_c_binding
            character(kind=c_char), intent(in), dimension(*) :: traverse
            integer(kind=c_int), intent(in), value :: niters, ncell, nlayers, ncell_3d
            real(kind=c_double), intent(inout), dimension(*) :: lhs
            real(kind=c_double), intent(in), dimension(*) :: x, matrix
            integer(kind=c_int), intent(in), dimension(*) :: map_lhs, map_x
            integer(kind=c_int), value :: ndf_lhs, undf_lhs, ndf_x, undf_x
        end subroutine wrapper_c_psy_layer
    end interface    

contains

subroutine run_psy_layer( &
        ! benchmark parameters
        traverse, niters, ncell, nlayers, ncell_3d, &
        ! lhs
        lhs, map_lhs, ndf_lhs, undf_lhs, &
        ! matrix
        matrix, &
        ! x
        x, map_x, ndf_x, undf_x)

    character(len=*), intent(in) :: traverse
    integer, intent(in) :: niters, nlayers, ncell, ncell_3d
    integer, intent(in) :: ndf_lhs, undf_lhs, ndf_x, undf_x
    real(kind=r_def), dimension(undf_lhs),               intent(inout) :: lhs
    real(kind=r_def), dimension(undf_x),                 intent(in)    :: x
    real(kind=r_def), dimension(ndf_lhs,ndf_x,ncell_3d), intent(in)    :: matrix
    integer, intent(in), allocatable, target :: map_lhs(:,:)
    integer, intent(in), allocatable, target :: map_x(:,:)

    call wrapper_c_psy_layer( &
        traverse, niters, ncell, nlayers, ncell_3d, &
        lhs, map_lhs, ndf_lhs, undf_lhs, &
        matrix, x, map_x, ndf_x, undf_x)

end subroutine run_psy_layer
end module psy_layer
