module initial_conditions_mod
  use physical_params_mod
  use field_mod
  implicit none
  private

  !> Amplitude of initial oscillations in stream function
  !! Used by invoke_init_stream_fn_kernel()
  REAL(wp), PARAMETER :: A = 1.0D6
  !> 2PI/{m,n}
  REAL(wp), SAVE  :: di, dj

  PUBLIC init_initial_condition_params
  PUBLIC invoke_init_stream_fn_kernel
  PUBLIC init_pressure
  PUBLIC init_velocity_u
  PUBLIC init_velocity_v

CONTAINS

  !===================================================

  !> \brief Set-up parameters related to the model domain which
  !! are stored in this module. We could compute these on the
  !! fly in init_stream_fn_code() and init_pressure() and
  !! rely on compiler magic to make sure they're not
  !! recomputed for every grid point.
  SUBROUTINE init_initial_condition_params(pfld)
    IMPLICIT none
    type(r2d_field), intent(in) :: pfld

    di = TPI/pfld%internal%nx
    dj = TPI/pfld%internal%ny

  END SUBROUTINE init_initial_condition_params

  !===================================================

  subroutine invoke_init_stream_fn_kernel(psifld)
    implicit none
    type(r2d_field), intent(inout) :: psifld
    ! Locals
    integer :: idim1, idim2
    integer :: i, j

    idim1 = SIZE(psifld%data, 1)
    idim2 = SIZE(psifld%data, 2)

    ! Loop over 'columns'
    DO J=1, idim2
      DO I=1, idim1

        CALL init_stream_fn_code(i, j, psifld%data)

      END DO
    END DO

  CONTAINS

    SUBROUTINE init_stream_fn_code(i, j, psi)
      IMPLICIT none
      !> The grid point (column) to act on
      INTEGER,      INTENT(in)                  :: i, j
      !> Array holding the stream function values
      REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: psi

      ! di = 2Pi/(Extent of mesh in x)
      ! dj = 2Pi/(Extent of mesh in y)
      ! Original code:
      !      PSI(I,J) = A*SIN((I-.5d0)*DI)*SIN((J-.5d0)*DJ)
      ! Psi is on F points. Our first internal F point is at 3,3
      ! rather than 2,2 so shift appropriately.
      PSI(I,J) = A*SIN((I-1.5d0)*DI)*SIN((J-1.5d0)*DJ)

    END SUBROUTINE init_stream_fn_code

  END SUBROUTINE invoke_init_stream_fn_kernel

  !===================================================

  SUBROUTINE init_pressure(pfld)
    IMPLICIT none
    type(r2d_field), target, intent(inout) :: pfld
    REAL(KIND=wp), DIMENSION(:,:), pointer :: p
    ! Locals
    INTEGER :: i, j, idim1, idim2
    !> Extent in x of model domain
    REAL(wp) :: el
    !> Computed amplitude of initial oscillations in
    !! pressure field.
    REAL(wp) :: pcf

    p => pfld%data

    EL = pfld%internal%nx * pfld%grid%dx
    PCF = PI*PI*A*A/(EL*EL)
     
    idim1 = SIZE(pfld%data, 1)
    idim2 = SIZE(pfld%data, 2)

    ! di = 2Pi/(Extent of mesh in x) where extent is from namelist
    ! dj = 2Pi/(Extent of mesh in y)   "     "     "   "     "
    DO J=1,idim2
       DO I=1, idim1
! Original code:
!          P(I,J) = PCF*(COS(2.0d0*(I-1)*DI)   & 
!               +COS(2.0d0*(J-1)*DJ))+50000.d0
! Our first internal T pt is at 2,2 rather than 1,1 so we shift
! the expression appropriately:
          P(I,J) = PCF*(COS(2.0d0*(I-2)*DI)   & 
               +COS(2.0d0*(J-2)*DJ))+50000.d0
       END DO
    END DO

  END SUBROUTINE init_pressure

  !===================================================

  subroutine init_velocity_u(ufld, psifld)
    implicit none
    ! The horizontal velocity field to initialise
    type(r2d_field), intent(inout), target :: ufld
    ! The stream function used in the initialisation
    type(r2d_field), intent(in),    target :: psifld
    ! Locals
    real(kind=wp), pointer, dimension(:,:) :: u, psi
    integer  :: j
    real(wp) :: dy

    u => ufld%data
    psi => psifld%data

    ! dy is a property of the mesh
    dy = ufld%grid%dy

    do J=1,ufld%internal%ystop
       U(:,J) = -(PSI(:,j+1) - PSI(:,j))/dy
    end do

  end subroutine init_velocity_u

  !===================================================

  SUBROUTINE init_velocity_v(vfld, psifld)
    implicit none
    ! The vertical velocity field to initialise
    type(r2d_field), intent(inout), target :: vfld
    ! The stream function used in the initialisation
    type(r2d_field), intent(in),    target :: psifld
    ! Locals
    real(kind=wp), pointer, dimension(:,:) :: v, psi
    integer  :: I
    real(wp) :: dx

    v => vfld%data
    psi => psifld%data

    dx = vfld%grid%dx

    do I=1, vfld%internal%xstop
       v(I,:) = (psi(i+1,:) - psi(i,:))/dx
    end do

  end subroutine init_velocity_v

end module initial_conditions_mod
