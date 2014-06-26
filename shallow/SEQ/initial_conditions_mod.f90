module initial_conditions_mod
  use physical_params_mod
  use field_mod
  implicit none
  private

  !> Amplitude of initial oscillations in stream function
  !! Used by invoke_init_stream_fn_kernel()
  REAL(wp), PARAMETER :: A = 1.0E6
  !> 2PI/{m,n}
  REAL(wp)  :: di, dj

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
    type(r2d_field_type), intent(in) :: pfld

    di = TPI/pfld%internal%nx
    dj = TPI/pfld%internal%ny


  END SUBROUTINE init_initial_condition_params

  !===================================================

  subroutine invoke_init_stream_fn_kernel(psifld)
    implicit none
    type(r2d_field_type), intent(inout) :: psifld
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

      PSI(I,J) = A*SIN((I-.5)*DI)*SIN((J-.5)*DJ)

    END SUBROUTINE init_stream_fn_code

  END SUBROUTINE invoke_init_stream_fn_kernel

  !===================================================

  SUBROUTINE init_pressure(pfld)
    IMPLICIT none
    type(r2d_field_type), target, intent(inout) :: pfld
    REAL(KIND=wp), DIMENSION(:,:), pointer :: p
    ! Locals
    INTEGER :: idim1, idim2
    INTEGER :: i, j
    !> Extent in x of model domain
    REAL(wp) :: el
    !> Computed amplitude of initial osciallations in
    !! pressure field.
    REAL(wp) :: pcf

    p => pfld%data

    idim1 = pfld%grid%nx
    idim2 = pfld%grid%ny

    EL = pfld%internal%nx * pfld%grid%dx
    PCF = PI*PI*A*A/(EL*EL)

    ! di = 2Pi/(Extent of mesh in x)
    ! dj = 2Pi/(Extent of mesh in y)
    ! Fields have one extra row and column than the mesh
    ! extent. Presumably because e.g. velocity fields are
    ! offset from pressure but this not dealt with that
    ! explicitly. Only field(1:m,1:n) is sent to be
    ! written to the netcdf file.
    DO J=1,idim2
       DO I=1,idim1
          P(I,J) = PCF*(COS(2.*(I-1)*DI)   & 
               +COS(2.*(J-1)*DJ))+50000.
       END DO
    END DO

  END SUBROUTINE init_pressure

  !===================================================

  subroutine init_velocity_u(ufld, psifld, m, n)
    implicit none
    type(r2d_field_type), intent(inout), target :: ufld
    type(r2d_field_type), intent(in),    target :: psifld
    integer,      intent(in) :: m, n
    ! Locals
    real(kind=wp), pointer, dimension(:,:) :: u, psi
    integer  :: i, j
    real(wp) :: dy

    u => ufld%data
    psi => psifld%data

    ! dy is a property of the mesh
    dy = ufld%grid%dy

    do J=1,N
       do I=1,M+1
          U(I,J) = -(PSI(I,J+1)-PSI(I,J))/dy
       end do
    end do

  end subroutine init_velocity_u

  !===================================================

  SUBROUTINE init_velocity_v(vfld, psifld)
    implicit none
    type(r2d_field_type), intent(inout), target :: vfld
    type(r2d_field_type), intent(in),    target :: psifld
    ! Locals
    real(kind=wp), pointer, dimension(:,:) :: v, psi
    integer  :: I, J
    real(wp) :: dx

    v => vfld%data
    psi => psifld%data
    dx = vfld%grid%dx

    DO J=vfld%internal%ystart, vfld%internal%ystop
       DO I=vfld%internal%xstart, vfld%internal%xstop
          call init_velocity_v_code(i,j,dx,v,psi)
!          V(I,J) = (PSI(I+1,J)-PSI(I,J))/DX
       END DO
    END DO

  CONTAINS

    subroutine init_velocity_v_code(i, j, dx, v, psi)
      implicit none
      integer, intent(in) :: i, j
      real(kind=wp),                 intent(in)  :: dx
      real(kind=wp), dimension(:,:), intent(out) :: v
      real(kind=wp), dimension(:,:), intent(in)  :: psi

      V(I,J) = (PSI(I+1,J)-PSI(I,J))/DX

    end subroutine init_velocity_v_code

  end subroutine init_velocity_v

END MODULE initial_conditions_mod
