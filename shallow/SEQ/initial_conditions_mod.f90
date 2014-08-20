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

    write(*,*) 'psi field, internal region: (',psifld%internal%xstart, ':', &
                                               psifld%internal%xstop,  ',', &
                                               psifld%internal%ystart, ':', &
                                               psifld%internal%ystop,  ')'
    ! Loop over 'columns'
    DO J=1, idim2
      DO I=1, idim1
    !DO J=psifld%internal%ystart, psifld%internal%ystop
    !   DO I=psifld%internal%xstart, psifld%internal%xstop

        CALL init_stream_fn_code(i, j, &
                                 psifld%internal%xstart, & 
                                 psifld%internal%ystart, &
                                 psifld%data)

      END DO
    END DO

  CONTAINS

    SUBROUTINE init_stream_fn_code(i, j, istart, jstart, psi)
      IMPLICIT none
      !> The grid point (column) to act on
      INTEGER,      INTENT(in)                  :: i, j
      INTEGER,      INTENT(in)                  :: istart, jstart
      !> Array holding the stream function values
      REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: psi

      ! di = 2Pi/(Extent of mesh in x)
      ! dj = 2Pi/(Extent of mesh in y)

      PSI(I,J) = A*SIN((I-istart+1.5d0)*DI)*SIN((J-jstart+1.5d0)*DJ)

    END SUBROUTINE init_stream_fn_code

  END SUBROUTINE invoke_init_stream_fn_kernel

  !===================================================

  SUBROUTINE init_pressure(pfld)
    IMPLICIT none
    type(r2d_field_type), target, intent(inout) :: pfld
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

    ! di = 2Pi/(Extent of mesh in x)
    ! dj = 2Pi/(Extent of mesh in y)
!    DO J=pfld%internal%ystart, pfld%internal%ystop
!       DO I=pfld%internal%xstart, pfld%internal%xstop
    DO J=1,idim2
       DO I=1, idim1
          P(I,J) = PCF*(COS(2.0d0*(I-pfld%internal%xstart)*DI)   & 
               +COS(2.0d0*(J-pfld%internal%ystart)*DJ))+50000.d0
       END DO
    END DO

  END SUBROUTINE init_pressure

  !===================================================

  subroutine init_velocity_u(ufld, psifld)
    implicit none
    ! The horizontal velocity field to initialise
    type(r2d_field_type), intent(inout), target :: ufld
    ! The stream function used in the initialisation
    type(r2d_field_type), intent(in),    target :: psifld
    ! Locals
    real(kind=wp), pointer, dimension(:,:) :: u, psi
    integer  :: i, j, ipsi
    real(wp) :: dy

    u => ufld%data
    psi => psifld%data

    ! dy is a property of the mesh
    dy = ufld%grid%dy

    do J=ufld%internal%ystart,ufld%internal%ystop
       do I=ufld%internal%xstart,ufld%internal%xstop

          ! Have to shift i right by one because psi is defined on f points
          ! which have xstart=2. This means it is shifted right relative
          ! to U points in original shallow which already had a halo at
          ! x=1. This ensures initial conditions are identical to those
          ! in original 'shallow.'
          ipsi = i - ufld%internal%xstart + psifld%internal%xstart
          U(I,J) = -(PSI(ipsi,j+1)-PSI(ipsi,j))/dy
       end do
    end do

  end subroutine init_velocity_u

  !===================================================

  SUBROUTINE init_velocity_v(vfld, psifld)
    implicit none
    ! The vertical velocity field to initialise
    type(r2d_field_type), intent(inout), target :: vfld
    ! The stream function used in the initialisation
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
       END DO
    END DO

  CONTAINS

    subroutine init_velocity_v_code(i, j, dx, v, psi)
      implicit none
      integer, intent(in) :: i, j
      real(kind=wp),                 intent(in)  :: dx
      real(kind=wp), dimension(:,:), intent(out) :: v
      real(kind=wp), dimension(:,:), intent(in)  :: psi
      ! Locals
      integer :: jpsi

      ! Have to shift j up by one because psi is defined on f points
      ! which have ystart=2. This means it is shifted upwards relative
      ! to V points in original shallow which already had a halo at
      ! y=1. This ensures initial conditions are identical to those
      ! in original 'shallow.'
      jpsi = j + 1 !- ufld%internal%ystart + psifld%internal%ystart
      V(I,J) = (PSI(I+1,jpsi)-PSI(I,jpsi))/DX

    end subroutine init_velocity_v_code

  end subroutine init_velocity_v

END MODULE initial_conditions_mod
