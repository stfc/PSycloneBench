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

    write(*,*) 'psi field, internal region: (',psifld%internal%xstart, ':', &
                                               psifld%internal%xstop,  ',', &
                                               psifld%internal%ystart, ':', &
                                               psifld%internal%ystop,  ')'
    ! Loop over 'columns'
    DO J=1, idim2
      DO I=1, idim1

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
          P(I,J) = PCF*(COS(2.0d0*(I-pfld%internal%xstart)*DI)   & 
               +COS(2.0d0*(J-pfld%internal%ystart)*DJ))+50000.d0
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
    integer  :: i, j, ipsi, jpsi
    real(wp) :: dy

    u => ufld%data
    psi => psifld%data

    ! dy is a property of the mesh
    dy = ufld%grid%dy

    do J=ufld%internal%ystart,ufld%internal%ystop
       do I=ufld%internal%xstart,ufld%internal%xstop

          ! Psi is on F pts and thus starts at (2,2) while U is on U pts
          ! and thus starts at (2,1).
          ! In original code:
          !      DO J=1,N
          !         DO I=1,M
          !            U(I+1,J) = -(PSI(I+1,J+1)-PSI(I+1,J))/DY
          !      or U(2,1) = -(psi(2,2) - psi(2,1))/dy
          !      so U(BL) = F( psi(BL) - psi(B-1L) )
          ! Have to shift i right by one because psi is defined on f points
          ! which have xstart=2. This means it is shifted right relative
          ! to U points in original shallow which already had a halo at
          ! x=1. This ensures initial conditions are identical to those
          ! in original 'shallow.'
          ipsi = i - ufld%internal%xstart + psifld%internal%xstart
          jpsi = j - ufld%internal%ystart + psifld%internal%ystart
          U(I,J) = -(PSI(ipsi,jpsi)-PSI(ipsi,jpsi-1))/dy
       end do
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
    integer  :: I, J
    integer  :: ipsi, jpsi
    real(wp) :: dx

    dx = vfld%grid%dx

    DO J=vfld%internal%ystart, vfld%internal%ystop
       DO I=vfld%internal%xstart, vfld%internal%xstop

          ! We must ensure that initial conditions are identical to those
          ! in original 'shallow.'
          ! Have to shift j up by one because psi is defined on f points
          ! which start at (2,2) in the orig. Shallow. V is on v points
          ! which start at (1,2) in the orig. Shallow.
          !      DO J=1,N
          !         DO I=1,M
          !            V(I,J+1) = (PSI(I+1,J+1)-PSI(I,J+1))/DX
          !     V(1,2) = (psi(2,2) - psi(1,2))/dx
          !  or V(B,L) = F( PSI(B,L) - PSI(B,L-1) )
          ! ipsi == psi%start WHEN i == ufld%start
          ipsi = i - vfld%internal%xstart + psifld%internal%xstart
          jpsi = j - vfld%internal%ystart + psifld%internal%ystart
          vfld%data(I,J) = (psifld%data(ipsi,jpsi)-psifld%data(ipsi-1,jpsi))/dx
       END DO
    END DO

  end subroutine init_velocity_v

END MODULE initial_conditions_mod
