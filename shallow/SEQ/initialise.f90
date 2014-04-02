MODULE initialise
  IMPLICIT none

CONTAINS

  !===================================================

  SUBROUTINE init_stream_fn(psi, A, DI, DJ)
    IMPLICIT none
    REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: psi
    REAL(KIND=8), INTENT(in) :: A, DI, DJ
    ! Locals
    INTEGER :: idim1, idim2
    INTEGER :: I,J

    idim1 = SIZE(psi, 1)
    idim2 = SIZE(psi, 2)

    ! di = 2Pi/(Extent of mesh in x)
    ! dj = 2Pi/(Extent of mesh in y)

    DO J=1, idim2
       DO I=1, idim1
          PSI(I,J) = A*SIN((I-.5)*DI)*SIN((J-.5)*DJ)
       END DO
    END DO

  END SUBROUTINE init_stream_fn

  !===================================================

  SUBROUTINE init_pressure(p, pcf, di, dj)
    IMPLICIT none
    REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: p
    REAL(KIND=8), INTENT(in) :: pcf, di, dj
    INTEGER :: idim1, idim2
    INTEGER :: i, j

    idim1 = SIZE(p, 1)
    idim2 = SIZE(p, 2)
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

  SUBROUTINE init_velocity_u(u, psi, m, n, dy)
    IMPLICIT none
    REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: u
    REAL(KIND=8), INTENT(in),  DIMENSION(:,:) :: psi
    INTEGER,      INTENT(in) :: m, n
    REAL(KIND=8), INTENT(in) :: dy
    ! Locals
    INTEGER :: I, J

    ! dy is a property of the mesh
    DO J=1,N
       DO I=1,M+1
          U(I,J) = -(PSI(I,J+1)-PSI(I,J))/DY
       END DO
    END DO
  END SUBROUTINE init_velocity_u

  !===================================================

  SUBROUTINE init_velocity_v(v, psi, m, n, dx)
    IMPLICIT none
    REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: v
    REAL(KIND=8), INTENT(in),  DIMENSION(:,:) :: psi
    INTEGER, INTENT(in) :: m, n
    REAL(KIND=8), INTENT(in) :: dx
    ! Locals
    INTEGER :: I, J

    ! dx is a property of the mesh
    DO J=1,N+1
       DO I=1,M
          V(I,J) = (PSI(I+1,J)-PSI(I,J))/DX
       END DO
    END DO
  END SUBROUTINE init_velocity_v

END MODULE initialise
