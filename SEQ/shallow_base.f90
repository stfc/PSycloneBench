PROGRAM shallow

!     BENCHMARK WEATHER PREDICTION PROGRAM FOR COMPARING THE
!     PREFORMANCE OF CURRENT SUPERCOMPUTERS. THE MODEL IS
!     BASED OF THE PAPER - THE DYNAMICS OF FINITE-DIFFERENCE
!     MODELS OF THE SHALLOW-WATER EQUATIONS, BY ROBERT SADOURNY
!     J. ATM. SCIENCES, VOL 32, NO 4, APRIL 1975.
!     
!     CODE BY PAUL N. SWARZTRAUBER, NATIONAL CENTER FOR
!     ATMOSPHERIC RESEARCH, BOULDER, CO,  OCTOBER 1984.
!     Modified by Juliana Rew, NCAR, January 2006
!
!     In this version, shallow4.f, initial and calculated values
!     of U, V, and P are written to a netCDF file
!     for later use in visualizing the results. The netCDF data
!     management library is freely available from
!     http://www.unidata.ucar.edu/software/netcdf
!     This code is still serial but has been brought up to modern
!     Fortran constructs and uses portable intrinsic Fortran 90 timing routines
!     This can be compiled on the IBM SP using:
!     xlf90 -qmaxmem=-1 -g -o shallow4 -qfixed=132 -qsclk=micro \
!     -I/usr/local/include shallow4.f -L/usr/local/lib32/r4i4 -l netcdf
!     where the -L and -I point to local installation of netCDF
!     
!     Changes from shallow4.f (Annette Osprey, January 2010):
!     - Converted to free-form fortran 90.  
!     - Some tidying up of old commented-out code.   
!     - Explicit type declarations.
!     - Variables n, m, itmax and mprint read in from namelist. 
!     - Dynamic array allocation.
!     - Only write to netcdf at mprint timesteps.
!     - Don't write wrap-around points to NetCDF file.
!     - Use 8-byte reals. 
!
!     This version heavily modified as part of the GOcean-2D project
!     with the mantra "all computation must occur in a kernel."
!     Andrew Porter, April 2014

  USE shallow_IO
  USE manual_invoke_initialise
  USE timing
  USE model
  IMPLICIT NONE

  ! solution arrays
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) ::                & 
                             u, v, p, unew, vnew, pnew,       & 
                             uold, vold, pold, cu, cv, z, h, psi  

  REAL(KIND=8) :: dt, tdt, dx, dy, alpha, & 
                  tdts8, tdtsdx, tdtsdy, fsdx, fsdy

  !> Checksum used for each array
  REAL(KIND=8) :: csum

  !> Loop counter for time-stepping loop
  INTEGER :: ncycle
   
  !> Integer tags for timers
  INTEGER :: idxt0, idxt1

  !  ** Initialisations ** 
  CALL model_init()

  ! NOTE BELOW THAT TWO DELTA T (TDT) IS SET TO DT ON THE FIRST
  ! CYCLE AFTER WHICH IT IS RESET TO DT+DT.
  DT = 90.
  TDT = DT
 
  DX = 1.E5
  DY = 1.E5

  ! Parameter for time smoothing
  ALPHA = .001

  FSDX = 4./DX
  FSDY = 4./DY

  MP1 = M+1
  NP1 = N+1

  CALL invoke_init_model_params_kernel(DX, M, N)

  !     Set up arrays

  ALLOCATE( u(MP1,NP1), v(MP1,NP1), p(MP1,NP1) ) 
  ALLOCATE( unew(MP1,NP1), vnew(MP1,NP1), pnew(MP1,NP1) ) 
  ALLOCATE( uold(MP1,NP1), vold(MP1,NP1), pold(MP1,NP1) )
  ALLOCATE( cu(MP1,NP1), cv(MP1,NP1) ) 
  ALLOCATE( z(MP1,NP1), h(MP1,NP1), psi(MP1,NP1) ) 

  !     INITIAL VALUES OF THE STREAM FUNCTION AND P

  CALL invoke_init_stream_fn_kernel(PSI)
  CALL init_pressure(P)

  !     INITIALIZE VELOCITIES

  CALL init_velocity_u(u, psi, m, n, dy)
  CALL init_velocity_v(v, psi, m, n, dx)

  !     PERIODIC CONTINUATION
  CALL apply_bcs_u(U)
  CALL apply_bcs_v(V)

  ! Initialise fields that will hold data at previous time step
  CALL copy_field(U, UOLD)
  CALL copy_field(V, VOLD)
  CALL copy_field(P, POLD)
     
  !     PRINT INITIAL VALUES
  CALL print_initial_values(n,m,dx,dy,dt,alpha, p, u, v)

  ! Write intial values of p, u, and v into a netCDF file   
  CALL model_write(0, p, u, v)

  !     Start timer
  CALL timer_start('Time-stepping',idxt0)

  !  ** Start of time loop ** 
  DO ncycle=1,itmax
    
    ! COMPUTE CAPITAL U, CAPITAL V, Z AND H

    CALL timer_start('Compute c{u,v},z,h', idxt1)

    CALL compute_cu(CU, P, U)
    CALL compute_cv(CV, P, V)
    CALL compute_z(z, P, U, V, FSDX, FSDY)
    CALL compute_h(h, P, U, V)

    CALL timer_stop(idxt1)

    ! PERIODIC CONTINUATION

    CALL apply_bcs_u(CU)
    CALL apply_bcs_p(H)
    CALL apply_bcs_v(CV)
    CALL apply_bcs_z(Z)

    ! COMPUTE NEW VALUES U,V AND P
    TDTS8 = TDT/8.
    TDTSDX = TDT/DX
    TDTSDY = TDT/DY

    CALL timer_start('Compute new fields', idxt1)
     
    CALL compute_unew(unew, uold, z, cv, h, TDTS8, TDTSDX)
    CALL compute_vnew(vnew, vold, z, cu, h, TDTS8, TDTSDY)
    CALL compute_pnew(pnew, pold, cu, cv, TDTSDX, TDTSDY)

    CALL timer_stop(idxt1)

    ! PERIODIC CONTINUATION

    CALL apply_bcs_u(UNEW)
    CALL apply_bcs_v(VNEW)
    CALL apply_bcs_p(PNEW)

    ! Time is in seconds but we never actually need it
    !TIME = TIME + DT

    CALL model_write(ncycle, p, u, v)

    ! TIME SMOOTHING AND UPDATE FOR NEXT CYCLE
    IF(NCYCLE .GT. 1) then

      CALL timer_start('Time smoothing',idxt1)

      CALL time_smooth(U, UNEW, UOLD, ALPHA)
      CALL time_smooth(V, VNEW, VOLD, ALPHA)
      CALL time_smooth(P, PNEW, POLD, ALPHA)

      CALL timer_stop(idxt1)

      CALL copy_field(UNEW, U)
      CALL copy_field(VNEW, V)
      CALL copy_field(PNEW, P)

    ELSE ! ncycle == 1

      ! Make TDT actually = 2*DT
      TDT = TDT+TDT

      CALL copy_field(U, UOLD)
      CALL copy_field(V, VOLD)
      CALL copy_field(P, POLD)
         
      CALL copy_field(UNEW, U)
      CALL copy_field(VNEW, V)
      CALL copy_field(PNEW, P)

    ENDIF ! ncycle > 1

  END DO

  !  ** End of time loop ** 

  CALL timer_stop(idxt0)

  CALL compute_checksum(pnew, csum)
  CALL model_write_log("('P CHECKSUM after ',I6,' steps = ',E15.7)", &
                       itmax, csum)

  CALL compute_checksum(unew, csum)
  CALL model_write_log("('U CHECKSUM after ',I6,' steps = ',E15.7)", &
                       itmax, csum)

  CALL compute_checksum(vnew, csum)
  CALL model_write_log("('V CHECKSUM after ',I6,' steps = ',E15.7)", &
                       itmax, csum)

  CALL model_finalise()

  !> Free memory \todo Move to model_finalise()
  DEALLOCATE( u, v, p, unew, vnew, pnew, uold, vold, pold )
  DEALLOCATE( cu, cv, z, h, psi ) 

CONTAINS

  !===================================================

  SUBROUTINE compute_checksum(field, val)
    IMPLICIT none
    REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: field
    REAL(KIND=8), INTENT(out) :: val

    val = SUM(field)

  END SUBROUTINE compute_checksum

  !===================================================

  SUBROUTINE apply_bcs_p(field)
    IMPLICIT none
    REAL(KIND=8), INTENT(inout), DIMENSION(:,:) :: field
    INTEGER :: M, MP1, N, NP1

    MP1 = SIZE(field, 1)
    NP1 = SIZE(field, 2)
    M = MP1 - 1
    N = NP1 - 1

    ! Last col = first col
    field(MP1,1:N) = field(1,  1:N)
    ! Last row = first row
    field(1:MP1,NP1) = field(1:MP1,1)

  END SUBROUTINE apply_bcs_p

  !===================================================

      SUBROUTINE apply_bcs_u(field)
        IMPLICIT none
        REAL(KIND=8), INTENT(inout), DIMENSION(:,:) :: field
        INTEGER :: M, MP1, N, NP1

        MP1 = SIZE(field, 1)
        NP1 = SIZE(field, 2)
        M = MP1 - 1
        N = NP1 - 1

        ! First col = last col
        field(1,    1:N) = field(MP1,  1:N)
        ! Last row = first row
        field(1:MP1,NP1) = field(1:MP1,1)

      END SUBROUTINE apply_bcs_u

      !===================================================

      SUBROUTINE apply_bcs_v(field)
        IMPLICIT none
        REAL(KIND=8), INTENT(inout), DIMENSION(:,:) :: field
        INTEGER :: M, MP1, N, NP1

        MP1 = SIZE(field, 1)
        NP1 = SIZE(field, 2)
        M = MP1 - 1
        N = NP1 - 1

        ! First row = last row
        field(1:M,1    ) = field(1:M,NP1)
        ! Last col = first col
        field(MP1,1:NP1) = field(1,  1:NP1)

      END SUBROUTINE apply_bcs_v

      !===================================================

      SUBROUTINE apply_bcs_z(field)
        IMPLICIT none
        REAL(KIND=8), INTENT(inout), DIMENSION(:,:) :: field
        INTEGER :: M, MP1, N, NP1

        MP1 = SIZE(field, 1)
        NP1 = SIZE(field, 2)
        M = MP1 - 1
        N = NP1 - 1

        ! First col = last col
        field(1,    2:NP1) = field(MP1,  2:NP1)
        ! First row = last row
        field(1:MP1,1)     = field(1:MP1,NP1)

      END SUBROUTINE apply_bcs_z
 
      !===================================================

      SUBROUTINE copy_field(field_in, field_out)
        IMPLICIT none
        REAL(KIND=8), INTENT(in),  DIMENSION(:,:) :: field_in
        REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: field_out
        
        field_out(:,:) = field_in(:,:)
        
      END SUBROUTINE copy_field

      !===================================================

      SUBROUTINE time_smooth(field, field_new, field_old, ALPHA)
        IMPLICIT none
        REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: field
        REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: field_new
        REAL(KIND=8), INTENT(inout), DIMENSION(:,:) :: field_old
        REAL(KIND=8), INTENT(in) :: ALPHA
        ! Locals
        INTEGER :: i, j
        INTEGER :: idim1, idim2

        idim1 = SIZE(field, 1)
        idim2 = SIZE(field, 2)

        DO J=1,idim2
           DO I=1,idim1
              field_old(I,J) = field(I,J)+ &
                ALPHA*(field_new(I,J)-2.*field(I,J)+field_old(I,J))
           END DO
        END DO

      END SUBROUTINE time_smooth

      !===================================================

      SUBROUTINE compute_cu(cu, p, u)
        IMPLICIT none
        REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: cu
        REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: p
        REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: u
        ! Locals
        INTEGER :: I, J
        INTEGER :: idim1, idim2

        idim1 = SIZE(cu, 1) - 1
        idim2 = SIZE(cu, 2) - 1

         DO J=1,idim1
            DO I=2,idim2+1
               CU(I,J) = .5*(P(I,J)+P(I-1,J))*U(I,J)
            END DO
         END DO

      END SUBROUTINE compute_cu

      !===================================================

      SUBROUTINE compute_cv(cv, p, v)
        IMPLICIT none
        REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: cv
        REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: p
        REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: v
        ! Locals
        INTEGER :: I, J
        INTEGER :: idim1, idim2

        idim1 = SIZE(cv, 1) - 1
        idim2 = SIZE(cv, 2) - 1

        DO J=2,idim2+1
           DO I=1,idim1
              CV(I,J) = .5*(P(I,J)+P(I,J-1))*V(I,J)
           END DO
        END DO

      END SUBROUTINE compute_cv

      !===================================================

      SUBROUTINE compute_z(z, p, u, v, fsdx, fsdy)
        IMPLICIT none
        REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: z
        REAL(KIND=8), INTENT(in),  DIMENSION(:,:) :: p, u, v
        REAL(KIND=8), INTENT(in)                  :: fsdx, fsdy
        ! Locals
        INTEGER :: I, J
        INTEGER :: idim1, idim2

        idim1 = SIZE(z, 1) - 1
        idim2 = SIZE(z, 2) - 1

        DO J=2,idim2+1
           DO I=2,idim1+1
              Z(I,J) =(FSDX*(V(I,J)-V(I-1,J))-FSDY*(U(I,J) & 
                   -U(I,J-1)))/(P(I-1,J-1)+P(I,J-1)+P(I,J)+P(I-1,J))
           END DO
        END DO

      END SUBROUTINE compute_z

      !===================================================

      SUBROUTINE compute_h(h, p, u, v)
        IMPLICIT none
        REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: h
        REAL(KIND=8), INTENT(in),  DIMENSION(:,:) :: p, u, v
        ! Locals
        INTEGER :: I, J
        INTEGER :: idim1, idim2

        idim1 = SIZE(z, 1) - 1
        idim2 = SIZE(z, 2) - 1

        DO J=1,idim2
           DO I=1,idim1
              H(I,J) = P(I,J)+.25*(U(I+1,J)*U(I+1,J)+U(I,J)*U(I,J)     & 
                    +V(I,J+1)*V(I,J+1)+V(I,J)*V(I,J))
           END DO
        END DO

      END SUBROUTINE compute_h

      !===================================================

      SUBROUTINE compute_unew(unew, uold, z, cv, h, tdts8, tdtsdx)
        IMPLICIT none
        REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: unew
        REAL(KIND=8), INTENT(in),  DIMENSION(:,:) :: uold, z, cv, h
        REAL(KIND=8), INTENT(in)                  :: tdts8, tdtsdx
        ! Locals
        INTEGER :: I, J
        INTEGER :: idim1, idim2

        idim1 = SIZE(z, 1) - 1
        idim2 = SIZE(z, 2) - 1

        DO J=1,idim2
          DO I=2,idim1+1
            UNEW(I,J) = UOLD(I,J)+                                        &
                 TDTS8*(Z(I,J+1)+Z(I,J))*(CV(I,J+1)+CV(I-1,J+1)+CV(I-1,J) &
                        +CV(I,J))-TDTSDX*(H(I,J)-H(I-1,J))
          END DO
        END DO

      END SUBROUTINE compute_unew

      !===================================================

      SUBROUTINE compute_vnew(vnew, vold, z, cu, h, tdts8, tdtsdy)
        IMPLICIT none
        REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: vnew
        REAL(KIND=8), INTENT(in),  DIMENSION(:,:) :: vold, z, cu, h
        REAL(KIND=8), INTENT(in)                  :: tdts8, tdtsdy
        ! Locals
        INTEGER :: I, J
        INTEGER :: idim1, idim2

        idim1 = SIZE(z, 1) - 1
        idim2 = SIZE(z, 2) - 1

         DO J=2,idim2+1
            DO I=1,idim1
               VNEW(I,J) = VOLD(I,J)-TDTS8*(Z(I+1,J)+Z(I,J)) & 
                   *(CU(I+1,J)+CU(I,J)+CU(I,J-1)+CU(I+1,J-1))        & 
                   -TDTSDY*(H(I,J)-H(I,J-1))
            END DO
         END DO
       END SUBROUTINE compute_vnew

      !===================================================

      SUBROUTINE compute_pnew(pnew, pold, cu, cv, tdtsdx, tdtsdy)
        IMPLICIT none
        REAL(KIND=8), INTENT(out), DIMENSION(:,:) :: pnew
        REAL(KIND=8), INTENT(in),  DIMENSION(:,:) :: pold, cu, cv
        REAL(KIND=8), INTENT(in)                  :: tdtsdx, tdtsdy
        ! Locals
        INTEGER :: I, J
        INTEGER :: idim1, idim2

        idim1 = SIZE(z, 1) - 1
        idim2 = SIZE(z, 2) - 1

        DO J=1,idim2
          DO I=1,idim1
            PNEW(I,J) = POLD(I,J)-TDTSDX*(CU(I+1,J)-CU(I,J))   & 
                   -TDTSDY*(CV(I,J+1)-CV(I,J))
          END DO
        END DO

      END SUBROUTINE compute_pnew

    END PROGRAM shallow

    !===================================================

    ! Check error code
    subroutine check(status, text)
      implicit none
      
      integer, intent(in) :: status
      character (len=*)   :: text
    
      if (status /= 0) then
        write(6,*) "error ", status
        write(6,*) text
        stop 2
      endif

      end subroutine check
