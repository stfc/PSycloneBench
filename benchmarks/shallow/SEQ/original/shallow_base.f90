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


      USE dl_timer
      use shallow_io_mod, only: ascii_write
      IMPLICIT NONE

      INTEGER :: m, n    ! global domain size
      INTEGER :: itmax   ! number of timesteps
      INTEGER :: mprint  ! frequency of output    
      NAMELIST/global_domain/ m, n, itmax, mprint

      LOGICAL :: l_out   ! produce output  
      NAMELIST/io_control/ l_out

      INTEGER :: m_len, n_len       ! array sizes
      INTEGER :: mp1, np1           ! m+1 and n+1

      ! solution arrays
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) ::                & 
                             u, v, p, unew, vnew, pnew,           & 
                             uold, vold, pold, cu, cv, z, h, psi  

      REAL(KIND=8) :: dt, tdt, dx, dy, a, alpha, el, pi, tpi, di, dj, pcf, & 
                      tdts8, tdtsdx, tdtsdy, fsdx, fsdy
      INTEGER :: mnmin, ncycle
      INTEGER :: i, j
   
      ! timer variables 
      !REAL(KIND=8) :: mfs100, mfs200, mfs300, mfs310, & 
      !                t100, t200, t300, t310,         & 
      !                tstart, ctime, tcyc
      !INTEGER :: c1, c2, r, max
      REAL(KIND=8) :: time, ptime

      ! namelist input 
      CHARACTER (LEN=8) :: nml_name = "namelist" 
      INTEGER :: input_unit = 99
      INTEGER :: ierr

      !> Integer tags for timers
      INTEGER :: idxt0, cu_timer, pbcs_timer, unew_timer, tsmooth_timer

!  ** Initialisations ** 
      call timer_init()
      
!     Read in namelist 
      OPEN(unit=input_unit, file=nml_name, status='old',iostat=ierr)
      CALL check(ierr, "open "//nml_name)
      READ(unit=input_unit, nml=global_domain, iostat=ierr)
      CALL check(ierr, "read "//nml_name)
      READ(unit=input_unit, nml=io_control, iostat=ierr)
      CALL check(ierr, "read "//nml_name)

      call timer_register(idxt0,     label='Time-stepping', &
                          num_repeats=int(itmax,8))
      call timer_register(cu_timer,  label='Compute CU,CV,CZ,H')
      call timer_register(pbcs_timer,label='PBCs')
      call timer_register(unew_timer,label='Compute {U,V,P}NEW')
      call timer_register(tsmooth_timer,label='Time smooth')

!     Set up arrays
      m_len = m+1
      n_len = n+1

      ALLOCATE( u(m_len,n_len), v(m_len,n_len), p(m_len,n_len) ) 
      ALLOCATE( unew(m_len,n_len), vnew(m_len,n_len), pnew(m_len,n_len) ) 
      ALLOCATE( uold(m_len,n_len), vold(m_len,n_len), pold(m_len,n_len) )
      ALLOCATE( cu(m_len,n_len), cv(m_len,n_len) ) 
      ALLOCATE( z(m_len,n_len), h(m_len,n_len), psi(m_len,n_len) ) 
 
!     NOTE BELOW THAT TWO DELTA T (TDT) IS SET TO DT ON THE FIRST
!     CYCLE AFTER WHICH IT IS RESET TO DT+DT.
      DT = 90.
      TDT = DT
 
      DX = 1.0D5
      DY = 1.0D5
      A = 1.0D6
      ALPHA = 0.001d0

      MP1 = M+1
      NP1 = N+1
      EL = N*DX
      PI = 4.0d0*ATAN(1.0d0)
      TPI = PI+PI
      DI = TPI/M
      DJ = TPI/N
      PCF = PI*PI*A*A/(EL*EL)
     
!     INITIAL VALUES OF THE STREAM FUNCTION AND P
      DO J=1,NP1
         DO I=1,MP1
            PSI(I,J) = A*SIN((I-.5d0)*DI)*SIN((J-.5d0)*DJ)
            P(I,J) = PCF*(COS(2.0d0*(I-1)*DI)   & 
                 +COS(2.0d0*(J-1)*DJ))+50000.0d0
         END DO
      END DO

!     INITIALIZE VELOCITIES
      DO J=1,N
         DO I=1,M
            U(I+1,J) = -(PSI(I+1,J+1)-PSI(I+1,J))/DY
            V(I,J+1) = (PSI(I+1,J+1)-PSI(I,J+1))/DX
         END DO
      END DO
     
!     PERIODIC CONTINUATION
      DO J=1,N
         U(1,J) = U(M+1,J)
         V(M+1,J+1) = V(1,J+1)
      END DO
      DO I=1,M
         U(I+1,N+1) = U(I+1,1)
         V(I,1) = V(I,N+1)
      END DO
      U(1,N+1) = U(M+1,1)
      V(M+1,1) = V(1,N+1)
      DO J=1,NP1
         DO I=1,MP1
            UOLD(I,J) = U(I,J)
            VOLD(I,J) = V(I,J)
            POLD(I,J) = P(I,J)
         END DO
      END DO
     
!     PRINT INITIAL VALUES
      IF (l_out) THEN 
         WRITE(6,390) N,M,DX,DY,DT,ALPHA
 390     FORMAT(" NUMBER OF POINTS IN THE X DIRECTION",I8,/    & 
                " NUMBER OF POINTS IN THE Y DIRECTION",I8,/    & 
                " GRID SPACING IN THE X DIRECTION    ",F8.0,/  & 
                " GRID SPACING IN THE Y DIRECTION    ",F8.0,/  & 
                " TIME STEP                          ",F8.0,/  & 
                " TIME FILTER PARAMETER              ",F8.3)
         MNMIN = MIN0(M,N)
         !WRITE(6,391) (P(I,I),I=1,MNMIN)
! 391     FORMAT(/' INITIAL DIAGONAL ELEMENTS OF P ' //,(8E15.6))
         !WRITE(6,392) (U(I,I),I=1,MNMIN)
! 392     FORMAT(/' INITIAL DIAGONAL ELEMENTS OF U ' //,(8E15.6))
         !WRITE(6,393) (V(I,I),I=1,MNMIN)
! 393     FORMAT(/' INITIAL DIAGONAL ELEMENTS OF V ' //,(8E15.6))

         ! Generate and output checksums of initial fields
         write(*,"('psi initial CHECKSUM = ',E24.16)") &
              sum( abs(psi(2:m+1,2:n+1)) )
         write(*,"('P initial CHECKSUM = ',E24.16)")   &
              sum( abs(p(1:m,1:n)) )
         write(*,"('U initial CHECKSUM = ',E24.16)")   &
              sum( abs(u(2:M+1,1:N)) )
         write(*,"('V initial CHECKSUM = ',E24.16)")   &
              sum( abs(v(1:m,2:N+1)) )

!        Write intial values of p, u, and v into a netCDF file   
         call ascii_write(0,'u_array.dat'  ,u  ,m,n,2,1)
         call ascii_write(0,'cu_array.dat' ,cu ,m,n,2,1)
         call ascii_write(0,'v_array.dat'  ,v  ,m,n,1,2)
         call ascii_write(0,'p_array.dat'  ,p  ,m,n,1,1)
         call ascii_write(0,'psi_array.dat',psi,m,n,2,2)
         call ascii_write(0,'cv_array.dat' ,cv ,m,n,1,2)
      ENDIF

!     Start timer
!      call system_clock (count=c1, count_rate=r, count_max=max)
!      TSTART = c1
!      T300 = 1.
!      T310 = 1.
      TIME = 0.

      !     Start timer
      CALL timer_start(idxt0)

!  ** Start of time loop ** 
      DO ncycle=1,itmax
    
!        COMPUTE CAPITAL U, CAPITAL V, Z AND H
         FSDX = 4.0d0/DX
         FSDY = 4.0d0/DY

         !call system_clock(count=c1, count_rate=r,count_max=max)
         !T100 = c1
         CALL timer_start(cu_timer)

         DO J=1,N
            DO I=1,M
               CU(I+1,J) = .5d0*(P(I+1,J)+P(I,J))*U(I+1,J)
               CV(I,J+1) = .5d0*(P(I,J+1)+P(I,J))*V(I,J+1)
               Z(I+1,J+1) =(FSDX*(V(I+1,J+1)-V(I,J+1))-FSDY*(U(I+1,J+1) & 
                    -U(I+1,J)))/(P(I,J)+P(I+1,J)+P(I+1,J+1)+P(I,J+1))
               H(I,J) = P(I,J)+0.25d0*(U(I+1,J)*U(I+1,J)+U(I,J)*U(I,J)  & 
                    +V(I,J+1)*V(I,J+1)+V(I,J)*V(I,J))
            END DO
         END DO

         CALL timer_stop(cu_timer)
         !call system_clock(count=c2,count_rate=r,count_max=max)
         !T100 = dble(c2-T100)/dble(r)

         CALL timer_start(pbcs_timer)
!        PERIODIC CONTINUATION
         DO J=1,N
            CU(1,J) = CU(M+1,J)
            CV(M+1,J+1) = CV(1,J+1)
            Z(1,J+1) = Z(M+1,J+1)
            H(M+1,J) = H(1,J)
         END DO
         DO I=1,M
            CU(I+1,N+1) = CU(I+1,1)
            CV(I,1) = CV(I,N+1)
            Z(I+1,1) = Z(I+1,N+1)
            H(I,N+1) = H(I,1)
         END DO
         CU(1,N+1) = CU(M+1,1)
         CV(M+1,1) = CV(1,N+1)
         Z(1,1) = Z(M+1,N+1)
         H(M+1,N+1) = H(1,1)
     
         CALL timer_stop(pbcs_timer)

!        COMPUTE NEW VALUES U,V AND P
         TDTS8 = TDT/8.0d0
         TDTSDX = TDT/DX
         TDTSDY = TDT/DY

         !call system_clock(count=c1, count_rate=r, count_max=max)
         !T200 = c1
         CALL timer_start(unew_timer)

         DO J=1,N
            DO I=1,M
               UNEW(I+1,J) = UOLD(I+1,J)+                                     &
                   TDTS8*(Z(I+1,J+1)+Z(I+1,J))*(CV(I+1,J+1)+CV(I,J+1)+CV(I,J) &
                   +CV(I+1,J))-TDTSDX*(H(I+1,J)-H(I,J))                       
               VNEW(I,J+1) = VOLD(I,J+1)-TDTS8*(Z(I+1,J+1)+Z(I,J+1)) & 
                   *(CU(I+1,J+1)+CU(I,J+1)+CU(I,J)+CU(I+1,J))        & 
                   -TDTSDY*(H(I,J+1)-H(I,J))
               PNEW(I,J) = POLD(I,J)-TDTSDX*(CU(I+1,J)-CU(I,J))   & 
                   -TDTSDY*(CV(I,J+1)-CV(I,J))
            END DO
         END DO

         CALL timer_stop(unew_timer)
         !call system_clock(count=c2, count_rate=r, count_max=max)
         !T200 = dble(c2 -T200)/dble(r)

         CALL timer_start(pbcs_timer)
!        PERIODIC CONTINUATION
         DO J=1,N
            UNEW(1,J) = UNEW(M+1,J)
            VNEW(M+1,J+1) = VNEW(1,J+1)
            PNEW(M+1,J) = PNEW(1,J)
         END DO
         DO I=1,M
            UNEW(I+1,N+1) = UNEW(I+1,1)
            VNEW(I,1) = VNEW(I,N+1)
            PNEW(I,N+1) = PNEW(I,1)
         END DO
         UNEW(1,N+1) = UNEW(M+1,1)
         VNEW(M+1,1) = VNEW(1,N+1)
         PNEW(M+1,N+1) = PNEW(1,1)

         CALL timer_stop(pbcs_timer)

         TIME = TIME + DT

         IF( l_out .AND. (MOD(NCYCLE,MPRINT) .EQ. 0) ) then
            
            PTIME = TIME/3600.
            WRITE(6,350) NCYCLE,PTIME
 350        FORMAT(//' CYCLE NUMBER',I5,' MODEL TIME IN  HOURS', F6.2)
            !WRITE(6,355) (PNEW(I,I),I=1,MNMIN)
! 355        FORMAT(/' DIAGONAL ELEMENTS OF P ' //(8E15.6))
            !WRITE(6,360) (UNEW(I,I),I=1,MNMIN)
! 360        FORMAT(/' DIAGONAL ELEMENTS OF U ' //(8E15.6))
            !WRITE(6,365) (VNEW(I,I),I=1,MNMIN)
! 365        FORMAT(/' DIAGONAL ELEMENTS OF V ' //(8E15.6))

!           jr added MFS310--don't know what actual mult factor should be
!           jr changed divide by 1 million to divide by 100K since system_clock
!           jr resolution is millisec rather than cpu_time's 10 millisec
            !MFS310 = 0.0
            !MFS100 = 0.0
            !MFS200 = 0.0
            !MFS300 = 0.0
            !IF (T310 .GT. 0) MFS310 = 24.*M*N/T310/1.D5
            !IF (T100 .GT. 0) MFS100 = 24.*M*N/T100/1.D5
            !IF (T200 .GT. 0) MFS200 = 26.*M*N/T200/1.D5
            !IF (T300 .GT. 0) MFS300 = 15.*M*N/T300/1.D5

            !call system_clock(count=c2, count_rate=r,count_max=max)
            !CTIME = dble(c2 - TSTART)/dble(r)
            !TCYC = CTIME/FLOAT(NCYCLE)

            !WRITE(6,375) NCYCLE,CTIME,TCYC,T310,MFS310,T200,MFS200,T300,MFS300
! 375        FORMAT(' CYCLE NUMBER',I5,' TOTAL COMPUTER TIME', D15.6,   & 
!                   ' TIME PER CYCLE', D15.6, /                           & 
!                   ' TIME AND MEGAFLOPS FOR LOOP 310 ', D15.6,2x,D6.1/   & 
!                   ' TIME AND MEGAFLOPS FOR LOOP 200 ', D15.6,2x,D6.1/   & 
!                   ' TIME AND MEGAFLOPS FOR LOOP 300 ', D15.6,2x,D6.1/ )

            call ascii_write(ncycle,'u_array.dat'  ,u  ,m,n,2,1)
            call ascii_write(ncycle,'cu_array.dat' ,cu ,m,n,2,1)
            call ascii_write(ncycle,'cv_array.dat' ,cv ,m,n,1,2)
            call ascii_write(ncycle,'v_array.dat'  ,v  ,m,n,1,2)
            call ascii_write(ncycle,'p_array.dat'  ,p  ,m,n,1,1)
            call ascii_write(ncycle,'z_array.dat'  ,z  ,m,n,2,2)
            call ascii_write(ncycle,'h_array.dat'  ,h  ,m,n,1,1)

         endif

!        Write out time if last timestep
!         IF (ncycle .EQ. itmax) THEN 
!            call system_clock(count=c2, count_rate=r,count_max=max)
!            CTIME = dble(c2 - TSTART)/dble(r)
!            WRITE(6,376) ctime  
! 376        FORMAT('system_clock time ', F15.6)
!         ENDIF


!        TIME SMOOTHING AND UPDATE FOR NEXT CYCLE
         IF(NCYCLE .GT. 1) then

            !call system_clock(count=c1,count_rate=r,count_max=max)
            !T300 = c1
            CALL timer_start(tsmooth_timer)

            DO J=1,N
               DO I=1,M
                  UOLD(I,J) = U(I,J)+ALPHA*(UNEW(I,J)-2.0d0*U(I,J)+UOLD(I,J))
                  VOLD(I,J) = V(I,J)+ALPHA*(VNEW(I,J)-2.0d0*V(I,J)+VOLD(I,J))
                  POLD(I,J) = P(I,J)+ALPHA*(PNEW(I,J)-2.0d0*P(I,J)+POLD(I,J))
                  U(I,J) = UNEW(I,J)
                  V(I,J) = VNEW(I,J)
                  P(I,J) = PNEW(I,J)
               END DO
            END DO

            CALL timer_stop(tsmooth_timer)
            !call system_clock(count=c2,count_rate=r, count_max=max)
            !T300 = dble(c2 - T300)/dble(r)
    
!           PERIODIC CONTINUATION
            DO J=1,N
               UOLD(M+1,J) = UOLD(1,J)
               VOLD(M+1,J) = VOLD(1,J)
               POLD(M+1,J) = POLD(1,J)
               U(M+1,J) = U(1,J)
               V(M+1,J) = V(1,J)
               P(M+1,J) = P(1,J)
            END DO
            DO I=1,M
               UOLD(I,N+1) = UOLD(I,1)
               VOLD(I,N+1) = VOLD(I,1)
               POLD(I,N+1) = POLD(I,1)
               U(I,N+1) = U(I,1)
               V(I,N+1) = V(I,1)
               P(I,N+1) = P(I,1)
            END DO
            UOLD(M+1,N+1) = UOLD(1,1)
            VOLD(M+1,N+1) = VOLD(1,1)
            POLD(M+1,N+1) = POLD(1,1)
            U(M+1,N+1) = U(1,1)
            V(M+1,N+1) = V(1,1)
            P(M+1,N+1) = P(1,1)
         else

            TDT = TDT+TDT

            !call system_clock(count=c1, count_rate=r,count_max=max)
            !T310 = c1

            DO J=1,NP1
               DO I=1,MP1
                  UOLD(I,J) = U(I,J)
                  VOLD(I,J) = V(I,J)
                  POLD(I,J) = P(I,J)
                  U(I,J) = UNEW(I,J)
                  V(I,J) = VNEW(I,J)
                  P(I,J) = PNEW(I,J)
               END DO
            END DO

            !call system_clock(count=c2, count_rate=r, count_max=max)
            !T310 = dble(c2 - T310)/dble(r)

         endif

     End do

!  ** End of time loop ** 
      CALL timer_stop(idxt0)

      WRITE(6,"('P CHECKSUM after ',I6,' steps = ',E24.16)") &
           itmax, SUM( ABS(PNEW(1:m,1:n)) )
      WRITE(6,"('U CHECKSUM after ',I6,' steps = ',E24.16)") &
           itmax, SUM( ABS(UNEW(2:M+1,1:N)) )
      WRITE(6,"('V CHECKSUM after ',I6,' steps = ',E24.16)") &
           itmax, SUM( ABS(VNEW(1:m,2:N+1)) )

      CALL timer_report()

!     Free memory
      DEALLOCATE( u, v, p, unew, vnew, pnew, uold, vold, pold )
      DEALLOCATE( cu, cv, z, h, psi ) 

      END PROGRAM shallow

      
!     Check error code
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
