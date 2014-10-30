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
  !     Fortran constructs and uses portable intrinsic Fortran 90
  !     timing routines
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

  use timing_mod
!$ USE omp_lib
  IMPLICIT NONE

  !      INCLUDE 'netcdf.inc'

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

  REAL(KIND=8) :: time, ptime

  ! NetCDF variables
  INTEGER :: ncid, t_id, p_id, u_id, v_id, iret, t_val
  INTEGER, DIMENSION(3) :: istart, icount 
  CHARACTER (LEN=13) :: ncfile = "shallowdat.nc"

  ! namelist input 
  CHARACTER (LEN=8) :: nml_name = "namelist" 
  INTEGER :: input_unit = 99
  INTEGER :: ierr

  !> Integer tags for timers
  INTEGER :: idxt0, idxt1

  INTEGER :: nthreads, chunk_size, myid

  !  ** Initialisations ** 
  CALL timer_init()

  !     Read in namelist 
  OPEN(unit=input_unit, file=nml_name, status='old',iostat=ierr)
  CALL check(ierr, "open "//nml_name)
  READ(unit=input_unit, nml=global_domain, iostat=ierr)
  CALL check(ierr, "read "//nml_name)
  READ(unit=input_unit, nml=io_control, iostat=ierr)
  CALL check(ierr, "read "//nml_name)

  !     Set up arrays
  m_len = m+1
  n_len = n+1

  ! OpenMP set-up
  nthreads = 1;
!$ nthreads = omp_get_max_threads();
  chunk_size = ceiling( real(M) / real(nthreads) )

  ALLOCATE( u(m_len,n_len), v(m_len,n_len), p(m_len,n_len) ) 
  ALLOCATE( unew(m_len,n_len), vnew(m_len,n_len), pnew(m_len,n_len) ) 
  ALLOCATE( uold(m_len,n_len), vold(m_len,n_len), pold(m_len,n_len) )
  ALLOCATE( cu(m_len,n_len), cv(m_len,n_len) ) 
  ALLOCATE( z(m_len,n_len), h(m_len,n_len), psi(m_len,n_len) ) 

  !     Prepare netCDF file to receive model output data
  IF (l_out) THEN 
     call netcdf_setup(ncfile,m,n,ncid,t_id,p_id,u_id,v_id,istart,icount)
  ENDIF

  !     NOTE BELOW THAT TWO DELTA T (TDT) IS SET TO DT ON THE FIRST
  !     CYCLE AFTER WHICH IT IS RESET TO DT+DT.
  DT = 90.
  TDT = DT

  DX = 1.E5
  DY = 1.E5
  A = 1.E6
  ALPHA = .001

  MP1 = M+1
  NP1 = N+1
  EL = N*DX
  PI = 4.*ATAN(1.)
  TPI = PI+PI
  DI = TPI/M
  DJ = TPI/N
  PCF = PI*PI*A*A/(EL*EL)

  !     INITIAL VALUES OF THE STREAM FUNCTION AND P
!$OMP PARALLEL DO default(none), shared(NP1,MP1,psi,p,A,DI,DJ,PCF), &
!$OMP             private(i,j)
  DO J=1,NP1
     DO I=1,MP1
        PSI(I,J) = A*SIN((I-.5)*DI)*SIN((J-.5)*DJ)
        P(I,J) = PCF*(COS(2.*(I-1)*DI)   & 
             +COS(2.*(J-1)*DJ))+50000.
     END DO
  END DO
!$OMP END PARALLEL DO

  !     INITIALIZE VELOCITIES
!$OMP PARALLEL DO default(none), shared(N,M,psi,Dx,Dy,U,V), &
!$OMP             private(i,j)
  DO J=1,N
     DO I=1,M
        U(I+1,J) = -(PSI(I+1,J+1)-PSI(I+1,J))/DY
        V(I,J+1) = (PSI(I+1,J+1)-PSI(I,J+1))/DX
     END DO
  END DO
!$OMP END PARALLEL DO

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

!$OMP PARALLEL DO default(none), shared(NP1,MP1,UOLD,VOLD,POLD,U,V,P), &
!$OMP             private(i,j)
  DO J=1,NP1
     DO I=1,MP1
        UOLD(I,J) = U(I,J)
        VOLD(I,J) = V(I,J)
        POLD(I,J) = P(I,J)
     END DO
  END DO
!$OMP END PARALLEL DO

  !     PRINT INITIAL VALUES
  IF (l_out) THEN 
     WRITE(6,390) N,M,DX,DY,DT,ALPHA
390  FORMAT(" NUMBER OF POINTS IN THE X DIRECTION",I8,/    & 
          " NUMBER OF POINTS IN THE Y DIRECTION",I8,/    & 
          " GRID SPACING IN THE X DIRECTION    ",F8.0,/  & 
          " GRID SPACING IN THE Y DIRECTION    ",F8.0,/  & 
          " TIME STEP                          ",F8.0,/  & 
          " TIME FILTER PARAMETER              ",F8.3)
     MNMIN = MIN0(M,N)
     WRITE(6,391) (P(I,I),I=1,MNMIN)
391  FORMAT(/' INITIAL DIAGONAL ELEMENTS OF P ' //,(8E15.6))
     WRITE(6,392) (U(I,I),I=1,MNMIN)
392  FORMAT(/' INITIAL DIAGONAL ELEMENTS OF U ' //,(8E15.6))
     WRITE(6,393) (V(I,I),I=1,MNMIN)
393  FORMAT(/' INITIAL DIAGONAL ELEMENTS OF V ' //,(8E15.6))

     !        Write intial values of p, u, and v into a netCDF file   
     t_val = 0   
     call my_ncwrite(ncid,p_id,istart,icount,p(1:m,1:n),m,n,t_id,t_val)
     call my_ncwrite(ncid,u_id,istart,icount,u(1:m,1:n),m,n,t_id,t_val)
     call my_ncwrite(ncid,v_id,istart,icount,v(1:m,1:n),m,n,t_id,t_val)
  ENDIF

  FSDX = 4./DX
  FSDY = 4./DY

  !     Start timer
  TIME = 0.

  !     Start timer
  CALL timer_start('Time-stepping',idxt0)

!$OMP PARALLEL default(none), &
!$OMP         shared(M, N, MP1, NP1, itmax, fsdx, fsdy, dx, dy, cu, cv,       &
!$OMP                u, v, z, h, p, uold, vold, pold, unew, vnew, pnew,       &
!$OMP                time, dt, alpha, l_out, mprint, ptime, mnmin, istart,    &
!$oMP                t_val, t_id, p_id, u_id, v_id, ncid, icount, chunk_size),&
!$OMP          private(i, j, ncycle, tdts8, tdtsdx, tdtsdy, idxt1, myid),     &
!$OMP          firstprivate(tdt)

  myid = 0
!$  myid = omp_get_thread_num()

  !  ** Start of time loop ** 
  DO ncycle=1,itmax

     !        COMPUTE CAPITAL U, CAPITAL V, Z AND H

     CALL timer_start('Compute CU,CV,CZ,H',idxt1)

!$OMP DO SCHEDULE(static,chunk_size)
     DO J=1,N
        DO I=1,M
           CU(I+1,J) = .5*(P(I+1,J)+P(I,J))*U(I+1,J)
        END DO
     END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
     DO J=1,N
        DO I=1,M
           CV(I,J+1) = .5*(P(I,J+1)+P(I,J))*V(I,J+1)
        END DO
     END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
     DO J=1,N
        DO I=1,M
           Z(I+1,J+1) =(FSDX*(V(I+1,J+1)-V(I,J+1))-FSDY*(U(I+1,J+1) & 
                -U(I+1,J)))/(P(I,J)+P(I+1,J)+P(I+1,J+1)+P(I,J+1))
        END DO
     END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
     DO J=1,N
        DO I=1,M
           H(I,J) = P(I,J)+.25*(U(I+1,J)*U(I+1,J)+U(I,J)*U(I,J)     & 
                +V(I,J+1)*V(I,J+1)+V(I,J)*V(I,J))
        END DO
     END DO
!$OMP END DO nowait

     CALL timer_stop(idxt1)
!$OMP BARRIER

     !        PERIODIC CONTINUATION
     CALL timer_start('PBCs',idxt1)
!$OMP SINGLE
     DO J=1,N
        CU(1,J) = CU(M+1,J)
        CV(M+1,J+1) = CV(1,J+1)
        Z(1,J+1) = Z(M+1,J+1)
        H(M+1,J) = H(1,J)
     END DO
     CU(1,N+1) = CU(M+1,1)
     CV(M+1,1) = CV(1,N+1)
     Z(1,1) = Z(M+1,N+1)
     H(M+1,N+1) = H(1,1)
!$OMP END SINGLE

!$OMP DO SCHEDULE(static,chunk_size)
     DO I=1,M
        CU(I+1,N+1) = CU(I+1,1)
        CV(I,1) = CV(I,N+1)
        Z(I+1,1) = Z(I+1,N+1)
        H(I,N+1) = H(I,1)
     END DO
!$OMP END DO nowait

     CALL timer_stop(idxt1)
!$OMP BARRIER

     !        COMPUTE NEW VALUES U,V AND P
     TDTS8 = TDT/8.
     TDTSDX = TDT/DX
     TDTSDY = TDT/DY

     CALL timer_start('Compute {U,V,P}NEW',idxt1)

!$OMP DO SCHEDULE(static,chunk_size)
     DO J=1,N
        DO I=1,M
           UNEW(I+1,J) = UOLD(I+1,J)+                                     &
                TDTS8*(Z(I+1,J+1)+Z(I+1,J))*(CV(I+1,J+1)+CV(I,J+1)+CV(I,J) &
                +CV(I+1,J))-TDTSDX*(H(I+1,J)-H(I,J))                       
        END DO
     END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
     DO J=1,N
        DO I=1,M
          VNEW(I,J+1) = VOLD(I,J+1)-TDTS8*(Z(I+1,J+1)+Z(I,J+1)) & 
                *(CU(I+1,J+1)+CU(I,J+1)+CU(I,J)+CU(I+1,J))        & 
                -TDTSDY*(H(I,J+1)-H(I,J))
        END DO
     END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
     DO J=1,N
        DO I=1,M
           PNEW(I,J) = POLD(I,J)-TDTSDX*(CU(I+1,J)-CU(I,J))   & 
                -TDTSDY*(CV(I,J+1)-CV(I,J))
        END DO
     END DO
!$OMP END DO nowait

     CALL timer_stop(idxt1)
!$OMP BARRIER

     CALL timer_start('PBCs',idxt1)
     !        PERIODIC CONTINUATION
!$OMP SINGLE
     DO J=1,N
        UNEW(1,J) = UNEW(M+1,J)
        VNEW(M+1,J+1) = VNEW(1,J+1)
        PNEW(M+1,J) = PNEW(1,J)
     END DO
     UNEW(1,N+1) = UNEW(M+1,1)
     VNEW(M+1,1) = VNEW(1,N+1)
     PNEW(M+1,N+1) = PNEW(1,1)
!$OMP END SINGLE

!$OMP DO SCHEDULE(static,chunk_size)
     DO I=1,M
        UNEW(I+1,N+1) = UNEW(I+1,1)
        VNEW(I,1) = VNEW(I,N+1)
        PNEW(I,N+1) = PNEW(I,1)
     END DO
!$OMP END DO

     CALL timer_stop(idxt1)

!$OMP SINGLE
     TIME = TIME + DT
!$OMP END SINGLE

     !        TIME SMOOTHING AND UPDATE FOR NEXT CYCLE
     IF(NCYCLE .GT. 1) then

        CALL timer_start('Time smooth',idxt1)

!$OMP DO SCHEDULE(static,chunk_size)
        DO J=1,NP1
           DO I=1,MP1
              UOLD(I,J) = U(I,J)+ALPHA*(UNEW(I,J)-2.*U(I,J)+UOLD(I,J))
           END DO
        END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
        DO J=1,NP1
           DO I=1,MP1
              VOLD(I,J) = V(I,J)+ALPHA*(VNEW(I,J)-2.*V(I,J)+VOLD(I,J))
           END DO
        END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
        DO J=1,NP1
           DO I=1,MP1
              POLD(I,J) = P(I,J)+ALPHA*(PNEW(I,J)-2.*P(I,J)+POLD(I,J))
           END DO
        END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
        DO J=1,NP1
           DO I=1,MP1
              U(I,J) = UNEW(I,J)
           END DO
        END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
        DO J=1,NP1
           DO I=1,MP1
              V(I,J) = VNEW(I,J)
           END DO
        END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(static,chunk_size)
        DO J=1,NP1
           DO I=1,MP1
              P(I,J) = PNEW(I,J)
           END DO
        END DO
!$OMP END DO nowait

        CALL timer_stop(idxt1)
!$OMP BARRIER

     else

        TDT = TDT+TDT

!$OMP DO SCHEDULE(static,chunk_size)
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

     endif

     IF( l_out .AND. (MOD(NCYCLE,MPRINT) .EQ. 0) .AND. (myid .EQ. 0)) then


        PTIME = TIME/3600.
        WRITE(6,350) NCYCLE,PTIME
350     FORMAT(//' CYCLE NUMBER',I5,' MODEL TIME IN  HOURS', F6.2)
        WRITE(6,355) (PNEW(I,I),I=1,MNMIN)
355     FORMAT(/' DIAGONAL ELEMENTS OF P ' //(8E15.6))
        WRITE(6,360) (UNEW(I,I),I=1,MNMIN)
360     FORMAT(/' DIAGONAL ELEMENTS OF U ' //(8E15.6))
        WRITE(6,365) (VNEW(I,I),I=1,MNMIN)
365     FORMAT(/' DIAGONAL ELEMENTS OF V ' //(8E15.6))

        !           Append calculated values of p, u, and v to netCDF file
        istart(3) = ncycle/mprint + 1
        t_val = ncycle

        !           Shape of record to be written (one ncycle at a time)
        call my_ncwrite(ncid,p_id,istart,icount,p(1:m,1:n),m,n,t_id,t_val)
        call my_ncwrite(ncid,u_id,istart,icount,u(1:m,1:n),m,n,t_id,t_val)
        call my_ncwrite(ncid,v_id,istart,icount,v(1:m,1:n),m,n,t_id,t_val)

     endif

  End do

!$OMP END PARALLEL

  !  ** End of time loop ** 
  CALL timer_stop(idxt0)

  WRITE(6,"('P CHECKSUM after ',I6,' steps = ',E15.7)") &
       itmax, SUM(ABS(PNEW(1:M,1:N)))
  WRITE(6,"('U CHECKSUM after ',I6,' steps = ',E15.7)") &
       itmax,SUM(ABS(UNEW(2:M+1,1:N)))
  WRITE(6,"('V CHECKSUM after ',I6,' steps = ',E15.7)") &
       itmax,SUM(abs(VNEW(1:M,2:N+1)))

  !     Close the netCDF file

  IF (l_out) THEN 

     iret = 0
     !iret = nf_close(ncid)
     call check_err(iret)
  ENDIF

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


      subroutine check_err(iret)
      integer iret
!      include 'netcdf.inc'
!      if(iret .ne. NF_NOERR) then
      if(iret .ne. 0) then
         !print *, nf_strerror(iret)
         stop
      endif
      end subroutine


      subroutine netcdf_setup(file,m,n,ncid,t_id,p_id,u_id,v_id,istart,icount)
!     Input args: file, m, n
!     Output args: ncid,t_id,p_id,u_id,v_id,istart,icount)
      character(len=*) file
      integer m,n
!     declarations for netCDF library
!      include 'netcdf.inc'
!     error status return
      integer iret
!     netCDF id
      integer ncid
!     dimension ids
      integer m_dim
      integer n_dim
      integer time_dim      
!     variable ids
      integer t_id
      integer p_id
      integer u_id
      integer v_id
!     rank (number of dimensions) for each variable
      integer p_rank, u_rank, v_rank
      parameter (p_rank = 3)
      parameter (u_rank = 3)
      parameter (v_rank = 3)
!     variable shapes
      integer t_dims(1)
      integer p_dims(p_rank)
      integer u_dims(u_rank)
      integer v_dims(v_rank)
      integer istart(p_rank)
      integer icount(p_rank)
      
      iret = 0
!     enter define mode
      !iret = nf_create(file, NF_CLOBBER,ncid)
      call check_err(iret)
!     define dimensions
      !iret = nf_def_dim(ncid, 'm', m, m_dim)
      call check_err(iret)
      !iret = nf_def_dim(ncid, 'n', n, n_dim)
      call check_err(iret)
!     time is an unlimited dimension so that any number of
!     records can be added
      !iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)
!     define coordinate variable for time      
      t_dims(1) = time_dim
      !iret = nf_def_var(ncid, 'time', NF_INT, 1, t_dims, t_id)
      call check_err(iret)
!     define variables
      p_dims(1) = m_dim
      p_dims(2) = n_dim
      p_dims(3) = time_dim
      !iret = nf_def_var(ncid, 'p', NF_DOUBLE, p_rank, p_dims, p_id)
      call check_err(iret)
      u_dims(1) = m_dim
      u_dims(2) = n_dim
      u_dims(3) = time_dim
      !iret = nf_def_var(ncid, 'u', NF_DOUBLE, u_rank, u_dims, u_id)
      call check_err(iret)
      v_dims(1) = m_dim
      v_dims(2) = n_dim
      v_dims(3) = time_dim
      !iret = nf_def_var(ncid, 'v', NF_DOUBLE, v_rank, v_dims, v_id)
      call check_err(iret)
!     start netCDF write at the first (1,1,1) position of the array
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
!     shape of record to be written (one ncycle at a time)
      icount(1) = m
      icount(2) = n
      icount(3) = 1
      
!     leave define mode
      !iret = nf_enddef(ncid)
      call check_err(iret)
      
!     end of netCDF definitions
      end subroutine netcdf_setup
 
     
      subroutine my_ncwrite(id,varid,istart,icount,var,m,n,t_id,t_val)
!     Input args: id, varid,istart,icount,var,m,n,t_id,t_val
!     Write a whole array out to the netCDF file
      integer id,varid,iret
      integer icount(3)
      integer istart(3)
      integer m,n
      real(kind=8) var (m,n)
      integer t_id,t_val
      integer t_start(1), t_count(1)

      iret = 0
      !iret = nf_put_vara_double(id,varid,istart,icount,var)
      call check_err(iret)
      
      t_start(1) = istart(3) 
      t_count(1) = 1
      !iret = nf_put_vara_int(id,t_id,t_start,t_count,t_val)
      call check_err(iret)

      end subroutine my_ncwrite

