MODULE shallow_io
  IMPLICIT none

  INTEGER :: m, n      !< global domain size
  INTEGER :: mp1, np1  !< m+1 and n+1 == array extents

  INTEGER :: itmax   !< number of timesteps
  INTEGER :: mprint  !< frequency of output    
  LOGICAL :: l_out   !< produce output  

  NAMELIST/global_domain/ m, n, itmax, mprint
  NAMELIST/io_control/ l_out

  ! NetCDF variables
  INCLUDE 'netcdf.inc'

  INTEGER :: ncid, t_id, p_id, u_id, v_id, iret, t_val
  INTEGER, DIMENSION(3) :: istart, icount 
  CHARACTER (LEN=13) :: ncfile = "shallowdat.nc"

CONTAINS

  !> Reads the namelist file for user-specified control 
  !! parameters
  SUBROUTINE read_namelist
    IMPLICIT none

    ! namelist input 
    CHARACTER (LEN=8) :: nml_name = "namelist" 
    INTEGER :: input_unit = 99
    INTEGER :: ierr

    !     Read in namelist 
    OPEN(unit=input_unit, file=nml_name, status='old',iostat=ierr)
    CALL check(ierr, "open "//nml_name)
    READ(unit=input_unit, nml=global_domain, iostat=ierr)
    CALL check(ierr, "read "//nml_name)
    READ(unit=input_unit, nml=io_control, iostat=ierr)
    CALL check(ierr, "read "//nml_name)

  END SUBROUTINE read_namelist

  !===================================================

  SUBROUTINE print_initial_values(n,m,dx,dy,dt,alpha,p,u,v)
    IMPLICIT none
    INTEGER,      INTENT(in) :: n, m
    REAL(KIND=8), INTENT(in) :: dx, dy, dt, alpha
    REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: p, u, v

    WRITE(6,390) N,M,DX,DY,DT,ALPHA
 390     FORMAT(" NUMBER OF POINTS IN THE X DIRECTION",I8,/    & 
                " NUMBER OF POINTS IN THE Y DIRECTION",I8,/    & 
                " GRID SPACING IN THE X DIRECTION    ",F8.0,/  & 
                " GRID SPACING IN THE Y DIRECTION    ",F8.0,/  & 
                " TIME STEP                          ",F8.0,/  & 
                " TIME FILTER PARAMETER              ",F8.3)

    CALL print_diagonals(p, u, v)

  END SUBROUTINE print_initial_values

  !===================================================

  SUBROUTINE print_diagonals(p, u, v)
    IMPLICIT none
    REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: p, u, v
    ! Locals
    INTEGER :: i, mnmin

    MNMIN = MIN0(M,N)

    WRITE(6,391) (P(I,I),I=1,MNMIN)
391 FORMAT(/' INITIAL DIAGONAL ELEMENTS OF P ' //,(8E15.6))
    WRITE(6,392) (U(I,I),I=1,MNMIN)
392 FORMAT(/' INITIAL DIAGONAL ELEMENTS OF U ' //,(8E15.6))
    WRITE(6,393) (V(I,I),I=1,MNMIN)
393 FORMAT(/' INITIAL DIAGONAL ELEMENTS OF V ' //,(8E15.6))

  END SUBROUTINE print_diagonals

  !===================================================
  
  SUBROUTINE netcdf_setup(file,m,n,ncid,t_id,p_id,u_id,v_id,istart,icount)
    !     Input args: file, m, n
    !     Output args: ncid,t_id,p_id,u_id,v_id,istart,icount)
    character(len=*) file
    integer m,n
    !     declarations for netCDF library
    include 'netcdf.inc'
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
      
    !     enter define mode
    iret = nf_create(file, NF_CLOBBER,ncid)
    call check_err(iret)
    !     define dimensions
    iret = nf_def_dim(ncid, 'm', m, m_dim)
    call check_err(iret)
    iret = nf_def_dim(ncid, 'n', n, n_dim)
    call check_err(iret)
    !     time is an unlimited dimension so that any number of
    !     records can be added
    iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
    call check_err(iret)
    !     define coordinate variable for time      
    t_dims(1) = time_dim
    iret = nf_def_var(ncid, 'time', NF_INT, 1, t_dims, t_id)
    call check_err(iret)
    !     define variables
    p_dims(1) = m_dim
    p_dims(2) = n_dim
    p_dims(3) = time_dim
    iret = nf_def_var(ncid, 'p', NF_DOUBLE, p_rank, p_dims, p_id)
    call check_err(iret)
    u_dims(1) = m_dim
    u_dims(2) = n_dim
    u_dims(3) = time_dim
    iret = nf_def_var(ncid, 'u', NF_DOUBLE, u_rank, u_dims, u_id)
    call check_err(iret)
    v_dims(1) = m_dim
    v_dims(2) = n_dim
    v_dims(3) = time_dim
    iret = nf_def_var(ncid, 'v', NF_DOUBLE, v_rank, v_dims, v_id)
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
    iret = nf_enddef(ncid)
    call check_err(iret)
      
    !     end of netCDF definitions
  END SUBROUTINE netcdf_setup

  !===================================================

  SUBROUTINE check_err(iret)
    INTEGER, INTENT(in) :: iret

    if(iret .ne. NF_NOERR) then
       print *, nf_strerror(iret)
       stop
    endif
  END SUBROUTINE check_err

  !===================================================
     
  SUBROUTINE my_ncwrite(id,varid,istart,icount,var,m,n,t_id,t_val)
!     Input args: id, varid,istart,icount,var,m,n,t_id,t_val
!     Write a whole array out to the netCDF file
      integer id,varid,iret
      integer icount(3)
      integer istart(3)
      integer m,n
      real(kind=8) var (m,n)
      integer t_id,t_val
      integer t_start(1), t_count(1)

      iret = nf_put_vara_double(id,varid,istart,icount,var)
      call check_err(iret)
      
      t_start(1) = istart(3) 
      t_count(1) = 1
      iret = nf_put_vara_int(id,t_id,t_start,t_count,t_val)
      call check_err(iret)

      END SUBROUTINE my_ncwrite

END MODULE shallow_io
