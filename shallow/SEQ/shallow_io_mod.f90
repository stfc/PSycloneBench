module shallow_io_mod
  use kind_params_mod
  use field_mod
  implicit none

  private

  logical :: l_out   !< Whether or not to produce output  
  integer :: mprint  !< frequency of output    

  ! NetCDF variables
  include 'netcdf.inc'

  integer :: ncid, t_id, p_id, u_id, v_id, iret, t_val
  integer, dimension(3) :: istart, icount 
  character (len=13) :: ncfile = "shallowdat.nc"

  public read_namelist, print_initial_values, print_diagonals
  public model_write_init, model_write, model_write_finalise
  public ascii_write

contains

  !===================================================

  !> Reads the namelist file for user-specified control 
  !! parameters
  SUBROUTINE read_namelist(m, n, itmax)
    IMPLICIT none
    INTEGER, INTENT(out) :: m, n
    INTEGER, INTENT(out) :: itmax
    ! namelist input 
    CHARACTER (LEN=8) :: nml_name = "namelist" 
    INTEGER :: input_unit = 99
    INTEGER :: ierr

    NAMELIST/global_domain/ m, n, itmax, mprint
    NAMELIST/io_control/ l_out

    ! Initialise these vars to problem values as defense
    ! against failure to read them properly and also
    ! to squelch incorrect compiler warnings about them
    ! not being assigned to.
    m = -1
    n = -1
    itmax = 0

    !     Read in namelist 
    OPEN(unit=input_unit, file=nml_name, status='old',iostat=ierr)

    CALL check(ierr, "open "//nml_name)
    READ(unit=input_unit, nml=global_domain, iostat=ierr)
    CALL check(ierr, "read "//nml_name)
    READ(unit=input_unit, nml=io_control, iostat=ierr)
    CALL check(ierr, "read "//nml_name)

    CLOSE(unit=input_unit)

  END SUBROUTINE read_namelist

  !===================================================

  !> Log initial parameter values
  SUBROUTINE print_initial_values(m,n,dx,dy,dt,alpha)
    IMPLICIT none
    INTEGER,      INTENT(in) :: m, n
    REAL(KIND=8), INTENT(in) :: dx, dy, dt, alpha

    IF(l_out)THEN

       WRITE(6,390) N,M,DX,DY,DT,ALPHA
 390     FORMAT(" NUMBER OF POINTS IN THE X DIRECTION",I8,/    & 
                " NUMBER OF POINTS IN THE Y DIRECTION",I8,/    & 
                " GRID SPACING IN THE X DIRECTION    ",F8.0,/  & 
                " GRID SPACING IN THE Y DIRECTION    ",F8.0,/  & 
                " TIME STEP                          ",F8.0,/  & 
                " TIME FILTER PARAMETER              ",F8.3)
    END IF

  END SUBROUTINE print_initial_values

  !===================================================

  !> Print diagonal elements of solution arrays - diagnostic
  SUBROUTINE print_diagonals(p, u, v)
    IMPLICIT none
    REAL(KIND=8), INTENT(in), DIMENSION(:,:) :: p, u, v
    ! Locals
    INTEGER :: i, mnmin, m, n

    m = SIZE(p,1)-1
    n = SIZE(p,2)-1

    MNMIN = MIN0(m,n)

    WRITE(6,391) (P(I,I),I=1,MNMIN)
391 FORMAT(/' INITIAL DIAGONAL ELEMENTS OF P ' //,(8E15.6))
    WRITE(6,392) (U(I,I),I=1,MNMIN)
392 FORMAT(/' INITIAL DIAGONAL ELEMENTS OF U ' //,(8E15.6))
    WRITE(6,393) (V(I,I),I=1,MNMIN)
393 FORMAT(/' INITIAL DIAGONAL ELEMENTS OF V ' //,(8E15.6))

  END SUBROUTINE print_diagonals

  !===================================================
  
  SUBROUTINE model_write_init(m, n)
    IMPLICIT none
    INTEGER, INTENT(in) :: m, n

    ! Prepare netCDF file to receive model output data
    IF (l_out) call netcdf_setup(ncfile,m,n,ncid,t_id,p_id,u_id,v_id, &
                                 istart,icount)

  END SUBROUTINE model_write_init
 
  !===================================================

  !> Write data for the current time step
  SUBROUTINE model_write(ncycle, pfld, ufld, vfld)
    IMPLICIT none
    integer,                      intent(in) :: ncycle
    type(r2d_field), intent(in), target :: pfld, ufld, vfld
    ! Locals
    real(kind=wp), dimension(:,:), pointer :: p, u, v
    INTEGER :: m, n

    u => ufld%data(ufld%internal%xstart:ufld%internal%xstop, &
                   ufld%internal%ystart:ufld%internal%ystop)

    v => vfld%data(vfld%internal%xstart:vfld%internal%xstop, &
                   vfld%internal%ystart:vfld%internal%ystop)

    p => pfld%data(pfld%internal%xstart:pfld%internal%xstop, &
                   pfld%internal%ystart:pfld%internal%ystop)


    IF( l_out .AND. (MOD(NCYCLE,MPRINT) .EQ. 0) ) then

       call ascii_write(ncycle, 'ufld.dat', ufld%data, &
                        ufld%internal%nx, ufld%internal%ny, &
                        ! We output data as though first internal
                        ! U pt is at (3,2). This is equivalent to
                        ! (2,1) on the mesh arrangement in the
                        ! original version of shallow. We only have
                        ! to do this so that the output can be
                        ! compared with the original in a straightforward
                        ! fashion.
                        3, 2 )!ufld%internal%xstart, ufld%internal%ystart)

       call ascii_write(ncycle, 'vfld.dat', vfld%data,      &
                        vfld%internal%nx, vfld%internal%ny, &
                        ! We output data as though first internal
                        ! V pt is at (2,3). This is equivalent to
                        ! (1,2) on the mesh arrangement in the
                        ! original version of shallow. We only have
                        ! to do this so that the output can be
                        ! compared with the original in a straightforward
                        ! fashion.
                        2, 3 )!vfld%internal%xstart, vfld%internal%ystart)

       call ascii_write(ncycle, 'pfld.dat', pfld%data,      &
                        pfld%internal%nx, pfld%internal%ny, &
                        ! We output data as though first internal
                        ! T pt is at (2,2) (which, unlike all the other
                        ! grid-pt types, it is). This is equivalent to
                        ! (1,1) on the mesh arrangement in the
                        ! original version of shallow. We only have
                        ! to do this so that the output can be
                        ! compared with the original in a straightforward
                        ! fashion.
                        2, 2 )!pfld%internal%xstart, pfld%internal%ystart)

       !CALL print_diagonals(p, u, v)
    
       m = ufld%internal%nx
       n = ufld%internal%ny

       !           Append calculated values of p, u, and v to netCDF file
       istart(3) = ncycle/mprint + 1
       t_val = ncycle

       !           Shape of record to be written (one ncycle at a time)
       call my_ncwrite(ncid,p_id,istart,icount, &
                       p, &
                       m,n,t_id,t_val)
       call my_ncwrite(ncid,u_id,istart,icount, &
                       u, &
                       m,n,t_id,t_val)
       call my_ncwrite(ncid,v_id,istart,icount, &
                       v, &
                       m,n,t_id,t_val)
    END IF

  END SUBROUTINE model_write

  !===================================================

  SUBROUTINE model_write_finalise()
    IMPLICIT none
    INTEGER :: iret
    
    ! Close the netCDF file

    IF (l_out) THEN 
       iret = nf_close(ncid)
       call check_err(iret)
    ENDIF

  END SUBROUTINE model_write_finalise

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

    !===================================================

    !> Dump the supplied 2D field in ASCII form suitable for
    !! gnuplot splot.
    !! @param[in] tag Tag to incorporate in filename - usually time-step
    !! @param[in] fname Base of filename that will have tag appended
    !! @param[in] var 2D double-precision array of data to write out
    !! @param[in] m Extent in x dimension of INTERNAL region of field
    !! @param[in] n Extent in y dimension of INTERNAL region of field
    !! @param[in] xstart Grid loc of start of INTERNAL region
    !! @param[in] ystart Grid loc of start of INTERNAL region
    subroutine ascii_write(tag, fname, var, m, n, xstart, ystart)
      implicit none
      integer,                 intent(in) :: tag, m, n
      integer, optional,       intent(in) :: xstart, ystart
      real(8), dimension(:,:), intent(in) :: var
      character(len=*),        intent(in) :: fname
      ! Locals
      integer :: ji, jj
      integer, parameter :: iounit = 24
      character(len=100) :: fname_full
      character(len=10) :: tag_str

      write(tag_str, fmt="(I10)") tag
      write(fname_full, fmt="((A),'.',(A))") TRIM(ADJUSTL(fname)), &
                                             TRIM(ADJUSTL(tag_str))
      open(unit=iounit, file=TRIM(fname_full), status='unknown', &
           action='write', iostat=ji)
      if( ji /= 0 )then
         write (*,*) 'ascii_write: ERROR: failed to open file'
         return
      end if

      if(present(xstart) .AND. present(ystart))then
         do jj=ystart, ystart+n-1, 1
            do ji=xstart, xstart+m-1, 1
               write(iounit,"(I5,1x,I5,1x,E24.16)") &
                    ji-xstart+1, jj-ystart+1, var(ji,jj)
            end do
            write(iounit,*)
         end do
      else
         do jj = 1, size(var, 2)
            do ji = 1, size(var, 1)
               write(iounit,"(I5,1x,I5,1x,E24.16)") &
                    ji, jj, var(ji,jj)
            end do
            write(iounit,*)
         end do
      end if

      close(unit=iounit)

    end subroutine ascii_write

END MODULE shallow_io_mod
