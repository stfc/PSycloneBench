module boundary_conditions_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use physical_params_mod
  use grid_mod
  use field_mod
  implicit none

  private

  public invoke_bc_solid_u,   invoke_bc_solid_v
  public invoke_bc_flather_u, invoke_bc_flather_v
  public invoke_bc_ssh

  !=======================================

  type, extends(kernel_type) :: bc_ssh
     type(arg), dimension(2) :: meta_args =  &
          (/ arg(READ,       R, POINTWISE),  &
             arg(READWRITE, CT, POINTWISE)   &
           /)

     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS

     !> This kernel is written assuming that the arrays for
     !! each field type are set-up such that the internal
     !! region of each field starts at the same array index (for
     !! both dimensions). If this weren't the case then
     !! these shifts (which are relative to the indexing used
     !! for fields on T points) would be non-zero.
     integer :: u_index_shift(2) = (/ 0, 0 /)
     integer :: v_index_shift(2) = (/ 0, 0 /)
     integer :: f_index_shift(2) = (/ 0, 0 /)

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => bc_ssh_code
  end type bc_ssh

  !=======================================

  type, extends(kernel_type) :: bc_solid_u
     type(arg), dimension(1) :: meta_args =  &
          (/ arg(WRITE, CU, POINTWISE)   &
           /)

     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS

     !> This kernel is written assuming that the arrays for
     !! each field type are set-up such that the internal
     !! region of each field starts at the same array index (for
     !! both dimensions). If this weren't the case then
     !! these shifts (which are relative to the indexing used
     !! for fields on T points) would be non-zero.
     integer :: u_index_shift(2) = (/ 0, 0 /)
     integer :: v_index_shift(2) = (/ 0, 0 /)
     integer :: f_index_shift(2) = (/ 0, 0 /)

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => bc_solid_u_code
  end type bc_solid_u

  !=======================================

  type, extends(kernel_type) :: bc_solid_v
     type(arg), dimension(1) :: meta_args =  &
          (/ arg(WRITE, CV, POINTWISE)   &
           /)

     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS

     !> This kernel is written assuming that the arrays for
     !! each field type are set-up such that the internal
     !! region of each field starts at the same array index (for
     !! both dimensions). If this weren't the case then
     !! these shifts (which are relative to the indexing used
     !! for fields on T points) would be non-zero.
     integer :: u_index_shift(2) = (/ 0, 0 /)
     integer :: v_index_shift(2) = (/ 0, 0 /)
     integer :: f_index_shift(2) = (/ 0, 0 /)

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => bc_solid_v_code
  end type bc_solid_v

  !=======================================

  type, extends(kernel_type) :: bc_flather_u
     type(arg), dimension(3) :: meta_args =  &
          (/ arg(READWRITE, CU, POINTWISE),  & ! ua
             arg(READ,      CU, POINTWISE),  & ! hu
             arg(READ,      CU, POINTWISE)   & ! sshn_u
           /)

     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS

     !> This kernel is written assuming that the arrays for
     !! each field type are set-up such that the internal
     !! region of each field starts at the same array index (for
     !! both dimensions). If this weren't the case then
     !! these shifts (which are relative to the indexing used
     !! for fields on T points) would be non-zero.
     integer :: u_index_shift(2) = (/ 0, 0 /)
     integer :: v_index_shift(2) = (/ 0, 0 /)
     integer :: f_index_shift(2) = (/ 0, 0 /)

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => bc_flather_u_code
  end type bc_flather_u

  !=======================================

  type, extends(kernel_type) :: bc_flather_v
     type(arg), dimension(3) :: meta_args =  &
          (/ arg(READWRITE, CV, POINTWISE),  & ! va
             arg(READ,      CV, POINTWISE),  & ! hv
             arg(READ,      CV, POINTWISE)   & ! sshn_v
           /)

     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS

     !> This kernel is written assuming that the arrays for
     !! each field type are set-up such that the internal
     !! region of each field starts at the same array index (for
     !! both dimensions). If this weren't the case then
     !! these shifts (which are relative to the indexing used
     !! for fields on T points) would be non-zero.
     integer :: u_index_shift(2) = (/ 0, 0 /)
     integer :: v_index_shift(2) = (/ 0, 0 /)
     integer :: f_index_shift(2) = (/ 0, 0 /)

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = OFFSET_NE

  contains
    procedure, nopass :: code => bc_flather_v_code
  end type bc_flather_v

contains
  
  !================================================

  subroutine invoke_bc_ssh(rtime, ssha)
    implicit none
    real(wp),             intent(in)    :: rtime
    type(r2d_field_type), intent(inout) :: ssha
    ! Locals
    real(wp) :: amp_tide, omega_tide
    integer  :: ji, jj

    amp_tide   = 0.2_wp
    omega_tide = 2.0_wp * 3.14159_wp / (12.42_wp * 3600._wp)

    !        DO jj = 1, jpj  
    !          DO ji = 1, jpi 

    DO jj = ssha%internal%ystart, ssha%internal%ystop
       DO ji = ssha%internal%xstart, ssha%internal%xstop
          call bc_ssh_code(ji, jj, ssha%grid%tmask, &
                           rtime, ssha%data)
       END DO
    END DO

  end subroutine invoke_bc_ssh
  
  !================================================

  subroutine bc_ssh_code(ji, jj, tmask, rtime, ssha)
    implicit none
    integer, intent(in)  :: ji, jj
    integer, dimension(:,:),  intent(in)    :: tmask
    real(wp),                 intent(in)    :: rtime
    real(wp), dimension(:,:), intent(inout) :: ssha
    ! Locals
    real(wp) :: amp_tide, omega_tide

    amp_tide   = 0.2_wp
    omega_tide = 2.0_wp * 3.14159_wp / (12.42_wp * 3600._wp)

    if(tmask(ji,jj) <= 0) return
    IF     (tmask(ji,jj-1) < 0) THEN
       ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
    ELSE IF(tmask(ji,jj+1) < 0) THEN
       ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
    ELSE IF(tmask(ji+1,jj) < 0) THEN
       ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
    ELSE IF(tmask(ji-1,jj) < 0) THEN
       ssha(ji,jj) = amp_tide * sin(omega_tide * rtime)
    END IF

  end subroutine bc_ssh_code
  
  !================================================

  !> Manual version of code to invoke kernel that applies solid 
  !! boundary conditions for u-velocity
  subroutine invoke_bc_solid_u(ua)
    implicit none
    type(r2d_field_type), intent(inout) :: ua
    ! Locals
    integer  :: ji, jj

! Original loop was:
!            DO jj = 1, jpj
!              DO ji = 0, jpi
! In original code, tmask is declared with one more row and column than
! any other field. ji==jpi IS last column of u field.
! 1/ How do I determine the full range of array indices to loop over for ua?
! 2/ If I do that, is tmask(ji+1,jj) going to stay within bounds?
    do jj = ua%whole%ystart, ua%whole%ystop, 1
       do ji = ua%whole%xstart, ua%whole%xstop, 1
          call bc_solid_u_code(ji, jj, ua%grid%tmask, ua%data)
       end do
    end do

  end subroutine invoke_bc_solid_u

  !================================================
  
  !> Kernel to apply solid boundary conditions for u-velocity
  subroutine bc_solid_u_code(ji, jj, tmask, ua)
    implicit none
    integer,                  intent(in)    :: ji, jj
    integer,  dimension(:,:), intent(in)    :: tmask
    real(wp), dimension(:,:), intent(inout) :: ua

    if(tmask(ji,jj) * tmask(ji+1,jj) == 0)then
       ua(ji,jj) = 0._wp
    end if

  end subroutine bc_solid_u_code
  
  !================================================

  !> Manual version of code to invoke the kernel
  !! that applies the solid-bc to a field on V pts.
  subroutine invoke_bc_solid_v(va)
    implicit none
    type(r2d_field_type), intent(inout) :: va
    ! Locals
    integer  :: ji, jj

    do jj = va%whole%ystart, va%whole%ystop, 1
       do ji = va%whole%xstart, va%whole%xstop, 1
          call bc_solid_v_code(ji,jj,va%grid%tmask,va%data)
      end do
    end do

  end subroutine invoke_bc_solid_v
  
  !================================================

  !> Kernel to apply solid boundary conditions for v-velocity
  subroutine bc_solid_v_code(ji, jj, tmask, va)
    implicit none
    integer,                 intent(in)    :: ji, jj
    integer, dimension(:,:), intent(in)    :: tmask
    real(wp),dimension(:,:), intent(inout) :: va

    if(tmask(ji,jj) * tmask(ji,jj+1) == 0)then
       va(ji,jj) = 0._wp
    end if

  end subroutine bc_solid_v_code
  
  !================================================

  !>                                  Du                 Dssh
  !!Flather open boundary condition [---- = sqrt(g/H) * ------]
  !!                                  Dn                 Dn
  !! ua and va in du/dn should be the specified tidal forcing
  subroutine invoke_bc_flather_u(ua, hu, sshn_u)
    implicit none
    type(r2d_field_type), intent(inout) :: ua
    type(r2d_field_type), intent(in) :: hu, sshn_u
    ! Locals
    integer  :: ji, jj

    ! Original loop was:
    !            DO jj = 1, jpj
    !              DO ji = 0, jpi  
    DO jj = ua%whole%ystart, ua%whole%ystop, 1
       DO ji = ua%whole%xstart, ua%whole%xstop, 1
          call bc_flather_u_code(ji,jj,ua%grid%tmask, &
                                 ua%data, hu%data, sshn_u%data)
       END DO
    END DO
  
  end subroutine invoke_bc_flather_u
  
  !================================================

  !> Kernel to apply Flather condition to U
  subroutine bc_flather_u_code(ji, jj, tmask, ua, hu, sshn_u)
    implicit none
    integer,                  intent(in)    :: ji, jj
    integer,  dimension(:,:), intent(in)    :: tmask
    real(wp), dimension(:,:), intent(inout) :: ua
    real(wp), dimension(:,:), intent(in)    :: hu, sshn_u
    ! Locals
    integer  :: jiu

    !                                  Du                 Dssh
    !Flather open boundary condition [---- = sqrt(g/H) * ------]
    !                                  Dn                 Dn
    ! ua and va in du/dn should be the specified tidal forcing

    ! Check whether this point lies within the domain
    if(tmask(ji,jj) + tmask(ji+1,jj) <= -1) return

    if(tmask(ji,jj) < 0) then
       jiu = ji + 1
       ua(ji,jj) = ua(jiu,jj) + &
                   sqrt(g/hu(ji,jj)) * (sshn_u(ji,jj) - sshn_u(jiu,jj))
    else if(tmask(ji+1,jj )< 0) then
       jiu = ji - 1 
       ua(ji,jj) = ua(jiu,jj) + sqrt(g/hu(ji,jj)) * &
            (sshn_u(ji,jj) - sshn_u(jiu,jj))
    end if
  
  end subroutine bc_flather_u_code

  !================================================

  !> Manual version of code to invoke the kernel for applying the
  !! Flather boundary condition to the v component of velocity.
  subroutine invoke_bc_flather_v(va, hv, sshn_v)
    implicit none
    type(r2d_field_type), intent(inout) :: va
    type(r2d_field_type), intent(in)    :: hv, sshn_v
    ! Locals
    integer  :: ji, jj

    !kernel Flather v 

    DO jj = va%whole%ystart, va%whole%ystop, 1
       DO ji = va%whole%xstart, va%whole%xstop, 1
          call bc_flather_v_code(ji,jj,va%grid%tmask, &
                                 va%data, hv%data, sshn_v%data)
       END DO
    END DO

  end subroutine invoke_bc_flather_v

  !================================================

  !> Kernel to apply Flather boundary condition to v component
  !! of velocity
  subroutine bc_flather_v_code(ji, jj, tmask, va, hv, sshn_v)
    implicit none
    integer,                  intent(in) :: ji, jj
    integer,  dimension(:,:), intent(in) :: tmask
    real(wp), dimension(:,:), intent(inout) :: va
    real(wp), dimension(:,:), intent(in) :: hv, sshn_v
    ! Locals
    integer  :: jiv

    ! Check whether this point is inside the simulated domain
    IF(tmask(ji,jj) + tmask(ji,jj+1) <= -1) return
    
    IF(tmask(ji,jj) < 0) THEN
       jiv = jj + 1
       va(ji,jj) = va(ji,jiv) + SQRT(g/hv(ji,jj)) * &
                                    (sshn_v(ji,jj) - sshn_v(ji,jiv))
    ELSE IF(tmask(ji,jj+1) < 0) THEN
       jiv = jj - 1 
       va(ji,jj) = va(ji,jiv) + SQRT(g/hv(ji,jj)) * &
                                    (sshn_v(ji,jj) - sshn_v(ji,jiv))
    END IF

  end subroutine bc_flather_v_code

  !================================================

end module boundary_conditions_mod
