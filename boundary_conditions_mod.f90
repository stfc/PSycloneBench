module boundary_conditions_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use physical_params_mod
  use grid_mod
  use field_mod
  implicit none

  private

  public invoke_bc_solid_u, bc_v_solid
  public bc_u_flather, bc_v_flather
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

  subroutine invoke_bc_solid_u(ua)
    implicit none
    type(r2d_field_type), intent(inout) :: ua
    ! Locals
    integer  :: ji, jj

    ! solid boundary conditions for u-velocity

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

  subroutine bc_solid_u_code(ji, jj, tmask, ua)
    implicit none
    integer,                  intent(in)    :: ji, jj
    integer,  dimension(:,:), intent(in)    :: tmask
    real(wp), dimension(:,:), intent(inout) :: ua

    ! solid boundary conditions for u-velocity

    if(tmask(ji,jj) * tmask(ji+1,jj) == 0)then
       ua(ji,jj) = 0._wp
    end if

  end subroutine bc_solid_u_code
  
  !================================================

  subroutine bc_v_solid(va)
    implicit none
    type(r2d_field_type), intent(inout) :: va
    ! Locals
    integer  :: ji, jj

    ! solid boundary conditions for v-velocity

    do jj = va%whole%ystart, va%whole%ystop, 1
       do ji = va%whole%xstart, va%whole%xstop, 1

        if(va%grid%tmask(ji,jj) * va%grid%tmask(ji,jj+1) == 0)then
          va%data(ji,jj) = 0._wp
        end if
      end do
    end do

  end subroutine bc_v_solid
  
  !================================================

  subroutine bc_u_flather(ua, hu, sshn_u)
    implicit none
    type(r2d_field_type), intent(inout) :: ua
    type(r2d_field_type), intent(in) :: hu, sshn_u
    ! Locals
    integer  :: ji, jj, jiu

    !                                  Du                 Dssh
    !Flather open boundary condition [---- = sqrt(g/H) * ------]
    !                                  Dn                 Dn
    ! ua and va in du/dn should be the specified tidal forcing

    ! Original loop was:
    !            DO jj = 1, jpj
    !              DO ji = 0, jpi  

    ! Flather u 
    DO jj = ua%whole%ystart, ua%whole%ystop, 1
       DO ji = ua%whole%xstart, ua%whole%xstop, 1

          IF(ua%grid%tmask(ji,jj) + ua%grid%tmask(ji+1,jj) <= -1) CYCLE  ! not in the domain

          IF(ua%grid%tmask(ji,jj) < 0) THEN
             jiu = ji + 1
             ua%data(ji,jj) = ua%data(jiu,jj) + &
                              SQRT(g/hu%data(ji,jj)) * (sshn_u%data(ji,jj) - &
                              sshn_u%data(jiu,jj))
          ELSE IF(ua%grid%tmask(ji+1,jj )< 0) THEN
             jiu = ji - 1 
             ua%data(ji,jj) = ua%data(jiu,jj) + SQRT(g/hu%data(ji,jj)) * &
                                 (sshn_u%data(ji,jj) - sshn_u%data(jiu,jj))
          END IF
       END DO
    END DO
  
  end subroutine bc_u_flather

  !================================================

  subroutine bc_v_flather(va, hv, sshn_v)
    implicit none
    type(r2d_field_type), intent(inout) :: va
    type(r2d_field_type), intent(in) :: hv, sshn_v
    ! Locals
    integer  :: ji, jj, jiv

    !kernel Flather v 

    DO jj = va%whole%ystart, va%whole%ystop, 1
       DO ji = va%whole%xstart, va%whole%xstop, 1

          ! Check whether this point is inside the simulated domain
          IF(va%grid%tmask(ji,jj) + va%grid%tmask(ji,jj+1) <= -1) CYCLE

          IF(va%grid%tmask(ji,jj) < 0) THEN
             jiv = jj + 1
             va%data(ji,jj) = va%data(ji,jiv) + SQRT(g/hv%data(ji,jj)) * &
                                    (sshn_v%data(ji,jj) - sshn_v%data(ji,jiv))
          ELSE IF(va%grid%tmask(ji,jj+1) < 0) THEN
             jiv = jj - 1 
             va%data(ji,jj) = va%data(ji,jiv) + SQRT(g/hv%data(ji,jj)) * &
                                    (sshn_v%data(ji,jj) - sshn_v%data(ji,jiv))
          END IF
       END DO
    END DO

  end subroutine bc_v_flather

  !================================================

end module boundary_conditions_mod