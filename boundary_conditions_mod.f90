module boundary_conditions_mod
  use kind_params_mod
  use physical_params_mod
  use field_mod
  implicit none

  private

  public bc_u_solid, bc_v_solid
  public bc_u_flather, bc_v_flather
  public bc_ssh

contains
  
  !================================================

  subroutine bc_ssh(rtime, ssha)
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
          IF(ssha%grid%tmask(ji,jj) <= 0) CYCLE
          IF     (ssha%grid%tmask(ji,jj-1) < 0) THEN
             ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
          ELSE IF(ssha%grid%tmask(ji,jj+1) < 0) THEN
             ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
          ELSE IF(ssha%grid%tmask(ji+1,jj) < 0) THEN
             ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
          ELSE IF(ssha%grid%tmask(ji-1,jj) < 0) THEN
             ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
          END IF
       END DO
    END DO

  end subroutine bc_ssh
  
  !================================================

  subroutine bc_u_solid(ua)
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

          if(ua%grid%tmask(ji,jj) * ua%grid%tmask(ji+1,jj) == 0)then
            ua%data(ji,jj) = 0._wp
          end if
       end do
    end do

  end subroutine bc_u_solid
  
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
