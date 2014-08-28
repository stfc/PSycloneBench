module time_update_mod
  use kind_params_mod
  use field_mod
  implicit none

contains

  !+++++++++++++++++++++++++++++++++++

  subroutine next(grid)
    use kind_params_mod
    use grid_mod
    use model_mod, only: jpi, jpj, ua, va, un, vn, ssha, sshn
    use model_mod, only: sshn_v, sshn_u
    implicit none
    type(grid_type), intent(in) :: grid
    ! Locals
    integer :: ji, jj
    real(wp) :: rtmp1

    ! update the now-velocity and ssh

  
    ! kernel  un updating
    DO jj = 1, jpj
       DO ji = 0, jpi
          un(ji,jj)   = ua(ji,jj)
       END DO
    END DO
    ! end of kernel sshn_u updating.

    ! kernel vn updating
    DO jj = 0, jpj
       DO ji = 1, jpi
          vn(ji,jj)   = va(ji,jj)
       END DO
    END DO
    ! end kernel vn updating.

    ! kernel sshn updating
    DO jj = 1, jpj
       DO ji = 1, jpi
          sshn(ji,jj) = ssha(ji,jj)
       END DO
    END DO
    ! end kernel sshn_u updating.

    ! kernel sshn_u updating
    DO jj = 1, jpj
       DO ji = 0, jpi
          IF(grid%tmask(ji,jj) + grid%tmask(ji+1,jj) <= 0)  CYCLE                              !jump over non-computational domain
          IF(grid%tmask(ji,jj) * grid%tmask(ji+1,jj) > 0) THEN
             rtmp1 = grid%e12t(ji,jj) * sshn(ji,jj) + grid%e12t(ji+1,jj) * sshn(ji+1,jj)
             sshn_u(ji,jj) = 0.5_wp * rtmp1 / grid%e12u(ji,jj) 
          ELSE IF(grid%tmask(ji,jj) <= 0) THEN
             sshn_u(ji,jj) = sshn(ji+1,jj)
          ELSE IF(grid%tmask(ji+1,jj) <= 0) THEN
             sshn_u(ji,jj) = sshn(ji,jj)
          END IF
       END DO
    END DO
    ! end kernel sshn_u updating.

    ! kernel: sshn_v updating
    DO jj = 0, jpj
       DO ji = 1, jpi
          IF(grid%tmask(ji,jj) + grid%tmask(ji,jj+1) <= 0)  CYCLE                              !jump over non-computational domain
          IF(grid%tmask(ji,jj) * grid%tmask(ji,jj+1) > 0) THEN
             rtmp1 = grid%e12t(ji,jj) * sshn(ji,jj) + grid%e12t(ji,jj+1) * sshn(ji,jj+1)
             sshn_v(ji,jj) = 0.5_wp * rtmp1 / grid%e12v(ji,jj) 
          ELSE IF(grid%tmask(ji,jj) <= 0) THEN
             sshn_v(ji,jj) = sshn(ji,jj+1)
          ELSE IF(grid%tmask(ji,jj+1) <= 0) THEN
             sshn_v(ji,jj) = sshn(ji,jj)
          END If
       END DO
    END DO
    ! end kernel sshn_v updating.
            
  END SUBROUTINE next

  !================================================

  subroutine invoke_next_u(un, ua)
    implicit none
    type(r2d_field_type), intent(inout) :: un
    type(r2d_field_type), intent(in)    :: ua
    ! Locals
    integer :: ji, jj

    do jj = un%internal%ystart, un%internal%ystop, 1
      do ji = un%internal%xstart, un%internal%xstop, 1
        call next_u_code(ji,jj, un%data, ua%data)
      end do
    end do

  contains

    subroutine next_u_code(ji,jj,un,ua)
      implicit none
      integer, intent(in) :: ji, jj
      real(wp), dimension(:,:), intent(inout) :: un
      real(wp), dimension(:,:), intent(in)    :: ua

      un(ji,jj)   = ua(ji,jj)

    end subroutine next_u_code

  end subroutine invoke_next_u

  !================================================

  subroutine invoke_next_v(vn, va)
    implicit none
    type(r2d_field_type), intent(inout) :: vn
    type(r2d_field_type), intent(in)    :: va
    ! Locals
    integer :: ji, jj

    do jj = vn%internal%ystart, vn%internal%ystop, 1
      do ji = vn%internal%xstart, vn%internal%xstop, 1
        call next_v_code(ji,jj, vn%data, va%data)
      end do
    end do

  contains

    subroutine next_v_code(ji,jj,vn,va)
      implicit none
      integer, intent(in) :: ji, jj
      real(wp), dimension(:,:), intent(inout) :: vn
      real(wp), dimension(:,:), intent(in)    :: va

      vn(ji,jj)   = va(ji,jj)
    end subroutine next_v_code

  end subroutine invoke_next_v

  !================================================

  subroutine invoke_next_ssht(sshn, ssha)
    implicit none
    type(r2d_field_type), intent(inout) :: sshn
    type(r2d_field_type), intent(in)    :: ssha
    ! Locals
    integer :: ji, jj

    do jj = sshn%internal%ystart, sshn%internal%ystop, 1
      do ji = sshn%internal%xstart, sshn%internal%xstop, 1
        call next_ssht_code(ji,jj,sshn%data, ssha%data)
      end do
    end do

  contains

    subroutine next_ssht_code(ji, jj, sshn, ssha)
      implicit none
      integer, intent(in) :: ji, jj
      real(wp), dimension(:,:), intent(inout) :: sshn
      real(wp), dimension(:,:), intent(in)    :: ssha

      sshn(ji,jj) = ssha(ji,jj)

    end subroutine next_ssht_code

  end subroutine invoke_next_ssht

  !================================================

  subroutine invoke_next_sshu(sshn_u, sshn)
    implicit none
    type(r2d_field_type), intent(inout) :: sshn_u
    type(r2d_field_type), intent(in)    :: sshn
    ! Locals
    integer :: ji, jj

    do jj = sshn_u%internal%ystart, sshn_u%internal%ystop, 1
      do ji = sshn_u%internal%xstart, sshn_u%internal%xstop, 1

        call next_sshu_code(ji, jj, sshn_u%grid%tmask, &
                            sshn_u%grid%e12t, sshn_u%grid%e12u, &
                            sshn_u%data, sshn%data)
      end do
    end do

  contains

    subroutine next_sshu_code(ji,jj,tmask,e12t,e12u, &
                              sshn_u, sshn)
      implicit none
      integer,                  intent(in)    :: ji, jj
      integer,  dimension(:,:), intent(in)    :: tmask
      real(wp), dimension(:,:), intent(in)    :: e12t, e12u
      real(wp), dimension(:,:), intent(inout) :: sshn_u
      real(wp), dimension(:,:), intent(in)    :: sshn
      ! Locals
      real(wp) :: rtmp1

      if(tmask(ji,jj) + tmask(ji+1,jj) <= 0)  return   !jump over non-computational domain

      IF(tmask(ji,jj) * tmask(ji+1,jj) > 0) THEN
        rtmp1 = e12t(ji,jj) * sshn(ji,jj) + e12t(ji+1,jj) * sshn(ji+1,jj)
        sshn_u(ji,jj) = 0.5_wp * rtmp1 / e12u(ji,jj) 
      ELSE IF(tmask(ji,jj) <= 0) THEN
         sshn_u(ji,jj) = sshn(ji+1,jj)
      ELSE IF(tmask(ji+1,jj) <= 0) THEN
         sshn_u(ji,jj) = sshn(ji,jj)
      END IF
    end subroutine next_sshu_code

  end subroutine invoke_next_sshu

  !================================================

  subroutine invoke_next_sshv(sshn_v, sshn)
    implicit none
    type(r2d_field_type), intent(inout) :: sshn_v
    type(r2d_field_type), intent(in)    :: sshn
    ! Locals
    integer :: ji, jj

    do jj = sshn_v%internal%ystart, sshn_v%internal%ystop, 1
      do ji = sshn_v%internal%xstart, sshn_v%internal%xstop, 1

        call next_sshv_code(ji, jj, sshn_v%grid%tmask, &
                            sshn_v%grid%e12t, sshn_v%grid%e12v, &
                            sshn_v%data, sshn%data )
      end do
    end do

  contains

    subroutine next_sshv_code(ji, jj, tmask, e12t, e12v, &
                              sshn_v, sshn)
      implicit none
      integer,                  intent(in)    :: ji, jj
      integer,  dimension(:,:), intent(in)    :: tmask
      real(wp), dimension(:,:), intent(in)    :: e12t, e12v
      real(wp), dimension(:,:), intent(inout) :: sshn_v
      real(wp), dimension(:,:), intent(in)    :: sshn
      ! Locals
      real(wp) :: rtmp1

      IF(tmask(ji,jj) + tmask(ji,jj+1) <= 0)  return !jump over non-computational domain
      IF(tmask(ji,jj) * tmask(ji,jj+1) > 0) THEN
         rtmp1 = e12t(ji,jj) * sshn(ji,jj) + e12t(ji,jj+1) * sshn(ji,jj+1)
         sshn_v(ji,jj) = 0.5_wp * rtmp1 / e12v(ji,jj) 
      ELSE IF(tmask(ji,jj) <= 0) THEN
         sshn_v(ji,jj) = sshn(ji,jj+1)
      ELSE IF(tmask(ji,jj+1) <= 0) THEN
         sshn_v(ji,jj) = sshn(ji,jj)
      END If

    end subroutine next_sshv_code

  end subroutine invoke_next_sshv

  !================================================

end module time_update_mod
