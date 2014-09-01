module time_update_mod
  use kind_params_mod
  use field_mod
  implicit none

contains

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
