module time_update_mod
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

  subroutine invoke_next_un
    implicit none

  contains

    subroutine next_un_code(ji,jj)
      implicit none
      integer, intent(in) :: ji, jj
      DO jj = 1, jpj
         DO ji = 0, jpi
            un(ji,jj)   = ua(ji,jj)
         END DO
      END DO
    end subroutine next_un_code

  end subroutine invoke_next_un

  !================================================

  subroutine invoke_next_vn
    implicit none

    DO jj = 0, jpj
       DO ji = 1, jpi
          vn(ji,jj)   = va(ji,jj)
       END DO
    END DO

  end subroutine invoke_next_vn

  !================================================

  subroutine invoke_next_sshn_t

    DO jj = 1, jpj
       DO ji = 1, jpi
          sshn(ji,jj) = ssha(ji,jj)
       END DO
    END DO

  end subroutine invoke_next_sshn_t

  !================================================

  subroutine invoke_next_sshn_u

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

  end subroutine invoke_next_sshn_u

  !================================================

  subroutine invoke_next_sshn_v(sshn_fld, sshn_v_fld)
    implicit none

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

  end subroutine invoke_next_sshn_v

  !================================================

end module time_update_mod
