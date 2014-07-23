module continuity_mod
  implicit none

  type, extends(kernel_type) :: continuity
     type(arg), dimension(3) :: meta_args =    &
          (/ arg(WRITE, CT, POINTWISE),        & ! ssha
             arg(READ,  CT, POINTWISE),        & ! 
             arg(READ,  CV, POINTWISE)         & ! 
           /)
     !> We only have one value per grid point and that means
     !! we have a single DOF per grid point.
     integer :: ITERATES_OVER = DOFS
  contains
    procedure, nopass :: code => continuity_code
  end type continuity

contains

  !===================================================

  subroutine invoke_continuity(ssha, sshn_u, sshn_v, hu, hv, un, vn)
    implicit none

    do jj = 1, jpj
       do ji = 1, jpi
          call continuity_code(ji, jj, ssha, sshn_u, sshn_v, hu, hv, un, vn)
       end do
    end do

  end subroutine invoke_continuity

  !===================================================

  subroutine continuity_code(ji, jj, ssha, sshn_u, sshn_v, hu, hv, un, vn)
    implicit none
    integer,                  intent(in)  :: ji, jj
    real(wp), dimension(:,:), intent(out) :: ssha
    real(wp), dimension(:,:), intent(in)  :: sshn_u, sshn_v, hu, hv, un, vn
    ! Locals
    real(wp) :: rtmp1, rtmp2, rtmp3, rtmp4

    rtmp1 = (sshn_u(ji  ,jj  ) + hu(ji  ,jj  )) * un(ji  ,jj  )
    rtmp2 = (sshn_u(ji-1,jj  ) + hu(ji-1,jj  )) * un(ji-1,jj  )
    rtmp3 = (sshn_v(ji  ,jj  ) + hv(ji  ,jj  )) * vn(ji  ,jj  )
    rtmp4 = (sshn_v(ji  ,jj-1) + hv(ji  ,jj-1)) * vn(ji  ,jj-1)

    ssha(ji,jj) = sshn(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * rdt / e12t(ji,jj)

  end subroutine continuity_code

end module continuity_mod
