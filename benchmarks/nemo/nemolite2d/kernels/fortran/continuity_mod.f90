module continuity_mod
  use kind_params_mod
  use kernel_mod
  use argument_mod
  use grid_mod
  use field_mod
  implicit none

  type, extends(kernel_type) :: continuity
     type(go_arg), dimension(9) :: meta_args =         &
           (/ go_arg(GO_WRITE, GO_CT, GO_POINTWISE),            & ! ssha
              go_arg(GO_READ,  GO_CT, GO_POINTWISE),            & ! sshn
              go_arg(GO_READ,  GO_CU, GO_STENCIL(000,110,000)), & ! sshn_u
              go_arg(GO_READ,  GO_CV, GO_STENCIL(000,010,010)), & ! sshn_v
              go_arg(GO_READ,  GO_CU, GO_STENCIL(000,110,000)), & ! hu
              go_arg(GO_READ,  GO_CV, GO_STENCIL(000,010,010)), & ! hv
              go_arg(GO_READ,  GO_CU, GO_STENCIL(000,110,000)), & ! un
              go_arg(GO_READ,  GO_CV, GO_STENCIL(000,010,010)), & ! vn
              go_arg(GO_READ,  GO_GRID_AREA_T)           &
           /)
     !> This kernel updates only internal points of the simulation
     !! domain
     integer :: ITERATES_OVER = GO_INTERNAL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = GO_OFFSET_NE

  contains
    procedure, nopass :: code => continuity_code
  end type continuity

contains

  !===================================================

  subroutine invoke_continuity(ssha, sshn, sshn_u, sshn_v, hu, hv, un, vn)
    implicit none
    type(r2d_field),     intent(inout) :: ssha
    type(r2d_field),     intent(in) :: sshn, sshn_u, sshn_v
    type(r2d_field),     intent(in) :: hu, hv, un, vn
    ! Locals
    integer :: ji, jj

    do jj = ssha%internal%ystart, ssha%internal%ystop, 1
       do ji = ssha%internal%xstart, ssha%internal%xstop, 1

          call continuity_code(ji, jj,                      &
                               ssha%data, sshn%data,        &
                               sshn_u%data, sshn_v%data,    &
                               hu%data, hv%data, un%data, vn%data, &
                               ssha%grid%area_t)
       end do
    end do

  end subroutine invoke_continuity

  !=================================================================

  subroutine invoke_continuity_arrays(nx, ny, M, N, rdt, ssha, &
                                      sshn_t, sshn_u, sshn_v, &
                                      hu, hv, un, vn, area_t)
    use kind_params_mod
    use dl_timer, only: timer_start, timer_stop, i_def64
    implicit none
    integer, intent(in) :: nx, ny, M, N
    real(go_wp), intent(in) :: rdt
    real(go_wp), intent(out) :: ssha(nx,ny)
    real(go_wp), intent(in)  :: sshn_u(nx,ny), sshn_v(nx,ny), sshn_t(nx,ny)
    real(go_wp), intent(in)  :: un(nx,ny), vn(nx,ny)
    real(go_wp), intent(in)  :: hu(nx,ny), hv(nx,ny), area_t(nx,ny)
    ! Locals
    integer :: jj, ji
    real(go_wp) :: rtmp1, rtmp2, rtmp3, rtmp4
    !> For timing
    integer, save :: idxt
    integer :: ic
    integer(i_def64) :: nrepeat
    integer, parameter :: ALIGNMENT = 4
!DIR$ ASSUME (MOD(NX,ALIGNMENT) .EQ. 0)
!DIR$ ASSUME (MOD(M,ALIGNMENT) .EQ. 0)
!DIR$ ASSUME_ALIGNED ssha:64, sshn_u:64, sshn_v:64, sshn_t:64
!DIR$ ASSUME_ALIGNED un:64, vn:64, hu:64, hv:64, area_t:64


    ! Runtime check
    if( mod(M, ALIGNMENT) .ne. 0 ) then
        write(*,*) "This PSy-layer is compiled expecting a ", ALIGNMENT, &
            "-wide alignment."
        write(*,*) "Use DL_ESM_ALIGNMENT environment variable to match", &
           " this requierement."
        stop
    endif

    !nrepeat = 100 
    nrepeat = 1
   
    call timer_start(idxt, label='Continuity', num_repeats=nrepeat)
    !call likwid_markerStartRegion('Continuity')
!DIR$ VECTOR ALIGNED
    !do ic = 1, nrepeat, 1
    do jj = 2, N, 1

      ! Explicit peel loop
      do ji = 2, ALIGNMENT
         rtmp1 = (sshn_u(ji  ,jj ) + hu(ji  ,jj  ))*un(ji  ,jj)
         rtmp2 = (sshn_u(ji-1,jj ) + hu(ji-1,jj  ))*un(ji-1,jj)
         rtmp3 = (sshn_v(ji ,jj )  + hv(ji  ,jj  ))*vn(ji ,jj)
         rtmp4 = (sshn_v(ji ,jj-1) + hv(ji  ,jj-1))*vn(ji,jj-1)
         ssha(ji,jj) = sshn_t(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * &
               rdt / area_t(ji,jj)
      end do
      
      do ji = ALIGNMENT+1, M, 1

         rtmp1 = (sshn_u(ji  ,jj ) + hu(ji  ,jj  ))*un(ji  ,jj)
         rtmp2 = (sshn_u(ji-1,jj ) + hu(ji-1,jj))*un(ji-1,jj)
         rtmp3 = (sshn_v(ji ,jj )  + hv(ji  ,jj  ))*vn(ji ,jj)
         rtmp4 = (sshn_v(ji ,jj-1) + hv(ji ,jj-1))*vn(ji,jj-1)
         ssha(ji,jj) = sshn_t(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * &
                       rdt / area_t(ji,jj)
      end do
      
    end do
    !end do
    !call likwid_markerStopRegion('Continuity')

    call timer_stop(idxt)

  end subroutine invoke_continuity_arrays

  !=================================================================

  subroutine invoke_continuity_arrays_basic(nx, ny, M, N, rdt, ssha, &
                                            sshn_t, sshn_u, sshn_v, &
                                            hu, hv, un, vn, area_t)
    use kind_params_mod
    use dl_timer, only: timer_start, timer_stop, i_def64
    implicit none
    integer, intent(in) :: nx, ny, M, N
    real(go_wp), intent(in) :: rdt
    real(go_wp), intent(out) :: ssha(nx,ny)
    real(go_wp), intent(in)  :: sshn_u(nx,ny), sshn_v(nx,ny), sshn_t(nx,ny)
    real(go_wp), intent(in)  :: un(nx,ny), vn(nx,ny)
    real(go_wp), intent(in)  :: hu(nx,ny), hv(nx,ny), area_t(nx,ny)
    ! Locals
    integer :: jj, ji
    real(go_wp) :: rtmp1, rtmp2, rtmp3, rtmp4
    !> For timing
    integer, save :: idxt
    integer(i_def64) :: nrepeat

    !nrepeat = 100 
    nrepeat = 1
   
    call timer_start(idxt, label='Continuity', num_repeats=nrepeat)

    do jj = 2, N, 1
      do ji = 2, M, 2

         rtmp1 = (sshn_u(ji  ,jj ) + hu(ji  ,jj  ))*un(ji  ,jj)
         rtmp2 = (sshn_u(ji-1,jj ) + hu(ji-1,jj))*un(ji-1,jj)
         rtmp3 = (sshn_v(ji ,jj )  + hv(ji  ,jj  ))*vn(ji ,jj)
         rtmp4 = (sshn_v(ji ,jj-1) + hv(ji ,jj-1))*vn(ji,jj-1)
         ssha(ji,jj) = sshn_t(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * &
                       rdt / area_t(ji,jj)

         rtmp1 = (sshn_u(ji+1  ,jj ) + hu(ji+1  ,jj  ))*un(ji+1  ,jj)
         rtmp2 = (sshn_u(ji,jj ) + hu(ji,jj))*un(ji,jj)
         rtmp3 = (sshn_v(ji+1 ,jj )  + hv(ji+1  ,jj  ))*vn(ji+1 ,jj)
         rtmp4 = (sshn_v(ji+1 ,jj-1) + hv(ji+1 ,jj-1))*vn(ji+1,jj-1)
         ssha(ji+1,jj) = sshn_t(ji+1,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * &
                       rdt / area_t(ji+1,jj)
      end do
      
    end do

    call timer_stop(idxt)

  end subroutine invoke_continuity_arrays_basic

  !===================================================

  subroutine continuity_code(ji, jj,                     &
                             ssha, sshn, sshn_u, sshn_v, &
                             hu, hv, un, vn, e12t)
    use model_mod, only: rdt
    implicit none
    integer,                  intent(in)  :: ji, jj
    real(go_wp), dimension(:,:), intent(in)  :: e12t
    real(go_wp), dimension(:,:), intent(out) :: ssha
    real(go_wp), dimension(:,:), intent(in)  :: sshn, sshn_u, sshn_v, &
                                             hu, hv, un, vn
    ! Locals
    real(go_wp) :: rtmp1, rtmp2, rtmp3, rtmp4

    rtmp1 = (sshn_u(ji  ,jj  ) + hu(ji  ,jj  )) * un(ji  ,jj  )
    rtmp2 = (sshn_u(ji-1,jj  ) + hu(ji-1,jj  )) * un(ji-1,jj  )
    rtmp3 = (sshn_v(ji  ,jj  ) + hv(ji  ,jj  )) * vn(ji  ,jj  )
    rtmp4 = (sshn_v(ji  ,jj-1) + hv(ji  ,jj-1)) * vn(ji  ,jj-1)

    ssha(ji,jj) = sshn(ji,jj) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) * &
                    rdt / e12t(ji,jj)

  end subroutine continuity_code

end module continuity_mod
