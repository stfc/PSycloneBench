module boundary_conditions_mod
  use kind_params_mod, only: go_wp
  use argument_mod, only: GO_READ, GO_READWRITE, GO_I_SCALAR, &
       GO_ARG, GO_R_SCALAR, GO_CU, GO_CV, GO_CT, GO_GRID_MASK_T, &
       GO_STENCIL
  use kernel_mod, only: kernel_type, GO_POINTWISE, GO_DOFS, &
      GO_ALL_PTS, GO_INTERNAL_PTS
  use physical_params_mod
  use grid_mod
  use field_mod
  implicit none

  private

  public bc_ssh, bc_solid_u, bc_solid_v, bc_flather_u, bc_flather_v
  public invoke_bc_solid_u,   invoke_bc_solid_v
  public invoke_bc_flather_u, invoke_bc_flather_v
  public invoke_bc_ssh, setup_vmask_code
  public bc_ssh_code, bc_solid_u_code, bc_solid_v_code
  public bc_flather_u_code, bc_flather_v_code

  !=======================================

  type, extends(kernel_type) :: bc_ssh
     type(go_arg), dimension(3) :: meta_args =        &
          (/ go_arg(GO_READ,      GO_I_SCALAR, GO_POINTWISE),  &
             go_arg(GO_READWRITE, GO_CT,       GO_POINTWISE),  &
             go_arg(GO_READ,      GO_GRID_MASK_T)           &
           /)

     !> Although this is a boundary-conditions kernel, it only
     !! acts on the internal points of the domain
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
    procedure, nopass :: code => bc_ssh_code
  end type bc_ssh

  !=======================================

  type, extends(kernel_type) :: bc_solid_u
     type(go_arg), dimension(2) :: meta_args =  &
          (/ go_arg(GO_READWRITE, GO_CU, GO_POINTWISE),  &
             go_arg(GO_READ,      GO_GRID_MASK_T)     &
           /)

     !> This is a boundary-conditions kernel and therefore
     !! acts on all points of the domain rather than just
     !! those that are internal
     integer :: ITERATES_OVER = GO_ALL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = GO_OFFSET_NE

  contains
    procedure, nopass :: code => bc_solid_u_code
  end type bc_solid_u

  !=======================================

  type, extends(kernel_type) :: bc_solid_v
     type(go_arg), dimension(2) :: meta_args =  &
          (/ go_arg(GO_READWRITE, GO_CV, GO_POINTWISE),  &
             go_arg(GO_READ,      GO_GRID_MASK_T)        &
           /)

     !> This is a boundary-conditions kernel and therefore
     !! acts on all points of the domain rather than just
     !! those that are internal
     integer :: ITERATES_OVER = GO_ALL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = GO_OFFSET_NE

  contains
    procedure, nopass :: code => bc_solid_v_code
  end type bc_solid_v

  !=======================================

  type, extends(kernel_type) :: bc_flather_u
     type(go_arg), dimension(4) :: meta_args =  &
          (/ go_arg(GO_READWRITE, GO_CU, GO_STENCIL(000,111,000)),  & ! ua
             go_arg(GO_READ,      GO_CU, GO_POINTWISE),  & ! hu
             go_arg(GO_READ,      GO_CU,  GO_STENCIL(000,111,000)),  & ! sshn_u
             go_arg(GO_READ,      GO_GRID_MASK_T)     &
           /)

     !> This is a boundary-conditions kernel and therefore
     !! acts on all points of the domain rather than just
     !! those that are internal
     integer :: ITERATES_OVER = GO_ALL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = GO_OFFSET_NE

  contains
    procedure, nopass :: code => bc_flather_u_code
  end type bc_flather_u

  !=======================================

  type, extends(kernel_type) :: bc_flather_v
     type(go_arg), dimension(4) :: meta_args =  &
          (/ go_arg(GO_READWRITE, GO_CV, GO_STENCIL(010,010,010)),  & ! va
             go_arg(GO_READ,      GO_CV, GO_POINTWISE),  & ! hv
             go_arg(GO_READ,      GO_CV, GO_STENCIL(010,010,010)),  & ! sshn_v
             go_arg(GO_READ,      GO_GRID_MASK_T)     &
           /)

     !> This is a boundary-conditions kernel and therefore
     !! acts on all points of the domain rather than just
     !! those that are internal
     integer :: ITERATES_OVER = GO_ALL_PTS

     !> Although the staggering of variables used in an Arakawa
     !! C grid is well defined, the way in which they are indexed is
     !! an implementation choice. This can be thought of as choosing
     !! which grid-point types have the same (i,j) index as a T
     !! point. This kernel assumes that the U,V and F points that
     !! share the same index as a given T point are those immediately
     !! to the North and East of it.
     integer :: index_offset = GO_OFFSET_NE

  contains
    procedure, nopass :: code => bc_flather_v_code
  end type bc_flather_v

contains
  
  !================================================

  subroutine invoke_bc_ssh(istep, ssha)
    implicit none
    integer,            intent(in)    :: istep
    type(r2d_field),    intent(inout) :: ssha
    ! Locals
    integer  :: ji, jj
    
    DO jj = ssha%internal%ystart, ssha%internal%ystop
       DO ji = ssha%internal%xstart, ssha%internal%xstop
          call bc_ssh_code(ji, jj, &
                           istep, ssha%data, ssha%grid%tmask)
       END DO
    END DO

  end subroutine invoke_bc_ssh
  
  !================================================

  subroutine bc_ssh_code(ji, jj, istep, ssha, tmask)
    use model_mod, only: rdt
    implicit none
    integer, intent(in)  :: ji, jj
    integer, dimension(:,:),  intent(in)    :: tmask
    integer,                  intent(in)    :: istep
    real(go_wp), dimension(:,:), intent(inout) :: ssha
    ! Locals
    real(go_wp) :: amp_tide, omega_tide, rtime

    amp_tide   = 0.2_go_wp
    omega_tide = 2.0_go_wp * 3.14159_go_wp / (12.42_go_wp * 3600._go_wp)
    rtime = real(istep) * rdt

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
    type(r2d_field), intent(inout) :: ua
    ! Locals
    integer  :: ji, jj

    do jj = ua%whole%ystart, ua%whole%ystop, 1
       do ji = ua%whole%xstart, ua%whole%xstop, 1
          call bc_solid_u_code(ji, jj, ua%data, ua%grid%tmask)
       end do
    end do

  end subroutine invoke_bc_solid_u

  !================================================
  
  !> Kernel to apply solid boundary conditions for u-velocity
  subroutine bc_solid_u_code(ji, jj, ua, tmask)
    implicit none
    integer,                  intent(in)    :: ji, jj
    integer,  dimension(:,:), intent(in)    :: tmask
    real(go_wp), dimension(:,:), intent(inout) :: ua

    if(tmask(ji,jj) * tmask(ji+1,jj) == 0)then
       ua(ji,jj) = 0._go_wp
    end if

  end subroutine bc_solid_u_code
  
  !================================================

  !> Manual version of code to invoke the kernel
  !! that applies the solid-bc to a field on V pts.
  subroutine invoke_bc_solid_v(va)
    implicit none
    type(r2d_field), intent(inout) :: va
    ! Locals
    integer  :: ji, jj

    do jj = va%whole%ystart, va%whole%ystop, 1
       do ji = va%whole%xstart, va%whole%xstop, 1
          call bc_solid_v_code(ji,jj,va%data,va%grid%tmask)
      end do
    end do

  end subroutine invoke_bc_solid_v
  
  !================================================

  !> Kernel to apply solid boundary conditions for v-velocity
  subroutine bc_solid_v_code(ji, jj, va, tmask)
    implicit none
    integer,                 intent(in)    :: ji, jj
    integer, dimension(:,:), intent(in)    :: tmask
    real(go_wp),dimension(:,:), intent(inout) :: va

    if(tmask(ji,jj) * tmask(ji,jj+1) == 0)then
       va(ji,jj) = 0._go_wp
    end if

  end subroutine bc_solid_v_code
  
  !================================================

  !>                                  Du                 Dssh
  !!Flather open boundary condition [---- = sqrt(g/H) * ------]
  !!                                  Dn                 Dn
  !! ua and va in du/dn should be the specified tidal forcing
  subroutine invoke_bc_flather_u(ua, hu, sshn_u)
    implicit none
    type(r2d_field), intent(inout) :: ua
    type(r2d_field), intent(in) :: hu, sshn_u
    ! Locals
    integer  :: ji, jj

    ! Original loop was:
    !            DO jj = 1, jpj
    !              DO ji = 0, jpi  
    DO jj = ua%whole%ystart, ua%whole%ystop, 1
       DO ji = ua%whole%xstart, ua%whole%xstop, 1
          call bc_flather_u_code(ji,jj, &
                                 ua%data, hu%data, sshn_u%data, &
                                 ua%grid%tmask)
       END DO
    END DO
  
  end subroutine invoke_bc_flather_u
  
  !================================================

  !> Kernel to apply Flather condition to U
  subroutine bc_flather_u_code(ji, jj, ua, hu, sshn_u, tmask)
    use physical_params_mod, only: g
    implicit none
    integer,                  intent(in)    :: ji, jj
    integer,  dimension(:,:), intent(in)    :: tmask
    real(go_wp), dimension(:,:), intent(inout) :: ua
    real(go_wp), dimension(:,:), intent(in)    :: hu, sshn_u
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
    type(r2d_field), intent(inout) :: va
    type(r2d_field), intent(in)    :: hv, sshn_v
    ! Locals
    integer  :: ji, jj

    DO jj = va%whole%ystart, va%whole%ystop, 1
       DO ji = va%whole%xstart, va%whole%xstop, 1
          call bc_flather_v_code(ji,jj, &
                                 va%data, hv%data, sshn_v%data, &
                                 va%grid%tmask)
       END DO
    END DO

  end subroutine invoke_bc_flather_v

  !================================================

  !> Kernel to apply Flather boundary condition to v component
  !! of velocity
  subroutine bc_flather_v_code(ji, jj, va, hv, sshn_v, tmask)
    use physical_params_mod, only: g
    implicit none
    integer,                  intent(in) :: ji, jj
    integer,  dimension(:,:), intent(in) :: tmask
    real(go_wp), dimension(:,:), intent(inout) :: va
    real(go_wp), dimension(:,:), intent(in) :: hv, sshn_v
    ! Locals
    integer  :: jiv

    ! Check whether this point is inside the simulated domain
    !\todo I could set-up a V-mask using exactly the same code structure
    !! as below. Could then apply the BC and multiply by V-mask and thus
    !! remove conditionals => get vectorisation.
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

  subroutine setup_vmask_code(ji, jj, vmask, tmask)
    integer, intent(in) :: ji, jj
    integer,  dimension(:,:), intent(inout) :: vmask
    integer,  dimension(:,:), intent(in) :: tmask

    vmask(ji,jj) = 0

    IF(tmask(ji,jj) + tmask(ji,jj+1) <= -1) return
    
    IF(tmask(ji,jj) < 0) THEN
       vmask(ji, jj + 1) = 1
    ELSE IF(tmask(ji,jj+1) < 0) THEN
       vmask(ji, jj - 1) = 1
    END IF

  end subroutine setup_vmask_code

  !================================================

end module boundary_conditions_mod
