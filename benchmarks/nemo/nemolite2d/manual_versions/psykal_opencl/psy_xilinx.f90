  MODULE psy_gocean2d
    USE field_mod
    USE kind_params_mod
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE continuity_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, &
&area_t, rdt)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, area_t
      REAL(KIND=go_wp), intent(in), target :: rdt
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the continuity_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(ssha_t), C_LOC(ssha_t))
      CALL check_status('clSetKernelArg: arg 4 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 5 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 6 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 7 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(hu), C_LOC(hu))
      CALL check_status('clSetKernelArg: arg 8 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 9, C_SIZEOF(hv), C_LOC(hv))
      CALL check_status('clSetKernelArg: arg 9 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 10, C_SIZEOF(un), C_LOC(un))
      CALL check_status('clSetKernelArg: arg 10 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 11, C_SIZEOF(vn), C_LOC(vn))
      CALL check_status('clSetKernelArg: arg 11 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 12, C_SIZEOF(area_t), C_LOC(area_t))
      CALL check_status('clSetKernelArg: arg 12 of continuity_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 13, C_SIZEOF(rdt), C_LOC(rdt))
      CALL check_status('clSetKernelArg: arg 13 of continuity_code', ierr)
    END SUBROUTINE continuity_code_set_args
    SUBROUTINE momentum_u_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, ua, un, vn, hu, hv, ht, ssha_u, sshn_t, sshn_u, &
&sshn_v, tmask, dx_u, dx_v, dx_t, dy_u, dy_t, area_u, gphiu, omega, d2r, g, rdt, cbfr, visc)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ua, un, vn, hu, hv, ht, ssha_u, sshn_t, sshn_u, sshn_v, tmask, dx_u, dx_v, &
&dx_t, dy_u, dy_t, area_u, gphiu
      REAL(KIND=go_wp), intent(in), target :: omega
      REAL(KIND=go_wp), intent(in), target :: d2r
      REAL(KIND=go_wp), intent(in), target :: g
      REAL(KIND=go_wp), intent(in), target :: rdt
      REAL(KIND=go_wp), intent(in), target :: cbfr
      REAL(KIND=go_wp), intent(in), target :: visc
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the momentum_u_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(ua), C_LOC(ua))
      CALL check_status('clSetKernelArg: arg 4 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(un), C_LOC(un))
      CALL check_status('clSetKernelArg: arg 5 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(vn), C_LOC(vn))
      CALL check_status('clSetKernelArg: arg 6 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(hu), C_LOC(hu))
      CALL check_status('clSetKernelArg: arg 7 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(hv), C_LOC(hv))
      CALL check_status('clSetKernelArg: arg 8 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 9, C_SIZEOF(ht), C_LOC(ht))
      CALL check_status('clSetKernelArg: arg 9 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 10, C_SIZEOF(ssha_u), C_LOC(ssha_u))
      CALL check_status('clSetKernelArg: arg 10 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 11, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 11 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 12, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 12 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 13, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 13 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 14, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 14 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 15, C_SIZEOF(dx_u), C_LOC(dx_u))
      CALL check_status('clSetKernelArg: arg 15 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 16, C_SIZEOF(dx_v), C_LOC(dx_v))
      CALL check_status('clSetKernelArg: arg 16 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 17, C_SIZEOF(dx_t), C_LOC(dx_t))
      CALL check_status('clSetKernelArg: arg 17 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 18, C_SIZEOF(dy_u), C_LOC(dy_u))
      CALL check_status('clSetKernelArg: arg 18 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 19, C_SIZEOF(dy_t), C_LOC(dy_t))
      CALL check_status('clSetKernelArg: arg 19 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 20, C_SIZEOF(area_u), C_LOC(area_u))
      CALL check_status('clSetKernelArg: arg 20 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 21, C_SIZEOF(gphiu), C_LOC(gphiu))
      CALL check_status('clSetKernelArg: arg 21 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 22, C_SIZEOF(omega), C_LOC(omega))
      CALL check_status('clSetKernelArg: arg 22 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 23, C_SIZEOF(d2r), C_LOC(d2r))
      CALL check_status('clSetKernelArg: arg 23 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 24, C_SIZEOF(g), C_LOC(g))
      CALL check_status('clSetKernelArg: arg 24 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 25, C_SIZEOF(rdt), C_LOC(rdt))
      CALL check_status('clSetKernelArg: arg 25 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 26, C_SIZEOF(cbfr), C_LOC(cbfr))
      CALL check_status('clSetKernelArg: arg 26 of momentum_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 27, C_SIZEOF(visc), C_LOC(visc))
      CALL check_status('clSetKernelArg: arg 27 of momentum_u_code', ierr)
    END SUBROUTINE momentum_u_code_set_args
    SUBROUTINE momentum_v_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, va, un, vn, hu, hv, ht, ssha_v, sshn_t, sshn_u, &
&sshn_v, tmask, dx_v, dx_t, dy_u, dy_v, dy_t, area_v, gphiv, omega, d2r, g, rdt, cbfr, visc)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: va, un, vn, hu, hv, ht, ssha_v, sshn_t, sshn_u, sshn_v, tmask, dx_v, dx_t, &
&dy_u, dy_v, dy_t, area_v, gphiv
      REAL(KIND=go_wp), intent(in), target :: omega
      REAL(KIND=go_wp), intent(in), target :: d2r
      REAL(KIND=go_wp), intent(in), target :: g
      REAL(KIND=go_wp), intent(in), target :: rdt
      REAL(KIND=go_wp), intent(in), target :: cbfr
      REAL(KIND=go_wp), intent(in), target :: visc
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the momentum_v_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(va), C_LOC(va))
      CALL check_status('clSetKernelArg: arg 4 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(un), C_LOC(un))
      CALL check_status('clSetKernelArg: arg 5 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(vn), C_LOC(vn))
      CALL check_status('clSetKernelArg: arg 6 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(hu), C_LOC(hu))
      CALL check_status('clSetKernelArg: arg 7 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(hv), C_LOC(hv))
      CALL check_status('clSetKernelArg: arg 8 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 9, C_SIZEOF(ht), C_LOC(ht))
      CALL check_status('clSetKernelArg: arg 9 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 10, C_SIZEOF(ssha_v), C_LOC(ssha_v))
      CALL check_status('clSetKernelArg: arg 10 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 11, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 11 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 12, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 12 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 13, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 13 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 14, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 14 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 15, C_SIZEOF(dx_v), C_LOC(dx_v))
      CALL check_status('clSetKernelArg: arg 15 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 16, C_SIZEOF(dx_t), C_LOC(dx_t))
      CALL check_status('clSetKernelArg: arg 16 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 17, C_SIZEOF(dy_u), C_LOC(dy_u))
      CALL check_status('clSetKernelArg: arg 17 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 18, C_SIZEOF(dy_v), C_LOC(dy_v))
      CALL check_status('clSetKernelArg: arg 18 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 19, C_SIZEOF(dy_t), C_LOC(dy_t))
      CALL check_status('clSetKernelArg: arg 19 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 20, C_SIZEOF(area_v), C_LOC(area_v))
      CALL check_status('clSetKernelArg: arg 20 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 21, C_SIZEOF(gphiv), C_LOC(gphiv))
      CALL check_status('clSetKernelArg: arg 21 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 22, C_SIZEOF(omega), C_LOC(omega))
      CALL check_status('clSetKernelArg: arg 22 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 23, C_SIZEOF(d2r), C_LOC(d2r))
      CALL check_status('clSetKernelArg: arg 23 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 24, C_SIZEOF(g), C_LOC(g))
      CALL check_status('clSetKernelArg: arg 24 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 25, C_SIZEOF(rdt), C_LOC(rdt))
      CALL check_status('clSetKernelArg: arg 25 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 26, C_SIZEOF(cbfr), C_LOC(cbfr))
      CALL check_status('clSetKernelArg: arg 26 of momentum_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 27, C_SIZEOF(visc), C_LOC(visc))
      CALL check_status('clSetKernelArg: arg 27 of momentum_v_code', ierr)
    END SUBROUTINE momentum_v_code_set_args
    SUBROUTINE bc_ssh_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, istp, ssha_t, tmask, rdt)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ssha_t, tmask
      INTEGER, intent(in), target :: istp
      REAL(KIND=go_wp), intent(in), target :: rdt
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_ssh_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(istp), C_LOC(istp))
      CALL check_status('clSetKernelArg: arg 4 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(ssha_t), C_LOC(ssha_t))
      CALL check_status('clSetKernelArg: arg 5 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 6 of bc_ssh_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(rdt), C_LOC(rdt))
      CALL check_status('clSetKernelArg: arg 7 of bc_ssh_code', ierr)
    END SUBROUTINE bc_ssh_code_set_args
    SUBROUTINE bc_solid_u_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, ua, tmask)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ua, tmask
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_solid_u_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(ua), C_LOC(ua))
      CALL check_status('clSetKernelArg: arg 4 of bc_solid_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 5 of bc_solid_u_code', ierr)
    END SUBROUTINE bc_solid_u_code_set_args
    SUBROUTINE bc_solid_v_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, va, tmask)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: va, tmask
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_solid_v_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(va), C_LOC(va))
      CALL check_status('clSetKernelArg: arg 4 of bc_solid_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 5 of bc_solid_v_code', ierr)
    END SUBROUTINE bc_solid_v_code_set_args
    SUBROUTINE bc_flather_u_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, ua, hu, sshn_u, tmask, g)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: ua, hu, sshn_u, tmask
      REAL(KIND=go_wp), intent(in), target :: g
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_flather_u_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(ua), C_LOC(ua))
      CALL check_status('clSetKernelArg: arg 4 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(hu), C_LOC(hu))
      CALL check_status('clSetKernelArg: arg 5 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 6 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 7 of bc_flather_u_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(g), C_LOC(g))
      CALL check_status('clSetKernelArg: arg 8 of bc_flather_u_code', ierr)
    END SUBROUTINE bc_flather_u_code_set_args
    SUBROUTINE bc_flather_v_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, va, hv, sshn_v, tmask, g)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: va, hv, sshn_v, tmask
      REAL(KIND=go_wp), intent(in), target :: g
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the bc_flather_v_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(va), C_LOC(va))
      CALL check_status('clSetKernelArg: arg 4 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(hv), C_LOC(hv))
      CALL check_status('clSetKernelArg: arg 5 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 6 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 7 of bc_flather_v_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(g), C_LOC(g))
      CALL check_status('clSetKernelArg: arg 8 of bc_flather_v_code', ierr)
    END SUBROUTINE bc_flather_v_code_set_args
    SUBROUTINE field_copy_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, un, ua)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: un, ua
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the field_copy_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(un), C_LOC(un))
      CALL check_status('clSetKernelArg: arg 4 of field_copy_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(ua), C_LOC(ua))
      CALL check_status('clSetKernelArg: arg 5 of field_copy_code', ierr)
    END SUBROUTINE field_copy_code_set_args
    SUBROUTINE next_sshu_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, sshn_u, sshn_t, tmask, area_t, area_u)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: sshn_u, sshn_t, tmask, area_t, area_u
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the next_sshu_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(sshn_u), C_LOC(sshn_u))
      CALL check_status('clSetKernelArg: arg 4 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 5 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 6 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(area_t), C_LOC(area_t))
      CALL check_status('clSetKernelArg: arg 7 of next_sshu_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(area_u), C_LOC(area_u))
      CALL check_status('clSetKernelArg: arg 8 of next_sshu_code', ierr)
    END SUBROUTINE next_sshu_code_set_args
    SUBROUTINE next_sshv_code_set_args(kernel_obj, xstart, xstop, ystart, ystop, sshn_v, sshn_t, tmask, area_t, area_v)
      USE clfortran, ONLY: clSetKernelArg
      USE iso_c_binding, ONLY: c_sizeof, c_loc, c_intptr_t
      USE ocl_utils_mod, ONLY: check_status
      INTEGER, intent(in), target :: xstart, xstop, ystart, ystop
      INTEGER(KIND=c_intptr_t), intent(in), target :: sshn_v, sshn_t, tmask, area_t, area_v
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), target :: kernel_obj
      ! Set the arguments for the next_sshv_code OpenCL Kernel
      ierr = clSetKernelArg(kernel_obj, 0, C_SIZEOF(xstart), C_LOC(xstart))
      CALL check_status('clSetKernelArg: arg 0 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 1, C_SIZEOF(xstop), C_LOC(xstop))
      CALL check_status('clSetKernelArg: arg 1 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 2, C_SIZEOF(ystart), C_LOC(ystart))
      CALL check_status('clSetKernelArg: arg 2 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 3, C_SIZEOF(ystop), C_LOC(ystop))
      CALL check_status('clSetKernelArg: arg 3 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 4, C_SIZEOF(sshn_v), C_LOC(sshn_v))
      CALL check_status('clSetKernelArg: arg 4 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 5, C_SIZEOF(sshn_t), C_LOC(sshn_t))
      CALL check_status('clSetKernelArg: arg 5 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 6, C_SIZEOF(tmask), C_LOC(tmask))
      CALL check_status('clSetKernelArg: arg 6 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 7, C_SIZEOF(area_t), C_LOC(area_t))
      CALL check_status('clSetKernelArg: arg 7 of next_sshv_code', ierr)
      ierr = clSetKernelArg(kernel_obj, 8, C_SIZEOF(area_v), C_LOC(area_v))
      CALL check_status('clSetKernelArg: arg 8 of next_sshv_code', ierr)
    END SUBROUTINE next_sshv_code_set_args
    SUBROUTINE invoke_0(ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v, istp)
      USE fortcl, ONLY: create_rw_buffer
      USE fortcl, ONLY: get_num_cmd_queues, get_cmd_queues, get_kernel_by_name
      USE clfortran
      USE iso_c_binding
      USE physical_params_mod, ONLY: omega, d2r, g
      USE model_mod, ONLY: rdt, cbfr, visc
      TYPE(r2d_field), intent(inout), target :: ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, un, vn, ua, ht, ssha_u, va, ssha_v
      INTEGER, intent(inout) :: istp
      INTEGER(KIND=c_size_t), target :: localsize_12(2)
      INTEGER(KIND=c_size_t), target :: globalsize_12(2)
      INTEGER(KIND=c_size_t), target :: localsize_11(2)
      INTEGER(KIND=c_size_t), target :: globalsize_11(2)
      INTEGER(KIND=c_size_t), target :: localsize_10(2)
      INTEGER(KIND=c_size_t), target :: globalsize_10(2)
      INTEGER(KIND=c_size_t), target :: localsize_9(2)
      INTEGER(KIND=c_size_t), target :: globalsize_9(2)
      INTEGER(KIND=c_size_t), target :: localsize_8(2)
      INTEGER(KIND=c_size_t), target :: globalsize_8(2)
      INTEGER(KIND=c_size_t), target :: localsize_7(2)
      INTEGER(KIND=c_size_t), target :: globalsize_7(2)
      INTEGER(KIND=c_size_t), target :: localsize_6(2)
      INTEGER(KIND=c_size_t), target :: globalsize_6(2)
      INTEGER(KIND=c_size_t), target :: localsize_5(2)
      INTEGER(KIND=c_size_t), target :: globalsize_5(2)
      INTEGER(KIND=c_size_t), target :: localsize_4(2)
      INTEGER(KIND=c_size_t), target :: globalsize_4(2)
      INTEGER(KIND=c_size_t), target :: localsize_3(2)
      INTEGER(KIND=c_size_t), target :: globalsize_3(2)
      INTEGER(KIND=c_size_t), target :: localsize_2(2)
      INTEGER(KIND=c_size_t), target :: globalsize_2(2)
      INTEGER(KIND=c_size_t), target :: localsize_1(2)
      INTEGER(KIND=c_size_t), target :: globalsize_1(2)
      INTEGER(KIND=c_size_t), target :: localsize(2)
      INTEGER(KIND=c_size_t), target :: globalsize(2)
      INTEGER(KIND=c_intptr_t), target :: write_event
      INTEGER(KIND=c_size_t) size_in_bytes
      INTEGER(KIND=c_intptr_t), target, save :: kernel_next_sshv_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_next_sshu_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_field_copy_code_2
      INTEGER(KIND=c_intptr_t), target, save :: kernel_field_copy_code_1
      INTEGER(KIND=c_intptr_t), target, save :: kernel_field_copy_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_flather_v_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_flather_u_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_solid_v_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_solid_u_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_bc_ssh_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_momentum_v_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_momentum_u_code
      INTEGER(KIND=c_intptr_t), target, save :: kernel_continuity_code
      LOGICAL, save :: first_time=.true.
      INTEGER ierr
      INTEGER(KIND=c_intptr_t), pointer, save :: cmd_queues(:)
      INTEGER, save :: num_cmd_queues
      IF (first_time) THEN
        first_time = .false.
        ! Ensure OpenCL run-time is initialised for this PSy-layer module
        CALL psy_init
        num_cmd_queues = get_num_cmd_queues()
        cmd_queues => get_cmd_queues()
        kernel_continuity_code = get_kernel_by_name("continuity_code")
        kernel_momentum_u_code = get_kernel_by_name("momentum_u_code")
        kernel_momentum_v_code = get_kernel_by_name("momentum_v_code")
        kernel_bc_ssh_code = get_kernel_by_name("bc_ssh_code")
        kernel_bc_solid_u_code = get_kernel_by_name("bc_solid_u_code")
        kernel_bc_solid_v_code = get_kernel_by_name("bc_solid_v_code")
        kernel_bc_flather_u_code = get_kernel_by_name("bc_flather_u_code")
        kernel_bc_flather_v_code = get_kernel_by_name("bc_flather_v_code")
        kernel_field_copy_code = get_kernel_by_name("field_copy_code")
        kernel_field_copy_code_1 = get_kernel_by_name("field_copy_code")
        kernel_field_copy_code_2 = get_kernel_by_name("field_copy_code")
        kernel_next_sshu_code = get_kernel_by_name("next_sshu_code")
        kernel_next_sshv_code = get_kernel_by_name("next_sshv_code")

        ! Create field device_ptr
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(ssha_t%data(1,1))
        ssha_t%device_ptr = create_rw_buffer(size_in_bytes)
        sshn_t%device_ptr = create_rw_buffer(size_in_bytes)
        sshn_u%device_ptr = create_rw_buffer(size_in_bytes)
        sshn_v%device_ptr = create_rw_buffer(size_in_bytes)
        hu%device_ptr = create_rw_buffer(size_in_bytes)
        hv%device_ptr = create_rw_buffer(size_in_bytes)
        un%device_ptr = create_rw_buffer(size_in_bytes)
        vn%device_ptr = create_rw_buffer(size_in_bytes)
        ua%device_ptr = create_rw_buffer(size_in_bytes)
        ssha_u%device_ptr = create_rw_buffer(size_in_bytes)
        va%device_ptr = create_rw_buffer(size_in_bytes)
        ssha_v%device_ptr = create_rw_buffer(size_in_bytes)
        ht%device_ptr = create_rw_buffer(size_in_bytes)

        ! Create grid device_ptr
        sshn_t%grid%dx_u_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%dx_v_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%dx_t_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%dy_u_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%dy_v_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%dy_t_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%area_u_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%area_v_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%area_t_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%gphiu_device = create_rw_buffer(size_in_bytes)
        sshn_t%grid%gphiv_device = create_rw_buffer(size_in_bytes)

        size_in_bytes = int(ssha_t%grid%nx*ssha_t%grid%ny, 8)*c_sizeof(ssha_t%grid%tmask(1,1))
        sshn_t%grid%tmask_device = create_rw_buffer(size_in_bytes)
        ierr = clFinish(cmd_queues(1))

      END IF

        CALL continuity_code_set_args(kernel_continuity_code, ssha_t%internal%xstart - 1, ssha_t%internal%xstop - 1, &
            &ssha_t%internal%ystart - 1, ssha_t%internal%ystop - 1, ssha_t%device_ptr, sshn_t%device_ptr, sshn_u%device_ptr, &
            &sshn_v%device_ptr, hu%device_ptr, hv%device_ptr, un%device_ptr, vn%device_ptr, sshn_t%grid%area_t_device, rdt)
      CALL momentum_u_code_set_args(kernel_momentum_u_code, ua%internal%xstart - 1, ua%internal%xstop - 1, ua%internal%ystart - 1, &
&ua%internal%ystop - 1, ua%device_ptr, un%device_ptr, vn%device_ptr, hu%device_ptr, hv%device_ptr, ht%device_ptr, &
&ssha_u%device_ptr, sshn_t%device_ptr, sshn_u%device_ptr, sshn_v%device_ptr, un%grid%tmask_device, un%grid%dx_u_device, &
&un%grid%dx_v_device, un%grid%dx_t_device, un%grid%dy_u_device, un%grid%dy_t_device, un%grid%area_u_device, un%grid%gphiu_device, &
&omega, d2r, g, rdt, cbfr, visc)
      CALL momentum_v_code_set_args(kernel_momentum_v_code, va%internal%xstart - 1, va%internal%xstop - 1, va%internal%ystart - 1, &
&va%internal%ystop - 1, va%device_ptr, un%device_ptr, vn%device_ptr, hu%device_ptr, hv%device_ptr, ht%device_ptr, &
&ssha_v%device_ptr, sshn_t%device_ptr, sshn_u%device_ptr, sshn_v%device_ptr, un%grid%tmask_device, un%grid%dx_v_device, &
&un%grid%dx_t_device, un%grid%dy_u_device, un%grid%dy_v_device, un%grid%dy_t_device, un%grid%area_v_device, un%grid%gphiv_device, &
&omega, d2r, g, rdt, cbfr, visc)
      CALL bc_ssh_code_set_args(kernel_bc_ssh_code, ssha_t%internal%xstart - 1, ssha_t%internal%xstop - 1, &
&ssha_t%internal%ystart - 1, ssha_t%internal%ystop - 1, istp, ssha_t%device_ptr, ssha_t%grid%tmask_device, rdt)
      CALL bc_solid_u_code_set_args(kernel_bc_solid_u_code, ua%whole%xstart - 1, ua%whole%xstop - 1, ua%whole%ystart - 1, &
&ua%whole%ystop - 1, ua%device_ptr, ua%grid%tmask_device)
      CALL bc_solid_v_code_set_args(kernel_bc_solid_v_code, va%whole%xstart - 1, va%whole%xstop - 1, va%whole%ystart - 1, &
&va%whole%ystop - 1, va%device_ptr, va%grid%tmask_device)
      CALL bc_flather_u_code_set_args(kernel_bc_flather_u_code, ua%whole%xstart - 1, ua%whole%xstop - 1, ua%whole%ystart - 1, &
&ua%whole%ystop - 1, ua%device_ptr, hu%device_ptr, sshn_u%device_ptr, hu%grid%tmask_device, g)
      CALL bc_flather_v_code_set_args(kernel_bc_flather_v_code, va%whole%xstart - 1, va%whole%xstop - 1, va%whole%ystart - 1, &
&va%whole%ystop - 1, va%device_ptr, hv%device_ptr, sshn_v%device_ptr, hv%grid%tmask_device, g)
      CALL field_copy_code_set_args(kernel_field_copy_code, 1 - 1, SIZE(un%data, 1) - 1, 1 - 1, SIZE(un%data, 2) - 1, &
&un%device_ptr, ua%device_ptr)
      CALL field_copy_code_set_args(kernel_field_copy_code, 1 - 1, SIZE(vn%data, 1) - 1, 1 - 1, SIZE(vn%data, 2) - 1, &
&vn%device_ptr, va%device_ptr)
      CALL field_copy_code_set_args(kernel_field_copy_code, 1 - 1, SIZE(sshn_t%data, 1) - 1, 1 - 1, SIZE(sshn_t%data, 2) - 1, &
&sshn_t%device_ptr, ssha_t%device_ptr)
      CALL next_sshu_code_set_args(kernel_next_sshu_code, sshn_u%internal%xstart - 1, sshn_u%internal%xstop - 1, &
&sshn_u%internal%ystart - 1, sshn_u%internal%ystop - 1, sshn_u%device_ptr, sshn_t%device_ptr, sshn_t%grid%tmask_device, &
&sshn_t%grid%area_t_device, sshn_t%grid%area_u_device)
      CALL next_sshv_code_set_args(kernel_next_sshv_code, sshn_v%internal%xstart - 1, sshn_v%internal%xstop - 1, &
&sshn_v%internal%ystart - 1, sshn_v%internal%ystop - 1, sshn_v%device_ptr, sshn_t%device_ptr, sshn_t%grid%tmask_device, &
&sshn_t%grid%area_t_device, sshn_t%grid%area_v_device)

      ! Ensure field data is on device
      IF (.NOT. ssha_t%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(ssha_t%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_t%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ssha_t%data_on_device = .true.
        ssha_t%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_t%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_t%data_on_device = .true.
        sshn_t%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_u%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_u%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_u%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_u%data_on_device = .true.
        sshn_u%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_v%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_v%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_v%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_v%data_on_device = .true.
        sshn_v%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. hu%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(hu%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hu%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hu%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        hu%data_on_device = .true.
        hu%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. hv%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(hv%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hv%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hv%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        hv%data_on_device = .true.
        hv%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. un%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(un%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        un%data_on_device = .true.
        un%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. vn%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(vn%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), vn%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(vn%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        vn%data_on_device = .true.
        vn%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (sshn_t%grid%area_t_device == 0) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%area_t(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%area_t_device, CL_TRUE, 0_8, size_in_bytes, &
&C_LOC(sshn_t%grid%area_t), 0, C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize = (/sshn_t%grid%nx, sshn_t%grid%ny/)
      localsize = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_continuity_code, 2, C_NULL_PTR, C_LOC(globalsize), C_LOC(localsize), 0, &
&C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. ua%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(ua%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        ua%data_on_device = .true.
        ua%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. un%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        un%data_on_device = .true.
        un%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. vn%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(vn%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), vn%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(vn%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        vn%data_on_device = .true.
        vn%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. hu%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(hu%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hu%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hu%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        hu%data_on_device = .true.
        hu%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. hv%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(hv%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hv%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hv%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        hv%data_on_device = .true.
        hv%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. ht%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(ht%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ht%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ht%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        ht%data_on_device = .true.
        ht%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. ssha_u%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(ssha_u%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_u%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_u%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ssha_u%data_on_device = .true.
        ssha_u%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_t%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_t%data_on_device = .true.
        sshn_t%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_u%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(sshn_u%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_u%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_u%data_on_device = .true.
        sshn_u%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_v%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(sshn_v%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_v%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_v%data_on_device = .true.
        sshn_v%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%tmask_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%tmask(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%tmask), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dx_u_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_u(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_u_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_u), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dx_v_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_v(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_v_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_v), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dx_t_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_t(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_t_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_t), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dy_u_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_u(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_u_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_u), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dy_t_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_t(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_t_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_t), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%area_u_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%area_u(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%area_u_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%area_u), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%gphiu_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%gphiu(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%gphiu_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%gphiu), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_1 = (/un%grid%nx, un%grid%ny/)
      localsize_1 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_momentum_u_code, 2, C_NULL_PTR, C_LOC(globalsize_1), C_LOC(localsize_1), &
&0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. va%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(va%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), va%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        va%data_on_device = .true.
        va%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. un%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        un%data_on_device = .true.
        un%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. vn%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(vn%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), vn%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(vn%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        vn%data_on_device = .true.
        vn%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. hu%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(hu%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hu%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hu%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        hu%data_on_device = .true.
        hu%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. hv%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(hv%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hv%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hv%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        hv%data_on_device = .true.
        hv%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. ht%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(ht%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ht%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ht%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        ht%data_on_device = .true.
        ht%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. ssha_v%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(ssha_v%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_v%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ssha_v%data_on_device = .true.
        ssha_v%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_t%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_t%data_on_device = .true.
        sshn_t%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_u%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(sshn_u%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_u%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_u%data_on_device = .true.
        sshn_u%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_v%data_on_device) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(sshn_v%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_v%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_v%data_on_device = .true.
        sshn_v%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%tmask_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%tmask(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%tmask), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dx_v_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_v(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_v_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_v), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dx_t_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dx_t(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dx_t_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dx_t), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dy_u_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_u(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_u_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_u), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dy_v_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_v(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_v_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_v), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%dy_t_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%dy_t(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%dy_t_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%dy_t), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%area_v_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%area_v(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%area_v_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%area_v), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (un%grid%gphiv_device == 0) THEN
        size_in_bytes = int(un%grid%nx*un%grid%ny, 8)*c_sizeof(un%grid%gphiv(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%grid%gphiv_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%grid%gphiv), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_2 = (/un%grid%nx, un%grid%ny/)
      localsize_2 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_momentum_v_code, 2, C_NULL_PTR, C_LOC(globalsize_2), C_LOC(localsize_2), &
&0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. ssha_t%data_on_device) THEN
        size_in_bytes = int(ssha_t%grid%nx*ssha_t%grid%ny, 8)*c_sizeof(ssha_t%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_t%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ssha_t%data_on_device = .true.
        ssha_t%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (ssha_t%grid%tmask_device == 0) THEN
        size_in_bytes = int(ssha_t%grid%nx*ssha_t%grid%ny, 8)*c_sizeof(ssha_t%grid%tmask(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_t%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, &
&C_LOC(ssha_t%grid%tmask), 0, C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_3 = (/ssha_t%grid%nx, ssha_t%grid%ny/)
      localsize_3 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_ssh_code, 2, C_NULL_PTR, C_LOC(globalsize_3), C_LOC(localsize_3), 0, &
&C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. ua%data_on_device) THEN
        size_in_bytes = int(ua%grid%nx*ua%grid%ny, 8)*c_sizeof(ua%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        ua%data_on_device = .true.
        ua%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (ua%grid%tmask_device == 0) THEN
        size_in_bytes = int(ua%grid%nx*ua%grid%ny, 8)*c_sizeof(ua%grid%tmask(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%grid%tmask), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_4 = (/ua%grid%nx, ua%grid%ny/)
      localsize_4 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_solid_u_code, 2, C_NULL_PTR, C_LOC(globalsize_4), C_LOC(localsize_4), &
&0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. va%data_on_device) THEN
        size_in_bytes = int(va%grid%nx*va%grid%ny, 8)*c_sizeof(va%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), va%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        va%data_on_device = .true.
        va%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (va%grid%tmask_device == 0) THEN
        size_in_bytes = int(va%grid%nx*va%grid%ny, 8)*c_sizeof(va%grid%tmask(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), va%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%grid%tmask), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_5 = (/va%grid%nx, va%grid%ny/)
      localsize_5 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_solid_v_code, 2, C_NULL_PTR, C_LOC(globalsize_5), C_LOC(localsize_5), &
&0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. ua%data_on_device) THEN
        size_in_bytes = int(hu%grid%nx*hu%grid%ny, 8)*c_sizeof(ua%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        ua%data_on_device = .true.
        ua%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. hu%data_on_device) THEN
        size_in_bytes = int(hu%grid%nx*hu%grid%ny, 8)*c_sizeof(hu%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hu%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hu%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        hu%data_on_device = .true.
        hu%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_u%data_on_device) THEN
        size_in_bytes = int(hu%grid%nx*hu%grid%ny, 8)*c_sizeof(sshn_u%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_u%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_u%data_on_device = .true.
        sshn_u%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (hu%grid%tmask_device == 0) THEN
        size_in_bytes = int(hu%grid%nx*hu%grid%ny, 8)*c_sizeof(hu%grid%tmask(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hu%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(hu%grid%tmask), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_6 = (/hu%grid%nx, hu%grid%ny/)
      localsize_6 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_flather_u_code, 2, C_NULL_PTR, C_LOC(globalsize_6), &
&C_LOC(localsize_6), 0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. va%data_on_device) THEN
        size_in_bytes = int(hv%grid%nx*hv%grid%ny, 8)*c_sizeof(va%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), va%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        va%data_on_device = .true.
        va%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. hv%data_on_device) THEN
        size_in_bytes = int(hv%grid%nx*hv%grid%ny, 8)*c_sizeof(hv%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hv%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(hv%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        hv%data_on_device = .true.
        hv%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_v%data_on_device) THEN
        size_in_bytes = int(hv%grid%nx*hv%grid%ny, 8)*c_sizeof(sshn_v%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_v%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_v%data_on_device = .true.
        sshn_v%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (hv%grid%tmask_device == 0) THEN
        size_in_bytes = int(hv%grid%nx*hv%grid%ny, 8)*c_sizeof(hv%grid%tmask(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), hv%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, C_LOC(hv%grid%tmask), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_7 = (/hv%grid%nx, hv%grid%ny/)
      localsize_7 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_bc_flather_v_code, 2, C_NULL_PTR, C_LOC(globalsize_7), &
&C_LOC(localsize_7), 0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. un%data_on_device) THEN
        size_in_bytes = int(ua%grid%nx*ua%grid%ny, 8)*c_sizeof(un%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), un%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(un%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        un%data_on_device = .true.
        un%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. ua%data_on_device) THEN
        size_in_bytes = int(ua%grid%nx*ua%grid%ny, 8)*c_sizeof(ua%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ua%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ua%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        ua%data_on_device = .true.
        ua%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_8 = (/ua%grid%nx, ua%grid%ny/)
      localsize_8 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_field_copy_code, 2, C_NULL_PTR, C_LOC(globalsize_8), C_LOC(localsize_8), &
&0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. vn%data_on_device) THEN
        size_in_bytes = int(va%grid%nx*va%grid%ny, 8)*c_sizeof(vn%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), vn%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(vn%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        vn%data_on_device = .true.
        vn%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. va%data_on_device) THEN
        size_in_bytes = int(va%grid%nx*va%grid%ny, 8)*c_sizeof(va%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), va%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(va%data), 0, C_NULL_PTR, &
&C_LOC(write_event))
        va%data_on_device = .true.
        va%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_9 = (/va%grid%nx, va%grid%ny/)
      localsize_9 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_field_copy_code, 2, C_NULL_PTR, C_LOC(globalsize_9), C_LOC(localsize_9), &
&0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. sshn_t%data_on_device) THEN
        size_in_bytes = int(ssha_t%grid%nx*ssha_t%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_t%data_on_device = .true.
        sshn_t%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. ssha_t%data_on_device) THEN
        size_in_bytes = int(ssha_t%grid%nx*ssha_t%grid%ny, 8)*c_sizeof(ssha_t%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), ssha_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(ssha_t%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        ssha_t%data_on_device = .true.
        ssha_t%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_10 = (/ssha_t%grid%nx, ssha_t%grid%ny/)
      localsize_10 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_field_copy_code, 2, C_NULL_PTR, C_LOC(globalsize_10), &
&C_LOC(localsize_10), 0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. sshn_u%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_u%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_u%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_u%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_u%data_on_device = .true.
        sshn_u%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_t%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_t%data_on_device = .true.
        sshn_t%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (sshn_t%grid%tmask_device == 0) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%tmask(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, &
&C_LOC(sshn_t%grid%tmask), 0, C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (sshn_t%grid%area_t_device == 0) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%area_t(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%area_t_device, CL_TRUE, 0_8, size_in_bytes, &
&C_LOC(sshn_t%grid%area_t), 0, C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (sshn_t%grid%area_u_device == 0) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%area_u(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%area_u_device, CL_TRUE, 0_8, size_in_bytes, &
&C_LOC(sshn_t%grid%area_u), 0, C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_11 = (/sshn_t%grid%nx, sshn_t%grid%ny/)
      localsize_11 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_next_sshu_code, 2, C_NULL_PTR, C_LOC(globalsize_11), &
&C_LOC(localsize_11), 0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Ensure field data is on device
      IF (.NOT. sshn_v%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_v%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_v%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_v%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_v%data_on_device = .true.
        sshn_v%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (.NOT. sshn_t%data_on_device) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%data(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%device_ptr, CL_TRUE, 0_8, size_in_bytes, C_LOC(sshn_t%data), 0, &
&C_NULL_PTR, C_LOC(write_event))
        sshn_t%data_on_device = .true.
        sshn_t%read_from_device_f => read_from_device
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (sshn_t%grid%tmask_device == 0) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%tmask(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%tmask_device, CL_TRUE, 0_8, size_in_bytes, &
&C_LOC(sshn_t%grid%tmask), 0, C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (sshn_t%grid%area_t_device == 0) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%area_t(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%area_t_device, CL_TRUE, 0_8, size_in_bytes, &
&C_LOC(sshn_t%grid%area_t), 0, C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      IF (sshn_t%grid%area_v_device == 0) THEN
        size_in_bytes = int(sshn_t%grid%nx*sshn_t%grid%ny, 8)*c_sizeof(sshn_t%grid%area_v(1,1))
        ! Create buffer on device
        ierr = clEnqueueWriteBuffer(cmd_queues(1), sshn_t%grid%area_v_device, CL_TRUE, 0_8, size_in_bytes, &
&C_LOC(sshn_t%grid%area_v), 0, C_NULL_PTR, C_LOC(write_event))
        ! Block until data copies have finished
        ierr = clFinish(cmd_queues(1))
      END IF
      globalsize_12 = (/sshn_t%grid%nx, sshn_t%grid%ny/)
      localsize_12 = (/64, 1/)
      ! Launch the kernel
      ierr = clEnqueueNDRangeKernel(cmd_queues(1), kernel_next_sshv_code, 2, C_NULL_PTR, C_LOC(globalsize_12), &
&C_LOC(localsize_12), 0, C_NULL_PTR, C_NULL_PTR)
      !
      ! Block until all kernels have finished
      ierr = clFinish(cmd_queues(1))
    END SUBROUTINE invoke_0
    SUBROUTINE read_from_device(from, to, nx, ny, width)
      USE iso_c_binding, ONLY: c_intptr_t
      USE fortcl, ONLY: read_buffer
      INTEGER(KIND=c_intptr_t), intent(in) :: from
      REAL(KIND=go_wp), intent(inout), dimension(:,:) :: to
      INTEGER, intent(in) :: nx, ny, width
      CALL read_buffer(from, to, int(width*ny, kind=8))
    END SUBROUTINE read_from_device
    SUBROUTINE psy_init()
      USE fortcl, ONLY: ocl_env_init, add_kernels
      CHARACTER(LEN=30) kernel_names(11)
      INTEGER :: ocl_device_num=1
      LOGICAL, save :: initialised=.False.
      ! Check to make sure we only execute this routine once
      IF (.not. initialised) THEN
        initialised = .True.
        ! Initialise the OpenCL environment/device
        CALL ocl_env_init(1, ocl_device_num, .False., .False.)
        ! The kernels this PSy layer module requires
        kernel_names(1) = "continuity_code"
        kernel_names(2) = "momentum_u_code"
        kernel_names(3) = "momentum_v_code"
        kernel_names(4) = "bc_ssh_code"
        kernel_names(5) = "bc_solid_u_code"
        kernel_names(6) = "bc_solid_v_code"
        kernel_names(7) = "bc_flather_u_code"
        kernel_names(8) = "bc_flather_v_code"
        kernel_names(9) = "field_copy_code"
        kernel_names(10) = "next_sshu_code"
        kernel_names(11) = "next_sshv_code"
        ! Create the OpenCL kernel objects. Expects to find all of the compiled
        ! kernels in FORTCL_KERNELS_FILE.
        CALL add_kernels(11, kernel_names)
      END IF
    END SUBROUTINE psy_init
  END MODULE psy_gocean2d
