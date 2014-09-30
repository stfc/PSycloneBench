module time_step_mod
  use kind_params_mod
  use field_mod, only: copy_field
  use topology_mod, only: M, N, mp1, np1
  implicit none

contains

  subroutine invoke_time_step(cufld, cvfld, ufld, unew, uold, &
                              vfld, vnew, vold, &
                              pfld, pnew, pold, &
                              hfld, zfld, tdt)
    use compute_cu_mod,   only: compute_cu_code
    use compute_cv_mod,   only: compute_cv_code
    use compute_z_mod,    only: compute_z_code
    use compute_h_mod,    only: compute_h_code
    use compute_unew_mod, only: compute_unew_code
    use compute_vnew_mod, only: compute_vnew_code
    use compute_pnew_mod, only: compute_pnew_code
    use time_smooth_mod,  only: time_smooth_code
    use timing_mod,       only: timer_start, timer_stop
    implicit none
    real(wp), dimension(mp1,np1), intent(inout) :: cufld, cvfld
    real(wp), dimension(mp1,np1), intent(inout) :: unew, vnew, pnew
    real(wp), dimension(mp1,np1), intent(inout) :: hfld, zfld, pfld, &
                                                   ufld, vfld
    real(wp), dimension(mp1,np1), intent(inout) :: uold, vold, pold
    real(wp),                     intent(in) :: tdt
    ! Locals
    integer :: I, J

    !============================================
    ! COMPUTE CAPITAL U, CAPITAL V, Z AND H

    !call invoke_compute_fluxes(CU, CV, z, h, P, U, V)

    !CALL invoke_compute_cu(CU, P, U)
    ! M/N obtained from topology look-up
    do J= 1, N, 1
       do I = 1, M, 1

          call compute_cu_code(i+1, j, cufld, pfld, ufld)

          call compute_cv_code(i, j+1, cvfld, pfld, vfld)

          call compute_z_code(i+1, j+1, zfld, pfld, ufld, vfld)

          CALL compute_h_code(i, j, hfld, pfld, ufld, vfld)
       end do
    end do

    !============================================
    ! PERIODIC CONTINUATION

    !CALL invoke_apply_bcs_uvtf(cufld, cvfld, hfld, Zfld)
    !call invoke(periodic_bc(cu), periodic_bc(cv), ....)

    !call invoke_apply_bcs_cu(ufield)
    ! Ultimately, this can be generated by PSyclone but in the
    ! absence of that we implement it manually here...
    ! First col = last col
!!$    cufld(1, 1:N) = cufld(M+1,1:N)
!!$    ! Last col = first col
!!$    hfld(M+1,1:N) = hfld(1,  1:N)
!!$    ! First col = last col
!!$    zfld(1, 2:N+1) = zfld(M+1,  2:N+1)
!!$
!!$    do i = 1, M, 1
!!$
!!$       ! Last row = first row
!!$       cufld(i,N+1) = cufld(i,1)
!!$
!!$       !call invoke_apply_bcs_ct(tfield)
!!$       ! Last row = first row
!!$       hfld(i,N+1) = hfld(i,1)
!!$
!!$       !call invoke_apply_bcs_cv(vfield)
!!$       cvfld(i,1 ) = cvfld(i,N+1)
!!$
!!$       ! First row = last row
!!$       zfld(i,1)   = zfld(i,N+1)
!!$
!!$    end do
!!$    cufld(M+1,N+1) = cufld(M+1,1)
!!$    hfld(M+1,N+1)  = hfld(M+1,1)
!!$    zfld(M+1,1)    = zfld(M+1,N+1)
!!$
!!$    ! Last col = first col
!!$    cvfld(M+1,1:N+1) = cvfld(1,  1:N+1)

         DO J=1,N
            CUfld(1,J) = CUfld(M+1,J)
            CVfld(M+1,J+1) = CVfld(1,J+1)
            Zfld(1,J+1) = Zfld(M+1,J+1)
            Hfld(M+1,J) = Hfld(1,J)
         END DO
         DO I=1,M
            CUfld(I+1,N+1) = CUfld(I+1,1)
            CVfld(I,1) = CVfld(I,N+1)
            Zfld(I+1,1) = Zfld(I+1,N+1)
            Hfld(I,N+1) = Hfld(I,1)
         END DO
         CUfld(1,N+1) = CUfld(M+1,1)
         CVfld(M+1,1) = CVfld(1,N+1)
         Zfld(1,1) = Zfld(M+1,N+1)
         Hfld(M+1,N+1) = Hfld(1,1)

    !============================================
    ! COMPUTE NEW VALUES U,V AND P

    !CALL manual_invoke_compute_unew(unew, uold,  z, cv, h, tdt)
    DO J=1, N, 1
       DO I= 1, M, 1

          CALL compute_unew_code(i+1, j, unew, uold, &
                                 zfld, cvfld, hfld, tdt)

          CALL compute_vnew_code(i, j+1, vnew, vold, &
                                 zfld, cufld, hfld, tdt)

          CALL compute_pnew_code(i, j, pnew, pold, &
                                 cufld, cvfld, tdt)
       END DO
    END DO

    !============================================
    ! PERIODIC CONTINUATION
    !CALL invoke_apply_bcs_uvt(UNEW, VNEW, PNEW)

    !call invoke_apply_bcs_cu(unew)
    ! Ultimately, this can be generated by PSyclone but in the
    ! absence of that we implement it manually here...
    ! First col = last col
!!$    unew(1,  1:N) = unew(M+1,1:N)
!!$    ! Last col = first col
!!$    pnew(M+1,1:N) = pnew(1,  1:N)
!!$    do i=1,M
!!$       ! Last row = first row
!!$       unew(i, N+1) = unew(i, 1)
!!$       ! Last row = first row
!!$       pnew(i, N+1) = pnew(i, 1)
!!$
!!$       ! First row = last row
!!$       vnew(i, 1  ) = vnew(i, N+1)
!!$    end do
!!$    unew(M+1,N+1) = unew(M+1,1)
!!$    pnew(M+1,N+1) = pnew(M+1,1)
!!$
!!$    ! Last col = first col
!!$    vnew(M+1,1:N+1) = vnew(1,  1:N+1)

         DO J=1,N
            UNEW(1,J) = UNEW(M+1,J)
            VNEW(M+1,J+1) = VNEW(1,J+1)
            PNEW(M+1,J) = PNEW(1,J)
         END DO
         DO I=1,M
            UNEW(I+1,N+1) = UNEW(I+1,1)
            VNEW(I,1) = VNEW(I,N+1)
            PNEW(I,N+1) = PNEW(I,1)
         END DO
         UNEW(1,N+1) = UNEW(M+1,1)
         VNEW(M+1,1) = VNEW(1,N+1)
         PNEW(M+1,N+1) = PNEW(1,1)

    !============================================
    ! The time-smoothing is applied to a field at *every* grid point
    ! This updates the 'old' fields...

    ! Loop over 'columns'
    DO J=1,N+1 !idim2
      DO I=1,M+1 !idim1
        CALL time_smooth_code(i,j,ufld,unew,uold)
        CALL time_smooth_code(i,j,vfld,vnew,vold)
        CALL time_smooth_code(i,j,pfld,pnew,pold)
        Ufld(i,j) = UNEW(i,j)
        Vfld(i,j) = VNEW(i,j)
        Pfld(i,j) = PNEW(i,j)
      END DO
    END DO

    !============================================
    ! Update for next step
    !CALL copy_field(UNEW, Ufld)
    !CALL copy_field(VNEW, Vfld)
    !CALL copy_field(PNEW, Pfld)
    !Ufld(:,:) = UNEW(:,:)
    !Vfld(:,:) = VNEW(:,:)
    !Pfld(:,:) = PNEW(:,:)

  end subroutine invoke_time_step

end module time_step_mod
