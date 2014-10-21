module time_step_mod
  use kind_params_mod
  use topology_mod, only: M, N
  implicit none

contains

  subroutine invoke_time_step(cufld, cvfld, ufld, unew, uold, &
                              vfld, vnew, vold, &
                              pfld, pnew, pold, &
                              hfld, zfld, tdt)
    implicit none
    real(wp), dimension(M+1,N+1), intent(inout) :: cufld, cvfld
    real(wp), dimension(M+1,N+1), intent(inout) :: unew, vnew, pnew
    real(wp), dimension(M+1,N+1), intent(inout) :: hfld, zfld, pfld, &
                                               ufld, vfld
    !> \todo Do we need these 'old' arrays as args or are they
    !! local workspace really?
    real(wp), dimension(M+1,N+1), intent(inout) :: uold, vold, pold
    real(wp),                     intent(in) :: tdt
    ! Locals
    integer :: I, J

    !============================================
    ! COMPUTE CAPITAL U, CAPITAL V, Z AND H

    ! M/N obtained from topology look-up
    do J= 1, N, 1
       do I = 1, M, 1

          call compute_cu_code(i+1, j, cufld, pfld, ufld)
!       end do
!    end do

 !   do J= 1, N, 1
 !      do I= 1, M, 1

          call compute_cv_code(i, j+1, cvfld, pfld, vfld)
       end do
    end do

    do J= 1, N, 1
       do I= 1, M, 1

          call compute_z_code(i+1, j+1, zfld, pfld, ufld, vfld)
!       end do
!    end do

!    DO J= 1, N, 1
!       DO I= 1, M, 1

          CALL compute_h_code(i, j, hfld, pfld, ufld, vfld)
       END DO
    END DO

    !============================================
    ! PERIODIC CONTINUATION

    !CALL invoke_apply_bcs_uvtf(cufld, cvfld, hfld, Zfld)
    !call invoke(periodic_bc(cu), periodic_bc(cv), ....)

    !call invoke_apply_bcs_cu(ufield)
    ! Ultimately, this can be generated by PSyclone but in the
    ! absence of that we implement it manually here...
    ! First col = last col
    cufld(1,    1:N) = cufld(M+1,  1:N)
    ! Last row = first row
    cufld(1:M+1,N+1) = cufld(1:M+1,1)

    !call invoke_apply_bcs_ct(tfield)
    ! Last col = first col
    hfld(M+1,1:N) = hfld(1,  1:N)
    ! Last row = first row
    hfld(1:M+1,N+1) = hfld(1:M+1,1)

    !call invoke_apply_bcs_cv(vfield)
    cvfld(1:M,1    ) = cvfld(1:M,N+1)
    ! Last col = first col
    cvfld(M+1,1:N+1) = cvfld(1,  1:N+1)

    !call invoke_apply_bcs_cf(ffield)
    ! First col = last col
    zfld(1,    2:N+1) = zfld(M+1,  2:N+1)
    ! First row = last row
    zfld(1:M+1,1)     = zfld(1:M+1,N+1)

    !============================================
    ! COMPUTE NEW VALUES U,V AND P

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
    unew(1,    1:N) = unew(M+1,  1:N)
    ! Last row = first row
    unew(1:M+1,N+1) = unew(1:M+1,1)

    !call invoke_apply_bcs_cv(vnew)
    ! First row = last row
    vnew(1:M,1    ) = vnew(1:M,N+1)
    ! Last col = first col
    vnew(M+1,1:N+1) = vnew(1,  1:N+1)

    !call invoke_apply_bcs_ct(pnew)
    ! Last col = first col
    pnew(M+1,1:N) = pnew(1,  1:N)
    ! Last row = first row
    pnew(1:M+1,N+1) = pnew(1:M+1,1)

    !============================================
    ! The time-smoothing is applied to a field at *every* grid point
    
    ! Loop over 'columns'
    DO J=1,N+1 !idim2
      DO I=1,M+1 !idim1
        CALL time_smooth_code(i,j,ufld,unew,uold)
!      END DO
!    END DO

    ! Loop over 'columns'
!    DO J=1,N+1 ! idim2
!      DO I=1,M+1 ! idim1
         CALL time_smooth_code(i,j,vfld,vnew,vold)
!      END DO
!    END DO

    ! Loop over 'columns'
!    DO J=1,N+1 ! idim2
!      DO I=1,M+1 ! idim1
         CALL time_smooth_code(i,j,pfld,pnew,pold)
      END DO
    END DO

    !============================================
    ! Update for next step
    CALL copy_field(UNEW, Ufld)
    CALL copy_field(VNEW, Vfld)
    CALL copy_field(PNEW, Pfld)

  end subroutine invoke_time_step

  !===================================================

  !> Compute the mass flux in the x direction at point (i,j)
  subroutine compute_cu_code(i, j, cu, p, u)
    implicit none
    integer,  intent(in) :: I, J
    real(wp), intent(out), dimension(:,:) :: cu
    real(wp), intent(in),  dimension(:,:) :: p, u

    CU(I,J) = .5*(P(I,J)+P(I-1,J))*U(I,J)

  end subroutine compute_cu_code

  !===================================================

  !> Compute the mass flux in the y direction at point (i,j)
  subroutine compute_cv_code(i, j, cv, p, v)
    implicit none
    integer,  intent(in) :: I, J
    real(wp), intent(out), dimension(:,:) :: cv
    real(wp), intent(in),  dimension(:,:) :: p, v

    CV(I,J) = .5*(P(I,J)+P(I,J-1))*V(I,J)

  end subroutine compute_cv_code

  !===================================================

  !> Compute the potential vorticity on the grid point (i,j)
  subroutine compute_z_code(i, j, z, p, u, v)
    use mesh_mod, only: fsdx, fsdy
    implicit none
    integer,  intent(in) :: I, J
    real(wp), intent(out), dimension(:,:) :: z
    real(wp), intent(in),  dimension(:,:) :: p, u, v

    Z(I,J) =(FSDX*(V(I,J)-V(I-1,J))-FSDY*(U(I,J)-U(I,J-1)))/ &
                 (P(I-1,J-1)+P(I,J-1)+P(I,J)+P(I-1,J))

  end subroutine compute_z_code

  !===================================================

  SUBROUTINE compute_h_code(i, j, h, p, u, v)
    IMPLICIT none
    INTEGER, INTENT(in) :: I, J
    REAL(wp), INTENT(out), DIMENSION(:,:) :: h
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: p, u, v

    H(I,J) = P(I,J)+.25*(U(I+1,J)*U(I+1,J)+U(I,J)*U(I,J)     & 
                        +V(I,J+1)*V(I,J+1)+V(I,J)*V(I,J))

  END SUBROUTINE compute_h_code

  !===================================================

  SUBROUTINE compute_unew_code(i, j, unew, uold, z, cv, h, tdt)
    USE model_mod, ONLY: dx
    IMPLICIT none
    INTEGER, INTENT(in) :: I, J
    REAL(wp), INTENT(out), DIMENSION(:,:) :: unew
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: uold, z, cv, h
    REAL(wp), INTENT(in) :: tdt
    ! Locals
    REAL(wp) :: tdts8, tdtsdx
   
    !> These quantities are computed here because tdt is not
    !! constant. (It is == dt for first time step, 2xdt for
    !! all remaining time steps.)
    tdts8 = tdt/8.0d0
    tdtsdx = tdt/dx

    UNEW(I,J) = UOLD(I,J) +                                 &
                TDTS8*(Z(I,J+1)+Z(I,J)) *                   &
                (CV(I,J+1)+CV(I-1,J+1)+CV(I-1,J)+CV(I,J)) - &
                TDTSDX*(H(I,J)-H(I-1,J))

  END SUBROUTINE compute_unew_code

  !===================================================

  SUBROUTINE compute_vnew_code(i, j, vnew, vold, z, cu, h, tdt)
    USE model_mod, ONLY: dy
    IMPLICIT none
    INTEGER, INTENT(in) :: I, J
    REAL(wp), INTENT(out), DIMENSION(:,:) :: vnew
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: vold, z, cu, h
    REAL(wp), INTENT(in) :: tdt
    ! Locals
    REAL(wp) :: tdts8, tdtsdy
   
    !> These quantities are computed here because tdt is not
    !! constant. (It is == dt for first time step, 2xdt for
    !! all remaining time steps.)
    tdts8 = tdt/8.0d0
    tdtsdy = tdt/dy

    VNEW(I,J) = VOLD(I,J)-TDTS8*(Z(I+1,J)+Z(I,J))                 & 
                *(CU(I+1,J)+CU(I,J)+CU(I,J-1)+CU(I+1,J-1))        & 
                 -TDTSDY*(H(I,J)-H(I,J-1))

  END SUBROUTINE compute_vnew_code

  !===================================================

  SUBROUTINE compute_pnew_code(i, j, pnew, pold, cu, cv, tdt)
    USE model_mod, ONLY: dx, dy
    IMPLICIT none
    INTEGER, INTENT(in) :: I, J
    REAL(wp), INTENT(out), DIMENSION(:,:) :: pnew
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: pold, cu, cv
    REAL(wp), INTENT(in) :: tdt
    ! Locals
    REAL(wp) :: tdtsdx, tdtsdy
   
    !> These quantities are computed here because tdt is not
    !! constant. (It is == dt for first time step, 2xdt for
    !! all remaining time steps.)
    tdtsdx = tdt/dx
    tdtsdy = tdt/dy

    PNEW(I,J) = POLD(I,J)-TDTSDX*(CU(I+1,J)-CU(I,J))   & 
                         -TDTSDY*(CV(I,J+1)-CV(I,J))

  END SUBROUTINE compute_pnew_code

  !===================================================

  !> Kernel to smooth supplied field in time
  SUBROUTINE time_smooth_code(i, j, field, field_new, field_old)
    use time_smooth_mod, only: alpha
    IMPLICIT none
    INTEGER,  INTENT(in)                    :: i, j
    REAL(wp), INTENT(in),    DIMENSION(:,:) :: field
    REAL(wp), INTENT(in),    DIMENSION(:,:) :: field_new
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: field_old

    field_old(i,j) = field(i,j) + &
         alpha*(field_new(i,j) - 2.*field(i,j) + field_old(i,j))

  END SUBROUTINE time_smooth_code

  !===================================================

  SUBROUTINE copy_field(field_in, field_out)
    IMPLICIT none
    REAL(wp), INTENT(in),  DIMENSION(:,:) :: field_in
    REAL(wp), INTENT(out), DIMENSION(:,:) :: field_out
        
    field_out(:,:) = field_in(:,:)
        
  END SUBROUTINE copy_field

end module time_step_mod
