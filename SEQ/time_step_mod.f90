module time_step_mod
  use kind_params_mod
!  use topology_mod, only: M, N
  use timing_mod
  implicit none

contains

  subroutine invoke_time_step(cufld, cvfld, ufld, unew, uold, &
                              vfld, vnew, vold, &
                              pfld, pnew, pold, &
                              hfld, zfld, tdt,  &
                              m, n)
    use time_smooth_mod, only: alpha
    use mesh_mod, only: dx, dy, fsdx, fsdy ! Properties of the grid (4.0/d{x,y})
    implicit none
    real(wp), dimension(M+1,N+1), intent(inout) :: cufld, cvfld
    real(wp), dimension(M+1,N+1), intent(inout) :: unew, vnew, pnew
    real(wp), dimension(M+1,N+1), intent(inout) :: hfld, zfld, pfld, &
                                               ufld, vfld
    !> \todo Do we need these 'old' arrays as args or are they
    !! local workspace really?
    real(wp), dimension(M+1,N+1), intent(inout) :: uold, vold, pold
    real(wp),                     intent(in) :: tdt
    ! We pass in the array extents and loop bounds so that the compiler
    ! can 'see' that they are held constant throughout the time-stepping
    ! loop
    integer,  intent(in) :: m, n
    ! Locals
    integer :: I, J, idxt
    REAL(wp) :: tdts8, tdtsdy, tdtsdx

    !============================================
    ! COMPUTE CAPITAL U, CAPITAL V, Z AND H
    call timer_start('Compute CU,CV,CZ,H',idxt)
         
    ! M/N obtained from topology look-up
    do J= 1, N, 1
       do I = 1, M, 1

          !call compute_cu_code(i+1, j, cufld, pfld, ufld)
          CUfld(I+1,J) = .5*(Pfld(I+1,J)+Pfld(I,J))*Ufld(I+1,J)

          !call compute_cv_code(i, j+1, cvfld, pfld, vfld)
          CVfld(I,J+1) = .5*(Pfld(I,J+1)+Pfld(I,J))*Vfld(I,J+1)

          !call compute_z_code(i+1, j+1, zfld, pfld, ufld, vfld)
          Zfld(I+1,J+1) =(FSDX*(Vfld(I+1,J+1)-Vfld(I,J+1))- &
                          FSDY*(Ufld(I+1,J+1)-Ufld(I+1,J)))/ &
                   (pfld(I,J)+Pfld(I+1,J)+Pfld(I+1,J+1)+Pfld(I,J+1))

          !call compute_h_code(i, j, hfld, pfld, ufld, vfld)
          Hfld(I,J) = Pfld(I,J) + &
                      .25*(Ufld(I+1,J)*Ufld(I+1,J)+Ufld(I,J)*Ufld(I,J) & 
                       + Vfld(I,J+1)*Vfld(I,J+1)+Vfld(I,J)*Vfld(I,J))

       END DO
    END DO

    call timer_stop(idxt)

    !============================================
    ! PERIODIC CONTINUATION

    ! Ultimately, this can be generated by PSyclone but in the
    ! absence of that we implement it manually here...
    ! First col = last col
    do j = 1, N
       cufld(1,  j  ) = cufld(M+1,  j)
       zfld(1,   j+1) = zfld(M+1,  j+1)
       hfld(M+1, j  ) = hfld(1, j)
       cvfld(M+1,j+1) = cvfld(1, j+1)
    end do

    cufld(1:M+1,N+1) = cufld(1:M+1,1)

    hfld(1:M+1,N+1) = hfld(1:M+1,1)

    cvfld(1:M,1   ) = cvfld(1:M,N+1)
    cvfld(M+1, 1  ) = cvfld(1, 1)
    zfld(1:M+1,1)   = zfld(1:M+1,N+1)

    !============================================
    ! COMPUTE NEW VALUES U,V AND P
    call timer_start('Compute {U,V,P}NEW',idxt)

    tdtsdx = tdt/dx
    tdtsdy = tdt/dy
    tdts8 = tdt/8.0d0

    DO J=1, N, 1
       DO I= 1, M, 1

          !CALL compute_unew_code(i+1, j, unew, uold, &
          !                       zfld, cvfld, hfld, tdt)
          UNEW(I,J) = UOLD(I,J) +                                 &
                      TDTS8*(Zfld(I,J+1)+Zfld(I,J)) *                   &
                      (CVfld(I,J+1)+CVfld(I-1,J+1)+CVfld(I-1,J)+CVfld(I,J)) - &
                       TDTSDX*(Hfld(I,J)-Hfld(I-1,J))

          !CALL compute_vnew_code(i, j+1, vnew, vold, &
          !                       zfld, cufld, hfld, tdt)
          VNEW(I,J) = VOLD(I,J)-TDTS8*(Zfld(I+1,J)+Zfld(I,J))           &
                      *(CUfld(I+1,J)+CUfld(I,J)+CUfld(I,J-1)+CUfld(I+1,J-1)) &
                      -TDTSDY*(Hfld(I,J)-Hfld(I,J-1))

          !CALL compute_pnew_code(i, j, pnew, pold, &
          !                       cufld, cvfld, tdt)
          PNEW(I,J) = POLD(I,J)-TDTSDX*(CUfld(I+1,J)-CUfld(I,J))   & 
                      -TDTSDY*(CVfld(I,J+1)-CVfld(I,J))
       END DO
    END DO

    call timer_stop(idxt)

    !============================================
    ! PERIODIC CONTINUATION
    !CALL invoke_apply_bcs_uvt(UNEW, VNEW, PNEW)

    ! Ultimately, this can be generated by PSyclone but in the
    ! absence of that we implement it manually here...
    ! First col = last col
    do j = 1, N
       unew(1, j) = unew(M+1, j)
       vnew(M+1,j+1) = vnew(1, j+1)
       pnew(M+1,j) = pnew(1, j)
    end do

    unew(1:M+1,N+1) = unew(1:M+1,1)
    vnew(1:M,1    ) = vnew(1:M,N+1)
    vnew(M+1,1)     = vnew(1,1)
    pnew(1:M+1,N+1) = pnew(1:M+1,1)

    !============================================
    ! The time-smoothing is applied to a field at *every* grid point
    call timer_start('Time smooth',idxt)

    ! Loop over 'columns'
    DO J=1,N+1
      DO I=1,M+1
        uold(i,j) = ufld(i,j) + &
              alpha*(unew(i,j) - 2.*ufld(i,j) + uold(i,j))

        vold(i,j) = vfld(i,j) + &
              alpha*(vnew(i,j) - 2.*vfld(i,j) + vold(i,j))

        pold(i,j) = pfld(i,j) + &
              alpha*(pnew(i,j) - 2.*pfld(i,j) + pold(i,j))

        Ufld(I,J) = UNEW(I,J)
        Vfld(I,J) = VNEW(I,J)
        Pfld(I,J) = PNEW(I,J)
      END DO
    END DO

    call timer_stop(idxt)

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

end module time_step_mod
