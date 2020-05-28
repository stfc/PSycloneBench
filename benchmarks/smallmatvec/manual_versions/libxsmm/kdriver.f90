program kdriver
  use constants_mod, only : i_def, r_def, str_max_filename
  use matrix_vector_kernel_mod
  use omp_lib
  implicit none
  ! the arguments and array sizes needed to read/write the arrays.
  integer(kind=i_def) :: ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, ndf_any_space_2_x, undf_any_space_2_x
  integer(kind=i_def) :: ncell, nlayers, ncell_3d, ncolour
  integer(kind=i_def), allocatable :: ncp_colour(:), cmap(:,:)
  integer(kind=i_def), allocatable,target :: map_any_space_1_theta_adv_term(:,:)
  integer(kind=i_def), allocatable,target :: map_any_space_2_x(:,:)
  real(kind=r_def), allocatable :: theta_adv_term_data(:), x_data(:), ans_data(:)
  real(kind=r_def), allocatable :: ptheta_2_local_stencil(:,:,:)
  real(kind=r_def), allocatable :: nlayers_first(:,:,:,:)
  integer(kind=i_def) :: df, df2, lp, nthreads, memstart,memend
  real(kind=r_def) :: diff, err
  real(kind=r_def) :: start, end

  integer(kind=i_def), pointer :: map1(:) =>null(), map2(:)=>null()
  integer(kind=i_def) :: colour, cell, idx, tile, r, row
  integer(kind=i_def) :: c1_index, c2_index, c3_index, c4_index

  integer(kind=i_def) :: Nreps, nsize, el, i1, i2, i
  character(len=32) :: arg, arg2, initialization, method, traverse, colouringkind 
  logical :: new_data_layout

  real(kind=r_def), dimension(144) :: SFC12, SFC12i

  integer(kind = OMP_lock_kind), allocatable :: lock(:,:)

  traverse = "colouring" ! Default traversing
  ! cmdline argument parsing
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)

    select case (arg)
    case ('-g', '--generate')
        initialization="generate"
        if (i + 2 <= command_argument_count()) then
            call get_command_argument(i+1, arg2)
            read(arg2,*)  nsize
            call get_command_argument(i+2, arg2)
            read(arg2,*)  nlayers
            write(*,*) "Generating grid with size ", nsize, "and ", nlayers," layers"
        else
            write(*,*) "Error parsing -gen [size][layers] parameter"
            call exit(-1)
        endif
    end select
  end do

  if (.TRUE.) then
      ! Generate scalars
      ncell = nsize * nsize
      ncolour = 4
      ndf_any_space_1_theta_adv_term = 8
      undf_any_space_1_theta_adv_term = (nsize+1)*(nsize+1)*(nlayers+1)
      ndf_any_space_2_x = 6
      undf_any_space_2_x = (2*(nsize+1)*(nsize)*nlayers) + (nsize*nsize*(nlayers+1))
      ncell_3d = nlayers * ncell

      ! Allocate arrays
      allocate( map_any_space_1_theta_adv_term(ndf_any_space_1_theta_adv_term,ncell))
      allocate( map_any_space_2_x(ndf_any_space_2_x,ncell) )
      allocate( theta_adv_term_data(undf_any_space_1_theta_adv_term) )
      allocate( x_data(undf_any_space_2_x) )
      allocate( ptheta_2_local_stencil( ndf_any_space_1_theta_adv_term, ndf_any_space_2_x, ncell_3d) )
      allocate( nlayers_first(nlayers, ndf_any_space_2_x, ndf_any_space_1_theta_adv_term, ncell) )
      allocate( ans_data(undf_any_space_1_theta_adv_term) )
      allocate( ncp_colour(ncolour) )

      ! Populate arrays
      do i2 = 0, nsize - 1
        do i1 = 0, nsize - 1
          ! Populate  map_any_space_1_theta_adv_term (8 vertices)
          do el = 1, 2
            map_any_space_1_theta_adv_term( 1 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = (i2 * (nsize+1) + i1) * (nlayers + 1) + el
            map_any_space_1_theta_adv_term( 2 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = (i2 * (nsize+1) + i1+1) * (nlayers + 1) + el
            map_any_space_1_theta_adv_term( 3 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = ((i2+1) * (nsize+1) + i1+1) * (nlayers + 1) + el
            map_any_space_1_theta_adv_term( 4 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = ((i2+1) * (nsize+1) + i1) * (nlayers + 1) + el
          enddo

          !Populate  map_any_space_2_x (6 surfaces)
          map_any_space_2_x(1, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + 1
          map_any_space_2_x(2, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + nlayers + 1
          if (i1.eq.nsize-1) then
            map_any_space_2_x(3, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + (i1+1) * (3*nlayers + 1) + 1
          else
            map_any_space_2_x(3, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + (i1+1) * (3*nlayers + 1) + nlayers + 1
          endif
          if (i2.eq.nsize-1) then
            map_any_space_2_x(4, (i2 * nsize) + i1 + 1 ) = (i2+1) * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * nlayers + 1
          else
            map_any_space_2_x(4, (i2 * nsize) + i1 + 1 ) = (i2+1) * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + 1
          endif
          map_any_space_2_x(5, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + 2*nlayers + 1
          map_any_space_2_x(6, (i2 * nsize) + i1 + 1 ) = i2 * ((3 * nsize * nlayers) + nsize + nlayers) + i1 * (3*nlayers + 1) + 2*nlayers + 2
        enddo
      enddo

      
      if (traverse.eq."colouring") then
        ! Number of cells per color = num times per row * num times per column
        ncp_colour(1) = ((nsize/2)+mod(nsize,2)) * ((nsize/2)+mod(nsize,2))
        ncp_colour(2) = ( nsize/2 )              * ((nsize/2)+mod(nsize,2))
        ncp_colour(3) = ((nsize/2)+mod(nsize,2)) * (nsize/2)
        ncp_colour(4) = ( nsize/2 )              * (nsize/2)

        allocate(cmap(ncolour,maxval(ncp_colour))) !color map (colour, cell number)

        c1_index = 1
        c2_index = 1
        c3_index = 1
        c4_index = 1
        do i2 = 0, nsize - 1
            do i1 = 0, nsize - 1
                if (mod(i1,2).eq.0) then
                    if (mod(i2,2).eq.0) then
                        cmap(1,c1_index) = (i2 * nsize) + i1 + 1
                        c1_index = c1_index + 1
                    else
                        cmap(2,c2_index) = (i2 * nsize) + i1 + 1
                        c2_index = c2_index + 1
                    endif
                else
                    if (mod(i2,2).eq.0) then
                        cmap(3,c3_index) = (i2 * nsize) + i1 + 1
                        c3_index = c3_index + 1
                    else
                        cmap(4,c4_index) = (i2 * nsize) + i1 + 1
                        c4_index = c4_index + 1
                    endif
                endif
            enddo
        enddo
      end if
 
      ! Random data
      call random_number(theta_adv_term_data)
      call random_number(x_data)
      call random_number(ptheta_2_local_stencil)
      theta_adv_term_data = theta_adv_term_data * 100000
      x_data = x_data * 100000
      ptheta_2_local_stencil = ptheta_2_local_stencil * 100000
      ans_data =1
  endif

  ! Starting computation phase
  write(*,*) "Starting mv with:"
  write(*,*) "Colors=", ncolour, " -> ", ncp_colour
  write(*,*) "Layers=", nlayers
  write(*,*) "Horitzontal Cells=", ncell_3d/nlayers, ncell
  write(*,*) "3D Cells=", ncell_3d
  write(*,*) "ndf1=", ndf_any_space_1_theta_adv_term, ", undf1=", undf_any_space_1_theta_adv_term
  write(*,*) "ndf2=", ndf_any_space_2_x, ", undf2=", undf_any_space_2_x


  start =  omp_get_wtime()
      do lp = 1, 10000
        do colour = 1, ncolour
          !$omp parallel do default(shared), private(cell, map1, map2), schedule(runtime)
          do cell = 1, ncp_colour(colour)
             map1=>map_any_space_1_theta_adv_term(:,cmap(colour,cell))
             map2=>map_any_space_2_x(:,cmap(colour,cell))
             call matrix_vector_kernel_code(cmap(colour,cell), nlayers, &
                 theta_adv_term_data, x_data, ncell_3d, ptheta_2_local_stencil, &
                 ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                 map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
          end do
          !$omp end parallel do
        end do
      end do

  end = omp_get_wtime()
  write(*,*) "T=",nthreads, "Loop time=",end - start, " s"

  if (.TRUE.) then
      ! print out reduction sum
      write(*,*) "Reduction value:", sum(theta_adv_term_data)
  endif

  deallocate( map_any_space_1_theta_adv_term,  map_any_space_2_x )
  deallocate( theta_adv_term_data, & 
       x_data,                           &
       ptheta_2_local_stencil )
end program kdriver
