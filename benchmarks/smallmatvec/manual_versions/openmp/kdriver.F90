#ifndef VERSION
#   define VERSION original
#endif

#define PASTER(v) matrix_vector_code_ ## v
#define MVFUNCTION(v) PASTER(v)

#define STR_(X) #X
#define STR(X) STR_(X)


program kdriver
  use constants_mod, only : i_def, r_def, str_max_filename
  use proflib_io_mod, only : dino_type
  use matrix_vector_kernel_mod
  use utils
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
  type(dino_type) :: dino

  integer(kind=i_def) :: Nreps, nsize, el, i1, i2, i
  character(len=32) :: arg, arg2, initialization, method, traverse, colouringkind 
  character(len=32) :: dotoutput
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
    case('-Nreps')
        if (i + 2 <= command_argument_count()) then
            call get_command_argument(i+1, arg2)
            read(arg2,*)  Nreps
            write(*,*) "Executing with", Nreps, "repetitions of the inner loop"
        endif
    case('-t','--traverse')
        if (i + 1 <= command_argument_count()) then
            call get_command_argument(i+1, arg2)
            select case(arg2)
            case ("linear")
                traverse = "linear"
            case ("omp-locking")
                traverse = "omp-locking"
            case ("ompall")
                traverse = "ompall"
            case ("colouring")
                traverse = "colouring"
            case ("colouring2")
                traverse = "colouring2"
            case ("colouring-rows")
                traverse = "colouring-rows"
            endselect
        endif
    case('-dg','--dotgraph')
        dotoutput = "write"
    case('-h','--help')
        write(*,*) "LFRic Matrix-Vector Multiplication isolated Kernel"
    end select
  end do


  call system_mem_usage(memstart)
  !$omp parallel
  nthreads=omp_get_num_threads()
  !$omp end parallel
  start =  omp_get_wtime()
  dino = dino_type()

  if (.not.(initialization.eq."generate")) then
      ! Get parameters
      call dino%input_scalar(ndf_any_space_1_theta_adv_term)
      call dino%input_scalar(undf_any_space_1_theta_adv_term)
      call dino%input_scalar(ndf_any_space_2_x) 
      call dino%input_scalar(undf_any_space_2_x)
      call dino%input_scalar(nlayers)
      call dino%input_scalar(ncell_3d)
      call dino%input_scalar(ncell)
      call dino%input_scalar(ncolour)

      ! Allocate and populate arrays
      allocate( map_any_space_1_theta_adv_term(ndf_any_space_1_theta_adv_term,ncell))
      allocate( map_any_space_2_x(ndf_any_space_2_x,ncell) )
      allocate( theta_adv_term_data(undf_any_space_1_theta_adv_term))
      allocate( x_data(undf_any_space_2_x) )
      allocate(ptheta_2_local_stencil( ndf_any_space_1_theta_adv_term, ndf_any_space_2_x, ncell_3d) )
      allocate(nlayers_first(nlayers, ndf_any_space_2_x, ndf_any_space_1_theta_adv_term, ncell))
      allocate(ans_data(undf_any_space_1_theta_adv_term))
      allocate(ncp_colour(ncolour))

      call dino%input_array(map_any_space_1_theta_adv_term,ndf_any_space_1_theta_adv_term,ncell)
      call dino%input_array(map_any_space_2_x,ndf_any_space_2_x,ncell)
      call dino%input_array(ncp_colour,ncolour)
      allocate(cmap(ncolour,maxval(ncp_colour)))
      call dino%input_array(cmap,ncolour,maxval(ncp_colour))
      call dino%input_array(theta_adv_term_data,undf_any_space_1_theta_adv_term)
      call dino%input_array(x_data,undf_any_space_2_x)
      call dino%input_array(ptheta_2_local_stencil, ndf_any_space_1_theta_adv_term, ndf_any_space_2_x, ncell_3d)
      call dino%input_array(ans_data,undf_any_space_1_theta_adv_term)
      call dino%io_close()
  else
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

      allocate(lock(nsize+1, nsize+1))

      do i1 = 1, nsize+1
        do i2 = 1, nsize+1
          call OMP_init_lock(lock(i1,i2))
        enddo
      enddo

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
      else if (traverse.eq."colouring2") then
        if (.not.(IAND (nsize, nsize-1) == 0)) then
          write(*,*) "Error colouring2 strategy just allowed for nsize is a power of 2"
          call exit(0)
        endif
      ! Number of cells per color = num times per row * num times per column
        ncp_colour(1) = ((nsize/2)+mod(nsize,2)) * ((nsize/2)+mod(nsize,2)) / 2
        ncp_colour(2) = ( nsize/2 )              * ((nsize/2)+mod(nsize,2)) / 2
        ncp_colour(3) = ((nsize/2)+mod(nsize,2)) * (nsize/2) / 2
        ncp_colour(4) = ( nsize/2 )              * (nsize/2) / 2
        write(*,*) ncp_colour

        allocate(cmap(ncolour,maxval(ncp_colour))) !color map (colour, cell number)

        c1_index = 1
        c2_index = 1
        c3_index = 1
        c4_index = 1
        do i2 = 0, nsize - 1
            do i1 = 0, nsize - 1
                if (mod(i1,4).eq.0) then
                    if (mod(i2,2).eq.0) then
                        cmap(1,c1_index) = (i2 * nsize) + i1 + 1
                        c1_index = c1_index + 1
                    else
                        cmap(2,c2_index) = (i2 * nsize) + i1 + 1
                        c2_index = c2_index + 1
                    endif
                else if (mod(i1,2).eq.0) then
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

      else if (traverse.eq."colouring-rows") then
          ! In this version ncp_colour contains the number of rows per color and
          ! cmap the first cell of each coloured row

        ncp_colour(1) = (nsize / 2) + mod(nsize,2)
        ncp_colour(2) = (nsize / 2)
        ncp_colour(3) = 0
        ncp_colour(4) = 0
        allocate(cmap(ncolour,maxval(ncp_colour))) !color map (colour, cell number)

        c1_index = 1
        c2_index = 1
            do i1 = 0, nsize - 1
                if (mod(i1,2).eq.0) then
                        cmap(1,c1_index) = (i1 * nsize) + 1
                        c1_index = c1_index + 1
                else
                        cmap(2,c2_index) = (i1 * nsize) + 1
                        c2_index = c2_index + 1
                endif
            enddo

      else
          write(*,*) "Error: No colouring strategy selected"
      endif
      
      ! Random data
      call random_number(theta_adv_term_data)
      call random_number(x_data)
      call random_number(ptheta_2_local_stencil)
      theta_adv_term_data = theta_adv_term_data * 100000
      x_data = x_data * 100000
      ptheta_2_local_stencil = ptheta_2_local_stencil * 100000
      if (.false.) then
              do i1 = 1, undf_any_space_1_theta_adv_term
                theta_adv_term_data(i1) = 1
              enddo
              do i1=1, undf_any_space_2_x
                x_data(i1) = 1
              enddo
              ptheta_2_local_stencil =1
      endif
      ans_data =1
  endif

  call system_mem_usage(memend)

  end = omp_get_wtime()
  write(*,*) "Allocation and initialization time =" , end - start, " s"
  write(*,*) "Datastructure RSS memory allocated =" , memend - memstart, " KB"

  ! Starting computation phase
  write(*,*) "Starting mv with:"
  write(*,*) "Colors=", ncolour, " -> ", ncp_colour
  write(*,*) "Layers=", nlayers
  write(*,*) "Horitzontal Cells=", ncell_3d/nlayers, ncell
  write(*,*) "3D Cells=", ncell_3d
  write(*,*) "ndf1=", ndf_any_space_1_theta_adv_term, ", undf1=", undf_any_space_1_theta_adv_term
  write(*,*) "ndf2=", ndf_any_space_2_x, ", undf2=", undf_any_space_2_x


  if (dotoutput.eq."write")then
      open(1234, file="2dgrid.dot", action="write")
      write(1234,*) "graph graphme{"
      do cell = 1, ncp_colour(1)
        write(1234,"(A,I3.3,A,I3.3,A,I1.1,A)") "N", cmap(1,cell), " [label=""", cmap(1,cell), """,color=blue];"
      enddo
      do cell = 1, ncp_colour(2)
        write(1234,"(A,I3.3,A,I3.3,A,I1.1,A)") "N", cmap(2,cell), " [label=""", cmap(2,cell), """,color=green];"
      enddo
      do cell = 1, ncp_colour(3)
        write(1234,"(A,I3.3,A,I3.3,A,I1.1,A)") "N", cmap(3,cell), " [label=""", cmap(3,cell), """,color=red];"
      enddo
      do cell = 1, ncp_colour(4)
        write(1234,"(A,I3.3,A,I3.3,A,I1.1,A)") "N", cmap(4,cell), " [label=""", cmap(4,cell), """,color=yellow];"
      enddo
      !do cell = 1, ncp_colour(5)
      !  write(1234,"(A,I3.3,A,I3.3,A,I1.1,A)") "N", cmap(5,cell), " [label=""", cmap(5,cell), """,color=purple];"
      !enddo
      !do cell = 1, ncp_colour(6)
      !  write(1234,"(A,I3.3,A,I3.3,A,I1.1,A)") "N", cmap(6,cell), " [label=""", cmap(6,cell), """,color=brown];"
      !enddo
      do colour = 1, ncolour
        do cell = 1, ncp_colour(colour)
            map1=>map_any_space_1_theta_adv_term(:,cmap(colour,cell))
            map2=>map_any_space_2_x(:,cmap(colour,cell))
            do lp = 1, ndf_any_space_1_theta_adv_term
                write(1234,"(A,I3.3,A,I4.4,A)") "N", cmap(colour,cell), " -- L", map1(lp)/nlayers, ";"
            enddo
            !do lp = 1, ndf_any_space_2_x
            !    write(1234,"(A,I3.3,A,I4.4,A)") "N", cmap(colour,cell), " -- R", map2(lp)/nlayers, ";"
            !enddo

        enddo
      enddo
      
      write(1234,*) "}"
      close(1234)
  endif

  if ((STR(VERSION) .eq. "nlayersf2") .or. &
     (STR(VERSION) .eq. "nlayersf") .or. &
     (STR(VERSION) .eq. "nlayersf_atomics") .or. &
     (STR(VERSION) .eq. "nlayersf_moreops") .or. &
     (STR(VERSION) .eq. "nlayersf_split")) then
    new_data_layout = .true.
  else
     new_data_layout = .false.
  endif

  write(*,*) " ----   ",STR(VERSION),"     ----- "

  if (new_data_layout) then
      start =  omp_get_wtime()
      do cell = 0, ncell-1
        do lp = 1, nlayers
            do df =  1, ndf_any_space_1_theta_adv_term
                do df2 =  1, ndf_any_space_2_x
                    nlayers_first(lp,df2,df,cell+1) = ptheta_2_local_stencil(df,df2,(cell*nlayers)+lp)
                end do
            end do
        end do
      end do
      end = omp_get_wtime()
      write(*,*) "Convert datastructure time=", end - start, " s"
  endif

  start =  omp_get_wtime()
  if (traverse.eq."linear") then  
      write(*,*) "Lineal traversing Version"
      do lp = 1, 1000
        do cell = 1, ncell
          map1=>map_any_space_1_theta_adv_term(:,cell)
          map2=>map_any_space_2_x(:,cell)
          if (new_data_layout) then
             call MVFUNCTION(VERSION)(cell, nlayers, theta_adv_term_data, &
             x_data, ncell_3d, nlayers_first, ndf_any_space_1_theta_adv_term, &
             undf_any_space_1_theta_adv_term,map1, &
             ndf_any_space_2_x, undf_any_space_2_x, map2 )
          else
            call MVFUNCTION(VERSION)(cell, nlayers, &
            theta_adv_term_data, x_data, ncell_3d, ptheta_2_local_stencil, &
            ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
            map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
          endif
        enddo
      enddo

  else if (traverse.eq."ompall") then  
      write(*,*) "Single OMP region (no colours)"
      do lp = 1, 1000
        !$omp parallel do default(shared), private(cell, map1, map2), schedule(runtime)
        do cell = 1, ncell
              map1=>map_any_space_1_theta_adv_term(:,cell)
              map2=>map_any_space_2_x(:,cell)
              ! Just implemented with nlayersf for now
              !call MVFUNCTION(VERSION)_atomics(cell, nlayers,  &
              !   theta_adv_term_data, x_data, ncell_3d, &
              !   nlayers_first, ndf_any_space_1_theta_adv_term, &
              !   undf_any_space_1_theta_adv_term, &
              !   map1, &
              !   ndf_any_space_2_x, undf_any_space_2_x,&
              !   map2 )
        enddo
      enddo

  else if (traverse.eq."tiling") then  
      write(*,*) "Tiling traversing Version"
      do lp = 1, 1000
        do tile = 1, ncell, 6
          do cell = tile, tile + 6
            map1=>map_any_space_1_theta_adv_term(:,cell)
            map2=>map_any_space_2_x(:,cell)
            if (new_data_layout) then
                call MVFUNCTION(VERSION)(cell, nlayers, theta_adv_term_data, &
                x_data, ncell_3d, nlayers_first, ndf_any_space_1_theta_adv_term, &
                undf_any_space_1_theta_adv_term, map1, ndf_any_space_2_x, &
                undf_any_space_2_x, map2 )
            else
                call MVFUNCTION(VERSION)(cell, nlayers, &
                theta_adv_term_data, x_data, ncell_3d, ptheta_2_local_stencil, &
                ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
            endif
          enddo
        enddo
      enddo
  else if (traverse.eq."SFC12") then  
      write(*,*) "Space Filling Curve Version"
  SFC12 = (/ 1, 2, 14, 13, 25, 37, 38, 26, 27, 39, 40, 28, 16, 15, 3, 4, &
  5, 17,18,6,7,8,20,19,31,32,44,43,42,30,29,41, &
  53,65,66,54,55,56,68,67,79,80,92,91,90,78,77,89, &
  88,87,75,76,64,52,51,63,62,50,49,61,73,74,86,85, &
  97,109,110,98,99,100,112,111,123,124,136,135,134,122,121,133, &
  140,128,127,139,138,137,125,126,114,113,101,102,103,115,116,104, &
  105,117,118,106,107,108,120,119,131,132,144,143,142,130,129,141, &
  96,84,83,95,94,93,81,82,70,69,57,58,59,71,72,60, &
  48,36,35,47,46,45,33,34,22,21,9,10,11,23,24,12 /)

  SFC12i = SFC12 - 1


      do lp = 1, 1000
        do tile = 1, ncell, 144
            do idx = 1, 144
                cell = tile + SFC12i(idx)
                map1=>map_any_space_1_theta_adv_term(:,cell)
                map2=>map_any_space_2_x(:,cell)
                if (new_data_layout) then
                    call MVFUNCTION(VERSION)(cell, nlayers, theta_adv_term_data, &
                    x_data, ncell_3d, nlayers_first, ndf_any_space_1_theta_adv_term, &
                    undf_any_space_1_theta_adv_term, map1, &
                    ndf_any_space_2_x, undf_any_space_2_x, map2 )
                else
                    call MVFUNCTION(VERSION)(cell, nlayers, &
                    theta_adv_term_data, x_data, ncell_3d, ptheta_2_local_stencil, &
                    ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                    map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
                endif
            enddo
        enddo
      enddo

  else if (traverse.eq."colouring") then
      write(*,*) "Coloring Version"
      do lp = 1, 1000
        do colour = 1, ncolour
          !$omp parallel do default(shared), private(cell, map1, map2), schedule(runtime)
          do cell = 1, ncp_colour(colour)
             map1=>map_any_space_1_theta_adv_term(:,cmap(colour,cell))
             map2=>map_any_space_2_x(:,cmap(colour,cell))
             ! All input parameters but 'theta_adv_term_data'
             if (new_data_layout) then
                call MVFUNCTION(VERSION)(cmap(colour,cell), nlayers, &
                    theta_adv_term_data, x_data, ncell_3d, nlayers_first, &
                    ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                    map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
             else
                call MVFUNCTION(VERSION)(cmap(colour,cell), nlayers, &
                    theta_adv_term_data, x_data, ncell_3d, ptheta_2_local_stencil, &
                    ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                    map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
              endif
          end do
          !$omp end parallel do
        end do
      end do

  else if (traverse.eq."colouring2") then
      write(*,*) "Coloring Version join 2 cells"
      do lp = 1, 1000
        do colour = 1, ncolour
          !$omp parallel do default(shared), private(cell, map1, map2), schedule(runtime)
          do cell = 1, ncp_colour(colour)
             !write(*,*) colour, cell, cmap(colour, cell)
             map1=>map_any_space_1_theta_adv_term(:,cmap(colour,cell))
             map2=>map_any_space_2_x(:,cmap(colour,cell))
             if (new_data_layout) then
                call MVFUNCTION(VERSION)(cmap(colour,cell), nlayers, &
                theta_adv_term_data, x_data, ncell_3d, nlayers_first, &
                ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
             else
                call MVFUNCTION(VERSION)(cmap(colour,cell), nlayers, &
                    theta_adv_term_data, x_data, ncell_3d, ptheta_2_local_stencil, &
                    ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                    map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
             endif
             map1=>map_any_space_1_theta_adv_term(:,cmap(colour,cell)+1)
             map2=>map_any_space_2_x(:,cmap(colour,cell)+1)
             if (new_data_layout) then
                call MVFUNCTION(VERSION)(cmap(colour,cell), nlayers, &
                theta_adv_term_data, x_data, ncell_3d, nlayers_first, &
                ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
             else
                call MVFUNCTION(VERSION)(cmap(colour,cell), nlayers, &
                    theta_adv_term_data, x_data, ncell_3d, ptheta_2_local_stencil, &
                    ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                    map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
             endif
          end do
          !$omp end parallel do
        end do
      end do

  else if (traverse.eq."colouring-rows") then
      write(*,*) "Coloring Version Rows"
      do lp = 1, 1000
        do colour = 1, 2
          !$omp parallel do default(shared), private(i1, row, map1, map2), schedule(runtime)
          do row = 1, ncp_colour(colour) ! for each row of this colour

            do i1= 0, nsize-1 ! For each element in the row
              !write(*,*) colour, row, cmap(colour, row) + i1
              map1=>map_any_space_1_theta_adv_term(:,cmap(colour,row) + i1 )
              map2=>map_any_space_2_x(:,cmap(colour,row) + i1)

              ! All input parameters but 'theta_adv_term_data'
              if (new_data_layout) then
                call MVFUNCTION(VERSION)(cmap(colour,row) + i1, nlayers, &
                theta_adv_term_data, x_data, ncell_3d, nlayers_first, &
                ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
              else
                call MVFUNCTION(VERSION)(cmap(colour,row) + i1, nlayers, &
                    theta_adv_term_data, x_data, ncell_3d, ptheta_2_local_stencil, &
                    ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                    map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
              endif
            enddo
          end do
          !$omp end parallel do
        end do
        !call exit(0)
      end do

   else if (traverse.eq."omp-locking") then
      write(*,*) "OMP locking traversing"
      do lp = 1, 1000
        !$omp parallel do default(shared), private(i1,i2, cell, map1, map2), schedule(runtime)
        do cell = 1, ncell
            !do r =1, Nreps
            i1 = ((cell-1) /nsize) + 1
            i2 = mod((cell-1), nsize) + 1
            !write(*,*) "Cell: ", cell, i1, i2
            call OMP_set_lock(lock(i1,i2))
            call OMP_set_lock(lock(i1+1,i2))
            call OMP_set_lock(lock(i1,i2+1))
            call OMP_set_lock(lock(i1+1,i2+1))
            map1=>map_any_space_1_theta_adv_term(:,cell)
            map2=>map_any_space_2_x(:,cell)
            if (new_data_layout) then
                call MVFUNCTION(VERSION)(cell, nlayers, theta_adv_term_data, &
                x_data, ncell_3d, nlayers_first, ndf_any_space_1_theta_adv_term, &
                undf_any_space_1_theta_adv_term, map1, ndf_any_space_2_x, &
                undf_any_space_2_x, map2 )
            else
                call MVFUNCTION(VERSION)(cell, nlayers, &
                    theta_adv_term_data, x_data, ncell_3d, ptheta_2_local_stencil, &
                    ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
                    map1, ndf_any_space_2_x, undf_any_space_2_x, map2 )
            endif
            call OMP_unset_lock(lock(i1,i2))
            call OMP_unset_lock(lock(i1+1,i2))
            call OMP_unset_lock(lock(i1,i2+1))
            call OMP_unset_lock(lock(i1+1,i2+1))
            !enddo
        enddo
      enddo

  else
    write(*,*) "No traversing method specified"
    call exit(-1)
  endif

  end = omp_get_wtime()
  write(*,*) "T=",nthreads, "Loop time=",end - start, " s"


  if (new_data_layout) then
      start =  omp_get_wtime()
      do cell = 0, ncell-1
        do lp = 1, nlayers
            do df =  1, ndf_any_space_1_theta_adv_term
                do df2 =  1, ndf_any_space_2_x
                    ptheta_2_local_stencil(df,df2,(cell*nlayers)+lp) = nlayers_first(lp,df2,df,cell+1)
                end do
            end do
        end do
      end do
      end = omp_get_wtime()
      write(*,*) "Convert back datastructure time=", end - start, " s"
  endif




  if (.not.(initialization.eq."generate")) then
      write(*,*) "dino_dump:Checking the answer ..."
      err = 0.01
      do df=1, undf_any_space_1_theta_adv_term
         diff=(theta_adv_term_data( df )-ans_data(df))/(theta_adv_term_data( df )+ans_data(df))
         if(abs(diff)>=err) then
            write(*,*) df,theta_adv_term_data( df ), ans_data(df), diff
            call exit(0)      
         end if
      end do
    else
        ! print out reduction sum
        write(*,*) "Reduction value:", sum(theta_adv_term_data)
    endif

  ! DUMP SOLUTION TO FILE, USE CAREFULLY
  ! call dino%output_array(theta_adv_term_data, undf_any_space_1_theta_adv_term)
  ! call dino%io_close()
  
  deallocate( map_any_space_1_theta_adv_term,  map_any_space_2_x )
  deallocate( theta_adv_term_data, & 
       x_data,                           &
       ptheta_2_local_stencil )
end program kdriver
