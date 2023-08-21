program kdriver
    use psy_layer, only: run_psy_layer

    implicit none

    ! Constants
    integer, parameter :: i_def        = 4       !< Default integer kind for application.
    integer, parameter :: r_def        = 8      !< Default real kind for application.

    call main()
contains
    
subroutine main()
    ! the arguments and array sizes needed to read/write the arrays.
    integer(kind=i_def) :: ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term
    integer(kind=i_def) :: ndf_any_space_2_x, undf_any_space_2_x
    integer(kind=i_def) :: ncell, nlayers, ncell_3d, ncolour
    integer(kind=i_def) :: ncp_colour(4)
    integer(kind=i_def), allocatable :: cmap(:,:)
    integer(kind=i_def), allocatable,target :: map_any_space_1_theta_adv_term(:,:)
    integer(kind=i_def), allocatable,target :: map_any_space_2_x(:,:)
    real(kind=r_def), allocatable :: theta_adv_term_data(:), x_data(:)
    real(kind=r_def), allocatable :: ptheta_2_local_stencil(:,:,:)
    integer(kind=i_def) :: memstart, memend, memmaps, memmatrix, memvectors, memcmap, memtotal
    real(kind=r_def) :: start, end, totaltime
    integer(kind=i_def) :: niters, nsize, i
    character(len=32) :: arg, arg2, traverse 
    ! real(kind=r_def), allocatable :: nlayers_first(:,:,:,:)
    ! logical :: new_data_layout = .false.

    ! Default values
    nsize = 64
    nlayers = 64
    niters = 100
    traverse = "linear"
    ! traverse = "colouring"
  
    ! cmdline argument parsing
    do i = 1, command_argument_count(), 2
        call get_command_argument(i, arg)
    
        select case (arg)
        case ('-s', '--nsize')
            if (i + 1 <= command_argument_count()) then
                call get_command_argument(i+1, arg2)
                read(arg2,*)  nsize
            else
                call print_help_and_exit()
            endif
        case ('-l', '--nlayers')
            if (i + 1 <= command_argument_count()) then
                call get_command_argument(i+1, arg2)
                read(arg2,*)  nlayers
            else
                call print_help_and_exit()
            endif
        case('-n', '--niters')
            if (i + 1 <= command_argument_count()) then
                call get_command_argument(i+1, arg2)
                read(arg2,*)  niters
            else
                call print_help_and_exit()
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
            else
                call print_help_and_exit()
            endif
        case('-h','--help')
            call print_help_and_exit()
        case DEFAULT
            call print_help_and_exit()
        end select
    end do

    ! Set up domain scalars
    ncell = nsize * nsize
    ncell_3d = nlayers * ncell
    ncolour = 4

    ! Set up fields sizes
    ndf_any_space_1_theta_adv_term = 8
    undf_any_space_1_theta_adv_term = (nsize+1)*(nsize+1)*(nlayers+1)
    ndf_any_space_2_x = 6
    undf_any_space_2_x = (2*(nsize+1)*(nsize)*nlayers) + (nsize*nsize*(nlayers+1)) ! FIXME: This looks wrong

    ! Print benchmark info
    write(*,*) "LFRic Matrix-Vector Multiplication (lhs = matrix * x) with:"
    write(*,*) " - Horitzonal size:                   ", ncell, " column-cells"
    write(*,*) " - Vertical size:                     ", nlayers," layers"
    write(*,*) " - Total 3d cells:                    ", ncell_3d," cells"
    write(*,*) " - lhs (vertices: 8 points per cell): ", undf_any_space_1_theta_adv_term, " points"
    write(*,*) " - x (faces: 6 points per cell):      ", undf_any_space_2_x, " points"
    write(*,*) " - matrix (8x6 per cell):             ", ndf_any_space_1_theta_adv_term * ndf_any_space_2_x * ncell_3d, " points"
    write(*,*) " - Time-steps:                        ", niters, " iterations"
    write(*,*) " - Each step does:                    ", ncell_3d, " small MV (8x6) * (6) multiplications"

    call cpu_time(start)
    call system_mem_usage(memstart)
    ! call random_init(.true., .true.) ! make random_number generation repeatable

    ! Allocate arrays
    allocate( map_any_space_1_theta_adv_term(ndf_any_space_1_theta_adv_term,ncell) ) 
    allocate( map_any_space_2_x(ndf_any_space_2_x,ncell) )
    call populate_maps(nsize, nlayers, map_any_space_1_theta_adv_term, map_any_space_2_x)
    call system_mem_usage(memend)
    memmaps = memend - memstart
    
    allocate( theta_adv_term_data(undf_any_space_1_theta_adv_term) )
    allocate( x_data(undf_any_space_2_x) )
    call ascending_init(theta_adv_term_data)
    call ascending_init(x_data)
    ! call random_number(theta_adv_term_data)
    ! call random_number(x_data)
    call system_mem_usage(memend)
    memvectors = memend - memmaps

    allocate( ptheta_2_local_stencil( ndf_any_space_1_theta_adv_term, ndf_any_space_2_x, ncell_3d) )
    !allocate( nlayers_first(nlayers, ndf_any_space_2_x, ndf_any_space_1_theta_adv_term, ncell) )
    call ascending_init_matrix(ptheta_2_local_stencil)
    ! call random_number(ptheta_2_local_stencil)
    call system_mem_usage(memend)
    memmatrix = memend - memvectors

    if (traverse.eq."colouring") then
        ! Number of cells per color = num times per row * num times per column
        ncp_colour(1) = ((nsize/2)+mod(nsize,2)) * ((nsize/2)+mod(nsize,2))
        ncp_colour(2) = ( nsize/2 )              * ((nsize/2)+mod(nsize,2))
        ncp_colour(3) = ((nsize/2)+mod(nsize,2)) * (nsize/2)
        ncp_colour(4) = ( nsize/2 )              * (nsize/2)

        allocate(cmap(4,maxval(ncp_colour))) !color map (colour, cell number)
        call populate_strided_cmap(nsize, cmap)
    endif
    call system_mem_usage(memend)
    memcmap = memend - memmatrix
    memtotal = memend - memstart
    call cpu_time(end)
    write(*,*) " - Allocation and init time:          " , ceiling((end - start) * 1000), " ms"
    write(*,*) " - Total memory allocated:            " , memtotal/1000, " MB  (", &
        memmaps/1000, "/", memvectors/1000, "/", memmatrix/1000, "/", memcmap/1000, ")"
    write(*,*) ""

    ! Starting computation phase
    call cpu_time(start)
    call run_psy_layer( &
        ! benchmark parameters
        traverse, niters, ncell, nlayers, ncell_3d, &
        ! lhs
        theta_adv_term_data, map_any_space_1_theta_adv_term, &
        ndf_any_space_1_theta_adv_term, undf_any_space_1_theta_adv_term, &
        ! matrix
        ptheta_2_local_stencil, &
        ! x
        x_data, map_any_space_2_x, ndf_any_space_2_x, undf_any_space_2_x &
    )
    call cpu_time(end)
    totaltime = end-start

    ! print out results
    write(*,*) "Reduction value:   ", sum(theta_adv_term_data)
    write(*,*) "Comutation time:   ", totaltime, " s"
    write(*,*) "Time/step:         ", totaltime/niters, " s"
    write(*,*) "Computation speed: ", (ncell_3d / 1000 * 6 * (2 * 8 - 1) * niters / totaltime) / 1000000 , " GFLOPs"
    write(*,*) "Mem bandwith usage:", memtotal * niters / 1000000 / totaltime , " GB/s"

end subroutine main

subroutine print_help_and_exit()
    write(*,*) "LFRic Matrix-Vector Multiplication Kernel Benchmark"
    write(*,*) "Usage: ./mv [-s NSIZE] [-l NLAYERS] [-n NITERS] [-t TRAVERSAL] [-h]"
    write(*,*) " -s/--size: The side size of an square horitzonal domain."
    write(*,*) " -l/--layers: The number of layers in the vertical domain."
    write(*,*) " -n/--niters: The number of iterations to execute."
    write(*,*) " -t/--traversal: The traversal strategy for the horitzonal domain."
    call exit(-1)
end subroutine print_help_and_exit

subroutine system_mem_usage(valueRSS)
    character(len=80) :: line
    integer, intent(out) :: valueRSS
    integer ::  ios, fu
    valueRSS = -1   

    open(newunit=fu, file='/proc/self/status', action='read')
    do
        read(fu, '(a)',iostat=ios ) line
        if(ios /=0) exit
        if(line(1:6) == 'VmRSS:') then
            read(line(7:), *) valueRSS
            exit
        endif
    enddo
    close(fu)
end subroutine system_mem_usage

 
subroutine populate_maps(nsize, nlayers, map_any_space_1_theta_adv_term, map_any_space_2_x)
  integer, intent(in) :: nsize, nlayers
  integer, intent(out), dimension(:,:) :: map_any_space_1_theta_adv_term
  integer, intent(out), dimension(:,:) :: map_any_space_2_x
  integer :: i2, i1, el

  ! Populate arrays
  do i2 = 0, nsize - 1
    do i1 = 0, nsize - 1
      ! Populate  map_any_space_1_theta_adv_term (8 vertices)
      do el = 1, 2
        map_any_space_1_theta_adv_term(1 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = (i2 * (nsize+1) + i1) * (nlayers + 1) + el
        map_any_space_1_theta_adv_term(2 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = (i2 * (nsize+1) + i1+1) * (nlayers + 1) + el
        map_any_space_1_theta_adv_term(3 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = ((i2+1) * (nsize+1) + i1+1) * (nlayers + 1) + el
        map_any_space_1_theta_adv_term(4 + 4*(el-1), (i2 * nsize) + i1 + 1 ) = ((i2+1) * (nsize+1) + i1) * (nlayers + 1) + el
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
end subroutine populate_maps

subroutine populate_strided_cmap(nsize, cmap)
    integer, intent(in) :: nsize
    integer, intent(out), dimension(:,:) :: cmap
    integer :: c1_index, c2_index, c3_index, c4_index, i1, i2

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

  ! if (traverse.eq."colouring2") then
  !   if (.not.(IAND (nsize, nsize-1) == 0)) then
  !     write(*,*) "Error colouring2 strategy just allowed for nsize is a power of 2"
  !     call exit(0)
  !   endif
  ! ! Number of cells per color = num times per row * num times per column
  !   ncp_colour(1) = ((nsize/2)+mod(nsize,2)) * ((nsize/2)+mod(nsize,2)) / 2
  !   ncp_colour(2) = ( nsize/2 )              * ((nsize/2)+mod(nsize,2)) / 2
  !   ncp_colour(3) = ((nsize/2)+mod(nsize,2)) * (nsize/2) / 2
  !   ncp_colour(4) = ( nsize/2 )              * (nsize/2) / 2
  !   write(*,*) ncp_colour

  !   allocate(cmap(ncolour,maxval(ncp_colour))) !color map (colour, cell number)

  !   c1_index = 1
  !   c2_index = 1
  !   c3_index = 1
  !   c4_index = 1
  !   do i2 = 0, nsize - 1
  !       do i1 = 0, nsize - 1
  !           if (mod(i1,4).eq.0) then
  !               if (mod(i2,2).eq.0) then
  !                   cmap(1,c1_index) = (i2 * nsize) + i1 + 1
  !                   c1_index = c1_index + 1
  !               else
  !                   cmap(2,c2_index) = (i2 * nsize) + i1 + 1
  !                   c2_index = c2_index + 1
  !               endif
  !           else if (mod(i1,2).eq.0) then
  !               if (mod(i2,2).eq.0) then
  !                   cmap(3,c3_index) = (i2 * nsize) + i1 + 1
  !                   c3_index = c3_index + 1
  !               else
  !                   cmap(4,c4_index) = (i2 * nsize) + i1 + 1
  !                   c4_index = c4_index + 1
  !               endif
  !           endif
  !       enddo
  !   enddo

  ! else if (traverse.eq."colouring-rows") then
  !     ! In this version ncp_colour contains the number of rows per color and
  !     ! cmap the first cell of each coloured row

  !   ncp_colour(1) = (nsize / 2) + mod(nsize,2)
  !   ncp_colour(2) = (nsize / 2)
  !   ncp_colour(3) = 0
  !   ncp_colour(4) = 0
  !   allocate(cmap(ncolour,maxval(ncp_colour))) !color map (colour, cell number)

  !   c1_index = 1
  !   c2_index = 1
  !       do i1 = 0, nsize - 1
  !           if (mod(i1,2).eq.0) then
  !                   cmap(1,c1_index) = (i1 * nsize) + 1
  !                   c1_index = c1_index + 1
  !           else
  !                   cmap(2,c2_index) = (i1 * nsize) + 1
  !                   c2_index = c2_index + 1
  !           endif
  !       enddo

  ! else
  !     write(*,*) "Error: No colouring strategy selected"
  ! endif
  
  ! ! Random data
end subroutine populate_strided_cmap

subroutine ascending_init(array)
    real(kind=r_def), dimension(:) :: array
    integer :: i

    do i = 1, size(array)
        array(i) = i
    enddo
end subroutine ascending_init
subroutine ascending_init_matrix(array)
    real(kind=r_def), dimension(:,:,:) :: array
    integer :: i, j, k

    do i = 1, size(array, 1)
        do j = 1, size(array, 2)
            do k = 1, size(array, 3)
                array(i, j, k) = i + j * size(array, 1) + k * size(array, 1) * size(array, 2)
            enddo
        enddo
    enddo
end subroutine ascending_init_matrix
end program kdriver
