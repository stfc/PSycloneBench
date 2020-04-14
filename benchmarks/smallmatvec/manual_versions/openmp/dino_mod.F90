!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.MetOffcie which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief IO library for the PSKE to inject into PSy layer for kernel extraction,
!! DO NOT use in general.
!>
!> @details Creates a Dynamo input-output object (dino). Can read-write scalars 
!! and arrays. Procedures are overloaded so code generation in PSy layer is 
!! simpler.
!! THIS IS SERIAL ONLY - Multiple MPI ranks will write to the same file.
!> @note dino_type uses a fixed file handle for simplicity. This could cause 
!! problems if more than one dino_type is instanciated at once, hence the 
!! close subroutine in  addition to the destructor.
!> 

module dino_mod
  use constants_mod, only : r_def, i_def, str_max_filename
  implicit none
  private
  type, public :: dino_type
     private
     integer :: file_handle
     logical :: file_is_open = .false.
   contains
     generic   :: output_scalar => output_integer, output_real
     generic   :: input_scalar => input_integer, input_real
     generic   :: output_array => output_2d_integer_array, output_1d_real_array, &
                                  output_3d_real_array, output_1d_integer_array
     generic   :: input_array => input_2d_integer_array, input_1d_real_array, &
                                 input_3d_real_array, input_1d_integer_array
     procedure :: output_integer
     procedure :: output_real
     procedure :: input_integer
     procedure :: input_real
     procedure :: io_close
     procedure :: output_2d_integer_array
     procedure :: input_2d_integer_array
     procedure :: output_1d_integer_array
     procedure :: output_1d_real_array
     procedure :: input_1d_real_array
     procedure :: output_3d_real_array
     procedure :: input_3d_real_array
     procedure :: input_1d_integer_array

  end type dino_type

  interface dino_type
     module procedure output_constructor
  end interface
contains

  type(dino_type)  function output_constructor() result(self)
    implicit none
    ! open a file
    integer(kind=i_def) :: ierr
    character(str_max_filename) :: fname
    character(str_max_filename) :: iomsg
    self%file_handle = 789
    
    write(fname,'(A)') "dinodump.dat"
    open(unit=self%file_handle,file=fname, status="unknown",iostat=ierr,iomsg=iomsg)
    if(ierr/=0) then
       write(*,'(A," ",A)') "output_constructor:cannot open file ",trim(fname), &
            iomsg
       stop
    end if
    self%file_is_open = .true.
  end function output_constructor

  subroutine output_3d_real_array(self, array, dim1, dim2, dim3)
    implicit none
    class(dino_type),                            intent(in) :: self
    integer(kind=i_def),                         intent(in) :: dim1
    integer(kind=i_def),                         intent(in) :: dim2
    integer(kind=i_def),                         intent(in) :: dim3
    real(kind=r_def), dimension(dim1,dim2,dim3), intent(in) :: array
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    write(self%file_handle,*) array

  end subroutine output_3d_real_array

  subroutine input_3d_real_array(self, array, dim1, dim2, dim3)
    implicit none
    class(dino_type),                            intent(in)  :: self
    integer(kind=i_def),                         intent(in)  :: dim1
    integer(kind=i_def),                         intent(in)  :: dim2
    integer(kind=i_def),                         intent(in)  :: dim3
    real(kind=r_def), dimension(dim1,dim2,dim3), intent(out) :: array
    integer :: i,j,k
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    read(self%file_handle,*) array

  end subroutine input_3d_real_array

  subroutine output_1d_integer_array(self, array, dim)
    implicit none
    class(dino_type),                 intent(in) :: self
    integer(kind=i_def),              intent(in) :: dim
    integer(kind=i_def), dimension(dim), intent(in) :: array
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    write(self%file_handle,*) array
    
  end subroutine output_1d_integer_array

  subroutine input_1d_integer_array(self, array, dim)
    implicit none
    class(dino_type),                    intent(in)  :: self
    integer(kind=i_def),                 intent(in)  :: dim
    integer(kind=i_def), dimension(dim), intent(out) :: array
    if(.not.self%file_is_open) then
       write(*,'(A)') "input_integer:Shriek, file not open"
       stop
    end if
    read(self%file_handle,*) array
    
  end subroutine input_1d_integer_array

  subroutine output_1d_real_array(self, array, dim)
    implicit none
    class(dino_type),                 intent(in) :: self
    integer(kind=i_def),              intent(in) :: dim
    real(kind=r_def), dimension(dim), intent(in) :: array
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    write(self%file_handle,*) array

  end subroutine output_1d_real_array

  subroutine input_1d_real_array(self, array, dim)
    implicit none
    class(dino_type),                  intent(in)  :: self
    integer(kind=i_def),              intent(in)  :: dim
    real(kind=r_def), dimension(dim), intent(out) :: array
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    read(self%file_handle,*) array

  end subroutine input_1d_real_array

  subroutine output_2d_integer_array(self, array, dim1, dim2)
    implicit none
    class(dino_type),                           intent(in) :: self
    integer(kind=i_def),                       intent(in) :: dim1
    integer(kind=i_def),                       intent(in) :: dim2
    integer(kind=i_def), dimension(dim1,dim2), intent(in) :: array
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    write(self%file_handle,*) array

  end subroutine output_2d_integer_array

  subroutine input_2d_integer_array(self, array, dim1, dim2)
    implicit none
    class(dino_type),                           intent(in)  :: self
    integer(kind=i_def),                       intent(in)  :: dim1
    integer(kind=i_def),                       intent(in)  :: dim2
    integer(kind=i_def), dimension(dim1,dim2), intent(out) :: array
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    read(self%file_handle,*) array

  end subroutine input_2d_integer_array
  
  subroutine output_integer(self,scalar)
    implicit none
    class(dino_type),  intent(inout) :: self
    integer(kind=i_def),            intent(in)    :: scalar
!    character(str_max_filename), intent(in)    :: label

    ! check the file is open
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    write(self%file_handle,*) scalar
    
  end subroutine output_integer

  subroutine output_real(self,scalar)
    implicit none
    class(dino_type), intent(inout) :: self
    real(kind=r_def),               intent(in)    :: scalar
!    character(str_max_filename), intent(in)    :: label

    ! check the file is open
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    write(self%file_handle,*) scalar
    
  end subroutine output_real

  subroutine input_integer(self,scalar)
    implicit none
    class(dino_type),  intent(inout) :: self
    integer(kind=i_def),            intent(out)    :: scalar

    ! check the file is open
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    read(self%file_handle,*) scalar
  end subroutine input_integer

  subroutine input_real(self,scalar)
    implicit none
    class(dino_type),  intent(inout) :: self
    real(kind=r_def),             intent(out)   :: scalar

    ! check the file is open
    if(.not.self%file_is_open) then
       write(*,'(A)') "output_integer:Shriek, file not open"
       stop
    end if
    read(self%file_handle,*) scalar
  end subroutine input_real

  subroutine io_close(self)
    implicit none
    class(dino_type), intent(inout) :: self

    if(self%file_is_open) then
       close(self%file_handle)
       self%file_is_open=.false.
    end if
  end subroutine io_close
  
end module dino_mod

