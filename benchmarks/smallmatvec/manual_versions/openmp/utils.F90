module utils
    implicit none
    contains

    ! This function has only been tested on Linux OS with Intel Compiler.
    ! It reads the UNIX /proc files, so it is only portable to OSs supporting
    ! them. getpid() can be used on gfortran without importing any module.
    subroutine system_mem_usage(valueRSS)
        use ifport ! Import getpid() on Intel compiler
        implicit none
        integer, intent(out) :: valueRSS

        character(len=200):: filename=' '
        character(len=80) :: line
        character(len=8)  :: pid_char=' '
        integer :: pid
        logical :: ifxst

        valueRSS=-1    ! return negative number if not found

        !--- get process ID
        pid=getpid()
        write(pid_char,'(I8)') pid
        filename='/proc/'//trim(adjustl(pid_char))//'/status'

        !--- read system file
        inquire(file=filename, exist=ifxst)
        if (.not.ifxst) then
          write (*,*) 'system file does not exist'
          return
        endif

        open(unit=100, file=filename, action='read')
        do
          read (100,'(a)',end=120) line
          if (line(1:6).eq.'VmRSS:') then
             read (line(7:),*) valueRSS
             exit
          endif
        enddo
        120 continue
        close(100)

        return
    end subroutine system_mem_usage

end module utils
