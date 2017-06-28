subroutine my_memory()
    ! Print out VSS, RSS by the current program and free memory.
    use tools, only : get_free_memory
    implicit none
    
    integer :: pid
    real :: vss, rss, peak, free
    character(len = 200) :: filename
    character(len = 80) :: line
    character(len = 8) :: pid_str
    logical :: file_exist

    pid = getpid()
    write(pid_str, '(I8)') pid
    filename = '/proc/' // trim(adjustl(pid_str)) // '/status'

    vss = -1
    rss = -1
    peak = -1
    
    inquire(file = filename, exist = file_exist)
    if (.not. file_exist) then
        write (*, *) 'Error: Memory file does not exist.'
        return
    endif

    open(unit = 100, file = filename, action = 'read')
    do
        read(100, '(A)') line
        if (line(1:6) .eq. 'VmRSS:') then
            read(line(7:), *) rss
            exit
        elseif (line(1:7) .eq. 'VmSize:') then
            read(line(8:), *) vss
        elseif (line(1:7) .eq. 'VmPeak:') then
            read(line(8:), *) peak
        endif
    enddo
    close(100)

    vss = vss / 1e3
    rss = rss / 1e3
    peak = peak / 1e3
    ! free = get_free_memory() / 1e6
    free = -1

    write (*, '(A, 4F12.2)') 'MEMORY USAGE VSS / RSS / PEAK / FREE (MB): ', vss, rss, peak, free
    call flush(6)

end subroutine my_memory
