  subroutine my_second (n,title)
! Written by Cyrus Umrigar
! Prints out cpu and wall-clock time,
! both since beginning of run and since last call.

! I used to use etime for cpu time and time for wall-clock time.
! I now use cpu_time and system_clock because they are standard.
! I tested out time and cpu_time for cpu time and
! system_clock and secnds for wall-clock time (system_clock is standard)
! to see if they are more accurate.  All four gave variations of
! a few tenths of a second for a job that ran less than 10 secs.
! cpu_time and time were the same using g77 in Linux

!**** Edited by AR[7/23/13]: for now the job is done just by master_core, may be extended to act globally later 

  use types, only: rk, i8b
  
  use mpi_routines, only: master_core
  
  implicit none
  integer, intent(in) :: n
  character(len=*), intent(in) :: title

  real(rk) etim,etim1,etim2,etimtot,etimlast, avail,avail2,avail3
! integer itim,itim1, itim2, itimtot, itimlast, icall, icount_rate, icountmax
  integer(i8b) itim,itim1, itim2, itimtot, itimlast, icall, icount_rate, icountmax

  save icall,itim1,itim2,etim1,etim2
  data icall/0/

! etim=etime(i)        ! unix
! etim=etime(tarray)   ! linux
! etim=mclock()*1.d-2  ! aix
! etim=mclock()*1.d-6

  call cpu_time(etim)  ! standard but use pgf90, not pgf77

! itim=time()
! itim=time(0)       ! aix

  call system_clock(itim,icount_rate,icountmax) ! standard but but use pgf90, not pgf77

  if(icall.eq.0) then
    icall=1
    itim1=itim
    itim2=itim
    etim1=etim
    etim2=etim
  endif

  if(itim.ge.itim1) then
    itimtot=nint(real(itim-itim1,rk)/icount_rate)
  else
    itimtot=nint(real(itim-itim1+icountmax,rk)/icount_rate)
  endif
  if(itim.ge.itim2) then
    itimlast=nint(real(itim-itim2,rk)/icount_rate)
  else
    itimlast=nint(real(itim-itim2+icountmax,rk)/icount_rate)
  endif
  itim2=itim

  etimtot=etim-etim1
  etimlast=etim-etim2
  etim2=etim

  call mem_avail(avail,avail2,avail3)
! call mem_avail

  if(n.eq.1) write (6,'(/,''BEG OF '',a,'' CP, REAL TIME '', t66,2f11.2,2i9,'', Mem avail'',9f11.2)') title,etimtot,etimlast,itimtot,itimlast, avail,avail2,avail3
  if(n.eq.2) write (6,'(''END OF '',a,'' CP, REAL TIME '', t66,2f11.2,2i9,'', Mem avail'',9f11.2)') title,etimtot,etimlast,itimtot,itimlast, avail,avail2,avail3

#ifdef DEBUG
  call my_memory() ! Also print out memory usage when in debug mode.
#endif

  call flush(6)

  return
  end
!-----------------------------------------------------------------------------------------

subroutine mem_avail(avail,avail2,avail3)
!subroutine mem_avail
! Returns 3 values of the available memory in MB.
! avail  : does not take into account allocated memory, but just the assigned memory
! avail2 : takes into account allocated memory
! avail3 : gives the smallest memory that was available at any point in the run.
!
! Memory used:
! rss (resident set size)   : does not take into account allocated memory, but just the assigned memory
! vsz (virtual memory size) : takes into account allocated memory
! peak (peak virtual memory size) : peak memory use, takes into account allocated memory
! vsz, rss and peak seem to come out more nearly the same than avail, when run several times.
! On newer kernels, more /proc/meminfo gives something like:
! MemTotal:        6106348 kB
! MemFree:         1132248 kB
! MemAvailable:    5537604 kB
! ...
! MemAvailable is what we need in order to know what is the largest array we can allocate
! at any point in the program (it changes as the program runs).
! However, it is not available on older kernels.  So, this program calculates it.
! Note: although it is called mem_avail it actually gives MemAvailable
! Note that it uses kB to mean 2^10=1024 bytes rather than 1000 bytes, so we convert to decimal.
! Written by Bastian Schaefer, Junhao Li and Cyrus Umrigar, Jan 2017

  implicit none
  real(8), intent(out) :: avail, avail2,avail3
  integer ierr
  character(len=80) line
  character(len=64) keyword
  real(8) val, memfree, activefile, inactivefile, SReclaimable, low
  logical file_open
  
  character(len=128) filename
  character(len=8) pid_str
  integer pid
  real(8) vsz, rss, peak
  logical file_exist
! save file_open
#ifdef MAC
!MJO MacOS does not have the same files as linux, so we cannot do mem_avail in this way.
  return 
#endif
  memfree=-1.d0 ; activefile=-1.d0 ; inactivefile=-1.d0 ; SReclaimable=-1.d0
  
  inquire(unit=33, opened=file_open)
  if(file_open) then
    rewind(33)
  else
   !write(6,'(''opening file 33'')')
    open(33,file='/proc/zoneinfo')
  endif
  low=0.d0
  do
    read(33,'(a)',iostat=ierr)line
    if(ierr/=0) exit
    read(line,*)keyword
    if(keyword=='low')then
      read(line,*)keyword,val
      low=low+val
    endif
  enddo
! close(33)

  inquire(unit=34, opened=file_open)
  if(file_open) then
    rewind(34)
  else
   !write(6,'(''opening file 34'')')
    open(34,file='/proc/meminfo')
  endif

  do
    read(34,'(a)',iostat=ierr)line
    if(ierr/=0) exit
    read(line,*)keyword,val
    select case (keyword)
       case('MemFree:')
         memfree=val
       case('Active(file):')
         activefile=val
       case('Inactive(file):')
         inactivefile=val
       case('SReclaimable:')
         SReclaimable=val
    end select
  enddo
! close(34)
  
  if(memfree < 0d0 .or. activefile < 0.d0 .or. &
     inactivefile < 0.d0 .or. SReclaimable < 0.d0) stop 'value in meminfo missing'
  
  avail = memfree + activefile + inactivefile + SReclaimable - 12.d0*low
  ! Divide by 8 to get double precision words
  !write(6,'(''memfree, activefile, inactivefile, SReclaimable, low, avail'',9es12.4)') &
  ! memfree*1024/8, activefile*1024/8, inactivefile*1024/8, SReclaimable*1024/8, low*1024/8, avail*1024/8
  avail = 1024e-6*avail !convert to MB
! write(6,'(''Memory avail (MB)='',f9.2)') avail

! avail does not take into account memory that is allocated, if it has not been assigned.
! So, calculate avail2, which does take this into account, and avail3 which shows the least memory avail at any point till now.

! If we keep the file open then it does not update, even though this is not a problem for /proc/meminfo and /proc/zoneinfo
! inquire(unit=35, opened=file_open)
! if(file_open) then
!   rewind(35)
! else
    pid = getpid()
    write(pid_str, '(i8)') pid
    filename = '/proc/' // trim(adjustl(pid_str)) // '/status'
!   write(6,*) filename
    
    inquire(file=filename, exist=file_exist)
    if (.not. file_exist) then
      write (6,'(''Error: Memory file for this pid does not exist.'')')
      return
    endif
    open(35,file=filename) ! do not use action='read', because then it will not update
! endif

  vsz = -1 ; rss = -1 ; peak = -1

  do
    read(35,'(a)') line
    if(line(1:6) .eq. 'VmRSS:') then
      read(line(7:),*) rss
      exit
    elseif(line(1:7) .eq. 'VmSize:') then
      read(line(8:),*) vsz
    elseif(line(1:7) .eq. 'VmPeak:') then
      read(line(8:),*) peak
    endif
  enddo
  close(35)

  vsz = 1024e-6*vsz ; rss = 1024e-6*rss ; peak = 1024e-6*peak !convert to MB

  avail2=avail+rss-vsz
  avail3=avail+rss-peak
! write(6,'(''Mem used (vsz, rss, peak), Mem avail, Total (vsz, rss) (MB)'',9f12.2)') vsz, rss, peak, avail, avail+vsz, avail+rss
! write(6,'(''Memory avail (MB)='',9f9.2)') avail, avail2, avail3
! call flush(6)

end subroutine mem_avail
!--------------------------------------------------------------------
