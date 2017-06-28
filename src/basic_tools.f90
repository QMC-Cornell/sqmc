module basic_tools

  contains

! ================================================================================
  subroutine get_date (date_nice)
! --------------------------------------------------------------------------------
! Description : return system current date in nice format
!
! Created     : J. Toulouse, 10 Jul 2007
! --------------------------------------------------------------------------------
  implicit none

! output
  character (len=*), intent(out) :: date_nice

! local:
  character (len=32) :: date, time, zone

! begin
  call date_and_time (date, time, zone)
  write (date_nice,'(a4,a,a2,a,a2,x,a2,a,a2,a,a2,a,a5,a)') date(1:4),'-',date(5:6),'-',date(7:8),time(1:2),':',time(3:4),':',time(4:6),' (',zone(1:5),')'

  end subroutine get_date

end module basic_tools
