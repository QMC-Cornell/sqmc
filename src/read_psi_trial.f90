module read_psi_trial

  use types, only: rk
  use common_ham, only: ndet
  use common_psi_t, only: ndet_psi_t, iwdet_psi_t, cdet_psi_t, psi_g
  implicit none

  contains

! ==============================================================================
  subroutine read_psi_t
! ------------------------------------------------------------------------------
! Description   : Read trial wavefunction
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------

  use common_run, only: ipr
  real(rk) norm_psi_t_sq
  integer istat,i

  allocate(psi_g(ndet),stat=istat)
  if(istat.ne.0) stop 'failed to allocate psi_g'

  read(5,*) ndet_psi_t
  write(6,'(/,''ndet_psi_t='',i5)') ndet_psi_t
  allocate(iwdet_psi_t(ndet_psi_t),stat=istat)
  if(istat.ne.0) stop 'failed to allocate iwdet_psi_t'
  allocate(cdet_psi_t(ndet_psi_t),stat=istat)
  if(istat.ne.0) stop 'failed to allocate cdet_psi_t'

  read(5,*) iwdet_psi_t
  write(6,'(''iwdet_psi_t='',1000i5)') iwdet_psi_t
  read(5,*) cdet_psi_t
  if(abs(minval(cdet_psi_t)).gt.maxval(cdet_psi_t,1)) cdet_psi_t=-cdet_psi_t ! flip sign if largest absolute value element is negative
  write(6,'(''cdet_psi_t='',1000d12.4)') cdet_psi_t

! Define psi_g
  norm_psi_t_sq=sum(cdet_psi_t(:)**2)
  write(6,'(/,''norm_psi_t_sq='',f8.5)') norm_psi_t_sq
  psi_g=sqrt((1-norm_psi_t_sq)/(ndet-ndet_psi_t))
  do i=1,ndet_psi_t
    psi_g(iwdet_psi_t(i))=cdet_psi_t(i)
  enddo
  write(6,'(''ipr,ndet'',9i5)') ipr,ndet
  if((ipr.ge.0.and.ndet.le.1000) .or. (ipr.ge.1.and.ndet.le.10000)) write(6,'(''psi_g='',1000d12.4)') psi_g

! Warning: tmp
! cdet_psi_t=sign(abs(cdet_psi_t)**.95,cdet_psi_t)

  cdet_psi_t=cdet_psi_t/sqrt(dot_product(cdet_psi_t,cdet_psi_t)) ! normalize psi_t
  write(6,'(/,''normalized cdet_psi_t='',1000d12.4)') cdet_psi_t

! Warning: tmp
! do i=1,ndet_psi_t
!   psi_g(iwdet_psi_t(i))=cdet_psi_t(i)
! enddo
! write(6,'(''psi_g='',1000d12.4)') psi_g

  end subroutine read_psi_t

end module read_psi_trial
