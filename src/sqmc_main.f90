program sqmc_main

! SQMC (Semistochastic Quantum Monte Carlo) is a program to do semistochastic projection, to find
! an extremal eigenvalue and expectation values for the corresponding eigenfunction of large matrices.
! The main authors of SQMC are:
! Cyrus J, Umrigar, Frank R. Petruzielo, Hitesh J. Changlani, Adam A. Holmes and Alessandro Roggero.
! The idea of semistochastic projection was conceived in a conversation between Peter Nightingale and Cyrus Umrigar.
! The method is described in:
! "Semistochastic Projector Monte Carlo Method",
! F. R. Petruzielo, A. A. Holmes, Hitesh J. Changlani, M. P. Nightingale, and C. J. Umrigar, PRL 109, 230201 (2012).
! This work was motivated by the i-FCIQMC method of Ali Alavi's group at Cambridge.
! Their method is described in:
! G. H. Booth, A. J. W. Thom, and A. Alavi, J. Chem.  Phys. 131, 054106 (2009).
! D. Cleland, G. H. Booth, and A. Alavi, J. Chem. Phys.  132, 041103 (2010).
! and subsequent papers.

!**** Edited by AR[7/23/13]: added MPI support
  use mpi_routines, only : mpi_barr, cluster_init, cluster_finalize, master_core
  use do_walk
  use basic_tools, only : get_date

  character (32) :: date = ''
  character (16) :: hostname = ''
  character (16) :: user = ''

  call cluster_init()

  write(6,*)
  write(6,'(8x,''**************************************************************'')')
  write(6,'(8x,''** SQMC (Semistochastic Quantum Monte Carlo)                **'')')
  write(6,'(8x,''** SHCI (Semistochastic Heat-bath Configuration Interation) **'')')
  write(6,'(8x,''**************************************************************'',/)')

  call get_date (date)
  call get_environment_variable ("HOSTNAME", hostname)
  call get_environment_variable ("USER", user)
  write(6,'(6a)') 'Executed by ',trim(user),' on ',trim(date), ' on master host ',trim(hostname)

  call my_second(1,'sqmc')
 
  call mpi_barr()
  flush(6) 
  call read_input

  call my_second(2,'sqmc')

  call cluster_finalize()
  
end program sqmc_main
