module hamiltonian_mod

  use types, only: rk, i4b, i8b, i16b, ik, ik_vec
#ifdef NUM_ORBITALS_GT_127
  use overload
#endif
  use read_psi_trial, only: read_psi_t
  use heg, only: read_heg, energy_pieces_heg, off_diagonal_move_heg,system_setup_heg
  use chemistry, only: read_chem, energy_pieces_chem, off_diagonal_move_chem_efficient_heatbath, off_diagonal_move_chem,system_setup_chem,off_diagonal_move_chem_cauchySchwarz
  use hubbard, only : nmultidet,coeffs,dets_up_multi,dets_dn_multi,wf_type, nsites,read_hubbard, energy_pieces_hubbard,energy_pieces_hubbard_k,off_diagonal_move_hubbard,off_diagonal_move_hubbard_k,off_diagonal_move_hubbard_dm,energy_exact_hubbard,ndet_hubbard,system_setup_hubbard
  use semistoch, only : generate_space_iterate,perform_selected_ci,perform_truncated_lanczos
  use common_run, only: tau, tau_multiplier, run_type
! use commons/common_ham, only: nelec,norb
  use common_ham, only: hamiltonian_type, nelec, nup, ndn, norb, n_core_orb, ndet, diagonalize_ham
  use common_psi_t, only : trial_wf_iters, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf, ndet_psi_t, iwdet_psi_t, dets_up_psi_t, dets_dn_psi_t, cdet_psi_t, ndet_psi_t_in, dets_up_psi_t_in, dets_dn_psi_t_in, psi_g
  use common_selected_ci, only : norb_det_sel,n_sym_uniq_det_det_sel,ndet_det_sel,dets_up_det_sel,dets_dn_det_sel,cdet_det_sel,det_sel_iters,lanczos_iters,lanczos_initiators,lanczos_truncate

  implicit none

! integer ndet
! real(rk) :: tau_multiplier, tau, energy_exact
  real(rk) :: energy_exact
  real(rk), allocatable :: ham(:,:),wavefn_exact(:)
! integer, private :: nelec,norb,bosonic,ipr
  integer, private :: bosonic,ipr
  real(rk),private :: spectrum_coef, spectrum_power, ham_diag_fluc, ham_offdiag_fluc

  contains

! ==============================================================================
  subroutine hamiltonian
! ------------------------------------------------------------------------------
! Description   : Create Hamiltonian matrix at beginning of run
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
!*** Edited by AR [7/24/13]: changed I/O to handle parallelism [just for chemistry right now]
    use mpi_routines, only:master_core,mpi_bsend
!***

  use tools, only: n_choose_k
  use common_ham, only: diagonal_ham_lowest, diagonal_ham_highest
  use common_walk, only: n_connected_dets, m_connected_dets
  use common_imp, only: semistochastic
  use common_psi_t, only : ndet_psi_t_in,dets_up_psi_t_in,dets_dn_psi_t_in,cdet_psi_t_in,use_psit_con_in,dtm_energy
  use chemistry, only : diagonalize_sparse_hamiltonian_chem,setup_efficient_heatbath
  use hubbard, only: diagonalize_sparse_hamiltonian_hubbard
  use heg, only : setup_efficient_heatbath_heg
  use mpi_routines, only : whoami,mpi_barr
  use common_selected_ci, only : cdets, tdets, eps_var, eps_pt, eps_pt_big, target_error, n_states
  use hci, only :  perform_hci
 !use hci, only : perform_cisdtq, perform_cisd_pt, perform_hci
! local
  integer istat
  integer  l_x, l_y, i, j, b, n_alpha, n_beta, num_alpha_configs, num_beta_configs
  integer  nelec_val, nup_val, ndn_val, norb_val, m_connected_dets_1excit, m_connected_dets_2excit
  real(rk) ndet_estim
 !integer(i16b) ndet_estim
  logical  pbc,neel_up_only
! real(rk) t, U, diagonal_ham_highest, diagonal_ham_lowest
  real(rk) t, U
  integer, allocatable :: bra_config(:),ket_config(:)
! integer signmod
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec) :: hf_up,hf_dn
#else
  integer(ik) :: hf_up,hf_dn
#endif

! character*16 hamiltonian_type

!***b
  if(master_core) then
    read(5,*) hamiltonian_type,ipr
    write(6,'(''hamiltonian_type= '',a)') trim(hamiltonian_type)
  endif
  call mpi_bsend(hamiltonian_type)
  call mpi_bsend(ipr)
!  read(5,*) hamiltonian_type,ipr
!  write(6,'(''hamiltonian_type= '',a)') trim(hamiltonian_type)
!***e
  if(hamiltonian_type.eq.'fictitious') then
    read(5,*) bosonic, spectrum_coef, spectrum_power, ham_diag_fluc, ham_offdiag_fluc ! bosonic=0 => sign problem, bosonic=1 => no sign problem
    write(6,'(''bosonic, spectrum_coef, spectrum_power, ham_diag_fluc, ham_offdiag_fluc ='',i2,9f8.3)') &
 &  bosonic,spectrum_coef, spectrum_power, ham_diag_fluc, ham_offdiag_fluc

!   read(5,*) nelec,norb ! At present this is just used to set ndet (number of determinants (basis functions))
!   ndet=fact(norb)/(fact(nelec)*fact(norb-nelec))
!   write(6,'(''nelec,norb,ndet='',9i5)') nelec,norb,ndet
    read(5,*) ndet
    write(6,'(''ndet='',9i5)') ndet
    call flush(6)

    allocate(ham(ndet,ndet),stat=istat)
    if(istat.ne.0) stop 'failed to allocate ham'

    call hamiltonian_fictitious
  elseif(hamiltonian_type.eq.'heg') then
    call read_heg
    write (6,'(''trial_wf_iters, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf:'',2i4,2i9)') &
   &             trial_wf_iters, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf
    call flush(6)

    if (run_type .eq. 'hci') then
      call setup_efficient_heatbath_heg()
      call my_second(1, 'HCI_HEG')
      call perform_hci(eps_pt, eps_pt_big, target_error, n_states=n_states)
      call my_second(2, 'HCI_HEG')
      call mpi_barr()
      stop
    endif

    if (ndet_psi_t_in>0) then
      call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=ndet_psi_t,dets_up=dets_up_psi_t,dets_dn=dets_dn_psi_t,values=cdet_psi_t,init_up=dets_up_psi_t_in,init_dn=dets_dn_psi_t_in,init_wts=cdet_psi_t_in,dtm_energy=dtm_energy) ! Generate psi trial here now instead of in hamiltonian
    elseif (.not.use_psit_con_in) then
      call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=ndet_psi_t,dets_up=dets_up_psi_t,dets_dn=dets_dn_psi_t,values=cdet_psi_t,dtm_energy=dtm_energy) ! Generate psi trial here now instead of in hamiltonian
    endif
  elseif(hamiltonian_type.eq.'chem') then

    call read_chem
    call system_setup_chem(hf_up,hf_dn)
    if(run_type.ne.'hci') call read_psit

    ! Set up efficient heatbath sampling probabilities (does nothing if proposal_method is not fast_heatbath)
    call setup_efficient_heatbath()

    ! Perform CISDTQ, then exit
    if (run_type.eq.'cisdtq') then
      call my_second(1,'CISDTQ')
      stop "CISDTQ no longer working"
     !call perform_cisdtq(eps_var)
      call my_second(2,'CISDTQ')
      stop
    endif

    ! Perform CISD+2PT, then exit
    if (run_type.eq.'cisd_pt') then
      call my_second(1,'CISD+2PT')
      stop "CISD+2PT no longer working"
     !call perform_cisd_pt(eps_pt)
      call my_second(2,'CISD+2PT')
      stop
    endif

    ! Perform Heat-bath CI, then exit
    if (run_type.eq.'hci') then
      call my_second(1,'HCI')
      if (hf_up==0_ik) then
        call perform_hci(eps_pt, eps_pt_big, target_error, n_states=n_states)
      else
        call perform_hci(eps_pt, eps_pt_big, target_error, n_states, hf_up, hf_dn)
      endif
      call my_second(2,'HCI')
      call mpi_barr()
      stop
    endif

    ! Edited by AAH on 1 May 2012
    if (run_type.eq.'selected_ci') then
      ! Deterministically selected subspace
!***b
      if(master_core) then
        call my_second(1,'generate deterministically selected subspace')
        call flush(6)
        write (6,*) "det_sel_iters=",det_sel_iters
        write (6,*) "norb_det_sel=",norb_det_sel
        write (6,*) "n_sym_uniq_det_det_sel=",n_sym_uniq_det_det_sel
        call flush(6)
      endif
!      call my_second(1,'generate deterministically selected subspace')
!      call flush(6)
!      write (6,*) "det_sel_iters=",det_sel_iters
!      write (6,*) "norb_det_sel=",norb_det_sel
!      write (6,*) "n_sym_uniq_det_det_sel=",n_sym_uniq_det_det_sel
!      call flush(6)
!***e
      call perform_selected_ci(det_sel_iters,norb_det_sel,n_sym_uniq_det_det_sel)
      stop
    endif
    if (run_type.eq.'trunc_lanc') then
      call my_second(1,'truncated lanczos')
      call perform_truncated_lanczos(lanczos_iters,lanczos_initiators,lanczos_truncate)
      call my_second(2,'truncated lanczos')
      stop
    endif
    if (.not.use_psit_con_in) then
    call my_second(1,'generate psi trial')
!***b
      if(master_core) then
          write (6,*) "trial_wf_iters:, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf:"
          write (6,*) trial_wf_iters, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf
          call flush(6)
      endif
    endif
!    write (6,*) "trial_wf_iters:, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf:"
!    write (6,*) trial_wf_iters, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf
!    call flush(6)
!***e
    if (ndet_psi_t_in>0) then
      call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=ndet_psi_t,dets_up=dets_up_psi_t,dets_dn=dets_dn_psi_t,values=cdet_psi_t,init_up=dets_up_psi_t_in,init_dn=dets_dn_psi_t_in,init_wts=cdet_psi_t_in,dtm_energy=dtm_energy) ! Generate psi trial here now instead of in hamiltonian
    elseif (.not.use_psit_con_in) then
      call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=ndet_psi_t,dets_up=dets_up_psi_t,dets_dn=dets_dn_psi_t,values=cdet_psi_t,dtm_energy=dtm_energy) ! Generate psi trial here now instead of in hamiltonian
    endif
   !call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=ndet_psi_t,dets_up=dets_up_psi_t,dets_dn=dets_dn_psi_t,values=cdet_psi_t) ! Generate psi trial here now instead of in hamiltonian
    call my_second(2,'generate psi trial')
!***b
    if(master_core) call flush(6)
!    call flush(6)
!***E
    ! End of edit


  elseif(hamiltonian_type.eq.'hubbard2') then
    call read_hubbard
    call read_psit
    if (wf_type .eq. 'gutz_multi' .or. wf_type .eq. 'cgutz_multi') then
        write (6,*) 'wf type is gutz_multi or cgutz_multi'
        call flush(6)
        hamiltonian_type='hubbardk'
        write (6,*) 'Calling generate_space_iterate - hubbard2/k'
        call flush(6)
        if (ndet_psi_t_in>0) then
          call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=nmultidet,dets_up=dets_up_multi,dets_dn=dets_dn_multi,values=coeffs,init_up=dets_up_psi_t_in,init_dn=dets_dn_psi_t_in,init_wts=cdet_psi_t_in)
        elseif (.not.use_psit_con_in) then
          call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=nmultidet,dets_up=dets_up_multi,dets_dn=dets_dn_multi,values=coeffs)
        endif
        write (6,*) 'Out of Calling generate_space_iterate - hubbard2/k'
        write (6,*) "Coefficients for multi determinantal wavefunction in Hubbard 2"
        call flush(6)
        write (6,*) coeffs
        call flush(6)
        hamiltonian_type='hubbard2'
    endif
    call system_setup_hubbard
    energy_exact=energy_exact_hubbard
    ndet=int(ndet_hubbard,i4b)
  elseif(hamiltonian_type.eq.'hubbardk') then
    call read_hubbard
    call system_setup_hubbard
    energy_exact=energy_exact_hubbard
    ndet=int(ndet_hubbard,i4b)
    call read_psit
    ! Edited by AAH on 1 May 2012
    if (run_type.eq.'selected_ci') then
      ! Deterministically selected subspace
      call my_second(1,'generate deterministically selected subspace')
      call flush(6)
      call perform_selected_ci(det_sel_iters,norb_det_sel,n_sym_uniq_det_det_sel)
      stop
    endif
    if (run_type.eq.'trunc_lanc') then
      call my_second(1,'truncated lanczos')
      call perform_truncated_lanczos(lanczos_iters,lanczos_initiators,lanczos_truncate)
      call my_second(2,'truncated lanczos')
      stop
    endif
    call my_second(1,'generate psi trial')
    write (6,'(''trial_wf_iters, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf:'',2i4,2i9)') &
   &             trial_wf_iters, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf
    call flush(6)
    if (ndet_psi_t_in>0) then
      call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=ndet_psi_t,dets_up=dets_up_psi_t,dets_dn=dets_dn_psi_t,values=cdet_psi_t,init_up=dets_up_psi_t_in,init_dn=dets_dn_psi_t_in,init_wts=cdet_psi_t_in,dtm_energy=dtm_energy) ! Generate psi trial here now instead of in hamiltonian
    elseif (.not.use_psit_con_in) then
      call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=ndet_psi_t,dets_up=dets_up_psi_t,dets_dn=dets_dn_psi_t,values=cdet_psi_t,dtm_energy=dtm_energy) ! Generate psi trial here now instead of in hamiltonian
    endif
   !call generate_space_iterate(trial_wf_iters,norb_trial_wf,n_initiators_trial_wf,n_truncate_trial_wf,n_det=ndet_psi_t,dets_up=dets_up_psi_t,dets_dn=dets_dn_psi_t,values=cdet_psi_t) ! Generate psi trial here now instead of in hamiltonian
    call my_second(2,'generate psi trial')
    call flush(6)
    ! End of edit
  elseif(hamiltonian_type.eq.'hubbarddm') then
    call read_hubbard
    ndet=int(ndet_hubbard,i4b)
  elseif(hamiltonian_type.eq.'hubbard') then
    read(5,*) l_x,l_y
    write(6,'(''l_x, l_y='', 9i4)') l_x,l_y
    read(5,*) pbc, neel_up_only
    write(6,'(''pbc, neel_up_only='',2l3)') pbc, neel_up_only
    read(5,*) t,U
    write(6,'(''t, U='', 9f8.2)') t,U
    read(5,*) n_alpha,n_beta
    write(6,'(''n_alpha, n_beta, filling factors='', 2i4, 9f8.2)') n_alpha,n_beta,real(n_alpha,rk)/(l_x*l_y),real(n_beta,rk)/(l_x*l_y)

    allocate(bra_config(2*l_x*l_y),stat=istat)
    allocate(ket_config(2*l_x*l_y),stat=istat)

    t=-t ! To make notation of t consistent with hubbard2 (HJC added this on April 7 2011)
    call binomial(10,4,b)
    if(ipr.ge.5) print *,"Binomial(10,4) = ",b

    call binomial(l_x*l_y,n_alpha,num_alpha_configs)
    call binomial(l_x*l_y,n_beta,num_beta_configs)

    write(6,'(''num_alpha_configs, num_beta_configs, total_configs='',2i8,i16)') &
 &   num_alpha_configs, num_beta_configs, num_alpha_configs*num_beta_configs
    call flush(6)

    ndet=num_alpha_configs*num_beta_configs
    allocate(ham(ndet,ndet),stat=istat)
    if(istat.ne.0) write(6,'(''Failed to allocate ham'')')
    call flush(6)

    call make_hubbard_matrix_2d(l_x,l_y, pbc, n_alpha, n_beta, t, U, num_alpha_configs*num_beta_configs,ham,ipr)

    if((ipr.ge.0.and.size(ham,1).le.100) .or. (ipr.ge.2.and.size(ham,1).le.10000)) then
      print *,"Hamiltonian"
      do i=1,size(ham,1)
        write (6,'(i5,8x,10000f10.5)') i,ham(i,:)
      enddo
    endif

    deallocate(bra_config)
    deallocate(ket_config)

!   print *,"S(6,8)=",signmod(6,8)
!   print *,"S(-6,8)=",signmod(-6,8)
!   print *,"S(-2,8)=",signmod(-2,8)
!   print *,"S(-1,8)=",signmod(-1,8)
!   print *,"S(2,3)=",signmod(2,3)
!   print *,"S(-2,3)=",signmod(-2,3)

   elseif(hamiltonian_type.eq.'read') then
    open(2,file='hamiltonian',form='formatted')
    read(2,*) ndet
    write(6,'(''Reading Hamiltonian of dimension'', i6)') ndet
    allocate(ham(ndet,ndet),stat=istat)
    if(istat.ne.0) write(6,'(''Failed to allocate ham'')')
    do i=1,ndet
      read(2,*) ham(i,:)
      if(ipr.ge.2) write(6,'(10000f9.5)') ham(i,:)
    enddo
    close(2)
   else
    stop 'hamiltonian_type must be fictitious or hubbard or read'

  endif

  ! Diagonalize ham for Hubbard 2 read while performing read_hubbard
  if ((hamiltonian_type .eq. 'read') .or. (hamiltonian_type .eq. 'fictitious')) then
      read(5,*) diagonalize_ham
      write(6,'(/,''diagonalize_ham='',i3)') diagonalize_ham ; call flush(6)
  endif

! Calculate the maximum number of connected determinants (including the one to itself)
  if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read') then
    m_connected_dets=0
    do j=1,ndet
      n_connected_dets=0
      do i=1,ndet
        if(ham(i,j).ne.0._rk) n_connected_dets=n_connected_dets+1
      enddo
      m_connected_dets=max(m_connected_dets,n_connected_dets)
    enddo
   elseif(hamiltonian_type.eq.'hubbard') then
    m_connected_dets=l_x*l_y*6 ! 6 for 3-dim
   elseif(hamiltonian_type.eq.'heg' .or. hamiltonian_type.eq.'chem') then
!   m_connected_dets=1+(nelec*(nelec-1)/2) * ((2*norb-nelec)/2) ! for heg?
! The calculation below is for chem.  It is an overestimate if there are orbitals of different symmetry.
    nelec_val=nelec-2*n_core_orb
    write (6,*) "nelec_val=",nelec_val
    nup_val=nup-n_core_orb
    write (6,*) "nup_val=",nup_val
    ndn_val=ndn-n_core_orb
    write (6,*) "ndn_val=",ndn_val
    norb_val=norb-n_core_orb
    write (6,*) "norb_val=",norb_val
    m_connected_dets_1excit=nelec_val*max(norb_val-nup_val,norb_val-ndn_val) ! single excit
    write (6,*) "m_connected_dets_1excit=",m_connected_dets_1excit
    m_connected_dets_2excit=(nup_val*(nup_val-1)*(norb_val-nup_val)*(norb_val-nup_val-1) + ndn_val*(ndn_val-1)*(norb_val-ndn_val)*(norb_val-ndn_val-1))/4 & ! par-spin double excit
 &  + nup_val*ndn_val*(2*norb_val-nelec_val)*max(norb_val-nup_val,norb_val-ndn_val)/2 ! antipar-spin double excit
    write (6,*) "m_connected_dets_2excit=",m_connected_dets_2excit
    m_connected_dets=1+m_connected_dets_1excit+m_connected_dets_2excit
    write (6,*) "m_connected_dets=",m_connected_dets
    ndet_estim=exp(log(real(n_choose_k(norb,nup)))+log(real(n_choose_k(norb,ndn))))
    write(6,'(/,''nelec, norb, m_connected_dets_1excit, m_connected_dets_2excit, m_connected_dets, ndet_estim='',2i4,i5,2i7,i15,'' ='',es11.4)') nelec, norb, m_connected_dets_1excit, m_connected_dets_2excit, m_connected_dets, int(ndet_estim,i8b), ndet_estim
   elseif(hamiltonian_type.eq.'hubbard2') then
    m_connected_dets=(nelec*4)+1
   elseif(hamiltonian_type.eq.'hubbardk') then
    m_connected_dets=(nelec*nelec*nsites)
    print *,"nsites=",nsites
    print *,"nelec=",nelec
   elseif(hamiltonian_type.eq.'hubbarddm') then
    m_connected_dets=1000
   endif

!***b
  if(master_core) then
    write(6,'(''m_connected_dets='',i6)') m_connected_dets
    call flush(6)
  endif
!  write(6,'(''m_connected_dets='',i6)') m_connected_dets
!  call flush(6)
!***e

  if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'hubbard' .or. hamiltonian_type.eq.'read') then
! Read Psi_trial
    call read_psi_t
    if(maxval(iwdet_psi_t,1).gt.ndet) stop 'iwdet_psi_t > ndet'

! Diagonalize Hamiltonian
    if(diagonalize_ham.eq.1) then
      call hamiltonian_diagonalize
     else
      diagonal_ham_highest=-1.e99_rk
      diagonal_ham_lowest=1.e99_rk
      do i=1,ndet
        diagonal_ham_highest=max(diagonal_ham_highest,ham(i,i))
        diagonal_ham_lowest=min(diagonal_ham_lowest,ham(i,i))
      enddo
      if(tau.eq.0._rk) then
        tau=tau_multiplier/(diagonal_ham_highest-diagonal_ham_lowest)
!***b
        if(master_core) write(6,'(/,''diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau='',9f12.6)') diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau
       else
        if(master_core) write(6,'(/,''diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, deduced tau='',9f12.6)') diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau_multiplier/(diagonal_ham_highest-diagonal_ham_lowest)
        if(master_core) write(6,'(''Input value (actually used) of tau='',9f12.6)') tau
      endif
      if(master_core) write(6,'(''tau='',f10.5)') tau
!        write(6,'(/,''diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau='',9f12.6)') diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau
!       else
!        write(6,'(/,''diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, deduced tau='',9f12.6)') diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau_multiplier/(diagonal_ham_highest-diagonal_ham_lowest)
!        write(6,'(''Input value (actually used) of tau='',9f12.6)') tau
!      endif
!      write(6,'(''tau='',f10.5)') tau
!***e
    endif
  endif

  end subroutine hamiltonian

! ==============================================================================
  subroutine hamiltonian_fictitious
! ------------------------------------------------------------------------------
! Description   : Create fictitious Hamiltonian matrix at beginning of run
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  implicit none

! local
  integer i,j
  real(rk) :: min_diag, rannyu
! real(rk), allocatable :: work(:),eigval(:),eigvec(:,:),propag(:,:),vec(:),vec_tmp(:)

! Construct hamiltonian matrix
! Diagonal has first element zero and the rest positive
! If bosonic.ge.1 then all the off-diagonal elements are negative and the eigenstate is entirely of one sign.
  do i=2,ndet
    ham(i,i)=spectrum_coef*(i-1)**spectrum_power + (rannyu()-0.5_rk)*ham_diag_fluc
    do j=1,i-1
      if(bosonic.ge.1) then
        ham(i,j)=-abs(rannyu()-0.5_rk)*ham_offdiag_fluc/abs(i-j)**spectrum_power
       else
        ham(i,j)=(rannyu()-0.5_rk)*ham_offdiag_fluc/abs(i-j)**spectrum_power
      endif
      ham(j,i)=ham(i,j)
    enddo
  enddo

! Shift the diagonal so that the smallest element is zero
  min_diag=1.d99
  do i=1,ndet
    min_diag=min(min_diag,ham(i,i))
  enddo
  do i=1,ndet
    ham(i,i)=ham(i,i)-min_diag
  enddo

  if(ipr.ge.1) then
    write(6,*)
    do i=1,ndet
      write(6,'(''ham='',100f10.5)') (ham(i,j),j=1,ndet)
    enddo
  endif

end subroutine hamiltonian_fictitious

! ==============================================================================
  subroutine hamiltonian_diagonalize
! ------------------------------------------------------------------------------
! use read_psi_trial
  use common_psi_t, only: e_trial
  implicit none

! local
  integer i,j,lwork,info,istat, dominant_det
  real(rk) :: tau_optimal_deterministic, tau_optimal_stochastic, dominant_component, top, bot, overlap_with_psi_t, sum_abs_eigvec_comp
! real(rk), allocatable :: work(:),eigval(:),eigvec(:,:),propag(:,:),vec(:),vec_tmp(:)
  real(rk), allocatable :: work(:),eigval(:),eigvec(:,:),propag(:,:)
  integer, allocatable :: index(:)

! Diagonalize Hamiltonian
  write(6,'(''ndet='',i15)') ndet
  call flush(6)
  allocate(eigvec(ndet,ndet),stat=istat)
  if(istat.ne.0) stop 'failed to allocate eigvec'
  allocate(wavefn_exact(ndet),stat=istat)
  if(istat.ne.0) stop 'failed to allocate wavefn_exact'
  eigvec=ham
  lwork=3*ndet-1
  write(6,'(''lwork='',i15)') lwork
  call flush(6)
  allocate(work(lwork),stat=istat)
  if(istat.ne.0) stop 'failed to allocate work'
  allocate(eigval(ndet),stat=istat)
  if(istat.ne.0) stop 'failed to allocate eigval'
  write(6,'(''calling dsyev'')')
  call flush(6)
  call dsyev('V','U',ndet,eigvec,ndet,eigval,work,lwork,info)
  dominant_det=maxloc(abs(eigvec(:,1)),1)
  dominant_component=maxval(abs(eigvec(:,1)),1)
  if(eigvec(dominant_det,1).lt.0._rk) eigvec(:,1)=-eigvec(:,1) ! Flip sign of lowest eigenvector if dominant component is negative
  write(6,'(/,''eigval='',1000f10.5)') eigval
! Check that it is normalized
  write(6,'(''Normalization of lowest eigenvector='',f9.6)') dot_product(eigvec(:,1),eigvec(:,1))
  if(ipr.ge.1) then
    write(6,*)
    do i=1,ndet
      write(6,'(''eigvec='',100f10.5)') (eigvec(i,j),j=1,ndet)
    enddo
  endif
  write(6,'(/,''diagonal of Ham='',10000f10.6)') (ham(i,i),i=1,ndet)
  write(6,'(/,''lowest eigvec  ='',10000f10.6)') (eigvec(i,1),i=1,ndet)
  write(6,'(''Largest components of lowest eigenvector='')',advance="no")
  do i=1,ndet
    if(abs(abs(eigvec(i,1))-dominant_component).lt.1.d-3) write(6,'(i6,f10.6)',advance="no") i,eigvec(i,1)
  enddo
  write(6,*)
! write(6,'(/,''The largest component of this eigenvector is'',i5,'' with value'',f10.5)') dominant_det, dominant_component

! Calculate overlap of exact eigenvector with psi_t
  overlap_with_psi_t=0
  do i=1,ndet_psi_t
    overlap_with_psi_t=overlap_with_psi_t+cdet_psi_t(i)*eigvec(iwdet_psi_t(i),1)
  enddo
  write(6,'(/,''Overlap of psi_t with psi_exact='',f9.6)') overlap_with_psi_t

! Calculate fraction of population on dominant determinant with perfect importance sampling and with no importance sampling
  sum_abs_eigvec_comp=sum(abs(eigvec(:,1)))
  write(6,'(/,''sum_abs_eigvec_comp='',f10.6)') sum_abs_eigvec_comp
  write(6,'(''Fraction of walkers on dominant det with perfect importance sampling and without any='',2f10.6)') dominant_component**2, dominant_component/sum_abs_eigvec_comp

! Set e_trial for use in projector
  write(6,'(''Setting e_trial (from diagonalizing Hamiltonian) ='',f10.5)') eigval(1)
  e_trial=eigval(1)

! Calculate optimal tau's for deterministic propagation and stochastic propagation if there is no sign problem.
  energy_exact=eigval(1)
  !wavefn_exact(:)=eigvec(:,1)*sign(1._rk,eigvec(1,1))
  wavefn_exact(:)=eigvec(:,1)
  write(6,'(/,''Sum of abs wavefn_exact components, ratio of largest component to sum='',10000f10.5)') sum(abs(wavefn_exact(1:ndet))),maxval(wavefn_exact,1)/sum(abs(wavefn_exact(1:ndet)))
  if(tau.eq.0._rk) then
    tau_optimal_deterministic=2/(eigval(ndet)+eigval(2)-2*eigval(1))
    tau_optimal_stochastic=1/(eigval(ndet)-eigval(1))
!   tau=tau_multiplier*tau_optimal_deterministic
    tau=tau_multiplier*tau_optimal_stochastic
    write(6,'(/,''tau_multiplier, tau_optimal_deterministic, tau_optimal_stochastic, tau='',9f10.6)') &
 &  tau_multiplier, tau_optimal_deterministic, tau_optimal_stochastic, tau
  endif
  call flush(6)

! Sort elements of lowest eigvec in order of descending absolute value.  Note this destroys original order.
  allocate(index(ndet),stat=istat)
  if(istat.ne.0) stop 'failed to allocate index'
  do i=1,ndet
    index(i)=i
  enddo
  call shell_sort_eigvec(eigvec(1,1),index,ndet)
  write(6,'(/,''Eigenvector components in order of descending absolute value:'')')
  write(6,'(''index ='',10000i12)') index
  write(6,'(''eigvec='',10000f12.8)') eigvec(:,1)
  deallocate(index)

! deallocate(eigvec)

! Calculate propagator and sign condition number.
! Verify that power method gives smallest eigenvector for small tau and highest eigenvector for large tau
  if(ipr.ge.0) then
    allocate(propag(ndet,ndet),stat=istat)
    if(istat.ne.0) stop 'failed to allocate propag'

! Calculate propagator
    propag=-tau*ham
    do i=1,ndet
      propag(i,i)=propag(i,i)+1+tau*energy_exact
    enddo

    if(ipr.ge.1) then
      write(6,*)
      do i=1,ndet
        write(6,'(''propag='',100f10.5)') (propag(i,j),j=1,ndet)
      enddo
    endif
    write(6,'(/,''diagonal of propag='',10000f10.5)') (propag(i,i),i=1,ndet)
    call flush(6)

  write(6,'(/,''Sum of elements of columns of propag'', 10000f9.4)') (sum(propag(:,i)),i=1,ndet)
  write(6,'(/,''Sum of absolute values of elements of columns of propag'', 10000f9.4)') (sum(abs(propag(:,i))),i=1,ndet)

! Calculate condition number of propagator
    bot=0
    do i=1,ndet
      bot=bot+dot_product(abs(propag(i,:)),abs(wavefn_exact(:)))
    enddo
    top=sum(abs(wavefn_exact(:)))
    write(6,'(''tau, top, bot, top/bot, condition_number(old)='',9f10.6)') tau,top,bot,top/bot,1-(1-top/bot)/tau
    write(6,'(''tau, top, bot, top/bot, condition_number(new)='',9f10.6)') tau,top,bot,top/bot,(top/bot)**(1/tau)

! Warning tmp:
! Calculate highest eigenvector of propagator (should be same as lowest eigenvec of Ham, provided that tau is not too large).
    eigvec=propag
    call dsyev('V','U',ndet,eigvec,ndet,eigval,work,lwork,info)
    dominant_det=maxloc(abs(eigvec(:,ndet)),1)
    dominant_component=maxval(abs(eigvec(:,ndet)),1)
    write(6,'(/,''eigval='',1000f10.5)') eigval
    if(ipr.ge.1) then
      write(6,*)
      do i=1,ndet
        write(6,'(''eigvec='',100f10.5)') (eigvec(i,j),j=1,ndet)
      enddo
    endif
!   write(6,'(/,''lowest eigvec  ='',10000f10.5)') (eigvec(i,1),i=1,ndet)
    write(6,'(/,''highest eigvec  ='',10000f10.5)') (eigvec(i,ndet),i=1,ndet)
    write(6,'(/,''The largest component of highest eigenvector of propagator is'',i5,'' with value'',f10.5)') &
 &  dominant_det, dominant_component
    call flush(6)

! Check power method gives the same eigenvector
!   allocate(vec(ndet),stat=istat)
!   allocate(vec_tmp(ndet),stat=istat)
!   vec=0
!   vec(1)=1
!   vec_tmp=0
!!  Project
!   do j=1,nint(100/tau_multiplier)
!     do i=1,ndet
!       vec_tmp(i)=dot_product(propag(i,:),vec(:))
!     enddo
!     vec=vec_tmp
!   enddo
!   vec=vec/vec(1)
!   if(ipr.ge.1) write(6,'(''vec='',1000f20.5)') vec
!!  write(6,'(''vec='',1000f20.5)') vec/vec(ndet)
!!  Project again to make sure first projection was long enough.
!   do j=1,nint(100/tau_multiplier)
!     do i=1,ndet
!       vec_tmp(i)=dot_product(propag(i,:),vec(:))
!     enddo
!     vec=vec_tmp
!   enddo
!   vec=vec/vec(1)
!   if(ipr.ge.1) write(6,'(''vec='',1000f20.5)') vec
!!  write(6,'(''vec='',1000f20.5)') vec/vec(ndet)
!   call flush(6)

! Calculate psi_0 importance-sampled propagator
    do i=1,ndet
      do j=1,ndet
        propag(i,j)=wavefn_exact(i)*propag(i,j)/wavefn_exact(j)
      enddo
    enddo
    if(ipr.ge.1) then
      write(6,*)
      do i=1,ndet
        write(6,'(''psi_0 imp sampl propag='',100f10.5)') (propag(i,j),j=1,ndet)
      enddo
    endif
    write(6,'(/,''diagonal of psi_0 imp sampl propag='',10000f10.5)') (propag(i,i),i=1,ndet)
    call flush(6)

! Write sum of columns of psi_0 importance-sampled propagator
  write(6,'(/,''Sum of elements of columns of psi_0 imp sampl propag (should be 1 for perfect importance sampling)'')')
  write(6,'(''Sum of elements of columns of psi_0 imp sampl propag'', 10000f9.4)') (sum(propag(:,i)),i=1,ndet)
  write(6,'(/,''Sum of absolute values of elements of columns of psi_0 imp sampl propag (should be 1 for perfect importance &
 & sampling and no sign problem)'')')
  write(6,'(''Sum of absolute values of elements of columns of psi_0 imp sampl propag'', 10000f9.4)') &
 & (sum(abs(propag(:,i))),i=1,ndet)

! Calculate sign condition number of psi_0 importance-sampled propagator
    bot=0
    do i=1,ndet
      bot=bot+dot_product(abs(propag(i,:)),wavefn_exact(:)**2)
    enddo

    top=sum(wavefn_exact(:)**2)

    write(6,'(''tau, top, bot, top/bot, psi_0 imp sampl condition_number(old)='',9f10.6)') tau,top,bot,top/bot,1-(1-top/bot)/tau
    write(6,'(''tau, top, bot, top/bot, psi_0 imp sampl condition_number(new)='',9f10.6)') tau,top,bot,top/bot,(top/bot)**(1/tau)
    call flush(6)

! Warning tmp:
! Use dgeev rather than dsyev since after importance sampling the propagator is no longer symmetric
! eigvec=propag
! call dsyev('V','U',ndet,eigvec,ndet,eigval,work,lwork,info)
!--
! call alloc ('work', work, lwork)
! call dgeev('V','V', orb_tot_nb, mat_a, orb_tot_nb, mat_wr, mat_wi, mat_vl, orb_tot_nb, mat_vr, orb_tot_nb, work, -1, info)
! if (info /= 0) then
!  call die (here, 'problem in dgeev')
! endif
! lwork = work(1)
! call alloc ('work', work, lwork)
! mat_a (:,:) = kappa (:,:)
! call dgeev('V','V', orb_tot_nb, mat_a, orb_tot_nb, mat_wr, mat_wi, mat_vl, orb_tot_nb, mat_vr, orb_tot_nb, work, lwork, info)
!--
! call dgeev('N','V',ndet,eigvec,ndet,eigval,eigvali,vl,1,eigvec,ndet,work,lwork,info)
! if (info /= 0) then
!  call die (here, 'problem in dgeev')
! endif
! lwork = work(1)
! call alloc ('work', work, lwork)
! call dgeev('N','V',ndet,eigvec,ndet,eigval,eigvali,vl,1,eigvec,ndet,work,lwork,info)

!     SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
!     SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
!    $                  LDVR, WORK, LWORK, INFO )

!  JOBVL   (input) CHARACTER*1
!          = 'N': left eigenvectors of A are not computed;
!          = 'V': left eigenvectors of A are computed.
!
!  JOBVR   (input) CHARACTER*1
!          = 'N': right eigenvectors of A are not computed;
!          = 'V': right eigenvectors of A are computed.
!
!  N       (input) INTEGER
!          The order of the matrix A. N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N-by-N matrix A.
!          On exit, A has been overwritten.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  WR      (output) DOUBLE PRECISION array, dimension (N)
!  WI      (output) DOUBLE PRECISION array, dimension (N)
!          WR and WI contain the real and imaginary parts,
!          respectively, of the computed eigenvalues.  Complex
!          conjugate pairs of eigenvalues appear consecutively
!          with the eigenvalue having the positive imaginary part
!          first.
!
!  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
!          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!          after another in the columns of VL, in the same order
!          as their eigenvalues.
!          If JOBVL = 'N', VL is not referenced.
!          If the j-th eigenvalue is real, then u(j) = VL(:,j),
!          the j-th column of VL.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
!          u(j+1) = VL(:,j) - i*VL(:,j+1).
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= 1; if
!          JOBVL = 'V', LDVL >= N.
!
!  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
!          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!          after another in the columns of VR, in the same order
!          as their eigenvalues.
!          If JOBVR = 'N', VR is not referenced.
!          If the j-th eigenvalue is real, then v(j) = VR(:,j),
!          the j-th column of VR.
!          If the j-th and (j+1)-st eigenvalues form a complex
!          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
!          v(j+1) = VR(:,j) - i*VR(:,j+1).
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= 1; if
!          JOBVR = 'V', LDVR >= N.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,3*N), and
!          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
!          performance, LWORK must generally be larger.

! dominant_det=maxloc(abs(eigvec(:,1)))
! dominant_component=maxval(abs(eigvec(:,1)))
! write(6,'(/,''eigval='',1000f10.5)') eigval
! if(ipr.ge.1) then
!   write(6,*)
!   do i=1,ndet
!     write(6,'(''eigvec='',100f10.5)') (eigvec(i,j),j=1,ndet)
!   enddo
! endif
! write(6,'(/,''lowest eigvec  ='',10000f10.5)') (eigvec(i,1),i=1,ndet)
! write(6,'(/,''largest eigvec  ='',10000f10.5)') (eigvec(i,ndet),i=1,ndet)
! write(6,'(/,''The largest component of this eigenvector is'',i5,'' with value'',f10.5)') dominant_det, dominant_component

!   vec=0
!   vec(1)=1
!   vec_tmp=0
!!  Project
!   do j=1,nint(100/tau_multiplier)
!     do i=1,ndet
!       vec_tmp(i)=dot_product(propag(i,:),vec(:))
!     enddo
!     vec=vec_tmp
!   enddo
!   vec=vec/vec(1)
!   if(ipr.ge.1) write(6,'(''vec='',1000f20.5)') vec
!!  write(6,'(''vec='',1000f20.5)') vec/vec(ndet)
!!  Project again to make sure first projection was long enough.
!   do j=1,nint(100/tau_multiplier)
!     do i=1,ndet
!       vec_tmp(i)=dot_product(propag(i,:),vec(:))
!     enddo
!     vec=vec_tmp
!   enddo
!   vec=vec/vec(1)
!   if(ipr.ge.1) write(6,'(''vec='',1000f20.5)') vec
!!  write(6,'(''vec='',1000f20.5)') vec/vec(ndet)
!   call flush(6)

! Calculate propagator again
    propag=-tau*ham
    do i=1,ndet
      propag(i,i)=propag(i,i)+1+tau*energy_exact
    enddo

! Calculate psi_g importance-sampled propagator
    do i=1,ndet
      do j=1,ndet
        propag(i,j)=psi_g(i)*propag(i,j)/psi_g(j)
      enddo
    enddo
    if(ipr.ge.1) then
      write(6,*)
      do i=1,ndet
        write(6,'(''psi_g imp sampl propag='',100f10.5)') (propag(i,j),j=1,ndet)
      enddo
    endif
    write(6,'(/,''diagonal of psi_g imp sampl propag='',10000f10.5)') (propag(i,i),i=1,ndet)
    call flush(6)

! Write sum of columns of psi_g importance-sampled propagator
  write(6,'(/,''Sum of elements of columns of psi_g imp sampl propag (should be 1 for perfect importance sampling)'')')
  write(6,'(''Sum of elements of columns of psi_g imp sampl propag'', 10000f9.4)') (sum(propag(:,i)),i=1,ndet)
  write(6,'(/,''Sum of absolute values of elements of columns of psi_g imp sampl propag (should be 1 for perfect importance &
 & sampling and no sign problem)'')')
  write(6,'(''Sum of absolute values of elements of columns of psi_g imp sampl propag'', 10000f9.4)') &
 &(sum(abs(propag(:,i))),i=1,ndet)

! Calculate sign condition number of psi_g importance-sampled propagator
    bot=0
    do i=1,ndet
      bot=bot+dot_product(abs(propag(i,:)),abs(wavefn_exact(:)*psi_g(:)))
    enddo

    top=sum(abs(wavefn_exact(:)*psi_g(:)))

    write(6,'(''tau, top, bot, top/bot, psi_g imp sampl condition_number(old)='',9f10.6)') tau,top,bot,top/bot,1-(1-top/bot)/tau
    write(6,'(''tau, top, bot, top/bot, psi_g imp sampl condition_number(new)='',9f10.6)') tau,top,bot,top/bot,(top/bot)**(1/tau)
    call flush(6)

!   deallocate(vec) ; deallocate(vec_tmp)
    deallocate(eigvec); deallocate(propag)
  endif ! ipr

  deallocate(work) ; deallocate(eigval)
  call my_second(2,'hamiltonian_diagonalize')
  call flush(6)

  end subroutine hamiltonian_diagonalize

! ==============================================================================
  real(rk) function fact(n)
! ------------------------------------------------------------------------------
! Description   : Factorial function
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  integer i,n
  fact = 1
  do i=1,n
    fact = fact*i
  enddo
  end function fact

! ===================================================================================================
    subroutine binomial(n1,n2,b)
! ----------------------------------------------------------------------------------------------------
! Description   : Calculate binomial/combinatorial function
! Author        : H.J. Changlani
! ----------------------------------------------------------------------------------------------------
    implicit none
    integer n1,n2,b
    integer i,denom

    if (n1.lt.n2) then
        b=0
        return
    endif

    b=1
    denom=1

    do i=1,n2
        b=b*(n1-n2+i)
    enddo

    do i=1,n2
        denom=denom*i
    enddo

    b=b/denom

    return

    end subroutine binomial

! ===================================================================================================
    integer function signmod(a,b)
! ----------------------------------------------------------------------------------------------------
! Description   : Get signed mod
! Author        : H.J. Changlani
! ----------------------------------------------------------------------------------------------------

        implicit none
        integer a,b
        integer tmp1,tmp2

        if (b.eq. 0 ) then
            signmod=0
            return
        endif

        tmp1=modulo(a,b)
        tmp2=modulo(-a,b)

        signmod=min(tmp1,tmp2)
        return
    end function

! ===================================================================================================
  subroutine make_hubbard_matrix_2d(l_x,l_y, pbc, n_alpha, n_beta, t, U, total_configs, ham, ipr)
! ----------------------------------------------------------------------------------------------------
! Description   : Generate full Hamiltonian for 1D lattice given filling
! Author        : H.J. Changlani
! ----------------------------------------------------------------------------------------------------
    implicit none

    !integer, parameter ::        rk = kind(1._rk)
    logical, intent(in) ::       pbc
    integer, intent(in) ::       l_x,l_y,n_alpha,n_beta,ipr, total_configs
    real(rk), intent(in) ::      t,U
    integer                      i,j,k,l,istat,tmp,tmp2,num_zeros,num_ones,flg,na,nb
    integer                      num_alpha_configs,num_beta_configs,non_zero
    real(rk)                     mat_elt
    integer, allocatable ::      alpha_config(:),beta_config(:),config(:)
    integer, allocatable ::      alpha_configs(:,:),beta_configs(:,:),all_configs(:,:)
    real(rk)                     ham(total_configs,total_configs)

!   integer lwork, info, dominant_det(1)
!   real(rk), allocatable :: work(:), eigval(:), eigvec(:,:), dominant_component

    if(ipr.ge.6) print *,"Calling Binomial"
    call binomial(l_x*l_y,n_alpha,num_alpha_configs)
    call binomial(l_x*l_y,n_beta,num_beta_configs)

    if(ipr.ge.6) print *,"Allocating alpha, beta and config arrays"
    allocate(alpha_config(l_x*l_y),stat=istat)
    allocate(beta_config(l_x*l_y), stat=istat)
    allocate(config(2*l_x*l_y),    stat=istat)

    if(ipr.ge.6) print *,"Allocating all alpha, beta and config arrays"
    allocate(alpha_configs(num_alpha_configs,l_x*l_y),stat=istat)
    allocate(beta_configs(num_beta_configs,l_x*l_y), stat=istat)
    allocate(all_configs(num_alpha_configs*num_beta_configs,2*l_x*l_y), stat=istat)

    na=0
    do i=1,2**(l_x*l_y)

        tmp=i-1
        num_zeros=0
        num_ones=0
        flg=0

        do j=1,(l_x*l_y)
            tmp2=mod(tmp,2)
            if(ipr.ge.6) print *,"i=",i
            if(ipr.ge.6) print *,"tmp2=",tmp2
            alpha_config(l_x*l_y+1-j)=tmp2
            tmp=tmp/2

            if (tmp2 .eq. 0) then
                num_zeros=num_zeros+1
            else
                num_ones=num_ones+1
            endif

            if ((num_zeros .gt. (l_x*l_y-n_alpha)) .or. (num_ones .gt. n_alpha) ) then
                flg=1
                exit
            endif

       enddo

        if (flg.eq.0) then
            na=na+1
            if(ipr.ge.6) print *,"na=",na
            alpha_configs(na,:)=alpha_config(:)
            if(ipr.ge.6) print *,alpha_config
        endif

    enddo

    nb=0
    do i=1,2**(l_x*l_y)

        tmp=i-1
        num_zeros=0
        num_ones=0
        flg=0

        do j=1,(l_x*l_y)
            tmp2=mod(tmp,2)
            beta_config(l_x*l_y+1-j)=tmp2
            tmp=tmp/2

            if (tmp2 .eq. 0) then
                num_zeros=num_zeros+1
            else
                num_ones=num_ones+1
            endif

            if ((num_zeros .gt. (l_x*l_y-n_beta)) .or. (num_ones .gt. n_beta) ) then
                flg=1
                exit
            endif

       enddo

        if (flg.eq.0) then
            nb=nb+1
            if(ipr.ge.6) print *,"nb=",nb
            beta_configs(nb,:)=beta_config(:)
        endif

    enddo


    if(ipr.ge.6) print *,"Fusing alpha and beta configs"

    do i=1,num_alpha_configs
        do j=1,l_x*l_y
            config((2*j)-1)=alpha_configs(i,j)
        enddo

        do k=1,num_beta_configs
            do l=1,l_x*l_y
                config(2*l)=beta_configs(k,l)
            enddo

            all_configs(((i-1)*num_beta_configs)+k,:)=config
            if(ipr.ge.2) write(6,'(i6,2x,1000i2)') ((i-1)*num_beta_configs)+k,config
        enddo

    enddo

    print *,""
    print *,""

    print *,"Size(all_configs,1)=",size(all_configs,1)
    print *,"Size(all_configs,2)=",size(all_configs,2)

    do i=1,size(all_configs,1)
        non_zero=0
        do j=i,size(all_configs,1)
            call hubbard_matrix_element_2d(l_x,l_y,pbc,all_configs(i,:),all_configs(j,:),t,U,mat_elt,ipr)
            if(ipr.ge.6) print *,all_configs(i,:)
            if(ipr.ge.6) print *,all_configs(j,:)
            ham(i,j)=mat_elt
            ham(j,i)=mat_elt
            if (mat_elt .ne. 0) then
                if(ipr.ge.6) print *,all_configs(i,:)
                if(ipr.ge.6) print *,all_configs(j,:)
                if(ipr.ge.6) print *,"Matrix element = ",mat_elt
                if(ipr.ge.6) print *,""
                if(ipr.ge.6) print *,""
                if(ipr.ge.6) print *,""
                non_zero=non_zero+1
            endif
        enddo
        if(ipr.ge.2) then
          write(6,'(/,''Number of non zero elements='',i5,'' for config'')') non_zero
          write(6,'(1000i2)') all_configs(i,:)
        endif
    enddo

!   allocate(eigvec(total_configs,total_configs),stat=istat)
!   if(istat.ne.0) stop 'failed to allocate eigvec'
!   allocate(eigval(total_configs),stat=istat)
!   if(istat.ne.0) stop 'failed to allocate eigval'
!   lwork=3*total_configs-1
!   allocate(work(lwork),stat=istat)
!   if(istat.ne.0) stop 'failed to allocate work'

!   eigvec=ham
!   call dsyev('V','U',total_configs,eigvec,total_configs,eigval,work,lwork,info)
!   dominant_det=maxloc(abs(eigvec(:,1)))
!   dominant_component=maxval(abs(eigvec(:,1)))
!   write(6,'(/,''eigval='',1000f10.5)') eigval
!   if(ipr.ge.1) then
!     write(6,*)
!     do i=1,total_configs
!       write(6,'(''eigvec='',100f10.5)') (eigvec(i,j),j=1,total_configs)
!     enddo
!   endif
!   write(6,'(/,''diagonal of Ham='',10000f10.5)') (ham(i,i),i=1,total_configs)
!   write(6,'(/,''lowest eigvec  ='',10000f10.5)') (eigvec(i,1),i=1,total_configs)
!   write(6,'(/,''lowest eigvec  ='',10000d12.4)') (eigvec(i,1),i=1,total_configs)
!   write(6,'(/,''The largest component of this eigenvector is'',i5,'' with value'',f10.5)') dominant_det, dominant_component

!   deallocate(eigvec) ; deallocate(eigval) ; deallocate(work)

    deallocate(alpha_config) ; deallocate(beta_config) ; deallocate(config)
    deallocate(alpha_configs) ; deallocate(beta_configs) ; deallocate(all_configs)

  end subroutine make_hubbard_matrix_2d


! ===================================================================================================
  subroutine hubbard_matrix_element_2d(l_x,l_y, pbc, bra_config, ket_config, t, U, mat_elt, ipr)
! ----------------------------------------------------------------------------------------------------
! Description   : Generate Hamiltonian matrix element for 2D lattice given bra and ket configurations
! Note          : Fermion signs are important
! Author        : H.J. Changlani
! ----------------------------------------------------------------------------------------------------

! Example of 2D numbering for 4x4 lattice

! Odd site = alpha (up electron)  Even site = beta (down electron)

! 25  26  27  28  29  30  31  32
! 17  18  19  20  21  22  23  24
! 9   10  11  12  13  14  15  16
! 1   2   3   4   5   6   7    8

! FORTRAN notation starts numbering with 1 instead of 0 (C++)

! I use 0,1 to specify occupancy for each "site"
! Rule for matrix element : The bra_config and the ket_config must be different in exactly
!                           2 positions or no positions. Moreover positions must be nearest
!                           neighbors (I have accounted for pbc vs obc as well)

    implicit none
    !integer, parameter ::  rk = kind(1._rk)
    integer, intent(in) :: l_x,l_y,ipr
    integer, intent(in) :: bra_config(2*l_x*l_y),ket_config(2*l_x*l_y)
    real(rk), intent(in) :: t,U
    real(rk), intent(inout) :: mat_elt
    logical, intent(in) :: pbc
    integer     i,pos_1,pos_2,diff,flg,doubly_occupied,num_between,config_size,x_1,x_2,y_1,y_2
    real(rk)  fermion_sign
    logical  allowed

    config_size=l_x*l_y*2
    diff=0
    flg=0
    doubly_occupied=0

    do i=1,config_size
        if (mod(i,2).ne.0) then
            if (ket_config(i).eq.1 .and. ket_config(i+1).eq.1) then
                doubly_occupied=doubly_occupied+1;
            endif
        endif

        if (ket_config(i).ne.bra_config(i)) then
            diff=diff+1;
            if(ipr.ge.6) print *,"Diff = ",diff
            if (diff.eq.1) then
                pos_1=i
            elseif (diff.eq.2) then
                pos_2=i
                if ((mod((pos_2-pos_1),2) .ne. 0) .or. ((ket_config(pos_1)+ket_config(pos_2)) .ne. 1) ) then
                    flg=1
                    exit
                endif
            else
                flg=1
                exit
            endif
        endif
    enddo

! By construction pos_1<pos_2

if(ipr.ge.6) print *,"Flag = ",flg

! If 1 differing position or disallowed
    if (flg.eq.1 .or. diff.eq.1 ) then
        mat_elt=0
        return
    endif
! If both the configs are the same then U * number of doubly occupied sites
    if (diff.eq.0) then
        mat_elt=U*real(doubly_occupied)
        return
    endif

! If the configs are different in exactly 2 sites

    if (diff.eq.2) then

        ! Getting 2-D coordinates
        if(ipr.ge.6) print *,"Pos 1 = ",pos_1
        y_1=((pos_1-1)/(2*l_x))+1
        x_1=pos_1-((y_1-1)*2*l_x)

        if(ipr.ge.6) print *,"x 1 = ",x_1
        if(ipr.ge.6) print *,"y 1 = ",y_1

        if(ipr.ge.6) print *,"Pos 2 = ",pos_2

        y_2=((pos_2-1)/(2*l_x))+1
        x_2=pos_2-((y_2-1)*2*l_x)

        if(ipr.ge.6) print *,"x 2 = ",x_2
        if(ipr.ge.6) print *,"y 2 = ",y_2

        ! Checking neighbors on rectangular lattice

        if (pbc.eqv. .false.) then
            if ((x_1.eq. x_2) .and. abs(y_1-y_2) .eq. 1) then
                allowed=.true.
            elseif ((y_1.eq. y_2) .and. abs(x_1-x_2) .eq. 2) then
                allowed=.true.
            else
                allowed=.false.
            endif
        else
            if ((x_1.eq. x_2) .and. (signmod(y_1-y_2,l_y) .eq. 1) ) then
                allowed=.true.
            elseif ((y_1.eq. y_2) .and. (signmod(x_1-x_2,2*l_x) .eq. 2)) then
                allowed=.true.
            else
                allowed=.false.
            endif
        endif

        num_between=0
        if (allowed .eqv. .true.) then
            do i=pos_1+2,pos_2-2,2
               num_between=num_between+ket_config(i)
            enddo

            if (mod(num_between,2) .eq. 0) then
                fermion_sign=1.0
            else
                fermion_sign=-1.0
            endif

            mat_elt=t*fermion_sign
            return

        else
            mat_elt=0
            return

        endif
    endif

end subroutine hubbard_matrix_element_2d

  subroutine shell_sort_eigvec(eigvec,index,ndet)
! ==============================================================================
! Description   : Shell-Metzger sort abs(eigvec) in descending order and give corresponding index.
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: ndet
  integer, intent(inout) :: index(ndet)
  real(rk), intent(inout) :: eigvec(ndet)
  integer it,nn,lognb2,m,k,j,i,l
  real(rk) t
  lognb2=int(log(real(ndet,rk))/log(2._rk)+10._rk**(-14))
  m=ndet
  do 20 nn=1,lognb2
    m=m/2
    k=ndet-m
    do 20 j=1,k
      do 10 i=j,1,-m
        l=i+m
        if (abs(eigvec(l)).lt.abs(eigvec(i))) goto 20
        t=eigvec(i)
        eigvec(i)=eigvec(l)
        eigvec(l)=t
        it=index(i)
        index(i)=index(l)
        index(l)=it
10    continue
20  continue
  return
  end subroutine shell_sort_eigvec

subroutine read_psit
  ! Read in a trial wave function
  ! Same format as that of the wave functions that are dumped out when psit_out is true
  ! A Holmes, 22 Jul 2013

  use common_psi_t, only : use_psit_in,use_psit_out,use_psit_con_in,use_psit_con_out,print_psit_wo_sqmc,psit_in_file,psit_out_file,psit_con_in_file,psit_con_out_file,ndet_psi_t_in,dets_up_psi_t_in,dets_dn_psi_t_in,cdet_psi_t_in,use_elems_in,use_elems_out,dtm_elems_in_file,dtm_elems_out_file
  use common_ham, only : n_core_orb
  use mpi_routines, only : master_core,whoami,mpi_bsend
  integer :: iostatus,i,j
  integer,allocatable :: tmp_reader(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec) :: filled_core
#else
  integer(ik) :: filled_core
#endif


  if (master_core) then
    read(5,*,iostat=iostatus) use_psit_in,use_psit_out,print_psit_wo_sqmc
    if (.not.iostatus<0) then
      if (use_psit_in.and.use_psit_out) then
        read(5,*) psit_in_file,psit_out_file
      elseif (use_psit_in.and..not.use_psit_out) then
        read(5,*) psit_in_file
      elseif (.not.use_psit_in.and.use_psit_out) then
        read(5,*) psit_out_file
      endif
    endif
  endif

  call mpi_bsend(use_psit_in)
  call mpi_bsend(use_psit_out)
  call mpi_bsend(print_psit_wo_sqmc)
  call mpi_bsend(psit_in_file)
  call mpi_bsend(psit_out_file)

  write(6,'(''run_type='',a)') run_type
  if (master_core .and. run_type.ne.'hci') then
    if (use_psit_in) then
      write (6,*) "Using input psi trial read from file ",psit_in_file,"whoami",whoami; call flush(6)
    else
      !write (6,*) "NOT using input psi trial, whoami=",whoami; call flush(6)
      write (6,*) "NOT using input psi trial"; call flush(6)
    endif
    if (use_psit_out) then
      write (6,*) "Dumping generated psi trial into file ",psit_out_file,"whoami",whoami; call flush(6)
    else
      write (6,*) "NOT dumping generated psi trial"; call flush(6)
    endif
    write (6,*) ""

    read(5,*,iostat=iostatus) use_psit_con_in,use_psit_con_out
    if (.not.iostatus<0) then
      if (use_psit_con_in.and.use_psit_con_out) then
        read(5,*) psit_con_in_file,psit_con_out_file
      elseif (use_psit_con_in.and..not.use_psit_con_out) then
        read(5,*) psit_con_in_file
      elseif (.not.use_psit_con_in.and.use_psit_con_out) then
        read(5,*) psit_con_out_file
      endif
    endif
    read(5,*,iostat=iostatus) use_elems_in,use_elems_out
    if (.not.iostatus<0) then
      if (use_elems_in.and.use_elems_out) then
        read(5,*) dtm_elems_in_file,dtm_elems_out_file
      elseif (use_elems_in.and..not.use_elems_out) then
        read(5,*) dtm_elems_in_file
      elseif (.not.use_elems_in.and.use_elems_out) then
        read(5,*) dtm_elems_out_file
      endif
    endif
  endif

  call mpi_bsend(use_psit_con_in)
  call mpi_bsend(use_psit_con_out)
  call mpi_bsend(psit_con_in_file)
  call mpi_bsend(psit_con_out_file)

  if (master_core) then
    if (use_psit_con_in) then
      write (6,*) "Using input psi trial connections read from file ",psit_con_in_file,"whoami",whoami; call flush(6)
    else
      write (6,*) "NOT using input psi trial connections"; call flush(6)
    endif
    if (use_psit_con_out) then
      write (6,*) "Dumping generated psi trial connections into file ",psit_con_out_file,"whoami",whoami; call flush(6)
    else
      write (6,*) "NOT dumping generated psi trial connections",whoami; call flush(6)
    endif
    write (6,*) ""
  endif

  call mpi_bsend(use_elems_in)
  call mpi_bsend(use_elems_out)
  call mpi_bsend(dtm_elems_in_file)
  call mpi_bsend(dtm_elems_out_file)

  if (master_core) then
    if (use_elems_in) then
      write (6,*) "Using input deterministic matrix elements read from file ",dtm_elems_in_file,"whoami",whoami; call flush(6)
    else
      write (6,*) "NOT using input deterministic matrix elements"; call flush(6)
    endif
    if (use_elems_out) then
      write (6,*) "Dumping generated deterministic matrix elements into file ",dtm_elems_out_file,"whoami",whoami; call flush(6)
    else
      write (6,*) "NOT dumping generated deterministic matrix elements"; call flush(6)
    endif
    write (6,*) ""
    if (print_psit_wo_sqmc .and. run_type.ne.'hci') then
      write (6,'(''NOT performing sqmc run'')') ; call flush(6)
    else
      write (6,'(''Performing sqmc run, process='',i6)') whoami ; call flush(6)
    endif

    close(5)
  endif

! Read in connections to Psi_T
  if (use_psit_in) then
    if (master_core) then
      open(55,file=psit_in_file,status='old')
      read(55,*) ndet_psi_t_in
      read(55,*)
      allocate(dets_up_psi_t_in(ndet_psi_t_in))
      allocate(dets_dn_psi_t_in(ndet_psi_t_in))
      allocate(cdet_psi_t_in(ndet_psi_t_in))
      allocate(tmp_reader(nup+ndn-2*n_core_orb))
      filled_core = 0_ik
      if (n_core_orb>0) then
        do i=1,n_core_orb
          filled_core = ibset(filled_core,i-1)
        enddo
      endif
      do i=1,ndet_psi_t_in
        read(55,*) tmp_reader ,cdet_psi_t_in(i)
        dets_up_psi_t_in(i) = filled_core
        do j=1,nup-n_core_orb
          dets_up_psi_t_in(i) = ibset(dets_up_psi_t_in(i),tmp_reader(j)+n_core_orb-1)
        enddo
        dets_dn_psi_t_in(i) = filled_core
        do j=1,ndn-n_core_orb
          dets_dn_psi_t_in(i) = ibset(dets_dn_psi_t_in(i),tmp_reader(nup-n_core_orb+j)+n_core_orb-1)
        enddo
      enddo
      close(55)
      write (6,*) "Read in a trial wave function composed of",ndet_psi_t_in,"determinants."
      call my_second(2,'reading Psi_T'); call flush(6)
    endif
    call mpi_bsend(ndet_psi_t_in)
    if (.not. master_core) then
      allocate(dets_up_psi_t_in(ndet_psi_t_in))
      allocate(dets_dn_psi_t_in(ndet_psi_t_in))
      allocate(cdet_psi_t_in(ndet_psi_t_in))
    endif
    call mpi_bsend(dets_up_psi_t_in)
    call mpi_bsend(dets_dn_psi_t_in)
    call mpi_bsend(cdet_psi_t_in)
  elseif (.not.use_psit_con_in) then
    write (6,*) "No input trial wave function given. Generating trial wave function iteratively starting from HF."
    ndet_psi_t_in = 0 ! flag for no input psi trial
  endif

end subroutine read_psit

end module hamiltonian_mod
