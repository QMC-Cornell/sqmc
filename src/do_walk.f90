module do_walk

  use types, only: i1b, i2b, i4b, i8b, i16b, ik, ik_vec, rk
#ifdef NUM_ORBITALS_GT_127
  use overload
#endif
  use common_walk, only: MWALK, nwalk, walk_det, walk_dets_up, walk_dets_dn, walk_wt, walk_wt_fn, e_num_walker,e_den_walker, n_connected_dets, m_connected_dets, connected_dets, connected_dets_sign, matrix_elements, imp_distance
  use common_psi_t, only: e_trial,ndet_psi_t_connected,psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den,trial_wf_iters,norb_trial_wf
  use common_run, only: run_type, proposal_method, ipr, partial_node_eps,importance_sampling, connected_dets_up, connected_dets_dn,use_efficient_heatbath, input_copy_unit
  use common_ham, only: norb
  use common_imp, only: norb_imp,n_imp_initiators,n_imp_truncate, n_imp, minus_tau_H_indices, minus_tau_H_nonzero_elements, minus_tau_H_values, semistochastic, imp_up, imp_dn, imp_iters, diff_from_psi_t, size_deterministic
  use hubbard, only: l_x, l_y
  use more_tools,only: print_real_matrix,gen_hist,add_to_hist
  use tools, only : do_n_times

!**** Edited by AR[7/23/13]: added module for parallelization
  use mpi_routines, only: snd_cnt, mpi_push_nwalk, mpi_bsend, mpi_bsend_diffroot, master_core, whoami, ncores, mpi_allred, init_snd_buffers, mpi_sendwalks, init_hash_owners, get_det_owner, mpi_red, mpi_sendnewwalks, mpi_barr, mpi_redscatt_real_dparray, mpi_stop, mpi_gath, mpi_red_max, mpi_red_min
!**** end of edit

  implicit none

! r_initiator = Almost the same as n_a in Cleland-Booth-Alavi-2010
!             = -1, No initiator
!             =  0, Initiator is not used, except that the permanent initiator must be occupied with at least one walker of the initial sign.
!             >  0, Initiator used, and the permanent initiator must be occupied with at least one walker of the initial sign.
! n_permanent_initiator = The number of permanent initiators.  The permanent initiators are presently chosen to be those dets. that have the largest absolute coeffs.
! sign_permanent_initiator = Sign of the permanent initiator coeffs. in the trial wavefn.  The trial wavefn. typically has more dets. than just the permanent initiators.
! initiator() = 3 permanent initiator provided r_initiator >= 0
!             = 2 initiator
!             = 1 spawned from initiator, or spawned from noninitiator and its weight exceeds the initiator threshold.  At the next generation it will become an initiator if its weight still exceeds threshold
!             = 0 spawned from noninitiator and its weight is less than the initiator threshold.  It is killed before it can contribute to expectation values.
! One difference from Cleland-Booth-Alavi-2010 is that it takes a determinant 2 MC steps to go from initiator=0 to 2, or from 2 to 0.
! imp_distance() = 0 Determinant is in deterministic space
!                =-1 Determinant was spawned by a determinant in the deterministic space
!                =-2 Determinant is connected to Psi_T.  When hf_to_psit=true, it is treated semideterministically, i.e. moves to and from 1st state are deterministic but all other moves are stochastic
!                > 0 Number of moves since last visit to deterministic space.  Note this can be less than number of moves needed to get to this det from deterministic space.
!
! run_type    =  no_fixed_node (exact, except for initiator)
!             =  none          (exact, except for initiator)
!             =  fixed_node1 (usual discrete-space FN)
!             =  fixed_node2
!             =  fixed_node3
!             =  fixed_node4
!             =  partial_node
!             =  release_node
!             =  selected_ci
!             =  trunc_lanc
!             =  sr          (stochastic reconfiguration)
!             =  vmc
!             =  hci (Heat-bath configuration interaction (includes perturbation theory))

! hamiltonian_type = read (read Hamiltonian matrix)
!                  = fictitious
!                  = hubbard2 (real-space Hubbard)
!                  = hubbardk (momentum-space Hubbard)
!                  = hubbarddm (Hitesh's density matrix from a patch idea for Hubbard)
!                  = heg (homogeneous electron gas in 2 or 3D)
!                  = chem (atoms, molecules, ... (read either integrals.dat or FCIDUMP))

  integer, parameter :: M_PERMANENT_INITIATOR=10000, M_WARNING=1000
  integer w_abs_gen_begin, w_abs_gen_target, nstep, nblk, nblk_eq, n_permanent_initiator, n_warning, initiator_power, initiator_min_distance, reached_w_abs_gen
  integer, dimension(4,2) :: irand_seed
  real(rk) :: population_control_exponent, e_trial_initial, e_est, e_est_gen, reweight_factor_inv, min_wt, r_initiator
  real(rk) :: H_off_multiplier
  integer(i1b), allocatable :: initiator(:), sign_permanent_initiator(:)
  integer, allocatable ::iorder(:)                 ! When key is sorted, iorder is sorted with it and then used to order all auxilliary arrays.
  integer, allocatable ::temp_i_1(:), temp_i_2(:)  ! For walk_det and iorder. Needed for merge_sort
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable ::temp_i16(:)           ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2
  type(ik_vec), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#else
  integer(ik), allocatable ::temp_i16(:)           ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2
  integer(ik), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#endif
! Warning: in above 2 lines we are temporarily declaring both temp_i16 and temp_i16_up, temp_i16_dn, but we really need only one or the other.
  real(rk), allocatable :: cumulative_probability(:),cumulative_probability_fn(:),diag_elems(:)
  real(rk) :: always_spawn_cutoff_wt ! always spawn if abs wt > always_spawn_cutoff_wt
  logical :: c_t_initiator ! true if all states in C(T) are permanent initiators; false otherwise
  integer :: leave_out,num_to_move ! for insertion sort of newly spawned walkers into list outside of ct
  integer,allocatable :: empty_outside_ct(:),indices_new_dets(:),leave_out_list(:) ! for merge_original_with_spawned3
  ! Edited by AAH on 31 Jan 2013. Make iblk, istep global variables, so that we don't have to keep passing them into subroutines.
  integer :: istep,iblk,iblkk
  ! End of edit
  real(rk) cdet_psi_t_abs_max !** Added by AR, variable needed even in 'walk' routine
  real(rk) :: log_n_con
  logical :: use_exponential_projector
  integer :: nbins
  real(rk),allocatable :: lbounds(:)
  integer(i8b),allocatable :: bins(:)
  integer(i8b),allocatable :: bins_single(:)
  integer(i8b),allocatable :: bins_double(:)
  real(rk),allocatable :: lbounds_hf(:)
  integer(i8b),allocatable :: bins_hf(:)
  real(rk),allocatable :: lbounds_new(:)
  integer(i8b),allocatable :: bins_new(:)
  real(rk),allocatable :: lbounds_sing(:)
  integer(i8b),allocatable :: bins_sing(:)
  real(rk),allocatable :: lbounds_doub(:)
  integer(i8b),allocatable :: bins_doub(:)
  integer,allocatable :: occ_up(:),occ_dn(:)
  real(rk) :: max_wt_ratio, min_wt_ratio, all_max_wt_ratio, all_min_wt_ratio
  integer max_wt_ratio_excit_level(1), recvcount
  integer,allocatable :: all_max_wt_ratio_excit_level(:)
  logical :: off_diag_histograms = .true.


!-------------------------------------------------------------------------------------
!*** Added by AR [27/11/13]
! using procedure pointers in order to avoid repeated string comparisons
    abstract interface
      subroutine off_diagonal_move_template(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j, importance_sampling,proposal_prob_inverse)
        use types, only : ik,ik_vec,rk
#ifdef NUM_ORBITALS_GT_127
          type(ik_vec), intent(in) :: det_i_up, det_i_dn
          type(ik_vec), intent(out) :: det_j_up, det_j_dn
#else
          integer(ik), intent(in) :: det_i_up, det_i_dn
          integer(ik), intent(out) :: det_j_up, det_j_dn
#endif
          real(rk), intent(out) :: weight_j
          integer,optional,intent(in) :: importance_sampling
          real(rk),optional,intent(out) :: proposal_prob_inverse !needed for hubbard_k
      end subroutine off_diagonal_move_template
    end interface

    procedure (off_diagonal_move_template), pointer :: off_diagonal_move => null ()

    abstract interface
        subroutine move_template(iwalk)
          integer, intent(in) :: iwalk
        end subroutine move_template
    end interface

    procedure (move_template), pointer :: move => null ()
!-------------------------------------------------------------------------------------



  contains

  subroutine prepare_namelist_copy
    implicit none

    integer :: input_size
    character(:), allocatable :: input_str
    character(1024) :: line
    integer :: io_state

    open(unit=input_copy_unit, status='scratch')
    do
      read(5, '(A)', iostat=io_state) line
      if (io_state == 0 .and. trim(line).ne.'end') then
        write(input_copy_unit, '(A)') line
      else
        rewind(5)
        return
      endif
    enddo
  end subroutine prepare_namelist_copy

! ==============================================================================
  subroutine read_input
! ------------------------------------------------------------------------------
! Description   : Create Hamiltonian matrix at beginning of run
! Author        : Cyrus Umrigar
! Comments      : Edited by HJC on May 13-14,2011. Please look for bypasses
!                 for hubbard2 code. Make sure it is not interfering with other
!                 code. Added option of run_type = vmc, fixed_node1, fixed_node2,
!                 release_node, or no_fixed_node/none
! ------------------------------------------------------------------------------

  use hamiltonian_mod
  use read_psi_trial
  use common_psi_t, only: ndet_psi_t, dets_up_psi_t, dets_dn_psi_t, cdet_psi_t,psi_g_energy,psi_g_epsilon,hf_to_psit, ndet_psi_t_in,dets_up_psi_t_in,dets_dn_psi_t_in,cdet_psi_t_in,use_psit_con_in,psit_con_in_file,psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den,use_elems_in,use_elems_out,print_psit_wo_sqmc,dtm_elems_in_file,dtm_elems_out_file,dtm_energy
  use common_run, only: reweight_factor_inv_max_multiplier, reweight_factor_inv_max, initiator_rescale_power,connected_dets_up,connected_dets_dn,connected_matrix_elements
  use common_ham, only: hamiltonian_type,nup,ndn
  use chemistry, only : counter_two_body,hamiltonian_chem, time_sym, hamiltonian_chem_time_sym, find_connected_dets_chem,setup_orb_by_symm, get_ind
  use hubbard,only: hubbard_ipr,hamiltonian_hubbard_k,hamiltonian_hubbard,space_sym,hamiltonian_hubbard_k_space_sym, find_connected_dets_hubbard_k
  use heg, only : hamiltonian_heg, find_connected_dets_heg
  use semistoch, only : generate_psi_t_connected_e_loc,generate_space_iterate,perform_selected_ci,perform_truncated_lanczos
  use more_tools, only : binary_search
  use types, only : num_words, bits_per_word
  use mpi_routines, only : mpi_snd, master_core, mpi_bsend, mpi_barr, whoami,ncores
  use common_selected_ci, only : hf_det, get_auto_hf, hf_symmetry, up, dn, irreps, irrep_occs_up, irrep_occs_dn, n_irrep, lz, g, u, natorb, get_natorbs, use_pt, greens_function, get_greens_function, n_w, w_min, w_max, &
      & selected_ci, eps_var_sched, eps_var, eps_pt, eps_pt_big, target_error, dump_wf_var, n_max_connections, n_states, n_mc, eps_pt_big_energy, &
      & active_space, n_var_e_up, n_var_e_dn, n_var_orbs, var_orbs, use_hash_generation, &
      & n_energy_batch

! local
! integer istat,iwalk,idet,i,i_psit,j,k
  integer istat,iwalk,idet,i,i_psit,j
  real(rk) cdet_psi_t_abs_sum,matrix_element
  !type(ik_vec),dimension(207184) :: final_up,final_dn
  !real(rk),dimension(207184) :: final_coeff
  !real(rk) :: abs_wt,abs_wt_tot
! integer, allocatable :: iorder_imp(:)
  integer :: nnzero
  integer,allocatable :: tmp_reader(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec) :: filled_core,tmp_up,tmp_dn
  type(ik_vec),allocatable :: my_up(:),my_dn(:)
  type(ik_vec) :: tmp_det
#else
  integer(ik) :: filled_core,tmp_up,tmp_dn
  integer(ik),allocatable :: my_up(:),my_dn(:)
  integer(ik) :: tmp_det
#endif
  real(rk) :: tmp_num,tmp_den
  real(rk),allocatable :: my_num(:),my_den(:)
  integer :: target_core,coreid
  character(len=16) :: fmt
  integer,allocatable :: n(:),n_old(:)
  integer :: ibatch,n_batches,istart,iend,ind
  integer :: read_batch_size_per_core
  integer :: ndet_psit_tmp
  integer,allocatable :: orb_up(:),orb_dn(:)
  integer :: n_nonzero_elems_upper_triangular,read_batch_size
  integer :: counter
  character(len=120) :: buf
  integer :: io_state
  integer :: a,b,old,new

  if(master_core) then
     !MJO FIXME prepare_namelist_copy does not work with MPI
     !          enabled, due to the rewind on stdin call
     call prepare_namelist_copy()
  endif
! read and set 2 random number seeds.
! The first for the the Alias method in SHCI and for the fictitious Hamiltonian in SQMC
! The second for the SQMC walk
  if(master_core) read(5,'(4i4,x,4i4)') irand_seed
  call mpi_bsend(irand_seed)
  irand_seed(:,1) = irand_seed(:,1)! ! MJO This seed should be the same for all cores for the stochastic pt alias routine because every core draws the same sample and keeps what belongs to it
  irand_seed(:,2) = irand_seed(:,2) + whoami ! gives every core a different seed for the SQMC walk
  write(6,'(/,''random number seeds [master_core]'',t35,4i4,x,4i4)') irand_seed

  ! Set random number seed for SHCI
  call setrn(irand_seed(1,1))

  if(master_core) read(5,*) run_type ! run_type can be vmc, fixed_node1, fixed_node2, fixed_node3, partial_node, release_node, no_fixed_node, sr, SR, none, trunc_lanc, or selected_ci, hci, but some have not been implemented yet.
  call mpi_bsend(run_type)
  write(6,'(/,''run_type= '',a)') trim(run_type)

  ! You must enter one of the run_types listed below
  if((run_type .ne.'fixed_node1') .and. (run_type .ne. 'fixed_node2') .and. (run_type .ne. 'fixed_node3') .and. (run_type .ne. 'fixed_node4') .and. (run_type .ne. 'partial_node') .and. (run_type .ne. 'release_node') .and. (run_type .ne. 'vmc') .and. (run_type .ne. 'no_fixed_node') .and. (run_type .ne. 'none') .and. (run_type .ne. 'selected_ci') .and. (run_type .ne. 'trunc_lanc') .and. (run_type .ne. 'sr') .and. (run_type .ne. 'SR') .and. (run_type.ne.'cisdtq') .and. (run_type.ne.'cisd_pt') .and. (run_type.ne.'hci')) then
    stop 'I could not recognize the run_type option (vmc,fixed_node1,fixed_node2,fixed_node3,partial_node,release_node,no_fixed_node/none,sr,SR,hci) I am ending the program.'
  endif

  if(run_type.ne.'hci') then

    if(master_core) read(5,*) nstep,nblk,nblk_eq,ipr
    call mpi_bsend(nstep)
    call mpi_bsend(nblk)
    call mpi_bsend(nblk_eq)
    call mpi_bsend(ipr)
    write(6, '(''nstep, nblk, nblk_eq, ipr='',2i8,9i5)') nstep, nblk, nblk_eq, ipr
    call flush(6)

    hubbard_ipr=ipr

    if(master_core) read(5,*) w_abs_gen_begin,w_abs_gen_target, MWALK ! Set starting and target numbers of walkers
    call mpi_bsend(w_abs_gen_begin)
    call mpi_bsend(w_abs_gen_target)
    call mpi_bsend(MWALK)
    write(6,'(''w_abs_gen_begin, w_abs_gen_target, MWALK='',i8,i9,i10)') w_abs_gen_begin, w_abs_gen_target, MWALK

    if(w_abs_gen_begin.gt.w_abs_gen_target .and. w_abs_gen_target.ge.0) then ! w_abs_gen_target<0 is used to stop the run before the MC steps.
      write(6,'(''w_abs_gen_begin must be <= w_abs_gen_target'')')
      stop 'w_abs_gen_begin must be <= w_abs_gen_target'
    endif

    if(master_core) read(5,*) tau_multiplier, tau ! tau_multiplier is almost always < 0.5, typically 0.025 to 0.2 (used only if tau=0).
    call mpi_bsend(tau_multiplier)
    call mpi_bsend(tau)
    write(6,'(''tau_multiplier, tau='',9f13.8)') tau_multiplier, tau
    if(master_core) read(5,*) reweight_factor_inv_max_multiplier, reweight_factor_inv_max ! 2. is a good choice for reweight_factor_inv_max_multiplier (used only if reweight_factor_inv_max_mult=0).
    call mpi_bsend(reweight_factor_inv_max_multiplier)
    call mpi_bsend(reweight_factor_inv_max)
    write(6,'(''reweight_factor_inv_max_multiplier, reweight_factor_inv_max='',9f10.5)') reweight_factor_inv_max_multiplier, reweight_factor_inv_max
    if(master_core) read(5,*) population_control_exponent, e_trial_initial, min_wt
    call mpi_bsend(population_control_exponent)
    call mpi_bsend(e_trial_initial)
    call mpi_bsend(min_wt)
    write(6,'(''population_control_exponent, e_trial_initial, min_wt='',2f11.5,f6.2)') population_control_exponent, e_trial_initial, min_wt
!   proposal method can be uniform, uniform2, CauchySchwarz, heat_bath, heat_bath2 or heat_bath3
!   0=no imp sampl, 1=imp. sampling.
!   initiator_power=0 => fixed initiator, initiator_power=1 => linearly increasing  initiator
!   initiator_rescale_power is used to rescale r_initiator during equilibration when both tau and r_initiator are made larger than their final values
    if(master_core) read(5,*) proposal_method, importance_sampling, r_initiator, initiator_power, initiator_rescale_power !, H_off_multiplier
    call mpi_bsend(proposal_method)
    if (proposal_method=='fast_heatbath') then
      use_efficient_heatbath = .true.
    else
      use_efficient_heatbath = .false.
    endif
    call mpi_bsend(use_efficient_heatbath)
    call mpi_bsend(importance_sampling)
    call mpi_bsend(r_initiator)
    call mpi_bsend(initiator_power)
    call mpi_bsend(initiator_rescale_power)

    if(proposal_method.ne.'uniform' .and. proposal_method.ne.'uniform2' .and. proposal_method.ne.'CauchySchwarz' .and. proposal_method .ne. 'fast_heatbath' .and. proposal_method.ne.'heat_bath' .and. proposal_method.ne.'heat_bath2' .and. proposal_method.ne.'heat_bath3') then
      call mpi_stop('proposal_method must be one of uniform, uniform2, CauchySchwarz, fast_heatbath, heat_bath, heat_bath2, heat_bath3')
    endif

    initiator_min_distance = 0

    write(6,'(''proposal_method, importance_sampling, r_initiator, initiator_power, initiator_min_distance, initiator_rescale_power= '',a,i3,f6.2,2i3,f6.3)') trim(proposal_method), importance_sampling, r_initiator, initiator_power, initiator_min_distance, initiator_rescale_power

    if(index(run_type,'fixed_node').ne.0 .and. r_initiator.gt.0) then
      write(6,'(/,''fixed_node runs must have r_initiator=0'')')
      stop 'fixed_node runs must have r_initiator=0'
    endif

    if (run_type.eq.'cisdtq') then
      if (master_core)  read(5,*) eps_var
    elseif (run_type.eq.'cisd_pt') then
      if (master_core)  read(5,*) eps_pt
    endif

    partial_node_eps=0._rk
    if ((run_type .eq. 'partial_node') .or. (run_type .eq. 'sr')) then
      if(master_core) read (5,*) partial_node_eps                                  ! partial node eps = 0 is the true projector
      call mpi_bsend(partial_node_eps)
      write(6,'(/,''partial_node_eps ='',f6.2)') partial_node_eps    ! partial node eps = 1 is the fixed node projector
    endif
    if (run_type .eq. 'sr') then
      if (importance_sampling .eq. 0) then
          write(6,'(/,''Stochastic reconfiguration does not work without importance sampling'')')
          stop
      endif
    endif

    if(master_core) read(5,*) semistochastic,use_exponential_projector ! semi-stochastic moves (true or false)
    call mpi_bsend(semistochastic)
    call mpi_bsend(use_exponential_projector)
    if(semistochastic) then
      write(6,'(/,''semistochastic run'')')
      if(master_core) read(5,*) diff_from_psi_t !different from psi trial?
      call mpi_bsend(diff_from_psi_t)
      if (diff_from_psi_t) then
        if(master_core) read(5,*) imp_iters
         call mpi_bsend(imp_iters)
         allocate(norb_imp(imp_iters))
         allocate(n_imp_initiators(imp_iters))
         allocate(n_imp_truncate(imp_iters))
         if(master_core) then
           read(5,*) norb_imp
           read(5,*) n_imp_initiators
           read(5,*) n_imp_truncate
         endif
         call mpi_bsend(norb_imp)
         call mpi_bsend(n_imp_initiators)
         call mpi_bsend(n_imp_truncate)

        if(maxval(norb_imp).gt.norb) then
          write(6,'(''norb_imp must be <= norb; norb_imp, norb='',9i5)') norb_imp, norb
          stop 'norb_imp must be <= norb'
        endif
      else
        write (6,'(''Generating trial wave function and important space simultaneously'')')
        if(master_core) read(5,*) size_deterministic
        call mpi_bsend(size_deterministic)
      endif
    else
      write(6,'(/,''not semistochastic run'')')
    endif

    if (master_core) then
      if (use_exponential_projector) then
        write (6,'(''Using Exponential Projector'')')
      else
        write (6,'(''Using Linear Projector'')')
      endif
      call flush(6)
    endif

    if(master_core .and. semistochastic) read(5,*) hf_to_psit, c_t_initiator, always_spawn_cutoff_wt ! replace HF with psi_T for 1st state?, psi_T connections (C(T)) are permanent initiators?, always spawn if abs wt > always_spawn_cutoff_wt
    call mpi_bsend(hf_to_psit)
    call mpi_bsend(c_t_initiator)
    call mpi_bsend(always_spawn_cutoff_wt)

    if (hf_to_psit) then
      write (6,*) "Replacing HF state with trial wave function"
    else
      write (6,*) "NOT replacing HF state with trial wave function"
    endif

    if (c_t_initiator) then
      write (6,*) "Connections of Psi_T are permanent initiators"
    else
      write (6,*) "Connections of Psi_T are NOT permanent initiators"
    endif

    write(6,'(''Cutoff weight for always proposing off-diagonal move='',es10.2)') always_spawn_cutoff_wt

    if (importance_sampling==1) then
      if(master_core) read(5,*) psi_g_energy,psi_g_epsilon
      call mpi_bsend(psi_g_energy)
      call mpi_bsend(psi_g_epsilon)
      write (6,*) "psi_g_energy=",psi_g_energy,"psi_g_epsilon=",psi_g_epsilon
    endif

    write(6,'(''energy_exact='',9f10.5)') energy_exact ; call flush(6)

  elseif (run_type.eq.'hci') then

    nup=1000
    ndn=1000
    allocate(up(nup))
    allocate(dn(ndn))
    n_irrep=1000
    allocate(irreps(n_irrep))
    allocate(irrep_occs_up(n_irrep))
    allocate(irrep_occs_dn(n_irrep))
    allocate(var_orbs(1000))
    if (master_core) then
      read(5,*) eps_var, eps_pt, target_error, n_states
      read(5,*) dump_wf_var
      write(6,'(''eps_var, eps_pt, target_error, n_states='',3es8.1,i3)') eps_var, eps_pt, target_error, n_states
      write(6,'(''dump_wf_var='',l)') dump_wf_var

      rewind(input_copy_unit)
      read(unit=input_copy_unit, nml=selected_ci, iostat=io_state)
      eps_var_sched=max(eps_var_sched,eps_var) ! Set those elements of eps_var_sched that are not set in the input
     !write(6,*) ; write(6,nml=selected_ci)
      write(6,'(/,''&SELECTED_CI'',/,''eps_var_sched='',30es9.2,/, &
      &''eps_var, eps_pt, eps_pt_big, target_error='',4es10.2,/, &
      &''dump_wf_var='',l,'' n_max_connections='',es12.4,'' n_states='',i3, '' n_mc='',i6,/,&
      &''n_energy_batch='',i2,'' use_hash_generation='',l)') &
      &eps_var_sched, eps_var, eps_pt, eps_pt_big, target_error, dump_wf_var, n_max_connections, n_states, n_mc, n_energy_batch, use_hash_generation

#ifdef MAC
      !MJO We need to error out if we are on a mac and use the flags
      ! for automatic memory, since mem_avail does not work on MacOS
      if (eps_pt_big.le.0 .or. n_mc.le.0) then
        write(6,*) "ERROR! MacOS requires explicit eps_pt_big and n_mc parameters"
        stop "ERROR! MacOS requires explicit eps_pt_big and n_mc parameters"
      endif
#endif
      ! Automatically choose HF det
      rewind(input_copy_unit)
      get_auto_hf = .false.
      hf_symmetry = 999
      up(:) = 0
      dn(:) = 0
      irrep_occs_up(:) = 0
      irrep_occs_dn(:) = 0
      read(unit=input_copy_unit, nml=hf_det, iostat=io_state)
      do i=1,nup
        if (up(i)==0) then
          nup = i-1
          exit
        endif
        if(i.eq.nup) stop 'up not dimensioned large enough'
      enddo
      do i=1,ndn
        if (dn(i)==0) then
          ndn = i-1
          exit
        endif
        if(i.eq.ndn) stop 'dn not dimensioned large enough'
      enddo
      do i=1,n_irrep
        if (irreps(i)==0) then
          n_irrep = i-1
          exit
        endif
        if(i.eq.n_irrep) stop 'n_irrep not dimensioned large enough'
      enddo
      if (minval(irreps(1:n_irrep))<0) then ! there are negative indices, so assume that the indices given are Sandeep's. This converts them to my indices.
        do i=1,n_irrep
          if (irreps(i).ne.1.and.irreps(i).ne.2) then
            old = irreps(i)
            a = abs(old)/2
            b = (abs(old)+1)/2
            new = a+3*b-8
            if (old<0)  new=new+2
            irreps(i) = new
          endif
        enddo
      endif
      write(6,'(/,''&HF_DET'')')
      if (hf_symmetry==999.and.(g.or.u)) then
        if (g) then
          write (6,*) "Lz=",lz,"g"
          call get_ind(lz,0,hf_symmetry)
        elseif (u) then
          write (6,*) "Lz=",lz,"u"
          call get_ind(lz,1,hf_symmetry)
        endif
      endif
      write(6,'(''hf_symmetry='',i2)') hf_symmetry
      write (6,*) "HF up electron occupancies (from input file)=",up(1:nup)
      write (6,*) "HF dn electron occupancies (from input file)=",dn(1:ndn)
      if (hf_symmetry.ne.999) then
        get_auto_hf = .true.
        write (6,*) "Automatically choosing HF det of total symmetry",hf_symmetry; call flush(6)
      elseif (up(1)>0) then
        write (6,*) "Starting with user-specified HF det"; call flush(6)
      endif

      ! HCI Variational Active Space
      rewind(input_copy_unit)
      var_orbs(:) = 0
      read(unit=input_copy_unit, nml=active_space, iostat=io_state)
      if (var_orbs(1)>0) then
        do i=1,n_var_orbs
          if (var_orbs(i)==0) then
            n_var_orbs = i-1
            exit
          endif
          if(i.eq.n_var_orbs) stop 'var_orbs not dimensioned large enough'
        enddo
      endif
      write (6,*) "&ACTIVE_SPACE"
      write (6,*) "n_var_e_up, n_var_e_dn, n_var_orbs=",n_var_e_up,n_var_e_dn,n_var_orbs
     !write (6,*) "var_orbs=",var_orbs(1:n_var_orbs)

      ! HCI Natural Orbitals
      rewind(input_copy_unit)
      read(unit=input_copy_unit, nml=natorb, iostat=io_state) ! HCI Natural Orbitals
     !write(6,*) ; write(6,nml=natorb)
      write(6,'(/,''&NATORB'',/,''get_natorbs, use_pt='',2l)') get_natorbs, use_pt

      if(get_natorbs .and. ncores.gt.1) then
        write(6,'(''At present FCIDUMP for NOs should be done on 1 core to avoid excessive memory'')')
!        stop 'At present FCIDUMP for NOs should be done on 1 core to avoid excessive memory'
      endif

      ! HCI Green's Function
      rewind(input_copy_unit)
      read(unit=input_copy_unit, nml=greens_function, iostat=io_state)
     !write(6,*) ; write(6,nml=greens_function)
      write(6,'(/,''&GREENS_FUNCTION'',/,''get_greens_function, n_w, w_min, w_max='',l,i4,9es11.4)') get_greens_function, n_w, w_min, w_max

    endif ! master_core

    ! Broadcast variables
    call mpi_bsend(eps_var_sched)
    call mpi_bsend(eps_var)
    call mpi_bsend(eps_pt)
    call mpi_bsend(eps_pt_big)
    call mpi_bsend(target_error)
    call mpi_bsend(dump_wf_var)
    call mpi_bsend(n_max_connections)
    call mpi_bsend(n_states)
    call mpi_bsend(n_mc)
    call mpi_bsend(eps_pt_big_energy)
    call mpi_bsend(n_energy_batch)
    call mpi_bsend(use_hash_generation)

    call mpi_bsend(get_auto_hf)
    call mpi_bsend(hf_symmetry)
    call mpi_bsend(lz)
    call mpi_bsend(g)
    call mpi_bsend(u)
    call mpi_bsend(up)
    call mpi_bsend(dn)

    call mpi_bsend(n_var_e_up)
    call mpi_bsend(n_var_e_dn)
    call mpi_bsend(n_var_orbs)
    call mpi_bsend(var_orbs)

    call mpi_bsend(get_natorbs)
    call mpi_bsend(use_pt)

    call mpi_bsend(get_greens_function)
    call mpi_bsend(n_w)
    call mpi_bsend(w_min)
    call mpi_bsend(w_max)

  endif ! run_type == hci

  if (ik==i1b) then
    write (6,*) "Using 8 bit integers to encode orbital occupations"
  elseif (ik==i2b) then
    write (6,*) "Using 16 bit integers to encode orbital occupations"
  elseif (ik==i4b) then
    write (6,*) "Using 32 bit integers to encode orbital occupations"
  elseif (ik==i8b) then
    write (6,*) "Using 64 bit integers to encode orbital occupations"
  elseif (ik==i16b) then
    write (6,*) "Using 128 bit integers to encode orbital occupations"
  else
    stop "Error: unknown integer type! Set ik to one of i1b, i2b, i4b, i8b, or i16b in types.f90"
  endif

  write (6,*) "num_words=",num_words,"bits_per_word=",bits_per_word,"max number of orbitals=",num_words*bits_per_word

! Set random number seed for the Hamiltonian (if it is random)
  call setrn(irand_seed(1,1))

 ! Read parameters specific to system, set up Hamiltonian matrix, possibly diagonalize it, find m_connected_dets (the maximum number of connected dets), do hci if it is an hci run.
  call hamiltonian

! Use and #ifndef here rather than an if(run_type.eq.'hci') to avoid compiling the rest of the file.
#ifndef SQMC
  stop
  end
#else

  allocate(connected_dets(m_connected_dets),stat=istat)
  if(istat.ne.0) stop 'failed to allocate connected_dets'
  allocate(connected_dets_sign(m_connected_dets),stat=istat)
  if(istat.ne.0) stop 'failed to allocate connected_dets_sign'
  allocate(cumulative_probability(m_connected_dets),stat=istat)
  if(istat.ne.0) stop 'failed to allocate cumulative_probability'
  if (run_type .eq. 'sr') then
    allocate(cumulative_probability_fn(m_connected_dets),stat=istat)
    if(istat.ne.0) stop 'failed to allocate cumulative_probability_fn'
  endif
  write (6,*) "Done with Hamiltonian setup..."

!*** added by AR [27/11/2013] for using interfaces
! call init_off_diagonal(hamiltonian_type)
  call init_move
!*** end added by AR

  if (run_type .eq. 'sr') then
    if (hamiltonian_type .ne. 'hubbard2') then
        write(6,'(''Stochastic reconfiguration only works for the hubbard2 hamiltonian'')')
        stop
    endif
  endif

  call flush(6)
  if(e_trial.eq. 0._rk) then
      write(6,'(''In do_walk.f90: Setting e_trial (from input) ='',f10.5)') e_trial_initial
      call flush(6)
    e_trial=e_trial_initial
  endif

  if ((hamiltonian_type .ne. 'hubbard2') .and. (hamiltonian_type .ne. 'heg')) then ! Hack by Hitesh - 28 th Aug 2012
      if(run_type.eq.'chemistry') then
        if(maxval(norb_trial_wf).gt.norb) then
            write(6,'(''norb_trial_wf must be <= norb; norb_trial_wf, norb='',9i5)') norb_trial_wf, norb
            call flush(6)
        endif
        stop 'norb_imp must be <= norb'
      endif
  endif

! Set random number seed for the walk
  call setrn(irand_seed(1,2))

  !***initialize owners hash table in mod_cluster.f90
  call init_hash_owners
!*** end Added by AR
  ! Set MWALK
  if (semistochastic) then
    if (hf_to_psit) then
      MWALK=int(3*(w_abs_gen_target/min_wt+ndet_psi_t_connected))/ncores
      !MWALK=max(MWALK,int(3*(w_abs_gen_target/min_wt+ndet_psi_t_connected)))
      write (6,*) "w_abs_gen_target,min_wt,ndet_psi_t_connected=",w_abs_gen_target,min_wt,ndet_psi_t_connected; call flush(6)
      write(6,'(''Setting MWALK=3*(w_abs_gen_target/min_wt+ndet_psi_t_connected)/ncores='',i10,'' ='',es9.2)') MWALK, float(MWALK)
      !write(6,'(''Setting MWALK=3*(w_abs_gen_target/min_wt+ndet_psi_t_connected)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
    else
      !MWALK=int(4*(w_abs_gen_target/min_wt+n_imp))/ncores
      if (use_efficient_heatbath) then
        MWALK=max(MWALK,int(3.5*(w_abs_gen_target/min_wt+n_imp)))
        write(6,'(''1Setting MWALK=3.5*(w_abs_gen_target/min_wt+n_imp)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
      else
        MWALK=max(MWALK,int(4*(w_abs_gen_target/min_wt+n_imp)))
        write(6,'(''1Setting MWALK=4*(w_abs_gen_target/min_wt+n_imp)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
      endif
    endif
  else
    !MWALK=int(4*(w_abs_gen_target/min_wt))/ncores
    MWALK=max(MWALK,int(4*(w_abs_gen_target/min_wt))/ncores)
    write(6,'(''2Setting MWALK=4*(w_abs_gen_target/min_wt)/ncores='',i10,'' ='',es9.2)') MWALK, float(MWALK)
    !write(6,'(''2Setting MWALK=4*(w_abs_gen_target/min_wt)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
  endif

  !***initialize send buffers in mod_cluster.f9C
  if(ncores>1) then
   !call init_snd_buffers(int(1.5*MWALK/ncores,8)) !we need less space for spawning dets
    call init_snd_buffers(int(0.75*MWALK,8)) !we need less space for spawning dets
  else
    write (6,*) "initializing send buffers; MWALK=",MWALK; call flush(6)
    call init_snd_buffers(MWALK)
  endif

  ! Edited by AAH on 20 Jan 2012
  ! Hack by Hitesh - generate_psi_t will not work for hubbard2, heg
  if (hamiltonian_type .ne. 'hubbard2' .and. hamiltonian_type .ne. 'heg' .and. hamiltonian_type .ne. 'fictitious' .and. hamiltonian_type .ne. 'read') then
    if (use_psit_con_in) then
      read_batch_size_per_core = 1000000/ncores
     !read_batch_size_per_core = 10000000/ncores ! This works, but the above is slightly faster
      write (6,*) "Reading in the local energies from the file ",psit_con_in_file
      call my_second(1,'Read local energies'); call flush(6)
      allocate(tmp_reader(nup+ndn-2*n_core_orb))
      filled_core = 0_ik
      if (n_core_orb>0) then
        do i=1,n_core_orb
          filled_core = ibset(filled_core,i-1)
        enddo
      endif

      if (master_core) then
        open(56,file=psit_con_in_file,status='old')
        read(56,*) ndet_psi_t,ndet_psi_t_connected
        read(56,*)
      endif
      call mpi_bsend(ndet_psi_t)
      call mpi_bsend(ndet_psi_t_connected)

      if (master_core) then
        allocate(dets_up_psi_t(ndet_psi_t))
        allocate(dets_dn_psi_t(ndet_psi_t))
        allocate(cdet_psi_t(ndet_psi_t))
      endif

      ! Arrays on each core that will be reallocated at end
      allocate(psi_t_connected_dets_up(int(1.5*ndet_psi_t_connected/ncores)))
      allocate(psi_t_connected_dets_dn(int(1.5*ndet_psi_t_connected/ncores)))
      allocate(psi_t_connected_e_loc_num(int(1.5*ndet_psi_t_connected/ncores)))
      allocate(psi_t_connected_e_loc_den(int(1.5*ndet_psi_t_connected/ncores)))

      ! Decide how many batches we need
      n_batches=(ndet_psi_t_connected-1)/(read_batch_size_per_core)+1
     !n_batches=(ndet_psi_t_connected-1)/(ncores*read_batch_size_per_core)+1

      ! Temporary buffer arrays
      allocate(my_up(int(1.5*ncores*read_batch_size_per_core)))
      allocate(my_dn(int(1.5*ncores*read_batch_size_per_core)))
      allocate(my_num(int(1.5*ncores*read_batch_size_per_core)))
      allocate(my_den(int(1.5*ncores*read_batch_size_per_core)))

      write (6,*) "About to read in psit of length",ndet_psi_t,"with",ndet_psi_t_connected,"connections."; call flush(6)
      allocate(n(ncores))
      allocate(n_old(ncores))
      ndet_psit_tmp=0
      n(:)=0
      n_old(:)=0
      call mpi_bsend(n_old)
      do ibatch = 1,n_batches
        if (master_core) then
          do i=(ibatch-1)*read_batch_size_per_core+1,min(ndet_psi_t_connected,ibatch*read_batch_size_per_core)
            read(56,*) tmp_reader, tmp_num, tmp_den
            tmp_up = filled_core
            do j=1,nup-n_core_orb
              tmp_up = ibset(tmp_up,tmp_reader(j)+n_core_orb-1)
            enddo
            tmp_dn = filled_core
            do j=1,ndn-n_core_orb
              tmp_dn = ibset(tmp_dn,tmp_reader(nup-n_core_orb+j)+n_core_orb-1)
            enddo
            if (abs(tmp_den)>1.e-12_rk) then
              ndet_psit_tmp=ndet_psit_tmp+1
              dets_up_psi_t(ndet_psit_tmp) = tmp_up
              dets_dn_psi_t(ndet_psit_tmp) = tmp_dn
              cdet_psi_t(ndet_psit_tmp) = tmp_den
            endif
            ! Now divide up the connected dets among the processors
            coreid=get_det_owner(tmp_up,tmp_dn)
            do target_core=0,ncores-1
              if (coreid==target_core) then
                ! add to list
                n(coreid+1) = n(coreid+1)+1
                ind = coreid*read_batch_size_per_core+n(coreid+1)-n_old(coreid+1)
                my_up(ind) = tmp_up
                my_dn(ind) = tmp_dn
                my_num(ind) = tmp_num
                my_den(ind) = tmp_den
              endif
            enddo
          enddo
        endif ! master_core
        call mpi_bsend(n)

        ! Now send them to respective cores
        do target_core=0,ncores-1
          istart = target_core*read_batch_size_per_core+1!n_old(target_core+1)+1
          iend = target_core*read_batch_size_per_core+n(target_core+1)-n_old(target_core+1)
#ifdef MPI
          if (target_core .ne. 0) then
            call mpi_snd(my_up(istart:iend),target_core)
            call mpi_snd(my_dn(istart:iend),target_core)
            call mpi_snd(my_num(istart:iend),target_core)
            call mpi_snd(my_den(istart:iend),target_core)
          endif
#endif
          ! Now move stuff from my_ arrays to psit_con arrays
          if (whoami==target_core) then
            psi_t_connected_dets_up(n_old(target_core+1)+1:n(target_core+1))=my_up(istart:iend)
            psi_t_connected_dets_dn(n_old(target_core+1)+1:n(target_core+1))=my_dn(istart:iend)
            psi_t_connected_e_loc_num(n_old(target_core+1)+1:n(target_core+1))=my_num(istart:iend)
            psi_t_connected_e_loc_den(n_old(target_core+1)+1:n(target_core+1))=my_den(istart:iend)
          endif
        enddo
        n_old = n
      enddo

      if (master_core)  close(56)
      call my_second(2,'reading connections to Psi_T'); call flush(6)
      if (.not.master_core) then
       !deallocate(dets_up_psi_t,dets_dn_psi_t,cdet_psi_t)
        allocate(dets_up_psi_t(ndet_psi_t))
        allocate(dets_dn_psi_t(ndet_psi_t))
        allocate(cdet_psi_t(ndet_psi_t))
      endif
      ! Send all psi trial dets to all processors
      call mpi_bsend(dets_up_psi_t)
      call mpi_bsend(dets_dn_psi_t)
      call mpi_bsend(cdet_psi_t)

      ! Finally reallocate
      if (size(my_up)<n(whoami+1)) then
        deallocate(my_up,my_dn,my_num,my_den)
        allocate(my_up(n(whoami+1)))
        allocate(my_dn(n(whoami+1)))
        allocate(my_num(n(whoami+1)))
        allocate(my_den(n(whoami+1)))
      endif
      my_up(1:n(whoami+1))=psi_t_connected_dets_up(1:n(whoami+1))
      my_dn(1:n(whoami+1))=psi_t_connected_dets_dn(1:n(whoami+1))
      my_num(1:n(whoami+1))=psi_t_connected_e_loc_num(1:n(whoami+1))
      my_den(1:n(whoami+1))=psi_t_connected_e_loc_den(1:n(whoami+1))
      deallocate(psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den)
      allocate(psi_t_connected_dets_up(n(whoami+1)))
      allocate(psi_t_connected_dets_dn(n(whoami+1)))
      allocate(psi_t_connected_e_loc_num(n(whoami+1)))
      allocate(psi_t_connected_e_loc_den(n(whoami+1)))
      psi_t_connected_dets_up(1:n(whoami+1))=my_up(1:n(whoami+1))
      psi_t_connected_dets_dn(1:n(whoami+1))=my_dn(1:n(whoami+1))
      psi_t_connected_e_loc_num(1:n(whoami+1))=my_num(1:n(whoami+1))
      psi_t_connected_e_loc_den(1:n(whoami+1))=my_den(1:n(whoami+1))
      deallocate(my_up,my_dn,my_num,my_den)

      write (6,*) "Finished reading",ndet_psi_t_connected,"determinants connected to trial wavefunction."
      call my_second(2,'Read local energies'); call flush(6)

      write (6,'(''Number of saved local energies, ndet_psi_t_connected='',i10,'' ='',es11.4)') ndet_psi_t_connected, real(ndet_psi_t_connected)
      call flush(6)

      write (6,*) "First 1000 (sorted by determinant label):"
      write (6,*) "  Det up     Det dn      E_Num        E_Den           E_L"
      write (fmt, '(i2)') 2*num_words
      do i=1,min(size(psi_t_connected_e_loc_num),1000)
        write(6,'(' // trim(fmt) // 'i10 f15.5 f15.5 f15.5 i3)') psi_t_connected_dets_up(i),psi_t_connected_dets_dn(i),psi_t_connected_e_loc_num(i),psi_t_connected_e_loc_den(i),psi_t_connected_e_loc_num(i)/psi_t_connected_e_loc_den(i), whoami
      enddo

    else ! use_psit_con_in
      ! Set MWALK
      if (semistochastic) then
        if (hf_to_psit) then
          !MWALK=int(3*(w_abs_gen_target/min_wt+ndet_psi_t_connected))/ncores
          MWALK=max(MWALK,int(3*(w_abs_gen_target/min_wt+ndet_psi_t_connected))/ncores)
          !write (6,*) "w_abs_gen_target,min_wt,ndet_psi_t_connected=",w_abs_gen_target,min_wt,ndet_psi_t_connected; call flush(6)
          write(6,'(''Setting MWALK=3*(w_abs_gen_target/min_wt+ndet_psi_t_connected)/ncores='',i10,'' ='',es9.2)') MWALK, float(MWALK)
          !write(6,'(''Setting MWALK=3*(w_abs_gen_target/min_wt+ndet_psi_t_connected)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
        else
          !MWALK=int(4*(w_abs_gen_target/min_wt+n_imp))/ncores
          MWALK=max(MWALK,int(4*(w_abs_gen_target/min_wt+n_imp))/ncores)
          write(6,'(''3Setting MWALK=4*(w_abs_gen_target/min_wt+n_imp)/ncores='',i10,'' ='',es9.2)') MWALK, float(MWALK)
          !write(6,'(''3Setting MWALK=4*(w_abs_gen_target/min_wt+n_imp)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
        endif
      else
        !MWALK=int(4*(w_abs_gen_target/min_wt))/ncores
        MWALK=max(MWALK,int(4*(w_abs_gen_target/min_wt))/ncores)
        write(6,'(''4Setting MWALK=4*(w_abs_gen_target/min_wt)/ncores='',i10,'' ='',es9.2)') MWALK, float(MWALK)
        !write(6,'(''4Setting MWALK=4*(w_abs_gen_target/min_wt)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
      endif

      call mpi_bsend(MWALK)
      !***initialize send buffers in mod_cluster.f90
     !if(ncores>1) then
     !  call init_snd_buffers(int(0.75*MWALK,8)) !we need less space for spawning dets
     !else
     !  call init_snd_buffers(MWALK)
     !endif
      call generate_psi_t_connected_e_loc(ndet_psi_t_connected,psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den)
      ! TODO: Send connections out in this case too
    endif
  else
      write (6,*) "Heg or hubbard2 or fictitious chosen, so generate_psi_t_eloc not being called"
      call flush(6)
  endif ! hamiltonian_type

  if(semistochastic) then
    if (use_elems_in) then

      read_batch_size = 10000000
      write (6,*) "Reading in the matrix elements from the file ",dtm_elems_in_file
      call my_second(1,'read deterministic matrix elements'); call flush(6)
      if (allocated(tmp_reader))  deallocate(tmp_reader)
      allocate(tmp_reader(nup+ndn-2*n_core_orb))
      filled_core = 0_ik
      if (n_core_orb>0) then
        do i=1,n_core_orb
          filled_core = ibset(filled_core,i-1)
        enddo
      endif

      if (master_core) then
        open(57,file=dtm_elems_in_file,status='old')
        read(57,*) n_imp,n_nonzero_elems_upper_triangular
      endif
      call mpi_bsend(n_imp)
      call mpi_bsend(n_nonzero_elems_upper_triangular)

      allocate(minus_tau_H_nonzero_elements(n_imp))
      allocate(minus_tau_H_values(n_nonzero_elems_upper_triangular))
      allocate(minus_tau_H_indices(n_nonzero_elems_upper_triangular))
      allocate(imp_up(n_imp))
      allocate(imp_dn(n_imp))

      if (master_core)  read(57,*) minus_tau_H_nonzero_elements(1:n_imp)

      ! Decide how many batches we need
      n_batches=(n_nonzero_elems_upper_triangular-1)/(read_batch_size)+1

      if (master_core) then
        write (6,*) "About to read in deterministic space of size",n_imp,"with",n_nonzero_elems_upper_triangular,"nonzero elements in its upper triangular form."; call flush(6)
        do ibatch = 1,n_batches
          do i=(ibatch-1)*read_batch_size+1,min(n_imp,ibatch*read_batch_size)
            read(57,*) ind, tmp_reader
            tmp_up = filled_core
            do j=1,nup-n_core_orb
              tmp_up = ibset(tmp_up,tmp_reader(j)+n_core_orb-1)
            enddo
            tmp_dn = filled_core
            do j=1,ndn-n_core_orb
              tmp_dn = ibset(tmp_dn,tmp_reader(nup-n_core_orb+j)+n_core_orb-1)
            enddo
            imp_up(ind) = tmp_up
            imp_dn(ind) = tmp_dn
          enddo
        enddo
        counter=1
        do ibatch = 1,n_batches
            do i=(ibatch-1)*read_batch_size+1,min(n_nonzero_elems_upper_triangular,ibatch*read_batch_size)
              read(57,*) minus_tau_H_indices(counter),minus_tau_H_values(counter)
              counter=counter+1
            enddo
        enddo
        close(57)
      endif ! master_core
      ! Send all elems to all processors
      minus_tau_H_values=-tau*minus_tau_H_values
      call mpi_bsend(minus_tau_H_nonzero_elements)
      call mpi_bsend(minus_tau_H_values)
      call mpi_bsend(minus_tau_H_indices)
      call mpi_bsend(imp_up)
      call mpi_bsend(imp_dn)


      write (6,*) "Finished reading",n_nonzero_elems_upper_triangular,"nonzero elements in deterministic Hamiltonian."
      call my_second(2,'read deterministic matrix elements')
    else
      call my_second(1,'generate important space')
      if (diff_from_psi_t.or.use_psit_con_in) then
        if (ndet_psi_t_in>0.or.use_psit_con_in) then
          call generate_space_iterate(imp_iters,norb_imp,n_imp_initiators,n_imp_truncate,n_det=n_imp,dets_up=imp_up,dets_dn=imp_dn,values=minus_tau_H_values,H_indices=minus_tau_H_indices,H_nonzero_elements=minus_tau_H_nonzero_elements,init_up=dets_up_psi_t_in,init_dn=dets_dn_psi_t_in,init_wts=cdet_psi_t_in,dtm_energy=dtm_energy)
        else
          call generate_space_iterate(imp_iters,norb_imp,n_imp_initiators,n_imp_truncate,n_det=n_imp,dets_up=imp_up,dets_dn=imp_dn,values=minus_tau_H_values,H_indices=minus_tau_H_indices,H_nonzero_elements=minus_tau_H_nonzero_elements,dtm_energy=dtm_energy)
        endif
      endif
      write (6,'(''Deterministic space Hamiltonian has'',i12,'' ='',es9.2'' nonzero values, n_imp='',i9,'' ='',es9.2)') size(minus_tau_H_values), real(size(minus_tau_H_values)), n_imp, real(n_imp)
      call my_second(2,'generate important space')

      if (use_elems_out) then
        if (master_core) then
          allocate(orb_up(nup))
          allocate(orb_dn(ndn))
          write (6,*) "Dumping the deterministic matrix elements into file ",dtm_elems_out_file
          call my_second(1,'dump deterministic matrix elements'); call flush(6)
          open(8, file=dtm_elems_out_file,status='new')
          write(8,*) n_imp, size(minus_tau_H_values), dtm_energy, "number of deterministic dets, number of nonzero deterministic Hamiltonian elements, ground state energy within deterministic space"
          write (fmt, '(i15)') n_imp
          write(8, '(' // trim(fmt) // 'i8)') minus_tau_H_nonzero_elements

          write (fmt, '(i5)') nup+ndn-2*n_core_orb-1
          do i = 1, n_imp
#ifdef NUM_ORBITALS_GT_127
            tmp_det = imp_up(i)-maskr_vec(n_core_orb)
#else
            tmp_det = imp_up(i)-maskr(n_core_orb,ik)
#endif

            do j=n_core_orb+1,nup
              orb_up(j) = trailz(tmp_det)+1
              tmp_det = ibclr(tmp_det,orb_up(j)-1)
            enddo
#ifdef NUM_ORBITALS_GT_127
            tmp_det = imp_dn(i)-maskr_vec(n_core_orb)
#else
            tmp_det = imp_dn(i)-maskr(n_core_orb,ik)
#endif

            do j=n_core_orb+1,ndn
              orb_dn(j) = trailz(tmp_det)+1
              tmp_det = ibclr(tmp_det,orb_dn(j)-1)
            enddo
            do j=n_core_orb+1,nup
              orb_up(j) = orb_up(j) - n_core_orb
            enddo
            do j=n_core_orb+1,ndn
              orb_dn(j) = orb_dn(j) - n_core_orb
            enddo
            write(8,*) i, orb_up(n_core_orb+1:nup), orb_dn(n_core_orb+1:ndn)
           !write(8,'(i5,' // trim(fmt) // 'i4)') i, orb_up(n_core_orb+1:nup), orb_dn(n_core_orb+1:ndn)
          enddo
          do i=1,size(minus_tau_H_values)
            write (8,*) minus_tau_H_indices(i),-minus_tau_H_values(i)/tau
           !write (8,'(i8,f15.21)') minus_tau_H_indices(i),-minus_tau_H_values(i)/tau
          enddo
          call my_second(2,'dump deterministic matrix elements'); call flush(6)
        endif
        close(8)
      endif ! use_elems_out

      if (print_psit_wo_sqmc) then
        write (6,*) "Printed deterministic matrix elements in ",dtm_elems_out_file
        write (6,*) "Terminating SQMC"
        call flush(6)
        stop "Printed deterministic matrix elements. Terminating SQMC."
      endif

    endif ! use_elems_in
    ! At this point, the deterministic projector is done being read/generated.

    ! Edited by AAH on 13 May 2014. Check that all connections to HF are in deterministic space Hamiltonian
   !if (hamiltonian_type.ne.'hubbard2') then
   !  call my_second(1,"Checking that all connections to HF are in deterministic space"); call flush(6)
   !  write (6,*) "HF=",imp_up(1),imp_dn(1); call flush(6)

   !  if(hamiltonian_type.eq.'heg') then
   !    call find_connected_dets_heg(imp_up(1), imp_dn(1), n_connected_dets, connected_dets_up, connected_dets_dn)
   !  elseif(hamiltonian_type.eq.'chem') then
   !    call find_connected_dets_chem(imp_up(1), imp_dn(1), n_connected_dets, connected_dets_up, connected_dets_dn, norb)
   !  elseif(hamiltonian_type.eq.'hubbardk') then
   !    call find_connected_dets_hubbard_k(imp_up(1), imp_dn(1), n_connected_dets, connected_dets_up, connected_dets_dn, nsites)
   !  else
   !    stop 'Error: hamiltonian_type must be one of heg, chem, hubbard2, hubbardk'
   !  endif

   !  do i=1,n_connected_dets
   !    call binary_search(connected_dets_up(i),connected_dets_dn(i),imp_up,imp_dn,k)
   !    if (k==0) then
   !      write (6,*) "This connection to HF not found in deterministic space:",connected_dets_up(i),connected_dets_dn(i); call flush(6)
   !      stop "Connections to HF not included in deterministic projector!"
   !    endif
   !  enddo

   !  call my_second(2,"Checking that all connections to HF are in deterministic space"); call flush(6)
   !else
   !  write (6,*) "Real space hubbard used, so skip checking whether connections to HF lie in deterministic space"; call flush(6)
   !endif

    ! Moved here by AAH on 23 Jan 2013 so that MWALK can be set using n_imp and to ndet_psi_t_connected
    ! If MWALK is clearly too small, reset it.  If it is smaller than the default, then just print warning.  If it is too large, stop.
    if(MWALK.lt.0) stop 'MWALK < 0'
    if (hf_to_psit .and. MWALK.lt.w_abs_gen_target+n_imp+ndet_psi_t_connected) then
      !MWALK=int(3*(w_abs_gen_target/min_wt+ndet_psi_t_connected))/ncores
      MWALK=max(MWALK,int(3*(w_abs_gen_target/min_wt+ndet_psi_t_connected))/ncores)
      write (6,*) "w_abs_gen_target,min_wt,ndet_psi_t_connected=",w_abs_gen_target,min_wt,ndet_psi_t_connected; call flush(6)
      write(6,'(''Setting MWALK=3*(w_abs_gen_target/min_wt+ndet_psi_t_connected)/ncores='',i10,'' ='',es9.2)') MWALK, float(MWALK)
      !write(6,'(''Setting MWALK=3*(w_abs_gen_target/min_wt+ndet_psi_t_connected)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
    elseif (.not.hf_to_psit .and. MWALK.lt.w_abs_gen_target+n_imp) then
      !MWALK=int(4*(w_abs_gen_target/min_wt+n_imp))/ncores
      MWALK=max(MWALK,int(4*(w_abs_gen_target/min_wt+n_imp))/ncores)
      write(6,'(''5Setting MWALK=4*(w_abs_gen_target/min_wt+n_imp)/ncores='',i10,'' ='',es9.2)') MWALK, float(MWALK)
      !write(6,'(''5Setting MWALK=4*(w_abs_gen_target/min_wt+n_imp)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
    endif

    if (hf_to_psit) then
      allocate(empty_outside_ct(max(w_abs_gen_target,ndet_psi_t_connected)))
      allocate(indices_new_dets(max(w_abs_gen_target,ndet_psi_t_connected)))
      allocate(leave_out_list(max(2*w_abs_gen_target,ndet_psi_t_connected)))
    endif

    if ((hf_to_psit .and. MWALK.lt.int(2.5*(w_abs_gen_target+n_imp+ndet_psi_t_connected))) .or. (.not.hf_to_psit .and. MWALK.lt.3*(w_abs_gen_target+n_imp))) then
        write (6,*) "MWALK=",MWALK
        write(6,'(''Warning: MWALK may be too small.  Set it to zero to get default MWALK'')')
        if (hf_to_psit) then
          MWALK=max(MWALK,nint(2.5*(w_abs_gen_target+n_imp+ndet_psi_t_connected),i8b)/ncores)
          write (6,*) "Setting MWALK=2.5*(w_abs_gen_target+n_imp+ndet_psi_t_connected)/ncores=",MWALK
        endif
        call flush(6)
    endif
    if(MWALK.gt.2147483647 .or. MWALK.lt.0) then
      write(6,'(''At present MWALK cannot be > 2^31 because iorder is 4-byte'')')
      stop 'At present MWALK cannot be > 2^31 because iorder is 4-byte'
    endif
    ! End of edit

    ! Now store diagonal elements for transformed projector (those in connections to psi trial but outside of deterministic space)
    if (hf_to_psit) then
      allocate(diag_elems(ndet_psi_t_connected))
      diag_elems(:) = 0._rk
      do i=1,ndet_psi_t_connected
        call binary_search(psi_t_connected_dets_up(i),psi_t_connected_dets_dn(i),imp_up(1:n_imp),imp_dn(1:n_imp),j)
        if (j==0) then
          if(hamiltonian_type.eq.'heg') then
            call hamiltonian_heg(psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), matrix_element)
          elseif(hamiltonian_type.eq.'chem') then
            if (time_sym) then
              call hamiltonian_chem_time_sym(psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), matrix_element)
            else
              call hamiltonian_chem(psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), 0, matrix_element)
            endif
          elseif(hamiltonian_type.eq.'hubbard2') then
            call hamiltonian_hubbard(psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), matrix_element)
          elseif(hamiltonian_type.eq.'hubbardk') then
            if (space_sym) then
              call hamiltonian_hubbard_k_space_sym(psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), matrix_element,nnzero)
            else
              call hamiltonian_hubbard_k(psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), psi_t_connected_dets_up(i), psi_t_connected_dets_dn(i), matrix_element)
            endif
          endif
          diag_elems(i) = matrix_element
        endif
      enddo
    endif
  else ! not semistochastic
    !MWALK=int(4*(w_abs_gen_target/min_wt))/ncores
    MWALK=max(MWALK,int(4*(w_abs_gen_target/min_wt))/ncores)
    write(6,'(''6Setting MWALK=4*(w_abs_gen_target/min_wt)/ncores='',i10,'' ='',es9.2)') MWALK, float(MWALK)
    !write(6,'(''6Setting MWALK=4*(w_abs_gen_target/min_wt)='',i10,'' ='',es9.2)') MWALK, float(MWALK)
  endif
  ! End of edit
  !write(6,'(/,''minus_tau_H number of nonzero elements per row'',100i9,/)') minus_tau_H_nonzero_elements

  allocate(walk_wt(MWALK),stat=istat)
  if(istat.ne.0) stop 'failed to allocate walk_wt'
  if (run_type .eq. 'sr') then ! If stochastic reconfiguration being used, store the fixed node weights separately
    allocate(walk_wt_fn(MWALK),stat=istat)
    if(istat.ne.0) stop 'failed to allocate walk_wt_fn'
  endif
  allocate(matrix_elements(MWALK),stat=istat)
  if(istat.ne.0) stop 'failed to allocate matrix_elements'
  matrix_elements(:) = 1.e51_rk
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    allocate(e_num_walker(MWALK),stat=istat)
    if(istat.ne.0) stop 'failed to allocate e_num_walker'
    e_num_walker(:)=1.e51_rk
    allocate(e_den_walker(MWALK),stat=istat)
    if(istat.ne.0) stop 'failed to allocate e_den_walker'
    e_den_walker(:)=1.e51_rk
  endif
  allocate(initiator(MWALK),stat=istat)
  if(istat.ne.0) stop 'failed to allocate initiator'
  allocate(iorder(MWALK),stat=istat)
  if(istat.ne.0) stop 'failed to allocate iorder'
  allocate(temp_i_2((MWALK+1)/2),stat=istat)
  if(istat.ne.0) stop 'failed to allocate temp_i_2'
  allocate(sign_permanent_initiator(M_PERMANENT_INITIATOR),stat=istat)
  if(istat.ne.0) stop 'failed to allocate sign_permanent_initiator'

  allocate(imp_distance(MWALK),stat=istat)
  if(istat.ne.0) stop 'failed to allocate imp_distance'
  imp_distance(:) = 1_i1b

  if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
    allocate(walk_det(MWALK),stat=istat)
    if(istat.ne.0) stop 'failed to allocate walk_det'
    allocate(temp_i_1((MWALK+1)/2),stat=istat)        ! for merge_sort, walk_det
    if(istat.ne.0) stop 'failed to allocate temp_i_1'

    call shell_sort_realwt(iwdet_psi_t,cdet_psi_t,ndet_psi_t)

    cdet_psi_t_abs_max=maxval(abs(cdet_psi_t(1:ndet_psi_t)))
    cdet_psi_t_abs_sum=sum(abs(cdet_psi_t(1:ndet_psi_t)))
    iwalk=0
    do idet=1,ndet_psi_t
      iwalk=iwalk+1
      walk_det(iwalk)=iwdet_psi_t(idet)
      if(hf_to_psit) then
        if(iwalk.eq.1) then
          walk_wt(iwalk)=w_abs_gen_begin
        else
          walk_wt(iwalk)=0
        endif
      else
        walk_wt(iwalk)=w_abs_gen_begin*cdet_psi_t(idet)/cdet_psi_t_abs_sum ! Note if we use importance sampling, it would be more consistent to use the squares in both the numerator and denominator
      endif
      if(abs(abs(cdet_psi_t(idet))-cdet_psi_t_abs_max).lt.1.d-3) then ! Only those that are nearly equal and of the same sign are deemed permanent initiators
        n_permanent_initiator=n_permanent_initiator+1
        if(n_permanent_initiator > M_PERMANENT_INITIATOR) then
          write(6,'(''n_permanent_initiator > M_PERMANENT_INITIATOR'')')
          stop 'n_permanent_initiator > M_PERMANENT_INITIATOR'
        endif
        !sign_permanent_initiator(n_permanent_initiator)=int(sign(1._rk,cdet_psi_t(idet)))            ! Correction - sign is 1 if importance sampling is used!
        if (importance_sampling .eq. 0) then
          sign_permanent_initiator(n_permanent_initiator)=int(sign(1._rk,cdet_psi_t(idet)),i1b)
        else
          sign_permanent_initiator(n_permanent_initiator)=1                                           ! Correction - sign is 1 if importance sampling is used!
        endif
        initiator(iwalk)=3
        if(semistochastic) imp_distance(iwalk)=0
!       if(.not.semistochastic) n_imp=n_imp+1
      else
        initiator(iwalk)=2
      endif
    enddo

    if(iwalk.eq.0) then
      iwalk=1
      idet=maxloc(abs(cdet_psi_t(:)),1)
      walk_det(iwalk)=iwdet_psi_t(idet)
      walk_wt(iwalk)=1
      initiator(iwalk)=3
      imp_distance(iwalk)=0
      write(6,'(''It should never get here'')')
      stop 'It should never get here'
    endif

    nwalk=iwalk
    if(run_type.eq.'fixed_node1' .and. importance_sampling.eq.1)  walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(run_type.eq.'fixed_node2' .and. importance_sampling.eq.1)  walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(run_type.eq.'fixed_node3' .and. importance_sampling.eq.1)  walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(run_type.eq.'sr'           .and. importance_sampling.eq.1) walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(run_type.eq.'vmc') walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(importance_sampling.eq.1) walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))

    write(6,'(''Initial nwalk, w_gen, w_abs_gen='',i8,9f8.1)') nwalk, sum(walk_wt(1:nwalk)), sum(abs(walk_wt(1:nwalk)))
    write(6,'(''Initial nwalk,walk_det ='',i5,'': '',1000i6)') nwalk,(walk_det(iwalk),iwalk=1,nwalk)
    write(6,'(''Initial nwalk,walk_wt  ='',i5,'': '',1000f6.1)') nwalk,(walk_wt(iwalk),iwalk=1,nwalk)
    write(6,'(''Initial nwalk,initiator='',i5,'': '',1000i6)') nwalk,(initiator(iwalk),iwalk=1,nwalk)
    write(6,'(''sign_permanent_initiator ='',100i3)') sign_permanent_initiator(1:n_permanent_initiator)

  elseif(hamiltonian_type.eq.'heg' .or. hamiltonian_type.eq.'chem' .or. hamiltonian_type .eq. 'hubbard2' .or. hamiltonian_type .eq. 'hubbardk' .or. hamiltonian_type .eq. 'hubbarddm') then

    allocate(walk_dets_up(MWALK),stat=istat)
    if(istat.ne.0) stop 'failed to allocate det_up'
    walk_dets_up(:) = 0_ik
    allocate(walk_dets_dn(MWALK),stat=istat)
    if(istat.ne.0) stop 'failed to allocate det_dn'
    walk_dets_dn(:) = 0_ik
!   allocate(imp_distance(MWALK),stat=istat)
!   if(istat.ne.0) stop 'failed to allocate imp_distance'
!   imp_distance(:) = 1_i1b
!   allocate(temp_i16((MWALK+1)/2),stat=istat)        ! for merge_sort2, walk_dets_up and walk_dets_dn
!   if(istat.ne.0) stop 'failed to allocate temp_i16'
    allocate(temp_i16_up((MWALK+1)/2),stat=istat)        ! for merge_sort2_up_dn, walk_dets_up and walk_dets_dn
    if(istat.ne.0) stop 'failed to allocate temp_i16_up'
    allocate(temp_i16_dn((MWALK+1)/2),stat=istat)        ! for merge_sort2_up_dn, walk_dets_up and walk_dets_dn
    if(istat.ne.0) stop 'failed to allocate temp_i16_dn'
! Warning: in above lines we are temporarily allocating both temp_i16 and temp_i16_up, temp_i16_dn, but we really need only one or the other.

! Populate all deterministic dets at start of run.
    iwalk=0
    do idet=1,n_imp
      iwalk=iwalk+1
      walk_dets_up(iwalk) = imp_up(idet)
      walk_dets_dn(iwalk) = imp_dn(idet)
      walk_wt(iwalk) = 0._rk
      initiator(iwalk) = 2
      imp_distance(iwalk)=0
    enddo

! Populate all dets in psi_T at start of run.  There are common dets in the deterministic space and in psi_T, but duplicates will be merged later.

    ! Sort psi trial by label, since it is currently sorted by abs wt
    call shell_sort2_realwt(dets_up_psi_t,dets_dn_psi_t,cdet_psi_t,ndet_psi_t)

    n_permanent_initiator=0

    if (.true.) then !.not.use_psit_con_in) then
      cdet_psi_t_abs_max=maxval(abs(cdet_psi_t(1:ndet_psi_t)))
      cdet_psi_t_abs_sum=sum(abs(cdet_psi_t(1:ndet_psi_t)))

      ! Edited by AAH on 8 Jan 2013
      if (hf_to_psit) then

        i_psit = 0
        do idet=1,ndet_psi_t_connected
          iwalk=iwalk+1
          walk_dets_up(iwalk) = psi_t_connected_dets_up(idet)
          walk_dets_dn(iwalk) = psi_t_connected_dets_dn(idet)
          if (idet==1) then
            i_psit=i_psit+1
            walk_wt(iwalk)=w_abs_gen_begin ! Note if we use importance sampling, it would be more consistent to use the squares in both the numerator and denominator
            if(abs(abs(cdet_psi_t(i_psit))-cdet_psi_t_abs_max).lt.1.d-3) then ! Only those that are nearly equal and of the same sign are deemed permanent initiators
              n_permanent_initiator=n_permanent_initiator+1
              if(n_permanent_initiator > M_PERMANENT_INITIATOR) then
                write(6,'(''n_permanent_initiator > M_PERMANENT_INITIATOR'')')
                stop 'n_permanent_initiator > M_PERMANENT_INITIATOR'
              endif
              !sign_permanent_initiator(n_permanent_initiator)=int(sign(1._rk,cdet_psi_t(idet)))            ! Correction - HJC - Nov 4
              if (importance_sampling .eq. 0) then
                sign_permanent_initiator(n_permanent_initiator)=int(sign(1._rk,cdet_psi_t(i_psit)),i1b)
              else
                sign_permanent_initiator(n_permanent_initiator)=1                                           ! Correction - sign is 1 if importance sampling is used!
              endif
              initiator(iwalk)=3
              imp_distance(iwalk)=-2
            else
              initiator(iwalk)=2 ! It makes no difference that small-wt walkers are being labeled as initiators, since that will be fixed in merge_original_with_spawned2
              imp_distance(iwalk)=-2
            endif
          else
            walk_wt(iwalk) = 0._rk
            initiator(iwalk) = 2
            imp_distance(iwalk)=-2
          endif
        enddo
      else

        do idet=1,ndet_psi_t
          iwalk=iwalk+1
          walk_dets_up(iwalk)=dets_up_psi_t(idet)
          walk_dets_dn(iwalk)=dets_dn_psi_t(idet)
          walk_wt(iwalk)=w_abs_gen_begin*cdet_psi_t(idet)/cdet_psi_t_abs_sum ! Note if we use importance sampling, it would be more consistent to use the squares in both the numerator and denominator
          if(abs(abs(cdet_psi_t(idet))-cdet_psi_t_abs_max).lt.1.d-3) then ! Only those that are nearly equal and of the same sign are deemed permanent initiators
            n_permanent_initiator=n_permanent_initiator+1
            if(n_permanent_initiator > M_PERMANENT_INITIATOR) then
              write(6,'(''n_permanent_initiator > M_PERMANENT_INITIATOR'')')
              stop 'n_permanent_initiator > M_PERMANENT_INITIATOR'
            endif
            !sign_permanent_initiator(n_permanent_initiator)=int(sign(1._rk,cdet_psi_t(idet)))            ! Correction - HJC - Nov 4
            if (importance_sampling .eq. 0) then
              sign_permanent_initiator(n_permanent_initiator)=int(sign(1._rk,cdet_psi_t(idet)),i1b)
            else
              sign_permanent_initiator(n_permanent_initiator)=1                                           ! Correction - sign is 1 if importance sampling is used!
            endif
            initiator(iwalk)=3
            imp_distance(iwalk)=0
           if(.not.semistochastic) n_imp=n_imp+1
          else
            initiator(iwalk)=2 ! It makes no difference that small-wt walkers are being labeled as initiators, since that will be fixed in merge_original_with_spawned2
          endif
        enddo
      endif
      ! End of edit
    endif ! hf_to_psit

    nwalk=iwalk

! If the largest absolute walk_wt is less than 1, make it 1.
    do iwalk=1,nwalk
      walk_wt(iwalk)=walk_wt(iwalk)/min(w_abs_gen_begin*cdet_psi_t_abs_max/cdet_psi_t_abs_sum,1._rk)
    enddo
    if((ipr.ge.1).and.(master_core)) write(6,'(''walk_wt='',1000f10.6)') walk_wt(1:n_imp+ndet_psi_t)

    call my_second(1,'populate initial dets')
! Merge the walkers populated from important dets and those populated from psi_T.  (There are usually some in common.)
!   n_imp=count(imp_distance(1:nwalk).eq.0)
    do i=1,nwalk
      iorder(i)=i
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Sort the determinants and assign the objects associated with them the right order
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call merge_sort2_up_dn(walk_dets_up,walk_dets_dn, iorder, nwalk, temp_i16_up, temp_i16_dn, temp_i_2)
    walk_wt(1:nwalk)=walk_wt(iorder(1:nwalk))

    if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
      e_num_walker(1:nwalk)=e_num_walker(iorder(1:nwalk))
      e_den_walker(1:nwalk)=e_den_walker(iorder(1:nwalk))
    endif
    initiator(1:nwalk)=initiator(iorder(1:nwalk))
    imp_distance(1:nwalk)=imp_distance(iorder(1:nwalk))
    matrix_elements(1:nwalk)=matrix_elements(iorder(1:nwalk))


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Take the sorted determinants and merge them - that is, perform cancellations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call merge_original_with_spawned2(0,0,walk_dets_up,walk_dets_dn,walk_wt,nwalk)!,n_permanent_initiator)
    ! Edited by AAH on 12 Mar 2013. This is necessary because otherwise, the merge routine makes low-wt dets non-initiators
    if (hf_to_psit) then
      do i=1,nwalk
        if (initiator(i).ne.3)  initiator(i)=2
      enddo
    endif
    ! End of edit
    call my_second(2,'populate initial dets')

    n_imp=count(imp_distance(1:nwalk).eq.0)
      write(6,'(/,''At start of run: nwalk, n_imp, w_abs_gen='',2i9,f10.1)') nwalk, n_imp, sum(abs(walk_wt(1:nwalk)))
      call flush(6)

    if(nwalk.eq.0) then
      nwalk=1
      idet=maxloc(abs(cdet_psi_t(:)),1)
      walk_dets_up(nwalk)=dets_up_psi_t(idet)
      walk_dets_dn(nwalk)=dets_dn_psi_t(idet)
      walk_wt(nwalk)=1
      initiator(nwalk)=3
      imp_distance(nwalk)=0
      write(6,'(''It should never get here'')')
      stop 'It should never get here'
    endif

    if(run_type.eq.'fixed_node1' .and. importance_sampling.eq.1) walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(run_type.eq.'fixed_node2' .and. importance_sampling.eq.1) walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(run_type.eq.'fixed_node3' .and. importance_sampling.eq.1) walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(run_type.eq.'sr'          .and. importance_sampling.eq.1) walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(importance_sampling.eq.1) walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    if(run_type.eq.'vmc') then
      walk_wt(1:nwalk)=abs(walk_wt(1:nwalk))
    endif

    write(6,'(/,''Initial nwalk, w_gen, w_abs_gen='',i8,9f9.1)') nwalk, sum(walk_wt(1:nwalk)), sum(abs(walk_wt(1:nwalk)))
    if(ipr.ge.1 .or. (ipr.ge.0 .and. nwalk.le.1000)) then
        write(6,'(''Initial walk_dets_up/dn='',1000(2i7,2x))') nwalk,(walk_dets_up(iwalk),walk_dets_dn(iwalk),iwalk=1,nwalk)
        write(6,'(''Initial walk_wt        ='',1000f16.1)') (walk_wt(iwalk),iwalk=1,nwalk)
        write(6,'(''Initial initiator      ='',1000i16)') (initiator(iwalk),iwalk=1,nwalk)
        write(6,'(''sign_permanent_initiator ='',1000i16)') sign_permanent_initiator(1:n_permanent_initiator)
    endif
    call flush(6)
  else
    write(6,'(''hamiltonian_type must be one of: fictitious read hubbard heg  chem'')')
    stop 'hamiltonian_type must be one of: fictitious read hubbard heg  chem'
  endif

  counter_two_body = 0

  if(reweight_factor_inv_max.eq.0._rk) reweight_factor_inv_max=1+reweight_factor_inv_max_multiplier*tau
  write(6,'(/,''reweight_factor_inv_max='',f12.8)') reweight_factor_inv_max
  log_n_con = log(real(ndet_psi_t_connected))/log(2.0)

  ! Create histograms of the weight generated by off-diagonal moves for all dets and for HF det only
  if (off_diag_histograms) then
    nbins = 10001
    call gen_hist(lo=0.0_rk,hi=10000._rk,nbins=nbins,lbounds=lbounds,bins=bins)
    call gen_hist(lo=0.0_rk,hi=10000._rk,nbins=nbins,lbounds=lbounds,bins=bins_single)
    call gen_hist(lo=0.0_rk,hi=10000._rk,nbins=nbins,lbounds=lbounds,bins=bins_double)
    call gen_hist(lo=0.0_rk,hi=10000._rk,nbins=nbins,lbounds=lbounds_hf,bins=bins_hf)
    ! Also create histograms of the weight after the move is performed, (i.e., including the weight of the walker on the initial det)
    call gen_hist(lo=0.0_rk,hi=10000._rk,nbins=nbins,lbounds=lbounds_new,bins=bins_new)
  endif

  call walk

  end subroutine read_input
!-------------------------------------------------------------------------------------

! ==============================================================================
! subroutine find_abs_largest_elems(ndet_psi_t, cdet_psi_t, n_abs_largest_elems)

! integer idet
! real(rk) eps=1.e-6_rk
! real(rk) cdet_abs_largest

! cdet_abs_largest=0
! do idet=1,ndet_psi_t
!   cdet_abs_largest=max(cdet_abs_largest,abs(cdet_psi_t(idet)))
! enddo

! iwalk=0
! do idet=1,ndet_psi_t
!   if(abs(cdet_abs_largest-abs(cdet_psi_t(idet))) .lt. eps) then
!     iwalk=iwalk+1
!     walk_det(iwalk)=iwdet_psi_t(idet)
!     walk_wt(iwalk)=nint(sign(1._rk,cdet_psi_t(idet)))
!     initiator(iwalk)=3
!   endif
! enddo

! end subroutine find_abs_largest_elems
!-------------------------------------------------------------------------------------

! ==============================================================================
  subroutine walk
! ------------------------------------------------------------------------------
! Description   : Do the walk
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use types, only: ik, rk, i8b
  use hamiltonian_mod, only: energy_exact, wavefn_exact, ham ! variables
  use hamiltonian_mod, only: shell_sort_eigvec               ! routines
  use heg, only: energy_pieces_heg
  use chemistry, only: energy_pieces_chem
  use hubbard, only: energy_pieces_hubbard,energy_pieces_hubbard_k,energy_pieces_hubbard_dm
  use read_psi_trial
  use common_ham, only: hamiltonian_type, diagonalize_ham, diagonal_ham_lowest, diagonal_ham_highest
! use common_psi_t, only: ndet_psi_t, dets_up_psi_t, dets_dn_psi_t, cdet_psi_t, n_excite_psi_t, trial_wf_iters, norb_trial_wf
  use common_psi_t, only: ndet_psi_t, dets_up_psi_t, dets_dn_psi_t, cdet_psi_t, trial_wf_iters, norb_trial_wf, hf_to_psit, ndet_psi_t_connected, use_psit_con_in
  use common_run, only: tau, reweight_factor_inv_max_multiplier, reweight_factor_inv_max, initiator_rescale_power
  use hubbard, only: wf_type
  use semistoch
  use tools, only : do_once
  use more_tools, only : sparse_matrix_multiply, fast_sparse_matrix_multiply, fast_sparse_matrix_multiply_upper_triangular, linear_search_list_and_update, sparse_fast_sparse_matrix_multiply_upper_triangular,reduce_matrix_indices,binary_search_list_and_update,fast_sparse_matrix_multiply_local_band, from_upper_triangular_to_band,fast_sparse_matrix_multiply_local_band_block,from_upper_triangular_to_band_shuffle,from_upper_triangular_to_band_shuffle_old
  use common_imp, only: norb_imp
  use mpi_routines, only : mpi_barr,whoami,mpi_bsend
! use common_walk, only: m_connected_dets ! presently used at top of module
  implicit none

! local
  character*22 write_string

  integer, parameter :: M_DOMINANT_DET=10000, max_imp_distance=15
  integer n_sign_flips
! integer(i8b) n_dominant_det
! integer(i8b) n_dominant_det, dominant_det_most, dominant_det, dominant_det_list(M_DOMINANT_DET)
  integer n_dominant_det
  integer(i8b) dominant_det_most, dominant_det, dominant_det_list(M_DOMINANT_DET)
  real(rk) dominant_det_wt_tmp, dominant_det_abs_wt_tmp, dominant_det_wt(M_DOMINANT_DET), dominant_abs_det_wt(M_DOMINANT_DET)

! integer(i8b) w_perm_initiator_gen, w_abs_gen, w_gen, w2_gen, w_blk, w2_blk
  ! Edited by AAH on 31 Jan 2013. Make iblk, istep global variables, so that we don't have to keep passing them into subroutines.
  integer(i8b) :: nwalk_before_merge
  integer iwalk,istat,idet,idet_begin!, changed_sign
 !integer iwalk,istat,nwalk_before_merge,idet,idet_begin!, changed_sign
 !integer iwalk,nwalk_old,istat,nwalk_before_merge,idet,idet_begin!, changed_sign
 !integer iwalk,nwalk_old,istep,iblk,iblkk,istat,nwalk_before_merge,idet,idet_begin!, changed_sign
  ! End of edit
  integer :: which_perm_initiator
  integer ntimes_nblk_eq, i_equilibration_run, equilibrated_w_abs_blk, equilibrated_e_blk, equilibrated_w_perm_initiator_blk, equilibrated_all
  integer, allocatable :: index(:)
  real(rk) :: tau_sav, tau_previous, tau_ratio, r_initiator_sav
  real(rk) :: proposal_prob_inv, norm_of_wavefn, overlap!, e_mix_numerator, e_mix_denominator!, cdet_psi_t_this_walk

! integer(i8b) nwalk_cum, nwalk_before_merge_cum, w_perm_initiator_cum, w_cum, w2_cum, w_gen2_cum, w_blk2_cum, w_abs_cum, w_genabs_cum, w_blkabs_cum, w_abs_before_merge_cum
  integer(i8b) nwalk_cum, nwalk_before_merge_cum
  real(rk) w_perm_initiator_gen, w_abs_gen_prev, w_abs_gen, w_abs_gen_imp, w_gen, w2_gen, w_blk, w2_blk, w_abs_blk, w_abs_blk_imp, w_abs_blk_prev, w_abs_genlast_cum
  real(rk) w_perm_initiator_blk, w_perm_initiator_blk_prev
  real(rk) w_perm_initiator_cum, w_cum, w2_cum, w_gen2_cum, w_blk2_cum, w_abs_cum, w_genabs_cum, w_blkabs_cum, w_abs_before_merge_cum
  real(rk) &
 &n_eff, n_abs_eff, n_genabs_eff, n_blkabs_eff, &
 &e_den_eff, e_den_abs_eff, e_den_genabs_eff, e_den_blkabs_eff, &
!&w_gen, w2_gen, w_blk, w2_blk, &
 &e_num, e_num_gen, e_num_blk, &
 &e_den, e_den_gen, e_den_blk, &
 &e_num_e_den, e_num_gen_e_den_gen, e_num_blk_e_den_blk, &
 &e_gen, e_genp, e_gen_gen1, &
 &e_num_cum, e_num2_cum, e_num_gen2_cum, e_num_blk2_cum, e_num_abs_cum, e_num_genabs_cum, e_num_blkabs_cum, &
 &e_den_cum, e_den2_cum, e_den_gen2_cum, e_den_blk2_cum, e_den_abs_cum, e_den_genabs_cum, e_den_blkabs_cum, &
 &e_gen_cum, e_gen2_cum, e_genp_e_genpp_cum, e_genp_cum, e_genpp_cum, &
 &e_num_e_den_cum, e_num_gen_e_den_gen_cum, e_num_blk_e_den_blk_cum, &
 &e_num_av, e_num_abs_av, e_num_genabs_av, e_num_blkabs_av, &
 &e_den_av, e_den_abs_av, e_den_genabs_av, e_den_blkabs_av, &
 &e_num_e_den_av, e_num_gen_e_den_gen_av, e_num_blk_e_den_blk_av, &
 &e_gen_av, e_gen2_av, e_genp_e_genpp_av, e_genp_av, e_genpp_av, &
 &e_num_sigma, e_num_abs_sigma, e_num_genabs_sigma, e_num_blkabs_sigma, &
 &e_den_sigma, e_den_abs_sigma, e_den_genabs_sigma, e_den_blkabs_sigma, &
 &e_num_e_den_var, e_num_gen_e_den_gen_var, e_num_blk_e_den_blk_var, &
 &e_num_gen_ave, e_num_gen_vn1, e_num_blk_ave, e_num_blk_vn1, &
 &e_den_gen_ave, e_den_gen_vn1, e_den_blk_ave, e_den_blk_vn1, &
 &e_num_genabs, e_num_blkabs, &
 &e_den_genabs, e_den_blkabs, &
 &e_num_genabs_del, e_num_genabs_ave, e_num_genabs_vn1, e_num_genabs_var, &
 &e_num_blkabs_del, e_num_blkabs_ave, e_num_blkabs_vn1, e_num_blkabs_var, &
 &e_den_genabs_del, e_den_genabs_ave, e_den_genabs_vn1, e_den_genabs_var, &
 &e_den_blkabs_del, e_den_blkabs_ave, e_den_blkabs_vn1, e_den_blkabs_var, &
 &e_genabs_ave, e_blkabs_ave, &
 &e_num_genabs_e_den_genabs_vn1, e_num_genabs_e_den_genabs_var, e_num_genabs_e_den_genabs_sigma, &
 &e_num_blkabs_e_den_blkabs_vn1, e_num_blkabs_e_den_blkabs_var, e_num_blkabs_e_den_blkabs_sigma, &
 &e_gen_del, e_gen_ave, e_gen_vn1, e_gen_var, &
 &e_genp_del, &
 &e_genpp_del, e_genpp_ave, &
 &e_genp_e_genpp_vn1, e_genp_e_genpp_var, &
 &e_blk, e_blk_prev, &
 &e_av, e_abs_av, e_genabs_av, e_blkabs_av, &
 &e_err, e_abs_err, e_genabs_err, e_blkabs_err, e_blkabs_corrected_err, &
 &tau_corr_nonint, t_corr_nonint, t_corr, overlap_with_psi_t, overlap_with_psi_t_sav, passes, &
 &imp_distance_fraction(-2:max_imp_distance), &
 &p_eigval

! Tmp ovlp only for hubbard 2
  real(rk) ::tmp_ovlp

! integer, parameter :: M_DOMINANT_DET=10000
! integer(i8b) n_dominant_det, dominant_det_most, dominant_det, dominant_det_list(M_DOMINANT_DET)
! real(rk) dominant_det_wt_tmp, dominant_det_abs_wt_tmp, dominant_det_wt(M_DOMINANT_DET), dominant_abs_det_wt(M_DOMINANT_DET)

! integer(i8b), allocatable :: wavefn(:)
  real(rk), allocatable :: wavefn(:)

  integer :: i,j,k,n
  real(rk) :: imp_wt(n_imp)!,delta_wt(n_imp) ! for deterministic projection
!  real(rk) :: delta_wt2(ndet_psi_t_connected) ! for HF -> Psi_T
!  real(rk) :: delta_wt3(ndet_psi_t) ! for HF -> Psi_T

  integer :: i_imp ! for energy pieces calculation
  integer(i8b) :: locations_of_imp_dets(n_imp)
  integer :: locations_of_psit(ndet_psi_t)
  integer :: i_permanent_initiator ! not really needed
  logical :: discard_walker
  logical,parameter :: use_improved_cancellation = .false.

!*** Added by AR
  real(rk) :: my_imp_wt(n_imp), my_w_abs_gen
  real(rk),allocatable :: deltaw(:)
  integer :: my_locations_of_imp_dets(n_imp)
  integer :: my_locations_of_psit(ndet_psi_t)

  logical :: iown_psit1, iown_first
  integer :: psit1_owner, first_owner, own_first_arr(1)

  integer :: my_nimp,my_ndet_psit,my_ndet_psi_t_connected,my_ndet_outside_ct
  integer(i8b), allocatable :: imp_mask(:)
  integer, allocatable :: ndet_psit_mask(:)

  integer :: my_nwlk_1
  integer :: my_new_walks,my_n_perm_init

  integer :: my_nwalk
  integer :: my_skip
 !integer :: my_nwalk,my_skip
  integer  :: my_cnt,icore
  integer(i8b) :: imp_stride


  integer(i8b) :: imp_core_mask(ncores)

  real(rk) ::tmp_real,wt1
  real(rk), allocatable ::  tmp_real_array(:)
  integer(i8b) :: disp_nwalk ! Just to display the total number of occupied determinants from all cores

  !*** needed for statistics
  real(rk) :: my_e_den2_cum,my_e_den_abs_cum,my_e_num2_cum,my_e_num_abs_cum,my_e_num_e_den_cum

  real(rk) :: my_w2_cum,my_w_abs_before_merge_cum,my_w_abs_blk_imp,my_w_perm_initiator_blk,my_w_perm_initiator_cum

  integer(i8b) :: my_nwalk_before_merge_cum,my_nwalk_cum
  integer :: hf_core

  logical :: growth_estim=.false. ! decide whether mixed (uses Psi_T) or growth estimator (uses Psi_G) is used to calculate e_trial in projector

! Warning: tmp
! open(11,file='wts_single')
! open(12,file='wts_double')
! open(13,file='wts_all')
! open(14,file='wts_large')

  if(growth_estim) then
    write(6,'(/,''Using growth estimator for population control'')')
  else
    write(6,'(/,''Using mixed estimator for population control'')')
  endif
  call flush(6)

  allocate(all_max_wt_ratio_excit_level(ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate max_wt_ratio_excit_level'

  if (hamiltonian_type.ne.'fictitious' .and. hamiltonian_type.ne.'read' .and. hamiltonian_type.ne.'hubbard'  .and. hamiltonian_type.ne.'hubbard2') hf_core = get_det_owner(walk_dets_up(1),walk_dets_dn(1))

  my_ndet_outside_ct=0
  my_e_den2_cum=0
  my_e_den_abs_cum=0
  my_e_num2_cum=0
  my_e_num_abs_cum=0
  my_e_num_e_den_cum=0

  allocate(imp_mask(n_imp),stat=istat)
  if(istat.ne.0) stop 'failed to allocate imp_mask'
  allocate(ndet_psit_mask(ndet_psi_t),stat=istat)
  if(istat.ne.0) stop 'failed to allocate ndet_psit_mask'

  allocate(tmp_real_array(15),stat=istat)
  if(istat.ne.0) stop 'failed to allocate tmp_real_array'

!  !***initialize send buffers in mod_cluster.f90
!  if(ncores>1) then
!    call init_snd_buffers(int(0.75*MWALK,8)) !we need less space for spowning dets
!  else
!    call init_snd_buffers(MWALK)
!  endif
!  !***initialize owners hash table in mod_cluster.f90
!  call init_hash_owners
!!*** end Added by AR

  ! Added by AR 16 Dec 2013
  ! reduces the dimension of diag_elems states list neglecting states not needed locally
  if(hf_to_psit .and. (ncores>1)) call reduce_connections_diag(psi_t_connected_dets_up,psi_t_connected_dets_dn,diag_elems)


!*** end Added by AR

  n_warning=0

  if(ipr.ge.0) then

!***b
    if(master_core) then
        open(1,file='walkalize',form='formatted')
        write(1,'(''#      i          f       w_abs_gen       e_gen         nwalk'')')
    endif
!---
!    open(1,file='walkalize',form='formatted')
!    write(1,'(''#      i          f         w_gen     e_gen_av    nwalk'')')
!***e
  endif

! We can calculate wavefn only when the full space is not large.
  if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
    allocate(wavefn(ndet),stat=istat)
    write(6,'(''ndet='',i8)') ndet
    if(istat.ne.0) stop 'failed to allocate wavefn'
  else
    allocate(wavefn(1),stat=istat) ! Done only to get rid of warning msg.
    if(istat.ne.0) stop 'failed to allocate wavefn'
  endif

  ! Edited by AAH on 20 Feb 2013. Skip this very slow and unnecessary step if we have many occupied dets, e.g., when we are using hf_to_psit!
! Calculate local energy pieces just to check.  Not needed for run
!      do iwalk=1,min(nwalk,ndet) ! For fictitious, read and hubbard we take det=walk.  For heg, chem we take walkers created at start.
!!       idet=(maxloc(iwdet_psi_t,1,iwdet_psi_t.eq.walk_det(iwalk)))
!        e_den=0 ; e_num=0
!        if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
!          do idet=1,ndet_psi_t
!            if(iwdet_psi_t(idet).eq.walk_det(iwalk)) then
!              if(importance_sampling.eq.0) then
!                e_den=e_den+cdet_psi_t(idet)
!              else
!                e_den=e_den+cdet_psi_t(idet)/psi_g(iwdet_psi_t(idet))
!              endif
!            endif
!            if(importance_sampling.eq.0) then
!              e_num=e_num+ham(walk_det(iwalk),iwdet_psi_t(idet)) * cdet_psi_t(idet)
!            else
!              e_num=e_num+ham(walk_det(iwalk),iwdet_psi_t(idet)) * (cdet_psi_t(idet) / psi_g(walk_det(iwalk)))
!            endif
!          enddo ! idet
!          !write(6,'(''iwalk, psi_g(walk_det(iwalk)), e_num, e_den, e_loc='',i5,9f12.5)') iwalk, psi_g(walk_det(iwalk)), e_num, e_den, e_num/e_den
!        elseif((hamiltonian_type.eq.'heg') .or. (hamiltonian_type.eq.'chem') .or. (hamiltonian_type .eq. 'hubbard2') .or. (hamiltonian_type .eq. 'hubbardk') .or. (hamiltonian_type .eq. 'hubbarddm')) then
!          if(hamiltonian_type.eq.'heg') then
!            call energy_pieces_heg(walk_dets_up(iwalk), walk_dets_dn(iwalk), e_mix_numerator, e_mix_denominator)
!          elseif(hamiltonian_type.eq.'chem') then
!            call energy_pieces_chem(walk_dets_up(iwalk), walk_dets_dn(iwalk), e_mix_numerator, e_mix_denominator)
!          elseif(hamiltonian_type.eq.'hubbard2') then
!           call energy_pieces_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), e_mix_numerator, e_mix_denominator, importance_sampling)
!          elseif(hamiltonian_type.eq.'hubbardk' .and. importance_sampling .eq. 0) then
!            call energy_pieces_hubbard_k(walk_dets_up(iwalk),walk_dets_dn(iwalk),e_mix_numerator,e_mix_denominator,importance_sampling)
!            ! Don't need to worry about retrieving saved energy pieces here because this is not during run.
!          elseif(hamiltonian_type.eq.'hubbardk' .and. importance_sampling .eq. 1) then
!            call energy_pieces_hubbard_k(walk_dets_up(iwalk), walk_dets_dn(iwalk), e_mix_numerator, e_mix_denominator,importance_sampling)
!          elseif(hamiltonian_type.eq.'hubbarddm') then
!            call energy_pieces_hubbard_dm(walk_dets_up(iwalk), walk_dets_dn(iwalk), e_mix_numerator, e_mix_denominator)
!          endif
!
!          if(importance_sampling.eq.0) then
!            e_den=e_den+e_mix_denominator
!            e_num=e_num+e_mix_numerator
!            tmp_ovlp=tmp_ovlp+e_den
!          elseif ((importance_sampling.eq.1) .and. (hamiltonian_type .eq. 'hubbard2' .or. hamiltonian_type .eq. 'hubbardk')) then
!            e_den=e_den+e_mix_denominator
!            e_num=e_num+e_mix_numerator
!            tmp_ovlp=tmp_ovlp+e_den
!          else
!            stop 'Error: importance sampling not yet implemented for heg, chem, hubbard2, hubbardk'
!          endif ! importance_sampling
!        endif ! hamiltonian_type
!        ! Line below bypassed for Hubbard2 by HJC - Seg fault if psi_g not specified ?
!        if ((hamiltonian_type .ne.'hubbard2') .and. (hamiltonian_type .ne. 'hubbardk') .and. (hamiltonian_type .ne. 'hubbarddm')) then
!          write(6,'(''iwalk, psi_g(walk_det(iwalk)), e_num, e_den, e_loc='',i5,9f12.5)') iwalk, psi_g(walk_det(iwalk)), e_num, e_den, e_num/e_den
!          call flush(6)
!        endif
!     enddo ! iwalk
!
!   Calculate overlap of initial walkers with psi_t
    !=================================================================================
    ! Comment on April 28 2011 by Hitesh J. Changlani
    ! For the Gutzwiller HF calculations HJC believes overlap_with_psi_t=e_den
    ! Hence HJC bypasses the steps below for this case
    !=================================================================================

  call my_second(1,'assign locations')
!***b
  ! Moved here by AAH on 21 Feb 2013.
  if(semistochastic) then
    ! Store locations of deterministic dets in list of occupied dets.
    !and shuffles the rows of hamiltonian to make it consistent with the division of imp_dets in different cores
    ! imp_core_mask(1:ncore) = number of deterministic states on each core
    ! imp_mask(1:my_nimp) = indices of the deterministic states on this core
    i=0
    imp_core_mask=0
    imp_mask=0
    my_nimp=0
    do icore=0,ncores-1
      n=0
      do iwalk=1,int(MWALK) !ndet_psi_t_connected
!write (6,*) "imp_dist=",imp_distance(iwalk); call flush(6)
        if (imp_distance(iwalk).eq.0) then
          n=n+1
          if(icore==get_det_owner(walk_dets_up(iwalk),walk_dets_dn(iwalk))) then
            i=i+1
            locations_of_imp_dets(i) = n
            imp_core_mask(icore+1) = imp_core_mask(icore+1) + 1
            if(whoami==icore) then
              my_nimp = my_nimp + 1
              imp_mask(my_nimp) = n
            endif
          endif
          if (n.ge.n_imp) then
            exit
          endif
        endif
      enddo
    enddo
    write(6,'(/,''Distribution of deterministic states on cores, (imp_core_mask)='',10000i6)') imp_core_mask(1:ncores)
    !write(6,'(''Indices of deterministic states on this core (imp_mask)='',10000i6)') imp_mask(1:n_imp)
    if (n.ne.n_imp) then
      write (6,*) "1 n=",n,"n_imp=",n_imp
      call flush(6)
      stop "locations of imp broken! 1"
    endif
    write(6,*) "n_imp: ",n_imp,"my_nimp: ", my_nimp,"whoami=",whoami; call flush(6)
    if (minval(imp_core_mask)==0) then
      call mpi_stop("Too few deterministic dets for this number of cores. Try running with a larger deterministic space!")
    endif


    !*** Resize Hamiltonian discarding dets not present locally and reorders rows for enhanced locality between cores
    if (ncores>1) call from_upper_triangular_to_band_shuffle(whoami,imp_core_mask,imp_mask,locations_of_imp_dets,n_imp,my_nimp,minus_tau_H_indices,minus_tau_H_nonzero_elements,minus_tau_H_values)
    deallocate(imp_mask) !not needed anymore
    imp_stride=sum(imp_core_mask(1:whoami))
  endif


  !***Added search for psi_t_connected dets since they can be distributed [30/7/13]
  !**** needed only when hf_to_psit = t
  if(hf_to_psit) then
    first_owner=0
    iown_first=.false.
    my_ndet_psi_t_connected=0
    own_first_arr(1)=10
    do iwalk=1,ndet_psi_t_connected
      if(whoami==get_det_owner(walk_dets_up(iwalk),walk_dets_dn(iwalk))) then
        my_ndet_psi_t_connected = my_ndet_psi_t_connected + 1
        if(iwalk==1) then
          iown_first=.true.
          first_owner=whoami
          own_first_arr(1)=1
        endif
      endif
    enddo
  endif
  call mpi_allred(first_owner)

  ! Added by AR 10 Dec 2013
  ! reduces the dimension of connected states list neglecting states not needed
  !** add a flag to exclude this when we will read these lists from checkpoint file
  if (.not.use_psit_con_in .and. hamiltonian_type .ne. 'hubbard2' .and.  hamiltonian_type .ne. 'fictitious' .and. hamiltonian_type .ne. 'read') then
    call reduce_connections(psi_t_connected_e_loc_num,psi_t_connected_e_loc_den,psi_t_connected_dets_up,psi_t_connected_dets_dn)
  endif

! At present MPI works for only a few hamiltonian_type.  Hence the following "if then"
  if (hamiltonian_type .ne. 'fictitious' .and. hamiltonian_type .ne. 'read') then

    ! Edited by AAH on 19 Dec 2012
    ! Store locations of psi_t in list of occupied dets,
    ! for replacing HF with Psi T.
    ! For now, assume that the deterministic space is contained in the connections to psi trial.
    iown_psit1=.false.
    psit1_owner=0
    my_ndet_psit=0
    n=0
    my_n_perm_init=0
    do while(n < ndet_psi_t)
      do iwalk=1,int(MWALK)
          if (walk_dets_up(iwalk)==dets_up_psi_t(n+1).and.walk_dets_dn(iwalk)==dets_dn_psi_t(n+1)) then
              locations_of_psit(n+1) = iwalk
              n=n+1
              if(whoami==get_det_owner(walk_dets_up(iwalk),walk_dets_dn(iwalk))) then
                  my_ndet_psit = my_ndet_psit + 1
                  ndet_psit_mask(my_ndet_psit) = n
                  if(n==1) then
                      iown_psit1 = .true.
                      psit1_owner = whoami
                  endif
                  if(.not.hf_to_psit) then
                    if(abs(abs(cdet_psi_t(n))-cdet_psi_t_abs_max).lt.1.d-3.and..not.use_psit_con_in) then ! Only those that are nearly equal and of the same sign are deemed permanent initiators
                      my_n_perm_init=my_n_perm_init+1
!                      if(n_permanent_initiator > M_PERMANENT_INITIATOR) then
!                          write(6,'(''n_permanent_initiator > M_PERMANENT_INITIATOR'')')
!                          stop 'n_permanent_initiator > M_PERMANENT_INITIATOR'
!                      endif
!                      !sign_permanent_initiator(n_permanent_initiator)=int(sign(1._rk,cdet_psi_t(idet)))            ! Correction - HJC - Nov 4
!                      if (importance_sampling .eq. 0) then
!                          sign_permanent_initiator(n_permanent_initiator)=int(sign(1._rk,cdet_psi_t(n)),i1b)
!                      else
!                          sign_permanent_initiator(n_permanent_initiator)=1                                           ! Correction - sign is 1 if importance sampling is used!
!                      endif

                      write(6,*) "found perm initiator: ", whoami,n_permanent_initiator,sign_permanent_initiator(n_permanent_initiator),walk_wt(iwalk)
                    endif
                  endif
              endif
              if (n.ge.ndet_psi_t) then
                  exit
              endif
          endif
      enddo
    enddo
    call mpi_allred(psit1_owner)
    if (n.ne.ndet_psi_t) then
      write (6,*) "n=",n,"ndet_psi_t=",ndet_psi_t,"nwalk=",nwalk,", ndet_psi_t_connected=",ndet_psi_t_connected,"1"
      call flush(6)
      stop "locations of psi t broken!"
    endif
    ! End of edit
    ! End of move

!***   added by AR: allocating deltaw buffer
!   if(ncores>1) then
     allocate(deltaw(n_imp+my_ndet_psi_t_connected+my_ndet_psit+3),stat=istat)
     if(istat.ne.0) stop 'failed to allocate deltaw'
!   else
!     allocate(deltaw(n_imp+ndet_psi_t_connected+ndet_psi_t),stat=istat)
!     if(istat.ne.0) stop 'failed to allocate deltaw'
!   endif

!***b   master core distributes population
    my_nwlk_1=1
    my_cnt=0
    if(master_core) my_cnt=nwalk
    my_nwalk=0
    my_skip=0
    if (my_cnt>MWALK) then
      write (6,*) "STOP: Need to set MWALK to at least",my_cnt; call flush(6)
      stop "Need to set MWALK higher!"
    endif
    call mpi_sendwalks(my_skip,my_nwalk,my_cnt,my_nwlk_1,initiator,e_num_walker,e_den_walker,matrix_elements)
    if(ncores == 1) my_nwalk=nwalk !**needed if compiled without MPI libraries
!***e master core distributes population

    !***b generate lists
    if(semistochastic) then
        n=0
        do iwalk=1,my_nwalk
            if (imp_distance(iwalk).eq.0) then
                n=n+1
                my_locations_of_imp_dets(n) = iwalk !local address
                my_imp_wt(n) = walk_wt(iwalk)
                if (n.ge.my_nimp) then
                    exit
                endif
            endif
        enddo
        if (n.ne.my_nimp) then
            write (6,*) "1 n=",n,"my_nimp=",my_nimp,"coreid=",whoami
            call flush(6)
            call mpi_stop("locations of my imp broken! 2")
        endif
        call mpi_allred(n)
        if (n.ne.n_imp) then
            write (6,*) "1 n=",n,"my_nimp=",my_nimp,"n_imp=",n_imp
            call flush(6)
            stop "locations of imp broken! 3"
        endif
    endif
    !***psi_t_connected states will be the first my_ndet_psi_t_connected(whoami+1) ones in each core
    !***--> no need to save locations
    !***should be possible with just one loop since they are already sorted [TODO]
    do n=1,my_ndet_psit
        do iwalk=1,my_nwalk
            if (walk_dets_up(iwalk)==dets_up_psi_t(ndet_psit_mask(n)).and.walk_dets_dn(iwalk)==dets_dn_psi_t(ndet_psit_mask(n))) then
                my_locations_of_psit(n) = iwalk
                exit
            endif
        enddo
    enddo
    n=my_ndet_psit
    call mpi_allred(n)
    if (n.ne.ndet_psi_t) then
        write (6,*) "nwalk=",nwalk,", ndet_psi_t_connected=",ndet_psi_t_connected
        call flush(6)
        stop "locations of psi t broken!"
    endif
    !***e generate lists


!--------------------------
!    ! Moved here by AAH on 21 Feb 2013.
!    if(semistochastic) then
!      ! Store locations of deterministic dets in list of occupied dets.
!      n=0
!      do iwalk=1,int(MWALK) !ndet_psi_t_connected
!        if (imp_distance(iwalk).eq.0) then
!          n=n+1
!          locations_of_imp_dets(n)=iwalk
!          if (n.ge.n_imp) then
!            exit
!          endif
!        endif
!      enddo
!      if (n.ne.n_imp) then
!        write (6,*) "1 n=",n,"n_imp=",n_imp
!        call flush(6)
!        stop "locations of imp broken!"
!      endif
!    endif
!    ! Edited by AAH on 19 Dec 2012
!    ! Store locations of psi_t in list of occupied dets,
!    ! for replacing HF with Psi T.
!    ! For now, assume that the deterministic space is contained in the connections to psi trial.
!    n=0
!    do iwalk=1,ndet_psi_t_connected
!          if (walk_dets_up(iwalk)==dets_up_psi_t(n+1).and.walk_dets_dn(iwalk)==dets_dn_psi_t(n+1)) then
!              locations_of_psit(n+1)=iwalk
!              n=n+1
!              if (n.ge.ndet_psi_t) then
!                  exit
!              endif
!          endif
!   enddo
!   if (n.ne.ndet_psi_t) then
!     write (6,*) "nwalk=",nwalk,", ndet_psi_t_connected=",ndet_psi_t_connected,"n=",n
!     call flush(6)
!     stop "locations of psi t broken!"
!   endif
!   ! End of edit
!   ! End of move
!***e

  else
    my_nwalk=nwalk
  endif ! hamiltonian_type

  call my_second(1,'calculate overlap with psi t'); call flush(6)

  overlap_with_psi_t=0
!***b
  if ((hamiltonian_type .eq. 'hubbard2') .and. wf_type .eq. 'gutz') then            !If - else statement by HJC on May 06 2011
      overlap_with_psi_t=tmp_ovlp
  else
      if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
      !***WARNING: for now parallelization works just for chem [TODO]
        if(ncores > 1) stop "for now parallelization works just for chem"
        k=1
        do iwalk=1,nwalk
          do idet=k,ndet_psi_t
            if(walk_det(iwalk).eq.iwdet_psi_t(idet)) then
              overlap_with_psi_t=overlap_with_psi_t+walk_wt(iwalk)*cdet_psi_t(idet)
            endif
          enddo
        enddo
      elseif((hamiltonian_type.eq.'heg') .or. (hamiltonian_type.eq.'chem') .or. (hamiltonian_type .eq. 'hubbard2') .or. (hamiltonian_type .eq. 'hubbardk') .or. (hamiltonian_type .eq. 'hubbarddm')) then
        ! Edited by AAH on 22 Feb 2013. Removed double for-loop.

        do idet=1,my_ndet_psit
          if(dets_up_psi_t(ndet_psit_mask(idet)).eq.walk_dets_up(my_locations_of_psit(idet)) .and. dets_dn_psi_t(ndet_psit_mask(idet)).eq.walk_dets_dn(my_locations_of_psit(idet))) then
            overlap_with_psi_t=overlap_with_psi_t+walk_wt(my_locations_of_psit(idet))*cdet_psi_t(ndet_psit_mask(idet))
          endif
        enddo
        call mpi_allred(overlap_with_psi_t)
      endif
      tmp_real=sum(real(walk_wt(1:my_nwalk),rk)**2)
!write(6,'(''my_nwalk, nwalk, ndet_psi_t, tmp_real'',3i5,es12.4)') my_nwalk, nwalk, ndet_psi_t, tmp_real
      call mpi_allred(tmp_real)
      overlap_with_psi_t=overlap_with_psi_t/sqrt(tmp_real)
  endif
!---
!  if ((hamiltonian_type .eq. 'hubbard2') .and. wf_type .eq. 'gutz') then            !If - else statement by HJC on May 06 2011
!      overlap_with_psi_t=tmp_ovlp
!  else
!      if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
!        k=1
!        do iwalk=1,nwalk
!          do idet=k,ndet_psi_t
!            if(walk_det(iwalk).eq.iwdet_psi_t(idet)) then
!              overlap_with_psi_t=overlap_with_psi_t+walk_wt(iwalk)*cdet_psi_t(idet)
!            endif
!          enddo
!        enddo
!      elseif((hamiltonian_type.eq.'heg') .or. (hamiltonian_type.eq.'chem') .or. (hamiltonian_type .eq. 'hubbard2') .or. (hamiltonian_type .eq. 'hubbardk') .or. (hamiltonian_type .eq. 'hubbarddm')) then
!        ! Edited by AAH on 22 Feb 2013. Removed double for-loop.
!        do idet=1,ndet_psi_t
!          if(dets_up_psi_t(idet).eq.walk_dets_up(locations_of_psit(idet)) .and. dets_dn_psi_t(idet).eq.walk_dets_dn(locations_of_psit(idet))) then
!            overlap_with_psi_t=overlap_with_psi_t+walk_wt(locations_of_psit(idet))*cdet_psi_t(idet)
!          endif
!        enddo
!      endif
!      overlap_with_psi_t=overlap_with_psi_t/sqrt(sum(real(walk_wt(1:nwalk),rk)**2))
!  endif
!***e


  write(6,'(''Overlap of initial walkers with psi_t='',f9.6,/)') overlap_with_psi_t
  call my_second(2,'calculate overlap with psi t')
  call flush(6)

  if(overlap_with_psi_t.lt.0 .and. (run_type.ne.'fixed_node2') .and. (run_type .ne. 'sr')) then
    walk_wt=-walk_wt
    overlap_with_psi_t=-overlap_with_psi_t
  endif
  overlap_with_psi_t_sav=overlap_with_psi_t

! Walk the walk
  n_sign_flips=0
  n_dominant_det=0
! walk_wt_projection_det_sav(1)=minval(walk_wt, mask=walk_det.eq.projection_det)
  if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
    wavefn=0
  endif

! Save quantities for using a larger tau and r_initiator during most of the equilibration
  tau_sav=tau
  tau_previous=tau
  r_initiator_sav=r_initiator

  e_est=e_trial
  proposal_prob_inv=ndet-1
  reweight_factor_inv=1
  iblkk=1
  ntimes_nblk_eq=1
  reached_w_abs_gen=0
  equilibrated_w_abs_blk=0
  equilibrated_w_perm_initiator_blk=0
  equilibrated_e_blk=0
  equilibrated_all=0
  e_blk_prev=1.d300

!***b
    w_abs_gen=sum(abs(walk_wt(1:my_nwalk)))
    call mpi_allred(w_abs_gen)
!    w_abs_gen=sum(abs(walk_wt(1:nwalk)))
!***e

    call my_second(1,'1st equilibration run')
    call flush(6)

    do while (iblkk <= ntimes_nblk_eq*nblk_eq+nblk)
! At present we reset e_est after each step but we keep e_trial fixed after the equilibration blocks are done.
! When e_trial is kept fixed, then the reweight factor corrects for inaccuracy in e_trial, but when it is changed it must depend only on the population.
! Zero out quantities after each set of equilibration blocks
        if((nblk_eq.eq.1 .or. mod(iblkk,nblk_eq).eq.1) .and. iblkk.le.ntimes_nblk_eq*nblk_eq+1) then
            if((nblk_eq.eq.1 .or. mod(iblkk,nblk_eq).eq.1) .and. iblkk.gt.nblk_eq) then
                i_equilibration_run=iblkk/nblk_eq
!***b
                write(6,'(/,''During equilibration run'',i3,'', n_sign_flips, n_sign_flips per unit time='',i6,f10.6)') i_equilibration_run, n_sign_flips, real(n_sign_flips,rk)/(real(nstep,rk)*nblk_eq*tau)
                write(write_string,'(''equilibration run'',i4)') i_equilibration_run
!----
!                write(6,'(/,''During equilibration run'',i3,'', n_sign_flips, n_sign_flips per unit time='',i6,f10.6)') i_equilibration_run, n_sign_flips, real(n_sign_flips,rk)/(real(nstep,rk)*nblk_eq*tau)
!                write(write_string,'(''equilibration run'',i4)') i_equilibration_run
!***e
                call my_second(2,write_string) ; write(6,*)
            endif

!            if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') write(6,'(''Dominant det is:'',i10)') dominant_det_most

            n_sign_flips=0
            n_dominant_det=0
            if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') wavefn=0
! Changing e_trial is now done near the end of the istep loop, so comment out next few lines
!           if(e_den_blkabs_cum.ne.0) then
!               e_trial=e_num_blkabs_cum/e_den_blkabs_cum
!           endif
            nwalk_cum=0 ; nwalk_before_merge_cum=0
            w_perm_initiator_cum=0 ; w_cum=0     ; w2_cum=0     ; w_gen2_cum=0     ; w_blk2_cum=0     ; w_abs_cum=0     ; w_genabs_cum=0     ; w_blkabs_cum=0     ; w_abs_before_merge_cum=0
            w_abs_genlast_cum=0
            e_num_cum=0 ; e_num2_cum=0 ; e_num_gen2_cum=0 ; e_num_blk2_cum=0 ; e_num_abs_cum=0 ; e_num_genabs_cum=0 ; e_num_blkabs_cum=0
            e_den_cum=0 ; e_den2_cum=0 ; e_den_gen2_cum=0 ; e_den_blk2_cum=0 ; e_den_abs_cum=0 ; e_den_genabs_cum=0 ; e_den_blkabs_cum=0
            e_num_e_den_cum=0 ; e_num_gen_e_den_gen_cum=0 ; e_num_blk_e_den_blk_cum=0
            e_gen_cum=0 ; e_gen2_cum=0 ; e_genp_e_genpp_cum=0 ; e_genp_cum=0 ; e_genpp_cum=0
            e_num_gen_ave=0 ; e_num_gen_vn1=0 ; e_num_blk_ave=0 ; e_num_blk_vn1=0
            e_den_gen_ave=0 ; e_den_gen_vn1=0 ; e_den_blk_ave=0 ; e_den_blk_vn1=0
            e_num_genabs_ave=0 ; e_num_genabs_vn1=0 ; e_num_blkabs_ave=0 ; e_num_blkabs_vn1=0
            e_den_genabs_ave=0 ; e_den_genabs_vn1=0 ; e_den_blkabs_ave=0 ; e_den_blkabs_vn1=0
            e_num_genabs_e_den_genabs_vn1 =0 ; e_num_blkabs_e_den_blkabs_vn1=0
            e_gen_ave=0 ; e_gen_vn1=0 ; e_genpp_ave=0 ; e_genp_e_genpp_vn1=0
            e_est=0 ! needed if growth estimator is used for e_est
            imp_distance_fraction(-2:max_imp_distance)=0
        endif

        if(equilibrated_all.eq.2 .and. iblkk.gt.ntimes_nblk_eq*nblk_eq) then
            iblk=iblkk-ntimes_nblk_eq*nblk_eq
        else
            iblk=iblkk-(ntimes_nblk_eq-1)*nblk_eq
        endif

        w_perm_initiator_blk=0 ; w_abs_blk=0 ; w_abs_blk_imp=0 ; w_blk=0 ; w2_blk=0 ; e_num_blk=0 ; e_den_blk=0 ; e_num_blk_e_den_blk=0


!** edit by AR [28/11/13] additional auxiliary variables for storing LOCAL estimates
        my_w2_cum=0d0; my_e_den2_cum=0d0; my_e_den_abs_cum=0d0; my_e_num_abs_cum=0d0; my_e_num2_cum=0d0; my_e_num_e_den_cum=0d0; my_w_abs_before_merge_cum=0d0
        my_nwalk_before_merge_cum = 0; my_nwalk_cum=0; my_w_perm_initiator_blk=0; my_w_perm_initiator_cum=0; my_w_abs_blk_imp=0
!*** end of edit


        do istep=1,nstep
            w_abs_gen_prev=w_abs_gen


            if(reached_w_abs_gen.eq.0) then
                tau=tau_sav*(1+log(w_abs_gen_target/w_abs_gen))
                tau_ratio=tau/tau_previous
                r_initiator=r_initiator_sav*(1+log(w_abs_gen_target/w_abs_gen))**initiator_rescale_power
                if(semistochastic) minus_tau_H_values=minus_tau_H_values*tau_ratio
!***b
                if(iblkk.eq.1 .and. istep.eq.1 .and. master_core)  write(6,'(/,''Variable tau and initiator used for equilibration.  At first step, tau, r_initiator='',2f10.6,/)') tau, r_initiator
!                if(iblkk.eq.1 .and. istep.eq.1) write(6,'(/,''Variable tau and initiator used for equilibration.  At first step, tau, r_initiator='',2f10.6,/)') tau, r_initiator
!***e
            endif


!***b
            if(semistochastic.and..not.hf_to_psit) then
            ! Store locations of deterministic dets in list of occupied dets.
                n=0
                do iwalk=1,my_nwalk
                    if (imp_distance(iwalk).eq.0) then
                        n=n+1
                        my_locations_of_imp_dets(n) = iwalk !local address
                        my_imp_wt(n) = walk_wt(iwalk)
                        if (n.ge.my_nimp) then
                            exit
                        endif
                    endif
                enddo
                if (n.ne.my_nimp) then
                    write (6,*) "1 n=",n,"my_nimp=",my_nimp,"coreid=",whoami
                    call flush(6)
                   stop "locations of my imp broken! 4"
                endif
!                call mpi_allred(n) !avoid this for performance later [TODO]
!                if (n.ne.n_imp) then
!                    write (6,*) "1 n=",n,"my_nimp=",my_nimp,"n_imp=",n_imp
!                    call flush(6)
!                    stop "locations of imp broken!"
!                endif
            endif


! Propagate the walkers
      if (use_improved_cancellation) then
        if(ncores>1) stop "improved_cancellation not tested yet!"
        call move_uniform_improved_cancellation(my_nwalk)
      else
        snd_cnt=0
        ! No need to spawn with HF det stochastically
        if (whoami==hf_core .and. hamiltonian_type.ne.'fictitious' .and. hamiltonian_type.ne.'read' .and. hamiltonian_type.ne.'hubbard' .and. hamiltonian_type.ne.'hubbard2') then
          do iwalk=1,my_nwalk
         !do iwalk=2,my_nwalk
            call move(iwalk)
          enddo ! iwalk
        else
          do iwalk=1,my_nwalk
            call move(iwalk)
          enddo ! iwalk
        endif
      endif
      my_new_walks=sum(snd_cnt)


      if(ncores>1) then
          call mpi_sendnewwalks(my_nwalk,initiator,e_num_walker,e_den_walker,matrix_elements)
      else
          my_nwalk = my_nwalk + my_new_walks
      endif

      if(ipr.ge.1) call write_walker_info('Before sorting,',my_nwalk)


!     Sort
      if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
        !***WARNING not tested
        if(ncores>1) stop 'just chem for now'
!       call shell_sort(walk_det,walk_wt,nwalk)
!       call merge_sort(walk_det,walk_wt,initiator,nwalk,temp_i_1,temp_i_2,temp_i1)
        call sort_walkers(my_nwalk)
      else
!        ! Edited by AAH on Dec 21, 2011. This merges the original with the deterministically "spawned" walkers.
!        ! If important space is being used, do the important space projection deterministically.
        if (.not.use_exponential_projector) then
            if(semistochastic) then
                if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem') then
                    if(ncores >1) then
                        call fast_sparse_matrix_multiply_local_band(n_imp,matrix_indices=minus_tau_H_indices,nelem_nonzero=minus_tau_H_nonzero_elements,matrix_values=minus_tau_H_values,vector=walk_wt(my_locations_of_imp_dets(1:my_nimp)),answer=deltaw)
                        call mpi_redscatt_real_dparray(deltaw(1:n_imp),int(imp_core_mask))
                    else
                        call fast_sparse_matrix_multiply_upper_triangular(n_imp,matrix_indices=minus_tau_H_indices,nelem_nonzero=minus_tau_H_nonzero_elements,matrix_values=minus_tau_H_values,vector=walk_wt(my_locations_of_imp_dets(1:n_imp)),answer=deltaw)
                    endif
                else
                    !***WARNING not tested
                    if(ncores>1) stop "just 'chem' hamiltonian for now"
                    !***WARNING not tested
                    call fast_sparse_matrix_multiply(n_imp,minus_tau_H_indices,minus_tau_H_nonzero_elements,minus_tau_H_values,imp_wt,deltaw)
                endif

                ! Edited by AAH on 4 Jan 2013. Multiply by first row and column of transformed projector.
                if (hf_to_psit) then
                    if(ncores>1) then
                        if(iown_first) wt1=walk_wt(1)
                        call mpi_bsend_diffroot(wt1,first_owner)
                        call sparse_fast_sparse_matrix_multiply_upper_triangular(own_first_arr,ndet_psi_t_connected,my_ndet_psi_t_connected,matrix_values=psi_t_connected_e_loc_num,vector=walk_wt(1:my_ndet_psi_t_connected),answer=deltaw(n_imp+1:n_imp+my_ndet_psi_t_connected+1),diag_elems=diag_elems,den1=psi_t_connected_e_loc_den(1),v1=wt1)
                        call mpi_red(deltaw(n_imp+1),first_owner)
                        if(iown_first) deltaw(n_imp+2)=deltaw(n_imp+1)

                        if(iown_psit1) wt1=walk_wt(my_locations_of_psit(1))
                        call mpi_bsend_diffroot(wt1,psit1_owner)
                        call sparse_fast_sparse_matrix_multiply_upper_triangular(ndet_psit_mask,ndet_psi_t,my_ndet_psit,matrix_values=(1._rk+tau*e_trial)*cdet_psi_t,vector=walk_wt(my_locations_of_psit(1:my_ndet_psit)),answer=deltaw(n_imp+my_ndet_psi_t_connected+2:n_imp+my_ndet_psi_t_connected+my_ndet_psit+2),v1=wt1)
                        call mpi_red(deltaw(n_imp+my_ndet_psi_t_connected+2),psit1_owner)
                        if(iown_psit1) deltaw(n_imp+my_ndet_psi_t_connected+3)=deltaw(n_imp+my_ndet_psi_t_connected+2)
                    else
                        call fast_sparse_matrix_multiply_upper_triangular(ndet_psi_t_connected,matrix_values=psi_t_connected_e_loc_num,vector=walk_wt(1:ndet_psi_t_connected),answer=deltaw(n_imp+1:n_imp+ndet_psi_t_connected),diag_elems=diag_elems,den1=psi_t_connected_e_loc_den(1))
                        call fast_sparse_matrix_multiply_upper_triangular(ndet_psi_t,matrix_values=(1._rk+tau*e_trial)*cdet_psi_t,vector=walk_wt(my_locations_of_psit(1:ndet_psi_t)),answer=deltaw(n_imp+ndet_psi_t_connected+1:n_imp+ndet_psi_t_connected+ndet_psi_t))
                    endif
                else
                     deltaw(1:my_nimp) = deltaw(1:my_nimp)+e_trial*tau*my_imp_wt(1:my_nimp)
                endif
            ! End of edit

            endif
        endif
!        if(ncores>1) then
!            call mpi_sendnewwalks(my_nwalk,initiator)
!        else
!            my_nwalk = my_nwalk + my_new_walks
!        endif
        if (.not.use_exponential_projector) then
            if(semistochastic) then
!                if(ncores>1) call mpi_redscatt_real_dparray(deltaw(1:n_imp),imp_core_mask)
                if (hf_to_psit) then
                    if(ncores > 1) then
                        do i=1,my_ndet_psi_t_connected
                            walk_wt(i) = walk_wt(i) - tau*deltaw(n_imp+i+1) + tau*e_trial*walk_wt(i)
                        enddo
                        do i=1,my_ndet_psit
                            walk_wt(my_locations_of_psit(i)) = walk_wt(my_locations_of_psit(i)) + deltaw(n_imp+my_ndet_psi_t_connected+i+2)
                        enddo
                    else
                        do i=1,my_ndet_psi_t_connected
                            walk_wt(i) = walk_wt(i) - tau*deltaw(n_imp+i) + tau*e_trial*walk_wt(i)
                        enddo
                        do i=1,my_ndet_psit
                            walk_wt(my_locations_of_psit(i)) = walk_wt(my_locations_of_psit(i)) + deltaw(n_imp+ndet_psi_t_connected+ndet_psit_mask(i))
                        enddo
                    endif
                endif
                do i=1,my_nimp
                    walk_wt(my_locations_of_imp_dets(i)) = walk_wt(my_locations_of_imp_dets(i)) + deltaw(i)
                enddo
            endif
        endif
        ! End of edit

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Sort the determinants and assign the objects associated with them the right order
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (hf_to_psit) then
           call sort_my_walkers3_up_dn(my_nwalk,my_ndet_psi_t_connected,my_ndet_outside_ct)
        else
          call sort_my_walkers2_up_dn(my_nwalk)
        endif

      endif

      if(ipr.ge.1) call write_walker_info('After sorting,',my_nwalk)

      my_w_abs_before_merge_cum=my_w_abs_before_merge_cum+sum(abs(walk_wt(1:my_nwalk)))
!     w_abs_before_merge_cum=w_abs_before_merge_cum+sum(abs(walk_wt(1:nwalk)))

      !nwalk=my_nwalk
      nwalk_before_merge=my_nwalk

      ! Merge walkers that are on the same determinant
      ! Then join small weight walkers that are not on the same determinant
      if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
        !***b WARNING: for now parallelization works just for chem [TODO]
        if(ncores>1) stop 'for now parallelization works just for chem'
        !***e
        !walk_wt(1:nwalk) = walk_wt(1:nwalk) * reweight_factor_inv ! New comment HJC: reweighting for toy systems made consistent with non-toy systems.
        call merge_original_with_spawned(iblk,istep,walk_det,walk_wt,my_nwalk)
!write(6,'(''3 my_nwalk, nwalk='',9i5)') my_nwalk, nwalk
        call join_walker(my_nwalk)
        nwalk=my_nwalk
!write(6,'(''4 my_nwalk, nwalk='',9i5)') my_nwalk, nwalk
      else  ! For non - toy Hamiltonians i.e. heg,hubbard2,hubbardk,hubbarddm,chemistry
       !walk_wt(1:nwalk) = walk_wt(1:nwalk) * reweight_factor_inv ! Warning: Adam changed where mult be reweight_factor_inv is done.  Make sure it is still consistent for toy version.

        ! Do different sequence of things if SR or not - HJC, September 29,2012
        if (run_type .ne. 'sr') then ! If not stochastic reconfiguration
            ! Take the sorted determinants and merge them - that is, perform cancellations
            if (hf_to_psit) then
!***b
               call merge_my_original_with_spawned3(my_nwalk, my_ndet_psi_t_connected, my_ndet_outside_ct)
!              call merge_original_with_spawned3
!***e
            else
!***b
              call merge_original_with_spawned2(iblk,istep,walk_dets_up(1:my_nwalk),walk_dets_dn(1:my_nwalk),walk_wt(1:my_nwalk),my_nwalk)!,my_n_perm_init)
!              call merge_original_with_spawned2(iblk,istep,walk_dets_up,walk_dets_dn,walk_wt,nwalk)
!***e
            endif
        else
            !***b WARNING: for now parallelization may not work for sr [TODO]
            if(ncores>1) stop 'run_type=sr not implemented yet'
            !***e
            ! Take the sorted determinants and merge them - that is, perform cancellations
            call merge_original_with_spawned2(iblk,istep,walk_dets_up,walk_dets_dn,walk_wt,my_nwalk)!,n_permanent_initiator)

            ! Perform stochastic reconfiguration of walkers - Keeps number of walkers fixed
            call stochastic_reconfiguration

            ! Walkers have new weights - now set initiators
            ! No merging needed , but initiator rules of merge_original_with_spawned2 being used - might be inefficient
            ! Seems ok as a temporary measure
            call merge_original_with_spawned2(iblk,istep,walk_dets_up,walk_dets_dn,walk_wt,my_nwalk)!,n_permanent_initiator)
        endif

        ! Multiply by T^-1 transpose and T^-1 if hf_to_psit is true, after the merge but before the join
        if (hf_to_psit) then

           !n=1
           !do iwalk=1,int(MWALK)
           !  if (walk_dets_up(iwalk)==dets_up_psi_t(n).and.walk_dets_dn(iwalk)==dets_dn_psi_t(n)) then
           !    psit_wt(n) = walk_wt(iwalk)
           !    locations_of_psit(n) = iwalk
           !    if (n.ge.ndet_psi_t) then
           !      exit
           !    endif
           !    n=n+1
           !  endif
           !enddo
           !if (n.ne.ndet_psi_t) stop "locations of psi t broken!"
!****b
          ! Multiply by T^-1
          tmp_real=0
          if((.not.iown_psit1).and.(my_ndet_psit>0)) tmp_real=cdet_psi_t(ndet_psit_mask(1))*walk_wt(my_locations_of_psit(1))
          do i=2,my_ndet_psit
            tmp_real=tmp_real+cdet_psi_t(ndet_psit_mask(i))*walk_wt(my_locations_of_psit(i))
          enddo
          call mpi_red(tmp_real,psit1_owner)
          if(iown_psit1) then
               walk_wt(my_locations_of_psit(1)) = walk_wt(my_locations_of_psit(1)) - tmp_real
               walk_wt(my_locations_of_psit(1)) = walk_wt(my_locations_of_psit(1)) / cdet_psi_t(ndet_psit_mask(1))
          endif

          ! Multiply by transpose(T^-1)
          if(iown_psit1) then
            walk_wt(my_locations_of_psit(1)) = walk_wt(my_locations_of_psit(1)) / cdet_psi_t(ndet_psit_mask(1))
            tmp_real=walk_wt(my_locations_of_psit(1))
          endif
          call mpi_bsend_diffroot(tmp_real,psit1_owner)

          if((.not.iown_psit1).and.(my_ndet_psit>0)) walk_wt(my_locations_of_psit(1)) = walk_wt(my_locations_of_psit(1)) - cdet_psi_t(ndet_psit_mask(1)) * tmp_real
          do i=2,my_ndet_psit
            walk_wt(my_locations_of_psit(i)) = walk_wt(my_locations_of_psit(i)) - cdet_psi_t(ndet_psit_mask(i)) * tmp_real
          enddo
!---
!          do i=2,ndet_psi_t
!            walk_wt(locations_of_psit(1)) = walk_wt(locations_of_psit(1)) - cdet_psi_t(i) * walk_wt(locations_of_psit(i))
!          enddo
!          walk_wt(locations_of_psit(1)) = walk_wt(locations_of_psit(1)) / cdet_psi_t(1)
!          walk_wt(locations_of_psit(1)) = walk_wt(locations_of_psit(1)) / cdet_psi_t(1)
!          do i=2,ndet_psi_t
!            walk_wt(locations_of_psit(i)) = walk_wt(locations_of_psit(i)) - cdet_psi_t(i) * walk_wt(locations_of_psit(1))
!          enddo
!****e
        endif

        if (hf_to_psit) then

          i_permanent_initiator = 0
          if (.not.c_t_initiator) then
!***b
            do i=1,my_ndet_psi_t_connected
              call check_initiator(i,discard_walker,i_permanent_initiator) ! Don't discard walkers here because we keep all dets in ct.
            enddo
!---
!            do i=1,ndet_psi_t_connected
!              call check_initiator(i,discard_walker,i_permanent_initiator) ! Don't discard walkers here because we keep all dets in ct.
!            enddo
!***e
          else ! This is necessary only to allow us to set w_perm_initiator to 1 if it is too small.
            do i=1,n_permanent_initiator
              call check_initiator(i,discard_walker,i_permanent_initiator)
            enddo
          endif

        endif

        if(ipr.ge.1) then
          call write_walker_info('Before reduce_my_walker or join_walker2,',my_nwalk)
          if(iblkk.le.ntimes_nblk_eq*nblk_eq) write(6,'(''nwalk_before_merge, my_nwalk, ratio = '',2i8,f8.4)') nwalk_before_merge, my_nwalk, real(my_nwalk,rk)/nwalk_before_merge
        endif

! Either join walkers or reduce them (similar to integerization but with real wt.)
        if(semistochastic) then
         !call join_walker_semistoch(my_nwalk)
          call reduce_my_walker(my_nwalk,my_ndet_psi_t_connected,my_ndet_outside_ct)
        else
          call join_walker2(my_nwalk)
        endif
        nwalk=my_nwalk

        if(ipr.ge.1) then
          call write_walker_info('After merging, multiplying by S^-1, and joining,',my_nwalk)
          if(iblkk.le.ntimes_nblk_eq*nblk_eq) write(6,'(''nwalk_before_merge, my_nwalk, ratio = '',2i8,f8.4)') nwalk_before_merge, my_nwalk, real(my_nwalk,rk)/nwalk_before_merge
        endif

      endif ! hamiltonian_type

!     walk_wt(1:nwalk) = walk_wt(1:nwalk) * reweight_factor_inv
      walk_wt(1:my_nwalk) = walk_wt(1:my_nwalk) * reweight_factor_inv

      if(my_nwalk.eq.0) then
        write(6,'(''stop my_nwalk=0, istep, iblk, e_est, e_trial, reweight_factor=''2i9,9d12.4)') istep, iblk, e_est, e_trial, 1._rk/reweight_factor_inv
        call flush(6)
        stop 'my_nwalk=0'
      endif

! Accumulate the wavefunction.  This we can do only because we do not have a huge space.
! Determine the dominant basis function for this generation.  This we can do even in a huge space.
      if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
        dominant_det_abs_wt_tmp=0
        do iwalk=1,my_nwalk
          wavefn(walk_det(iwalk))=wavefn(walk_det(iwalk))+walk_wt(iwalk)
          if(abs(walk_wt(iwalk)).gt.dominant_det_abs_wt_tmp) then
            dominant_det=walk_det(iwalk)
            dominant_det_wt_tmp=walk_wt(iwalk)
            dominant_det_abs_wt_tmp=abs(walk_wt(iwalk))
          endif
        enddo

        do iwalk=1,n_dominant_det
          if(dominant_det.eq.dominant_det_list(iwalk)) goto 10
        enddo
        n_dominant_det=n_dominant_det+1
        if(n_dominant_det.gt.M_DOMINANT_DET) then
          write(6,'(''n_dominant_det > M_DOMINANT_DET'')')
          stop 'n_dominant_det > M_DOMINANT_DET'
        endif
        dominant_det_list(n_dominant_det)=dominant_det
        dominant_det_wt(n_dominant_det)=dominant_det_wt_tmp
        dominant_abs_det_wt(n_dominant_det)=abs(dominant_det_wt_tmp)
        goto 20
   10   dominant_det_wt(iwalk)=dominant_det_wt(iwalk)+dominant_det_wt_tmp
        dominant_abs_det_wt(iwalk)=dominant_abs_det_wt(iwalk)+abs(dominant_det_wt_tmp)
   20   continue

! Select the most dominant determinant.  At present I am doing it based on the wt rather than on the absolute wt.
        dominant_det_abs_wt_tmp=0
        do iwalk=1,n_dominant_det
          if(abs(dominant_det_wt(iwalk)).gt.dominant_det_abs_wt_tmp) then
            dominant_det_most=dominant_det_list(iwalk)
            dominant_det_abs_wt_tmp=abs(dominant_det_wt(iwalk))
          endif
        enddo

        if(ipr.ge.4) write(6,'(''dominant_det_most, dominant_det_abs_wt_tmp='',i10,f10.1)') dominant_det_most, dominant_det_abs_wt_tmp
        if(ipr.ge.4) write(6,'(''dominant_det_list  = '',100i10)')   (dominant_det_list(iwalk),iwalk=1,n_dominant_det)
        if(ipr.ge.4) write(6,'(''dominant_det_wt    = '',100f10.1)') (dominant_det_wt(iwalk),iwalk=1,n_dominant_det)
        if(ipr.ge.4) write(6,'(''dominant_abs_det_wt= '',100f10.1)') (dominant_abs_det_wt(iwalk),iwalk=1,n_dominant_det)

! Calculate norm of averaged wavefn
        if(ipr.ge.1) then
          norm_of_wavefn=0
          do idet=1,ndet
            !norm_of_wavefn=norm_of_wavefn+(dfloat(wavefn(idet)))**2
            norm_of_wavefn=norm_of_wavefn+(wavefn(idet))**2
          enddo
          norm_of_wavefn=sqrt(abs(norm_of_wavefn))
! Calculate overlap of averaged wavefn with exact one
          if(diagonalize_ham.eq.1) then
            overlap=0
            do idet=1,ndet
              overlap=overlap+wavefn_exact(idet)*wavefn(idet)
            enddo
            overlap=sqrt(abs(overlap)/norm_of_wavefn)
            if(ipr.ge.2) write(6,'(''wavefn_exact='',10000i5)') nint(100*sqrt(dfloat(ndet))*wavefn_exact)
          endif
          if(ipr.ge.2) write(6,'(''wavefn_proj ='',10000i5)') nint(100*sqrt(dfloat(ndet))*(wavefn/norm_of_wavefn))
        endif
      endif ! hamiltonian_type

        ! A temporary check to make sure number of important dets has not changed.
        !if (semistochastic) then  ! This if added by Hitesh needs to talk to Adam
        !    n_imp2=0
        !      do iwalk=1,nwalk
        !        if(imp_distance(iwalk).eq.0) n_imp2=n_imp2+1
        !      enddo
        !      if(n_imp2.ne.n_imp) then
        !        write(6,'(''before join: n_imp, n_imp2'',9i8)') n_imp, n_imp2
        !        stop 'before join: n_imp2 \= n_imp'
        !      endif
        !      call flush(6)
        !endif

! Zero out accumulators for this generation
      w_perm_initiator_gen=0 ; w_gen=0 ; w2_gen=0 ; w_abs_gen=0 ; e_num_gen=0 ; e_den_gen=0 ; e_num_gen_e_den_gen=0; w_abs_gen_imp=0

      my_nwalk_cum=my_nwalk_cum+my_nwalk
      my_nwalk_before_merge_cum=my_nwalk_before_merge_cum+nwalk_before_merge

      ! Added by HJC - tmp_ovlp

      tmp_ovlp=0._rk
      which_perm_initiator=0
      i_imp=0
      idet_begin=1

      if(ipr.ge.4) write(6,'(''my_nwalk,nwalk='',9i8)') my_nwalk,nwalk
      if(ipr.ge.4) write(6,'(''e_num_walker='',1000es14.6)') e_num_walker(1:my_nwalk)
      if(ipr.ge.4) write(6,'(''e_den_walker='',1000es14.6)') e_den_walker(1:my_nwalk)

      do iwalk=1,my_nwalk
!       idet=(maxloc(iwdet_psi_t,1,iwdet_psi_t.eq.walk_det(iwalk)))
!       if(initiator(iwalk).eq.3) w_perm_initiator_gen=w_perm_initiator_gen+walk_wt(iwalk)      ! Prior to Oct 20 2011 CJU
        if(initiator(iwalk).eq.3) then
            which_perm_initiator=which_perm_initiator+1                                         ! Assumption of sorted list used here
            w_perm_initiator_gen=w_perm_initiator_gen+walk_wt(iwalk)*sign_permanent_initiator(which_perm_initiator) ! Corrected by HJC on Oct 20 2011 - walk_wt*sign
        endif
        w_gen=w_gen+walk_wt(iwalk)
        w2_gen=w2_gen+walk_wt(iwalk)**2
        w_abs_gen=w_abs_gen+abs(walk_wt(iwalk))
        if (imp_distance(iwalk).eq.0.or.(imp_distance(iwalk).eq.-2.and.c_t_initiator))  w_abs_gen_imp=w_abs_gen_imp+abs(walk_wt(iwalk))
        e_den=0 ; e_num=0
        if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
            !***b WARNING not tested
            if(ncores>1) stop "MPI not yet tested for hamiltonian_type fictitious, read and hubbard"
            !***e WARNING not tested

            do idet=1,ndet_psi_t
                if(iwdet_psi_t(idet).eq.walk_det(iwalk)) then
                    if(importance_sampling.eq.0) then
                        e_den=e_den+walk_wt(iwalk)*cdet_psi_t(idet)
                    else
!                       if(hamiltonian_type.eq.'hubbard' .and. run_type.eq.'fixed_node2') then
!                           e_den=e_den+walk_wt(iwalk)*abs(cdet_psi_t(idet)/psi_g(iwdet_psi_t(idet)))
!                       else
                            e_den=e_den+walk_wt(iwalk)*cdet_psi_t(idet)/psi_g(iwdet_psi_t(idet))
!                       endif
                    endif
                endif
                if(importance_sampling.eq.0) then
                    e_num=e_num+ham(walk_det(iwalk),iwdet_psi_t(idet)) * cdet_psi_t(idet) * walk_wt(iwalk)
                else
!             if(hamiltonian_type.eq.'hubbard' .and. run_type.eq.'fixed_node2' .and. walk_det(iwalk).ne.iwdet_psi_t(idet)) then
!               e_num=e_num-abs(ham(walk_det(iwalk),iwdet_psi_t(idet)) * (cdet_psi_t(idet) / psi_g(walk_det(iwalk)))) * walk_wt(iwalk)
!              else
                e_num=e_num+ham(walk_det(iwalk),iwdet_psi_t(idet)) * (cdet_psi_t(idet) / psi_g(walk_det(iwalk))) * walk_wt(iwalk)
!             endif
            endif
          enddo ! idet
          if(ipr.ge.1) write(6,'(''iwalk, psi_g(walk_det(iwalk)), e_num, e_den, e_loc='',i5,9f12.5)') iwalk, psi_g(walk_det(iwalk)), e_num, e_den, e_num/e_den
        elseif((hamiltonian_type.eq.'heg') .or. (hamiltonian_type.eq.'chem') .or. (hamiltonian_type .eq. 'hubbard2') .or. (hamiltonian_type .eq. 'hubbardk') .or. (hamiltonian_type .eq. 'hubbarddm')) then

          !if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'chem') then
          if (hamiltonian_type.eq.'hubbard2') then
                !print *,"Here for hubbard2 energy pieces"
                !***b WARNING not tested
                if(ncores>1) stop "MPI not yet tested for hamiltonian_type hubbard2"
                !***e WARNING not tested
                !if ( (e_num_walker(iwalk) .ge. 1.e50_rk) .or. (e_den_walker(iwalk) .ge. 1.e50_rk)) then
                    call energy_pieces_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), e_num, e_den, importance_sampling)
                    e_num_walker(iwalk)=e_num
                    e_den_walker(iwalk)=e_den
                    !write (6,*) "Computed and now stored....",e_num_walker(iwalk),e_den_walker(iwalk)
                !else  ! Energy pieces already stored in list from previous computations, hence one can save time
                    !write (6,*) "Stored already...",e_num_walker(iwalk),e_den_walker(iwalk)
                    !call energy_pieces_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), e_num, e_den, importance_sampling)
                    !write (6,*) "Verifying ...",e_num, e_den
                !    e_num=e_num_walker(iwalk)
                !    e_den=e_den_walker(iwalk)
                !endif
                e_num=e_num*walk_wt(iwalk)
                e_den=e_den*walk_wt(iwalk)
                !print *,"e_num =",e_num
                !print *,"e_den =",e_den
                if(ipr.ge.1) write(6,'(''iwalk, e_num, e_den, e_loc='',i5,9f12.5)') iwalk, e_num, e_den, e_num/e_den
               !      Old style of doing energy pieces of heg commented out during Zach's A exam (Sept 2012). This part will be eventually removed
               !      elseif (hamiltonian_type .eq. 'heg') then
               !      call energy_pieces_heg(walk_dets_up(iwalk), walk_dets_dn(iwalk), e_num, e_den)
               !      e_num=e_num*walk_wt(iwalk)
               !      e_den=e_den*walk_wt(iwalk)
               !      if(ipr.ge.1) write(6,'(''iwalk, e_num, e_den, e_loc='',i5,9f12.5)') iwalk, e_num, e_den, e_num/e_den
           endif
          if(ipr.ge.1) write(6,'(''iwalk, e_num, e_den, e_loc='',i5,9f12.5)') iwalk, e_num, e_den, e_num/e_den
        endif ! hamiltonian_type
      ! Edited by AAH on 28 Mar 2012. Find local energies corresponding to each walker (so we can store more local energies than the size of the important space).
        if (.not.(hamiltonian_type.eq.'chem'.or.hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type .eq. 'heg')) then
          !***b WARNING not tested
          if(ncores>1) stop "MPI not yet tested for hamiltonian_type other than chem, hubbardk and heg"
          !***e WARNING not tested
          e_den_gen=e_den_gen+e_den
!*** edit by AR [28/11/13] -> updated at the end of block, use local variable for accumulation
          my_e_den2_cum=my_e_den2_cum+e_den**2
!           e_den2_cum=e_den2_cum+e_den**2
          my_e_den_abs_cum=my_e_den_abs_cum+abs(e_den)
!           e_den_abs_cum=e_den_abs_cum+abs(e_den)
          my_e_num_abs_cum=my_e_num_abs_cum+e_num*sign(1._rk,e_den)
!           e_num_abs_cum=e_num_abs_cum+e_num*sign(1._rk,e_den)
          my_e_num2_cum=my_e_num2_cum+e_num**2
!           e_num2_cum=e_num2_cum+e_num**2
          e_num_e_den=e_num*e_den
          my_e_num_e_den_cum=my_e_num_e_den_cum+e_num_e_den
!           e_num_e_den_cum=e_num_e_den_cum+e_num_e_den
!*** end of edit
          e_num_gen=e_num_gen+e_num
        endif
      enddo ! iwalk

!***b
      my_w_abs_gen=w_abs_gen
      if(ncores > 1) then
        tmp_real_array(1)=w_gen
        tmp_real_array(2)=w_abs_gen
      endif
      !e_num_gen=0  ! HJC comment May 2014
      !e_den_gen=0  ! HJC comment May 2014

!***e

      if ((hamiltonian_type.eq.'chem').or.(hamiltonian_type.eq.'hubbardk') .or. (hamiltonian_type .eq.'heg')) then     ! Wont work for wfs with amplitude on everything!!!
        e_num_gen=0  ! HJC comment May 2014 This was outside the "if" and hence overruling hubbard2 calcs
        e_den_gen=0  ! HJC comment May 2014
        ! Edited by AAH on 9 Apr 2013. No need for a linear search when we already know the locations of the psi_t connections
        if (hf_to_psit) then
!***b
            do j=1,my_ndet_psi_t_connected
                if (j==1.and.iown_first) then
                    e_num=psi_t_connected_e_loc_num(1)/psi_t_connected_e_loc_den(1)*walk_wt(1)
                    e_den=walk_wt(1)
                else
                    e_num=psi_t_connected_e_loc_num(j)*walk_wt(j)
                    e_den=psi_t_connected_e_loc_den(j)*walk_wt(j)
                endif
                if (abs(e_den).lt.1.d-22)  e_den=abs(e_den) ! This is necessary because the sign of e_num depends on the sign of e_den, which can be either +0 or -0
                e_den_gen=e_den_gen+e_den
                my_e_den2_cum=my_e_den2_cum+e_den**2
                my_e_den_abs_cum=my_e_den_abs_cum+abs(e_den)
                e_num_gen=e_num_gen+e_num
                my_e_num2_cum=my_e_num2_cum+e_num**2
                my_e_num_abs_cum=my_e_num_abs_cum+e_num*sign(1._rk,e_den)
                my_e_num_e_den_cum=my_e_num_e_den_cum+e_num*e_den
            enddo
          if(ncores>1) then
            tmp_real_array(3)=e_den_gen
            tmp_real_array(4)=e_num_gen
            tmp_real_array(5)=w_perm_initiator_gen
            tmp_real_array(6)=nwalk
            tmp_real_array(7)=w_abs_gen_imp
          endif
!------
!          do j=1,ndet_psi_t_connected
!            if (hf_to_psit.and.j==1) then
!              e_num=psi_t_connected_e_loc_num(1)/psi_t_connected_e_loc_den(1)*walk_wt(1)
!              e_den=walk_wt(1)
!            else
!              e_num=psi_t_connected_e_loc_num(j)*walk_wt(j)
!              e_den=psi_t_connected_e_loc_den(j)*walk_wt(j)
!            endif
!            if (abs(e_den).lt.1.d-22)  e_den=abs(e_den) ! This is necessary because the sign of e_num depends on the sign of e_den, which can be either +0 or -0
!            e_den_gen=e_den_gen+e_den
!            e_den2_cum=e_den2_cum+e_den**2
!            e_den_abs_cum=e_den_abs_cum+abs(e_den)
!            e_num_gen=e_num_gen+e_num
!            e_num2_cum=e_num2_cum+e_num**2
!            e_num_abs_cum=e_num_abs_cum+e_num*sign(1._rk,e_den)
!            e_num_e_den_cum=e_num_e_den_cum+e_num*e_den
!          enddo
!***e
         !e_den_gen = dot_product(walk_wt(1:ndet_psi_t_connected),psi_t_connected_e_loc_den)
         !e_den2_cum = dot_product(walk_wt(1:ndet_psi_t_connected)**2,psi_t_connected_e_loc_den**2)
         !e_den_abs_cum = dot_product(abs(walk_wt(1:ndet_psi_t_connected)),abs(psi_t_connected_e_loc_den))
         !e_num_gen = dot_product(walk_wt(1:ndet_psi_t_connected),psi_t_connected_e_loc_num)
         !e_num2_cum = dot_product(walk_wt(1:ndet_psi_t_connected)**2,psi_t_connected_e_loc_num**2)
         !e_num_abs_cum = dot_product(psi_t_connected_e_loc_num,sign(walk_wt(1:ndet_psi_t_connected),psi_t_connected_e_loc_den))
         !e_num_e_den_cum = dot_product(walk_wt(1:ndet_psi_t_connected)**2,psi_t_connected_e_loc_den*psi_t_connected_e_loc_num)
        else
!***b
          if (my_nwalk*log_n_con<my_nwalk+ndet_psi_t_connected) then
            call binary_search_list_and_update(walk_dets_up(1:my_nwalk),walk_dets_dn(1:my_nwalk),walk_wt(1:my_nwalk),psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den,e_num_gen,e_den_gen,my_e_num2_cum,my_e_den2_cum,my_e_num_abs_cum,my_e_den_abs_cum,my_e_num_e_den_cum,e_num_walker,e_den_walker)
          else
            call linear_search_list_and_update(walk_dets_up(1:my_nwalk),walk_dets_dn(1:my_nwalk),walk_wt(1:my_nwalk),psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den,e_num_gen,e_den_gen,my_e_num2_cum,my_e_den2_cum,my_e_num_abs_cum,my_e_den_abs_cum,my_e_num_e_den_cum,e_num_walker,e_den_walker)
          endif
          if(ncores>1) then
            tmp_real_array(3)=e_den_gen
            tmp_real_array(4)=e_num_gen
            tmp_real_array(5)=w_perm_initiator_gen
            tmp_real_array(6)=nwalk
            tmp_real_array(7)=w_abs_gen_imp
          endif
!----
!          call linear_search_list(walk_dets_up(1:nwalk),walk_dets_dn(1:nwalk),walk_wt(1:nwalk),psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den,e_num_gen,e_den_gen,e_num2_cum,e_den2_cum,e_num_abs_cum,e_den_abs_cum,e_num_e_den_cum)
!***e
        endif ! hf_to_psit
        ! End of edit
      endif ! hamiltonian_type
      ! End of edit

!***b

      if(ncores>1) then
        call mpi_allred(tmp_real_array(1:7))
        w_gen=tmp_real_array(1)
        w_abs_gen=tmp_real_array(2)
        e_den_gen=tmp_real_array(3)
        e_num_gen=tmp_real_array(4)
        w_perm_initiator_gen=tmp_real_array(5)
        disp_nwalk=nint(tmp_real_array(6))
        w_abs_gen_imp=tmp_real_array(7)
      else
        disp_nwalk=nwalk
      endif
!***e

      call flush(6)
      passes=(iblk-1)*nstep+istep
      w_abs_blk=w_abs_blk+w_abs_gen
      w_blk=w_blk+w_gen ; e_num_blk=e_num_blk+e_num_gen ; e_den_blk=e_den_blk+e_den_gen
      if(ipr.ge.1) write(6,'(''iblk,istep,e_num_gen,e_den_gen='',2i4,9f9.4)') iblk,istep,e_num_gen,e_den_gen

!*** edit by AR [28/11/13] -> updated at the end of block, use local variable for accumulation
!     my_w_abs_blk_imp=my_w_abs_blk_imp+w_abs_gen_imp
      w_abs_blk_imp=w_abs_blk_imp+w_abs_gen_imp
!     my_w_perm_initiator_blk=my_w_perm_initiator_blk+w_perm_initiator_gen
!     my_w_perm_initiator_cum=my_w_perm_initiator_cum+w_perm_initiator_gen
      w_perm_initiator_blk=w_perm_initiator_blk+w_perm_initiator_gen
      w_perm_initiator_cum=w_perm_initiator_cum+w_perm_initiator_gen
      my_w2_cum=my_w2_cum + w2_gen
!     w2_cum=w2_cum+w2_gen
!*** end of edit
      w_abs_cum=w_abs_cum+w_abs_gen
      w_gen2_cum=w_gen2_cum+w_gen**2
      e_den_gen2_cum=e_den_gen2_cum+e_den_gen**2
      e_num_gen2_cum=e_num_gen2_cum+e_num_gen**2
      w_genabs_cum=w_genabs_cum+abs(w_gen)
      e_den_genabs_cum=e_den_genabs_cum+abs(e_den_gen)
      e_num_genabs_cum=e_num_genabs_cum+e_num_gen*sign(1._rk,e_den_gen)
      e_num_gen_e_den_gen_cum=e_num_gen_e_den_gen_cum+e_num_gen*e_den_gen

      e_num_genabs=e_num_gen*sign(1._rk,e_den_gen)
      e_num_genabs_del=e_num_genabs-e_num_genabs_ave
      e_num_genabs_ave=e_num_genabs_ave+e_num_genabs_del/passes
      e_num_genabs_vn1=e_num_genabs_vn1+(e_num_genabs-e_num_genabs_ave)*e_num_genabs_del

      e_den_genabs=abs(e_den_gen)
      e_den_genabs_del=e_den_genabs-e_den_genabs_ave
      e_den_genabs_ave=e_den_genabs_ave+e_den_genabs_del/passes
      e_den_genabs_vn1=e_den_genabs_vn1+(e_den_genabs-e_den_genabs_ave)*e_den_genabs_del

      e_num_genabs_e_den_genabs_vn1=e_num_genabs_e_den_genabs_vn1+(e_num_genabs-e_num_genabs_ave)*e_den_genabs_del

      e_num_genabs_var=e_num_genabs_vn1/(passes-1)
      e_den_genabs_var=e_den_genabs_vn1/(passes-1)
      e_num_genabs_e_den_genabs_var=e_num_genabs_e_den_genabs_vn1/(passes-1)
      e_num_genabs_sigma=sqrt(e_num_genabs_var)
      e_den_genabs_sigma=sqrt(e_den_genabs_var)
      e_num_genabs_e_den_genabs_sigma=sqrt(e_num_genabs_e_den_genabs_var)
      e_genabs_ave=(e_num_genabs_ave/e_den_genabs_ave)
      if(passes.ge.2) then
        e_genabs_ave=e_genabs_ave/(1+(e_den_genabs_var/e_den_genabs_ave**2-(e_num_genabs_e_den_genabs_var/(e_num_genabs_ave*e_den_genabs_ave)))/passes)
        e_genabs_err=abs(e_num_genabs_ave/e_den_genabs_ave)*sqrt((e_num_genabs_var/e_num_genabs_ave**2 + e_den_genabs_var/e_den_genabs_ave**2 -2*e_num_genabs_e_den_genabs_var/(e_num_genabs_ave*e_den_genabs_ave))/passes)
      else
        e_genabs_err=0
      endif

      e_gen=e_num_gen/e_den_gen ! This and the next 31 lines are for calculating the non-integrated autocorrelation time
      e_gen_cum=e_gen_cum+e_gen
      e_gen2_cum=e_gen2_cum+e_gen**2
      e_gen_del=e_gen-e_gen_ave
      e_gen_ave=e_gen_ave+e_gen_del/passes
      e_gen_vn1=e_gen_vn1+(e_gen-e_gen_ave)*e_gen_del
      e_gen_var=e_gen_vn1/(passes-1)
!     write(6,'(/,''e_gen, e_gen_ave, e_gen_del, e_gen_vn1'',9f16.8)') e_gen, e_gen_ave, e_gen_del, e_gen_vn1

      if(passes.ge.2) then
        e_genpp_del=e_gen-e_genpp_ave
        e_genpp_ave=(passes*e_gen_ave-e_gen_gen1)/(passes-1)
        e_genp_e_genpp_cum=e_genp_e_genpp_cum+e_gen*e_genp
        e_genp_cum=e_gen_cum-e_gen
        e_genpp_cum=e_gen_cum-e_gen_gen1
        e_gen_av=e_gen_cum/((iblk-1)*passes)
        e_gen2_av=e_gen2_cum/((iblk-1)*passes)
        e_genp_e_genpp_av=e_genp_e_genpp_cum/((iblk-1)*passes-1)
        e_genp_av=e_genp_cum/((iblk-1)*passes-1)
        e_genpp_av=e_genpp_cum/((iblk-1)*passes-1)
        if(passes.ge.3) then
          e_genp_e_genpp_vn1=e_genp_e_genpp_vn1+e_genp_del*e_genpp_del
          e_genp_e_genpp_var=e_genp_e_genpp_vn1/(passes-2)
        endif
!       tau_corr_nonint=-1/log((e_genp_e_genpp_av-e_genp_av*e_genpp_av)/(e_gen2_av-e_gen_av**2))
        tau_corr_nonint=-1/log(e_genp_e_genpp_var/e_gen_var)
!       tau_corr_nonint=min(max(0._rk,tau_corr_nonint),1.d300)
        t_corr_nonint=1+2*tau_corr_nonint
!***b
        if((e_genp_e_genpp_var.le.0._rk .or. e_gen_var.le.0._rk .or. e_genp_e_genpp_var.gt.e_gen_var) .and. master_core) write(6,'(''e_genp_e_genpp_var, e_gen_var='',9d12.4)') e_genp_e_genpp_var, e_gen_var
!        if(e_genp_e_genpp_var.le.0._rk .or. e_gen_var.le.0._rk .or. e_genp_e_genpp_var.gt.e_gen_var) write(6,'(''e_genp_e_genpp_var, e_gen_var='',9d12.4)') e_genp_e_genpp_var, e_gen_var
!***e
      else
        e_gen_gen1=e_gen
      endif
      e_genp_del=e_gen_del

! Calculate e_est in preparation for calculating e_trial and reweight_factor for population control
      if(growth_estim) then
          p_eigval=w_abs_gen/(reweight_factor_inv*w_abs_gen_prev)
          e_est_gen=e_trial+(1-p_eigval)/tau
          e_est=e_est+(e_est_gen-e_est)/passes
      else
        if(e_den_blkabs_cum+abs(e_den_blk).ne.0) then
          e_est=(e_num_blkabs_cum+e_num_blk*sign(1._rk,e_den_blk))/(e_den_blkabs_cum+abs(e_den_blk))
        endif
      endif

! After equilibration blocks are done, use e_trial and e_est for calculating reweight_factor but before that use just the current and the target populations.
! e_trial is updated during the equilibration runs, but then held fixed.
! It is used for the projector and its reweight_factor after one block of the main run is done, but not for the reweight factor before that.
!     reweight_factor_inv_prev=reweight_factor_inv
      if(iblkk.le.ntimes_nblk_eq*nblk_eq) then
        e_trial=e_trial+sign(min(abs(e_est-e_trial),1._rk),e_est-e_trial) ! Warning: when we remove population control error, we should not update e_trial after first equilibration is done

        reweight_factor_inv = min(2._rk, max(0.5_rk,                             (real(w_abs_gen_target,rk)/w_abs_gen)**min(1._rk,tau*population_control_exponent)))
      else
        reweight_factor_inv = min(2._rk, max(0.5_rk, (1/(1+tau*(e_trial-e_est)))*(real(w_abs_gen_target,rk)/w_abs_gen)**min(1._rk,tau*population_control_exponent)))
      endif
      reweight_factor_inv = min(reweight_factor_inv, reweight_factor_inv_max)

!     if(ipr.ge.1 .or. (ipr.ge.0.and.iblkk.le.ntimes_nblk_eq*nblk_eq)) write(6,'(''iblkk, istep, w_perm_initiator_gen, nwalk, my_w_abs_gen, w_abs_gen_imp, e_gen, e_est, e_trial, reweight_factor_inv, t_c='',2i6,f8.1,i9,f11.1,f10.1,f10.4,2f11.5,f10.5,f8.1)') &
!&    iblkk, istep, w_perm_initiator_gen, nwalk, my_w_abs_gen, w_abs_gen_imp, e_num_gen/e_den_gen, e_est, e_trial, reweight_factor_inv, t_corr_nonint
      if(ipr.ge.1 .or. (ipr.ge.0.and.iblkk.le.ntimes_nblk_eq*nblk_eq)) write(6,'(''iblkk, istep, w_perm_initiator_gen, nwalk, w_abs_gen, w_abs_gen_imp, e_gen, e_est, e_trial, reweight_factor_inv, t_c='',2i6,f8.1,i9,f11.1,f10.1,f10.4,2f11.5,f10.5,f8.1)') &
 &    iblkk, istep, w_perm_initiator_gen, disp_nwalk, w_abs_gen, w_abs_gen_imp, e_num_gen/e_den_gen, e_est, e_trial, reweight_factor_inv, t_corr_nonint
      call flush(6)

      !if(.not.(iblkk==1 .and. istep ==1)) call stochastic_reconfiguration2(my_nwalk,walk_wt,e_num_walker,e_den_walker)
      if (run_type .eq. 'SR') call stochastic_reconfiguration2(my_nwalk,walk_wt,e_num_walker,e_den_walker)

! Set tau and r_initiator back to original values when w_abs_gen_target is reached for the 1st time.
      if(reached_w_abs_gen.eq.0 .and. w_abs_gen.ge.w_abs_gen_target) then
        reached_w_abs_gen=2
        tau_ratio=tau_sav/tau
        tau=tau_sav
        r_initiator=r_initiator_sav
        if (allocated(minus_tau_H_values)) then                     ! Added by HJC if minus tau H not used!!!!
          minus_tau_H_values=minus_tau_H_values*tau_ratio
        endif
        write(6,'(/,''w_abs_gen_target='',i9,'' reached at iblkk, istep='',2i6,'', tau, r_initiator reset to actual tau='',2f10.6)') w_abs_gen_target, iblkk, istep, tau, r_initiator
        call my_second(2,'w_abs_gen_target reached') ; write(6,*)
      endif

      tau_previous=tau

      if(ipr.ge.2) write(6,'(''reweight_factor='',f10.5)') 1._rk/reweight_factor_inv
      if(ipr.ge.2) write(6,'(''nwalk, nwalk_before_merge, w_abs_gen'', i6,i7,f10.1)') nwalk, nwalk_before_merge, w_abs_gen
!     if(ipr.ge.0) write(1,'(i10,i6,2i8,2f15.3)') nstep*(iblkk-1)+istep, nwalk, nint(w_gen), w_abs_gen, e_den_gen, 1._rk/reweight_factor_inv
      if(ipr.ge.0) write(1,'(i10,f12.6,es13.6,f19.12,i9)') nstep*(iblkk-1)+istep, 1/reweight_factor_inv, w_abs_gen, e_num_gen/e_den_gen, nwalk
      call flush(6)
      if(ipr.ge.1) call flush(1)

    enddo ! istep

!*** edit by AR [28/11/13] : update block quantities
    if(ncores>1) then
      tmp_real_array(1)=my_w2_cum
      tmp_real_array(2)=my_e_den2_cum
      tmp_real_array(3)=my_e_den_abs_cum
      tmp_real_array(4)=my_e_num2_cum
      tmp_real_array(5)=my_e_num_abs_cum
      tmp_real_array(6)=my_e_num_e_den_cum
      tmp_real_array(7)=my_w_abs_before_merge_cum
      tmp_real_array(8)=my_nwalk_before_merge_cum
      tmp_real_array(9)=my_nwalk_cum
      tmp_real_array(10)=my_w_perm_initiator_blk
      tmp_real_array(11)=my_w_perm_initiator_cum
      tmp_real_array(12)=my_w_abs_blk_imp
     !tmp_real_array(13)=nwalk
      call mpi_allred(tmp_real_array(1:12))
      w2_cum=w2_cum+tmp_real_array(1)
      e_den2_cum = e_den2_cum + tmp_real_array(2)
      e_den_abs_cum = e_den_abs_cum + tmp_real_array(3)
      e_num2_cum = e_num2_cum + tmp_real_array(4)
      e_num_abs_cum = e_num_abs_cum + tmp_real_array(5)
      e_num_e_den_cum = e_num_e_den_cum + tmp_real_array(6)
      w_abs_before_merge_cum=w_abs_before_merge_cum+tmp_real_array(7)
      nwalk_before_merge_cum=nwalk_before_merge_cum+int(tmp_real_array(8),i8b)
      nwalk_cum=nwalk_cum+int(tmp_real_array(9),i8b)
      w_perm_initiator_blk=w_perm_initiator_blk+tmp_real_array(10)
      w_perm_initiator_cum=w_perm_initiator_cum+tmp_real_array(11)
      w_abs_blk_imp=w_abs_blk_imp+tmp_real_array(12)
     !disp_nwalk=nint(tmp_real_array(13))
    else
      e_den2_cum = e_den2_cum + my_e_den2_cum
      e_den_abs_cum = e_den_abs_cum + my_e_den_abs_cum
      e_num2_cum = e_num2_cum + my_e_num2_cum
      e_num_abs_cum = e_num_abs_cum + my_e_num_abs_cum
      e_num_e_den_cum = e_num_e_den_cum + my_e_num_e_den_cum
      w_abs_before_merge_cum=w_abs_before_merge_cum+my_w_abs_before_merge_cum
      nwalk_before_merge_cum=nwalk_before_merge_cum+my_nwalk_before_merge_cum
      nwalk_cum=nwalk_cum+my_nwalk_cum
      w_perm_initiator_blk=w_perm_initiator_blk+my_w_perm_initiator_blk
      w_perm_initiator_cum=w_perm_initiator_cum+my_w_perm_initiator_cum
      w_abs_blk_imp=w_abs_blk_imp+my_w_abs_blk_imp
    endif
!*** end of edit
!    return

    w_cum=w_cum+w_blk ; e_num_cum=e_num_cum+e_num_blk ; e_den_cum=e_den_cum+e_den_blk
    w_blk2_cum=w_blk2_cum+w_blk**2
    e_den_blk2_cum=e_den_blk2_cum+e_den_blk**2
    e_num_blk2_cum=e_num_blk2_cum+e_num_blk**2
    w_blkabs_cum=w_blkabs_cum+abs(w_blk)
    e_den_blkabs_cum=e_den_blkabs_cum+abs(e_den_blk)
    e_num_blkabs_cum=e_num_blkabs_cum+e_num_blk*sign(1._rk,e_den_blk)
    e_num_blk_e_den_blk_cum=e_num_blk_e_den_blk_cum+e_num_blk*e_den_blk

    e_den_av=e_den_cum/w_abs_cum
    e_den_abs_av=e_den_abs_cum/w_abs_cum
    e_den_genabs_av=e_den_genabs_cum/(nstep*iblk)
    e_den_blkabs_av=e_den_blkabs_cum/iblk
    e_num_av=e_num_cum/w_abs_cum
    e_num_abs_av=e_num_abs_cum/w_abs_cum
    e_num_genabs_av=e_num_genabs_cum/(nstep*iblk)
    e_num_blkabs_av=e_num_blkabs_cum/iblk

    e_den_sigma=sqrt(e_den2_cum/w_abs_cum -e_den_av**2)
    e_den_abs_sigma=sqrt(e_den2_cum/w_abs_cum -e_den_abs_av**2)
    e_den_genabs_sigma=sqrt(e_den_gen2_cum/(nstep*iblk)-e_den_genabs_av**2)
    e_den_blkabs_sigma=sqrt(e_den_blk2_cum/iblk-e_den_blkabs_av**2)
    e_num_sigma=sqrt(e_num2_cum/w_abs_cum -e_num_av**2)
    e_num_abs_sigma=sqrt(e_num2_cum/w_abs_cum -e_num_abs_av**2)
    e_num_genabs_sigma=sqrt(e_num_gen2_cum/(nstep*iblk)-e_num_genabs_av**2)
    e_num_blkabs_sigma=sqrt(e_num_blk2_cum/iblk-e_num_blkabs_av**2)
    e_num_e_den_var=e_num_e_den_cum/w_abs_cum-e_num_abs_av*e_den_abs_av
    e_num_gen_e_den_gen_var=e_num_gen_e_den_gen_cum/(nstep*iblk)-e_num_genabs_av*e_den_genabs_av
    if(iblk.ge.2) e_num_blk_e_den_blk_var=(e_num_blk_e_den_blk_cum/iblk-e_num_blkabs_av*e_den_blkabs_av)*iblk/(iblk-1)

    e_num_e_den_av=e_num_e_den_cum/w_abs_cum
    e_num_gen_e_den_gen_av=e_num_gen_e_den_gen_cum/(nstep*iblk)
    e_num_blk_e_den_blk_av=e_num_blk_e_den_blk_cum/iblk

    e_blk=e_num_blk/e_den_blk

    e_num_blkabs=e_num_blk*sign(1._rk,e_den_blk)
    e_num_blkabs_del=e_num_blkabs-e_num_blkabs_ave
    e_num_blkabs_ave=e_num_blkabs_ave+e_num_blkabs_del/iblk
    e_num_blkabs_vn1=e_num_blkabs_vn1+(e_num_blkabs-e_num_blkabs_ave)*e_num_blkabs_del

    e_den_blkabs=abs(e_den_blk)
    e_den_blkabs_del=e_den_blkabs-e_den_blkabs_ave
    e_den_blkabs_ave=e_den_blkabs_ave+e_den_blkabs_del/iblk
    e_den_blkabs_vn1=e_den_blkabs_vn1+(e_den_blkabs-e_den_blkabs_ave)*e_den_blkabs_del

    e_num_blkabs_e_den_blkabs_vn1=e_num_blkabs_e_den_blkabs_vn1+(e_num_blkabs-e_num_blkabs_ave)*e_den_blkabs_del

    if(iblk.ge.2) then
      e_num_blkabs_var=e_num_blkabs_vn1/(iblk-1)
      e_den_blkabs_var=e_den_blkabs_vn1/(iblk-1)
      e_num_blkabs_e_den_blkabs_var=e_num_blkabs_e_den_blkabs_vn1/(iblk-1)
    endif
    e_num_blkabs_sigma=sqrt(e_num_blkabs_var)
    e_den_blkabs_sigma=sqrt(e_den_blkabs_var)
    e_num_blkabs_e_den_blkabs_sigma=sqrt(e_num_blkabs_e_den_blkabs_var)
    e_blkabs_ave=(e_num_blkabs_ave/e_den_blkabs_ave)
    if(iblk.ge.2) then
      e_blkabs_ave=e_blkabs_ave/(1+(e_den_blkabs_var/e_den_blkabs_ave**2-(e_num_blkabs_e_den_blkabs_var/(e_num_blkabs_ave*e_den_blkabs_ave)))/iblk) ! Remove 1st-order bias in ratio of averages
      e_blkabs_err=abs(e_num_blkabs_ave/e_den_blkabs_ave)*sqrt((e_num_blkabs_var/e_num_blkabs_ave**2 + e_den_blkabs_var/e_den_blkabs_ave**2 -2*e_num_blkabs_e_den_blkabs_var/(e_num_blkabs_ave*e_den_blkabs_ave))/iblk)
      t_corr=(e_blkabs_err/e_genabs_err)**2
    else
      e_blkabs_err=0
      t_corr=0
    endif
    if(t_corr.gt.1000._rk) then
      write(6,'(''e_blkabs_err/e_genabs_err'',9d12.4)') e_blkabs_err/e_genabs_err
      write(6,'(''e_blkabs_err, e_genabs_err, e_num_blkabs_var, e_den_blkabs_var, e_num_blkabs_e_den_blkabs_var='',9d12.4)') e_blkabs_err, e_genabs_err, e_num_blkabs_var, e_den_blkabs_var, e_num_blkabs_e_den_blkabs_var
      write(6,'(''abs(e_num_blkabs_ave/e_den_blkabs_ave), e_num_blkabs_var/e_num_blkabs_ave**2, e_den_blkabs_var/e_den_blkabs_ave**2, e_num_blkabs_e_den_blkabs_var/(e_num_blkabs_ave*e_den_blkabs_ave)'',9es12.4)') abs(e_num_blkabs_ave/e_den_blkabs_ave), e_num_blkabs_var/e_num_blkabs_ave**2, e_den_blkabs_var/e_den_blkabs_ave**2, e_num_blkabs_e_den_blkabs_var/(e_num_blkabs_ave*e_den_blkabs_ave)
      write(6,'(''abs(e_num_genabs_ave/e_den_genabs_ave), e_num_genabs_var/e_num_genabs_ave**2, e_den_genabs_var/e_den_genabs_ave**2, e_num_genabs_e_den_genabs_var/(e_num_genabs_ave*e_den_genabs_ave)'',9es12.4)') abs(e_num_genabs_ave/e_den_genabs_ave), e_num_genabs_var/e_num_genabs_ave**2, e_den_genabs_var/e_den_genabs_ave**2, e_num_genabs_e_den_genabs_var/(e_num_genabs_ave*e_den_genabs_ave)
    endif

    e_av=e_num_cum/e_den_cum
    e_abs_av=e_num_abs_cum/e_den_abs_cum
    e_genabs_av=e_num_genabs_cum/e_den_genabs_cum
    e_blkabs_av=e_num_blkabs_cum/e_den_blkabs_cum
!   corrected=e_blkabs_av/(1+((e_den_blkabs_sigma/e_den_blkabs_av)**2 - e_num_blk_e_den_blk_var/(e_num_blkabs_av*e_den_blkabs_av))/iblk)
    if(iblk.ge.2) e_blkabs_av=e_blkabs_av/(1+((e_den_blkabs_sigma/e_den_blkabs_av)**2 - e_num_blkabs_e_den_blkabs_var/(e_num_blkabs_av*e_den_blkabs_av))/iblk)  ! Remove 1st-order bias in ratio of averages
!   write(6,'(''e_num_blkabs_av, e_num_blkabs_ave, e_num_blkabs_av-e_num_blkabs_ave, e_den_blkabs_av, e_den_blkabs_ave, e_den_blkabs_av-e_den_blkabs_ave, &
!  & e_num_blk_e_den_blk_var, e_num_blkabs_e_den_blkabs_var, e_num_blk_e_den_blk_var-e_num_blkabs_e_den_blkabs_var'',6f16.6,9f18.3)') e_num_blkabs_av, e_num_blkabs_ave, e_num_blkabs_av-e_num_blkabs_ave, e_den_blkabs_av, e_den_blkabs_ave, e_den_blkabs_av-e_den_blkabs_ave, &
!  & e_num_blk_e_den_blk_var, e_num_blkabs_e_den_blkabs_var, e_num_blk_e_den_blk_var-e_num_blkabs_e_den_blkabs_var

    e_err=abs(e_av)*sqrt(((e_num_sigma/e_num_av)**2+(e_den_sigma/e_den_av)**2-2*e_num_e_den_var/(e_num_av*e_den_av))/w_abs_cum)
    e_abs_err=abs(e_abs_av)*sqrt(((e_num_abs_sigma/e_num_abs_av)**2 + (e_den_abs_sigma/e_den_abs_av)**2 &
 &  -2*e_num_e_den_var/(e_num_abs_av*e_den_abs_av))/w_abs_cum)

!   e_genabs_err=abs(e_genabs_av)*sqrt(((e_num_genabs_sigma/e_num_genabs_av)**2 + (e_den_genabs_sigma/e_den_genabs_av)**2 &
!&  -2*e_num_gen_e_den_gen_var/(e_num_genabs_av*e_den_genabs_av))/(nstep*iblk))
!   e_blkabs_err=abs(e_blkabs_av)*sqrt(((e_num_blkabs_sigma/e_num_blkabs_av)**2 + (e_den_blkabs_sigma/e_den_blkabs_av)**2 &
!&  -2*e_num_blk_e_den_blk_var/(e_num_blkabs_av*e_den_blkabs_av))/iblk)

    n_eff=real(w_cum,rk)**2/real(w2_cum,rk)
    n_abs_eff=real(w_abs_cum,rk)**2/real(w2_cum,rk)
    n_genabs_eff=real(w_genabs_cum,rk)**2/real(w_gen2_cum,rk)
    n_blkabs_eff=real(w_blkabs_cum,rk)**2/real(w_blk2_cum,rk)
    e_den_eff=e_den_cum**2/e_den2_cum
    e_den_abs_eff=e_den_abs_cum**2/e_den2_cum
    e_den_genabs_eff=e_den_genabs_cum**2/e_den_gen2_cum
    e_den_blkabs_eff=e_den_blkabs_cum**2/e_den_blk2_cum

!   if(ipr.ge.1) call write_debug

! This takes into account that if nstep is not >> t_corr then t_corr may be greatly underestimated.  In that case, t_corr_nonint may be a better estimate than t_corr.
    !write(6,'(''t_corr_nonint,t_corr='',9es12.4)') t_corr_nonint,t_corr
    e_blkabs_corrected_err=e_blkabs_err*sqrt(max(1._rk,t_corr_nonint/t_corr))


!   write(6,'(''e_blkabs_av, corrected, correction'',2(f11.6,''(''i6,'')''), 9d12.4)')  e_blkabs_av, nint(1000000*e_blkabs_err), corrected, nint(1000000*e_blkabs_err), corrected-e_blkabs_av, e_den_blkabs_sigma, e_num_blk_e_den_blk_var

!   write(6,'(''iblk, w_perm_initiator, nwalk, w_abs, w_abs_imp='',i6,f9.1,i8,i9,i9,'' e_blk='',f10.4,'' e='',4(f12.6,''(''i6,'')''), '' e_trial, rew_fac_inv, t_c='',2f11.5,f8.1)') iblk, w_perm_initiator_blk/nstep, nwalk, nint(w_abs_blk/nstep), nint(w_abs_blk_imp/nstep), &
!&  e_num_blk/e_den_blk, e_abs_av, nint(1000000*e_abs_err), e_genabs_av, nint(1000000*e_genabs_err), e_blkabs_av, nint(1000000*e_blkabs_err), e_av, nint(1000000*e_err), e_trial, reweight_factor_inv, t_corr_nonint
!   write(6,'(''iblk, w_perm_initiator, nwalk, w_abs, w_abs_imp='',i6,f9.1,i8,i9,i9,'' e_blk='',f10.4,'' e='',4(f14.8,''(''i8,'')''), '' e_trial, rew_fac_inv, t_c='',2f11.5,f8.1)') iblk, w_perm_initiator_blk/nstep, nwalk, nint(w_abs_blk/nstep), nint(w_abs_blk_imp/nstep), &
!&  e_num_blk/e_den_blk, e_abs_av, nint(100000000*e_abs_err), e_genabs_av, nint(100000000*e_genabs_err), e_blkabs_ave, nint(100000000*e_blkabs_err), e_av, nint(100000000*e_err), e_trial, reweight_factor_inv, t_corr_nonint

    call mpi_barr()
    write(6,'(''iblk, w_perm_initiator, nwalk, w_abs, w_abs_imp='',i6,f9.1,2i9,i8,'' e_blk='',f10.4,'' e='',3(f14.8,''(''i8,'')''), '' e_trial, rew_fac_inv, t_c='',2f11.5,f8.1)') iblk, w_perm_initiator_blk/nstep, disp_nwalk, nint(w_abs_blk/nstep), nint(w_abs_blk_imp/nstep), &
 &  e_num_blk/e_den_blk, e_genabs_av, nint(100000000*e_genabs_err), e_blkabs_ave, nint(100000000*e_blkabs_err), e_blkabs_ave, nint(100000000*e_blkabs_corrected_err), e_trial, reweight_factor_inv, t_corr_nonint
    !call my_second (2,'block')
    call flush(6)

!   if(ipr.ge.1 .and. (hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard')) then
!     write(6,'(''dominant_det_most, dominant_det_abs_wt_tmp='',i10,f10.1)') dominant_det_most, dominant_det_abs_wt_tmp
!     write(6,'(''dominant_det_list  = '',100i10)')   (dominant_det_list(iwalk),iwalk=1,n_dominant_det)
!     write(6,'(''dominant_det_wt    = '',100f10.1)') (dominant_det_wt(iwalk),iwalk=1,n_dominant_det)
!     write(6,'(''dominant_abs_det_wt= '',100f10.1)') (dominant_abs_det_wt(iwalk),iwalk=1,n_dominant_det)
!   endif
!   call flush(6)

! If the run is not equilibrated, increase ntimes_nblk_eq
! The equilibration is done first at a larger tau and then at the actual tau to be used.
! Equilibration is judged by 4 criteria:
! 1) First w_abs_gen should have exceeded w_abs_gen_target
! Then
! 2) w_abs_blk has to stop growing and become smaller and then larger again
! 3) w_perm_initiator_blk has to stop growing and become smaller and then larger again
! 4) e_blk has to stop dropping and rise and drop again
! Finally for small systems where the energies and the populations are nearly constant, there is a special criterion for judging that everything is equilibrated
! Warning tmp changed lt and gt to le and ge so that it works when wavefn. is very good and there are no oscillations:
!   if(reached_w_abs_gen.eq.2 .and. equilibrated_w_abs_blk.eq.0 .and. w_abs_blk.lt.w_abs_blk_prev) equilibrated_w_abs_blk=1
!   if(equilibrated_w_abs_blk.eq.1 .and. w_abs_blk.gt.w_abs_blk_prev) then
    if(reached_w_abs_gen.eq.2 .and. equilibrated_w_abs_blk.eq.0 .and. w_abs_blk.le.w_abs_blk_prev) equilibrated_w_abs_blk=1
    if(equilibrated_w_abs_blk.eq.1 .and. w_abs_blk.ge.w_abs_blk_prev) then
      equilibrated_w_abs_blk=2 ! w_abs_gen is considered equilibrated after it increases, decreases and increases again
      write(6,'(''Equilibration of population achieved at iblkk='',i6,'', w_perm_initiator, w_abs_imp, e_blk_ave, e_blk='',2f9.1,2f14.6)') iblkk, w_perm_initiator_blk/nstep, w_abs_blk_imp/nstep, e_blkabs_ave, e_blk
      call flush(6)
    endif

    if (run_type .ne. 'fixed_node1' ) then
! Warning tmp changed lt and gt to le and ge so that it works when wavefn. is very good and there are no oscillations:
!       if(reached_w_abs_gen.eq.2 .and. equilibrated_w_perm_initiator_blk.eq.0 .and. w_perm_initiator_blk.lt.w_perm_initiator_blk_prev) equilibrated_w_perm_initiator_blk=1
!       if(equilibrated_w_perm_initiator_blk.eq.1 .and. w_perm_initiator_blk.gt.w_perm_initiator_blk_prev) then
        if(reached_w_abs_gen.eq.2 .and. equilibrated_w_perm_initiator_blk.eq.0 .and. w_perm_initiator_blk.le.w_perm_initiator_blk_prev) equilibrated_w_perm_initiator_blk=1
        if(equilibrated_w_perm_initiator_blk.eq.1 .and. w_perm_initiator_blk.ge.w_perm_initiator_blk_prev) then
          equilibrated_w_perm_initiator_blk=2 ! w_perm_initiator_blk is considered equilibrated after it increases, decreases and increases again
          write(6,'(''Equilibration of perm_initiator population achieved at iblkk='',i6,'', w_perm_initiator, w_abs_imp, e_blk_ave, e_blk='',2f9.1,2f14.6)') iblkk, w_perm_initiator_blk/nstep, w_abs_blk_imp/nstep, e_blkabs_ave, e_blk
          call flush(6)
        endif
    else
      equilibrated_w_perm_initiator_blk=2 ! fixed_node1 is often used when there is no state with a large weight, so the perm_initiator concept makes no sense.  So, consider it always equilibrated.
    endif

! Warning tmp changed lt and gt to le and ge so that it works when wavefn. is very good and there are no oscillations:
!   if(reached_w_abs_gen.eq.2 .and. equilibrated_e_blk.eq.0 .and. e_blk.gt.e_blk_prev) equilibrated_e_blk=1
!   if(equilibrated_e_blk.eq.1 .and. e_blk.lt.e_blk_prev) then
    if(reached_w_abs_gen.eq.2 .and. equilibrated_e_blk.eq.0 .and. e_blk.ge.e_blk_prev) equilibrated_e_blk=1
    if(equilibrated_e_blk.eq.1 .and. e_blk.le.e_blk_prev) then
      equilibrated_e_blk=2 ! e_blk is considered equilibrated after it increases, decreases and increases again
      write(6,'(''Equilibration of energy achieved at iblkk='',i6,'', w_perm_initiator, w_abs_imp, e_blk_ave, e_blk='',2f9.1,2f14.6)') iblkk, w_perm_initiator_blk/nstep, w_abs_blk_imp/nstep, e_blkabs_ave, e_blk
      call flush(6)
    endif
    if(equilibrated_w_abs_blk.eq.2 .and. equilibrated_w_perm_initiator_blk.eq.2 .and. equilibrated_e_blk.eq.2 .and. equilibrated_all.le.1) then
      equilibrated_all=2
      write(6,'(''Equilibration of everything achieved at iblkk='',i6,'', w_perm_initiator, w_abs_imp, e_blk_ave, e_blk='',2f9.1,2f14.6)') iblkk, w_perm_initiator_blk/nstep, w_abs_blk_imp/nstep, e_blkabs_ave, e_blk
      call flush(6)
    endif
    if(equilibrated_all.le.1 .and. iblkk.eq.ntimes_nblk_eq*nblk_eq) ntimes_nblk_eq=ntimes_nblk_eq+1
    w_abs_blk_prev=w_abs_blk
    w_perm_initiator_blk_prev=w_perm_initiator_blk
    e_blk_prev=e_blk

! Print out distribution of distances from deterministic space after every nblk_eq equilibrium blocks and after the entire run.
! To save time, only the distances at the end of each block are averaged
    w_abs_genlast_cum=w_abs_genlast_cum+sum(abs(walk_wt(1:nwalk)))
    do iwalk=1,nwalk
      imp_distance_fraction(min(imp_distance(iwalk),max_imp_distance))=imp_distance_fraction(min(imp_distance(iwalk),max_imp_distance))+abs(walk_wt(iwalk))
     !imp_distance_fraction(max(0,min(imp_distance(iwalk),max_imp_distance)))=imp_distance_fraction(max(0,min(imp_distance(iwalk),max_imp_distance)))+abs(walk_wt(iwalk))  ! This includes -2 in the count for 0
    enddo

    if((mod(iblkk,nblk_eq).eq.0 .and. iblkk.le.ntimes_nblk_eq*nblk_eq) .or. iblkk.eq.ntimes_nblk_eq*nblk_eq+nblk) then
      write(6,'(/,''imp_distance         =           0          -2'',15i12  )') (i,i=1,max_imp_distance)
      write(6,'(  ''imp_distance_fraction='',17f12.8)') imp_distance_fraction(0)/w_abs_genlast_cum, imp_distance_fraction(-2)/w_abs_genlast_cum, imp_distance_fraction(1:max_imp_distance)/w_abs_genlast_cum
      call flush(6)
    endif

    !Reduce the max_wt_ratio and min_wt_ratio to find the max and min over all processes
    !Note that if one is using symmetry then the excitation level is not correct, because representatives can differ by more than 2 excitations.
    !call MPI_REDUCE(max_wt_ratio,all_max_wt_ratio,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    !call MPI_REDUCE(min_wt_ratio,all_min_wt_ratio,1,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    call mpi_red_max(max_wt_ratio,all_max_wt_ratio,0)
    call mpi_red_min(min_wt_ratio,all_min_wt_ratio,0)
    call mpi_gath(max_wt_ratio_excit_level,1,all_max_wt_ratio_excit_level,1,0)
    if(ipr.ge.1) write(6,'(''Max Ratio, Min Ratio, tau_multiplier, tau for spawning max wt 1='',9es12.4)') all_max_wt_ratio, all_min_wt_ratio, (diagonal_ham_highest-diagonal_ham_lowest)/all_max_wt_ratio, 1/all_max_wt_ratio
    if(ipr.ge.1) write(6,'(''all_max_wt_ratio_excit_level'',1000i3)') all_max_wt_ratio_excit_level(1:ncores)

    iblkk=iblkk+1

  enddo ! iblk

  write (6,*) "ndet=",ndet
  call flush(6)

  if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
  if(ipr.ge.1) write(6,'(''wavefn ='',10000es9.1)') (wavefn(idet),idet=1,ndet)
  if(ipr.ge.1) write(6,'(''psi_g ='',10000es9.1)') (psi_g(idet),idet=1,ndet)
  ! Divide by psi_g to get projected wavefn.
    if(importance_sampling.eq.1) then
      do idet=1,ndet
        if(psi_g(idet).ne.0._rk) then
          wavefn(idet)=wavefn(idet)/psi_g(idet)
        else
          wavefn(idet)=0
        endif
      enddo
    endif

  ! Calculate norm of averaged wavefn
    norm_of_wavefn=0
    do idet=1,ndet
      norm_of_wavefn=norm_of_wavefn+wavefn(idet)**2
    enddo
    norm_of_wavefn=sqrt(abs(norm_of_wavefn))
  endif

  if(ipr.lt.1 .and. (hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') .and. diagonalize_ham.eq.1)  then
!   write(6,'(''dominant_det_most, dominant_det_abs_wt_tmp='',i10,f10.1)') dominant_det_most, dominant_det_abs_wt_tmp
!   write(6,'(''dominant_det_list  = '',100i10)')   (dominant_det_list(iwalk),iwalk=1,n_dominant_det)
!   write(6,'(''dominant_det_wt    = '',100f10.1)') (dominant_det_wt(iwalk),iwalk=1,n_dominant_det)
!   write(6,'(''dominant_abs_det_wt= '',100f10.1)') (dominant_abs_det_wt(iwalk),iwalk=1,n_dominant_det)

! Calculate overlap of averaged wavefn with exact one
    if(diagonalize_ham.eq.1)  then
      overlap=0
      do idet=1,ndet
        overlap=overlap+wavefn_exact(idet)*wavefn(idet)
      enddo
      if(ipr.ge.2) write(6,'(''wavefn_exact='',10000i5)') nint(100*sqrt(real(ndet,rk))*wavefn_exact)
    endif
    if(ipr.ge.2) write(6,'(''wavefn_proj ='',10000i5)') nint(100*sqrt(real(ndet,rk))*(wavefn/norm_of_wavefn))
    if(overlap.lt.0._rk) then
      write(6,'(''Changing sign of projected wavefn.'')')
      overlap=-overlap
      wavefn=-wavefn
    endif
    overlap=sqrt(overlap/norm_of_wavefn)
    write(6,'(/,''Overlap of exact and projected wavefns='',f9.6)') overlap
  endif

  if(nwalk.le.1000) then
    write(6,'(/,''Walkers for last generation'')')
    if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
      write(6,'(''walk_det='',10000i5)')   (walk_det(iwalk),iwalk=1,nwalk)
      write(6,'(''walk_wt= '',10000f5.1)') (walk_wt(iwalk),iwalk=1,nwalk)
    else
      write(6,'(''walk_dets_up/dn='',10000(2i5,2x))') (walk_dets_up(iwalk),walk_dets_dn(iwalk),iwalk=1,nwalk)
      write(6,'(''walk_wt      ='',10000f12.1)') (walk_wt(iwalk),iwalk=1,nwalk)
    endif
  endif
! write(6,'(''nwalk, nwalk_before_merge, w_abs_gen'',i6,9i9)') nwalk, nwalk_before_merge, w_abs_gen

  if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
    write(6,'(/,''Comparison of exact wavefn. with wavefn. averaged over run:)'')')
    write(6,'(''idet        ='',10000i5)') (idet,idet=1,ndet)
    if(diagonalize_ham.eq.1) write(6,'(''wavefn_exact='',10000i5)') nint(100*sqrt(real(ndet,rk))*wavefn_exact)
    write(6,'(''wavefn_proj ='',10000i5)') nint(100*sqrt(real(ndet,rk))*(wavefn/norm_of_wavefn))
    write(6,'(''w_cum, overlap, e_blkabs_ave, energy_exact='',9f12.6)') w_cum, overlap, e_blkabs_ave, energy_exact
! Sort elements of lowest eigvec in order of descending absolute value.  Note this destroys original order.
    write(6,'(''diagonalize_ham='',i5)') diagonalize_ham
    if(diagonalize_ham.eq.1) then
      allocate(index(ndet),stat=istat)
      if(istat.ne.0) stop 'failed to allocate index'
      do idet=1,ndet
        index(idet)=idet
      enddo
      call shell_sort_eigvec(wavefn_exact(1),index,ndet)
      write(6,'(/,''Comparison of components in order of descending absolute value:'')')
      write(6,'(''index       ='',10000i5)') index
      write(6,'(''wavefn_exact='',10000i5)') nint(100*sqrt(real(ndet,rk))*wavefn_exact)
      write(6,'(''wavefn_proj ='',10000i5)') (nint(100*sqrt(real(ndet,rk))*(wavefn(index(idet))/norm_of_wavefn)),idet=1,ndet)
      write(6,'(''wavefn_proj ='',10000es9.1)') ((100*sqrt(real(ndet,rk))*(wavefn(index(idet))/norm_of_wavefn)),idet=1,ndet)
      deallocate(index)
    endif
  endif

  passes=real(nstep,rk)*real(nblk,rk)
  write(6,'(''e_genabs_ave, e_genabs_err, e_blkabs_ave, e_blkabs_err='',9es12.4)') e_genabs_ave, e_genabs_err, e_blkabs_ave, e_blkabs_err
  t_corr=(e_blkabs_err/e_genabs_err)**2
  if(n_warning.ne.0) write(6,'(''n_warning='',i10)') n_warning
  write(6,'(/,''During main run, n_sign_flips, n_sign_flips per unit time='',i6,f10.6)') &
 &n_sign_flips, real(n_sign_flips,rk)/(real(nstep,rk)*nblk*tau)
  if(ipr.ge.0) write(1,'(2i6,i9,f12.6,f10.6'' nstep, nblk, w_abs_gen_target, e_trial, tau'')') nstep, nblk, w_abs_gen_target, e_trial, tau
  write(6,'(''w_perm_initiator_av, nwalk_av, nwalk_before_merge_av, ratio, w_av, w_abs_av, w_abs_before_merge_av, ratio, e_den_genabs_av, e_den_blkabs_av/nstep='',f8.1,2f10.1,f6.2,3f10.1,f6.2,9f10.1)') &
 &w_perm_initiator_cum/passes, nwalk_cum/passes, nwalk_before_merge_cum/passes, real(nwalk_cum,rk)/nwalk_before_merge_cum, w_cum/passes, w_abs_cum/passes, w_abs_before_merge_cum/passes, real(w_abs_cum,rk)/w_abs_before_merge_cum, e_den_genabs_av, e_den_blkabs_av/nstep

! write(6,'(''Energy='',4(f11.6,''(''i6,'')''),'' energy_exact='',f10.5,'' t_corr_nonint, t_corr, nstep='',2f8.2,i6)') &
!&  e_abs_av, nint(1000000*e_abs_err), e_genabs_av, nint(1000000*e_genabs_err), e_blkabs_av, nint(1000000*e_blkabs_err), &
!&  e_av, nint(1000000*e_err), energy_exact, t_corr_nonint, t_corr, nstep
! write(6,'(''Energy='',4(f13.8,''(''i8,'')''),'' energy_exact='',f10.5,'' t_corr_nonint, t_corr, nstep='',2f8.2,i6)') &
!&  e_abs_av, nint(100000000*e_abs_err), e_genabs_av, nint(100000000*e_genabs_err), e_blkabs_ave, nint(100000000*e_blkabs_err), &
!&  e_av, nint(100000000*e_err), energy_exact, t_corr_nonint, t_corr, nstep
  write(6,'(''Energy='',3(f14.8,''(''i8,'')''),'' energy_exact='',f10.5,'' t_corr_nonint, t_corr, nstep='',2f8.2,i6)') &
 &  e_genabs_av, nint(100000000*e_genabs_err), e_blkabs_ave, nint(100000000*e_blkabs_err), e_blkabs_ave, nint(100000000*e_blkabs_corrected_err), &
 &  energy_exact, t_corr_nonint, t_corr, nstep


  if (off_diag_histograms) then
    write (6,*) "Histogram of abs weight multipliers spawned by HF:"
    write (6,*) "Bin,Lowerbound,Count"
    do i=1,nbins
      write (6,'(i5,f10.3,9i11)') i,lbounds_hf(i),bins_hf(i)
    enddo
    write (6,*) "Total=",sum(bins_hf)
    call flush(6)

    write (6,*) "Histogram of abs weight multipliers spawned by any states:"
    write (6,*) "Bin,Lowerbound,Singles,Doubles,Total,Singles%,Doubles%,Total%"
    do i=1,nbins
      write (6,'(i5,f10.3,3i11,3f10.6)') i,lbounds(i),bins_single(i),bins_double(i),bins(i),float(bins_single(i))/float(sum(bins_single)),float(bins_double(i))/float(sum(bins_double)),float(bins(i))/float(sum(bins))
    enddo
    write (6,'(a,9x,9i11)') "Total=",sum(bins_single),sum(bins_double),sum(bins)
    call flush(6)
  endif

  !Reduce the max_wt_ratio and min_wt_ratio to find the max and min over all processes
  !Note that if one is using symmetry then the excitation level is not correct, because representatives can differ by more than 2 excitations.
  !call MPI_REDUCE(max_wt_ratio,all_max_wt_ratio,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  !call MPI_REDUCE(min_wt_ratio,all_min_wt_ratio,1,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  call mpi_red_max(max_wt_ratio,all_max_wt_ratio,0)
  call mpi_red_min(min_wt_ratio,all_min_wt_ratio,0)
  call mpi_gath(max_wt_ratio_excit_level,1,all_max_wt_ratio_excit_level,1,0)
  if(ipr.ge.1) write(6,'(''Max Ratio, Min Ratio, tau_multiplier, tau for spawning max wt 1='',9es12.4)') all_max_wt_ratio, all_min_wt_ratio, (diagonal_ham_highest-diagonal_ham_lowest)/all_max_wt_ratio, 1/all_max_wt_ratio
  !write(6,'(''all_max_wt_ratio_excit_level'',1000i10)') all_max_wt_ratio_excit_level(1:ncores)
  if(ipr.ge.1) write(6,'(''all_max_wt_ratio_excit_level'',1000i3)') all_max_wt_ratio_excit_level(1:ncores)

  end subroutine walk
!-------------------------------------------------------------------------------------

!======AR 10/12/13: make the list of connected states local
subroutine reduce_connections(bpsi_t_connected_e_loc_num,bpsi_t_connected_e_loc_den,bpsi_t_connected_dets_up,bpsi_t_connected_dets_dn)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable,intent(inout) :: bpsi_t_connected_dets_up(:),bpsi_t_connected_dets_dn(:)
#else
  integer(ik),allocatable,intent(inout) :: bpsi_t_connected_dets_up(:),bpsi_t_connected_dets_dn(:)
#endif
  real(rk),allocatable,intent(inout) :: bpsi_t_connected_e_loc_num(:),bpsi_t_connected_e_loc_den(:)

  integer  :: iw
  real(rk) :: tmp_psi_t_connected_e_loc_num(size(bpsi_t_connected_e_loc_num)),tmp_psi_t_connected_e_loc_den(size(bpsi_t_connected_e_loc_num))
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec) :: tmp_psi_t_connected_dets_up(size(bpsi_t_connected_e_loc_num)),tmp_psi_t_connected_dets_dn(size(bpsi_t_connected_e_loc_num))
#else
  integer(ik) :: tmp_psi_t_connected_dets_up(size(bpsi_t_connected_e_loc_num)),tmp_psi_t_connected_dets_dn(size(bpsi_t_connected_e_loc_num))
#endif
  integer :: my_psitc,istat

  my_psitc = 0

  do iw=1,size(bpsi_t_connected_e_loc_num)
    if(whoami==get_det_owner(bpsi_t_connected_dets_up(iw),bpsi_t_connected_dets_dn(iw))) then
      my_psitc=my_psitc+1
      tmp_psi_t_connected_e_loc_num(my_psitc)=bpsi_t_connected_e_loc_num(iw)
      tmp_psi_t_connected_e_loc_den(my_psitc)=bpsi_t_connected_e_loc_den(iw)
      tmp_psi_t_connected_dets_up(my_psitc)=bpsi_t_connected_dets_up(iw)
      tmp_psi_t_connected_dets_dn(my_psitc)=bpsi_t_connected_dets_dn(iw)
    endif
  enddo

  deallocate(bpsi_t_connected_e_loc_num,bpsi_t_connected_e_loc_den,bpsi_t_connected_dets_up,bpsi_t_connected_dets_dn)
  allocate(bpsi_t_connected_e_loc_num(my_psitc),stat=istat)
  if(istat.ne.0) stop 'failed to reallocate psi_t_connected_e_loc_num'
  allocate(bpsi_t_connected_e_loc_den(my_psitc),stat=istat)
  if(istat.ne.0) stop 'failed to reallocate psi_t_connected_e_loc_den'
  allocate(bpsi_t_connected_dets_up(my_psitc),stat=istat)
  if(istat.ne.0) stop 'failed to reallocate psi_t_connected_dets_up'
  allocate(bpsi_t_connected_dets_dn(my_psitc),stat=istat)
  if(istat.ne.0) stop 'failed to reallocate psi_t_connected_dets_dn'

  bpsi_t_connected_e_loc_num(1:my_psitc)=tmp_psi_t_connected_e_loc_num(1:my_psitc)
  bpsi_t_connected_e_loc_den(1:my_psitc)=tmp_psi_t_connected_e_loc_den(1:my_psitc)
  bpsi_t_connected_dets_up(1:my_psitc)=tmp_psi_t_connected_dets_up(1:my_psitc)
  bpsi_t_connected_dets_dn(1:my_psitc)=tmp_psi_t_connected_dets_dn(1:my_psitc)
end subroutine reduce_connections
!-------------------------------------------------------------------------------------

!======AR 10/12/13: make the diag_elems list local
subroutine reduce_connections_diag(bpsi_t_connected_dets_up,bpsi_t_connected_dets_dn,bdiag_elems)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable,intent(in) :: bpsi_t_connected_dets_up(:),bpsi_t_connected_dets_dn(:)
#else
  integer(ik),allocatable,intent(in) :: bpsi_t_connected_dets_up(:),bpsi_t_connected_dets_dn(:)
#endif
  real(rk),allocatable,intent(inout) :: bdiag_elems(:)

  integer  :: iw
  real(rk) :: tmp_diag_elems(size(bpsi_t_connected_dets_up))

  integer :: my_psitc,istat

  my_psitc = 0

  do iw=1,size(bpsi_t_connected_dets_up)
    if(whoami==get_det_owner(bpsi_t_connected_dets_up(iw),bpsi_t_connected_dets_dn(iw))) then
      my_psitc=my_psitc+1
      tmp_diag_elems(my_psitc)=bdiag_elems(iw)
    endif
  enddo

  deallocate(bdiag_elems)
  allocate(bdiag_elems(my_psitc),stat=istat)
  if(istat.ne.0) stop 'failed to reallocate diag_elems'

  bdiag_elems(1:my_psitc)=tmp_diag_elems(1:my_psitc)
end subroutine reduce_connections_diag
!-------------------------------------------------------------------------------------

! using procedure pointers in order to avoid repeated string comparisons
!-- [TODO] need to make off_diagonal_moves standard!!!
!   subroutine init_off_diagonal(ham_type)
!     use hamiltonian_mod
!     implicit none
!     character*16,intent(in) :: ham_type

!     if(ham_type.eq.'heg') then
!         off_diagonal_move => off_diagonal_move_heg
!     elseif(ham_type.eq.'chem') then
!         off_diagonal_move => off_diagonal_move_chem
!     elseif(ham_type.eq.'hubbard2') then
!         stop "need to extend off_diagonal_move_hubbard for compatibilty"
!     elseif(ham_type.eq.'hubbardk') then
!         off_diagonal_move => off_diagonal_move_hubbard_k
! !         stop "need to extend off_diagonal_move_hubbard_k for compatibilty"
!     elseif(ham_type.eq.'hubbarddm') then
!         stop "need to extend off_diagonal_move_hubbard_dm for compatibilty"
!     endif
!   end subroutine init_off_diagonal
!-------------------------------------------------------------------------------------

  subroutine init_move
! define pointers for subroutine move
    implicit none

    ! MJO -- Init max and min ratio
    max_wt_ratio = -1.d99
    min_wt_ratio = 1.d99

    if(proposal_method.eq.'uniform') then
      move => move_uniform
    elseif(proposal_method.eq.'uniform2' .or. proposal_method.eq.'CauchySchwarz' .or. proposal_method.eq.'fast_heatbath') then
      if (use_exponential_projector) then
        move => move_uniform_exponential_projector
      elseif (run_type .eq. 'vmc') then
        move => move_uniform_vmc
      else
        move => move_uniform2
        write (6,*) "Proposal method ",proposal_method,"selected"
      endif
    elseif(proposal_method.eq.'heat_bath') then
        move => move_heat_bath
    elseif(proposal_method.eq.'heat_bath2') then
        move => move_heat_bath2
    elseif(proposal_method.eq.'heat_bath3') then
        move => move_heat_bath3
        !stop 'WARNING: need to move Hack by HJC [Aug 29 2012] directly in move_heat_bath3 '
    else
        stop 'proposal method must be uniform1/2 or CauchySchwarz or fast_heatbath or heat_bath1/2/3'
    endif
  end subroutine init_move
!-------------------------------------------------------------------------------------


  subroutine move_uniform(iwalk)
! ==========================================================================================================
! Description   : Do off-diagonal move choosing uniformly from n_connected_dets elements of the propagator
!               : for small systems (entire Hamiltonian is stored) (hamiltonian_type fictitious, read, hubbard)
! Author        : Cyrus Umrigar
! ----------------------------------------------------------------------------------------------------------
  use hamiltonian_mod
  implicit none
  integer, intent(in) :: iwalk
  integer proposal_prob_inv, iw_allowed_move, walk_det_new, i, iwalk_child
  integer :: nwalk_child
! real(rk) random, rannyu, walk_wt_child
  real(rk) rannyu, walk_wt_tmp, walk_wt_child

! Count the number of non-zero off-diagonal elements and set proposal_prob_inv
  n_connected_dets=0
  do i=1,ndet
    if(i.ne.walk_det(iwalk) .and. ham(i,walk_det(iwalk)).ne.0._rk) then
      n_connected_dets=n_connected_dets+1
      connected_dets(n_connected_dets)=i
    endif
  enddo
  proposal_prob_inv=n_connected_dets

  nwalk_child=max(nint(abs(walk_wt(iwalk))),1)
  walk_wt_child=walk_wt(iwalk)/nwalk_child
! do iwalk_child=1,abs(walk_wt(iwalk))
  do iwalk_child=1,nwalk_child

! Decide which of these states to move to
    iw_allowed_move=int(n_connected_dets*rannyu())+1

    walk_det_new=connected_dets(iw_allowed_move)

! Calculate wt of this determinant
    if(importance_sampling.eq.0) then
!     random=sign(rannyu(),ham(walk_det_new,walk_det(iwalk)))
!     walk_wt_tmp=-int(reweight_factor_inv*tau*proposal_prob_inv*ham(walk_det_new,walk_det(iwalk))+random) * sign(1,walk_wt(iwalk))
      !walk_wt_tmp=-(reweight_factor_inv*tau*proposal_prob_inv*ham(walk_det_new,walk_det(iwalk))) * walk_wt_child
      walk_wt_tmp=-(tau*proposal_prob_inv*ham(walk_det_new,walk_det(iwalk))) * walk_wt_child ! (HJC Sep 29 2012) reweight factor removed to be consistent with non-toy systems
    else
!     random=sign(rannyu(),ham(walk_det_new,walk_det(iwalk))*psi_g(walk_det_new)/psi_g(walk_det(iwalk)))
!     walk_wt_tmp=-int(reweight_factor_inv*tau*proposal_prob_inv*ham(walk_det_new,walk_det(iwalk))*psi_g(walk_det_new)/psi_g(walk_det(iwalk))+random) * sign(1,walk_wt(iwalk))
      !walk_wt_tmp=-(reweight_factor_inv*tau*proposal_prob_inv*ham(walk_det_new,walk_det(iwalk))*psi_g(walk_det_new)/psi_g(walk_det(iwalk))) * walk_wt_child
      walk_wt_tmp=-(tau*proposal_prob_inv*ham(walk_det_new,walk_det(iwalk))*psi_g(walk_det_new)/psi_g(walk_det(iwalk))) * walk_wt_child  ! HJC(Sep 29 2012) reweight factor removed to be consistent with non-toy systems
    endif
    if(ipr.ge.1) write(6,'(''iwalk, walk_det(iwalk), walk_det_new, walk_wt_child, walk_wt_tmp'', 3i5,2f10.6)') iwalk, walk_det(iwalk), walk_det_new, walk_wt_child, walk_wt_tmp
! write(6, '(''reweight_factor, tau, proposal_prob_inv, ham(walk_det_new, walk_det(iwalk)), random'', 9d12.4)') &
!& 1._rk/reweight_factor_inv, tau, proposal_prob_inv, ham(walk_det_new, walk_det(iwalk)), random
!   if(walk_wt_tmp.ne.0) then
      nwalk=nwalk+1
      snd_cnt(1)=snd_cnt(1)+1
      if(nwalk.gt.MWALK) then
        write(6,'(''stop move_uniform: nwalk>MWALK, e_est, e_trial, reweight_factor=''9d12.4)') e_est, e_trial, 1._rk/reweight_factor_inv
        call flush(6)
        stop 'nwalk>MWALK'
      endif
      walk_det(nwalk)=walk_det_new
      walk_wt(nwalk)=walk_wt_tmp
      if(initiator(iwalk).ge.2) then ! If it is spawned by an initiator then 1 otherwise 0
        initiator(nwalk)=1
      else
        initiator(nwalk)=0
      endif
!   endif

  enddo ! iwalk_child

! random=sign(rannyu(),walk_wt(iwalk)*(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk)))))
! walk_wt(iwalk)=int(reweight_factor_inv*walk_wt(iwalk)*(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk))))+random) ! Diagonal move
  !walk_wt(iwalk)=reweight_factor_inv*walk_wt(iwalk)*(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk)))) ! Diagonal move
  walk_wt(iwalk)=walk_wt(iwalk)*(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk)))) ! Diagonal move
  ! NO REWEIGHT FACTOR TO BE CONSISTENT WITH NON-TOY SYSTEMS
  if(ipr.ge.1) write(6,'(''Diagonal: iwalk, walk_det(iwalk), walk_wt(iwalk)='',2i5,f10.6)') iwalk, walk_det(iwalk), walk_wt(iwalk)

  end subroutine move_uniform
!-------------------------------------------------------------------------------------

  subroutine move_uniform2(iwalk)
! ==========================================================================================================
! Description   : Do off-diagonal move choosing uniformly from n_connected_dets elements of the propagator
!               : for large systems (entire Hamiltonian cannot be stored)
! Author        : Cyrus Umrigar
! Modified      : A Holmes, 17 Jan 2013. Spawn low-wt walkers with probability ~ abs(wt).
! ----------------------------------------------------------------------------------------------------------
  use hamiltonian_mod
  use semistoch
  use common_psi_t, only : hf_to_psit
  use hubbard, only : hamiltonian_hubbard, hamiltonian_hubbard_k, hamiltonian_hubbard_dm, hamiltonian_hubbard_k_space_sym, space_sym
  use chemistry, only : hamiltonian_chem, hamiltonian_chem_time_sym, time_sym,excitation_level
  use heg, only : hamiltonian_heg
  use tools, only : count_excitations
  use types, only : num_words, bits_per_word

  implicit none
  integer, intent(in) :: iwalk
  integer iwalk_child, nwalk_child, count_diag_neg
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec) det_j_up, det_j_dn
  type(ik_vec) tmp_dets_up(2),tmp_dets_dn(2)
#else
  integer(ik) det_j_up, det_j_dn
  integer(ik) tmp_dets_up(2),tmp_dets_dn(2)
#endif
! real(rk) random, rannyu, diagonal_matrix_element, walk_wt_child
  real(rk) :: rannyu
  real(rk) diagonal_matrix_element, diagonal_factor, walk_wt_child, weight_j
  logical :: spawn,use_wt
  character(len=7) :: fmt
  integer :: nnzero, excite_level(2), i
  data count_diag_neg/0/
  real(rk) :: weight_j_list(2)
  integer :: n_new_dets

  if (.not.(hf_to_psit .and. walk_dets_up(iwalk)==dets_up_psi_t(1) .and. walk_dets_dn(iwalk)==dets_dn_psi_t(1))) then ! For first state, all moves are deterministic, so skip stochastic step, as well as diagonal move.

    ! Don't spawn low-wt walkers every time.
    if (abs(walk_wt(iwalk)) < always_spawn_cutoff_wt) then
      spawn = (rannyu() < abs(walk_wt(iwalk)/always_spawn_cutoff_wt))
      use_wt = .false.
    else
      spawn = .true.
      use_wt = .true.
    endif

    if (spawn) then

      if (use_wt) then
        nwalk_child=max(nint(abs(walk_wt(iwalk))),1)
        walk_wt_child=walk_wt(iwalk)/nwalk_child
      else
        nwalk_child=1
        walk_wt_child=sign(always_spawn_cutoff_wt,walk_wt(iwalk))
      endif
      do iwalk_child=1,nwalk_child

    ! Pick an off-diagonal move and figure out its weight
!*** edited by AR [27/11/2013]: using interfaces
!         call off_diagonal_move(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, importance_sampling)
        if(hamiltonian_type.eq.'heg') then
          call off_diagonal_move_heg(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j)
        elseif(hamiltonian_type.eq.'chem') then
          if (use_efficient_heatbath) then
            ! Removing choose_pairs for now, but can put it back in later (for no-time_sym only!)
            call off_diagonal_move_chem_efficient_heatbath(walk_dets_up(iwalk), walk_dets_dn(iwalk), tmp_dets_up, tmp_dets_dn, weight_j_list,iwalk_child,importance_sampling,n_new_dets,excite_level)
            if (n_new_dets==1) then
              call add_walker(tmp_dets_up(1),tmp_dets_dn(1),walk_wt_child,weight_j_list(1),tau,iwalk,excite_level(1))
              cycle
            else
              call add_walker(tmp_dets_up(1),tmp_dets_dn(1),walk_wt_child,weight_j_list(1),tau,iwalk,excite_level(1))
              call add_walker(tmp_dets_up(2),tmp_dets_dn(2),walk_wt_child,weight_j_list(2),tau,iwalk,excite_level(2))
              cycle
            endif
          elseif(proposal_method.eq.'CauchySchwarz') then
            call off_diagonal_move_chem_CauchySchwarz(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, importance_sampling, iwalk_child)
          else
            call off_diagonal_move_chem(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, importance_sampling, iwalk_child)
          endif
          ! This part only executed if use_efficient_heatbath=.false.
          if (off_diag_histograms) then
            call excitation_level(walk_dets_up(iwalk),walk_dets_dn(iwalk),det_j_up,det_j_dn,excite_level(1))
            if (excite_level(1)==1) then
              call add_to_hist(abs(weight_j)/tau,nbins,lbounds,bins_single)
            elseif (excite_level(1)==2) then
              call add_to_hist(abs(weight_j)/tau,nbins,lbounds,bins_double)
            endif
          endif
        elseif(hamiltonian_type.eq.'hubbard2') then
          call off_diagonal_move_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j)
        elseif(hamiltonian_type.eq.'hubbardk') then
          call off_diagonal_move_hubbard_k(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, importance_sampling)
        elseif(hamiltonian_type.eq.'hubbarddm') then
          call off_diagonal_move_hubbard_dm(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j)
        endif
!**** end edit AR
        if (off_diag_histograms) then
          call add_to_hist(abs(weight_j)/tau,nbins,lbounds,bins)
          if (iwalk==1) then
            call add_to_hist(abs(weight_j)/tau,nbins,lbounds_hf,bins_hf)
          endif
        endif

! Warning: tmp
! Save the max and min ratios
!       if(abs(weight_j).gt.tau*max_wt_ratio .or. abs(weight_j).gt.1000._rk) then
!!        write(6,*) 'walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, proposal_prob=', walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j/tau, proposal_prob
!         write(6,'(''walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, proposal_prob='',4i20,2es9.2)') walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j/tau ! , proposal_prob ! Warning: the format should depend on how many words are used for determinant labels, i.e., num_words in types.f90
!       endif

        if(abs(weight_j)/tau.gt.1e-9) then
           if(abs(weight_j).gt.tau*max_wt_ratio) then
             write(6,*) 'walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j=', walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j/tau
            !write(6,'(''walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j='',4i20,f9.2)') walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j/tau
             max_wt_ratio = abs(weight_j)/tau
             if(time_sym) then
               max_wt_ratio_excit_level=min(count_excitations(walk_dets_up(iwalk),det_j_up)+count_excitations(walk_dets_dn(iwalk),det_j_dn),count_excitations(walk_dets_up(iwalk),det_j_dn)+count_excitations(walk_dets_dn(iwalk),det_j_up))
             else
               max_wt_ratio_excit_level=count_excitations(walk_dets_up(iwalk),det_j_up)+count_excitations(walk_dets_dn(iwalk),det_j_dn)
             endif
           endif
           if(abs(weight_j).lt.tau*min_wt_ratio) min_wt_ratio = abs(weight_j)/tau
        endif

        if(importance_sampling.eq.0) then
          weight_j=walk_wt_child*weight_j
        else
          weight_j=walk_wt_child*weight_j!* psi_g(walk_det_new)/psi_g(walk_det(iwalk))
        endif
        if (off_diag_histograms) then
          call add_to_hist(abs(weight_j)/tau,nbins,lbounds_new,bins_new)
        endif

        if(ipr.ge.1) write(6,*) 'walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j=', walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j
       !if(ipr.ge.1) write(6,'(''walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j'',4i20,f8.5)') walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j

        ! Edited by AAH on 15 Jan 2013. Don't allow stochastic spawning onto first state.
        if (hf_to_psit .and. det_j_up==dets_up_psi_t(1) .and. det_j_dn==dets_dn_psi_t(1))  weight_j=0._rk
        ! End of edit

        if(weight_j.ne.0) then

          nwalk=nwalk+1
          if (ncores==1)  snd_cnt(1)=snd_cnt(1)+1

          if(nwalk.gt.MWALK) then
            write(6,'(''stop move_uniform2: nwalk>MWALK, e_est, e_trial, reweight_factor=''9d12.4)') e_est, e_trial, 1._rk/reweight_factor_inv
            write(6,'(''nwalk='',i9)') nwalk
            write(6,'(1000(es9.1,i2))') (walk_wt(i), imp_distance(i),i=1,nwalk)
            call flush(6)
            stop 'nwalk>MWALK'
          endif

          !****B building of new walker
          walk_dets_up(nwalk)=det_j_up
          walk_dets_dn(nwalk)=det_j_dn
          walk_wt(nwalk)=weight_j

    ! Warning: Impose spin-reversal symmetry assuming +ve relative sign
    !     if(walk_dets_up(nwalk).gt.walk_dets_dn(nwalk)) then
    !       det_j_up=walk_dets_up(nwalk)
    !       walk_dets_up(nwalk)=walk_dets_dn(nwalk)
    !       walk_dets_dn(nwalk)=det_j_up
    !     endif
          if (imp_distance(iwalk)==-2) then
            ! Edited by AAH on 18 Jan 2013. Treat imp_distance=-2 as effectively 1 if c_t_initiator is false
            if (c_t_initiator) then
              imp_distance(nwalk)=1_i1b
            else
              imp_distance(nwalk)=2_i1b
            endif
            ! End of edit
          else
            imp_distance(nwalk)=min(imp_distance(iwalk),huge(1_i1b)-1_i1b)+1_i1b ! Increment imp_distance making sure not to exceed largest i1b integer
          endif
          if(semistochastic) then
            ! Prevents stochastic spawning between important dets
            if (imp_distance(iwalk).eq.0)  imp_distance(nwalk)=-1
          endif
          if(initiator(iwalk).ge.2) then ! If it is spawned by an initiator then 1 otherwise 0
            initiator(nwalk)=1
          else
            initiator(nwalk)=0
          endif
          if (c_t_initiator.and.imp_distance(iwalk)==-2)  initiator(nwalk)=1 ! if c_t_initiator true, then always spawn from C(T)
          if(semistochastic) then
             if (imp_distance(iwalk).eq.0) initiator(nwalk)=1 !always spawn from important space FP
          endif
!Warning: Is this in the right place for a parallel run?
          e_num_walker(nwalk)=1.e51_rk
          e_den_walker(nwalk)=1.e51_rk
          matrix_elements(nwalk)=1.e51_rk
          if (ncores>1)  call mpi_push_nwalk(nwalk,initiator(nwalk),e_num_walker(nwalk),e_den_walker(nwalk),matrix_elements(nwalk))
          !****E building of new walker

    !   else
    !     stop 'move_uniform2: weight_j=0'
        endif

      enddo ! iwalk_child

    endif ! spawn

  ! Diagonal move
    if((.not.semistochastic) .or. imp_distance(iwalk).ge.1) then ! Do diagonal move only if walker is outside deterministic space (and psi_t_connection space if hf_to_psit is true)

      if (matrix_elements(iwalk).gt.1.e50_rk) then
        if(hamiltonian_type.eq.'heg') then
          call hamiltonian_heg(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_matrix_element)
        elseif(hamiltonian_type.eq.'chem') then
          if (time_sym) then
            call hamiltonian_chem_time_sym(walk_dets_up(iwalk),walk_dets_dn(iwalk),walk_dets_up(iwalk),walk_dets_dn(iwalk),diagonal_matrix_element)
          else
            call hamiltonian_chem(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), 0, diagonal_matrix_element)
          endif
        elseif(hamiltonian_type.eq.'hubbard2') then
          call hamiltonian_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_matrix_element)
        elseif(hamiltonian_type.eq.'hubbardk') then
          if (space_sym) then
            call hamiltonian_hubbard_k_space_sym(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_matrix_element,nnzero)
          else
            call hamiltonian_hubbard_k(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_matrix_element)
          endif
        endif
        matrix_elements(iwalk) = diagonal_matrix_element

      else

        diagonal_matrix_element = matrix_elements(iwalk)

      endif

      diagonal_factor=1._rk+tau*(e_trial-diagonal_matrix_element)

      if(diagonal_factor.lt.0) then ! It is OK for diagonal_factor to be < 0 during early equil when tau is large, but not later on.
        if(reached_w_abs_gen.le.1) then
          count_diag_neg=count_diag_neg+1
          if(count_diag_neg.le.100) then
            write (fmt, '(i2,''i'',i4)') (2*num_words), bits_per_word
            write(6,'(''Warning: For determinant walk_dets_up(iwalk), walk_dets_dn(iwalk)='',' // trim(fmt) // ','' walk_wt multiplier is negative'',f10.6,'' Reduce tau='',f10.5)') &
           &walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_factor, tau
          endif
          if(count_diag_neg.eq.100) then
            write(6,'(''Warning: No more walk_wt multiplier is negative msgs. will be written'')')
          endif
        else
          write (fmt, '(i2,''i'',i4)') (2*num_words), bits_per_word
          write(6,'(''Warning: For determinant walk_dets_up(iwalk), walk_dets_dn(iwalk)='',' // trim(fmt) // ','' walk_wt multiplier is negative'',f10.6,'' Reduce tau='',f10.5)') &
         &walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_factor, tau
          call mpi_stop('diagonal_factor<0 after target population has been reached')
        endif
        diagonal_factor=0
      endif

      walk_wt(iwalk)=walk_wt(iwalk)*diagonal_factor ! Diagonal move
      if(ipr.ge.1) write(6,*) 'Diagonal: iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)=', iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)
     !if(ipr.ge.1) write(6,'(''Diagonal: iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)='',i10,2i20,f10.6)') iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)
    endif

  endif

  end subroutine move_uniform2
  !-------------------------------------------------------------------------------------

  subroutine move_uniform_exponential_projector(iwalk)
  ! ==========================================================================================================
  ! Description   : For each walker on a given determinant iwalk, perform (off-diagonal) moves by sampling
  !               : tau's until you reach the measurement time step tau_e (currently equal to the global
  !               : variable tau), at which point the off-diagonal move that took you across the tau_e
  !               : threshold is replaced by a diagonal move up until that time.
  ! Author        : A Holmes, 28 Aug 2013.
  ! ----------------------------------------------------------------------------------------------------------
    use hamiltonian_mod
    use semistoch
    use common_psi_t, only : hf_to_psit
    use hubbard, only : hamiltonian_hubbard, hamiltonian_hubbard_k, hamiltonian_hubbard_dm, hamiltonian_hubbard_k_space_sym, space_sym,U,nsites,find_connected_dets_hubbard_k
    use chemistry, only : hamiltonian_chem, hamiltonian_chem_time_sym, time_sym
    use heg, only : hamiltonian_heg
    use tools, only : count_excitations
    use types, only : num_words, bits_per_word

    implicit none
    integer, intent(in) :: iwalk
    integer iwalk_child, nwalk_child
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec) det_j_up, det_j_dn, my_up, my_dn
#else
    integer(ik) det_j_up, det_j_dn, my_up, my_dn
#endif
    ! real(rk) random, rannyu, diagonal_matrix_element, walk_wt_child
    real(rk) :: rannyu
    real(rk) walk_wt_child, weight_j
    integer :: i
    real(rk) :: H_off,sampled_tau,my_tau,elem,my_wt
    logical :: diagonal_move
    logical, parameter :: print_stuff = .false.
    real(rk) :: xi,stoch_wt
    logical, parameter :: deterministic_diagonal = .true.
    integer(i8b) :: nwalk_old

    nwalk_old = nwalk
    H_off_multiplier = 1.0

    if (print_stuff) then
      write (6,*) ""
      write (6,*) "iwalk=",iwalk; call flush(6)
      write (6,*) ""
    endif
    if (print_stuff .and. iwalk==1)  write (6,*) "wt before:",walk_wt(1:nwalk_old)

    H_off = -U/real(nsites)*nup*ndn*(nsites-nup)*H_off_multiplier ! just for hubbard: abs(off-diag element) * number of off-diagonal elements

    if (deterministic_diagonal) then
      call hamiltonian_hubbard_k(walk_dets_up(iwalk),walk_dets_dn(iwalk),walk_dets_up(iwalk),walk_dets_dn(iwalk),elem)
      stoch_wt = walk_wt(iwalk) * (1._rk-exp(H_off*tau))
      walk_wt(iwalk) = walk_wt(iwalk) * exp(tau*(e_trial-elem))
      !write (6,*) "tau=",tau,"H_off=",H_off
      !write (6,*) "diag",exp(tau*(H_off)),"off_diag",(1._rk-exp(tau*H_off)); call flush(6)
      !write (6,*) "walk_wt=",walk_wt(iwalk)
      !write (6,*) "stoch_wt=",stoch_wt
    else
      stoch_wt = walk_wt(iwalk)
    endif

    nwalk_child=max(nint(abs(stoch_wt)),1)
    walk_wt_child=stoch_wt/nwalk_child
    xi = rannyu()

    if (print_stuff)  write (6,*) "nwalk_child=",nwalk_child,"walk_wt_child=",walk_wt_child; call flush(6)

    !write (6,*) "size(connected_dets_up)",size(connected_dets_up) ; call flush(6)
    !deallocate(connected_dets_up,connected_dets_dn)
    !if (.not.allocated(connected_dets_up)) then
    !  allocate(connected_dets_up(19))
    !  allocate(connected_dets_dn(19))
    !endif
    !call find_connected_dets_hubbard_k(walk_dets_up(iwalk),walk_dets_dn(iwalk),n_con,connected_dets_up,connected_dets_dn)
    !H_off = n_con*(-U/real(nsites))
    !e_trial = -0.286_rk
    !tau = 0.01_rk
    if (print_stuff)  write (6,*) "H_off=",H_off,"e_trial=",e_trial,"tau=",tau; call flush(6)

    do iwalk_child=1,nwalk_child
      my_tau = 0._rk
      my_up = walk_dets_up(iwalk)
      my_dn = walk_dets_dn(iwalk)
      my_wt = walk_wt_child
      diagonal_move = .false.
      if (.not.deterministic_diagonal)  xi = mod(xi+1._rk/nwalk_child,1._rk)

      do i=1,1000 ! should do forever
        if (i>900)  stop "Infinite loop!"

        ! First, pick a tau from the exponential distribution exp(-tau*H_off)
        if (deterministic_diagonal) then
          if (i==1) then
            sampled_tau = mod(log(rannyu())/H_off,tau)
          else
            sampled_tau = log(rannyu())/H_off
          endif
        else
          if (i==1) then
            sampled_tau = log(xi)/H_off
          else
            sampled_tau = log(rannyu())/H_off
          endif
        endif

        if (print_stuff)  write (6,*) "sampled_tau=",sampled_tau

        if (my_tau + sampled_tau > tau) then
          ! Update this piece of the wt and go to next walker
          call hamiltonian_hubbard_k(my_up,my_dn,my_up,my_dn,elem)
          if (print_stuff)  write (6,*) "H_ii=",elem; call flush(6)
          if (i==1) then
            if (deterministic_diagonal)  stop "Should not get here!"
            diagonal_move = .true.
            walk_wt(iwalk) = walk_wt(iwalk) + walk_wt_child * (exp((tau-my_tau)*(e_trial-elem-H_off)) - 1._rk)
            ! DELETE THIS:
            my_wt=walk_wt_child *exp((tau-my_tau)*(e_trial-elem-H_off))
          else ! This move is diagonal but an off-diagonal move has happened before
            my_wt = my_wt * exp((tau-my_tau)*(e_trial-elem-H_off))
            !walk_wt(iwalk) = walk_wt(iwalk) - walk_wt_child
            if (.not.deterministic_diagonal)  walk_wt(iwalk) = walk_wt(iwalk) - walk_wt_child
          endif
          exit
        endif

        ! Pick new state j with uniform probability (for now)

        !weight_j = 0._rk
        !do while (abs(weight_j)<1.e-10_rk)
        if(hamiltonian_type.eq.'heg') then
          call off_diagonal_move_heg(my_up, my_dn, det_j_up, det_j_dn, weight_j)
        elseif(hamiltonian_type.eq.'chem') then
          if (use_efficient_heatbath) then
            call mpi_stop('Efficient heatbath not yet implemented for this run mode!')
           !call off_diagonal_move_chem_efficient_heatbath(my_up, my_dn, det_j_up, det_j_dn, weight_j,iwalk_child,(nwalk_child>nelec),importance_sampling)
          elseif(proposal_method.eq.'CauchySchwarz') then
            call off_diagonal_move_chem_CauchySchwarz(my_up, my_dn, det_j_up, det_j_dn, weight_j, importance_sampling, iwalk_child)
          else
            call off_diagonal_move_chem(my_up, my_dn, det_j_up, det_j_dn, weight_j, importance_sampling, iwalk_child)
          endif
        elseif(hamiltonian_type.eq.'hubbard2') then
          call off_diagonal_move_hubbard(my_up, my_dn, det_j_up, det_j_dn, weight_j)
        elseif(hamiltonian_type.eq.'hubbardk') then
          call off_diagonal_move_hubbard_k(my_up, my_dn, det_j_up, det_j_dn, weight_j, importance_sampling)
        elseif(hamiltonian_type.eq.'hubbarddm') then
          call off_diagonal_move_hubbard_dm(my_up, my_dn, det_j_up, det_j_dn, weight_j)
        endif
        !enddo
        if (off_diag_histograms) then
          call add_to_hist(abs(weight_j)/tau,nbins,lbounds,bins)
          if (iwalk==1) then
            call add_to_hist(abs(weight_j)/tau,nbins,lbounds_hf,bins_hf)
          endif
        endif

        !Save the max and min ratios
        if(abs(weight_j)/tau.gt.1e-9) then
           if(abs(weight_j).gt.tau*max_wt_ratio) then
             max_wt_ratio = abs(weight_j)/tau
             if(time_sym) then
               max_wt_ratio_excit_level=min(count_excitations(walk_dets_up(iwalk),det_j_up)+count_excitations(walk_dets_dn(iwalk),det_j_dn),count_excitations(walk_dets_up(iwalk),det_j_dn)+count_excitations(walk_dets_dn(iwalk),det_j_up))
             else
               max_wt_ratio_excit_level=count_excitations(walk_dets_up(iwalk),det_j_up)+count_excitations(walk_dets_dn(iwalk),det_j_dn)
             endif
           endif
           if(abs(weight_j).lt.tau*min_wt_ratio) min_wt_ratio = abs(weight_j)/tau
        endif

        call hamiltonian_hubbard_k(my_up,my_dn,my_up,my_dn,elem)
        if (print_stuff)  write (6,*) "H_ii=",elem; call flush(6)
        if (print_stuff)  call hamiltonian_hubbard_k(my_up,my_dn,det_j_up,det_j_dn,elem)
        if (print_stuff)  write (6,*) "H_ij=",-weight_j/tau/(nup*ndn*(nsites-nup)),elem,U/real(nsites); call flush(6)

        ! Update tau counter
        my_tau = my_tau + sampled_tau

        ! Move to the new det
        my_up = det_j_up
        my_dn = det_j_dn
        my_wt=-my_wt*(-weight_j/tau/(nup*ndn*(nsites-nup)))*exp(sampled_tau*(e_trial-elem-H_off))/(H_off_multiplier*U/real(nsites))

      enddo

      if (print_stuff)  write (6,*) "diagonal_move=",diagonal_move,"my_up=",my_up,"my_dn=",my_dn,"my_wt=",my_wt; call flush(6)

      if(ipr.ge.1) write(6,'(''walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j'',4i20,f8.5)') walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j

      if(my_wt.ne.0.and..not.diagonal_move) then
        if (print_stuff)  write (6,*) "This is an off-diagonal move" ; call flush(6)
        nwalk=nwalk+1
        snd_cnt(1)=snd_cnt(1)+1
        if(nwalk.gt.MWALK) then
          write(6,'(''stop move_uniform_exponential_projector: nwalk>MWALK, e_est, e_trial, reweight_factor=''9d12.4)') e_est, e_trial, 1._rk/reweight_factor_inv
          call flush(6)
          stop 'nwalk>MWALK'
        endif
        walk_dets_up(nwalk)=my_up
        walk_dets_dn(nwalk)=my_dn
        walk_wt(nwalk)=my_wt

        ! Warning: Impose spin-reversal symmetry assuming +ve relative sign
        !     if(walk_dets_up(nwalk).gt.walk_dets_dn(nwalk)) then
        !       det_j_up=walk_dets_up(nwalk)
        !       walk_dets_up(nwalk)=walk_dets_dn(nwalk)
        !       walk_dets_dn(nwalk)=det_j_up
        !     endif
        if (imp_distance(iwalk)==-2) then
          ! Edited by AAH on 18 Jan 2013. Treat imp_distance=-2 as effectively 1 if c_t_initiator is false
          if (c_t_initiator) then
            imp_distance(nwalk)=1_i1b
          else
            imp_distance(nwalk)=2_i1b
          endif
          ! End of edit
        else
          imp_distance(nwalk)=min(imp_distance(iwalk),huge(1_i1b)-1_i1b)+1_i1b ! Increment imp_distance making sure not to exceed largest i1b integer
        endif
        !if(semistochastic) then
        !  ! Prevents stochastic spawning between important dets
        !  if (imp_distance(iwalk).eq.0)  imp_distance(nwalk)=-1
        !endif
        if(initiator(iwalk).ge.2) then ! If it is spawned by an initiator then 1 otherwise 0
          initiator(nwalk)=1
        else
          initiator(nwalk)=0
        endif
        if (c_t_initiator.and.imp_distance(iwalk)==-2)  initiator(nwalk)=1 ! if c_t_initiator true, then always spawn from C(T)
        if(semistochastic) then
        if (imp_distance(iwalk).eq.0) initiator(nwalk)=1 !always spawn from important space FP
        endif
!Warning: Is this in the right place for a parallel run?
        e_num_walker(nwalk)=1.e51_rk
        e_den_walker(nwalk)=1.e51_rk
        matrix_elements(nwalk)=1.e51_rk
      endif

    enddo ! iwalk_child

    if (print_stuff .and. iwalk==1)  write (6,*) "wt after:",walk_wt(1:nwalk_old)


  end subroutine move_uniform_exponential_projector
!-------------------------------------------------------------------------------------

  subroutine move_uniform_improved_cancellation(nwalk)
! ==========================================================================================================
! Description   : Move all walkers as follows:
!               :   1) Pick one occupied det at random with probability proportional to its absolute wt, make uniform move
!               :   2) Use bit-wise xor to identify all occupied dets that could be connected to target det
!               :   3) Calculate contributions from all those dets to target det
! Author        : A Holmes, 14 Aug 2013
! ----------------------------------------------------------------------------------------------------------
  use hamiltonian_mod
  use semistoch
  use common_psi_t, only : hf_to_psit
  use hubbard, only : hamiltonian_hubbard, hamiltonian_hubbard_k, hamiltonian_hubbard_dm, hamiltonian_hubbard_k_space_sym, space_sym,is_connected_hubbard_fast
  use chemistry, only : hamiltonian_chem, hamiltonian_chem_time_sym, time_sym
  use heg, only : hamiltonian_heg
  use types, only : num_words, bits_per_word
  use tools, only : count_excitations, random_int, choose_det_with_prob_prop_to_abs_wt

  implicit none
  integer, intent(inout) :: nwalk
  integer iwalk,i
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec) det_j_up, det_j_dn
#else
  integer(ik) det_j_up, det_j_dn
#endif
! real(rk) random, rannyu, diagonal_matrix_element, walk_wt_child
  real(rk) diagonal_matrix_element, diagonal_factor, weight_j,proposal_prob_inv,elem
  character(len=7) :: fmt
  integer :: nnzero, nwalk_old, count_diag_neg
  real(rk) :: counter,tot_abs_wt
  data count_diag_neg/0/

  nwalk_old = nwalk

    call choose_det_with_prob_prop_to_abs_wt(iwalk,tot_abs_wt)

    ! Pick an off-diagonal move and figure out its weight
    ! Automatically divides by probability to pick from iwalk's connections -> this is exactly what we want!
    if(hamiltonian_type.eq.'heg') then
      call off_diagonal_move_heg(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j)
    elseif(hamiltonian_type.eq.'chem') then
      if (use_efficient_heatbath) then
        call mpi_stop('Efficient heatbath not yet implemented for this run mode!')
       !call off_diagonal_move_chem_efficient_heatbath(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j,1,.false.,importance_sampling)
      elseif(proposal_method.eq.'CauchySchwarz') then
        call off_diagonal_move_chem_CauchySchwarz(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, importance_sampling, 1)
      else
        call off_diagonal_move_chem(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, importance_sampling, 1)
      endif
    elseif(hamiltonian_type.eq.'hubbard2') then
      call off_diagonal_move_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j)
    elseif(hamiltonian_type.eq.'hubbardk') then
      call off_diagonal_move_hubbard_k(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, importance_sampling,proposal_prob_inv)
    elseif(hamiltonian_type.eq.'hubbarddm') then
      call off_diagonal_move_hubbard_dm(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j)
    endif

    if (off_diag_histograms) then
      call add_to_hist(abs(weight_j)/tau,nbins,lbounds,bins)
      if (iwalk==1) then
        call add_to_hist(abs(weight_j)/tau,nbins,lbounds_hf,bins_hf)
      endif
    endif

    !Save the max and min ratios
    if(abs(weight_j)/tau.gt.1e-9) then
       if(abs(weight_j).gt.tau*max_wt_ratio) then
         max_wt_ratio = abs(weight_j)/tau
         if(time_sym) then
           max_wt_ratio_excit_level=min(count_excitations(walk_dets_up(iwalk),det_j_up)+count_excitations(walk_dets_dn(iwalk),det_j_dn),count_excitations(walk_dets_up(iwalk),det_j_dn)+count_excitations(walk_dets_dn(iwalk),det_j_up))
         else
           max_wt_ratio_excit_level=count_excitations(walk_dets_up(iwalk),det_j_up)+count_excitations(walk_dets_dn(iwalk),det_j_dn)
         endif
       endif
       if(abs(weight_j).lt.tau*min_wt_ratio) min_wt_ratio = abs(weight_j)/tau
    endif

    weight_j=weight_j*walk_wt(iwalk)
    counter = abs(walk_wt(iwalk))

    ! Now new det = det_j_up,det_j_dn,  new det weight = weight_j

    ! Now for each walker, figure out whether it could spawn to det_j

    if (walk_dets_up(iwalk).eq.det_j_up.and.walk_dets_dn(iwalk).eq.det_j_dn) then
      weight_j = 0._rk
    else
      do i=1,nwalk_old
        if (i==iwalk)  cycle ! already taken care of
        if (walk_dets_up(i)==det_j_up.and.walk_dets_dn(i)==det_j_dn) cycle ! No diagonal moves allowed in stochastic piece
        if (is_connected_hubbard_fast(walk_dets_up(i),walk_dets_dn(i),det_j_up,det_j_dn)) then   ! Note: could we speed up is_con_hub_fast with count_excit??
          call hamiltonian_hubbard_k(walk_dets_up(i),walk_dets_dn(i),det_j_up,det_j_dn,elem,.true.)
          weight_j = weight_j -tau*elem*walk_wt(i)*proposal_prob_inv
          counter = counter + abs(walk_wt(i))
        endif
      enddo
    endif
    weight_j = weight_j * tot_abs_wt/counter

    if(weight_j.ne.0) then
      nwalk=nwalk+1
      snd_cnt(1)=snd_cnt(1)+1
      if(nwalk.gt.MWALK) then
        write(6,'(''stop move_uniform_improved_cancellation: nwalk>MWALK, e_est, e_trial, reweight_factor=''9d12.4)') e_est, e_trial, 1._rk/reweight_factor_inv
        call flush(6)
        stop 'nwalk>MWALK'
      endif
      walk_dets_up(nwalk)=det_j_up
      walk_dets_dn(nwalk)=det_j_dn
      walk_wt(nwalk)=weight_j

! Warning: Impose spin-reversal symmetry assuming +ve relative sign
!     if(walk_dets_up(nwalk).gt.walk_dets_dn(nwalk)) then
!       det_j_up=walk_dets_up(nwalk)
!       walk_dets_up(nwalk)=walk_dets_dn(nwalk)
!       walk_dets_dn(nwalk)=det_j_up
!     endif
      if (imp_distance(iwalk)==-2) then
        ! Edited by AAH on 18 Jan 2013. Treat imp_distance=-2 as effectively 1 if c_t_initiator is false
        if (c_t_initiator) then
          imp_distance(nwalk)=1_i1b
        else
          imp_distance(nwalk)=2_i1b
        endif
        ! End of edit
      else
        imp_distance(nwalk)=min(imp_distance(iwalk),huge(1_i1b)-1_i1b)+1_i1b ! Increment imp_distance making sure not to exceed largest i1b integer
      endif
      if(semistochastic) then
        ! Prevents stochastic spawning between important dets
        if (imp_distance(iwalk).eq.0)  imp_distance(nwalk)=-1
      endif
      if(initiator(iwalk).ge.2) then ! If it is spawned by an initiator then 1 otherwise 0
        initiator(nwalk)=1
      else
        initiator(nwalk)=0
      endif
      if (c_t_initiator.and.imp_distance(iwalk)==-2)  initiator(nwalk)=1 ! if c_t_initiator true, then always spawn from C(T)
      if(semistochastic) then
         if (imp_distance(iwalk).eq.0) initiator(nwalk)=1 !always spawn from important space FP
      endif
!Warning: Is this in the right place for a parallel run?
      e_num_walker(nwalk)=1.e51_rk
      e_den_walker(nwalk)=1.e51_rk
      matrix_elements(nwalk)=1.e51_rk
    endif

  ! Diagonal moves
  do iwalk=1,nwalk_old
    if((.not.semistochastic) .or. imp_distance(iwalk).ge.1) then ! Do diagonal move only if walker is outside deterministic space (and psi_t_connection space if hf_to_psit is true)

      if (matrix_elements(iwalk).gt.1.e50_rk) then
        if(hamiltonian_type.eq.'heg') then
          call hamiltonian_heg(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_matrix_element)
        elseif(hamiltonian_type.eq.'chem') then
          if (time_sym) then
            call hamiltonian_chem_time_sym(walk_dets_up(iwalk),walk_dets_dn(iwalk),walk_dets_up(iwalk),walk_dets_dn(iwalk),diagonal_matrix_element)
          else
            call hamiltonian_chem(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), 0, diagonal_matrix_element)
          endif
        elseif(hamiltonian_type.eq.'hubbard2') then
          call hamiltonian_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_matrix_element)
        elseif(hamiltonian_type.eq.'hubbardk') then
          if (space_sym) then
            call hamiltonian_hubbard_k_space_sym(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_matrix_element,nnzero)
          else
            call hamiltonian_hubbard_k(walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_matrix_element)
          endif
        endif
        matrix_elements(iwalk) = diagonal_matrix_element

      else

        diagonal_matrix_element = matrix_elements(iwalk)

      endif

      diagonal_factor=1._rk+tau*(e_trial-diagonal_matrix_element)

      if(diagonal_factor.lt.0) then ! It is OK for diagonal_factor to be < 0 during early equil when tau is large, but not later on.
        if(reached_w_abs_gen.le.1) then
          count_diag_neg=count_diag_neg+1
          if(count_diag_neg.le.100) then
            write (fmt, '(i2,''i'',i4)') (2*num_words), bits_per_word
            write(6,'(''Warning: For determinant walk_dets_up(iwalk), walk_dets_dn(iwalk)='',' // trim(fmt) // ','' walk_wt multiplier is negative'',f10.6,'' Reduce tau='',f10.5)') &
           &walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_factor, tau
          endif
          if(count_diag_neg.eq.100) then
            write(6,'(''Warning: No more walk_wt multiplier is negative msgs. will be written'')')
          endif
        else
          write (fmt, '(i2,''i'',i4)') (2*num_words), bits_per_word
          write(6,'(''Warning: For determinant walk_dets_up(iwalk), walk_dets_dn(iwalk)='',' // trim(fmt) // ','' walk_wt multiplier is negative'',f10.6,'' Reduce tau='',f10.5)') &
         &walk_dets_up(iwalk), walk_dets_dn(iwalk), diagonal_factor, tau
          call mpi_stop('diagonal_factor<0 after target population has been reached')
        endif
        diagonal_factor=0
      endif

      walk_wt(iwalk)=walk_wt(iwalk)*diagonal_factor ! Diagonal move
      if(ipr.ge.1) write(6,'(''Diagonal: iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)='',i10,2i20,f10.6)') iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)
    endif
  enddo

  end subroutine move_uniform_improved_cancellation
!-------------------------------------------------------------------------------------

  subroutine move_uniform_vmc(iwalk)
! ==========================================================================================================
! Description   : Do off-diagonal move choosing uniformly from n_connected_dets elements of the propagator
!               : for VMC
! Notes         : Under testing and works ONLY in 1 walker mode at the moment
! Author        : Hitesh J Changlani May 28-31,2011
! ----------------------------------------------------------------------------------------------------------
  use hamiltonian_mod
  implicit none
  integer, intent(in) :: iwalk
! integer iwalk_child, weight_j
  integer iwalk_child
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec) det_j_up, det_j_dn
#else
  integer(ik) det_j_up, det_j_dn
#endif
  real(rk) rannyu, wf_ratio, ratio_square, weight_j

! Warning: tmp
! do iwalk_child=1,abs(walk_wt(iwalk))                                     ! Works ONLY in 1 walker mode
  do iwalk_child=1,1                                                       ! Works ONLY in 1 walker mode

! Pick an off-diagonal move and figure out its weight
    if(hamiltonian_type.eq.'heg') then
      !call off_diagonal_move_heg(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, diagonal_matrix_element)
      stop 'VMC not yet implemented in heg'
    elseif(hamiltonian_type.eq.'chem') then
      !call off_diagonal_move_chem(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, diagonal_matrix_element, iwalk_child)
      stop 'VMC not yet implemented in chem'
    elseif(hamiltonian_type.eq.'hubbard2') then
      call off_diagonal_move_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, wf_ratio,.true.)
     !elseif(hamiltonian_type.eq.'hubbardk') then
     ! call off_diagonal_move_hubbard_k(walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j, wf_ratio,.true.)
    endif

    if(ipr.ge.1) write(6,'(''walk_dets_up(iwalk),walk_dets_dn(iwalk),det_j_up, det_j_dn, weight_j'',4i20,f10.6)') walk_dets_up(iwalk),walk_dets_dn(iwalk),det_j_up, det_j_dn, weight_j

    ratio_square=(wf_ratio)**(2._rk)

    if (ratio_square>rannyu()) then ! accept according to metropolis
        walk_dets_up(nwalk)=det_j_up
        walk_dets_dn(nwalk)=det_j_dn
        walk_wt(nwalk)=1                                                  ! Works ONLY in 1 walker mode
    endif

   enddo ! iwalk_child

  end subroutine move_uniform_vmc
!-------------------------------------------------------------------------------------

  subroutine move_heat_bath(iwalk)
! ====================================================================================================================
! Description   : Do off-diagonal move using heat bath on the absolute value of the elements of the appropriate column of the propagator
!               : Version for problems small enough that Hamiltonian can be stored (hamiltonian_type fictitious, read, hubbard)
! Author        : Cyrus Umrigar
! -------------------------------------------------------------------------------------------------------------------
  use hamiltonian_mod
  implicit none
  integer, intent(in) :: iwalk
  integer walk_det_new, i, i_nonzero, iwalk_child, nwalk_child
  real(rk) random, rannyu, probability, sgn, walk_wt_tmp, walk_wt_child

! Count the number of non-zero off-diagonal elements and set cumulative probabilities
  n_connected_dets=0
  do i=1,ndet
    if(i.ne.walk_det(iwalk) .and. ham(i,walk_det(iwalk)).ne.0._rk) then
      n_connected_dets=n_connected_dets+1
      connected_dets(n_connected_dets)=i
      if(importance_sampling.eq.0) then
        probability=abs(ham(i,walk_det(iwalk)))
      else
        probability=abs(ham(i,walk_det(iwalk)) * psi_g(i)/psi_g(walk_det(iwalk)))
      endif
      if(n_connected_dets.eq.1) then
        cumulative_probability(n_connected_dets)=probability
      else
        cumulative_probability(n_connected_dets)=cumulative_probability(n_connected_dets-1)+probability
      endif
    endif
  enddo

! Warning: should be changed to:
! nwalk_child=max(nint(abs(walk_wt(iwalk))),1)
  nwalk_child=max(nint(abs(walk_wt(iwalk)*tau*cumulative_probability(n_connected_dets))),1)
  walk_wt_child=walk_wt(iwalk)/nwalk_child
! do iwalk_child=1,abs(walk_wt(iwalk))
  do iwalk_child=1,nwalk_child

! Decide which of these states to move to
    random=cumulative_probability(n_connected_dets)*rannyu()

    do i_nonzero=1,n_connected_dets
      if(cumulative_probability(i_nonzero).gt.random) then
        walk_det_new=connected_dets(i_nonzero)
        exit
      endif
    enddo

! Calculate wt of this determinant
    if(importance_sampling.eq.0) then
      sgn=-sign(1._rk,ham(walk_det_new,walk_det(iwalk))*walk_wt(iwalk))
    else
      sgn=-sign(1._rk,ham(walk_det_new,walk_det(iwalk))*walk_wt(iwalk)*psi_g(walk_det_new)/psi_g(walk_det(iwalk)))
    endif
!   walk_wt_tmp=sgn*int(reweight_factor_inv*tau*cumulative_probability(n_connected_dets)+rannyu())
!   walk_wt_tmp=sgn*(reweight_factor_inv*tau*cumulative_probability(n_connected_dets))*walk_wt_child      ! reweighting done in subroutine walk
    walk_wt_tmp=sgn*(tau*cumulative_probability(n_connected_dets))*walk_wt_child

! write(6, '(''reweight_factor, tau, proposal_prob_inv, ham(walk_det_new, walk_det(iwalk)), random'', 9d12.4)') &
!& 1._rk/reweight_factor_inv, tau, proposal_prob_inv, ham(walk_det_new, walk_det(iwalk)), random
    if(walk_wt_tmp.ne.0) then
      nwalk=nwalk+1
      snd_cnt(1)=snd_cnt(1)+1
      if(nwalk.gt.MWALK) then
        write(6,'(''stop move_heat_bath: nwalk>MWALK, e_est, e_trial, reweight_factor=''9d12.4)') e_est, e_trial, 1._rk/reweight_factor_inv
        call flush(6)
        stop 'nwalk>MWALK'
      endif
      walk_det(nwalk)=walk_det_new
      walk_wt(nwalk)=walk_wt_tmp
      if(initiator(iwalk).ge.2) then ! If it is spawned by an initiator then 1 otherwise 0
        initiator(nwalk)=1
      else
        initiator(nwalk)=0
      endif
!Warning: Is this in the right place for a parallel run?
      e_num_walker(nwalk)=1.e51_rk
    endif

  enddo ! iwalk_child

! random=sign(rannyu(),walk_wt(iwalk)*(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk)))))
! walk_wt(iwalk)=int(reweight_factor_inv*walk_wt(iwalk)*(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk))))+random) ! Diagonal move
! walk_wt(iwalk)=reweight_factor_inv*walk_wt(iwalk)*(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk)))) ! Diagonal move !reweighting done in subroutine walk
  walk_wt(iwalk)=walk_wt(iwalk)*(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk)))) ! Diagonal move
  if(ipr.ge.1) write(6,'(''Diagonal: iwalk, walk_det(iwalk), walk_wt(iwalk)='',2i5,f10.6)') iwalk, walk_det(iwalk), walk_wt(iwalk)

  end subroutine move_heat_bath
!-------------------------------------------------------------------------------------

  subroutine move_heat_bath2(iwalk)
! ==============================================================================
! Description   : Do move using heat bath on the absolute value of the elements of the appropriate column of the propagator
!               : Version for problems small enough that Hamiltonian can be stored (hamiltonian_type fictitious, read, hubbard)
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use hamiltonian_mod
  implicit none
  integer, intent(in) :: iwalk
  integer i, i_nonzero, nwalk_new, nwalk_this_det, nwalk_all_det
  real(rk) random, rannyu, probability
  integer diag_conn
  real(rk) probability_vec(size(connected_dets,1)), walk_wt_child, diagonal_dump

  if (run_type .ne. 'fixed_node1') then
    ! Count the number of non-zero off-diagonal elements and set cumulative probabilities
      n_connected_dets=0
      do i=1,ndet
        if(i.eq.walk_det(iwalk) .or. ham(i,walk_det(iwalk)).ne.0._rk) then
          n_connected_dets=n_connected_dets+1
          connected_dets(n_connected_dets)=i
          if(i.eq.walk_det(iwalk)) then
            probability=(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk))))
          else
            if(importance_sampling.eq.0) then
              probability=-tau*ham(i,walk_det(iwalk))
            else
              if(run_type.eq.'fixed_node2') then
                probability=abs(-tau*ham(i,walk_det(iwalk)) * psi_g(i)/psi_g(walk_det(iwalk)))
              else
                probability=-tau*ham(i,walk_det(iwalk)) * psi_g(i)/psi_g(walk_det(iwalk))
              endif
            endif
          endif
          !probability=reweight_factor_inv*walk_wt(iwalk)*probability    ! reweighting done in subroutine walk
          probability=walk_wt(iwalk)*probability
          connected_dets_sign(n_connected_dets)=int(sign(1._rk,probability))
          if(n_connected_dets.eq.1) then
            cumulative_probability(n_connected_dets)=abs(probability)
          else
            cumulative_probability(n_connected_dets)=cumulative_probability(n_connected_dets-1)+abs(probability)
          endif
        endif
      enddo
  else ! fixed_node1
      ! Count the number of non-zero off-diagonal elements and set cumulative probabilities
      n_connected_dets=0
      diagonal_dump=0
      do i=1,ndet
        if(i.eq.walk_det(iwalk) .or. ham(i,walk_det(iwalk)).ne.0._rk) then
          n_connected_dets=n_connected_dets+1
          connected_dets(n_connected_dets)=i
          if(i.eq.walk_det(iwalk)) then
            probability=(1+tau*(e_trial-ham(walk_det(iwalk),walk_det(iwalk))))
            diag_conn=n_connected_dets
          else
            if(importance_sampling.eq.0) then
              !probability=-tau*ham(i,walk_det(iwalk))*sign(1._rk,psi_g(i)*psi_g(walk_det(iwalk)))
              probability=-tau*ham(i,walk_det(iwalk))
            else
              probability=-tau*ham(i,walk_det(iwalk)) * psi_g(i)/psi_g(walk_det(iwalk))
            endif
            if(probability .lt. 0._rk) then
              diagonal_dump=diagonal_dump+probability
              probability=0._rk
            endif
          endif
          probability_vec(n_connected_dets)=probability
        endif
      enddo
      probability_vec(diag_conn)=probability_vec(diag_conn)+diagonal_dump
      !write(6,'(''diag_conn,probability_vec(diag_conn),diagonal_dump='',i5,9es12.4)') diag_conn,probability_vec(diag_conn),diagonal_dump
      if(probability_vec(diag_conn).lt.0) then
        write(6,'(''probability_vec(diag_conn)='',es12.4)') probability_vec(diag_conn)
        stop 'probability_vec(diag_conn).lt.0'
      endif

      do i=1,n_connected_dets
        !probability_vec(i)=reweight_factor_inv*walk_wt(iwalk)*probability_vec(i) !reweighting in subroutine walk
        probability_vec(i)=walk_wt(iwalk)*probability_vec(i)
        connected_dets_sign(i)=int(sign(1._rk,probability_vec(i)))
        if(i.eq.1) then
          cumulative_probability(i)=abs(probability_vec(i))
        else
          cumulative_probability(i)=cumulative_probability(i-1)+abs(probability_vec(i))
        endif
      enddo
  endif

  if(cumulative_probability(n_connected_dets).le.0._rk) then
    write(6,'(''Warning: cumulative_probability(n_connected_dets)='',es12.4)') cumulative_probability(n_connected_dets)
    stop 'cumulative_probability(n_connected_dets).le.0._rk'
  endif

! Decide on total absolute walker wts coming from this walker
! nwalk_new=int(cumulative_probability(n_connected_dets)+rannyu())
  nwalk_new=max(nint(cumulative_probability(n_connected_dets)),1)
  walk_wt_child=cumulative_probability(n_connected_dets)/nwalk_new

  if(ipr.ge.1) write(6,'(/,''iwalk, walk_det(iwalk), walk_wt(iwalk), n_connected_dets, nwalk_new='',2i5,f7.1,9i5)') &
 & iwalk, walk_det(iwalk), walk_wt(iwalk), n_connected_dets, nwalk_new
  if(ipr.ge.1) write(6,'(''reweight_factor'',f10.4)') 1._rk/reweight_factor_inv
  if(ipr.ge.1) write(6,'(''cumulative_probability='',1000f9.4)') cumulative_probability(1:n_connected_dets)
  if(ipr.ge.1) write(6,'(''connected_dets_sign='',1000i3)') connected_dets_sign(1:n_connected_dets)

  if(nwalk_new.eq.0) then
    walk_wt(iwalk)=0
   else
    random=rannyu()
!   random=mod(random,1._rk/nwalk_new) ! This is a random number between 0 and 1/nwalk_new
    nwalk_all_det=0
    do i_nonzero=1,n_connected_dets
!     nwalk_this_det=int((cumulative_probability(i_nonzero)/cumulative_probability(n_connected_dets)-random)*nwalk_new+1)-nwalk_all_det
      nwalk_this_det=int((cumulative_probability(i_nonzero)/cumulative_probability(n_connected_dets))*nwalk_new+random)-nwalk_all_det
      if(nwalk_this_det.gt.0) then
        nwalk_all_det=nwalk_all_det+nwalk_this_det
        if(connected_dets(i_nonzero).eq.walk_det(iwalk)) then ! diagonal move
          walk_wt(iwalk)=nwalk_this_det*connected_dets_sign(i_nonzero)*walk_wt_child
          if(ipr.ge.1) write(6,'(''iwalk,walk_wt(iwalk)='',i8,f10.6)') iwalk,walk_wt(iwalk)
        else ! off-diagonal move
          nwalk=nwalk+1
          snd_cnt(1)=snd_cnt(1)+1
          if(nwalk.gt.MWALK) then
            write(6,'(''stop nwalk>MWALK, e_trial, reweight_factor=''9d12.4)') e_trial, 1._rk/reweight_factor_inv
            call flush(6)
            stop 'nwalk>MWALK'
          endif
          walk_det(nwalk)=connected_dets(i_nonzero)
          walk_wt(nwalk)=nwalk_this_det*connected_dets_sign(i_nonzero)*walk_wt_child
          if(initiator(iwalk).ge.2) then ! If it is spawned by an initiator then 1 otherwise 0
            initiator(nwalk)=1
          else
            initiator(nwalk)=0
          endif
!Warning: Is this in the right place for a parallel run?
          e_num_walker(nwalk)=1.e51_rk
          if(ipr.ge.1) write(6,'(''nwalk,walk_det(nwalk),walk_wt(nwalk)='',2i8,f9.1)') nwalk,walk_det(nwalk),walk_wt(nwalk)
        endif ! end of diagonal and off-diagonal moves
      else
        if(connected_dets(i_nonzero).eq.walk_det(iwalk)) walk_wt(iwalk)=0 ! If nwalk_this_det=0 and diagonal move
      endif
    enddo
    if(nwalk_all_det.ne.nwalk_new) then
      write(6,'(''Warning: iwalk, nwalk_all_det, nwalk_new='',9i5)') iwalk, nwalk_all_det, nwalk_new
      stop 'In move_heat_bath2 nwalk_all_det.ne.nwalk_new'
    endif
  endif
  !write(6,'(''iwalk,nwalk,nwalk_all_det,nwalk_new,walk_det(nwalk-nwalk_all_det+1:nwalk)'',104i5)') iwalk,nwalk,nwalk_all_det,nwalk_new,walk_det(nwalk-nwalk_new+1:nwalk)
  !write(6,'(''iwalk,nwalk,nwalk_all_det,nwalk_new'',104i5)') iwalk,nwalk,nwalk_all_det,nwalk_new

  end subroutine move_heat_bath2
!-------------------------------------------------------------------------------------

  subroutine move_heat_bath3(iwalk)
! ==============================================================================
! Description   : Do move using heat bath on the absolute value of the elements of the appropriate column of the propagator
!               : Version for realistic problems where Hamiltonian cannot be stored
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use hamiltonian_mod
  use heg
  use chemistry, only : find_connected_dets_chem,norb
  use hubbard

  implicit none
  integer, intent(in) :: iwalk
  integer             :: i, nwalk_new
  integer             :: nwalk_this_det, nwalk_all_det, nnew
  real(rk)            :: random, rannyu, probability, walk_wt_child
 !real(rk)            :: probability_fn
  real(rk)            ::  norm

  logical             :: accepted

  accepted=.true.
  norm=0._rk

  if(hamiltonian_type.eq.'heg') then
    call find_connected_dets_heg(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements)
  elseif(hamiltonian_type.eq.'chem') then
    ! need to fix this if heat bath is to be used; currently, only the first n_connected_dets elements of connected_dets_up,dn and matrix_elements are used.
    call find_connected_dets_chem(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, norb, connected_matrix_elements)
  elseif((hamiltonian_type.eq.'hubbard2') .and. (importance_sampling .eq.0)) then
    call find_connected_dets_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements)
  ! Hubbard2 Importance sampling by HJC on Apr 8 2011
  elseif((hamiltonian_type.eq.'hubbard2') .and. (importance_sampling .eq.1)) then
      if (main_wf_type .eq. 'gutz') then
          if (run_type .eq. 'fixed_node1') then
              call find_connected_dets_imp_gutz_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, 1)
          elseif (run_type .eq. 'fixed_node2') then
              call find_connected_dets_imp_gutz_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, 2)
          elseif (run_type .eq. 'partial_node') then
              call find_connected_dets_imp_gutz_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, 3)
          elseif (run_type .eq. 'fixed_node4') then
              call find_connected_dets_imp_gutz_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, 4)
          elseif (run_type .eq. 'sr') then
              call find_connected_dets_imp_gutz_hubbard_sr(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, connected_matrix_elements_fn)
          else
              call find_connected_dets_imp_gutz_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements)
          endif
      else
          stop 'Importance sampling only implemented for gutzwiller on Hubbard'
      endif
  elseif(hamiltonian_type.eq.'hubbardk') then
    ! need to fix this if heat bath is to be used; currently, only the first n_connected_dets elements of connected_dets_up,dn and connected_matrix_elements are used.
    call find_connected_dets_hubbard_k(walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets, connected_dets_up, connected_dets_dn, nsites, connected_matrix_elements)
  else
    stop 'move_heat_bath3: hamiltonian_type must be one of heg, chem, hubbard2, hubbardk'
  endif

if(ipr.ge.4) write(6,'(/,''e_trial='',9f10.5)') e_trial
if(ipr.ge.4) write(6,'(''iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets='',9i6)') iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), n_connected_dets
if(ipr.ge.4) write(6,'(''connected_dets_up, connected_dets_dn, connected_matrix_elements'',/,(2i6,d12.4))') (connected_dets_up(i), connected_dets_dn(i), connected_matrix_elements(i),i=1,n_connected_dets)

! Set cumulative probabilities
  do i=1,n_connected_dets
    if(walk_dets_up(iwalk).eq.connected_dets_up(i) .and. walk_dets_dn(iwalk).eq.connected_dets_dn(i)) then
      probability=(1+tau*(e_trial-connected_matrix_elements(i)))
      !if (run_type .eq. 'sr') probability_fn=(1+tau*(e_trial-connected_matrix_elements_fn(i)))
      !write(6,'(''diagonal connected: i, probability='',i5,es12.4)') i, probability
    else
      if(importance_sampling.eq.0) then
        probability=-tau*connected_matrix_elements(i)
      ! Hubbard2 Importance sampling by HJC April 8 2011
      elseif((importance_sampling.eq.1) .and. (hamiltonian_type .eq. 'hubbard2')) then
        probability=-tau*connected_matrix_elements(i)
        !if (run_type .eq. 'sr') probability_fn=-tau*connected_matrix_elements_fn(i)
        if (run_type .eq. 'fixed_node3' .or. run_type .eq. 'vmc') then
            if (probability .lt. 0._rk) probability=0._rk
        endif
      elseif((importance_sampling.eq.1) .and. (hamiltonian_type .ne. 'hubbard2')) then  ! Guiding wavefunctions are not being used for any other Hamiltonian - so this part never accessed
        probability=-tau*connected_matrix_elements(i) * psi_g(i)/psi_g(walk_det(iwalk))
      else
        stop 'importance_sampling must be 0 or 1'
      endif
    endif
    if (run_type .eq. 'vmc') then
      norm=norm+probability
    endif

    !if (run_type .eq. 'vmc') then
    !  probability=reweight_factor_inv*walk_wt(iwalk)*probability
    !else
    probability=walk_wt(iwalk)*probability           ! removed reweight factor as  it is being accounted for in subroutine walk
    !if (run_type .eq. 'sr') probability_fn=reweight_factor_inv*walk_wt(iwalk)*probability_fn
    !endif

    connected_dets_sign(i)=int(sign(1._rk,probability))
    if(i.eq.1) then
      cumulative_probability(i)=abs(probability)
      !if (run_type .eq. 'sr') cumulative_probability_fn(i)=abs(probability_fn)
    else
      cumulative_probability(i)=cumulative_probability(i-1)+abs(probability)
      !if (run_type .eq. 'sr') cumulative_probability_fn(i)=cumulative_probability_fn(i-1)+abs(probability_fn)
    endif
  enddo

  if (run_type .eq. 'vmc') then
    do i=1,n_connected_dets
      cumulative_probability(i)=cumulative_probability(i)/norm
    enddo
  endif

! nwalk_new=int(cumulative_probability(n_connected_dets)+rannyu()) ! Total number of progeny from this walker
  nwalk_new=max(nint(cumulative_probability(n_connected_dets)),1) ! Total number of progeny from this walker
  !write(6,'(''nwalk_new, n_connected_dets, cumulative_probability'',2i12,es12.4)') nwalk_new, n_connected_dets, cumulative_probability(n_connected_dets)
  walk_wt_child=cumulative_probability(n_connected_dets)/nwalk_new
  !write(6,'(''cumulative_probability(i)='',1000es12.4)') (cumulative_probability(i), i=1,n_connected_dets)

  ! For VMC number of progeny must equal walk_wt(iwalk) because no NET walkers are created or destroyed - comment by HJC circa April 2011

  if(ipr.ge.1) then
    if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
      write(6,'(/,''iwalk, walk_det(iwalk), walk_wt(iwalk), n_connected_dets, nwalk_new='',2i8,f10.4,2i8)') &
 &    iwalk, walk_det(iwalk), walk_wt(iwalk), n_connected_dets, nwalk_new
    else
      write(6,'(/,''iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk), n_connected_dets, nwalk_new='',3i8,f10.4,2i8)') &
 &    iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk), n_connected_dets, nwalk_new
    endif
    write(6,'(''reweight_factor'',f10.4)') 1._rk/reweight_factor_inv
    write(6,'(''cumulative_probability='',1000f14.4)') cumulative_probability(1:n_connected_dets)
    write(6,'(''connected_dets_sign='',1000i3)') connected_dets_sign(1:n_connected_dets)
  endif
  call flush(6)

  if(nwalk_new.eq.0) then
    walk_wt(iwalk)=0
  else
    random=rannyu()
!   random=mod(random,1._rk/nwalk_new) ! This is a random number between 0 and 1/nwalk_new
    nwalk_all_det=0
    nnew=0
    do i=1,n_connected_dets
      if(cumulative_probability(n_connected_dets).ne.0._rk) then
        nwalk_this_det=int((cumulative_probability(i)/cumulative_probability(n_connected_dets))*nwalk_new+random)-nwalk_all_det
      else
        nwalk_this_det=0
      endif
      if(nwalk_this_det.eq.0) then                                                                                         ! If # of walkers on the current det should be 0
        if(walk_dets_up(iwalk).eq.connected_dets_up(i) .and. walk_dets_dn(iwalk).eq.connected_dets_dn(i)) then
           walk_wt(iwalk)=0 ! If diag. move and # of walkers on the current det should be 0, set wt=0
        endif
      else
        nwalk_all_det=nwalk_all_det+nwalk_this_det
        if(walk_dets_up(iwalk).eq.connected_dets_up(i) .and. walk_dets_dn(iwalk).eq.connected_dets_dn(i)) then  ! diagonal move
          walk_wt(iwalk)=nwalk_this_det*connected_dets_sign(i)*walk_wt_child
!         nnew=nnew+abs(walk_wt(iwalk))
          nnew=nnew+nwalk_this_det
          if(ipr.ge.1) write(6,'(''iwalk, nnew (at diag), walk_wt(iwalk),='',i8,f10.6,i8)') iwalk, walk_wt(iwalk), nnew
        else                                                                                                    ! off-diag move
          nwalk=nwalk+1
          snd_cnt(1)=snd_cnt(1)+1
          if(nwalk.gt.MWALK) then
            write(6,'(''stop nwalk>MWALK, e_trial, reweight_factor=''9d12.4)') e_trial, 1._rk/reweight_factor_inv
            call flush(6)
            stop 'nwalk>MWALK'
          endif
          walk_dets_up(nwalk)=connected_dets_up(i)
          walk_dets_dn(nwalk)=connected_dets_dn(i)
          walk_wt(nwalk)=nwalk_this_det*connected_dets_sign(i)*walk_wt_child

          if(initiator(iwalk).ge.2) then ! If it is spawned by an initiator then 1 otherwise 0
            initiator(nwalk)=1
          else
            initiator(nwalk)=0
          endif
!Warning: Is this in the right place for a parallel run?
          e_num_walker(nwalk)=1.e51_rk
          e_den_walker(nwalk)=1.e51_rk
          if(ipr.ge.1) write(6,'(''nwalk,walk_dets_up(nwalk),walk_dets_dn(nwalk),walk_wt(nwalk)='',3i8,f10.4)') nwalk,walk_dets_up(nwalk),walk_dets_dn(nwalk),walk_wt(nwalk)
!         nnew=nnew+abs(walk_wt(nwalk))
          nnew=nnew+nwalk_this_det
          ! Comment by HJC nnew=nnew+abs(walk_wt(nwalk)) Moved up
          if(ipr.ge.1) write(6,'(''nnew (at off diag)'',9i8)') nnew
        endif ! diagonal/off-diag
      endif   ! nwalk_this_det
    enddo

!!Warning: Temporarily turn this off because with SR walk_wt can be zero (I should remove those walkers but I am not at present) and so nnew=0 but nwalk_new=1
!   if (run_type .ne. 'fixed_node3') then                         ! Bypass this check for fixed node 3 - HJC May 26 2011
!       if(nnew.ne.nwalk_new) then
!           write(6,'(''Warning: iwalk, nnew, nwalk_new='',3i20)') iwalk, nnew, nwalk_new
!           stop 'nnew .ne. nwalk_new'
!       endif
!   endif

  endif
  call flush(6)

!  if (allocated(connected_dets_sign_real)) deallocate(connected_dets_sign_real)

  end subroutine move_heat_bath3
!------------------------------------------------------------------------------


  subroutine shell_sort(walk_det,walk_wt,nwalk)
! ==============================================================================
! Description   : Shell-Metzger sort in ascending order.
! Not presently being used.
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  implicit none
  integer, intent(inout) :: walk_det(nwalk)
  integer, intent(in) :: nwalk
  real(rk), intent(inout) :: walk_wt(nwalk)
  integer it,nn,lognb2,m,k,j,i,l
  real(rk) rt
  lognb2=int(log(real(nwalk,rk))/log(2._rk)+10._rk**(-14))
  m=nwalk
  do 20 nn=1,lognb2
    m=m/2
    k=nwalk-m
    do 20 j=1,k
      do 10 i=j,1,-m
        l=i+m
        if (abs(walk_det(l)).gt.abs(walk_det(i))) goto 20
        it=walk_det(i)
        walk_det(i)=walk_det(l)
        walk_det(l)=it
        rt=walk_wt(i)
        walk_wt(i)=walk_wt(l)
        walk_wt(l)=rt
10    continue
20  continue
  return
  end subroutine shell_sort
!-------------------------------------------------------------------------------------

  subroutine shell_sort2(walk_dets_up,walk_dets_dn,walk_wt,nwalk)
! ==============================================================================
! Description   : Shell-Metzger sort in ascending order.
! Not presently being used.
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(in) :: nwalk
  real(rk), intent(inout) :: walk_wt(nwalk)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
  type(ik_vec) it_ik
#else
  integer(ik), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
  integer(ik) it_ik
#endif
  integer nn,lognb2,m,k,j,i,l
  real(rk) rt
  lognb2=int(log(real(nwalk,rk))/log(2._rk)+10._rk**(-14))
  m=nwalk
  do 20 nn=1,lognb2
    m=m/2
    k=nwalk-m
    do 20 j=1,k
      do 10 i=j,1,-m
        l=i+m
        if (abs(walk_dets_up(l)) .gt. abs(walk_dets_up(i))) then
            goto 20
        elseif ((walk_dets_up(l) .eq. walk_dets_up(i)) .and. (abs(walk_dets_dn(l)) .gt. abs(walk_dets_dn(i)))) then
            goto 20
        endif
        it_ik=walk_dets_up(i)
        walk_dets_up(i)=walk_dets_up(l)
        walk_dets_up(l)=it_ik
        it_ik=walk_dets_dn(i)
        walk_dets_dn(i)=walk_dets_dn(l)
        walk_dets_dn(l)=it_ik
        rt=walk_wt(i)
        walk_wt(i)=walk_wt(l)
        walk_wt(l)=rt
10    continue
20  continue
  return
  end subroutine shell_sort2
!-------------------------------------------------------------------------------

!  subroutine shell_sort2(walk_dets_up,walk_dets_dn,walk_wt,nwalk)
!! ==============================================================================
!! Description   : Shell-Metzger sort in ascending order.
!! Author        : Cyrus Umrigar
!! ------------------------------------------------------------------------------
!  use types, only: ik
!  implicit none
!  type(ik_vec), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
!  integer, intent(in) :: nwalk
!  real(rk), intent(inout) :: walk_wt(nwalk)
!  integer nn,lognb2,m,k,j,i,l
!  type(ik_vec) it_ik
!  real(rk) rt
!  type(ik_vec), parameter :: two_to_63=9223372036854775808_ik
!! type(ik_vec), parameter :: two_to_31=2147483648_ik
!  lognb2=int(log(real(nwalk,rk))/log(2._rk)+10._rk**(-14))
!  m=nwalk
!  do 20 nn=1,lognb2
!    m=m/2
!    k=nwalk-m
!    do 20 j=1,k
!      do 10 i=j,1,-m
!        l=i+m
!        if (two_to_63*walk_dets_up(l)+walk_dets_dn(l) .gt. two_to_63*walk_dets_up(i)+walk_dets_dn(i)) goto 20
!!       if (two_to_31*walk_dets_up(l)+walk_dets_dn(l) .gt. two_to_31*walk_dets_up(i)+walk_dets_dn(i)) goto 20
!        it_ik=walk_dets_up(i)
!        walk_dets_up(i)=walk_dets_up(l)
!        walk_dets_up(l)=it_ik
!        it_ik=walk_dets_dn(i)
!        walk_dets_dn(i)=walk_dets_dn(l)
!        walk_dets_dn(l)=it_ik
!        rt=walk_wt(i)
!        walk_wt(i)=walk_wt(l)
!        walk_wt(l)=rt
!10    continue
!20  continue
!  return
!  end subroutine shell_sort2
!-------------------------------------------------------------------------------

  subroutine shell_sort_realwt(walk_det,walk_wt,nwalk)
! ==============================================================================
! Description   : Shell-Metzger sort in ascending order.
! Author        : Cyrus Umrigar/ HJC - just change walk_wt to real
!                 Will be merged eventually when we use real weights.
!                 This is a temporary hack
! Date           : Oct 21 2011
! ------------------------------------------------------------------------------
  implicit none
  integer, intent(inout) :: walk_det(nwalk)
  real(rk), intent(inout):: walk_wt(nwalk)
  integer, intent(in) :: nwalk
  integer  it,nn,lognb2,m,k,j,i,l
  real(rk) itr

  lognb2=int(log(real(nwalk,rk))/log(2._rk)+10._rk**(-14))
  m=nwalk
  do 20 nn=1,lognb2
    m=m/2
    k=nwalk-m
    do 20 j=1,k
      do 10 i=j,1,-m
        l=i+m
        if (abs(walk_det(l)).gt.abs(walk_det(i))) goto 20
        it=walk_det(i)
        walk_det(i)=walk_det(l)
        walk_det(l)=it
        itr=walk_wt(i)
        walk_wt(i)=walk_wt(l)
        walk_wt(l)=itr
10    continue
20  continue
  return
  end subroutine shell_sort_realwt
!-------------------------------------------------------------------------------------

  subroutine shell_sort2_realwt(walk_dets_up,walk_dets_dn,walk_wt,nwalk)
! ==============================================================================
! Description   : Shell-Metzger sort in ascending order.
! Author        : Cyrus Umrigar/ HJC - just change walk_wt to real
!                 Will be merged eventually when we use real weights.
!                 This is a temporary hack
! Date           : Oct 21 2011
! ------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  real(rk), intent(inout) :: walk_wt(nwalk)
  integer, intent(in) :: nwalk
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
  type(ik_vec) it_ik
#else
  integer(ik), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
  integer(ik) it_ik
#endif
  integer nn,lognb2,m,k,j,i,l
  real(rk) :: itr
  lognb2=int(log(real(nwalk,rk))/log(2._rk)+10._rk**(-14))
  m=nwalk
  do 20 nn=1,lognb2
    m=m/2
    k=nwalk-m
    do 20 j=1,k
      do 10 i=j,1,-m
        l=i+m
        if (abs(walk_dets_up(l)) .gt. abs(walk_dets_up(i))) then
            goto 20
        elseif ((walk_dets_up(l) .eq. walk_dets_up(i)) .and. (abs(walk_dets_dn(l)) .gt. abs(walk_dets_dn(i)))) then
            goto 20
        endif
        it_ik=walk_dets_up(i)
        walk_dets_up(i)=walk_dets_up(l)
        walk_dets_up(l)=it_ik
        it_ik=walk_dets_dn(i)
        walk_dets_dn(i)=walk_dets_dn(l)
        walk_dets_dn(l)=it_ik
        itr=walk_wt(i)
        walk_wt(i)=walk_wt(l)
        walk_wt(l)=itr
10    continue
20  continue
  return
  end subroutine shell_sort2_realwt
!-------------------------------------------------------------------------------

!  subroutine shell_sort2_realwt(walk_dets_up,walk_dets_dn,walk_wt,nwalk)
!! ==============================================================================
!! Description   : Shell-Metzger sort in ascending order.
!! Author        : Cyrus Umrigar/ HJC - just change walk_wt to real
!!                 Will be merged eventually when we use real weights.
!!                 This is a temporary hack
!! Date           : Oct 21 2011
!! ------------------------------------------------------------------------------
!  use types, only: ik
!  implicit none
!  type(ik_vec), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
!  real(rk), intent(inout) :: walk_wt(nwalk)
!  integer, intent(in) :: nwalk
!  integer nn,lognb2,m,k,j,i,l
!  type(ik_vec) it_ik
!  type(ik_vec), parameter :: two_to_63=9223372036854775808_ik
!  real(rk) :: itr
!! type(ik_vec), parameter :: two_to_31=2147483648_ik
!  lognb2=int(log(real(nwalk,rk))/log(2._rk)+10._rk**(-14))
!  m=nwalk
!  do 20 nn=1,lognb2
!    m=m/2
!    k=nwalk-m
!    do 20 j=1,k
!      do 10 i=j,1,-m
!        l=i+m
!        if (two_to_63*walk_dets_up(l)+walk_dets_dn(l) .gt. two_to_63*walk_dets_up(i)+walk_dets_dn(i)) goto 20
!!       if (two_to_31*walk_dets_up(l)+walk_dets_dn(l) .gt. two_to_31*walk_dets_up(i)+walk_dets_dn(i)) goto 20
!        it_ik=walk_dets_up(i)
!        walk_dets_up(i)=walk_dets_up(l)
!        walk_dets_up(l)=it_ik
!        it_ik=walk_dets_dn(i)
!        walk_dets_dn(i)=walk_dets_dn(l)
!        walk_dets_dn(l)=it_ik
!        itr=walk_wt(i)
!        walk_wt(i)=walk_wt(l)
!        walk_wt(l)=itr
!10    continue
!20  continue
!  return
!  end subroutine shell_sort2_realwt
!!-------------------------------------------------------------------------------

  subroutine sort_walkers(nwalk)
! ==============================================================================
! Description   : Call merge_sort to sort with key walk_det.
!               : Then use iorder to sort the auxilliary items, viz., walk_wt and initiator.
! Written by    : Cyrus Umrigar
! ------------------------------------------------------------------------------

  use types, only: ik
  use common_ham, only : hamiltonian_type
  implicit none
  integer, intent(in) :: nwalk
  integer i

  do i=1,nwalk
    iorder(i)=i
  enddo
  call merge_sort(walk_det, iorder, nwalk, temp_i_1, temp_i_2)
  walk_wt(1:nwalk)=walk_wt(iorder(1:nwalk))
  initiator(1:nwalk)=initiator(iorder(1:nwalk))
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(1:nwalk)=e_num_walker(iorder(1:nwalk))
    e_den_walker(1:nwalk)=e_den_walker(iorder(1:nwalk))
  endif
  matrix_elements(1:nwalk) = matrix_elements(iorder(1:nwalk))

  return
  end subroutine sort_walkers
!-------------------------------------------------------------------------------

  subroutine sort_walkers2
! ==============================================================================
! Description   : Call merge_sort2 to sort with primary key walk_dets_up and secondary key walk_dets_dn.
!               : Then use iorder to sort the auxilliary items, viz., walk_wt and initiator.
! Written by    : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use types, only: ik
  use common_ham, only : hamiltonian_type
  implicit none
! type(ik_vec), allocatable ::walk_dets_up(:), walk_dets_dn(:), temp_i16(:)
! integer, allocatable ::  walk_wt(:), iorder(:), temp_i_2(:)
  integer iwalk, jwalk, i, n_same_det_up

! read(5,*) (walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk), iwalk=1,nwalk)
! write(6,'(//,(3i6))') (walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk), iwalk=1,nwalk)
  call flush(6)

  do i=1,nwalk
    iorder(i)=i
  enddo
  call merge_sort2(walk_dets_up, iorder, nwalk, temp_i16, temp_i_2)
  walk_dets_dn(1:nwalk)=walk_dets_dn(iorder(1:nwalk))
  walk_wt(1:nwalk)=walk_wt(iorder(1:nwalk))
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(1:nwalk)=e_num_walker(iorder(1:nwalk))
    e_den_walker(1:nwalk)=e_den_walker(iorder(1:nwalk))
  endif
  initiator(1:nwalk)=initiator(iorder(1:nwalk))
  ! Edited by AAH on 8 Feb 2012.
  imp_distance(1:nwalk)=imp_distance(iorder(1:nwalk))
  ! End of edit
  matrix_elements(1:nwalk)=matrix_elements(iorder(1:nwalk))

  iwalk=1
  do while (iwalk <= nwalk)
    n_same_det_up=1
    do jwalk=iwalk+1,nwalk
      if(walk_dets_up(jwalk) == walk_dets_up(iwalk)) then                                  ! same det
        n_same_det_up=n_same_det_up+1
      endif
      if(walk_dets_up(jwalk) /= walk_dets_up(iwalk) .or. iwalk+n_same_det_up > nwalk) then ! either not same det or at end of list
        if(n_same_det_up >= 2) then
          do i=iwalk,iwalk+n_same_det_up-1
            iorder(i)=i
          enddo
          call merge_sort2(walk_dets_dn(iwalk), iorder(iwalk), n_same_det_up, temp_i16, temp_i_2)
          walk_dets_up(iwalk:iwalk+n_same_det_up-1)=walk_dets_up(iorder(iwalk:iwalk+n_same_det_up-1))
          walk_wt(iwalk:iwalk+n_same_det_up-1)=walk_wt(iorder(iwalk:iwalk+n_same_det_up-1))
          initiator(iwalk:iwalk+n_same_det_up-1)=initiator(iorder(iwalk:iwalk+n_same_det_up-1))
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk:iwalk+n_same_det_up-1)=e_num_walker(iorder(iwalk:iwalk+n_same_det_up-1))
            e_den_walker(iwalk:iwalk+n_same_det_up-1)=e_den_walker(iorder(iwalk:iwalk+n_same_det_up-1))
          endif
          ! Edited by AAH on 8 Feb 2012.
          imp_distance(iwalk:iwalk+n_same_det_up-1)=imp_distance(iorder(iwalk:iwalk+n_same_det_up-1))
          ! End of edit
          matrix_elements(iwalk:iwalk+n_same_det_up-1)=matrix_elements(iorder(iwalk:iwalk+n_same_det_up-1))
          exit ! jwalk
        endif
      endif
    enddo ! jwalk
    iwalk=iwalk+n_same_det_up
  enddo ! do while

! write(6,'(/,(4i6))') (walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk), iorder(iwalk), iwalk=1,nwalk)

  return
  end subroutine sort_walkers2
!-------------------------------------------------------------------------------

  subroutine sort_walkers3
! ==============================================================================
! Description   : Call merge_sort2 to sort with primary key walk_dets_up and secondary key walk_dets_dn.
!               : Then use iorder to sort the auxilliary items, viz., walk_wt and initiator.
! Written by    : Cyrus Umrigar
! Modified      : A Holmes, 28 Jan 2013. Sort walkers first by whether or not imp_distance<=0, then by label.
!               : Assumes that the first ndet_psi_t_connected dets are already sorted, as well as the next ndet_outside_ct dets.
! ------------------------------------------------------------------------------
  use types, only: ik
  use common_ham, only : hamiltonian_type
  use common_run, only : ndet_outside_ct
  implicit none
  integer iwalk, jwalk, i, n_same_det_up
  integer :: nwalk_eff

  if (nwalk<ndet_psi_t_connected+ndet_outside_ct) then
    write (6,*) "stopping because nwalk<ndet_psi_t_connected+ndet_outside_ct; nwalk=",nwalk,", ndet_psi_t_connected=",ndet_psi_t_connected,", ndet_outside_ct=",ndet_outside_ct
    call flush(6)
    stop 'in sort_walkers3: nwalk<ndet_psi_t_connected+ndet_outside_ct'
  endif
  if (nwalk==ndet_psi_t_connected+ndet_outside_ct)  return

  ! Sort by up det
  do i=ndet_psi_t_connected+ndet_outside_ct+1,nwalk
    iorder(i)=i
  enddo
  call merge_sort2(walk_dets_up(ndet_psi_t_connected+ndet_outside_ct+1:nwalk), iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk), nwalk_eff, temp_i16, temp_i_2)
  walk_dets_dn(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=walk_dets_dn(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
  walk_wt(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=walk_wt(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=e_num_walker(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
    e_den_walker(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=e_den_walker(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
  endif
  initiator(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=initiator(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
  imp_distance(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=imp_distance(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
  matrix_elements(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=matrix_elements(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))

  ! Within each up det, sort by dn det
  iwalk=ndet_psi_t_connected+ndet_outside_ct+1
  do while (iwalk <= nwalk)
    n_same_det_up=1
    do jwalk=iwalk+1,nwalk
      if(walk_dets_up(jwalk) == walk_dets_up(iwalk)) then                                  ! same det
        n_same_det_up=n_same_det_up+1
      endif
      if(walk_dets_up(jwalk) /= walk_dets_up(iwalk) .or. iwalk+n_same_det_up > nwalk) then ! either not same det or at end of list
        if(n_same_det_up >= 2) then
          do i=iwalk,iwalk+n_same_det_up-1
            iorder(i)=i
          enddo
          call merge_sort2(walk_dets_dn(iwalk), iorder(iwalk), n_same_det_up, temp_i16, temp_i_2)
          walk_dets_up(iwalk:iwalk+n_same_det_up-1)=walk_dets_up(iorder(iwalk:iwalk+n_same_det_up-1))
          walk_wt(iwalk:iwalk+n_same_det_up-1)=walk_wt(iorder(iwalk:iwalk+n_same_det_up-1))
          initiator(iwalk:iwalk+n_same_det_up-1)=initiator(iorder(iwalk:iwalk+n_same_det_up-1))
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk:iwalk+n_same_det_up-1)=e_num_walker(iorder(iwalk:iwalk+n_same_det_up-1))
            e_den_walker(iwalk:iwalk+n_same_det_up-1)=e_den_walker(iorder(iwalk:iwalk+n_same_det_up-1))
          endif
          ! Edited by AAH on 8 Feb 2012.
          imp_distance(iwalk:iwalk+n_same_det_up-1)=imp_distance(iorder(iwalk:iwalk+n_same_det_up-1))
          ! End of edit
          matrix_elements(iwalk:iwalk+n_same_det_up-1)=matrix_elements(iorder(iwalk:iwalk+n_same_det_up-1))
          exit ! jwalk
        endif
      endif
    enddo ! jwalk
    iwalk=iwalk+n_same_det_up
  enddo ! do while

  return
  end subroutine sort_walkers3
!-------------------------------------------------------------------------------

  subroutine sort_my_walkers2_up_dn(my_nwalk)
! ==============================================================================
! Description   : Call merge_sort2_up_dn to sort with primary key walk_dets_up and secondary key walk_dets_dn.
!               : Then use iorder to sort the auxilliary items, viz., walk_wt and initiator.
! Written by    : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use types, only: ik
  use common_ham, only : hamiltonian_type
  implicit none
  integer, intent(in) :: my_nwalk
  integer i

  do i=1,my_nwalk
    iorder(i)=i
  enddo
  call merge_sort2_up_dn(walk_dets_up, walk_dets_dn, iorder, my_nwalk, temp_i16_up, temp_i16_dn, temp_i_2)
  walk_wt(1:my_nwalk)=walk_wt(iorder(1:my_nwalk))
  initiator(1:my_nwalk)=initiator(iorder(1:my_nwalk))
  imp_distance(1:my_nwalk)=imp_distance(iorder(1:my_nwalk))
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(1:my_nwalk)=e_num_walker(iorder(1:my_nwalk))
    e_den_walker(1:my_nwalk)=e_den_walker(iorder(1:my_nwalk))
  endif
  matrix_elements(1:my_nwalk)=matrix_elements(iorder(1:my_nwalk))

! write(6,'(/,(4i6))') (walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk), iorder(iwalk), iwalk=1,my_nwalk)

  return
  end subroutine sort_my_walkers2_up_dn
!-------------------------------------------------------------------------------

  subroutine sort_walkers3_up_dn
! ==============================================================================
! Description   : Call merge_sort2_up_dn to sort with primary key walk_dets_up and secondary key walk_dets_dn.
!               : Then use iorder to sort the auxilliary items, viz., walk_wt and initiator.
! Written by    : Cyrus Umrigar
! Modified      : A Holmes, 29 Jan 2013. Only sort newly spawned dets.
! No longer used.  Replaced by sort_my_walkers3_up_dn
! ------------------------------------------------------------------------------
  use types, only: ik
  use common_ham, only : hamiltonian_type
  use common_run, only : ndet_outside_ct
  implicit none
  integer i,nwalk_eff

!***b
!  if (my_nwalk<my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)) then
!    write (6,*) "nwalk=",nwalk,"ndet_psi_t_connected=",ndet_psi_t_connected,"ndet_outside_ct=",ndet_outside_ct
!    call flush(6)
!    stop 'in sort_walkers3_up_dn [core=",whoami,"]: nwalk<ndet_psi_t_connected+ndet_outside_ct'
!  endif
!  if (my_nwalk==my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1))  return
!
!  nwalk_eff = my_nwalk-my_ndet_psi_t_connected(whoami+1)-my_ndet_outside_ct(whoami+1)
!
!  do i=my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1,my_nwalk
!    iorder(i)=i
!  enddo
!  call merge_sort2_up_dn(walk_dets_up(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk), walk_dets_dn(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk), iorder(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk), nwalk_eff, temp_i16_up, temp_i16_dn, temp_i_2)
!  walk_wt(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk)=walk_wt(iorder(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk))
!  initiator(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk)=initiator(iorder(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk))
!  imp_distance(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk)=imp_distance(iorder(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk))
!  if (hamiltonian_type.eq.'hubbard2'.or.run_type.eq.'sr') then!***WARNING check compatibility [TODO]
!    e_num_walker(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk)=e_num_walker(iorder(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk))
!    e_den_walker(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk)=e_den_walker(iorder(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk))
!  endif
!  matrix_elements(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk)=matrix_elements(iorder(my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1:my_nwalk))
!  my_nwalk = nwalk_eff + my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)
!-------
  if (nwalk<ndet_psi_t_connected+ndet_outside_ct) then
    write (6,*) "nwalk=",nwalk,"ndet_psi_t_connected=",ndet_psi_t_connected,"ndet_outside_ct=",ndet_outside_ct
    call flush(6)
    stop 'in sort_walkers3_up_dn: nwalk<ndet_psi_t_connected+ndet_outside_ct'
  endif
  if (nwalk==ndet_psi_t_connected+ndet_outside_ct)  return

  nwalk_eff = nwalk-ndet_psi_t_connected-ndet_outside_ct

  do i=ndet_psi_t_connected+ndet_outside_ct+1,nwalk
    iorder(i)=i
  enddo
  call merge_sort2_up_dn(walk_dets_up(ndet_psi_t_connected+ndet_outside_ct+1:nwalk), walk_dets_dn(ndet_psi_t_connected+ndet_outside_ct+1:nwalk), iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk), nwalk_eff, temp_i16_up, temp_i16_dn, temp_i_2)
  walk_wt(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=walk_wt(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
  initiator(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=initiator(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
  imp_distance(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=imp_distance(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=e_num_walker(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
    e_den_walker(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=e_den_walker(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))
  endif
  matrix_elements(ndet_psi_t_connected+ndet_outside_ct+1:nwalk)=matrix_elements(iorder(ndet_psi_t_connected+ndet_outside_ct+1:nwalk))

! write(6,'(/,(4i6))') (walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk), iorder(iwalk), iwalk=1,nwalk)

  nwalk = nwalk_eff + ndet_psi_t_connected + ndet_outside_ct
!***e
  return
  end subroutine sort_walkers3_up_dn
!-------------------------------------------------------------------------------

  subroutine sort_my_walkers3_up_dn(my_nwalk,my_ndet_psi_t_connected,my_ndet_outside_ct)!****REMOVE WHEN DONE
! ==============================================================================
! Description   : Call merge_sort2_up_dn to sort with primary key walk_dets_up and secondary key walk_dets_dn.
!               : Then use iorder to sort the auxilliary items, viz., walk_wt and initiator.
! Written by    : Cyrus Umrigar
! Modified      : A Holmes, 29 Jan 2013. Only sort newly spawned dets.
! ------------------------------------------------------------------------------
  use types, only: ik
  use common_ham, only : hamiltonian_type
  use common_run, only : ndet_outside_ct
  implicit none
  integer,intent(inout) :: my_nwalk,my_ndet_outside_ct
  integer,intent(in)    :: my_ndet_psi_t_connected
  integer i,nwalk_eff

!***b
  if (my_nwalk<my_ndet_psi_t_connected+my_ndet_outside_ct) then
    write (6,*) "nwalk=",nwalk,"ndet_psi_t_connected=",ndet_psi_t_connected,"ndet_outside_ct=",ndet_outside_ct
    call flush(6)
    stop 'in sort_my_walkers3_up_dn [core=",whoami,"]: nwalk<ndet_psi_t_connected+ndet_outside_ct'
  endif
  if (my_nwalk==my_ndet_psi_t_connected+my_ndet_outside_ct)  return

  nwalk_eff = my_nwalk-my_ndet_psi_t_connected-my_ndet_outside_ct

  do i=my_ndet_psi_t_connected+my_ndet_outside_ct+1,my_nwalk
    iorder(i)=i
  enddo
  call merge_sort2_up_dn(walk_dets_up(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk), walk_dets_dn(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk), iorder(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk), nwalk_eff, temp_i16_up, temp_i16_dn, temp_i_2)
  walk_wt(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk)=walk_wt(iorder(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk))
  initiator(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk)=initiator(iorder(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk))
  imp_distance(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk)=imp_distance(iorder(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk))
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then!***WARNING check compatibility [TODO]
    e_num_walker(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk)=e_num_walker(iorder(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk))
    e_den_walker(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk)=e_den_walker(iorder(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk))
  endif
  matrix_elements(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk)=matrix_elements(iorder(my_ndet_psi_t_connected+my_ndet_outside_ct+1:my_nwalk))
  my_nwalk = nwalk_eff + my_ndet_psi_t_connected+my_ndet_outside_ct
  return
  end subroutine sort_my_walkers3_up_dn
!-------------------------------------------------------------------------------

recursive subroutine merge_sort(key, iorder, nwalk, temp_i_1, temp_i_2)
! ==============================================================================
! Description   : Merge sort in ascending order.
!               : It sorts nwalk items in key.  Index is then used outside this routine to sort auxilliary items.
!               : temp_i_1 is an auxilliary array of length (nwalk+1)/2 of same type as key.
!               : temp_i_2 is an auxilliary array of length (nwalk+1)/2 of same type as iorder.
!               : In the present version key is type(ik_vec) and iorder is integer.
!               : temp_i_1 and temp_i_2 are passed in to avoid overhead of automatic arrays or allocations.
! Adapted by    : Cyrus Umrigar
! ------------------------------------------------------------------------------
  implicit none
  integer, intent(inout) :: key(nwalk), iorder(nwalk)
  integer, intent(in) :: nwalk
  integer,     dimension((nwalk+1)/2), intent (out) :: temp_i_1, temp_i_2 ! temporary arrays actually neither in nor out
! Local variables
  integer t0
  integer :: na,nb,i

  if (nwalk < 2) return
  if (nwalk == 2) then
    if (key(1) > key(2)) then
      t0 = key(1)
      key(1) = key(2)
      key(2) = t0
      i = iorder(1)
      iorder(1) = iorder(2)
      iorder(2) = i
    endif
    return
  endif

  na=(nwalk+1)/2
  nb=nwalk-na
  call merge_sort(key, iorder, na, temp_i_1, temp_i_2)
  call merge_sort(key(na+1), iorder(na+1), nb, temp_i_1, temp_i_2)

  if (key(na) > key(na+1)) then
    temp_i_1(1:na)=key(1:na)
    temp_i_2(1:na)=iorder(1:na)
    call merge(temp_i_1, temp_i_2, na, key(na+1), iorder(na+1), nb, key, iorder, nwalk)
  endif

  return
end subroutine merge_sort
!-------------------------------------------------------------------------------

recursive subroutine merge_sort2(key, iorder, nwalk, temp_i16, temp_i_2)
! ==============================================================================
! Description   : Merge sort in ascending order.
!               : It sorts nwalk items in key.  Index is then used outside this routine to sort auxilliary items.
!               : temp_i16 is an auxilliary array of length (nwalk+1)/2 of same type as key.
!               : temp_i_2 is an auxilliary array of length (nwalk+1)/2 of same type as iorder.
!               : In the present version key is type(ik_vec) and iorder is integer.
!               : temp_i16 and temp_i_2 are passed in to avoid overhead of automatic arrays or allocations.
! Adapted by    : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(inout) :: iorder(nwalk)
  integer, intent(in) :: nwalk
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: key(nwalk)
  type(ik_vec), dimension((nwalk+1)/2), intent (out) :: temp_i16 ! temporary array actually neither in nor out
  type(ik_vec) t0
#else
  integer(ik), intent(inout) :: key(nwalk)
  integer(ik), dimension((nwalk+1)/2), intent (out) :: temp_i16 ! temporary array actually neither in nor out
  integer(ik) t0
#endif
  integer,     dimension((nwalk+1)/2), intent (out) :: temp_i_2 ! temporary array actually neither in nor out
! Local variables
  integer :: na,nb,i

  if (nwalk < 2) return
  if (nwalk == 2) then
    if (key(1) > key(2)) then
      t0 = key(1)
      key(1) = key(2)
      key(2) = t0
      i = iorder(1)
      iorder(1) = iorder(2)
      iorder(2) = i
    endif
    return
  endif

  na=(nwalk+1)/2
  nb=nwalk-na
  call merge_sort2(key, iorder, na, temp_i16, temp_i_2)
  call merge_sort2(key(na+1), iorder(na+1), nb, temp_i16, temp_i_2)

  if (key(na) > key(na+1)) then
    temp_i16(1:na)=key(1:na)
    temp_i_2(1:na)=iorder(1:na)
    call merge2(temp_i16, temp_i_2, na, key(na+1), iorder(na+1), nb, key, iorder, nwalk)
  endif

  return
end subroutine merge_sort2
!-------------------------------------------------------------------------------

recursive subroutine merge_sort2_up_dn(key_up, key_dn, iorder, nwalk, temp_i16_up, temp_i16_dn, temp_i_2)
! ==============================================================================
! Description   : Merge sort in ascending order.
!               : It sorts nwalk items in key.  Index is then used outside this routine to sort auxilliary items.
!               : temp_i16 is an auxilliary array of length (nwalk+1)/2 of same type as key.
!               : temp_i_2 is an auxilliary array of length (nwalk+1)/2 of same type as iorder.
!               : In the present version key is type(ik_vec) and iorder is integer.
!               : temp_i16 and temp_i_2 are passed in to avoid overhead of automatic arrays or allocations.
! Adapted by    : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(inout) :: iorder(nwalk)
  integer, intent(in) :: nwalk
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: key_up(nwalk), key_dn(nwalk)
  type(ik_vec), dimension((nwalk+1)/2), intent (out) :: temp_i16_up, temp_i16_dn ! temporary array actually neither in nor out
  type(ik_vec) t0
#else
  integer(ik), intent(inout) :: key_up(nwalk), key_dn(nwalk)
  integer(ik), dimension((nwalk+1)/2), intent (out) :: temp_i16_up, temp_i16_dn ! temporary array actually neither in nor out
  integer(ik) t0
#endif
  integer,     dimension((nwalk+1)/2), intent (out) :: temp_i_2                 ! temporary array actually neither in nor out
! Local variables
  integer :: na,nb,i

  if (nwalk < 2) return
  if (nwalk == 2) then
    if (abs(key_up(1)) > abs(key_up(2))) then ! If key_up's in wrong order, sort on key_up
      t0 = key_up(1)
      key_up(1) = key_up(2)
      key_up(2) = t0
      t0 = key_dn(1)
      key_dn(1) = key_dn(2)
      key_dn(2) = t0
      i = iorder(1)
      iorder(1) = iorder(2)
      iorder(2) = i
    endif
    if (key_up(1) == key_up(2) .and. abs(key_dn(1)) > abs(key_dn(2))) then ! If key_up's are equal, sort on key_dn
      t0 = key_dn(1)
      key_dn(1) = key_dn(2)
      key_dn(2) = t0
      i = iorder(1)
      iorder(1) = iorder(2)
      iorder(2) = i
    endif
    return
  endif

  na=(nwalk+1)/2
  nb=nwalk-na
  call merge_sort2_up_dn(key_up, key_dn, iorder, na, temp_i16_up, temp_i16_dn, temp_i_2)
  call merge_sort2_up_dn(key_up(na+1), key_dn(na+1), iorder(na+1), nb, temp_i16_up, temp_i16_dn, temp_i_2)

  if (abs(key_up(na)) > abs(key_up(na+1)) .or. (key_up(na) == key_up(na+1) .and. abs(key_dn(na)) > abs(key_dn(na+1)))) then
    temp_i16_up(1:na)=key_up(1:na)
    temp_i16_dn(1:na)=key_dn(1:na)
    temp_i_2(1:na)=iorder(1:na)
    call merge2_up_dn(temp_i16_up, temp_i16_dn, temp_i_2, na, key_up(na+1), key_dn(na+1), iorder(na+1), nb, key_up, key_dn, iorder, nwalk)
  endif

  return
end subroutine merge_sort2_up_dn
!-------------------------------------------------------------------------------

subroutine merge(a,a2,na, b,b2,nb, c,c2,nc)
! ==============================================================================
! Description   : Called by merge_sort
! Adapted by    : Cyrus Umrigar
! ------------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: na,nb,nc               ! na+nb = nc
  integer    , intent(inout) :: a(na), c(nc)    ! b  overlays c(na+1:nc)
  integer    , intent(inout) :: a2(na), c2(nc)  ! b2 overlays c2(na+1:nc)
  integer    , intent(in)    :: b(nb)
  integer    , intent(in)    :: b2(nb)
  integer :: i,j,k

  i = 1; j = 1; k = 1;
  do while(i <= na .and. j <= nb)
    if (a(i) <= b(j)) then
      c(k) = a(i)
      c2(k) = a2(i)
      i = i+1
    else
      c(k) = b(j)
      c2(k) = b2(j)
      j = j+1
    endif
    k = k + 1
  enddo

  do while (i <= na)
    c(k) = a(i)
    c2(k) = a2(i)
    i = i + 1
    k = k + 1
  enddo

  return
end subroutine merge
!-------------------------------------------------------------------------------

subroutine merge2(a,a2,na, b,b2,nb, c,c2,nc)
! ==============================================================================
! Description   : Called by merge_sort2
! Adapted by    : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(in) :: na,nb,nc               ! na+nb = nc
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: a(na), c(nc)    ! b  overlays c(na+1:nc)
  type(ik_vec), intent(in)    :: b(nb)
#else
  integer(ik), intent(inout) :: a(na), c(nc)    ! b  overlays c(na+1:nc)
  integer(ik), intent(in)    :: b(nb)
#endif
  integer    , intent(inout) :: a2(na), c2(nc)  ! b2 overlays c2(na+1:nc)
  integer    , intent(in)    :: b2(nb)
  integer :: i,j,k

  i = 1; j = 1; k = 1;
  do while(i <= abs(na) .and. j <= abs(nb))
    if (abs(a(i)) <= abs(b(j))) then
      c(k) = a(i)
      c2(k) = a2(i)
      i = i+1
    else
      c(k) = b(j)
      c2(k) = b2(j)
      j = j+1
    endif
    k = k + 1
  enddo

  do while (i <= abs(na))
    c(k) = a(i)
    c2(k) = a2(i)
    i = i + 1
    k = k + 1
  enddo

  return
end subroutine merge2
!-------------------------------------------------------------------------------

subroutine merge2_up_dn(a_up,a_dn,a2,na, b_up,b_dn,b2,nb, c_up,c_dn,c2,nc)
! ==============================================================================
! Description   : Called by merge_sort2_up_dn
! Adapted by    : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(in) :: na,nb,nc                                      ! na+nb = nc
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: a_up(na), a_dn(na), c_up(nc), c_dn(nc) ! b_up overlays c_up(na+1:nc), and similarly for b_dn, c_dn
  type(ik_vec), intent(in)    :: b_up(nb), b_dn(nb)
  type(ik_vec) :: a_up_i, a_dn_i, b_up_j, b_dn_j
#else
  integer(ik), intent(inout) :: a_up(na), a_dn(na), c_up(nc), c_dn(nc) ! b_up overlays c_up(na+1:nc), and similarly for b_dn, c_dn
  integer(ik), intent(in)    :: b_up(nb), b_dn(nb)
  integer(ik) :: a_up_i, a_dn_i, b_up_j, b_dn_j
#endif
  integer    , intent(inout) :: a2(na), c2(nc)                         ! b2   overlays c2(na+1:nc)
  integer    , intent(in)    :: b2(nb)
  integer :: i,j,k
  integer :: na_abs, nb_abs

  na_abs = abs(na)
  nb_abs = abs(nb)

  i = 1; j = 1; k = 1;
  do while(i <= na_abs .and. j <= nb_abs)
    a_up_i = a_up(i)
    a_dn_i = a_dn(i)
    b_up_j = b_up(j)
    b_dn_j = b_dn(j)
    if (abs(a_up_i) < abs(b_up_j) .or. (a_up_i == b_up_j .and. abs(a_dn_i) <= abs(b_dn_j))) then
      c_up(k) = a_up_i
      c_dn(k) = a_dn_i
      c2(k) = a2(i)
      i = i+1
    else
      c_up(k) = b_up_j
      c_dn(k) = b_dn_j
      c2(k) = b2(j)
      j = j+1
    endif
    k = k + 1
  enddo

  do while (i <= na_abs)
    c_up(k) = a_up(i)
    c_dn(k) = a_dn(i)
    c2(k) = a2(i)
    i = i + 1
    k = k + 1
  enddo

  return
end subroutine merge2_up_dn
!-------------------------------------------------------------------------------

  subroutine merge_original_with_spawned(iblk,istep,walk_det,walk_wt,nwalk)
! ==============================================================================
! Description   : Merge walkers (hamiltonian_type fictitious, read, hubbard)
! Author        : Cyrus Umrigar
! Modified by Adam Holmes to put in semistochastic projection
! ------------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: iblk,istep
  integer, intent(inout) :: walk_det(nwalk), nwalk
  real(rk), intent(inout) :: walk_wt(nwalk)
  integer iwalk, nshift, i_permanent_initiator!, n_imp2

  nshift=0
  i_permanent_initiator=0
  do iwalk=2,nwalk
    if(walk_det(iwalk).eq.walk_det(iwalk-nshift-1)) then           ! This walker is on same det as previous walker
      nshift=nshift+1
!     The present scheme for combining walkers sequentially in pairs is a bit unsatisfactory because the same set of wts and initiators generated in a different order could result in a different final initiator.
!     An alternative, that cures this is to replace the next 13 lines by max(initiator(iwalk-nshift),initiator(iwalk)), but this could result in a high initiator value even though most of the wt comes from low initiator values.
      if(walk_wt(iwalk)*walk_wt(iwalk-nshift).gt.0) then           ! Same sign walkers being combined
        initiator(iwalk-nshift)=max(initiator(iwalk-nshift),initiator(iwalk))
        if(imp_distance(iwalk-nshift).ne.0) imp_distance(iwalk-nshift)=min(imp_distance(iwalk-nshift),abs(imp_distance(iwalk)))
      else                                                         ! Different sign
        if(abs(walk_wt(iwalk-nshift)).lt.abs(walk_wt(iwalk))) then ! Different sign and sum has same sign as the second of these walkers.
          if(initiator(iwalk-nshift).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to that of the second walker
            initiator(iwalk-nshift)=initiator(iwalk)
            if(imp_distance(iwalk-nshift).ne.0) imp_distance(iwalk-nshift)=abs(imp_distance(iwalk))
          endif
        elseif(abs(walk_wt(iwalk-nshift)).eq.abs(walk_wt(iwalk))) then  ! Different sign and they add to zero
          if(initiator(iwalk-nshift).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to the smaller of the two or equivalently to 0 since the wt is zero anyway
            initiator(iwalk-nshift)=0
          endif
        endif
      endif
      if (.not.((imp_distance(iwalk-nshift).eq.0.and.imp_distance(iwalk).eq.-1)))  walk_wt(iwalk-nshift)=walk_wt(iwalk-nshift)+walk_wt(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
      if(initiator(iwalk-nshift-1).eq.3 .and. r_initiator.ge.0) then    ! Make sure that permanent initiator has an abs wt of atleast 1 with the right sign
        i_permanent_initiator=i_permanent_initiator+1
        if(walk_wt(iwalk-nshift-1)*sign_permanent_initiator(i_permanent_initiator) < 1._rk) then
          if(n_warning.le.M_WARNING) write(6,'(''Warning: iblk, istep'',2i6,'' walk_wt of permanent initiator changed from'',f4.1,'' to'',i3)') iblk, istep, walk_wt(iwalk-nshift-1), sign_permanent_initiator(i_permanent_initiator)
          if(n_warning.eq.M_WARNING) write(6,'(''Warning: n_warning = M_WARNING, so no more warnings will be issued'')')
          n_warning=n_warning+1
          walk_wt(iwalk-nshift-1)=sign_permanent_initiator(i_permanent_initiator)
        endif
      elseif(initiator(iwalk-nshift-1).eq.2 .and. abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power)) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
     !elseif(initiator(iwalk-nshift-1).eq.2 .and. abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power)) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
          if(ipr.ge.1) write(6,'(''Remove from initiators1'',9i6)') walk_det(iwalk-nshift-1), initiator(iwalk-nshift-1)
!         initiator(iwalk-nshift-1)=min(initiator(iwalk-nshift-1),1)
          initiator(iwalk-nshift-1)=1
      elseif(initiator(iwalk-nshift-1).lt.2 .and. abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
     !elseif(initiator(iwalk-nshift-1).lt.2 .and. abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
          if(ipr.ge.1) write(6,'(''For det'',i8,'' incrementing initiator from'',i3,'' to'',i3)') walk_det(iwalk-nshift-1), initiator(iwalk-nshift-1), initiator(iwalk-nshift-1)+1
          initiator(iwalk-nshift-1)=initiator(iwalk-nshift-1)+1_i1b
      endif
      if(ipr.ge.1) write(6,'(''walk_wt(iwalk-nshift-1), initiator(iwalk-nshift-1), r_initiator, imp_distance(iwalk-nshift-1)'',f10.5,i5,f10.5,i5)') walk_wt(iwalk-nshift-1), initiator(iwalk-nshift-1), r_initiator, imp_distance(iwalk-nshift-1)

!     if(((walk_wt(iwalk-nshift-1).eq.0 .and. (initiator(iwalk-nshift-1).ne.3 .or. r_initiator.lt.0)) .or. initiator(iwalk-nshift-1).eq.0) .and. imp_distance(iwalk-nshift-1).ge.1) then
      if((walk_wt(iwalk-nshift-1).eq.0 .and. (initiator(iwalk-nshift-1).ne.3 .or. r_initiator.lt.0)) .or. (initiator(iwalk-nshift-1).eq.0 .and. imp_distance(iwalk-nshift-1).ge.1)) then
        nshift=nshift+1 ! If the sum of wts for the previous det is 0 and it is not a permanent initiator, or, if it was spawned from a noninitiator and the abs wt is below threshold then one can reuse this location
      endif
! Done working on previous det.  Now move position of current one.
      walk_det(iwalk-nshift)=walk_det(iwalk)
      walk_wt(iwalk-nshift)=walk_wt(iwalk)
      initiator(iwalk-nshift)=initiator(iwalk)
      imp_distance(iwalk-nshift)=imp_distance(iwalk)
      matrix_elements(iwalk-nshift)=matrix_elements(iwalk)
      if (imp_distance(iwalk).eq.-1) then
        imp_distance(iwalk-nshift)=1
!       write(6,'(''walk_dets_up(iwalk-nshift), walk_dets_dn(iwalk-nshift), walk_wt(iwalk-nshift), initiator(iwalk-nshift)'',2i8,f10.6,i3)') walk_dets_up(iwalk-nshift), walk_dets_dn(iwalk-nshift), walk_wt(iwalk-nshift), initiator(iwalk-nshift)
      endif
    endif
  enddo

! write(6,'(''nwalk, nshift='',2i6)') nwalk, nshift
! Fix up last det, just as we did above after "else".  This is needed when the last walker matches the previous walker so that the "else" above is not executed for this walker.
  if(initiator(nwalk-nshift).eq.3 .and. r_initiator.ge.0) then   ! Make sure that permanent initiator has an abs wt of atleast 1 with the right sign
    i_permanent_initiator=i_permanent_initiator+1
    if(walk_wt(nwalk-nshift)*sign_permanent_initiator(i_permanent_initiator) < 1._rk) then
      if(n_warning.le.M_WARNING) write(6,'(''Warning: iblk, istep'',2i6,'' walk_wt of permanent initiator changed from'',f4.1,'' to'',i3)') iblk, istep, walk_wt(nwalk-nshift), sign_permanent_initiator(i_permanent_initiator)
      if(n_warning.eq.M_WARNING) write(6,'(''Warning: n_warning = M_WARNING, so no more warnings will be issued'')')
      n_warning=n_warning+1
      walk_wt(nwalk-nshift)=sign_permanent_initiator(i_permanent_initiator)
    endif
  elseif(initiator(nwalk-nshift).eq.2 .and. abs(walk_wt(nwalk-nshift)).le.r_initiator*(max(0,imp_distance(nwalk-nshift)-initiator_min_distance)**initiator_power)) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
 !elseif(initiator(nwalk-nshift).eq.2 .and. abs(walk_wt(nwalk-nshift)).le.r_initiator*(imp_distance(nwalk-nshift)**initiator_power)) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
    if(ipr.ge.1) write(6,'(''Remove from initiators1'',9i6)') walk_det(nwalk-nshift), initiator(nwalk-nshift)
    initiator(nwalk-nshift)=min(initiator(nwalk-nshift),1_i1b)
  elseif(initiator(nwalk-nshift).lt.2 .and. abs(walk_wt(nwalk-nshift)).gt.r_initiator*(max(0,imp_distance(nwalk-nshift)-initiator_min_distance)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
 !elseif(initiator(nwalk-nshift).lt.2 .and. abs(walk_wt(nwalk-nshift)).gt.r_initiator*(imp_distance(nwalk-nshift)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
    if(ipr.ge.1) write(6,'(''For det'',i8,'' incrementing initiator from'',i3,'' to'',i3)') walk_det(nwalk-nshift), initiator(nwalk-nshift), initiator(nwalk-nshift)+1
    initiator(nwalk-nshift)=initiator(nwalk-nshift)+1_i1b
  endif
  if (imp_distance(iwalk-nshift).eq.-1) then ! Now set imp_distance=1 for states spawned from deterministic space to stochastic space
    imp_distance(iwalk-nshift)=1
  endif

! if(((walk_wt(nwalk-nshift).eq.0 .and. (initiator(nwalk-nshift).ne.3 .or. r_initiator.lt.0)) .or. initiator(nwalk-nshift).eq.0) .and. imp_distance(nwalk-nshift).ge.1) nshift=nshift+1 ! If the sum of wts for the previous det is 0 and it is not a permanent initiator, or, if it was spawned from a noninitiator and the abs wt is below threshold then one can reuse this location
  if((walk_wt(nwalk-nshift).eq.0 .and. (initiator(nwalk-nshift).ne.3 .or. r_initiator.lt.0)) .or. (initiator(nwalk-nshift).eq.0 .and. imp_distance(nwalk-nshift).ge.1)) nshift=nshift+1 ! If the sum of wts for the previous det is 0 and it is not a permanent initiator, or, if it was spawned from a noninitiator and the abs wt is below threshold then one can reuse this location

  nwalk=nwalk-nshift
! imp_distance(nwalk+1:MWALK)=1_i1b
  if(r_initiator.ge.0 .and. i_permanent_initiator/=n_permanent_initiator) then
    write(6,'(''i_permanent_initiator, n_permanent_initiator='',2i5)') i_permanent_initiator, n_permanent_initiator
    stop 'i_permanent_initiator /= n_permanent_initiator'
  endif

  !if(walk_wt(1).eq.0) then
  !  write(6,'(''nwalk='',i8)') nwalk
  !  write(6,'(''walk_wt='',1000f8.1)') (walk_wt(iwalk),iwalk=1,nwalk)
  !  write(6,'(''merge_original_with_spawned: no walkers left'')')
  !  stop 'merge_original_with_spawned: no walkers left'
  !endif

  if(ipr.ge.1) call write_walker_info('End of merge_original_with_spawned, ',nwalk)

  end subroutine merge_original_with_spawned
!-------------------------------------------------------------------------------------

  subroutine merge_original_with_spawned2_wo_initiator(iblk,istep,walk_dets_up,walk_dets_dn,walk_wt,nwalk)
! ==============================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified by Adam Holmes to put in semistochastic projection
! Not presently being used.
! ------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(in) :: iblk,istep
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
#else
  integer(ik), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
#endif
  integer, intent(inout) :: nwalk
  real(rk), intent(inout) :: walk_wt(nwalk)
  integer iwalk, nshift, i_permanent_initiator!, n_imp2

  nshift=0
  i_permanent_initiator=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
!     The present scheme for combining walkers sequentially in pairs is a bit unsatisfactory because the same set of wts and initiators generated in a different order could result in a different final initiator.
!     An alternative, that cures this is to replace the next 13 lines by max(initiator(iwalk-nshift),initiator(iwalk)), but this could result in a high initiator value even though most of the wt comes from low initiator values.
      if(walk_wt(iwalk)*walk_wt(iwalk-nshift).gt.0) then           ! Same sign walkers being combined
        initiator(iwalk-nshift)=max(initiator(iwalk-nshift),initiator(iwalk))
        if(imp_distance(iwalk-nshift).ne.0) imp_distance(iwalk-nshift)=min(imp_distance(iwalk-nshift),abs(imp_distance(iwalk)))
      else                                                         ! Different sign
        if(abs(walk_wt(iwalk-nshift)).lt.abs(walk_wt(iwalk))) then ! Different sign and sum has same sign as the second of these walkers.
          if(initiator(iwalk-nshift).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to that of the second walker
            initiator(iwalk-nshift)=initiator(iwalk)
            if(imp_distance(iwalk-nshift).ne.0) imp_distance(iwalk-nshift)=abs(imp_distance(iwalk))
          endif
        elseif(abs(walk_wt(iwalk-nshift)).eq.abs(walk_wt(iwalk))) then  ! Different sign and they add to zero
          if(initiator(iwalk-nshift).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to the smaller of the two or equivalently to 0 since the wt is zero anyway
            initiator(iwalk-nshift)=0
          endif
        endif
      endif
!     if(walk_dets_up(iwalk).eq.9007267974217728_ik) write(6,'(''1 iwalk, nshift, initiator(iwalk-nshift), imp_distance(iwalk-nshift)'',9i5)') iwalk, nshift, initiator(iwalk-nshift), imp_distance(iwalk-nshift)
      if (.not.((imp_distance(iwalk-nshift).eq.0.and.imp_distance(iwalk).eq.-1)))  walk_wt(iwalk-nshift)=walk_wt(iwalk-nshift)+walk_wt(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
      if(initiator(iwalk-nshift-1).eq.3 .and. r_initiator.ge.0) then    ! Make sure that permanent initiator has an abs wt of atleast 1 with the right sign
        i_permanent_initiator=i_permanent_initiator+1
        if(walk_wt(iwalk-nshift-1)*sign_permanent_initiator(i_permanent_initiator) < 1._rk) then
          if(n_warning.le.M_WARNING) write(6,'(''Warning: iblk, istep'',2i6,'' walk_wt of permanent initiator changed from'',f4.1,'' to'',i3)') iblk, istep, walk_wt(iwalk-nshift-1), sign_permanent_initiator(i_permanent_initiator)
          if(n_warning.eq.M_WARNING) write(6,'(''Warning: n_warning = M_WARNING, so no more warnings will be issued'')')
          n_warning=n_warning+1
          walk_wt(iwalk-nshift-1)=sign_permanent_initiator(i_permanent_initiator)
        endif
      elseif(initiator(iwalk-nshift-1).eq.2 .and. abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power)) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
     !elseif(initiator(iwalk-nshift-1).eq.2 .and. abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power)) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
          if(ipr.ge.1) write(6,'(''Remove from initiators1'',9i6)') walk_dets_up(iwalk-nshift-1), walk_dets_dn(iwalk-nshift-1), initiator(iwalk-nshift-1)
!         initiator(iwalk-nshift-1)=min(initiator(iwalk-nshift-1),1)
          initiator(iwalk-nshift-1)=1
      elseif(initiator(iwalk-nshift-1).lt.2 .and. abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
     !elseif(initiator(iwalk-nshift-1).lt.2 .and. abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
          if(ipr.ge.1) write(6,'(''For det'',2i8,'' incrementing initiator from'',i3,'' to'',i3)') walk_dets_up(iwalk-nshift-1), walk_dets_dn(iwalk-nshift-1), initiator(iwalk-nshift-1), initiator(iwalk-nshift-1)+1
          initiator(iwalk-nshift-1)=initiator(iwalk-nshift-1)+1_i1b
      endif
!     if(walk_dets_up(iwalk-nshift-1).eq.9007267974217728_ik) write(6,'(''iwalk, nshift, initiator(iwalk-nshift-1), imp_distance(iwalk-nshift-1)'',9i5)') iwalk, nshift, initiator(iwalk-nshift-1), imp_distance(iwalk-nshift-1)
      !if(((walk_wt(iwalk-nshift-1).eq.0 .and. (initiator(iwalk-nshift-1).ne.3 .or. r_initiator.lt.0)) .or. initiator(iwalk-nshift-1).eq.0) .and. imp_distance(iwalk-nshift-1).ge.1) then
      !  nshift=nshift+1 ! If the sum of wts for the previous det is 0 and it is not a permanent initiator, or, if it was spawned from a noninitiator and the abs wt is below threshold then one can reuse this location
      !endif
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      walk_wt(iwalk-nshift)=walk_wt(iwalk)
      initiator(iwalk-nshift)=initiator(iwalk)
      imp_distance(iwalk-nshift)=imp_distance(iwalk)
      matrix_elements(iwalk-nshift)=matrix_elements(iwalk)
!     if(walk_dets_up(iwalk-nshift).eq.9007267974217728_ik) write(6,'(''2 iwalk, nshift, initiator(iwalk-nshift), imp_distance(iwalk-nshift)'',9i5)') iwalk, nshift, initiator(iwalk-nshift), imp_distance(iwalk-nshift)
      if (imp_distance(iwalk).eq.-1) then
        imp_distance(iwalk-nshift)=1
!       write(6,'(''walk_dets_up(iwalk-nshift), walk_dets_dn(iwalk-nshift), walk_wt(iwalk-nshift), initiator(iwalk-nshift)'',2i8,f10.6,i3)') walk_dets_up(iwalk-nshift), walk_dets_dn(iwalk-nshift), walk_wt(iwalk-nshift), initiator(iwalk-nshift)
      endif
    endif
  enddo

! Fix up last det, just as we did above after "else".  This is needed when the last walker matches the previous walker so that the "else" above is not executed for this walker.
  if(initiator(nwalk-nshift).eq.3 .and. r_initiator.ge.0) then   ! Make sure that permanent initiator has an abs wt of atleast 1 with the right sign
    i_permanent_initiator=i_permanent_initiator+1
    if(walk_wt(nwalk-nshift)*sign_permanent_initiator(i_permanent_initiator) < 1._rk) then
      write(6,'(''Warning4:'',f24.17,es12.4)') walk_wt(iwalk-nshift-1),walk_wt(iwalk-nshift-1)-sign_permanent_initiator(i_permanent_initiator)
      if(n_warning.le.M_WARNING) write(6,'(''Warning: iblk, istep'',2i6,'' walk_wt of permanent initiator changed from'',f4.1,'' to'',i3)') iblk, istep, walk_wt(nwalk-nshift), sign_permanent_initiator(i_permanent_initiator)
      if(n_warning.eq.M_WARNING) write(6,'(''Warning: n_warning = M_WARNING, so no more warnings will be issued'')')
      n_warning=n_warning+1
      walk_wt(nwalk-nshift)=sign_permanent_initiator(i_permanent_initiator)
    endif
  elseif(initiator(nwalk-nshift).eq.2 .and. abs(walk_wt(nwalk-nshift)).le.r_initiator*(max(0,imp_distance(nwalk-nshift)-initiator_min_distance)**initiator_power)) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
 !elseif(initiator(nwalk-nshift).eq.2 .and. abs(walk_wt(nwalk-nshift)).le.r_initiator*(imp_distance(nwalk-nshift)**initiator_power)) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
    if(ipr.ge.1) write(6,'(''Remove from initiators1'',9i6)') walk_dets_up(nwalk-nshift), walk_dets_dn(nwalk-nshift), initiator(nwalk-nshift)
    initiator(nwalk-nshift)=min(initiator(nwalk-nshift),1_i1b)
  elseif(initiator(nwalk-nshift).lt.2 .and. abs(walk_wt(nwalk-nshift)).gt.r_initiator*(max(0,imp_distance(nwalk-nshift)-initiator_min_distance)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
 !elseif(initiator(nwalk-nshift).lt.2 .and. abs(walk_wt(nwalk-nshift)).gt.r_initiator*(imp_distance(nwalk-nshift)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
    if(ipr.ge.1) write(6,'(''For det'',2i8,'' incrementing initiator from'',i3,'' to'',i3)') walk_dets_up(nwalk-nshift), walk_dets_dn(nwalk-nshift), initiator(nwalk-nshift), initiator(nwalk-nshift)+1
    initiator(nwalk-nshift)=initiator(nwalk-nshift)+1_i1b
  endif
  if (imp_distance(iwalk-nshift).eq.-1) then ! Now set imp_distance=1 for states spawned from deterministic space to stochastic space
    imp_distance(iwalk-nshift)=1
  endif

  if(((walk_wt(nwalk-nshift).eq.0 .and. (initiator(nwalk-nshift).ne.3 .or. r_initiator.lt.0)) .or. initiator(nwalk-nshift).eq.0) .and. imp_distance(nwalk-nshift).ge.1) nshift=nshift+1 ! If the sum of wts for the previous det is 0 and it is not a permanent initiator, or, if it was spawned from a noninitiator and the abs wt is below threshold then one can reuse this location

  nwalk=nwalk-nshift
! imp_distance(nwalk+1:MWALK)=1_i1b
  if(r_initiator.ge.0 .and. i_permanent_initiator/=n_permanent_initiator) then
    write(6,'(''i_permanent_initiator, n_permanent_initiator='',2i5)') i_permanent_initiator, n_permanent_initiator
    write(6,'(''walk_dets_up='',1000i20)') walk_dets_up(1:nwalk)
    write(6,'(''walk_dets_dn='',1000i20)') walk_dets_dn(1:nwalk)
    write(6,'(''walk_wt     ='',1000f20.3)') walk_wt(1:nwalk)
    write(6,'(''initiator   ='',1000i20)') initiator(1:nwalk)
    stop 'i_permanent_initiator /= n_permanent_initiator'
  endif

  !Temporary comment - HJC Sept 25 2012
  !if(walk_wt(1).eq.0) then
  !  write(6,'(''nwalk='',i8)') nwalk
  !  write(6,'(''walk_wt='',1000f8.1)') (walk_wt(iwalk),iwalk=1,nwalk)
  !  write(6,'(''merge_original_with_spawned2: no walkers left'')')
  !  stop 'merge_original_with_spawned2: no walkers left'
  !endif
  !END Temporary comment - HJC Sept 25 2012

  if(ipr.ge.1) call write_walker_info('End of merge_original_with_spawned2_wo_initiator,',nwalk)

  end subroutine merge_original_with_spawned2_wo_initiator
!-------------------------------------------------------------------------------------

  subroutine merge_original_with_spawned2(iblk,istep,walk_dets_up,walk_dets_dn,walk_wt,nwalk)
!  subroutine merge_original_with_spawned2(iblk,istep,walk_dets_up,walk_dets_dn,walk_wt,nwalk,my_n_perm_init)
! ==============================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified by Adam Holmes to put in semistochastic projection
! Merges walk_dets_up/dn, walk_wt, initiator, imp_distance
! ------------------------------------------------------------------------------
  use types, only: ik
  use common_psi_t, only : hf_to_psit
  use common_ham, only : hamiltonian_type

  implicit none
  integer, intent(in) :: iblk,istep!,my_n_perm_init

#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
#else
  integer(ik), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
#endif
  integer, intent(inout) :: nwalk
  real(rk), intent(inout) :: walk_wt(nwalk)
  integer iwalk, nshift, i_permanent_initiator!, n_imp2


  nshift=0
  i_permanent_initiator=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
!     The present scheme for combining walkers sequentially in pairs is a bit unsatisfactory because the same set of wts and initiators generated in a different order could result in a different final initiator.
!     An alternative, that cures this is to replace the next 13 lines by max(initiator(iwalk-nshift),initiator(iwalk)), but this could result in a high initiator value even though most of the wt comes from low initiator values.
      if(walk_wt(iwalk)*walk_wt(iwalk-nshift).gt.0) then           ! Same sign walkers being combined
        initiator(iwalk-nshift)=max(initiator(iwalk-nshift),initiator(iwalk))
        if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
          e_num_walker(iwalk-nshift)=min(e_num_walker(iwalk-nshift),e_num_walker(iwalk))         ! HJC: If calculated energy, then use minimum value as default is positive and set to 1.e51_rk
          e_den_walker(iwalk-nshift)=min(e_den_walker(iwalk-nshift),e_den_walker(iwalk))
        endif
        matrix_elements(iwalk-nshift) = min(matrix_elements(iwalk-nshift),matrix_elements(iwalk))
        ! Edited by AAH on 8 Jan 2013
        if(imp_distance(iwalk-nshift).eq.-2) then
          if (imp_distance(iwalk).eq.0) then
            imp_distance(iwalk-nshift)=0
          endif
        elseif (imp_distance(iwalk).eq.-2) then
          if(imp_distance(iwalk-nshift).ne.0) then
            imp_distance(iwalk-nshift)=-2
          endif
        elseif(imp_distance(iwalk-nshift).ne.0 .and. imp_distance(iwalk-nshift).ne.-2) then
          imp_distance(iwalk-nshift)=min(imp_distance(iwalk-nshift),abs(imp_distance(iwalk)))
        endif
       !if(imp_distance(iwalk-nshift).ne.0) imp_distance(iwalk-nshift)=min(imp_distance(iwalk-nshift),abs(imp_distance(iwalk)))
        ! End of edit
      else                                                         ! Different sign
!Warning:  I don't think the next 4 lines are necessary since the 2 walkers are on the same state
!       if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
!         e_num_walker(iwalk-nshift)=min(e_num_walker(iwalk-nshift),e_num_walker(iwalk))         ! HJC : If calculated energy, then use minimum value as default is positive and set to 1.e51_rk
!         e_den_walker(iwalk-nshift)=min(e_den_walker(iwalk-nshift),e_den_walker(iwalk))
!       endif
        matrix_elements(iwalk-nshift) = min(matrix_elements(iwalk-nshift),matrix_elements(iwalk))
        ! Edited by AAH on 8 Jan 2013
        if(imp_distance(iwalk-nshift).eq.-2) then
          if (imp_distance(iwalk).eq.0) then
            imp_distance(iwalk-nshift)=0
          endif
        elseif (imp_distance(iwalk).eq.-2) then
          if(imp_distance(iwalk-nshift).ne.0) then
            imp_distance(iwalk-nshift)=-2
          endif
        elseif(imp_distance(iwalk-nshift).ne.0 .and. imp_distance(iwalk-nshift).ne.-2) then
          imp_distance(iwalk-nshift)=min(imp_distance(iwalk-nshift),abs(imp_distance(iwalk)))
        endif
        ! End of edit
        if(abs(walk_wt(iwalk-nshift)).lt.abs(walk_wt(iwalk))) then ! Different sign and sum has same sign as the second of these walkers.
          if(initiator(iwalk-nshift).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to that of the second walker
            initiator(iwalk-nshift)=initiator(iwalk)
          endif
        elseif(abs(walk_wt(iwalk-nshift)).eq.abs(walk_wt(iwalk))) then  ! Different sign and they add to zero
          if(initiator(iwalk-nshift).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to the smaller of the two or equivalently to 0 since the wt is zero anyway
            initiator(iwalk-nshift)=0
          endif
        endif
      endif
!     if(walk_dets_up(iwalk).eq.9007267974217728_ik) write(6,'(''1 iwalk, nshift, initiator(iwalk-nshift), imp_distance(iwalk-nshift)'',9i5)') iwalk, nshift, initiator(iwalk-nshift), imp_distance(iwalk-nshift)
      if (.not.((imp_distance(iwalk-nshift).eq.0.and.imp_distance(iwalk).eq.-1)))  walk_wt(iwalk-nshift)=walk_wt(iwalk-nshift)+walk_wt(iwalk) ! Combine wt only if parent and child det are not both in deterministic space
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
      if(initiator(iwalk-nshift-1).eq.3 .and. r_initiator.ge.0) then    ! Make sure that permanent initiator has an abs wt of atleast 1 with the right sign
        i_permanent_initiator=i_permanent_initiator+1
        if(walk_wt(iwalk-nshift-1)*sign_permanent_initiator(i_permanent_initiator) < 1._rk) then
          if(n_warning.le.M_WARNING) write(6,'(''Warning: iblk, istep'',2i6,'' walk_wt of permanent initiator changed from'',f4.1,'' to'',i3)') iblk, istep, walk_wt(iwalk-nshift-1), sign_permanent_initiator(i_permanent_initiator)
          if(n_warning.eq.M_WARNING) write(6,'(''Warning: n_warning = M_WARNING, so no more warnings will be issued'')')
          n_warning=n_warning+1
          walk_wt(iwalk-nshift-1)=sign_permanent_initiator(i_permanent_initiator)
        endif
      elseif(initiator(iwalk-nshift-1).eq.2 .and. ((abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power).and.imp_distance(iwalk-nshift-1)>0) .or. ((abs(walk_wt(iwalk-nshift-1)).le.r_initiator.and..not.c_t_initiator).and.imp_distance(iwalk-nshift-1)==-2))) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
     !elseif(initiator(iwalk-nshift-1).eq.2 .and. abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power) .and. imp_distance(iwalk-nshift-1)>=0) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
          if(ipr.ge.1) write(6,'(''Remove from initiators1'',9i6)') walk_dets_up(iwalk-nshift-1), walk_dets_dn(iwalk-nshift-1), initiator(iwalk-nshift-1)
          initiator(iwalk-nshift-1)=1
      elseif(initiator(iwalk-nshift-1).lt.2 .and. ((abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power).and.imp_distance(iwalk-nshift-1)>=0) .or. ((abs(walk_wt(iwalk-nshift-1)).gt.r_initiator.or.c_t_initiator).and.imp_distance(iwalk-nshift-1)==-2))) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
     !elseif(initiator(iwalk-nshift-1).lt.2 .and. abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
          if(ipr.ge.1) write(6,'(''For det'',2i8,'' incrementing initiator from'',i3,'' to'',i3)') walk_dets_up(iwalk-nshift-1), walk_dets_dn(iwalk-nshift-1), initiator(iwalk-nshift-1), initiator(iwalk-nshift-1)+1
          initiator(iwalk-nshift-1)=initiator(iwalk-nshift-1)+1_i1b
      endif
!     if(walk_dets_up(iwalk-nshift-1).eq.9007267974217728_ik) write(6,'(''iwalk, nshift, initiator(iwalk-nshift-1), imp_distance(iwalk-nshift-1)'',9i5)') iwalk, nshift, initiator(iwalk-nshift-1), imp_distance(iwalk-nshift-1)
      if(((walk_wt(iwalk-nshift-1).eq.0 .and. (initiator(iwalk-nshift-1).ne.3 .or. r_initiator.lt.0)) .or. initiator(iwalk-nshift-1).eq.0) .and. imp_distance(iwalk-nshift-1).ge.1) then
        nshift=nshift+1 ! If the sum of wts for the previous det is 0 and it is not a permanent initiator, or, if it was spawned from a noninitiator and the abs wt is below threshold then one can reuse this location
      endif
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      walk_wt(iwalk-nshift)=walk_wt(iwalk)
      initiator(iwalk-nshift)=initiator(iwalk)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(iwalk-nshift)=e_num_walker(iwalk)
        e_den_walker(iwalk-nshift)=e_den_walker(iwalk)
      endif
      imp_distance(iwalk-nshift)=imp_distance(iwalk)
      matrix_elements(iwalk-nshift)=matrix_elements(iwalk)
!     if(walk_dets_up(iwalk-nshift).eq.9007267974217728_ik) write(6,'(''2 iwalk, nshift, initiator(iwalk-nshift), imp_distance(iwalk-nshift)'',9i5)') iwalk, nshift, initiator(iwalk-nshift), imp_distance(iwalk-nshift)
      if (imp_distance(iwalk).eq.-1) then
        imp_distance(iwalk-nshift)=1
!       write(6,'(''walk_dets_up(iwalk-nshift), walk_dets_dn(iwalk-nshift), walk_wt(iwalk-nshift), initiator(iwalk-nshift)'',2i8,f10.6,i3)') walk_dets_up(iwalk-nshift), walk_dets_dn(iwalk-nshift), walk_wt(iwalk-nshift), initiator(iwalk-nshift)
      endif
    endif
  enddo

! Fix up last det, just as we did above after "else".  This is needed when the last walker matches the previous walker so that the "else" above is not executed for this walker.
  if(initiator(nwalk-nshift).eq.3 .and. r_initiator.ge.0) then   ! Make sure that permanent initiator has an abs wt of atleast 1 with the right sign
    i_permanent_initiator=i_permanent_initiator+1
    if(walk_wt(nwalk-nshift)*sign_permanent_initiator(i_permanent_initiator) < 1._rk) then
      write(6,'(''Warning4:'',f24.17,es12.4)') walk_wt(iwalk-nshift-1),walk_wt(iwalk-nshift-1)-sign_permanent_initiator(i_permanent_initiator)
      if(n_warning.le.M_WARNING) write(6,'(''Warning: iblk, istep'',2i6,'' walk_wt of permanent initiator changed from'',f4.1,'' to'',i3)') iblk, istep, walk_wt(nwalk-nshift), sign_permanent_initiator(i_permanent_initiator)
      if(n_warning.eq.M_WARNING) write(6,'(''Warning: n_warning = M_WARNING, so no more warnings will be issued'')')
      n_warning=n_warning+1
      walk_wt(nwalk-nshift)=sign_permanent_initiator(i_permanent_initiator)
    endif

  elseif(initiator(iwalk-nshift-1).eq.2 .and. ((abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power).and.imp_distance(iwalk-nshift-1)>0) .or. ((abs(walk_wt(iwalk-nshift-1)).le.r_initiator.and..not.c_t_initiator).and.imp_distance(iwalk-nshift-1)==-2))) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
 !elseif(initiator(iwalk-nshift-1).eq.2 .and. ((abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power).and.imp_distance(iwalk-nshift-1)>0) .or. (abs(walk_wt(iwalk-nshift-1)).le.r_initiator.and.imp_distance(iwalk-nshift-1)==-2))) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
 !elseif(initiator(iwalk-nshift-1).eq.2 .and. ((abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power).and.imp_distance(iwalk-nshift-1)>0) .or. (abs(walk_wt(iwalk-nshift-1)).le.r_initiator.and.imp_distance(iwalk-nshift-1)==-2))) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
  ! Edited by AAH on 18 Jan 2013. Makes psi_t connections permanent initiators if hf_to_psit is true only if c_t_initiator is true.
 !elseif(initiator(iwalk-nshift-1).eq.2 .and. abs(walk_wt(iwalk-nshift-1)).le.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power) .and. imp_distance(iwalk-nshift-1)>=0) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
  ! End of edit
      if(ipr.ge.1) write(6,'(''Remove from initiators1'',9i6)') walk_dets_up(iwalk-nshift-1), walk_dets_dn(iwalk-nshift-1), initiator(iwalk-nshift-1)
!     initiator(iwalk-nshift-1)=min(initiator(iwalk-nshift-1),1)
      initiator(iwalk-nshift-1)=1
  elseif(initiator(iwalk-nshift-1).lt.2 .and. ((abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power).and.imp_distance(iwalk-nshift-1)>=0) .or. ((c_t_initiator.or.abs(walk_wt(iwalk-nshift-1)).gt.r_initiator).and.imp_distance(iwalk-nshift-1)==-2))) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
 !elseif(initiator(iwalk-nshift-1).lt.2 .and. ((abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(max(0,imp_distance(iwalk-nshift-1)-initiator_min_distance)**initiator_power).and.imp_distance(iwalk-nshift-1)>=0) .or. (abs(walk_wt(iwalk-nshift-1)).gt.r_initiator.and.imp_distance(iwalk-nshift-1)==-2))) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
 !elseif(initiator(iwalk-nshift-1).lt.2 .and. ((abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power).and.imp_distance(iwalk-nshift-1)>=0) .or. (abs(walk_wt(iwalk-nshift-1)).gt.r_initiator.and.imp_distance(iwalk-nshift-1)==-2))) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
  ! Edited by AAH on 18 Jan 2013. Makes psi_t connections permanent initiators if hf_to_psit is true only if c_t_initiator is true.
 !elseif(initiator(iwalk-nshift-1).lt.2 .and. abs(walk_wt(iwalk-nshift-1)).gt.r_initiator*(imp_distance(iwalk-nshift-1)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1

! ! Edited by AAH on 8 Jan 2013
! elseif(initiator(nwalk-nshift).eq.2 .and. abs(walk_wt(nwalk-nshift)).le.r_initiator*(imp_distance(nwalk-nshift)**initiator_power) .and. imp_distance(nwalk-nshift)>=0) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
!!elseif(initiator(nwalk-nshift).eq.2 .and. abs(walk_wt(nwalk-nshift)).le.r_initiator*(imp_distance(nwalk-nshift)**initiator_power)) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
! ! End of edit
!   if(ipr.ge.1) write(6,'(''Remove from initiators1'',9i6)') walk_dets_up(nwalk-nshift), walk_dets_dn(nwalk-nshift), initiator(nwalk-nshift)
!   initiator(nwalk-nshift)=min(initiator(nwalk-nshift),1_i1b)
! elseif(initiator(nwalk-nshift).lt.2 .and. abs(walk_wt(nwalk-nshift)).gt.r_initiator*(imp_distance(nwalk-nshift)**initiator_power)) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1

  ! End of edit
    if(ipr.ge.1) write(6,'(''For det'',2i8,'' incrementing initiator from'',i3,'' to'',i3)') walk_dets_up(nwalk-nshift), walk_dets_dn(nwalk-nshift), initiator(nwalk-nshift), initiator(nwalk-nshift)+1
    initiator(nwalk-nshift)=initiator(nwalk-nshift)+1_i1b
  endif

!***b
  do iwalk=1,nwalk-nshift
    if (imp_distance(iwalk).eq.-1) then ! Now set imp_distance=1 for states spawned from deterministic space to stochastic space
      imp_distance(iwalk)=1
    endif
  enddo
!***e
  if(((walk_wt(nwalk-nshift).eq.0 .and. (initiator(nwalk-nshift).ne.3 .or. r_initiator.lt.0)) .or. initiator(nwalk-nshift).eq.0) .and. imp_distance(nwalk-nshift).ge.1) nshift=nshift+1 ! If the sum of wts for the previous det is 0 and it is not a permanent initiator, or, if it was spawned from a noninitiator and the abs wt is below threshold then one can reuse this location

  nwalk=nwalk-nshift
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(nwalk+1:MWALK)=1.e51_rk
    e_den_walker(nwalk+1:MWALK)=1.e51_rk
  endif
  matrix_elements(nwalk+1:MWALK)=1.e51_rk

!Temporary comment - AAH 12 Dec 2013
! if(r_initiator.ge.0 .and. i_permanent_initiator/=my_n_perm_init) then
!   write(6,*) "coreid =",whoami
!   write(6,'(''i_permanent_initiator, n_permanent_initiator='',2i5)') i_permanent_initiator, my_n_perm_init,n_permanent_initiator
!   write(6,'(''walk_dets_up='',1000i20)') walk_dets_up(1:nwalk)
!   write(6,'(''walk_dets_dn='',1000i20)') walk_dets_dn(1:nwalk)
!   write(6,'(''walk_wt     ='',1000f20.3)') walk_wt(1:nwalk)
!   write(6,'(''initiator   ='',1000i20)') initiator(1:nwalk)
!   call flush(6)
!   stop 'i_permanent_initiator /= my_n_perm_init'
! endif
! End temp comment

  !Temporary comment - HJC Sept 25 2012
  !if(walk_wt(1).eq.0) then
  !  write(6,'(''nwalk='',i8)') nwalk
  !  write(6,'(''walk_wt='',1000f8.1)') (walk_wt(iwalk),iwalk=1,nwalk)
  !  write(6,'(''merge_original_with_spawned2: no walkers left'')')
  !  stop 'merge_original_with_spawned2: no walkers left'
  !endif
  !END Temporary comment - HJC Sept 25 2012

  if(ipr.ge.1) call write_walker_info('End of merge_original_with_spawned2,',nwalk)

! n_imp2=0
! do iwalk=1,nwalk
!   if(imp_distance(iwalk).eq.0) n_imp2=n_imp2+1
!   if(imp_distance(iwalk).eq.-1) write(6,'(''nwalk, iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)='',2i4,2i20,f10.6)') nwalk, iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)
!   if(imp_distance(iwalk).eq.-1) write(6,'(''nwalk, iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)='',2i4,2b60,f10.6)') nwalk, iwalk, walk_dets_up(iwalk), walk_dets_dn(iwalk), walk_wt(iwalk)
! enddo
! if(n_imp2.ne.n_imp) then
!   write(6,'(''after merge: n_imp, n_imp2'',9i8)') n_imp, n_imp2
!   stop 'after merge: n_imp2 \= n_imp'
! endif
! call flush(6)

  end subroutine merge_original_with_spawned2
!-------------------------------------------------------------------------------------

! ==============================================================================
  subroutine merge_original_with_spawned3 !(walk_dets_up,walk_dets_dn,walk_wt,nwalk)
! ==============================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified by Adam Holmes to put in semistochastic projection
! Modified : A Holmes, 28 Jan 2013. First add the newly spawned dets into the connections to psi trial, then merges them with the list outside the connections from the previous iterations.
! ------------------------------------------------------------------------------
  use types, only: ik,i8b
  use common_psi_t, only : hf_to_psit
  use common_ham, only : hamiltonian_type
  use common_run, only : ndet_outside_ct

  implicit none
 !integer, intent(in) :: iblk,istep
 !type(ik_vec), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
 !integer, intent(inout) :: nwalk
 !real(rk), intent(inout) :: walk_wt(nwalk)
  integer iwalk, i_permanent_initiator
  integer :: i_ct,start_ct,i_outside_ct,start_outside_ct,start_empty_outside_ct
  logical :: found,discard_walker,merge_with_prev,reorder
  integer :: i,j
  integer(i8b) :: from,to,first,last,shift

 !write (6,*) "merge_original_with_spawned3 called; nwalk=",nwalk,", ndet_psi_t_connected=",ndet_psi_t_connected,", ndet_outside_ct=",ndet_outside_ct
 !call flush(6)
!***b
!  if (my_nwalk<my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1))  stop 'nwalk<ndet_psi_t_connected+ndet_outside_ct'
  if (nwalk<ndet_psi_t_connected+ndet_outside_ct) then
    stop 'nwalk<ndet_psi_t_connected+ndet_outside_ct'
  endif
!***e
  i_permanent_initiator=0
  leave_out=0 ! Number of dets in outside_ct that will be discarded
  num_to_move=0 ! Number of newly spawned dets to insert into the list outside of ct
  merge_with_prev=.false. ! Whether a walker that is not in ct or outside_ct needs to be merged with the one to its left

!***b
!  if (my_nwalk==my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)) then ! Already sorted, but still need to apply initiator criterion to all dets
  if (nwalk==ndet_psi_t_connected+ndet_outside_ct) then ! Already sorted, but still need to apply initiator criterion to all dets
!***e
   !if (.not.c_t_initiator) then
   !  do i_ct=1,ndet_psi_t_connected
   !    call check_initiator(i_ct,discard_walker,i_permanent_initiator) ! Don't discard walkers here because we keep all dets in ct.
   !  enddo
   !else ! This is necessary only to allow us to set w_perm_initiator to 1 if it is too small.
   !  do i_ct=1,n_permanent_initiator
   !    call check_initiator(i_ct,discard_walker,i_permanent_initiator)
   !  enddo
   !endif
!***b
!    do i_outside_ct=my_ndet_psi_t_connected(whoami+1)+1,my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)
    do i_outside_ct=ndet_psi_t_connected+1,ndet_psi_t_connected+ndet_outside_ct
!***e
      call check_initiator(i_outside_ct,discard_walker,i_permanent_initiator)
      if (discard_walker) then
        leave_out=leave_out+1
        leave_out_list(leave_out)=i_outside_ct
      endif
    enddo

  else

    start_ct=1
!***b
!    start_outside_ct=my_ndet_psi_t_connected(whoami+1)+1
!    do iwalk=my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1,my_nwalk ! Scan over list of (sorted) newly spawned walkers
    start_outside_ct=ndet_psi_t_connected+1
    do iwalk=ndet_psi_t_connected+ndet_outside_ct+1,nwalk ! Scan over list of (sorted) newly spawned walkers
!***e
      found=.false.
      reorder=.false.
      ! First check to see if it's in list of ct
!***b
!      do i_ct=start_ct,my_ndet_psi_t_connected(whoami+1)
      do i_ct=start_ct,ndet_psi_t_connected
!***e
        if(walk_dets_up(iwalk).eq.walk_dets_up(i_ct) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(i_ct)) then ! This walker is on same det as previous walker

          if(walk_wt(iwalk)*walk_wt(i_ct).gt.0) then           ! Same sign walkers being combined
            initiator(i_ct)=max(initiator(i_ct),initiator(iwalk))
          else                                                         ! Different sign
            if(abs(walk_wt(i_ct)).lt.abs(walk_wt(iwalk))) then ! Different sign and sum has same sign as the second of these walkers.
              if(initiator(i_ct).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to that of the second walker
                initiator(i_ct)=initiator(iwalk)
              endif
            elseif(abs(walk_wt(i_ct)).eq.abs(walk_wt(iwalk))) then  ! Different sign and they add to zero
              if(initiator(i_ct).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to the smaller of the two or equivalently to 0 since the wt is zero anyway
                initiator(i_ct)=0
              endif
            endif
          endif
          if (.not.(imp_distance(i_ct).eq.0.and.imp_distance(iwalk).eq.-1))  walk_wt(i_ct)=walk_wt(i_ct)+walk_wt(iwalk) ! Combine wt only if parent and child det are not both in deterministic space

          found=.true.
          start_ct=i_ct ! not i_ct + 1 because there could be duplicates in the newly spawned dets
          exit

        elseif (walk_dets_up(iwalk)>walk_dets_up(i_ct) .or. (walk_dets_up(iwalk)==walk_dets_up(i_ct).and.walk_dets_dn(iwalk)>walk_dets_dn(i_ct))) then ! This walker (i_ct) in ct is not being stochastically spawned onto anymore. Just check initiator criterion.

         !if ((.not.c_t_initiator).or.i_ct<n_permanent_initiator) then
         !  call check_initiator(i_ct,discard_walker,i_permanent_initiator) ! Don't discard walkers here because we keep all dets in ct.
         !endif
          start_ct=i_ct+1
          ! Would like to be able to calculate energy for psi_t_connections here, but we have to multiply by S^-1 first.

        else ! No reason to keep scanning ct, since there won't be any more matches.

          start_ct = i_ct
          exit

        endif

      enddo ! i_ct

      if (found)  cycle

      ! Now loop over the outside_ct dets from the previous MC iteration
!***b
!      do i_outside_ct=start_outside_ct,my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)
      do i_outside_ct=start_outside_ct,ndet_psi_t_connected+ndet_outside_ct
!***e
        if(walk_dets_up(iwalk).eq.walk_dets_up(i_outside_ct) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(i_outside_ct)) then ! This walker is on same det as previous walker

          if(walk_wt(iwalk)*walk_wt(i_outside_ct).gt.0) then           ! Same sign walkers being combined
            initiator(i_outside_ct)=max(initiator(i_outside_ct),initiator(iwalk))
            imp_distance(i_outside_ct)=min(imp_distance(i_outside_ct),abs(imp_distance(iwalk)))
          else                                                         ! Different sign
            imp_distance(i_outside_ct)=min(imp_distance(i_outside_ct),abs(imp_distance(iwalk)))
            if(abs(walk_wt(i_outside_ct)).lt.abs(walk_wt(iwalk))) then ! Different sign and sum has same sign as the second of these walkers.
              if(initiator(i_outside_ct).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to that of the second walker
                initiator(i_outside_ct)=initiator(iwalk)
              endif
            elseif(abs(walk_wt(i_outside_ct)).eq.abs(walk_wt(iwalk))) then  ! Different sign and they add to zero
              if(initiator(i_outside_ct).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to the smaller of the two or equivalently to 0 since the wt is zero anyway
                initiator(i_outside_ct)=0
              endif
            endif
          endif
          walk_wt(i_outside_ct)=walk_wt(i_outside_ct)+walk_wt(iwalk)

          found=.true.
          start_outside_ct=i_outside_ct ! not i_outside_ct + 1 because there could be duplicates
          exit

        elseif (walk_dets_up(i_outside_ct)<walk_dets_up(iwalk).or.(walk_dets_up(i_outside_ct)==walk_dets_up(iwalk).and.walk_dets_dn(i_outside_ct)<walk_dets_dn(iwalk))) then ! This walker (i_outside_ct) in outside ct is not being stochastically spawned onto anymore. Just check initiator criterion.

          call check_initiator(i_outside_ct,discard_walker,i_permanent_initiator)
          if (discard_walker) then
            leave_out=leave_out+1
            leave_out_list(leave_out)=i_outside_ct
          endif
          start_outside_ct = i_outside_ct+1

        else ! No reason to keep scanning outside_ct, since there won't be any more matches.
        ! We passed the slot where we want to insert det iwalk in the list of saved outside_ct dets
          start_outside_ct = i_outside_ct
          reorder = .true.
          exit

        endif

      enddo ! i_outside_ct

      if (.not.found) then

        if (merge_with_prev) then
          ! combine iwalk with iwalk-1, store at iwalk

          if(walk_wt(iwalk-1)*walk_wt(iwalk).gt.0) then           ! Same sign walkers being combined
            initiator(iwalk)=max(initiator(iwalk),initiator(iwalk-1))
            if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
              e_num_walker(iwalk)=min(e_num_walker(iwalk),e_num_walker(iwalk-1))         ! HJC: If calculated energy, then use minimum value as default is positive and set to 1.e51_rk
              e_den_walker(iwalk)=min(e_den_walker(iwalk),e_den_walker(iwalk-1))
            endif
            imp_distance(iwalk)=min(imp_distance(iwalk),abs(imp_distance(iwalk-1)))
            matrix_elements(iwalk)=min(matrix_elements(iwalk),(matrix_elements(iwalk-1)))
          else                                                         ! Different sign
            if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
              e_num_walker(iwalk)=min(e_num_walker(iwalk),e_num_walker(iwalk-1))         ! HJC : If calculated energy, then use minimum value as default is positive and set to 1.e51_rk
              e_den_walker(iwalk)=min(e_den_walker(iwalk),e_den_walker(iwalk-1))
            endif
            imp_distance(iwalk)=min(imp_distance(iwalk),abs(imp_distance(iwalk-1)))
            matrix_elements(iwalk)=min(matrix_elements(iwalk),(matrix_elements(iwalk-1)))
            if(abs(walk_wt(iwalk)).lt.abs(walk_wt(iwalk-1))) then ! Different sign and sum has same sign as the second of these walkers.
              if(initiator(iwalk).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to that of the second walker
                initiator(iwalk)=initiator(iwalk-1)
              endif
            elseif(abs(walk_wt(iwalk)).eq.abs(walk_wt(iwalk-1))) then  ! Different sign and they add to zero
              if(initiator(iwalk).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to the smaller of the two or equivalently to 0 since the wt is zero anyway
                initiator(iwalk)=0
              endif
            endif
          endif
          walk_wt(iwalk)=walk_wt(iwalk)+walk_wt(iwalk-1)
          merge_with_prev=.false.

        endif

!***b
!        if (iwalk<my_nwalk.and.walk_dets_up(iwalk)==walk_dets_up(iwalk+1).and.walk_dets_dn(iwalk)==walk_dets_dn(iwalk+1)) then ! Try to combine this with other newly spawned walkers
        if (iwalk<nwalk.and.walk_dets_up(iwalk)==walk_dets_up(iwalk+1).and.walk_dets_dn(iwalk)==walk_dets_dn(iwalk+1)) then ! Try to combine this with other newly spawned walkers
!***e
          merge_with_prev=.true.
        else
          ! Check initiator criterion on det iwalk.
          call check_initiator(iwalk,discard_walker,i_permanent_initiator)
          if (.not.discard_walker) then
            ! Add to list of new walkers that we keep
            num_to_move=num_to_move+1
            if (reorder) then
              empty_outside_ct(num_to_move)=i_outside_ct
            else
              ! This det is not in ct or outside_ct and should be added to the end.
!***b
!              empty_outside_ct(num_to_move)=my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)+1
              empty_outside_ct(num_to_move)=ndet_psi_t_connected+ndet_outside_ct+1
!***e
            endif
            indices_new_dets(num_to_move)=iwalk
          endif
        endif

      endif

    enddo ! iwalk

    ! Now apply initiator criterion to all remaining ct's and outside_ct's that we haven't scanned over yet
   !if (.not.c_t_initiator) then
   !  do i_ct=start_ct,ndet_psi_t_connected
   !    call check_initiator(i_ct,discard_walker,i_permanent_initiator)
   !  enddo
   !else
   !  do i_ct=1,n_permanent_initiator
   !    call check_initiator(i_ct,discard_walker,i_permanent_initiator)
   !  enddo
   !endif
!***b
!    do i_outside_ct=start_outside_ct,my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)
    do i_outside_ct=start_outside_ct,ndet_psi_t_connected+ndet_outside_ct
!***e
      call check_initiator(i_outside_ct,discard_walker,i_permanent_initiator)
      if (discard_walker) then
        leave_out=leave_out+1
        leave_out_list(leave_out)=i_outside_ct
      endif
    enddo

  endif

  ! First remove the discarded walkers from outside_ct. Ideally we would like to combine this with the next
  ! step (sorting outside_ct and newly spawned walkers), but keeping them separate is simplest.

  if (leave_out>0) then
   !write (6,*) "leave_out=",leave_out
    start_empty_outside_ct=1
    do i=1,leave_out
      first = leave_out_list(i)+1
      if (i<leave_out) then
        last = leave_out_list(i+1)-1
      else
!***b
!        last = my_ndet_outside_ct(whoami+1)
        last = ndet_outside_ct
!***e
      endif
      shift = -i
      walk_dets_up(first+shift:last+shift)=walk_dets_up(first:last)
      walk_dets_dn(first+shift:last+shift)=walk_dets_dn(first:last)
      walk_wt(first+shift:last+shift)=walk_wt(first:last)
      initiator(first+shift:last+shift)=initiator(first:last)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(first+shift:last+shift)=e_num_walker(first:last)
        e_den_walker(first+shift:last+shift)=e_den_walker(first:last)
      endif
      imp_distance(first+shift:last+shift)=imp_distance(first:last)
      matrix_elements(first+shift:last+shift)=matrix_elements(first:last)

      ! Update empty_outside_ct
      do j=start_empty_outside_ct,num_to_move
        if (empty_outside_ct(j)<first) then
          cycle
        elseif (empty_outside_ct(j)>=first.and.empty_outside_ct(j)<=last) then
          empty_outside_ct(j) = empty_outside_ct(j)-i
        else
          start_empty_outside_ct = j
          exit
        endif
      enddo
    enddo

    ! Update ndet_outside_ct
!***b
!    my_ndet_outside_ct(whoami+1)=my_ndet_outside_ct(whoami+1)-leave_out
    ndet_outside_ct=ndet_outside_ct-leave_out
!***e
  endif

  ! Now sort outside_ct and newly spawned dets into one list
  if (num_to_move>0) then
    ! First move all the newly spawned dets to end of list so we don't overwrite them
    do i=num_to_move,1,-1
      from = indices_new_dets(i)
      to = MWALK+i-num_to_move
      if (indices_new_dets(i).ge.MWALK+1-num_to_move) then
        write (6,*) "ndet_psi_t_connected=",ndet_psi_t_connected,"num_to_move=",num_to_move,"MWALK=",MWALK
        write (6,*) "STOP: Need to set MWALK higher!"; call flush(6)
        stop 'Need to set MWALK higher!'
      endif
      walk_dets_up(to)=walk_dets_up(from)
      walk_dets_dn(to)=walk_dets_dn(from)
      walk_wt(to)=walk_wt(from)
      initiator(to)=initiator(from)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(to)=e_num_walker(from)
        e_den_walker(to)=e_den_walker(from)
      endif
      imp_distance(to)=imp_distance(from)
      matrix_elements(to)=matrix_elements(from)
    enddo
    ! Next make room for these dets
    do i=num_to_move,1,-1
      shift = i
      first = empty_outside_ct(i)
      if (i==num_to_move) then
!***b
!        last = my_ndet_psi_t_connected(whoami+1)+my_ndet_outside_ct(whoami+1)
        last = ndet_psi_t_connected+ndet_outside_ct
!***e
      else
        last = empty_outside_ct(i+1)-1
      endif
      walk_dets_up(first+shift:last+shift)=walk_dets_up(first:last)
      walk_dets_dn(first+shift:last+shift)=walk_dets_dn(first:last)
      walk_wt(first+shift:last+shift)=walk_wt(first:last)
      initiator(first+shift:last+shift)=initiator(first:last)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(first+shift:last+shift)=e_num_walker(first:last)
        e_den_walker(first+shift:last+shift)=e_den_walker(first:last)
      endif
      imp_distance(first+shift:last+shift)=imp_distance(first:last)
      matrix_elements(first+shift:last+shift)=matrix_elements(first:last)
    enddo
    ! Finally insert the newly spawned dets into these empty spots
    do i=1,num_to_move
      from = MWALK+i-num_to_move
      to = empty_outside_ct(i)+i-1
      walk_dets_up(to)=walk_dets_up(from)
      walk_dets_dn(to)=walk_dets_dn(from)
      walk_wt(to)=walk_wt(from)
      initiator(to)=initiator(from)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(to)=e_num_walker(from)
        e_den_walker(to)=e_den_walker(from)
      endif
      imp_distance(to)=imp_distance(from)
      matrix_elements(to)=matrix_elements(from)
    enddo

  endif

!***b
!  my_ndet_outside_ct(whoami+1) = my_ndet_outside_ct(whoami+1) + num_to_move
!  my_nwalk = my_ndet_psi_t_connected(whoami+1) + my_ndet_outside_ct(whoami+1)
!  nwalk=my_nwalk
!  call allred(nwalk) !**update total population count
!  if (hamiltonian_type.eq.'hubbard2'.or.run_type.eq.'sr') then
!    e_num_walker(my_nwalk+1:MWALK)=1.e51_rk
!    e_den_walker(my_nwalk+1:MWALK)=1.e51_rk
!  endif
!  matrix_elements(my_nwalk+1:MWALK)=1.e51_rk
  ndet_outside_ct = ndet_outside_ct + num_to_move
  nwalk = ndet_psi_t_connected + ndet_outside_ct
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(nwalk+1:MWALK)=1.e51_rk
    e_den_walker(nwalk+1:MWALK)=1.e51_rk
  endif
  matrix_elements(nwalk+1:MWALK)=1.e51_rk
!***e

 !if(r_initiator.ge.0 .and. i_permanent_initiator/=n_permanent_initiator) then
 !  write(6,'(''i_permanent_initiator, n_permanent_initiator='',2i5)') i_permanent_initiator, n_permanent_initiator
 !  write(6,'(''walk_dets_up='',1000i20)') walk_dets_up(1:nwalk)
 !  write(6,'(''walk_dets_dn='',1000i20)') walk_dets_dn(1:nwalk)
 !  write(6,'(''walk_wt     ='',1000f20.3)') walk_wt(1:nwalk)
 !  write(6,'(''initiator   ='',1000i20)') initiator(1:nwalk)
 !  stop 'i_permanent_initiator /= n_permanent_initiator'
 !endif

  if(ipr.ge.1) call write_walker_info('End of merge_original_with_spawned3,',nwalk)

! if (my_ndet_outside_ct(whoami)<0)  stop 'ndet_outside_ct<0'
  if (ndet_outside_ct<0)  stop 'ndet_outside_ct<0'

  end subroutine merge_original_with_spawned3
!-------------------------------------------------------------------------------------

! ==============================================================================
  subroutine merge_my_original_with_spawned3 (my_nwalk,my_ndet_psi_t_connected,my_ndet_outside_ct)!(walk_dets_up,walk_dets_dn,walk_wt,nwalk)
! ==============================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified by Adam Holmes to put in semistochastic projection
! Modified : A Holmes, 28 Jan 2013. First add the newly spawned dets into the connections to psi trial, then merges them with the list outside the connections from the previous iterations.
! ------------------------------------------------------------------------------
  use types, only: ik,i8b
  use common_psi_t, only : hf_to_psit
  use common_ham, only : hamiltonian_type
  use common_run, only : ndet_outside_ct

  implicit none
  integer,intent(inout) :: my_nwalk,my_ndet_outside_ct
  integer,intent(in)    :: my_ndet_psi_t_connected


 !integer, intent(in) :: iblk,istep
 !type(ik_vec), intent(inout) :: walk_dets_up(nwalk),walk_dets_dn(nwalk)
 !integer, intent(inout) :: nwalk
 !real(rk), intent(inout) :: walk_wt(nwalk)
  integer iwalk, i_permanent_initiator
  integer :: i_ct,start_ct,i_outside_ct,start_outside_ct,start_empty_outside_ct
  logical :: found,discard_walker,merge_with_prev,reorder
  integer :: i,j
  integer(i8b) :: from,to,first,last,shift

 !write (6,*) "merge_original_with_spawned3 called; nwalk=",nwalk,", ndet_psi_t_connected=",ndet_psi_t_connected,", ndet_outside_ct=",ndet_outside_ct
 !call flush(6)
  if (my_nwalk<my_ndet_psi_t_connected+my_ndet_outside_ct) then
    stop 'nwalk<ndet_psi_t_connected+ndet_outside_ct'
  endif
  i_permanent_initiator=0
  leave_out=0 ! Number of dets in outside_ct that will be discarded
  num_to_move=0 ! Number of newly spawned dets to insert into the list outside of ct
  merge_with_prev=.false. ! Whether a walker that is not in ct or outside_ct needs to be merged with the one to its left

  if (my_nwalk==my_ndet_psi_t_connected+my_ndet_outside_ct) then ! Already sorted, but still need to apply initiator criterion to all dets
   !if (.not.c_t_initiator) then
   !  do i_ct=1,ndet_psi_t_connected
   !    call check_initiator(i_ct,discard_walker,i_permanent_initiator) ! Don't discard walkers here because we keep all dets in ct.
   !  enddo
   !else ! This is necessary only to allow us to set w_perm_initiator to 1 if it is too small.
   !  do i_ct=1,n_permanent_initiator
   !    call check_initiator(i_ct,discard_walker,i_permanent_initiator)
   !  enddo
   !endif
    do i_outside_ct=my_ndet_psi_t_connected+1,my_ndet_psi_t_connected+my_ndet_outside_ct
      call check_initiator(i_outside_ct,discard_walker,i_permanent_initiator)
      if (discard_walker) then
        leave_out=leave_out+1
        leave_out_list(leave_out)=i_outside_ct
      endif
    enddo

  else

    start_ct=1
    start_outside_ct=my_ndet_psi_t_connected+1
    do iwalk=my_ndet_psi_t_connected+my_ndet_outside_ct+1,my_nwalk ! Scan over list of (sorted) newly spawned walkers
      found=.false.
      reorder=.false.
      ! First check to see if it's in list of ct
      do i_ct=start_ct,my_ndet_psi_t_connected
        if(walk_dets_up(iwalk).eq.walk_dets_up(i_ct) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(i_ct)) then ! This walker is on same det as previous walker

          if(walk_wt(iwalk)*walk_wt(i_ct).gt.0) then           ! Same sign walkers being combined
            initiator(i_ct)=max(initiator(i_ct),initiator(iwalk))
          else                                                         ! Different sign
            if(abs(walk_wt(i_ct)).lt.abs(walk_wt(iwalk))) then ! Different sign and sum has same sign as the second of these walkers.
              if(initiator(i_ct).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to that of the second walker
                initiator(i_ct)=initiator(iwalk)
              endif
            elseif(abs(walk_wt(i_ct)).eq.abs(walk_wt(iwalk))) then  ! Different sign and they add to zero
              if(initiator(i_ct).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to the smaller of the two or equivalently to 0 since the wt is zero anyway
                initiator(i_ct)=0
              endif
            endif
          endif
          if (.not.(imp_distance(i_ct).eq.0.and.imp_distance(iwalk).eq.-1))  walk_wt(i_ct)=walk_wt(i_ct)+walk_wt(iwalk) ! Combine wt only if parent and child det are not both in deterministic space

          found=.true.
          start_ct=i_ct ! not i_ct + 1 because there could be duplicates in the newly spawned dets
          exit

        elseif (walk_dets_up(iwalk)>walk_dets_up(i_ct) .or. (walk_dets_up(iwalk)==walk_dets_up(i_ct).and.walk_dets_dn(iwalk)>walk_dets_dn(i_ct))) then ! This walker (i_ct) in ct is not being stochastically spawned onto anymore. Just check initiator criterion.

         !if ((.not.c_t_initiator).or.i_ct<n_permanent_initiator) then
         !  call check_initiator(i_ct,discard_walker,i_permanent_initiator) ! Don't discard walkers here because we keep all dets in ct.
         !endif
          start_ct=i_ct+1
          ! Would like to be able to calculate energy for psi_t_connections here, but we have to multiply by S^-1 first.

        else ! No reason to keep scanning ct, since there won't be any more matches.

          start_ct = i_ct
          exit

        endif

      enddo ! i_ct

      if (found)  cycle

      ! Now loop over the outside_ct dets from the previous MC iteration
      do i_outside_ct=start_outside_ct,my_ndet_psi_t_connected+my_ndet_outside_ct
        if(walk_dets_up(iwalk).eq.walk_dets_up(i_outside_ct) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(i_outside_ct)) then ! This walker is on same det as previous walker

          if(walk_wt(iwalk)*walk_wt(i_outside_ct).gt.0) then           ! Same sign walkers being combined
            initiator(i_outside_ct)=max(initiator(i_outside_ct),initiator(iwalk))
            imp_distance(i_outside_ct)=min(imp_distance(i_outside_ct),abs(imp_distance(iwalk)))
          else                                                         ! Different sign
            imp_distance(i_outside_ct)=min(imp_distance(i_outside_ct),abs(imp_distance(iwalk)))
            if(abs(walk_wt(i_outside_ct)).lt.abs(walk_wt(iwalk))) then ! Different sign and sum has same sign as the second of these walkers.
              if(initiator(i_outside_ct).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to that of the second walker
                initiator(i_outside_ct)=initiator(iwalk)
              endif
            elseif(abs(walk_wt(i_outside_ct)).eq.abs(walk_wt(iwalk))) then  ! Different sign and they add to zero
              if(initiator(i_outside_ct).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to the smaller of the two or equivalently to 0 since the wt is zero anyway
                initiator(i_outside_ct)=0
              endif
            endif
          endif
          walk_wt(i_outside_ct)=walk_wt(i_outside_ct)+walk_wt(iwalk)

          found=.true.
          start_outside_ct=i_outside_ct ! not i_outside_ct + 1 because there could be duplicates
          exit

        elseif (walk_dets_up(i_outside_ct)<walk_dets_up(iwalk).or.(walk_dets_up(i_outside_ct)==walk_dets_up(iwalk).and.walk_dets_dn(i_outside_ct)<walk_dets_dn(iwalk))) then ! This walker (i_outside_ct) in outside ct is not being stochastically spawned onto anymore. Just check initiator criterion.

          call check_initiator(i_outside_ct,discard_walker,i_permanent_initiator)
          if (discard_walker) then
            leave_out=leave_out+1
            leave_out_list(leave_out)=i_outside_ct
          endif
          start_outside_ct = i_outside_ct+1

        else ! No reason to keep scanning outside_ct, since there won't be any more matches.
        ! We passed the slot where we want to insert det iwalk in the list of saved outside_ct dets
          start_outside_ct = i_outside_ct
          reorder = .true.
          exit

        endif

      enddo ! i_outside_ct

      if (.not.found) then

        if (merge_with_prev) then
          ! combine iwalk with iwalk-1, store at iwalk

          if(walk_wt(iwalk-1)*walk_wt(iwalk).gt.0) then           ! Same sign walkers being combined
            initiator(iwalk)=max(initiator(iwalk),initiator(iwalk-1))
            if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
              e_num_walker(iwalk)=min(e_num_walker(iwalk),e_num_walker(iwalk-1))         ! HJC: If calculated energy, then use minimum value as default is positive and set to 1.e51_rk
              e_den_walker(iwalk)=min(e_den_walker(iwalk),e_den_walker(iwalk-1))
            endif
            imp_distance(iwalk)=min(imp_distance(iwalk),abs(imp_distance(iwalk-1)))
            matrix_elements(iwalk)=min(matrix_elements(iwalk),(matrix_elements(iwalk-1)))
          else                                                         ! Different sign
            if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
              e_num_walker(iwalk)=min(e_num_walker(iwalk),e_num_walker(iwalk-1))         ! HJC : If calculated energy, then use minimum value as default is positive and set to 1.e51_rk
              e_den_walker(iwalk)=min(e_den_walker(iwalk),e_den_walker(iwalk-1))
            endif
            imp_distance(iwalk)=min(imp_distance(iwalk),abs(imp_distance(iwalk-1)))
            matrix_elements(iwalk)=min(matrix_elements(iwalk),(matrix_elements(iwalk-1)))
            if(abs(walk_wt(iwalk)).lt.abs(walk_wt(iwalk-1))) then ! Different sign and sum has same sign as the second of these walkers.
              if(initiator(iwalk).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to that of the second walker
                initiator(iwalk)=initiator(iwalk-1)
              endif
            elseif(abs(walk_wt(iwalk)).eq.abs(walk_wt(iwalk-1))) then  ! Different sign and they add to zero
              if(initiator(iwalk).ne.3 .or. r_initiator.eq.-1._rk) then  ! Provided it is not a permanent initiator, set the initiator to the smaller of the two or equivalently to 0 since the wt is zero anyway
                initiator(iwalk)=0
              endif
            endif
          endif
          walk_wt(iwalk)=walk_wt(iwalk)+walk_wt(iwalk-1)
          merge_with_prev=.false.

        endif

        if (iwalk<my_nwalk.and.walk_dets_up(iwalk)==walk_dets_up(iwalk+1).and.walk_dets_dn(iwalk)==walk_dets_dn(iwalk+1)) then ! Try to combine this with other newly spawned walkers
          merge_with_prev=.true.
        else
          ! Check initiator criterion on det iwalk.
          call check_initiator(iwalk,discard_walker,i_permanent_initiator)
          if (.not.discard_walker) then
            ! Add to list of new walkers that we keep
            num_to_move=num_to_move+1
            if (reorder) then
              empty_outside_ct(num_to_move)=i_outside_ct
            else
              ! This det is not in ct or outside_ct and should be added to the end.
              empty_outside_ct(num_to_move)=my_ndet_psi_t_connected+my_ndet_outside_ct+1
            endif
            indices_new_dets(num_to_move)=iwalk
          endif
        endif

      endif

    enddo ! iwalk

    ! Now apply initiator criterion to all remaining ct's and outside_ct's that we haven't scanned over yet
   !if (.not.c_t_initiator) then
   !  do i_ct=start_ct,ndet_psi_t_connected
   !    call check_initiator(i_ct,discard_walker,i_permanent_initiator)
   !  enddo
   !else
   !  do i_ct=1,n_permanent_initiator
   !    call check_initiator(i_ct,discard_walker,i_permanent_initiator)
   !  enddo
   !endif
    do i_outside_ct=start_outside_ct,my_ndet_psi_t_connected+my_ndet_outside_ct
      call check_initiator(i_outside_ct,discard_walker,i_permanent_initiator)
      if (discard_walker) then
        leave_out=leave_out+1
        leave_out_list(leave_out)=i_outside_ct
      endif
    enddo

  endif

  ! First remove the discarded walkers from outside_ct. Ideally we would like to combine this with the next
  ! step (sorting outside_ct and newly spawned walkers), but keeping them separate is simplest.

  if (leave_out>0) then
   !write (6,*) "leave_out=",leave_out
    start_empty_outside_ct=1
    do i=1,leave_out
      first = leave_out_list(i)+1
      if (i<leave_out) then
        last = leave_out_list(i+1)-1
      else
        last = my_ndet_outside_ct
      endif
      shift = -i
      walk_dets_up(first+shift:last+shift)=walk_dets_up(first:last)
      walk_dets_dn(first+shift:last+shift)=walk_dets_dn(first:last)
      walk_wt(first+shift:last+shift)=walk_wt(first:last)
      initiator(first+shift:last+shift)=initiator(first:last)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(first+shift:last+shift)=e_num_walker(first:last)
        e_den_walker(first+shift:last+shift)=e_den_walker(first:last)
      endif
      imp_distance(first+shift:last+shift)=imp_distance(first:last)
      matrix_elements(first+shift:last+shift)=matrix_elements(first:last)

      ! Update empty_outside_ct
      do j=start_empty_outside_ct,num_to_move
        if (empty_outside_ct(j)<first) then
          cycle
        elseif (empty_outside_ct(j)>=first.and.empty_outside_ct(j)<=last) then
          empty_outside_ct(j) = empty_outside_ct(j)-i
        else
          start_empty_outside_ct = j
          exit
        endif
      enddo
    enddo

    ! Update ndet_outside_ct
    my_ndet_outside_ct=my_ndet_outside_ct-leave_out
    ndet_outside_ct=ndet_outside_ct-leave_out
  endif

  ! Now sort outside_ct and newly spawned dets into one list
  if (num_to_move>0) then
    ! First move all the newly spawned dets to end of list so we don't overwrite them
    do i=num_to_move,1,-1
      from = indices_new_dets(i)
      to = MWALK+i-num_to_move
      if (indices_new_dets(i).ge.MWALK+1-num_to_move) then
        write (6,*) "ndet_psi_t_connected=",ndet_psi_t_connected,"num_to_move=",num_to_move,"MWALK=",MWALK
        write (6,*) "STOP: Need to set MWALK higher!"; call flush(6)
        stop 'Need to set MWALK higher!'
      endif
      walk_dets_up(to)=walk_dets_up(from)
      walk_dets_dn(to)=walk_dets_dn(from)
      walk_wt(to)=walk_wt(from)
      initiator(to)=initiator(from)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(to)=e_num_walker(from)
        e_den_walker(to)=e_den_walker(from)
      endif
      imp_distance(to)=imp_distance(from)
      matrix_elements(to)=matrix_elements(from)
    enddo
    ! Next make room for these dets
    do i=num_to_move,1,-1
      shift = i
      first = empty_outside_ct(i)
      if (i==num_to_move) then
        last = my_ndet_psi_t_connected+my_ndet_outside_ct
      else
        last = empty_outside_ct(i+1)-1
      endif
      walk_dets_up(first+shift:last+shift)=walk_dets_up(first:last)
      walk_dets_dn(first+shift:last+shift)=walk_dets_dn(first:last)
      walk_wt(first+shift:last+shift)=walk_wt(first:last)
      initiator(first+shift:last+shift)=initiator(first:last)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(first+shift:last+shift)=e_num_walker(first:last)
        e_den_walker(first+shift:last+shift)=e_den_walker(first:last)
      endif
      imp_distance(first+shift:last+shift)=imp_distance(first:last)
      matrix_elements(first+shift:last+shift)=matrix_elements(first:last)
    enddo
    ! Finally insert the newly spawned dets into these empty spots
    do i=1,num_to_move
      from = MWALK+i-num_to_move
      to = empty_outside_ct(i)+i-1
      walk_dets_up(to)=walk_dets_up(from)
      walk_dets_dn(to)=walk_dets_dn(from)
      walk_wt(to)=walk_wt(from)
      initiator(to)=initiator(from)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(to)=e_num_walker(from)
        e_den_walker(to)=e_den_walker(from)
      endif
      imp_distance(to)=imp_distance(from)
      matrix_elements(to)=matrix_elements(from)
    enddo

  endif

  my_ndet_outside_ct = my_ndet_outside_ct + num_to_move
  my_nwalk = my_ndet_psi_t_connected + my_ndet_outside_ct
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(my_nwalk+1:MWALK)=1.e51_rk
    e_den_walker(my_nwalk+1:MWALK)=1.e51_rk
  endif
  matrix_elements(my_nwalk+1:MWALK)=1.e51_rk

 !if(r_initiator.ge.0 .and. i_permanent_initiator/=n_permanent_initiator) then
 !  write(6,'(''i_permanent_initiator, n_permanent_initiator='',2i5)') i_permanent_initiator, n_permanent_initiator
 !  write(6,'(''walk_dets_up='',1000i20)') walk_dets_up(1:nwalk)
 !  write(6,'(''walk_dets_dn='',1000i20)') walk_dets_dn(1:nwalk)
 !  write(6,'(''walk_wt     ='',1000f20.3)') walk_wt(1:nwalk)
 !  write(6,'(''initiator   ='',1000i20)') initiator(1:nwalk)
 !  stop 'i_permanent_initiator /= n_permanent_initiator'
 !endif

  if(ipr.ge.1) call write_walker_info('End of merge_my_original_with_spawned2,',my_nwalk)

  if (my_ndet_outside_ct<0)  stop 'ndet_outside_ct<0'

end subroutine merge_my_original_with_spawned3
!-------------------------------------------------------------------------------------


! ==============================================================================
  subroutine check_initiator(iwalk,discard_walker,i_permanent_initiator)
! ==============================================================================
! Description   : Check initiator criterion, i.e.:
!               : (1) Make sure permanent initiator has abs wt of at least 1 with right sign
!               : (2) Increase or decrease initiator value if wt crosses initiator threshold
!               : (3) Output discard_walker, whether or not to discard this walker because
!               :     it was spawned by a noninitiator and has abs wt below threshold
! Author        : A Holmes, 31 Jan 2013.
! ------------------------------------------------------------------------------

    integer,intent(in) :: iwalk
    logical,intent(out) :: discard_walker ! whether or not to discard this iwalk (if it was spawned by a noninitiator and has abs wt below threshold)
    integer,intent(inout) :: i_permanent_initiator ! update this count of number of permanent initiators (gives an error if not what we expect)

      discard_walker=.false.
      if(initiator(iwalk).eq.3 .and. r_initiator.ge.0) then    ! Make sure that permanent initiator has an abs wt of atleast 1 with the right sign
        i_permanent_initiator=i_permanent_initiator+1
        if(walk_wt(iwalk)*sign_permanent_initiator(i_permanent_initiator) < 1._rk) then
          if(n_warning.le.M_WARNING) write(6,'(''Warning: iblk, istep'',2i6,'' walk_wt of permanent initiator changed from'',f4.1,'' to'',i3)') iblk, istep, walk_wt(iwalk), sign_permanent_initiator(i_permanent_initiator)
          if(n_warning.eq.M_WARNING) write(6,'(''Warning: n_warning = M_WARNING, so no more warnings will be issued'')')
          n_warning=n_warning+1
          walk_wt(iwalk)=sign_permanent_initiator(i_permanent_initiator)
        endif
      elseif(initiator(iwalk).eq.2 .and. ((abs(walk_wt(iwalk)).le.r_initiator*(max(0,imp_distance(iwalk)-initiator_min_distance)**initiator_power).and.imp_distance(iwalk)>0) .or. ((abs(walk_wt(iwalk)).le.r_initiator.and..not.c_t_initiator).and.imp_distance(iwalk)==-2))) then ! If it was an initiator but not a permanent initiator and the abs wt goes below threshold make it a noninitiator
          if(ipr.ge.1) write(6,'(''Remove from initiators1'',9i6)') walk_dets_up(iwalk), walk_dets_dn(iwalk), initiator(iwalk)
          initiator(iwalk)=1
      elseif(initiator(iwalk).lt.2 .and. ((abs(walk_wt(iwalk)).gt.r_initiator*(max(0,imp_distance(iwalk)-initiator_min_distance)**initiator_power).and.imp_distance(iwalk)>=0) .or. ((abs(walk_wt(iwalk)).gt.r_initiator.or.c_t_initiator).and.imp_distance(iwalk)==-2))) then ! If it was not an initiator and the abs wt goes above threshold boost its initiator value by 1
          if(ipr.ge.1) write(6,'(''For det'',2i8,'' incrementing initiator from'',i3,'' to'',i3)') walk_dets_up(iwalk), walk_dets_dn(iwalk), initiator(iwalk), initiator(iwalk)+1
          initiator(iwalk)=initiator(iwalk)+1_i1b
      endif
      if(((walk_wt(iwalk).eq.0 .and. (initiator(iwalk).ne.3 .or. r_initiator.lt.0)) .or. initiator(iwalk).eq.0) .and. imp_distance(iwalk).ge.1) then
        discard_walker=.true. !nshift=nshift+1 ! If the sum of wts for the previous det is 0 and it is not a permanent initiator, or, if it was spawned from a noninitiator and the abs wt is below threshold then one can reuse this location
      endif

  end subroutine check_initiator
! ------------------------------------------------------------------------------


  subroutine join_walker(my_nwalk)
! ==============================================================================
! Description   : Join walkers with small absolute weights until the absolute wt is larger than min_wt
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use common_ham, only : hamiltonian_type
  implicit none
  integer , intent(inout) :: my_nwalk
  integer iwalk, iwalk2, ipair, nshift
  real(rk) rannyu, wttot

! if(ipr.ge.0) write(6,'(''walk_wt_in ='',1000f8.3)') (walk_wt(iwalk),iwalk=1,my_nwalk)

! Join low +ve wt walkers in pairs
  ipair=0
  do iwalk=1,my_nwalk
!   if(abs(walk_wt(iwalk)).lt.min_wt .and. initiator(iwalk).lt.3) then
    if(walk_wt(iwalk).gt.0._rk .and. abs(walk_wt(iwalk)).lt.min_wt .and. initiator(iwalk).lt.3) then
      if(ipair.eq.0) then
        ipair=1
        iwalk2=iwalk
      else
        wttot=abs(walk_wt(iwalk))+abs(walk_wt(iwalk2))
        if(rannyu().gt.(abs(walk_wt(iwalk))/wttot)) then
          walk_wt(iwalk2)=sign(wttot,walk_wt(iwalk2))
          walk_wt(iwalk)=0
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk)=1.e51_rk
            e_den_walker(iwalk)=1.e51_rk
          endif
          matrix_elements(my_nwalk+1:MWALK)=1.e51_rk
        else
          walk_wt(iwalk)=sign(wttot,walk_wt(iwalk))
          walk_wt(iwalk2)=0
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk2)=1.e51_rk
            e_den_walker(iwalk2)=1.e51_rk
          endif
          matrix_elements(my_nwalk+1:MWALK)=1.e51_rk
          iwalk2=iwalk ! used only if also wttot>min_wt
        endif
        if(wttot>min_wt) then
          ipair=0
        endif
      endif
    endif
  enddo

! Join low -ve wt walkers in pairs
  ipair=0
  do iwalk=1,my_nwalk
    if(walk_wt(iwalk).lt.0_rk .and. abs(walk_wt(iwalk)).lt.min_wt .and. initiator(iwalk).lt.3) then
      if(ipair.eq.0) then
        ipair=1
        iwalk2=iwalk
      else
        wttot=abs(walk_wt(iwalk))+abs(walk_wt(iwalk2))
        if(rannyu().gt.(abs(walk_wt(iwalk))/wttot)) then
          walk_wt(iwalk2)=sign(wttot,walk_wt(iwalk2))
          walk_wt(iwalk)=0
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk)=1.e51_rk
            e_den_walker(iwalk)=1.e51_rk
          endif
          matrix_elements(my_nwalk+1:MWALK)=1.e51_rk
        else
          walk_wt(iwalk)=sign(wttot,walk_wt(iwalk))
          walk_wt(iwalk2)=0
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk2)=1.e51_rk
            e_den_walker(iwalk2)=1.e51_rk
          endif
          matrix_elements(my_nwalk+1:MWALK)=1.e51_rk
          iwalk2=iwalk ! used only if also wttot>min_wt
        endif
        if(wttot>min_wt) then
          ipair=0
        endif
      endif
    endif
  enddo

! Eliminate zero-wt walkers
  nshift=0
  do iwalk=1,my_nwalk
    if(walk_wt(iwalk).eq.0_rk .and. imp_distance(iwalk).ge.1) then
      nshift=nshift+1
    else
      walk_wt(iwalk-nshift)=walk_wt(iwalk)
      walk_det(iwalk-nshift)=walk_det(iwalk)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(iwalk-nshift)=e_num_walker(iwalk)
        e_den_walker(iwalk-nshift)=e_den_walker(iwalk)
      endif
      initiator(iwalk-nshift)=initiator(iwalk)
      imp_distance(iwalk-nshift)=imp_distance(iwalk)
      matrix_elements(iwalk-nshift)=matrix_elements(iwalk)
    endif
  enddo
  my_nwalk=my_nwalk-nshift
! imp_distance(my_nwalk+1:MWALK)=1_i1b
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(my_nwalk+1:MWALK)=1.e51_rk
    e_den_walker(my_nwalk+1:MWALK)=1.e51_rk
  endif
  matrix_elements(my_nwalk+1:MWALK)=1.e51_rk

! if(ipr.ge.0) write(6,'(''walk_wt_out='',1000f8.3)') (walk_wt(iwalk),iwalk=1,my_nwalk)

  if(ipr.ge.1) call write_walker_info('End of join_walker,',my_nwalk)

  end subroutine join_walker
!-------------------------------------------------------------------------------------

  subroutine join_walker2(my_nwalk)
! ==============================================================================
! Description   : Join walkers with small absolute weights until the absolute wt is larger than min_wt
! Author        : Cyrus Umrigar
! ------------------------------------------------------------------------------
  use common_ham, only : hamiltonian_type
  implicit none
  integer , intent(inout) :: my_nwalk
  integer iwalk, iwalk2, ipair, nshift
  real(rk) rannyu, wttot

! if(ipr.ge.0) write(6,'(''walk_wt_in ='',1000f8.3)') (walk_wt(iwalk),iwalk=1,my_nwalk)

! Join low +ve wt walkers in pairs
  ipair=0
  do iwalk=1,my_nwalk
    if(walk_wt(iwalk).gt.0._rk .and. abs(walk_wt(iwalk)).lt.min_wt .and. initiator(iwalk).lt.3) then
      if(ipair.eq.0) then
        ipair=1
        iwalk2=iwalk
      else
        wttot=abs(walk_wt(iwalk))+abs(walk_wt(iwalk2))
        if(rannyu().gt.(abs(walk_wt(iwalk))/wttot)) then
          walk_wt(iwalk2)=sign(wttot,walk_wt(iwalk2))
          walk_wt(iwalk)=0
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk)=1.e51_rk
            e_den_walker(iwalk)=1.e51_rk
          endif
          matrix_elements(my_nwalk+1:MWALK)=1.e51_rk
        else
          walk_wt(iwalk)=sign(wttot,walk_wt(iwalk))
          walk_wt(iwalk2)=0
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk2)=1.e51_rk
            e_den_walker(iwalk2)=1.e51_rk
          endif
          matrix_elements(my_nwalk+1:MWALK)=1.e51_rk
          iwalk2=iwalk ! used only if also wttot>min_wt
        endif
        if(wttot>min_wt) then
          ipair=0
        endif
      endif
    endif
  enddo

! Join low -ve wt walkers in pairs
  ipair=0
  do iwalk=1,my_nwalk
    if(walk_wt(iwalk).lt.0._rk .and. abs(walk_wt(iwalk)).lt.min_wt .and. initiator(iwalk).lt.3) then
      if(ipair.eq.0) then
        ipair=1
        iwalk2=iwalk
      else
        wttot=abs(walk_wt(iwalk))+abs(walk_wt(iwalk2))
        if(rannyu().gt.(abs(walk_wt(iwalk))/wttot)) then
          walk_wt(iwalk2)=sign(wttot,walk_wt(iwalk2))
          walk_wt(iwalk)=0
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk)=1.e51_rk
            e_den_walker(iwalk)=1.e51_rk
          endif
          matrix_elements(my_nwalk+1:MWALK)=1.e51_rk
        else
          walk_wt(iwalk)=sign(wttot,walk_wt(iwalk))
          walk_wt(iwalk2)=0
          if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
            e_num_walker(iwalk2)=1.e51_rk
            e_den_walker(iwalk2)=1.e51_rk
          endif
          matrix_elements(my_nwalk+1:MWALK)=1.e51_rk
          iwalk2=iwalk ! used only if also wttot>min_wt
        endif
        if(wttot>min_wt) then
          ipair=0
        endif
      endif
    endif
  enddo

! Eliminate zero-wt walkers
  nshift=0
  do iwalk=1,my_nwalk
    !if(walk_wt(iwalk).eq.0_rk .and. imp_distance(iwalk).ge.1) then - Do we need imp_distance for stochastic??
    if(walk_wt(iwalk).eq.0_rk) then
      nshift=nshift+1
    else
      walk_wt(iwalk-nshift)=walk_wt(iwalk)
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      e_num_walker(iwalk-nshift)=e_num_walker(iwalk)
      e_den_walker(iwalk-nshift)=e_den_walker(iwalk)
      initiator(iwalk-nshift)=initiator(iwalk)
      if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
        e_num_walker(iwalk-nshift)=e_num_walker(iwalk)
        e_den_walker(iwalk-nshift)=e_den_walker(iwalk)
      endif
      ! I am uncommenting the following line because I don't see why it was commented out. - AAH, 18 Jun 2013.
      imp_distance(iwalk-nshift)=imp_distance(iwalk)
      matrix_elements(iwalk-nshift)=matrix_elements(iwalk)
    endif
  enddo
  my_nwalk=my_nwalk-nshift
! imp_distance(my_nwalk+1:MWALK)=1_i1b
  if (hamiltonian_type.eq.'hubbard2'.or.hamiltonian_type.eq.'hubbardk'.or.hamiltonian_type.eq.'chem'.or.run_type.eq.'sr') then
    e_num_walker(my_nwalk+1:MWALK)=1.e51_rk
    e_den_walker(my_nwalk+1:MWALK)=1.e51_rk
  endif
  matrix_elements(my_nwalk+1:MWALK)=1.e51_rk

! if(ipr.ge.0) write(6,'(''walk_wt_out='',1000f8.3)') (walk_wt(iwalk),iwalk=1,my_nwalk)

  end subroutine join_walker2
!-------------------------------------------------------------------------------------

  subroutine join_walker_semistoch(my_nwalk)
! ==============================================================================
! Description   : Join walkers with small absolute weights until the absolute wt is larger than min_wt
! Author        : Cyrus Umrigar
!               : Edited by AAH. Same as join_walker2, but avoids joining important walkers when using semi-stochastic moves
! ------------------------------------------------------------------------------

  use semistoch
  use common_run, only : ndet_outside_ct

  implicit none

  integer , intent(inout) :: my_nwalk
  integer iwalk, iwalk2, ipair, nshift
  real(rk) rannyu, wttot

! if(ipr.ge.0) write(6,'(''walk_wt_in ='',1000f8.3)') (walk_wt(iwalk),iwalk=1,my_nwalk)

! Join low +ve wt walkers in pairs
  ipair=0
  do iwalk=1,my_nwalk ! iwalk=ndet_psi_t_connected+1,my_nwalk
    if(walk_wt(iwalk).gt.0._rk .and. abs(walk_wt(iwalk)).lt.min_wt .and. initiator(iwalk).lt.3 .and. (imp_distance(iwalk).ge.1)) then
      if(ipair.eq.0) then
        ipair=1
        iwalk2=iwalk
      else
        wttot=abs(walk_wt(iwalk))+abs(walk_wt(iwalk2))
        if(rannyu().gt.(abs(walk_wt(iwalk))/wttot)) then
          walk_wt(iwalk2)=sign(wttot,walk_wt(iwalk2))
          walk_wt(iwalk)=0
        else
          walk_wt(iwalk)=sign(wttot,walk_wt(iwalk))
          walk_wt(iwalk2)=0
          iwalk2=iwalk ! used only if also wttot>min_wt
        endif
        if(wttot>min_wt) then
          ipair=0
        endif
      endif
    endif
  enddo

! Join low -ve wt walkers in pairs
  ipair=0
  do iwalk=1,my_nwalk
    if(walk_wt(iwalk).lt.0_rk .and. abs(walk_wt(iwalk)).lt.min_wt .and. initiator(iwalk).lt.3 .and. (imp_distance(iwalk).ge.1)) then
      if(ipair.eq.0) then
        ipair=1
        iwalk2=iwalk
      else
        wttot=abs(walk_wt(iwalk))+abs(walk_wt(iwalk2))
        if(rannyu().gt.(abs(walk_wt(iwalk))/wttot)) then
          walk_wt(iwalk2)=sign(wttot,walk_wt(iwalk2))
          walk_wt(iwalk)=0
        else
          walk_wt(iwalk)=sign(wttot,walk_wt(iwalk))
          walk_wt(iwalk2)=0
          iwalk2=iwalk ! used only if also wttot>min_wt
        endif
        if(wttot>min_wt) then
          ipair=0
        endif
      endif
    endif
  enddo

! Eliminate zero-wt walkers
  nshift=0
  do iwalk=1,my_nwalk
    if(walk_wt(iwalk).eq.0_rk .and. imp_distance(iwalk).ge.1) then
      nshift=nshift+1
    else
      walk_wt(iwalk-nshift)=walk_wt(iwalk)
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      e_num_walker(iwalk-nshift)=e_num_walker(iwalk)
      e_den_walker(iwalk-nshift)=e_den_walker(iwalk)
      initiator(iwalk-nshift)=initiator(iwalk)
      imp_distance(iwalk-nshift)=imp_distance(iwalk)
      matrix_elements(iwalk-nshift)=matrix_elements(iwalk)
    endif
  enddo
  my_nwalk=my_nwalk-nshift
  ndet_outside_ct = ndet_outside_ct-nshift

! if(ipr.ge.0) write(6,'(''walk_wt_out='',1000f8.3)') (walk_wt(iwalk),iwalk=1,my_nwalk)

  end subroutine join_walker_semistoch
!-------------------------------------------------------------------------------------

  subroutine reduce_my_walker(my_nwalk,my_ndet_psi_t_connected,my_ndet_outside_ct) !***this version is different just in that it takes a general arguments as input [MERGE WITH ORIGINAL: TODO]
! ==============================================================================
! Description   : Performs similar function to join (reduces number of low-wt walkers) but in a different way:
!               : Instead of combining low-wt walkers in pairs, each low-wt walker has a chance of being set
!               : to min_wt with probability proportional to its abs wt.
! Author        : A Holmes, 2 Apr 2013.
! ------------------------------------------------------------------------------

  use common_run, only : ndet_outside_ct
  use common_imp, only : n_imp
  use common_psi_t, only : hf_to_psit

  implicit none

  integer,intent(inout) :: my_nwalk,my_ndet_outside_ct
  integer,intent(in) :: my_ndet_psi_t_connected

  integer iwalk, nshift, i_start
  real(rk) rannyu

  if (hf_to_psit) then
    i_start = my_ndet_psi_t_connected+1
  else
    i_start = 1 ! Not n_imp+1 because the imp dets are not necessarily the first ones
  endif

  do iwalk=i_start,my_nwalk
    if (imp_distance(iwalk).ge.1) then
      if (abs(walk_wt(iwalk))<min_wt) then
        if (rannyu()<(abs(walk_wt(iwalk))/min_wt)) then
          walk_wt(iwalk) = sign(min_wt,walk_wt(iwalk))
        else
          walk_wt(iwalk) = 0._rk
        endif
      endif
    endif
  enddo

! Eliminate zero-wt walkers
  nshift=0
  do iwalk=i_start,my_nwalk
    if(walk_wt(iwalk).eq.0_rk .and. imp_distance(iwalk).ge.1) then
      nshift=nshift+1
    elseif (nshift.ne.0) then
      walk_wt(iwalk-nshift)=walk_wt(iwalk)
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      e_num_walker(iwalk-nshift)=e_num_walker(iwalk)
      e_den_walker(iwalk-nshift)=e_den_walker(iwalk)
      initiator(iwalk-nshift)=initiator(iwalk)
      imp_distance(iwalk-nshift)=imp_distance(iwalk)
      matrix_elements(iwalk-nshift)=matrix_elements(iwalk)
    endif
  enddo

  my_nwalk=my_nwalk-nshift
  my_ndet_outside_ct = my_ndet_outside_ct-nshift

  end subroutine reduce_my_walker
!-------------------------------------------------------------------------------------


subroutine stochastic_reconfiguration
! ====================================================================================

! Created by : Hitesh J. Changlani, September 17, 2012
! Purpose    : redistribute walker weights so as to make most of them positive keeping important expectation values unchanged
!              The number of such observables is p in Sorella's notation. If we keep the e_mix_numerator and e_mix denominator unchanged
!              then p = 2. This routine has been written specifically for p=2 but can be easily generalized
! ====================================================================================

use hubbard,only :energy_pieces_hubbard

implicit none

real(rk) :: avg_e_nums,avg_e_dens,avg_e_num_e_num,avg_e_den_e_den,avg_e_num_e_den,sum_fn_wts
real(rk) :: true_avg_e_nums,true_avg_e_dens,sum_wts
real(rk) :: bracket
integer  :: p=1
real(rk) :: a(1,1),diff_expec(1),alphas(1)
real(rk) :: det
integer  :: iwalk
real(rk) :: c

true_avg_e_nums=0._rk
true_avg_e_dens=0._rk
sum_wts=0._rk
avg_e_nums=0._rk
avg_e_dens=0._rk
sum_fn_wts=0._rk
avg_e_num_e_num=0._rk
avg_e_den_e_den=0._rk
avg_e_num_e_den=0._rk

! Compute avg_e_nums, avg_e_dens and more averages
if (nwalk .gt. 1) then
    do iwalk=1,nwalk
        ! Cyrus' suggestion - set walk_wt_fn(iwalk)=0                if walk_wt is negative
        !                    else walk_wt_fn(iwalk)=walk_wt(iwalk)
        if (walk_wt(iwalk) .lt. 0._rk) then
            walk_wt_fn(iwalk)=0._rk
        else
            walk_wt_fn(iwalk)=walk_wt(iwalk)
        endif
        if (e_num_walker(iwalk) .gt. 1.e50_rk .or. e_den_walker(iwalk) .gt. 1.e50_rk) then       ! The numerators and denominators may not have been computed
              call energy_pieces_hubbard(walk_dets_up(iwalk), walk_dets_dn(iwalk), e_num_walker(iwalk), e_den_walker(iwalk), importance_sampling)
        endif
        true_avg_e_nums=true_avg_e_nums+(walk_wt(iwalk)*e_num_walker(iwalk))
        true_avg_e_dens=true_avg_e_dens+(walk_wt(iwalk)*e_den_walker(iwalk))
        avg_e_nums=avg_e_nums+(walk_wt_fn(iwalk)*e_num_walker(iwalk))
        avg_e_num_e_num=avg_e_num_e_num+(walk_wt_fn(iwalk)*e_num_walker(iwalk)*e_num_walker(iwalk))
        avg_e_dens=avg_e_dens+(walk_wt_fn(iwalk)*e_den_walker(iwalk))
        avg_e_den_e_den=avg_e_den_e_den+(walk_wt_fn(iwalk)*e_den_walker(iwalk)*e_den_walker(iwalk))
        avg_e_num_e_den=avg_e_num_e_den+(walk_wt_fn(iwalk)*e_num_walker(iwalk)*e_den_walker(iwalk))
        sum_wts=sum_wts+walk_wt(iwalk)
        sum_fn_wts=sum_fn_wts+walk_wt_fn(iwalk)
    enddo

    !call print_real_matrix(1,nwalk,walk_wt)
    !call print_real_matrix(1,nwalk,walk_wt_fn)
    !write (6,*) e_num_walker(1:nwalk)
    !write (6,*) e_den_walker(1:nwalk)
    !write (6,*) "sum_wts=",sum_wts
    !write (6,*) "sum_fn_wts=",sum_fn_wts
    call flush(6)

    true_avg_e_nums=true_avg_e_nums/sum_wts
    true_avg_e_dens=true_avg_e_dens/sum_wts
    avg_e_nums=avg_e_nums/sum_fn_wts
    avg_e_dens=avg_e_dens/sum_fn_wts
    avg_e_num_e_num=avg_e_num_e_num/sum_fn_wts
    avg_e_den_e_den=avg_e_den_e_den/sum_fn_wts
    avg_e_num_e_den=avg_e_num_e_den/sum_fn_wts
    c=sum_wts/sum_fn_wts

    ! Compute - the A matrix first to get alpha(1),alpha(2)
    a(1,1)=avg_e_num_e_num-(avg_e_nums**2)
    !a(1,2)=avg_e_num_e_den-(avg_e_nums*avg_e_dens)
    !a(2,1)=a(1,2)
    !a(2,2)=avg_e_den_e_den-(avg_e_dens**2)
    diff_expec(1)=true_avg_e_nums-avg_e_nums
    !diff_expec(2)=true_avg_e_dens-avg_e_dens

    !write (6,*) "A before"
    !call print_real_matrix(p,p,a)
    ! To get these we have to solve a pxp problem
    ! Invert A: A is replaced by its inverse in matinv.f
    call matinv(a,p,det)
    ! multiply A_inv by difference of expectations
    alphas=matmul(a,diff_expec)
    !write (6,*) "A after"
    !call print_real_matrix(p,p,a)
    !write (6,*) "alphas",alphas(1),alphas(2)
    !call flush(6)

    ! Loops over all walkers and assign new weights
    do iwalk=1,nwalk
        bracket=1._rk+ (alphas(1)*(e_num_walker(iwalk)-avg_e_nums))
        !bracket=bracket+ alphas(2)*(e_den_walker(iwalk)-avg_e_dens)
        walk_wt(iwalk)=walk_wt_fn(iwalk)*bracket*c
    enddo
endif

end subroutine stochastic_reconfiguration
!-------------------------------------------------------------------------------------

subroutine stochastic_reconfiguration2(nwalk,wt,e_num,e_den)
! Do Sandro Sorella's stochastic reconfiguration keeping just the reconfigured energy denominator and numerator equal to the true value.
! This does not require solving a set of linear equations.  If we keep any more observables constant, then we do.
! In this version we are assuming that PsiG=1.  In order to not have a sign problem, we need to impose that
! wt_fn*e_num <= 0 (assuming that the expectation value of the energy is negative).
! We have 2 options here for setting wt_fn, controlled by the abs_wt logical variable.
! 1) set wt_fn=-wt*sign(e_num),                                               , if abs_wt=true
! 2) set wt_fn=-wt*sign(e_num), if it is not sign violating and zero otherwise, if abs_wt=false (better choice)
! e_den = PsiT/PsiG
! e_num = H PsiT/PsiG
! wt    = Psi0*PsiG
! wtt   = wt * e_den = Psi0*PsiT.  For a dense wavefn, if this is negative, walker is sign violating.
! This routine can be used either for
! 1) PsiG=PsiT, in which case e_den=1 (this happens for real-space Hubbard with importance sampling)
! 2) PsiG=1, in which case e_den=PsiT (this happens for orbital-space: momentum-space Hubbard, chemistry and HEG without importance sampling)
! More generally e_den=PsiT/PsiG and e_num=E_loc*PsiT/PsiG.  wt represents Psi0*PsiG and wtt represents Psi0*PsiT.
! For sparse PsiT when e_den=0, we replace it with -eps*e_num, so that E_loc is a large negative number. (Warning: Check if this is best)

 use types, only : rk
 implicit none

 integer,intent(in) :: nwalk
 real(rk),intent(in) :: e_num(nwalk), e_den(nwalk)
 real(rk),intent(inout) :: wt(nwalk)
 real(rk) e_num_tmp(nwalk), e_den_tmp(nwalk)

 integer iwalk
 real(rk) wt_fn(nwalk), wt_new(nwalk), wtt(nwalk), wtt_fn(nwalk), wtt_new(nwalk), e_loc(nwalk)
 real(rk) wtt_sum, wtt_fn_sum, wtt_new_sum, e_sum, e_fn_sum, e2_fn_sum, e_new_sum, wtt_av, wtt_fn_av, wtt_new_av, e_av, e_fn_av, e2_fn_av, e_new_av, & ! e2_sum, e2_av, &
& c, alpha, sign_cond, sign_cond_fn, sign_cond_new, min_ratio, max_ratio, p
 real(rk) :: eps=1.e-300_rk, eps2=1.e-6_rk
! false is better choice for abs_wt
 logical abs_wt /.false./
!logical use_new_wts

 if(nwalk.le.60) then
   write(6,'(''  wt imp_dist e_den    e_num'')')
   do iwalk=1,nwalk
     write(6,'(f10.6,i5,9f10.6)') wt(iwalk), imp_distance(iwalk), e_den(iwalk), e_num(iwalk)
   enddo
 endif

 e_num_tmp(1:nwalk)=e_num(1:nwalk)
 e_den_tmp(1:nwalk)=e_den(1:nwalk)

!write(6,'(''e_den='',100es9.1)') e_den(1:nwalk) ! For hubbard2 the wavefunction values in e_den (if imp. sampling is not used) are small.

! The net effect of this is to set the local energy to a large negative number if state is not in PsiT,
! regardless of whether it is in Con(PsiT) or not.
! I don't think it matters whether we use a large negative or positive number since all that does is to change
! the sign of alpha.
! If not in Con(PsiT) set e_num, then set ratio e_num/e_den=-1/eps2
 do iwalk=1,nwalk
   if(abs(e_den_tmp(iwalk)).lt.eps) then
     if(e_num_tmp(iwalk).eq.0._rk) then
       e_num_tmp(iwalk)=-eps*sign(1._rk,wt(iwalk))
     endif
     e_den_tmp(iwalk)=-eps2*e_num_tmp(iwalk) ! alpha changes sign, but it makes no difference to the final result if we use -eps2 or eps2.
   endif
 enddo

 wtt(1:nwalk)=wt(1:nwalk)*e_den_tmp(1:nwalk)
 wtt_sum=sum(wtt(1:nwalk))
 if(wtt_sum.le.0._rk) then ! Switch all the wts because wt is PsiG*Psi0 and wtt is PsiT*Psi0 and we want Psi0 to have +ve overlap with PsiT
   wt(1:nwalk)=-wt(1:nwalk)
   wtt(1:nwalk)=-wtt(1:nwalk)
   wtt_sum=-wtt_sum
 endif

 e_loc(1:nwalk)=e_num_tmp(1:nwalk)/e_den_tmp(1:nwalk)
 e_sum=dot_product(wtt(1:nwalk),e_loc(1:nwalk))
 e_av=e_sum/wtt_sum
!e2_sum=dot_product(abs(wtt(1:nwalk)),e_loc(1:nwalk)**2)
!e2_av=e2_sum/wtt_sum

!write(6,'(''e2_av, e_av**2, e2_sum, e_sum'',9es14.6)') e2_av, e_av**2, e2_sum, e_sum
!if(e2_av.lt.e_av**2) then
!  stop 'e2_av < e_av**2'
!endif

 do iwalk=1,nwalk
   if(e_den(iwalk).ne.0._rk) then                                       ! In PsiT
!    if(wt(iwalk)*e_den(iwalk).gt.0 .or. imp_distance(iwalk).eq.0) then ! Not sign violating or in deterministic space
     if(wtt(iwalk).gt.0) then                                           ! Not sign violating
       wt_fn(iwalk)=wt(iwalk)
!      wt_fn(iwalk)=wt(iwalk) * exp(-(e_loc(iwalk)-e_av)**2/(10*(e2_av-e_av**2)))
!      if((e_loc(iwalk)-e_av)**2.lt.10*(e2_av-e_av**2)+1) then
!        wt_fn(iwalk)=wt(iwalk)
!      else
!        wt_fn(iwalk)=wt(iwalk)/((e_loc(iwalk)-e_av)**2-10*(e2_av-e_av**2))
!      endif
!    else                                                               ! Sign violating and not in deterministic space
     else                                                               ! Sign violating
       if(abs_wt) then
         wt_fn(iwalk)=-wt(iwalk)
       else
         wt_fn(iwalk)=0
       endif
     endif
   elseif(e_num(iwalk).ne.0._rk .or. imp_distance(iwalk).eq.0) then     ! In Con(PsiT)
     if(wt(iwalk)*e_num(iwalk).lt.0) then                               ! Not sign violating or in deterministic space
       wt_fn(iwalk)=wt(iwalk)
     else                                                               ! Sign violating and not in deterministic space
       if(abs_wt) then
         wt_fn(iwalk)=-wt(iwalk)
       else
         wt_fn(iwalk)=0
       endif
     endif
   else                                                                 ! Not in PsiT or Con(PsiT), so don't know about sign, so leave it unchanged
     if(imp_distance(iwalk).eq.0) then                                  ! In deterministic space
       wt_fn(iwalk)=wt(iwalk)
     else                                                               ! Not in deterministic space
       wt_fn(iwalk)=wt(iwalk)*0.9
       !wt_fn(iwalk)=wt(iwalk)
     endif
   endif
 enddo

 wtt_fn(1:nwalk)=wt_fn(1:nwalk)*e_den_tmp(1:nwalk)
!e_loc(1:nwalk)=e_num_tmp(1:nwalk)/e_den_tmp(1:nwalk)

 wtt_fn_sum=sum(wtt_fn(1:nwalk))
!e_sum=dot_product(wtt(1:nwalk),e_loc(1:nwalk))
 e_fn_sum=dot_product(wtt_fn(1:nwalk),e_loc(1:nwalk))
 e2_fn_sum=dot_product(wtt_fn(1:nwalk),e_loc(1:nwalk)**2)

 wtt_av=wtt_sum/nwalk
 wtt_fn_av=wtt_fn_sum/nwalk
!e_av=e_sum/wtt_sum
 e_fn_av=e_fn_sum/wtt_fn_sum
 e2_fn_av=e2_fn_sum/wtt_fn_sum

 !write(6,'(/,''wtt_sum,wtt_fn_sum='',9es12.4)') wtt_sum,wtt_fn_sum
 c=wtt_sum/wtt_fn_sum
 if(e2_fn_av-e_fn_av**2.ne.0._rk) then
   alpha=(e_av-e_fn_av)/(e2_fn_av-e_fn_av**2)
 else
   alpha=0 ! This can be violated if there are few walkers, so in calculation of alpha, denominator is zero but numerator is nonzero.  Just use more starting walkers to cure this.
 endif
 write(6,'(/,''c, alpha='',f10.6,es14.6)') c, alpha

!write(6,'(''   wt        wt_fn     wt_new     e_den_tmp      e_num_tmp      e_loc  wt_new/wt  wt_new/wt_fn'')')
 do iwalk=1,nwalk
   wtt_new(iwalk)=c*wtt_fn(iwalk)*(1+alpha*(e_loc(iwalk)-e_fn_av))
   wt_new(iwalk)=wtt_new(iwalk)/e_den_tmp(iwalk)
   !!write(6,'(5f10.6,f11.6,9f10.6)') wt(iwalk),wt_fn(iwalk),wt_new(iwalk),e_den_tmp(iwalk),e_num_tmp(iwalk),e_loc(iwalk), wt_new(iwalk)/wt(iwalk), wt_new(iwalk)/wt_fn(iwalk)
   !write(6,'(3f10.6,2es12.4,es12.4,9f10.6)') wt(iwalk),wt_fn(iwalk),wt_new(iwalk),e_den_tmp(iwalk),e_num_tmp(iwalk),e_loc(iwalk), wt_new(iwalk)/wt(iwalk), wt_new(iwalk)/wt_fn(iwalk)
 enddo

 wtt_new_sum=sum(wtt_new(1:nwalk))
 e_new_sum=dot_product(wtt_new(1:nwalk),e_loc(1:nwalk))
 wtt_new_av=wtt_new_sum/nwalk
 e_new_av=e_new_sum/wtt_new_sum
!write(6,'(''wtt_sum, wtt_fn_sum, wtt_new_sum, e_sum, e_fn_sum, e_new_sum ='',9f12.6)') wtt_sum, wtt_fn_sum, wtt_new_sum, e_sum, e_fn_sum, e_new_sum
 write(6,'(''wtt_av, wtt_fn_av, wtt_new_av, e_av, e_fn_av, e_new_av, e-av-e_fn_av='',9f12.6)') wtt_av, wtt_fn_av, wtt_new_av, e_av, e_fn_av, e_new_av, e_av-e_fn_av
!write(6,'(''e2_av-e_av**2, e2_fn_av-e_fn_av**2'',9es12.4)') e2_av-e_av**2, e2_fn_av-e_fn_av**2

!Check that total wtt and average energy are unchanged
 if(abs((wtt_sum-wtt_new_sum)/wtt_new_sum).gt.1.e-6_rk) then
   write(6,'(''Warning: wtt_sum != wtt_new_sum'',2f15.6)') wtt_sum, wtt_new_sum
   stop 'abs((wtt_sum-wtt_new_sum)/wtt_new_sum).gt.1.e-6_rk'
 endif
 if(abs((e_sum-e_new_sum)/e_new_sum).gt.1.e-6_rk) then
   write(6,'(''Warning: e_sum != e_new_sum'',2f15.6)') e_sum, e_new_sum
   write(6,'(''e_av,e_fn_av,e2_fn_av,e_fn_av**2),(e_av-e_fn_av),(e2_fn_av-e_fn_av**2)'',9es14.6)') e_av,e_fn_av,e2_fn_av,e_fn_av**2,(e_av-e_fn_av),(e2_fn_av-e_fn_av**2)
   write(6,'(''abs((e_sum-e_new_sum)/e_new_sum).gt.1.e-6_rk) In calculation of alpha, denominator is zero but numerator is nonzero.  Just use more starting walkers to cure this.'')')
   stop 'abs((e_sum-e_new_sum)/e_new_sum).gt.1.e-6_rk) In calculation of alpha, denominator is zero but numerator is nonzero. Just use more starting walkers to cure this.'
 endif

 !!write(6,*) dot_product(wtt_new(1:nwalk),e_loc(1:nwalk)),dot_product(wt_new(1:nwalk),e_num_tmp(1:nwalk)) ! These should be identical

! Note sign_cond_new can be < 1 only if min_ratio is < 0.  So, I use p*wt_new+(1-p)*wt, which should always have sign_cond_new=1.
! However, what is being written here is the sign_cond_new if we used wt_new rather than p*wt_new+(1-p)*wt.
 sign_cond=sum(wtt(1:nwalk))/sum(abs(wtt(1:nwalk)))
 sign_cond_fn=sum(wtt_fn(1:nwalk))/sum(abs(wtt_fn(1:nwalk)))
 sign_cond_new=sum(wtt_new(1:nwalk))/sum(abs(wtt_new(1:nwalk)))
 write(6,'(''sign_cond, sign_cond_new='',9f16.12)') sign_cond, sign_cond_fn, sign_cond_new, sign_cond_fn-sign_cond, sign_cond_new-sign_cond
 if(sign_cond_fn.lt.0.999d0) then
   write(6,'(''wt,    wtt,    wtt_fn,    wtt_new,   eden_tmp,   enum_tmp'')')
   do iwalk=1,nwalk
     write(6,'(9f10.6)') wt(iwalk), wtt(iwalk), wtt_fn(iwalk), wtt_new(iwalk), e_den_tmp(iwalk), e_num_tmp(iwalk)
   enddo
 endif

 p=1
 min_ratio=minval(wtt_new(1:nwalk)/wtt(1:nwalk))
 max_ratio=maxval(wtt_new(1:nwalk)/wtt(1:nwalk))
 if(min_ratio.lt.0) p=1/(1-min_ratio)
 if(max_ratio.gt.2) p=min(p,1/(max_ratio-1))
 write(6,'(''min, max wt_new/wt, wt_new/wt_fn'',9f10.6)') min_ratio, max_ratio, minval(wtt_new(1:nwalk)/wtt_fn(1:nwalk)),maxval(wtt_new(1:nwalk)/wtt_fn(1:nwalk))

 if(min_ratio.lt.0 .or. max_ratio.gt.2) then
   if(nwalk.le.2000) then
!    write(6,'(''Warning: minval(wt_new(1:nwalk)/wt(1:nwalk)), maxval(wt_new(1:nwalk)/wt(1:nwalk)), p='',2es12.4, f6.3)') minval(wt_new(1:nwalk)/wt(1:nwalk)), maxval(wt_new(1:nwalk)/wt(1:nwalk)), p
     write(6,'(''Warning: min_ratio, max_ratio, p='',2es12.4, f6.3)') min_ratio, max_ratio, p
     write(6,'(''   wt        wt_fn     wt_new     e_den_tmp    e_num_tmp      e_loc  wt_new/wt  wt_new/wt_fn'')')
     do iwalk=1,nwalk
       write(6,'(3f10.6,2es12.4,f11.6,9f10.6)') wt(iwalk),wt_fn(iwalk),wt_new(iwalk),e_den_tmp(iwalk),e_num_tmp(iwalk),e_loc(iwalk), wt_new(iwalk)/wt(iwalk), wt_new(iwalk)/wt_fn(iwalk)
     enddo
   endif
 endif

!use_new_wts=.true.
!if(sign_cond_new.lt.sign_cond) then
!  use_new_wts=.false.
!  write(6,'(''Warning: sign_cond, sign_cond_new='',2f16.12)') sign_cond, sign_cond_new
!  if(nwalk.le.2000) then
!    write(6,'(''   wt        wt_fn     wt_new     e_den_tmp    e_num_tmp      e_loc  wt_new/wt  wt_new/wt_fn'')')
!    do iwalk=1,nwalk
!      write(6,'(3f10.6,2es12.4,f11.6,9f10.6)') wt(iwalk),wt_fn(iwalk),wt_new(iwalk),e_den_tmp(iwalk),e_num_tmp(iwalk),e_loc(iwalk), wt_new(iwalk)/wt(iwalk), wt_new(iwalk)/wt_fn(iwalk)
!    enddo
!  endif
!endif

!if(use_new_wts) wt(1:nwalk)=p*wt_new(1:nwalk)+(1-p)*wt(1:nwalk)

 wt(1:nwalk)=p*wt_new(1:nwalk)+(1-p)*wt(1:nwalk)

end subroutine stochastic_reconfiguration2
!-------------------------------------------------------------------------------------

  subroutine add_walker(det_j_up, det_j_dn, walk_wt_child, weight_j, tau,iwalk,excite_level)
    ! Add new walker to walker list
    ! A Holmes, 14 Sep 2015

    use common_psi_t, only : hf_to_psit,dets_up_psi_t,dets_dn_psi_t ! Only for HF->Psi_T
    use chemistry, only : time_sym,excitation_level
    use tools, only : count_excitations

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_j_up,det_j_dn
#else
    integer(ik),intent(in) :: det_j_up,det_j_dn
#endif
    real(rk),intent(in) :: walk_wt_child,weight_j,tau
    integer,intent(in) :: iwalk
    integer,intent(in) :: excite_level
    real(rk) :: weight_j_local


        if (off_diag_histograms) then
          call add_to_hist(abs(weight_j)/tau,nbins,lbounds,bins)
          if (iwalk==1) then
            call add_to_hist(abs(weight_j)/tau,nbins,lbounds_hf,bins_hf)
          endif
        endif
        !Save the max and min ratios
        if(abs(weight_j)/tau.gt.1e-9) then
           if(abs(weight_j).gt.tau*max_wt_ratio) then
             write(6,*) 'walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j=', walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j/tau
            !write(6,'(''walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j='',4i20,f9.2)') walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up, det_j_dn, weight_j/tau
             max_wt_ratio = abs(weight_j)/tau
             if(time_sym) then
               max_wt_ratio_excit_level=min(count_excitations(walk_dets_up(iwalk),det_j_up)+count_excitations(walk_dets_dn(iwalk),det_j_dn),count_excitations(walk_dets_up(iwalk),det_j_dn)+count_excitations(walk_dets_dn(iwalk),det_j_up))
             else
               max_wt_ratio_excit_level=count_excitations(walk_dets_up(iwalk),det_j_up)+count_excitations(walk_dets_dn(iwalk),det_j_dn)
             endif
           endif
           if(abs(weight_j).lt.tau*min_wt_ratio) min_wt_ratio = abs(weight_j)/tau
        endif

        if(importance_sampling.eq.0) then
          weight_j_local=walk_wt_child*weight_j
        else
          weight_j_local=walk_wt_child*weight_j!* psi_g(walk_det_new)/psi_g(walk_det(iwalk))
        endif
        if (off_diag_histograms) then
          call add_to_hist(abs(weight_j_local)/tau,nbins,lbounds_new,bins_new)
          if (excite_level==1) then
            call add_to_hist(abs(weight_j)/tau,nbins,lbounds,bins_single)
          elseif (excite_level==2) then
            call add_to_hist(abs(weight_j)/tau,nbins,lbounds,bins_double)
          endif
        endif

        if(ipr.ge.1) write(6,*) 'walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j=', walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j_local
       !if(ipr.ge.1) write(6,'(''walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j'',4i20,f8.5)') walk_dets_up(iwalk), walk_dets_dn(iwalk), det_j_up,  det_j_dn,  weight_j

        ! Edited by AAH on 15 Jan 2013. Don't allow stochastic spawning onto first state.
        if (hf_to_psit .and. det_j_up==dets_up_psi_t(1) .and. det_j_dn==dets_dn_psi_t(1))  weight_j_local=0._rk
        ! End of edit

        if(weight_j_local.ne.0) then

          nwalk=nwalk+1
          if (ncores==1)  snd_cnt(1)=snd_cnt(1)+1

          if(nwalk.gt.MWALK) then
            write(6,'(''stop nwalk>MWALK, e_est, e_trial, reweight_factor=''9d12.4)') e_est, e_trial, 1._rk/reweight_factor_inv
            call flush(6)
            stop 'nwalk>MWALK'
          endif

          walk_dets_up(nwalk)=det_j_up
          walk_dets_dn(nwalk)=det_j_dn
          walk_wt(nwalk)=weight_j_local

    ! Warning: Impose spin-reversal symmetry assuming +ve relative sign
    !     if(walk_dets_up(nwalk).gt.walk_dets_dn(nwalk)) then
    !       det_j_up=walk_dets_up(nwalk)
    !       walk_dets_up(nwalk)=walk_dets_dn(nwalk)
    !       walk_dets_dn(nwalk)=det_j_up
    !     endif
          if (imp_distance(iwalk)==-2) then
            ! Edited by AAH on 18 Jan 2013. Treat imp_distance=-2 as effectively 1 if c_t_initiator is false
            if (c_t_initiator) then
              imp_distance(nwalk)=1_i1b
            else
              imp_distance(nwalk)=2_i1b
            endif
            ! End of edit
          else
            imp_distance(nwalk)=min(imp_distance(iwalk),huge(1_i1b)-1_i1b)+1_i1b ! Increment imp_distance making sure not to exceed largest i1b integer
          endif
          if(semistochastic) then
            ! Prevents stochastic spawning between important dets
            if (imp_distance(iwalk).eq.0)  imp_distance(nwalk)=-1
          endif
          if(initiator(iwalk).ge.2) then ! If it is spawned by an initiator then 1 otherwise 0
            initiator(nwalk)=1
          else
            initiator(nwalk)=0
          endif
          if (c_t_initiator.and.imp_distance(iwalk)==-2)  initiator(nwalk)=1 ! if c_t_initiator true, then always spawn from C(T)
          if(semistochastic) then
             if (imp_distance(iwalk).eq.0) initiator(nwalk)=1 !always spawn from important space FP
          endif
          e_num_walker(nwalk)=1.e51_rk
          e_den_walker(nwalk)=1.e51_rk
          matrix_elements(nwalk) = 1e51_rk
          if (ncores>1)  call mpi_push_nwalk(nwalk,initiator(nwalk),e_num_walker(nwalk),e_den_walker(nwalk),matrix_elements(nwalk))

        endif

  end subroutine add_walker
!-------------------------------------------------------------------------------------

  subroutine write_walker_info(string,my_nwalk)

  use common_ham, only: hamiltonian_type
  character(len=*),intent(in) :: string
  integer,intent(in) :: my_nwalk
  integer iwalk
  write(6,'(/,a,'' my_nwalk='',i6)') string, my_nwalk
  if(hamiltonian_type.eq.'fictitious' .or. hamiltonian_type.eq.'read' .or. hamiltonian_type.eq.'hubbard') then
      write(6,'(''walk_det='',10000i6)')   (walk_det(iwalk),iwalk=1,my_nwalk)
      write(6,'(''walk_wt= '',10000f6.1)') (walk_wt(iwalk),iwalk=1,my_nwalk)
      write(6,'(''initiator='',10000i6)')   (initiator(iwalk),iwalk=1,my_nwalk)
  else
      !do iwalk=1,my_nwalk
      !    write (6,*) walk_dets_up(iwalk),walk_dets_dn(iwalk),walk_wt(iwalk)
      !enddo
      write(6,'(''walker number  ='',10000i14)') (iwalk,iwalk=1,my_nwalk)
      write(6,'(''walk_dets_up/dn='',10000(2x,2i6))') (walk_dets_up(iwalk),walk_dets_dn(iwalk),iwalk=1,my_nwalk)
      write(6,'(''walk_wt        ='',10000f14.6)') (walk_wt(iwalk),iwalk=1,my_nwalk)
      write(6,'(''e_num_walker   ='',10000es14.6)') (e_num_walker(iwalk),iwalk=1,my_nwalk)
      write(6,'(''e_den_walker   ='',10000es14.6)') (e_den_walker(iwalk),iwalk=1,my_nwalk)
      write(6,'(''initiator      ='',10000i14)') (initiator(iwalk),iwalk=1,my_nwalk)
      write(6,'(''imp_distance   ='',10000i14)') (imp_distance(iwalk),iwalk=1,my_nwalk)
  endif
  call flush(6)

  end subroutine write_walker_info

#endif ! SQMC
end module do_walk
