module hubbard
  use types, only          :    i4b,i8b,rk, ik, ik_vec
#ifdef NUM_ORBITALS_GT_127
  use overload
#endif
  use constants, only      :    pi
  use generic_sort, only   :    sort
  use tools, only          :    n_choose_k, random_int
  use more_tools, only     :    signmod, constrained_dets,fermionic_phase, get_nbr,check_nbr, is_in_list_mod,               &
                              & print_walker,print_walker2,print_real_matrix, print_complex_matrix, print_int_matrix,print_map_on_square, &
                              & print_real_array_with_index,sherman_morrison,print_bin,                                     &
                              & get_nbr_from_adj_list,get_adj_list_44_square,get_adj_list_12_cross,                         &
                              & choose_entry_by_its_weight,choose_entry_by_its_weight_rank2,convert_integer_to_nary,        &
                              & real_symmetric_diagonalize,real_symmetric_diagonalize_ow_ham,real_general_diagonalize_ow_ham, matrix_lanczos,                &
                              & real_sym_gen_eig,binary_search,binary_search_single,overlap_two_slater,gp_two_slater,       &
                              & create_kspace_sym_maps,create_rspace_sym_maps,generate_fourfold_k_configs,print_sym_configs,&
                              & get_rep_only,generate_fourfold_k_configs_efficient,generate_fourfold_k_configs_efficient_given_locations

  use common_run, only     :  tau,run_type,max_connected_dets,connected_dets_up,connected_dets_dn,connected_matrix_elements, connected_matrix_elements_fn, partial_node_eps
  use common_ham, only     :  nelec,nup, ndn,diagonalize_ham, hamiltonian_type, max_energy, energy_exact
  use common_psi_t, only   :  ndet_psi_t, dets_up_psi_t, dets_dn_psi_t, cdet_psi_t, e_trial,trial_wf_iters,norb_trial_wf,n_initiators_trial_wf, n_truncate_trial_wf,ndeg
  use common_selected_ci, only   :  det_sel_iters, norb_det_sel, n_sym_uniq_det_det_sel,lanczos_iters,lanczos_initiators,lanczos_truncate
  use common_imp, only     :  norb_imp, n_imp,imp_up,imp_dn
  implicit none
  save
  public  ::   read_hubbard, system_setup_hubbard,hamiltonian_hubbard, hamiltonian_hubbard_k,hamiltonian_hubbard_dm,       &
             & off_diagonal_move_hubbard, off_diagonal_move_hubbard_k, off_diagonal_move_hubbard_dm,                       &
             & find_connected_dets_hubbard, find_connected_dets_hubbard_k, find_connected_dets_imp_gutz_hubbard,           &
             & energy_pieces_hubbard,energy_pieces_hubbard_k, check_momentum_conservation,generate_sparse_ham_hubbardk,    &
             & generate_sparse_ham_hubbardk_upper_triangular,symmetry_reduce_and_replace_hubbardk                             !subroutines and functions

  ! Common to real,DM and k-space
  integer                   :: z,p                                                                       ! Time reversal and spatial reflection symmetry
  logical                   :: space_sym=.false.                                                         ! Use spatial symmetry of the square or not
  integer                   :: l_x,l_y,nsites                                                            ! Dimensions of the 2D lattice
  real(rk)                  :: t,U                                                                       ! Hopping and Coulomb repulsion
  logical                   :: pbc                                                                       ! Boundary conditions
  real(rk)                  :: g,ln_g                                                                    ! Gutzwiller factor and its logarithm
  real(rk)                  :: aniso=0._rk                                                               ! Anisotropy in t MAY get rid of certain symmetry nodes - breaks trans. inv.
  logical                   :: neel_up_only                                                              ! For the defect wavefunction, it determines if only the up- Neel or both Neel have to be considered
  character*16              :: wf_type,main_wf_type                                                      ! Wavefunction type = RHF, UHF, RHF + UHF, RHF + epsilon, UHF + epsilon
  logical                   :: ham_cross,ham_square                                                      ! Hamiltonian cross and Hamiltonian square turned or not
  real(rk)                  :: energy_exact_hubbard,energy_fn_exact_hubbard                              ! Exact energy from ED or Lanczos
  real(rk)                  :: tau_diag,tau_defect,tau_filling,tau_user                                  ! Various estimates of tau
  real(rk)                  :: w_min
  integer(ik)               :: ndet_hubbard                                                              ! Total number of determinants in the space
  integer                   :: n_det_up,n_det_dn                                                         ! Total number of up and down determinants without symmetry included
  integer                   :: hubbard_ipr                                                               ! Printing option for hubbard
  integer,allocatable       :: c4_map(:,:),reflection_map(:)                                             ! Spatial symmetries of the square lattice
  integer                   :: global_fn_option                                                          ! Fixed node option
  logical                   :: use_projector=.false.                                                     ! Use projector or Hamiltonian for Lanczos/Arnoldi used as a power method

! Specific to DM basis
  integer(ik),allocatable   :: inverse_map_up12(:),inverse_map_dn12(:)                                   ! Inverse maps needed for Lanczos
  real(rk),allocatable      :: lowest_vec_dmbasis(:)                                                     ! Eigenvector in DM basis
  real(rk),allocatable      :: dm_eigenvecs(:,:),dm_eigenvalues(:)                                       ! Linear combinations living on the blocks
  integer,allocatable       :: nup_dm(:),ndn_dm(:)                                                       ! Nup and Ndn corresponding to linear combinations
  integer                   :: ndm_start_end(0:4,0:4,2)                                                  ! Locations at where the eigenvectors are stored
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable   :: all_dets_up(:),all_dets_dn(:)
#else
  integer(ik),allocatable   :: all_dets_up(:),all_dets_dn(:)
#endif
  integer                   :: inverse_dm_map(256)
  character *16             :: dm_cluster                                                                ! Which cluster? square or cross?
  integer                   :: adj_list_44_square(16,4),adj_list_12_cross(12,4)                          ! Adjacency lists - square and cross
  real(rk),allocatable      :: ham_on_2x2(:,:)                                                           ! Hamiltonian on 2x2 patch
  integer(ik)               :: ordered_dets(256)
  real(rk)                  :: eye_256(256,256)                                                          ! 256x256 Identity matrix
  logical                   :: use_eye,use_eig,use_dm                                                    ! Work with which rotation matrix? identity, RDM eigenvectors of 2x2 block? or Eigenbasis of 2x2 block?
  real(rk),allocatable      :: table_block_matel(:,:,:,:,:,:,:)                                          ! Table of 2 block Hamiltonian matrix elements
  logical                   :: diag_44_done                                                              ! Has the 4x4 square lattice Hamiltonian previously been diagonalized?
  real(rk),allocatable      :: stored_guiding_wf(:)                                                      ! The guiding wavefunction can be stored for small systems
  real(rk),allocatable      :: h_psi_t(:)                                                                ! Action of H on psi_T
  logical                   :: guiding_wf_stored=.false.

! Specific to Real space code
  integer,allocatable       :: xycoord(:,:)                                                              ! X and Y coordinates of real space sites
  real(rk),allocatable      :: rhf_up_orbitals(:,:),rhf_dn_orbitals(:,:)                                 ! Up and down Rhf 1 particle orbitals (actually they are the same... so redundant)
  real(rk),allocatable      :: uhf_up_orbitals(:,:),uhf_dn_orbitals(:,:)                                 ! Up and down UHF 1 particle orbitals with a symmetry breaking field m - Refer HJC notes
  real(rk),allocatable      :: up_orbitals(:,:,:),dn_orbitals(:,:,:)                                     ! Up and down 1 particle orbitals for multidet wavefunction
  complex*16,allocatable:: c_up_orbitals(:,:,:),c_dn_orbitals(:,:,:)                                 ! Up and down Complex Rhf 1 particle orbitals (these are just exp(ikr) i.e. plane waves)
 !double complex,allocatable:: c_up_orbitals(:,:,:),c_dn_orbitals(:,:,:)                                 ! Up and down Complex Rhf 1 particle orbitals (these are just exp(ikr) i.e. plane waves)
  real(rk)                  :: cutoff                                                                    ! Wavefunction value cutoff for selecting configuration to start walk
  integer                   :: nmultidet                                                                 ! Number of determinants in multidet wavefunction
  integer                   :: ncsfs                                                                     ! Number of csfs user wants to keep
  real(rk),allocatable      :: coeffs(:)                                                                 ! Coefficients of each multidet wavefunction
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable   :: dets_up_multi(:),dets_dn_multi(:)                                         ! Coefficients of each multidet wavefunction
#else
  integer(ik),allocatable   :: dets_up_multi(:),dets_dn_multi(:)                                         ! Coefficients of each multidet wavefunction
#endif

! Specific to K space code
  integer,allocatable       :: k_vectors(:,:),kmap(:,:,:)                                                ! k vectors and k map denoting scattering map
  real(rk),allocatable      :: k_energies(:)                                                             ! for k-space hubbard model
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec)              :: k_hf_up,k_hf_dn                                                           ! The Hartree Fock determinant in k space
  type(ik_vec),allocatable  :: k_hf_con_up(:),k_hf_con_dn(:)                                             ! dets connected to HF ground state
  type(ik_vec),allocatable  :: rep_k_hf_con_up(:),rep_k_hf_con_dn(:)                                     ! dets connected to HF ground state
  type(ik_vec),allocatable  :: k_hf_deg_up(:),k_hf_deg_dn(:)                                             ! list of configurations with same energy and momentum as HF ground state
#else
  integer(ik)               :: k_hf_up,k_hf_dn                                                           ! The Hartree Fock determinant in k space
  integer(ik),allocatable   :: k_hf_con_up(:),k_hf_con_dn(:)                                             ! dets connected to HF ground state
  integer(ik),allocatable   :: rep_k_hf_con_up(:),rep_k_hf_con_dn(:)                                     ! dets connected to HF ground state
  integer(ik),allocatable   :: k_hf_deg_up(:),k_hf_deg_dn(:)                                             ! list of configurations with same energy and momentum as HF ground state
#endif
  integer,allocatable       :: up_hf_location(:),up_hf_empty_location(:),dn_hf_location(:)               ! List of occupied and unocciped sites in hartree fock state
  integer,dimension(2)      :: k_hf_up_mom,k_hf_dn_mom                                                   ! momentum of hf states up and down stored separately
  integer,dimension(2)      :: ktot                                                                      ! 2D Momentum of the Hartree Fock determinant
  integer                   :: n_k_hf_con,n_k_hf_con_sym                                                 ! number of dets connected to HF ground state
  real(rk)                  :: ubyn                                                                      ! U/number of sites
  real(rk),allocatable      :: hf_matrix_elements(:,:,:)                                                 ! 3 dimensional array to efficiently store connected determinants to hartree Fock and their matrix elements
  real(rk),allocatable      :: k_hf_matrix_elements(:)                                                   ! 1 d array eqvt to the above
  real(rk),allocatable      :: lowest_k_eigenvec(:)                                                      ! Lowest eigenvec in k space
  real(rk)                  :: k_hf_energy                                                               ! Hartree Fock energy
  real(rk),allocatable      :: c_sym_psi_t(:)                                                            ! Trial wavefunction
  integer,allocatable       :: degeneracy_list(:)                                                        ! list of degenerate open shell lattice sites
  integer                   :: nshells                                                                   ! number of individual electron k-shells that we consider for our trial wave function
  integer                   :: n_core_orb=0                                                              ! number of sites in the core. Set to 0 for now because we are not freezing the core in Hubbard.
  integer                   :: init_multiplier, trunc_multiplier

! Constrained Path Monte Carlo
  real(rk),allocatable      :: h_k(:,:),exp_h_k(:,:),cpmc_orbitals(:,:)
  real(rk),allocatable      :: phi_t_up(:,:),phi_t_dn(:,:)
  real(rk),allocatable      :: n_r_up(:),n_r_dn(:)
  real(rk),allocatable      :: k_eigenvalues(:)
  real(rk)                  :: e_t_cpmc,e_v_cpmc,e_k_cpmc
  real(rk)                  :: gamma_,fac_norm
  real(rk)                  :: aux_fld(2,2)
  real(rk)                  :: aux_fld_even,aux_fld_odd
  real(rk),allocatable      :: phi_up(:,:,:),phi_dn(:,:,:)
  real(rk),allocatable      :: walk_wt_cpmc(:),overlaps_cpmc(:)
  integer                   :: nwalk_cpmc

contains
  !===========================================================================
  subroutine read_hubbard
  !---------------------------------------------------------------------------
  ! Description : Read input relevant to the Hubbard.
  ! Created by  : Hitesh J. Changlani/ C.J.U  21 Nov 2010
  ! Edited by   : A.R. 3/12/13 added MPI compatibility
  !---------------------------------------------------------------------------

    use types, only: rk, ik, num_words, bits_per_word
!*** Edited by AR [3/12/13]: changed I/O to handle parallelism
    use mpi_routines, only: master_core, mpi_bsend
!***
    implicit none
    integer i

    if(master_core) then
      read(5,*) l_x,l_y
      write(6,'(''l_x, l_y='', 9i4)') l_x,l_y
    endif
    call mpi_bsend(l_x)
    call mpi_bsend(l_y)
    if((num_words*bits_per_word)<(l_x*l_y)) then
      write(6,'(''Integer word size specified by ik in types.f90 is not large enough for this lattice, (num_words*bits_per_word), (l_x*l_y)='',2i12)') (num_words*bits_per_word), (l_x*l_y)
      stop 'Integer word size specified by ik in types.f90 is not large enough for this lattice'
    endif
    if(master_core) then
      call flush(6)
      read(5,*) pbc, neel_up_only
      write(6,'(''pbc, neel_up_only='',2l3)') pbc, neel_up_only
      call flush(6)
      read(5,*) t,U
      write(6,'(''t, U='', 9f8.2)') t,U
      call flush(6)
      read(5,*) nup,ndn
    endif
    call mpi_bsend(pbc)
    call mpi_bsend(neel_up_only)
    call mpi_bsend(t)
    call mpi_bsend(U)
    call mpi_bsend(nup)
    call mpi_bsend(ndn)
    nelec=nup+ndn
    write(6,'(''nelec, nup, ndn='', 9i4)') nelec,nup,ndn
    call flush(6)
    ! Read spatial symm if hubbardk
    if (hamiltonian_type .eq. 'hubbardk') then
      if(master_core) read(5,*) space_sym         ! True or false
      call mpi_bsend(space_sym)
      if (space_sym) then
         if(master_core) then
           write(6,*) 'Spatial + time symmetries of hubbardk being used'
           read(5,*) z,p            ! Time reversal and reflection (parity)
           write(6,'(''z(time reversal), p(reflection about y-axis)='', 9i4)') z,p
         endif
         call mpi_bsend(z)
         call mpi_bsend(p)
         if (nup .ne. ndn) then
            write (6,*) "Time symmetry for hubbardk will not work when nup and ndn are not equal! Stopping the program...."
            stop 'Time symmetry for hubbardk will not work when nup and ndn are not equal! Stopping the program....'
         endif
         if(abs(z).ne.1) then
           write(6,'(''z(time reversal) must be either 1 or -1'')')
           stop 'z(time reversal) must be either 1 or -1'
         endif
         if(abs(p).ne.1) then
           write(6,'(''p(reflection about y-axis) must be either 1 or -1'')')
           stop 'p(reflection about y-axis) must be either 1 or -1'
         endif
      else
         write(6,*) 'NO symmetries of hubbardk being used'
      endif
    endif

    if(master_core) then
      read(5,*) wf_type
      write(6,'(''trial_wf_type='',a)') wf_type
    endif
    call mpi_bsend(wf_type)
    if (wf_type .eq. 'gutz_aniso') then
      if(master_core) then
        read(5,*) aniso
        write(6,'(''aniso='',f8.3)') aniso
      endif
      call mpi_bsend(aniso)
      wf_type ='gutz'
    endif
    if ((wf_type .eq. 'gutz_multi') .or. (wf_type .eq. 'cgutz_multi')) then
        if(master_core) read(5,*) trial_wf_iters
        call mpi_bsend(trial_wf_iters)
        allocate(norb_trial_wf(max(0,trial_wf_iters)))
        allocate(n_initiators_trial_wf(max(0,trial_wf_iters)))
        allocate(n_truncate_trial_wf(max(0,trial_wf_iters)))
        if(master_core) then
          read(5,*) norb_trial_wf
          read(5,*) n_initiators_trial_wf
          read(5,*) n_truncate_trial_wf
          write (6,*) "trial_wf_iters",trial_wf_iters
          write (6,*) "norb_trial_wf",norb_trial_wf
          write (6,*) "n_initiators_trial_wf",n_initiators_trial_wf
          write (6,*) "n_truncate_trial_wf",n_truncate_trial_wf
          call flush(6)
        endif
        call mpi_bsend(norb_trial_wf)
        call mpi_bsend(n_initiators_trial_wf)
        call mpi_bsend(n_truncate_trial_wf)
    endif
    call flush(6)
    if (wf_type .eq. 'gutz' .or. wf_type .eq. 'gutz_rhf' .or. wf_type .eq. 'gutz_uhf' .or. wf_type .eq. 'gutz_multi' .or. wf_type .eq. 'cgutz_multi'.or. wf_type .eq. 'gutz_rhf_epsilon' .or. wf_type .eq. 'gutz_uhf_epsilon') then
        if(master_core) then
          read(5,*) g
          write(6,'(''g='',f8.3)') g
          call flush(6)
          read(5,*) diagonalize_ham
          write(6,'(''diagonalize_ham (hubbard2/gutz) ='',i4)') diagonalize_ham
        endif
        call mpi_bsend(g)
        call mpi_bsend(diagonalize_ham)
    elseif (wf_type .eq. 'sym') then
      if(master_core) read(5,*) trial_wf_iters
      call mpi_bsend(trial_wf_iters)
      allocate(norb_trial_wf(trial_wf_iters))
      allocate(n_initiators_trial_wf(trial_wf_iters))
      allocate(n_truncate_trial_wf(trial_wf_iters))
      if(master_core) then
        read(5,*) norb_trial_wf
        read(5,*) n_initiators_trial_wf
        read(5,*) n_truncate_trial_wf
        write (6,*) "trial_wf_iters",trial_wf_iters
        write (6,*) "norb_trial_wf",norb_trial_wf
        write (6,*) "n_initiators_trial_wf",n_initiators_trial_wf
        write (6,*) "n_truncate_trial_wf",n_truncate_trial_wf
        read(5,*) diagonalize_ham
        write(6,'(''diagonalize_ham (hubbard) ='',i4)') diagonalize_ham
      endif
      call mpi_bsend(norb_trial_wf)
      call mpi_bsend(n_initiators_trial_wf)
      call mpi_bsend(n_truncate_trial_wf)
      call mpi_bsend(diagonalize_ham)
      w_min=0._rk
      !read(5,*) w_min
      !write(6,'(''w_min ='',f8.3)') w_min
      ! Deterministically selected subspace
      if (run_type.eq.'selected_ci') then
        if(master_core) read(5,*) det_sel_iters
        call mpi_bsend(det_sel_iters)
        allocate(norb_det_sel(det_sel_iters))
        allocate(n_sym_uniq_det_det_sel(det_sel_iters))
        if(master_core) then
          read(5,*) norb_det_sel
          read(5,*) n_sym_uniq_det_det_sel
          read(5,*) init_multiplier,trunc_multiplier
        endif
        call mpi_bsend(norb_det_sel)
        call mpi_bsend(n_sym_uniq_det_det_sel)
        call mpi_bsend(init_multiplier)
        call mpi_bsend(trunc_multiplier)
        norb_det_sel(:) = init_multiplier*norb_det_sel(:)
        n_sym_uniq_det_det_sel(:) = trunc_multiplier*n_sym_uniq_det_det_sel(:)
        if(master_core) then
          write (6,*) "det_sel_iters",det_sel_iters
          write (6,*) "norb_det_sel",norb_det_sel
          write (6,*) "n_sym_uniq_det_det_sel",n_sym_uniq_det_det_sel
        endif
      endif
      if (run_type.eq.'trunc_lanc') then
        if(master_core) then
          read(5,*) lanczos_iters,lanczos_initiators,lanczos_truncate
          write (6,*) "lanczos_iters",lanczos_iters
          write (6,*) "lanczos_initiators",lanczos_initiators
          write (6,*) "lanczos_truncate",lanczos_truncate
        endif
        call mpi_bsend(lanczos_iters)
        call mpi_bsend(lanczos_initiators)
        call mpi_bsend(lanczos_truncate)
      endif
      if(maxval(norb_trial_wf).gt.l_x*l_y) then
        write(6,'(''maxval(norb_trial_wf) must be <= norb = l_x*l_y; norb_trial_wf, norb='',9i5)') maxval(norb_trial_wf), l_x*l_y
        stop 'maxval(norb_imp) must be <= norb = l_x*l_y'
      endif
    endif

    if (wf_type .eq. 'gutz' .or. wf_type .eq. 'gutz_rhf' .or. wf_type .eq. 'gutz_uhf' .or. wf_type .eq. 'gutz_rhf_epsilon' .or. wf_type .eq. 'gutz_uhf_epsilon' .or. wf_type .eq. 'gutz_multi' .or. wf_type .eq. 'cgutz_multi') then
        if (wf_type .eq. 'gutz' .and. hamiltonian_type .eq. 'hubbard2') wf_type='gutz_rhf' ! RHF is the default for real space calculations and is well defined only for closed shells
        main_wf_type='gutz'
    else
        g=0.0
    endif
!   write(6,'(''diagonalize_ham='',i4)') diagonalize_ham
    if(master_core) call flush(6)

    nsites=l_x*l_y
    write (6,*) "Generating k (momentum) vectors"
    call flush(6)
    call generate_k_vectors()                                                                                         !Sets ktot to be momentum of Hartree Fock state and also gives the Hartree Fock state

    allocate (xycoord(2,nsites))
    do i=1,nsites
        xycoord(2,i)=((i-1)/(l_x))+1              ! Y coordinate
        xycoord(1,i)=i-((xycoord(2,i)-1)*l_x)     ! X coordinate
        xycoord(1,i)=xycoord(1,i)-1
        xycoord(2,i)=xycoord(2,i)-1
        write (6,*) "i,X,Y =",i,xycoord(1,i),xycoord(2,i)
    enddo

    do i=1,nsites
      write (6,*) "i,K_X,K_Y =",i,k_vectors(1,i),k_vectors(2,i)
    enddo

    if (hamiltonian_type .eq. 'hubbard2' .and. (wf_type .eq. 'gutz_multi' .or. wf_type .eq. 'cgutz_multi')) then
         ! allocate the connected dets vectors
         max_connected_dets = nup*ndn*min(nsites-nup,nsites-ndn)+1
         write (6,*) "max_connected_dets",max_connected_dets
         call flush(6)
         if (allocated(connected_dets_up))         deallocate(connected_dets_up)
         if (allocated(connected_dets_dn))         deallocate(connected_dets_dn)
         if (allocated(connected_matrix_elements)) deallocate(connected_matrix_elements)

         allocate(connected_dets_up(max_connected_dets))
         allocate(connected_dets_dn(max_connected_dets))
         allocate(connected_matrix_elements(max_connected_dets))
         ! Set maximum energy
         max_energy = 0._rk
         do i=1,nup
           max_energy = max_energy + k_energies(nsites-i+1)
         enddo
         do i=1,ndn
           max_energy = max_energy + k_energies(nsites-i+1)
         enddo
         write (6,*) "Max energy=",max_energy; call flush(6)
     endif

     global_fn_option=0
     if (run_type .eq. 'fixed_node1') then
        global_fn_option=1
     endif
     if (run_type .eq. 'fixed_node2') then
        global_fn_option=2
     endif
     if (run_type .eq. 'fixed_node3') then
        global_fn_option=3
     endif
     if (run_type .eq. 'fixed_node4') then
        global_fn_option=4
     endif

  end subroutine read_hubbard
  !=============================================================================

  subroutine system_setup_hubbard
  !---------------------------------------------------------------------------
  ! Description : Perform various checks
  !               Generate trial wavefunction (Hartree Fock if needed)
  ! Created     : H.J.Changlani
  ! Edited by   : A.R. 3/12/13 added MPI compatibility
  !---------------------------------------------------------------------------
    use common_run, only: tau_multiplier
!*** Edited by AR [3/12/13]: changed I/O to handle parallelism
    use mpi_routines, only: master_core
!***

    implicit none
    !dummy
    real(rk),allocatable         :: ham(:,:)
    !local
    real(rk)                     :: largest_diag_matel,smallest_diag_matel
    integer                      :: i,j,zsav,psav,excite_level
    integer                      :: n,ctr
    integer                      :: n_det
    real(rk), parameter          :: epsilon = 1.e-8_rk
    integer                      :: first(nup+ndn+1),last(nup+ndn+1)
    real(rk)                     :: kdotr

    ham_cross=.false.
    ham_square=.false.
    diag_44_done=.false.

    tau_user=tau
    tau_diag=0._rk
    tau_defect=0._rk

    ! tau_filling needs extra checking when U=0 or rather when t > U
    largest_diag_matel=real(min(nup,ndn))*U/abs(t)
    smallest_diag_matel=0._rk ! As long as we are at or below half filling
    if (nup .gt. (nsites/2) .and. ndn .gt. (nsites/2)) then
        smallest_diag_matel=real(nup+ndn-nsites)*U/abs(t)
    endif
    tau_filling=tau_multiplier/(largest_diag_matel-smallest_diag_matel)
    !tau=tau_filling  ! HJC has some doubt why this statement was added by him

    ! ndet_hubbard = Total number of real space configurations with the right particle number
    n_det_up=int(n_choose_k(nsites,nup),i4b)
    n_det_dn=int(n_choose_k(nsites,ndn),i4b)
    ndet_hubbard=n_choose_k(nsites,nup)*n_choose_k(nsites,ndn)

    ! Real space Hubbard Hartree Fock calculation
    if (hamiltonian_type .eq. 'hubbard2') then
        !write (6,*) "Generating k vectors (for Hubbard2)"
        !call flush(6)
        !call generate_k_vectors()                                                                         !Sets ktot to be momentum of Hartree Fock state and also gives the Hartree Fock state
        if (main_wf_type .eq. 'gutz') then
            write(6,*)
            write(6,*) "TRACE: Doing Hartree Fock calculation"
            call flush(6)
            call do_hartree_fock()
            write(6,*)
            call flush(6)

            ! RHF
            if ((wf_type .eq. 'gutz_rhf') .or.(wf_type .eq. 'gutz_multi') .or. (wf_type .eq. 'cgutz_multi') .or. (wf_type .eq. 'gutz_rhf_epsilon')) then
              if (l_x*l_y .le. 30) then
                  write(6,*) " RHF up orbitals"
                  call print_real_matrix(nsites,nsites,rhf_up_orbitals)
                  write(6,*)
                  write(6,*) " RHF down orbitals"
                  call print_real_matrix(nsites,nsites,rhf_dn_orbitals)
                  write(6,*)
              else
                  write(6,*) "TRACE: I am not printing RHF orbitals to save screen space"
              endif
              call flush(6)
            endif

            ! Multideterminant
            if (wf_type .eq. 'gutz_multi') then
                ! Set orbitals
                allocate (up_orbitals(nmultidet,nsites,nup))
                allocate (dn_orbitals(nmultidet,nsites,ndn))
                do n=1,nmultidet
                    ctr=0
                    do i=1,nsites
                        if (btest(dets_up_multi(n),i-1) .eqv. .true.) then
                            ctr=ctr+1
                            up_orbitals(n,:,ctr)=rhf_up_orbitals(:,i)
                        endif
                    enddo
                    ctr=0
                    do i=1,nsites
                        if (btest(dets_dn_multi(n),i-1) .eqv. .true.) then
                            ctr=ctr+1
                            dn_orbitals(n,:,ctr)=rhf_dn_orbitals(:,i)
                        endif
                    enddo
                    if(master_core) then
                      write (6,*) "n=",n
                      call print_real_matrix(nsites,nup,up_orbitals(n,:,:))
                      write (6,*)
                      call print_real_matrix(nsites,ndn,dn_orbitals(n,:,:))
                      call flush(6)
                    endif
                enddo
            endif

            ! Complex Multideterminant
            if (wf_type .eq. 'cgutz_multi') then
                ! Set orbitals
                allocate (c_up_orbitals(nmultidet,nsites,nup))
                allocate (c_dn_orbitals(nmultidet,nsites,ndn))
                do n=1,nmultidet
                    ctr=0
                    do i=1,nsites
                        if (btest(dets_up_multi(n),i-1) .eqv. .true.) then
                            ctr=ctr+1
                            do j=1,nsites
                                kdotr=(xycoord(1,j)*k_vectors(1,i))+(xycoord(2,j)*k_vectors(2,i))
                                kdotr=kdotr*(pi/real(l_x))               ! For square lattice only - not rectangular
                                c_up_orbitals(n,j,ctr)=(1._rk/sqrt(real(nsites)))*dcmplx(cos(kdotr),sin(kdotr))
                            enddo
                        endif
                    enddo
                    ctr=0
                    do i=1,nsites
                        if (btest(dets_dn_multi(n),i-1) .eqv. .true.) then
                            ctr=ctr+1
                            do j=1,nsites
                                kdotr=(xycoord(1,j)*k_vectors(1,i))+(xycoord(2,j)*k_vectors(2,i))
                                kdotr=kdotr*(pi/real(l_x))
                                c_dn_orbitals(n,j,ctr)=(1._rk/sqrt(real(nsites)))*dcmplx(cos(kdotr),sin(kdotr))
                            enddo
                        endif
                    enddo
                    if(master_core) then
                      write (6,*) "n=",n
                      !write(6,*) c_up_orbitals(n,:,:)
                      !write(6,*) c_dn_orbitals(n,:,:)
                      call print_complex_matrix(nsites,nup,c_up_orbitals(n,:,:))
                      write (6,*)
                      call print_complex_matrix(nsites,ndn,c_dn_orbitals(n,:,:))
                      call flush(6)
                    endif
                enddo
            endif

            ! UHF - introduction of staggered field to break sublattice symmetry
            if(master_core) then
              if ((wf_type .eq. 'gutz_uhf') .or. (wf_type .eq. 'gutz_uhf_epsilon')) then
                if (l_x*l_y .le. 30) then
                    write(6,*) " UHF up orbitals"
                    call print_real_matrix(nsites,nsites,uhf_up_orbitals)
                    write(6,*)
                    write(6,*) " UHF down orbitals"
                    call print_real_matrix(nsites,nsites,uhf_dn_orbitals)
                    write(6,*)
                else
                    write(6,*) "TRACE: I am not printing UHF orbitals to save screen space"
                endif
                call flush(6)
              endif
            endif

            if(master_core) then
              write(6,*)"TRACE: Finished Hartree Fock calculation"
              write(6,*)
              call flush(6)
              write(6,*)"TRACE: Setting Gutzwiller by performing a short VMC calculation"
              write(6,*)
            endif
            call set_gutzwiller()
            if(master_core) then
              write(6,*)
              write(6,*)"TRACE: Finished setting Gutzwiller"
              write(6,*)
              call flush(6)

              write(6,'(''  nup ='',i10)')   nup
              write(6,'(''  ndn ='',i10)')   ndn
              write(6,'(''  g   ='',f10.3)') g
              write(6,'(''  tau   ='',f10.3)') tau
              call flush(6)

              write(6,'(''Now comparing exact wavefunction and gutz-HF'')')
            endif
            if ((diagonalize_ham .eq. 1)) call compare_exact_and_gutz()
            if(master_core) then
              write(6,'(''done comparing exact wavefunction and gutz-HF'')')
              call flush(6)
            endif
        else
            if ( (nup .gt. l_x*l_y) .or. (ndn .gt. l_x*l_y) .or. (nelec .gt. 2*l_x*l_y) ) then
              if(master_core) write(6,*) " Too many electrons for Neel + Defect wavefunction"
            endif

            if ( (nup .eq. (l_x*l_y)/2) .and. (ndn .eq. (l_x*l_y)/2)) then
             call generate_trial_wf_half_filled_hubbard()
            else
             write (6,*)
             write (6,'(''Incorrect number of electrons for half filling wfs; nup, l_x*l_y/2, ndn, l_x*l_y/2'',4i6)') nup, l_x*l_y/2, ndn, l_x*l_y/2
             stop 'Incorrect number of electrons for half filling wfs..... Ending program'
            endif
        endif

     elseif (hamiltonian_type .eq. 'hubbardk') then                                                                      ! Momentum space Hubbard Hamiltonian
        if (g .eq. 0.)  g=1.                                                                                             ! Purely Hartree Fock wavefunction
        if(master_core) write(6,'(''  g   ='',f10.3)') g
        ln_g=log(g)

        !write (6,*) "Generating k vectors"
        !call flush(6)
        !call generate_k_vectors()                                                                                         !Sets ktot to be momentum of Hartree Fock state and also gives the Hartree Fock state
        ! Also sets the degenerate HF state, but not the trial wf any longer.
        ! Also create maps of rotational + reflection + time reversal symmetry

        if (hamiltonian_type .eq. 'hubbardk') then
            ! allocate the connected dets vectors
            max_connected_dets = nup*ndn*min(nsites-nup,nsites-ndn)+1
            write (6,*) "max_connected_dets",max_connected_dets
            call flush(6)
            allocate(connected_dets_up(max_connected_dets))
            allocate(connected_dets_dn(max_connected_dets))
            allocate(connected_matrix_elements(max_connected_dets))
            ! Set maximum energy
            max_energy = 0._rk
            do i=1,nup
              max_energy = max_energy + k_energies(nsites-i+1)
            enddo
            do i=1,ndn
              max_energy = max_energy + k_energies(nsites-i+1)
            enddo
            write (6,*) "Max energy=",max_energy; call flush(6)
        endif

        if (space_sym) then
            zsav=z
            psav=p
            z=1
            p=1
            if(master_core) then
              write(6,*) "z=1,p=1"
              call flush(6)
              call print_sym_configs(c4_map,reflection_map,z,p,15_ik,29_ik)
              call print_sym_configs(c4_map,reflection_map,z,p,23_ik,27_ik)
              call print_sym_configs(c4_map,reflection_map,z,p,31_ik,31_ik)
              call print_sym_configs(c4_map,reflection_map,z,p,157_ik,285_ik)
              call print_sym_configs(c4_map,reflection_map,z,p,1053_ik,155_ik)
            endif
            z=1
            p=-1
            if(master_core) then
              write(6,*) "z=1,p=-1"
              call flush(6)
              call print_sym_configs(c4_map,reflection_map,z,p,15_ik,29_ik)
              call print_sym_configs(c4_map,reflection_map,z,p,23_ik,27_ik)
              call print_sym_configs(c4_map,reflection_map,z,p,31_ik,31_ik)
              call print_sym_configs(c4_map,reflection_map,z,p,157_ik,285_ik)
              call print_sym_configs(c4_map,reflection_map,z,p,1053_ik,155_ik)
            endif
            z=zsav
            p=psav
        endif

        if (diagonalize_ham .eq. 1) then                                                                                  ! If diagonalization requested
            if (ndet_hubbard .le. 30000*(nsites)) then
                if(master_core) then
                  write (6,*) "Performing Exact calculation in k space"
                  call flush(6)
                endif
                call make_hubbard_hamiltonian_k(ktot,ham,all_dets_up,all_dets_dn,lowest_k_eigenvec,energy_exact_hubbard)  ! Matrix diagonalization
                if(master_core) then
                  write (6,*) "Lowest eigenvalue of full Hamiltonian:",energy_exact_hubbard
                  call flush(6)
                endif
            else
                if(master_core) then
                  write (6,*) "Performing Lanczos diagonalization in k space"
                  call flush(6)
                endif
                if (nsites .le. 16 ) then
                    call lanczos_hubbard(lowest_k_eigenvec,energy_exact_hubbard,all_dets_up,all_dets_dn)                      ! Lanczos diagonalization with indexing
                else
                    call lanczos_hubbard_binary_search(lowest_k_eigenvec,energy_exact_hubbard,all_dets_up,all_dets_dn)        ! Lanczos diagonalization with binary search
                endif
                if(master_core) then
                  write (6,*) "Size of lowest_k_eigenvec",size(lowest_k_eigenvec,1)
                  call flush(6)
                endif
                call sort(lowest_k_eigenvec, all_dets_up, all_dets_dn)
                if(master_core) then
                  write (6,*) "Exact Eigenvector in momentum space"
                  call flush(6)
                endif
                if (hubbard_ipr.ge.1) then
                  first(:)=size(all_dets_up,1)
                  last(:)=1
                endif
                do i=1,size(all_dets_up,1)
                  ! Edited by AAH on 17 Jul 2012. Added output of excite level, as well as at what point states at a new excite level are reached.
                  excite_level=popcnt(iand(all_dets_up(i),not(all_dets_up(1))))+popcnt(iand(all_dets_dn(i),not(all_dets_dn(1))))
                  if ((i.le.1000).and.master_core)  write(6,'(2i10,i3,6g14.6)') all_dets_up(i), all_dets_dn(i),excite_level,lowest_k_eigenvec(i)
                  if (hubbard_ipr.ge.1) then
                    if (i<first(excite_level+1))  first(excite_level+1)=i
                    if (i>last(excite_level+1))  last(excite_level+1)=i
                  elseif (i.ge.1000) then
                    exit
                  endif
                enddo
                if ((hubbard_ipr.ge.1).and.master_core) then
                  write (6,*) "excite_level  first_index  last_index"
                  do i=1,nup+ndn+1
                    if (last(i)>1.or.first(i)<size(all_dets_up,1))  write (6,*) i-1,first(i),last(i)
                  enddo
                  call flush(6)
                endif
            endif
            ! Following is for truncating and rediagonalizing to see how well selected CI works, based on hightest 100000 weights from full wavefunction (Commented out for now)
!           previous_value=lowest_k_eigenvec(100000)
!           do i=100001,size(lowest_k_eigenvec,1)
!             if ( abs(abs(previous_value) - abs(lowest_k_eigenvec(i))) > epsilon ) then
!               n_det = i-1
!               exit
!             endif
!           enddo
!           call generate_sparse_ham_hubbardk(all_dets_up(1:n_det),all_dets_dn(1:n_det),H_indices,H_nonzero_elements,H_values)
!           deallocate(lowest_k_eigenvec)
!           allocate(lowest_k_eigenvec(n_det))
!           !write(6,'(''Number of nonzero elements in sparse Hamiltonian='',i10)') nonzero_elements
!           call my_second(2, 'generate_sparse_ham_hubbard')
!           write (6,'(''Performing Lanczos using matrix_lanczos in more_tools, n_det='',i10)') n_det
!           call flush(6)
!           call matrix_lanczos(n_det,lowest_k_eigenvec,energy_exact_hubbard,H_indices,H_nonzero_elements,H_values)
!           call lanczos_hubbard(lowest_k_eigenvec,energy_exact_hubbard,all_dets_up,all_dets_dn)                      ! Lanczos diagonalization
!           write (6,*) "After truncation and rediagonalization, size of lowest_k_eigenvec",n_det
!           call flush(6)
!           call sort(lowest_k_eigenvec(1:n_det), all_dets_up(1:n_det), all_dets_dn(1:n_det))
!           write (6,*) "Exact Eigenvector in momentum space"
!           call flush(6)
            do i=1,n_det
                ! Edited by AAH on 17 Jul 2012. Added output of excite level, as well as at what point states at a new excite level are reached.
                excite_level=popcnt(iand(all_dets_up(i),not(all_dets_up(1))))+popcnt(iand(all_dets_dn(i),not(all_dets_dn(1))))
                if(master_core) write(6,'(2i10,i3,6g14.6)') all_dets_up(i), all_dets_dn(i),excite_level,lowest_k_eigenvec(i)
            enddo

        endif
        energy_exact=energy_exact_hubbard

        if(master_core) then
          write(6,'(''  nup ='',i10)') nup
          write(6,'(''  ndn ='',i10)') ndn
          call flush(6)
        endif
        if (wf_type .eq. 'gutz') then
            if (allocated(cdet_psi_t))    deallocate(cdet_psi_t)
            if (allocated(dets_up_psi_t)) deallocate(dets_up_psi_t)
            if (allocated(dets_dn_psi_t)) deallocate(dets_dn_psi_t)
            allocate(cdet_psi_t(1))
            allocate(dets_up_psi_t(1))
            allocate(dets_dn_psi_t(1))
            cdet_psi_t(1)=1._rk
            dets_up_psi_t(1)=k_hf_up
            dets_dn_psi_t(1)=k_hf_dn
            ndet_psi_t=1
        endif
    elseif (hamiltonian_type .eq. 'hubbarddm') then                                                     ! DM basis calculations
        call get_adj_list_44_square(pbc,adj_list_44_square)                                              ! Obtain adjacency
        call get_adj_list_12_cross(pbc,adj_list_12_cross)
        do i=1,256
            eye_256(i,i)=1._rk
        enddo
        if (real(nup)/real(nsites) .gt. 0.4 .and. real(ndn)/real(nsites) .gt. 0.4) then
            dm_cluster="cross"
            if(master_core) then
              write (6,*)
              write(6,'(/,''==============================================================='')')
              write(6,'(/,''Calculating Density Matrix of a 2x2 patch from a 12 site cross'')')
              write(6,'(/,''==============================================================='')')
            endif
        else
            dm_cluster="square"
            if(master_core) then
              write (6,*)
              write(6,'(/,''==============================================================='')')
              write(6,'(/,''Calculating Density Matrix of a 2x2 patch from a 16 site square'')')
              write(6,'(/,''==============================================================='')')
            endif
        endif
        if(master_core) call flush(6)
        use_eye=.false.
        !use_dm=.false.
        !use_eig=.true.
        use_eig=.false.
        use_dm=.true.

        if(master_core) then
          write (6,*) "Computing reduced density matrix of 2x2 block and its eigenvectors"
          call flush(6)
        endif
        call density_matrix_2by2(dm_cluster)
        if(master_core) then
          write (6,*) "FINISHED Computing reduced density matrix of 2x2 block and its eigenvectors"
          call flush(6)
          write (6,*) "Computing 2 block Hamiltonian tables"
          call flush(6)
        endif
        call make_hamiltonian_tables_two_blocks()
        if(master_core) then
          write (6,*) "FINISHED Computing 2 block Hamiltonian tables"
          call flush(6)
          write (6,*) "Testing Hamiltonian Hubbard DM"
          call flush(6)
        endif
        call test_hamiltonian_hubbard_dm()
        if(master_core) then
          write (6,*) "FINISHED Testing Hamiltonian Hubbard DM"
          call flush(6)
          write (6,*) "Constructing Hamiltonian in most probable basis"
          call flush(6)
        endif
        call ham_most_probable_basis()  ! Sets the trial wavefunction as well
        if(master_core) then
          write (6,*) "FINISHED Constructing Hamiltonian in most probable basis"
          call flush(6)
        endif
        !write (6,*)
        !call project_44wf_dmbasis()
        !write (6,*)
     endif

     ! Rules for tau selection
     ! Choose final tau based on real space considerations
     ! Tau_diag (exact) has first preference even that over user tau
     ! Tau_user has preference over Tau_defect and Tau_filling
     ! Tau_defect wins over Tau_filling

     if(master_core) then
       write (6,*)
       write(6,'(''tau_user    ='',f10.5)') tau_user
       write(6,'(''tau_diag    ='',f10.5)') tau_diag
       write(6,'(''tau_defect  ='',f10.5)') tau_defect
       write(6,'(''tau_filling ='',f10.5)') tau_filling
       call flush(6)
     endif

     if(tau_diag.ne.0._rk) then
            tau=tau_diag
            if(master_core) write(6,'(/,''WARNING: Tau exact might be overwriting user input tau, tau='',f10.5)') tau
     elseif (tau_user .eq.0.0_rk) then
        if (tau_defect .ne. 0.0_rk .and. tau_defect .lt. tau_filling) then
            tau=tau_defect
            if(master_core) write(6,'(/,''Tau from defect wf chosen, tau='',f10.5)') tau
        else
            tau=tau_filling
            if(master_core) write(6,'(/,''Tau from filling consideration chosen, tau='',f10.5)') tau
        endif
     else
           tau=tau_user
           if(master_core) write(6,'(/,''Tau from input chosen, tau='',f10.5)') tau
     endif
     if(master_core) write (6,*)

    end subroutine system_setup_hubbard
  !============================================================================

  !============================================================================
  subroutine do_hartree_fock()
  !============================================================================
  !---------------------------------------------------------------------------
  ! Description : Generates HF solution (for pbc - plane waves)
  ! Method:       Diagonalize One body Kinetic energy operator on lattice for pbc
  !               without worrying about self consistency
  ! Created     : H.J.Changlani
  !---------------------------------------------------------------------------
    implicit none

    integer                                        ::   i,j,k,l
    integer                                        ::   x,y
    integer                                        ::   nbr
    logical                                        ::   allowed
    real(rk), allocatable                          ::   eigenvalues_up(:), eigenvalues_dn(:)
    real(rk)                                       ::   rnd,rannyu
    real(rk)                                       ::   avgn,m
    real(rk),dimension(nsites,nsites)              ::   fock_up,fock_dn,new_fock_up,new_fock_dn
    real(rk)                                       ::   rhf_energy,uhf_energy
    real(rk),allocatable                           ::   uhf_ups(:,:,:),uhf_dns(:,:,:)
    real(rk),allocatable                           ::   h(:,:),evecs(:,:),overlap(:,:),eigs(:)
    integer                                        ::   num_wfs
    real(rk)                                       ::   gp_up(nsites,nsites),gp_dn(nsites,nsites),ov_up,ov_dn

    allocate(eigenvalues_up(nsites))
    allocate(eigenvalues_dn(nsites))
    allocate(rhf_up_orbitals(nsites,nsites))
    allocate(rhf_dn_orbitals(nsites,nsites))
    allocate(uhf_up_orbitals(nsites,nsites))
    allocate(uhf_dn_orbitals(nsites,nsites))

    ! Restricted Hartree Fock with periodic boundary conditions - No SCF required
    fock_up(:,:)=0._rk
    fock_dn(:,:)=0._rk

    do i=1,nsites
        do j=1,4
            call get_nbr(l_x,l_y,pbc,i,j,nbr,allowed)
            if (allowed .and. nbr>i) then
                rnd=rannyu()
                fock_up(i,nbr)=-1._rk+(aniso*rnd)
                fock_up(nbr,i)=-1._rk+(aniso*rnd)
                fock_dn(i,nbr)=-1._rk+(aniso*rnd)
                fock_dn(nbr,i)=-1._rk+(aniso*rnd)
            endif
         enddo
    enddo

    new_fock_up=fock_up
    new_fock_dn=fock_dn

    ! For the RHF case the symmetry breaking factor m is zero
    m=0.0
    avgn=real(nup+ndn)/real(nsites)
    write(6,'(''<n>, m='',2f10.5)') avgn,m
    call flush(6)

    do i=1,nsites
        y=((i-1)/(l_x))+1
        x=i-((y-1)*l_x)
        if (mod(x+y,2) .eq. 0) then
            new_fock_up(i,i)=new_fock_up(i,i)+0.5_rk*U*(avgn+m)
            new_fock_dn(i,i)=new_fock_dn(i,i)+0.5_rk*U*(avgn-m)
        else
            new_fock_up(i,i)=new_fock_up(i,i)+0.5_rk*U*(avgn-m)
            new_fock_dn(i,i)=new_fock_dn(i,i)+0.5_rk*U*(avgn+m)
        endif
    enddo

    call real_symmetric_diagonalize (nsites, new_fock_up,rhf_up_orbitals, eigenvalues_up)
    call real_symmetric_diagonalize (nsites, new_fock_dn,rhf_dn_orbitals, eigenvalues_dn)

    rhf_energy=sum(eigenvalues_up(1:nup))+sum(eigenvalues_dn(1:ndn))
    rhf_energy=rhf_energy-(U/4._rk)*(avgn)*(avgn)*(nsites)
    write(6,'(''RHF energy'',f10.5)') rhf_energy

    write(6,*) "Final RHF Eigenvalues UP electrons"
    write(6,'(100f10.5)') (eigenvalues_up(j),j=1,nsites)
    write(6,*) "Final RHF Eigenvalues DN electrons"
    write(6,'(100f10.5)') (eigenvalues_dn(j),j=1,nsites)
    call flush(6)

    num_wfs=nint(2._rk*0.6/(0.1))+1                                         ! Uniformly spaced "m" wavefunctions

    allocate(uhf_ups(num_wfs,nsites,nup))
    allocate(uhf_dns(num_wfs,nsites,ndn))
    allocate(h(num_wfs,num_wfs))
    allocate(evecs(num_wfs,num_wfs))
    allocate(overlap(num_wfs,num_wfs))
    allocate(eigs(num_wfs))

   avgn=real(nup+ndn)/real(nsites)

    do j=1,num_wfs
        ! m = <n up - ndn> so it can not exceed 1 or be below -1
        m=-0.6+((j-1)*0.1)
        !if (abs(m).le.avgn) then
            new_fock_up=fock_up
            new_fock_dn=fock_dn
            do i=1,nsites
                y=((i-1)/(l_x))+1
                x=i-((y-1)*l_x)
                if (mod(x+y,2) .eq. 0) then
                    new_fock_up(i,i)=new_fock_up(i,i)+0.5_rk*U*(avgn+m)
                    new_fock_dn(i,i)=new_fock_dn(i,i)+0.5_rk*U*(avgn-m)
                else
                    new_fock_up(i,i)=new_fock_up(i,i)+0.5_rk*U*(avgn-m)
                    new_fock_dn(i,i)=new_fock_dn(i,i)+0.5_rk*U*(avgn+m)
                endif
            enddo

            call real_symmetric_diagonalize (nsites, new_fock_up,uhf_up_orbitals, eigenvalues_up)
            call real_symmetric_diagonalize (nsites, new_fock_dn,uhf_dn_orbitals, eigenvalues_dn)

            uhf_ups(j,:,:)=uhf_up_orbitals(:,1:nup)
            uhf_dns(j,:,:)=uhf_dn_orbitals(:,1:ndn)

            uhf_energy=sum(eigenvalues_up(1:nup))+sum(eigenvalues_dn(1:ndn))
            uhf_energy=uhf_energy-(U/4._rk)*(avgn-m)*(avgn+m)*(nsites)
            ! Structure of RHF and UHF orbitals
            !      E1 E2 E3 E4........
            ! R1
            ! R2
            ! R3
            write(6,'(''m,UHF energy'',2f10.5)') m,uhf_energy
            call flush(6)
        !endif
    enddo

    ! Get the multi-determinant wavefunctions linar combination of UHF wavefunctions
    do i=1,num_wfs
       do j=1,num_wfs
            call overlap_two_slater(nsites,nup,uhf_ups(i,:,:),uhf_ups(j,:,:),ov_up)
            call overlap_two_slater(nsites,ndn,uhf_dns(i,:,:),uhf_dns(j,:,:),ov_dn)
            overlap(i,j)=ov_up*ov_dn
            call gp_two_slater(nsites,nup,uhf_ups(i,:,:),uhf_ups(j,:,:),gp_up)
            call gp_two_slater(nsites,ndn,uhf_dns(i,:,:),uhf_dns(j,:,:),gp_dn)
            h(i,j)=0._rk
            do k=1,nsites
                h(i,j)=h(i,j)+(U*gp_up(k,k)*gp_dn(k,k))
                do l=1,nsites
                    h(i,j)=h(i,j)+(fock_up(k,l)*gp_up(k,l))+(fock_dn(k,l)*gp_dn(k,l))
                enddo
            enddo
            h(i,j)=h(i,j)*overlap(i,j)
       enddo
    enddo

    write (6,*) "Hamiltonian in UHF basis"
    call print_real_matrix(num_wfs,num_wfs,h)
    write (6,*) "Overlap Matrix of UHF functions"
    call print_real_matrix(num_wfs,num_wfs,overlap)
    !call real_symmetric_diagonalize(num_wfs,overlap,evecs,eigs)
    !write (6,*) "Eigenvalues of Overlap matrix"
    !call print_real_matrix(num_wfs,1,eigs)
    !call matinv(overlap,num_wfs,ov_up)
    !write (6,*) "Finished inverse of overlap matrix"
    !write (6,*) "Determinant of overlap matrix is ",ov_up
    !h=matmul(overlap,h)
    !write (6,*) "Finished multiplying overlap^-1 times h"
    !stop
    !call real_sym_gen_eig(num_wfs,h,overlap,evecs,eigs)
    !write (6,*) "Transformed Hamiltonian in UHF basis"
    !call print_real_matrix(num_wfs,num_wfs,h)
    !call real_symmetric_diagonalize(num_wfs,h,evecs,eigs)
    !write (6,*) "Eigenvalues of Generalized Eigenproblem"
    !call print_real_matrix(num_wfs,1,eigs)
    !stop

    deallocate(eigenvalues_up)
    deallocate(eigenvalues_dn)

  end subroutine do_hartree_fock
  !=================================================================================================

  !=================================================================================================
  subroutine choose_random_electron(det_up,det_dn,chosen_site,chosen_spin)
  !=================================================================================================
  !-------------------------------------------------------------------------------------------------
  ! Description : Chooses a random electron on the lattice in a given determinant
  !               Return chosen site and chosen spin (0 for spin down, 1 for spin up)
  ! Created by :  Hitesh J Changlani May 27, 2011 for CJU's SQMC code
  !-------------------------------------------------------------------------------------------------
   implicit none
   ! Dummy
#ifdef NUM_ORBITALS_GT_127
   type(ik_vec),intent(in)    ::            det_up,det_dn
#else
   integer(ik),intent(in)    ::            det_up,det_dn
#endif
   integer,intent(out)       ::            chosen_site,chosen_spin
   ! Local
   logical                   ::            found

   found =.false.
   do
        ! pick a 'electronic' site i.e. divide real site into up+dn site at random
        chosen_site=random_int(2*l_x*l_y)
        if (mod(chosen_site,2) .eq. 0) then           ! spin down
            chosen_spin=0
            chosen_site=chosen_site/2
            if (btest(det_dn,chosen_site-1) .eqv. .true.) found=.true.
        else                                          ! spin up
            chosen_spin=1
            chosen_site=(chosen_site+1)/2
            if (btest(det_up,chosen_site-1) .eqv. .true.) found =.true.
      endif
      if (found) exit
  enddo

  end subroutine choose_random_electron

  !=======================================================================================================
  subroutine compute_vmc_energy_gutz_wf(g_in, vmc_energy_mean, vmc_energy_error, max_det_up, max_det_dn)
  !=======================================================================================================
  !-------------------------------------------------------------------------------------------------------
  ! Description : Computes VMC energy of Gutzwiller- HF wavefunction by performing a short run
  !               We have used the simple Metropolis algorithm here
  ! Created by :  Hitesh J Changlani May 1, 2011 for CJU's SQMC code
  !               Sherman-Morrison formula added by Adam Holmes, Jun 16, 2011
  !-------------------------------------------------------------------------------------------------------
  implicit none

  real(rk),intent(in)     :: g_in
  real(rk), intent(out)   :: vmc_energy_mean, vmc_energy_error
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(out) :: max_det_up,max_det_dn
  type(ik_vec)             :: det_up,det_dn,det_up_new,det_dn_new
#else
  integer(ik),intent(out) :: max_det_up,max_det_dn
  integer(ik)             :: det_up,det_dn,det_up_new,det_dn_new
#endif
  integer                 :: i,j,k
  integer                 :: nsamples=100000, equil=20000
  !integer                 :: nsamples=10000000, equil=20000
  real(rk)                :: local_energy,e_den,energy_tot
  real(rk)                :: rannyu
  real(rk)                :: gsav,ginv,wf,wf_ratio,ratio,weight
  real(rk)                :: up_detval,dn_detval,wf_detval
  real(rk)                :: up_inv(nup,nup),dn_inv(ndn,ndn),up_inv_new(nup,nup),dn_inv_new(ndn,ndn)
  integer                 :: hopelec,hopfrom,hopto,ud ! ud = 1(0) means hopping happens in up(dn) det
  real(rk)                :: oldrow_up(nup),newrow_up(nup),oldrow_dn(ndn),newrow_dn(ndn)
  integer                 :: double_occ,delta_dbl_occ
  integer                 :: e_types(4*(nup+ndn)),chosen_e_row(4*(nup+ndn)),chosen_site(4*(nup+ndn)),target_site(4*(nup+ndn)),diff_double_occs(4*(nup+ndn))

  vmc_energy_mean=0._rk
  vmc_energy_error=0._rk
  local_energy=0._rk
  e_den=0._rk
  energy_tot=0._rk
  wf=0._rk

  gsav=g
  g=g_in
  ginv=1.0/g

  do
    call generate_random_config(det_up,det_dn)
    call wf_calc(det_up,det_dn,wf)
    call flush(6)
    !write (6,*) "wf=",wf
    call flush(6)
    if (abs(wf) .gt. cutoff) exit
  enddo

  max_det_up=det_up
  max_det_dn=det_dn
  call get_connected_dets_and_info_hubbard(det_up,det_dn,chosen_e_row,chosen_site,target_site, &
                                         & e_types,up_inv,dn_inv,up_detval,dn_detval,wf_detval,&
                                         & double_occ,diff_double_occs)
  do i=1,nsamples
        call off_diagonal_move_hubbard(det_up,det_dn,det_up_new,det_dn_new,weight,wf_ratio,.true.)
        if ((wf_type .ne. 'gutz_multi') .and. (wf_type .ne. 'cgutz_multi')) then
            if (((det_up .eq. det_up_new) .and. (det_dn .eq. det_dn_new)) .eqv. .false.) then ! if at least one determinant changes...
              ! given binary strings det_up,det_dn,det_up_new,det_dn_new, find which number electron hopped, where from, where to
              ! hopelec is the electron number that has moved, from hopfrom to hopto.

              if (det_up .ne. det_up_new) ud = 1
              if (det_dn .ne. det_dn_new) ud = 0

              if (ud .eq. 1) then
                k = 0
                do j=1,nsites
                  if (btest(det_up,j-1) .eqv. .true.) then
                    k = k+1
                  endif
                  if ((btest(det_up,j-1) .eqv. btest(det_up_new,j-1 )) .eqv. .false.) then ! if electron position different
                    if (btest(det_up,j-1) .eqv. .true.) then
                      hopelec = k
                      hopfrom = j
                    endif
                    if (btest(det_up,j-1) .eqv. .false.) then
                      hopto = j
                    endif
                  endif
                enddo
              else
                k=0
                do j=1,nsites
                  if (btest(det_dn,j-1) .eqv. .true.) then
                    k = k+1
                  endif
                  if ((btest(det_dn,j-1) .eqv. btest(det_dn_new,j-1 )) .eqv. .false.) then
                    if (btest(det_dn,j-1) .eqv. .true.) then
                      hopelec = k
                      hopfrom = j
                    endif
                    if (btest(det_dn,j-1) .eqv. .false.) then
                      hopto = j
                    endif
                  endif
                enddo
              endif

              ! for sm code, newrow = HF orbital evaluated at hopto, rownum = hopelec,
              ! oldrow = HF orbital evaluated at hopfrom.
              ! (for now, we assume that the rows of the matrices are ordered)
              if (ud .eq. 1) then
                do j=1,nup
                  if (wf_type .eq. 'gutz_rhf') then
                      oldrow_up(j) = rhf_up_orbitals(hopfrom,j)
                      newrow_up(j) = rhf_up_orbitals(hopto,j)
                  endif
                  if (wf_type .eq. 'gutz_uhf') then
                      oldrow_up(j) = uhf_up_orbitals(hopfrom,j)
                      newrow_up(j) = uhf_up_orbitals(hopto,j)
                  endif
                enddo
                call sherman_morrison(nup,hopelec,oldrow_up,newrow_up,up_inv,up_inv_new,ratio)
              else
                do j=1,ndn
                  if (wf_type .eq. 'gutz_rhf') then
                    oldrow_dn(j) = rhf_dn_orbitals(hopfrom,j)
                    newrow_dn(j) = rhf_dn_orbitals(hopto,j)
                  endif
                  if (wf_type .eq. 'gutz_uhf') then
                    oldrow_dn(j) = uhf_dn_orbitals(hopfrom,j)
                    newrow_dn(j) = uhf_dn_orbitals(hopto,j)
                  endif
                enddo
                call sherman_morrison(ndn,hopelec,oldrow_dn,newrow_dn,dn_inv,dn_inv_new,ratio)
              endif
            else
                ratio=0._rk
            endif
        endif

        ! count change in number of doubly occupied sites
        delta_dbl_occ=0
        do j=1,nsites
          if (btest(det_up_new,j-1) .and. btest(det_dn_new,j-1)) then
            delta_dbl_occ=delta_dbl_occ+1
          endif
          if (btest(det_up,j-1) .and. btest(det_dn,j-1)) then
            delta_dbl_occ=delta_dbl_occ-1
          endif
        enddo

        if ((wf_type .ne. 'gutz_multi') .and. (wf_type .ne. 'cgutz_multi')) then
            if (delta_dbl_occ>0) ratio=ratio*g
            if (delta_dbl_occ<0) ratio=ratio*ginv
        else
            ratio=abs(wf_ratio)
            !write (6,*) "wf_ratio=",wf_ratio
            !write (6,*) "ratio=",ratio
            !call flush(6)
        endif

        if ((ratio*ratio)>rannyu()) then ! accept and make new measurement
          ! now, we need to shuffle the rows of up_mat,dn_mat (actually just care about the columns of up_inv,dn_inv),
          ! so that the electrons are in order of the site number they occupy
          ! note that permuting rows can only change the sign of the determinant, which doesn't matter for now, since we are squaring the ratio
          if (abs(ratio)>1._rk) then
            max_det_up=det_up_new
            max_det_dn=det_dn_new
          endif

          if ((wf_type .ne. 'gutz_multi') .and. (wf_type .ne. 'cgutz_multi')) then
              ! given binary representations det_up,det_up_new (or same for down),
              ! figure out proper order of columns of up_inv(dn_inv)
              if (ud .eq. 1) then
                up_inv = up_inv_new
                if (nup .gt. 1) call sort_electrons(nup,det_up,det_up_new,up_inv) ! sort columns of inverse if there is more than 1 column
                det_up=det_up_new
              else
                dn_inv = dn_inv_new
                if (ndn .gt. 1) call sort_electrons(ndn,det_dn,det_dn_new,dn_inv)
                det_dn=det_dn_new
              endif
          else
                det_up=det_up_new
                det_dn=det_dn_new
          endif

          if (i>equil) then
            call energy_pieces_hubbard(det_up,det_dn,local_energy,e_den,1)
          endif
       !else
           ! else rejected - do not do anything
        endif

        if (i>equil) energy_tot=energy_tot+local_energy
  enddo
  vmc_energy_mean=energy_tot/real(nsamples-equil)
  g=gsav

  end subroutine compute_vmc_energy_gutz_wf

  !===========================================================================================
  subroutine compare_exact_and_gutz()
  !===========================================================================================
  !-------------------------------------------------------------------------------------------
  ! Description : Compares Exact wavefunction with the Gutzwiller-HF wavefunction
  !               Objective is mainly to compare the quality of nodes and sign structure of wfs
  ! Created by :  Hitesh J Changlani May 1, 2011 for CJU's SQMC code
  !-------------------------------------------------------------------------------------------
   real(rk),allocatable     ::  ham(:,:),propag(:,:)
   real(rk),allocatable     ::  lowest_eigenvec(:),lowest_fn_eigenvec(:)
   integer                  ::  i,j
#ifdef NUM_ORBITALS_GT_127
   type(ik_vec)              ::  det_up,det_dn
#else
   integer(ik)              ::  det_up,det_dn
#endif
   real(rk)                 ::  up_detval,dn_detval
   real(rk)                 ::  wfi,wfj,wf_val,wfsign, lowest_eigenvec_previous
   real(rk)                 ::  top,bot,overlap,cross,norm
   real(rk)                 ::  e_mix_numerator,e_mix_denominator,local_energy
   logical                  ::  exact

   if (ndet_hubbard .le. 200) then
       exact=.true.
       call make_hubbard_matrix_2d(ham,lowest_eigenvec,energy_exact_hubbard,all_dets_up,all_dets_dn,.true.)
       allocate(propag(size(ham,1),size(ham,2)))
       propag=-tau_diag*ham
       do i=1,size(ham,1)
          propag(i,i)=propag(i,i)+1+tau_diag*energy_exact_hubbard
       enddo

       bot=0
       do i=1,size(ham,1)
          bot=bot+dot_product(abs(propag(i,:)),abs(lowest_eigenvec(:)))
       enddo
       top=sum(abs(lowest_eigenvec(:)))

       ! Importance sample projector
       do i=1,size(ham,1)
                call wf_calc(all_dets_up(i),all_dets_dn(i),wfi)
                do j=1,size(ham,1)
                    call wf_calc(all_dets_up(j),all_dets_dn(j),wfj)
                    propag(i,j)=0._rk
                    if (wfj .ne. 0) propag(i,j)=propag(i,j)*wfi/wfj
                enddo
       enddo
       deallocate(propag)
       write(6,'(''tau_diag, top, bot, top/bot, condition_number(new)='',9f10.6)') tau_diag,top,bot,top/bot,(top/bot)**(1/tau_diag)
       write(6,'(''tau_diag, top, bot, top/bot, condition_number(old)='',9f10.6)') tau_diag,top,bot,top/bot,1-(1-top/bot)/tau_diag
       !
    else
        !allocate(stored_guiding_wf(size(all_dets_up)))
        !stored_guiding_wf(:)=0._rk
        exact=.false.
        call flush(6)
        write(6,*) "Lanczos in real space....."
        call flush(6)
        call lanczos_hubbard(lowest_eigenvec,energy_exact_hubbard,all_dets_up,all_dets_dn)
        write(6,*) "Finished Lanczos....."
        call flush(6)

        if (global_fn_option .ne. -1) then
            write(6,*) "Arnoldi with fixed node Hamiltonian in real space....."
            call flush(6)
            call arnoldi_hubbard_binary_search(lowest_fn_eigenvec,energy_fn_exact_hubbard,all_dets_up,all_dets_dn)
            write(6,*) "Finished Arnoldi....."
            call flush(6)
            !stop 'Temporary stop statement'
        endif
   endif
   write(6,'(''  energy_exact ='', f10.6)') energy_exact_hubbard

   e_trial=energy_exact_hubbard

   write(6,*)
   print *,"-------------------------------------------------------------------------------------------------------------"
   print *,"             Comparing Exact wavefunction with Hartree Fock wavefunction                                     "
   write(6,'(''   det_up   det_dn    exact     approx wo J   approx w J            Hpsi        psi        Local energy'')')
   print *,"-------------------------------------------------------------------------------------------------------------"

   call sort(lowest_eigenvec, all_dets_up, all_dets_dn)

   lowest_eigenvec_previous=1
   overlap=0._rk
   norm=0._rk
   cross=0._rk
   wfsign=1._rk

   do i=1,int(ndet_hubbard,i4b)
        det_up=all_dets_up(i)
        det_dn=all_dets_dn(i)
        call wf_calc(det_up,det_dn,wf_val)
        if ( abs(wf_val)>cutoff .and. abs(lowest_eigenvec(i))>cutoff) then
            wfsign=lowest_eigenvec(i)/wf_val
            exit
        endif
   enddo
   write(6,*)

   do i=1,int(ndet_hubbard,i4b)
        det_up=all_dets_up(i)
        det_dn=all_dets_dn(i)
        call energy_pieces_hubbard(det_up, det_dn, e_mix_numerator, e_mix_denominator, 0)
        call wf_calc(det_up,det_dn,wf_val)
        cross=cross+(wf_val*lowest_eigenvec(i))
        norm=norm+(wf_val**(2._rk))
        local_energy=e_mix_numerator/e_mix_denominator
!       if (i==1) wfsign=lowest_eigenvec(i)/wf_val
        if(abs(wf_val*wfsign).lt.1.d-10) then
            write(6,'(2i10,6g14.6,'' zero in HF'')') det_up, det_dn, lowest_eigenvec(i), up_detval*dn_detval*wfsign, wf_val*wfsign,e_mix_numerator,e_mix_denominator,local_energy
        elseif(lowest_eigenvec(i)*wf_val*wfsign.le.-1.d-16) then
            write(6,'(2i10,6g14.6,'' sign flip'')') det_up, det_dn, lowest_eigenvec(i), up_detval*dn_detval*wfsign, wf_val*wfsign,e_mix_numerator,e_mix_denominator,local_energy
        !endif
        else
        !if(abs(abs(lowest_eigenvec(i))-lowest_eigenvec_previous).gt.1.d-12) then
            write(6,'(2i10,6g14.6)') det_up, det_dn, lowest_eigenvec(i), up_detval*dn_detval*wfsign, wf_val*wfsign,e_mix_numerator,e_mix_denominator,local_energy
        endif
        !endif
        lowest_eigenvec_previous=abs(lowest_eigenvec(i))
!       if(lowest_eigenvec(i)*wf_val*wfsign.ge.0) then
!         write(6,'(2i10,3g14.6)') all_dets_up(i), all_dets_dn(i), lowest_eigenvec(i), up_detval*dn_detval*wfsign, wf_val*wfsign
!        else
!         write(6,'(2i10,3f14.6,'' sign flip'')') all_dets_up(i), all_dets_dn(i), lowest_eigenvec(i), up_detval*dn_detval*wfsign, wf_val*wfsign
!       endif
   enddo
   overlap=cross/(norm**(0.5_rk))
   write(6,'(''Overlap ='', g14.6)') overlap
   allocate (lowest_k_eigenvec(size(lowest_eigenvec)))
   lowest_k_eigenvec=lowest_eigenvec
   deallocate(lowest_eigenvec)

  end subroutine compare_exact_and_gutz
  !=========================================================================================

  !=========================================================================================
  subroutine generate_random_config(det_up,det_dn)
  !-----------------------------------------------------------------------------------------
  ! Description : Generates random state, with correct nup,ndown
  ! Created by : Hitesh J Changlani May 1, 2011 for CJU's SQMC code
  !-----------------------------------------------------------------------------------------
  implicit none
  ! Dummy
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(out)  :: det_up,det_dn
#else
  integer(ik),intent(out)  :: det_up,det_dn
#endif
  ! Local
  integer                  :: nup_ctr,ndn_ctr
  integer                  :: rnd_int_1

  det_up=0_ik
  det_dn=0_ik
  nup_ctr=0
  ndn_ctr=0
  do
        if (nup_ctr .eq. nup) exit
        rnd_int_1=random_int(l_x*l_y)
        if ( btest(det_up,rnd_int_1-1) .eqv. .false.) then
            det_up=ibset(det_up,rnd_int_1-1)
            nup_ctr=nup_ctr+1
        endif
   enddo

  do
        if (ndn_ctr .eq. ndn) exit
        rnd_int_1=random_int(l_x*l_y)
        if ( btest(det_dn,rnd_int_1-1) .eqv. .false.) then
            det_dn=ibset(det_dn,rnd_int_1-1)
            ndn_ctr=ndn_ctr+1
        endif
  enddo

  end subroutine generate_random_config
  !============================================================================

  subroutine set_gutzwiller()
  !------------------------------------------------------------------------------------------
  ! Description : Sets the gutzwiller factor based on VMC or user specified (or table lookup)
  ! Created     : H.J.Changlani
  !------------------------------------------------------------------------------------------
    implicit none

    integer               :: i,j,istat
    real(rk)              :: gvals(20),energy_vals(20)
    real(rk)              :: lowest_energy,error
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec)           :: detpsiup,detpsidn,max_updet,max_dndet
#else
    integer(ik)           :: detpsiup,detpsidn,max_updet,max_dndet
#endif

    write(6,'(''TRACE: Presently setting cutoff '')')
    call flush(6)

    ! Calculate cutoff for wavefunction value
    cutoff=0._rk

    if (wf_type .eq. 'gutz_rhf' .or. wf_type .eq. 'gutz_multi' .or. wf_type .eq. 'cgutz_multi') then
        do i=1,max(nup,ndn)
            do j=1,max(nup,ndn)
                    if (abs(rhf_up_orbitals(i,j)) .gt. cutoff) then
                        cutoff=abs(rhf_up_orbitals(i,j))
                    endif
            enddo
        enddo
    endif

    if (wf_type .eq. 'gutz_uhf') then
        do i=1,max(nup,ndn)
            do j=1,max(nup,ndn)
                    if (abs(uhf_up_orbitals(i,j)) .gt. cutoff) then
                        cutoff=abs(uhf_up_orbitals(i,j))
                    endif
            enddo
        enddo
    endif

    cutoff=(cutoff/2._rk)**real(nup+ndn)
    write(6,'('' wf_value_cutoff ='',g10.3)') cutoff
    call flush(6)

    lowest_energy=0._rk

    ! Optimize Gutzwiller factor if set to 0 by user
    if (g .eq. 0.0_rk) then
        write (6,*) "The user has set a Gutzwiller factor of 0, so I am doing a line search for the best factor"
        do j=1,20
          gvals(j)=real(j)*0.05
          call compute_vmc_energy_gutz_wf(gvals(j),energy_vals(j),error,max_updet,max_dndet)
          write(6,'(''g, VMC energy='',2f10.5)') gvals(j), energy_vals(j)
          call flush(6)
          if ((j==1) .or. (energy_vals(j)<lowest_energy)) then
            lowest_energy=energy_vals(j)
            g=gvals(j)
            detpsiup=max_updet
            detpsidn=max_dndet
          endif
        enddo
        write(6,'(''Optimal g, VMC energy='',2f10.5)') g, lowest_energy
        call flush(6)
    else  ! If user has set a Gutzwiller factor simply use it
        write (6,*) "The user set a NON ZERO Gutzwiller factor so I am using it for the Monte Carlo simulation"
        call compute_vmc_energy_gutz_wf(g,lowest_energy,error,detpsiup,detpsidn)
        write(6,'(''g (read in), VMC energy='',2f10.5)') g, lowest_energy
    endif
    call flush(6)

    ! Note we can not store this wavefunction, however we initially populate using this wavefunction
    ndet_psi_t=1
    allocate(cdet_psi_t(ndet_psi_t),stat=istat)
    allocate(dets_up_psi_t(ndet_psi_t),stat=istat)
    allocate(dets_dn_psi_t(ndet_psi_t),stat=istat)

    do i=1,ndet_psi_t
        if (i==1) then
            dets_up_psi_t(i)=detpsiup
            dets_dn_psi_t(i)=detpsidn
            call wf_calc(dets_up_psi_t(i),dets_dn_psi_t(i),cdet_psi_t(i))
        else
            write(6,'(''i='',i5)') i
            write(6,'(''calling generate_random_config'')')
            call generate_random_config(dets_up_psi_t(i),dets_dn_psi_t(i))
            write(6,'(''Finished calling generate_random_config'')')
            call wf_calc(dets_up_psi_t(i),dets_dn_psi_t(i),cdet_psi_t(i))
        endif
    enddo

    write(6,'(''run_type (hubbard) ='',a)') run_type
    call flush(6)
    if ((run_type .eq. 'fixed_node1') .or. (run_type .eq. 'fixed_node2')) then
        do i=1,ndet_psi_t
            cdet_psi_t(i)=abs(cdet_psi_t(i))
        enddo
    endif

  end subroutine set_gutzwiller
  !=============================================================================================

  !=============================================================================================
  subroutine hamiltonian_hubbard(det_i_up, det_i_dn, det_j_up,det_j_dn, matrix_element)
  !=============================================================================================
  !---------------------------------------------------------------------------------------------
  ! Description   : Bit packed version to generate Hamiltonian matrix element for 2D lattice
  ! Note          : NOT checking for neighbors, assumes inherent connectedness
  !                 Check option
  ! Author        : H.J. Changlani, 21 Nov 2010
  !---------------------------------------------------------------------------------------------

  implicit none
  ! dummy variables
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in)      :: det_i_up,det_i_dn,det_j_up,det_j_dn
  type(ik_vec)                 :: tmp_config,tmp_config_1,tmp_config_2
#else
  integer(ik),intent(in)      :: det_i_up,det_i_dn,det_j_up,det_j_dn
  integer(ik)                 :: tmp_config,tmp_config_1,tmp_config_2
#endif

  !local variables
  real(rk),   intent(out)     :: matrix_element
  integer                     :: i,pos_1,pos_2,double_occupied,num_ones
  logical                     :: flag

  double_occupied=0
  pos_1=0
  pos_2=0
  tmp_config=0_ik
  tmp_config_1=0_ik
  tmp_config_2=0_ik
  flag=.false.

  ! Diagonal term
  if (det_i_up==det_j_up .and. det_i_dn==det_j_dn) then
    tmp_config_1=iand(det_i_up,det_i_dn)
    do i=1,l_x*l_y
      if (btest(tmp_config_1,i-1)) then
         double_occupied=double_occupied+1
      endif
    enddo
    matrix_element=U*double_occupied
    return
  endif

  ! Terms separated by 1 down electron hop
  if (det_i_up==det_j_up) then
   tmp_config_1=iand(det_i_dn,not(det_j_dn))
   tmp_config_2=iand(not(det_i_dn),det_j_dn)
   tmp_config=det_i_dn
   flag=.true.
  endif

  ! Terms separated by 1 up electron hop
  if (det_i_dn==det_j_dn) then
    tmp_config_1=iand(det_i_up,not(det_j_up))
    tmp_config_2=iand(not(det_i_up),det_j_up)
    tmp_config=det_i_up
    flag=.true.
  endif

  if (flag .eqv. .false.) then
    matrix_element=0._rk
    return
  endif

  num_ones=0
  ! Check 1st tmp config for position 1
  do i=1,l_x*l_y
    if (btest(tmp_config_1,i-1)) then
        num_ones=num_ones+1
        if (num_ones .eq. 1) then
            pos_1=i
            !print *,"pos_1=",pos_1
        else
            matrix_element=0._rk
            return
        endif
    endif
  enddo

  if (num_ones .eq. 0 ) then
     matrix_element=0._rk
     return
  endif

  num_ones=0
  ! Check 2nd tmp config for position 2
  do i=1,l_x*l_y
    if (btest(tmp_config_2,i-1)) then
        num_ones=num_ones+1
        if (num_ones .eq. 1) then
            pos_2=i
            !print *,"pos_2=",pos_2
        else
            matrix_element=0._rk
            return
        endif
    endif
  enddo

  if (num_ones .eq. 0 ) then
     matrix_element=0._rk
     return
  endif

  ! Matrix element = -t * fermion phase
  matrix_element=-t*fermionic_phase(tmp_config,min(pos_1,pos_2),max(pos_1,pos_2))

  return
  end subroutine hamiltonian_hubbard

! ============================================================================================
  subroutine make_hamiltonian_tables_two_blocks()
  ! ==========================================================================================
  ! ------------------------------------------------------------------------------------------
  ! Description   : Bit packed version to generate Hamiltonian matrix elements for 2D lattice
  ! Author        : H.J. Changlani, Feb 15 2012
  ! ------------------------------------------------------------------------------------------

  implicit none

  ! Local variables
  integer(ik)              :: i,j,k,l
  integer                  :: nstates_bra_block1,nstates_bra_block2
  integer                  :: nup1,nup2,ndn1,ndn2,nup1_new,nup2_new,ndn1_new,ndn2_new
  integer                  :: nbr_type
  integer                  :: shift1,shift2
  integer                  :: updn,lost,elost

  allocate (table_block_matel(2,2,4,256,256,36,36))
  table_block_matel(:,:,:,:,:,:,:)=0._rk

  do updn=1,2
      do lost=1,2
          if (lost .eq. 1) elost=-1
          if (lost .eq. 2) elost=1
          do nbr_type=1,4
              do i=1,256
                nup1=nup_dm(i)
                ndn1=ndn_dm(i)
                do j=1,256
                   nup2=nup_dm(j)
                   ndn2=ndn_dm(j)
                   if (updn .eq. 1) then
                    nup1_new=nup1+elost
                    nup2_new=nup2-elost
                    ndn1_new=ndn1
                    ndn2_new=ndn2
                   else
                    nup1_new=nup1
                    nup2_new=nup2
                    ndn1_new=ndn1+elost
                    ndn2_new=ndn2-elost
                   endif
                   if ((nup1_new .ge. 0) .and. (nup1_new .le.4) .and. (ndn1_new .ge. 0) .and. (ndn1_new .le. 4) .and. (nup2_new .ge. 0) .and. (nup2_new .le.4) .and. (ndn2_new .ge. 0) .and. (ndn2_new .le. 4)) then
                       shift1=ndm_start_end(nup1_new,ndn1_new,1)-1
                       shift2=ndm_start_end(nup2_new,ndn2_new,1)-1
                       nstates_bra_block1=ndm_start_end(nup1_new,ndn1_new,2)-ndm_start_end(nup1_new,ndn1_new,1)+1
                       nstates_bra_block2=ndm_start_end(nup2_new,ndn2_new,2)-ndm_start_end(nup2_new,ndn2_new,1)+1
                       do k=1,nstates_bra_block1
                            do l=1,nstates_bra_block2
                                call compute_hamiltonian_matrix_element_dm(i,j,k+shift1,l+shift2,nbr_type,table_block_matel(updn,lost,nbr_type,i,j,k,l))
                            enddo
                       enddo
                   endif
                enddo
              enddo
          enddo
     enddo
 enddo

 end subroutine make_hamiltonian_tables_two_blocks

! ========================================================================================================
  subroutine compute_hamiltonian_matrix_element_dm(ket_b1,ket_b2,bra_b1,bra_b2,nbr_type,matrix_element)
  ! ======================================================================================================
  ! ------------------------------------------------------------------------------------------------------
  ! Author        : H.J. Changlani, Feb 15 2012
  ! ------------------------------------------------------------------------------------------------------

  implicit none

  ! Dummy
  integer(ik),intent(in)              :: ket_b1,ket_b2,bra_b1,bra_b2
  integer,intent(in)                  :: nbr_type
  real(rk),intent(out)                :: matrix_element

  ! Local
  integer                             :: site_1,site_2,site_11,site_12,site_21,site_22
  real(rk)                            :: w1,w2,w3,w4,sign_w
  logical                             :: spr
  integer                             :: i,j,k,term
  integer                             :: n_tot,n_tot_j
  integer                             :: up_diff_1,up_diff_2,dn_diff_1,dn_diff_2
  integer                             :: updn_shift
  logical                             :: up,electron_lost
  integer(ik)                         :: ind1,ind2,ind3,ind4

!  spr=.true.
  spr=.false.

  up_diff_1=nup_dm(bra_b1)-nup_dm(ket_b1)
  up_diff_2=nup_dm(bra_b2)-nup_dm(ket_b2)
  dn_diff_1=ndn_dm(bra_b1)-ndn_dm(ket_b1)
  dn_diff_2=ndn_dm(bra_b2)-ndn_dm(ket_b2)

 if (spr) then
     write (6,*) "up_diff_1",up_diff_1
     write (6,*) "up_diff_2",up_diff_2
     write (6,*) "dn_diff_1",dn_diff_1
     write (6,*) "dn_diff_2",dn_diff_2
     call flush(6)
 endif

 if ((dn_diff_1 .eq. 0) .and. (dn_diff_2 .eq. 0)) then
     if ((up_diff_1 .eq. -1) .and. (up_diff_2 .eq. 1)) then
        up=.true.
        electron_lost=.true.
     elseif ((up_diff_1 .eq. 1) .and. (up_diff_2 .eq. -1)) then
        up=.true.
        electron_lost=.false.
     else
        matrix_element=0._rk
        return
     endif
 elseif ((up_diff_1 .eq. 0) .and. (up_diff_2 .eq. 0)) then
     if ((dn_diff_1 .eq. -1) .and. (dn_diff_2 .eq. 1)) then
        up=.false.
        electron_lost=.true.
     elseif ((dn_diff_1 .eq. 1) .and. (dn_diff_2 .eq. -1)) then
        up=.false.
        electron_lost=.false.
     else
        matrix_element=0._rk
        return
     endif
 else
    matrix_element=0._rk
    return
 endif

 if (spr) then
     write (6,*) "up=",up
     call flush(6)
 endif

 if (nbr_type .eq. 1) then     ! Left orientation of block2 wrt block 1
    site_11=1
    site_21=2
    site_12=3
    site_22=4
 endif

 if (nbr_type .eq. 2) then     ! Right orientation
    site_11=2
    site_21=1
    site_12=4
    site_22=3
 endif

 if (nbr_type .eq. 3) then     ! Up orientation
    site_11=1
    site_21=3
    site_12=2
    site_22=4
 endif

 if (nbr_type .eq. 4) then     ! Down orientation
    site_11=3
    site_21=1
    site_12=4
    site_22=2
 endif

 ! Evaluate matrix_element CAREFULLY !! Can be made more efficient!!!!!!
 updn_shift=0
 if (up) updn_shift=4

 do term=1,2
     if (term .eq. 1) then
        site_1=site_11
        site_2=site_21
     else
        site_1=site_12
        site_2=site_22
     endif

     if (spr) then
         write (6,*) "site_1 =",site_1
         write (6,*) "site_2 =",site_2
         write (6,*)
         call flush(6)
         write(6,*) "electron_lost from block (in ket) =",electron_lost
         call flush(6)
     endif

     if (electron_lost) then
         do i=ndm_start_end(nup_dm(ket_b1),ndn_dm(ket_b1),1),ndm_start_end(nup_dm(ket_b1),ndn_dm(ket_b1),2)
            n_tot=0
            ind1=ordered_dets(i)
            if (btest(ind1,updn_shift+site_1-1)) then                                              ! Electron present
                w1=dm_eigenvecs(i,ket_b1)
                !if (spr) then
                !    write(6,*) "ind1=",ind1
                !    write(6,*) "w1=",w1
                !    write(6,*) "i=",i
                !    write(6,*) "inverse_dm_map(ind1+1)... must be equal to i",inverse_dm_map(ind1+1)
                !    call flush(6)
                !endif
                if (abs(w1) .gt. 1.0e-12) then
                    ind3=ibclr(ind1,updn_shift+site_1-1)
                    w3=dm_eigenvecs(inverse_dm_map(ind3+1),bra_b1)
                    !if (spr) then
                    !    write(6,*) "ind3=",ind3
                    !    write(6,*) "w3=",w3
                    !    write(6,*) "inverse_dm_map(ind3+1)",inverse_dm_map(ind3+1)
                    !    call flush(6)
                    !endif
                    if (abs(w3) .gt. 1.0e-12) then
                        do k=site_1+1,4
                            if (btest(ind1,updn_shift+k-1)) n_tot=n_tot+1                                ! Electrons after site in block 1 for getting phase
                        enddo
                        do j=ndm_start_end(nup_dm(ket_b2),ndn_dm(ket_b2),1),ndm_start_end(nup_dm(ket_b2),ndn_dm(ket_b2),2)
                            n_tot_j=n_tot
                            ind2=ordered_dets(j)
                            if (btest(ind2,updn_shift+site_2-1) .eqv. .false.) then                     ! Electron missing
                                w2=dm_eigenvecs(j,ket_b2)
                                !if (spr) write(6,*) "ind2=",ind2
                                !if (spr) write(6,*) "w2=",w2
                                if (abs(w2).gt. 1.0e-12) then
                                    do k=1,site_2-1
                                        if (btest(ind2,updn_shift+k-1)) n_tot_j=n_tot_j+1               ! Electrons before site in block 2 for getting phase
                                    enddo
                                    if (mod(n_tot_j,2).ne.0) then
                                        sign_w=-1._rk
                                    else
                                        sign_w=1._rk
                                    endif
                                    ind4=ibset(ind2,updn_shift+site_2-1)
                                    w4=dm_eigenvecs(inverse_dm_map(ind4+1),bra_b2)
                                    !if (spr) write(6,*) "ind4=",ind4
                                    !if (spr) write(6,*) "w4=",w4
                                    matrix_element=matrix_element+(w1*w2*w3*w4*sign_w*(-t))
                                    !if (spr) write (6,*) "w1*w2*w3*w4*sign_w*(-t)=",(w1*w2*w3*w4*sign_w*(-t))
                                    !if (spr) write (6,*) "Present value of matrix element",matrix_element
                                    !if (spr) write (6,*)
                                    !if (spr) call flush(6)
                                endif
                            endif
                        enddo
                    endif
                endif
             endif
         enddo
     else
         do i=ndm_start_end(nup_dm(ket_b1),ndn_dm(ket_b1),1),ndm_start_end(nup_dm(ket_b1),ndn_dm(ket_b1),2)
            n_tot=0
            ind1=ordered_dets(i)
            if (btest(ind1,updn_shift+site_1-1) .eqv. .false.) then                                 ! Electron absent
                w1=dm_eigenvecs(i,ket_b1)
                !if (spr) then
                !    write(6,*) "ind1=",ind1
                !    write(6,*) "w1=",w1
                !    write(6,*) "i=",i
                !    write(6,*) "inverse_dm_map(ind1+1)... must be equal to i",inverse_dm_map(ind1+1)
                !endif
                if (abs(w1) .gt. 1.0e-12) then
                    ind3=ibset(ind1,updn_shift+site_1-1)
                    w3=dm_eigenvecs(inverse_dm_map(ind3+1),bra_b1)
                    !if (spr) then
                    !    write(6,*) "ind3=",ind3
                    !    write(6,*) "w3=",w3
                    !    write(6,*) "inverse_dm_map(ind3+1)",inverse_dm_map(ind3+1)
                    !endif
                    if (abs(w3) .gt. 1.0e-12) then
                        do k=site_1+1,4
                            if (btest(ind1,updn_shift+k-1)) n_tot=n_tot+1                                 ! Electrons after site in block 1 for getting phase
                        enddo
                        do j=ndm_start_end(nup_dm(ket_b2),ndn_dm(ket_b2),1),ndm_start_end(nup_dm(ket_b2),ndn_dm(ket_b2),2)
                            n_tot_j=n_tot
                            ind2=ordered_dets(j)
                            if (btest(ind2,updn_shift+site_2-1)) then                                   ! Electron present
                                w2=dm_eigenvecs(j,ket_b2)
                                !if (spr) write(6,*) "ind2=",ind2
                                !if (spr) write(6,*) "w2=",w2
                                if (abs(w2).gt. 1.0e-12) then
                                    do k=1,site_2-1
                                        if (btest(ind2,updn_shift+k-1)) n_tot_j=n_tot_j+1                ! Electrons before site in block 2 for getting phase
                                    enddo
                                    if (mod(n_tot_j,2).ne.0) then
                                        sign_w=-1._rk
                                    else
                                        sign_w=1._rk
                                    endif
                                    ind4=ibclr(ind2,updn_shift+site_2-1)
                                    w4=dm_eigenvecs(inverse_dm_map(ind4+1),bra_b2)
                                    !if (spr) write(6,*) "ind4=",ind4
                                    !if (spr) write(6,*) "w4=",w4
                                    !if (spr) call flush(6)
                                    matrix_element=matrix_element+(w1*w2*w3*w4*sign_w*(-t))
                                    !if (spr) write (6,*) "w1*w2*w3*w4*sign_w*(-t)=",(w1*w2*w3*w4*sign_w*(-t))
                                    !if (spr) write (6,*) "Present value of matrix element",matrix_element
                                    !if (spr) write (6,*)
                                    !if (spr) call flush(6)
                                endif
                            endif
                        enddo
                    endif
                endif
             endif
         enddo
     endif
 enddo

 end subroutine compute_hamiltonian_matrix_element_dm

! ======================================================================================================================
 subroutine get_matrix_element_from_table(ket_b1,ket_b2,bra_b1,bra_b2, nbr_type, electron_lost, up, matrix_element)
! =====================================================================================================================

 implicit none
 real(rk),intent(out)     :: matrix_element
 integer(ik),intent(in)   :: ket_b1,ket_b2,bra_b1,bra_b2
 integer,intent(in)       :: nbr_type
 logical,intent(in)       :: electron_lost,up

 ! Dummy
 integer                  :: shift1,shift2
 integer                  :: upint,elost

 shift1=ndm_start_end(nup_dm(bra_b1),ndn_dm(bra_b1),1)-1
 shift2=ndm_start_end(nup_dm(bra_b2),ndn_dm(bra_b2),1)-1

 if (up) then
     upint=1
 else
    upint=2
 endif

 if (electron_lost) then
    elost=1
 else
    elost=2
 endif

 matrix_element=table_block_matel(upint,elost,nbr_type,ket_b1,ket_b2,bra_b1-shift1,bra_b2-shift2)

 end subroutine get_matrix_element_from_table

! =======================================================================================================
  subroutine hamiltonian_hubbard_dm(det_ket_up, det_ket_dn, det_bra_up,det_bra_dn, matrix_element)
  ! =====================================================================================================
  ! ------------------------------------------------------------------------------------------
  ! Description   : Bit packed version to generate Hamiltonian matrix element for 2D lattice
  !                 One needs to interpret the dets in the DM basis
  ! Author        : H.J. Changlani, Nov/Dec 2011
  ! ------------------------------------------------------------------------------------------

  implicit none
  ! dummy variables
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in)      :: det_ket_up,det_ket_dn,det_bra_up,det_bra_dn
#else
  integer(ik),intent(in)      :: det_ket_up,det_ket_dn,det_bra_up,det_bra_dn
#endif
  real(rk),   intent(out)     :: matrix_element

  ! Local variables
  real(rk)                    :: matrix_element_tmp
  integer                     :: num_blocks
  integer                     :: difference,which_diff_1,which_diff_2
  integer(ik)                 :: ind_ket(nsites/4),ind_bra(nsites/4)
  integer                     :: ket_nup_block(nsites/4),ket_ndn_block(nsites/4)
  integer                     :: bra_nup_block(nsites/4),bra_ndn_block(nsites/4)
  integer                     :: i,j,n,num_types,nbr_types(4)
  integer                     :: n_tot
  integer                     :: up_diff_1,up_diff_2,dn_diff_1,dn_diff_2
  logical                     :: up,electron_lost
  logical                     :: spr

!  spr=.true.
  spr=.false.
  num_blocks=nsites/4
  matrix_element=0._rk

  do i=1,num_blocks
    ind_ket(i)=0_ik
    do j=1,4
        if (btest(det_ket_up,4*(i-1)+j-1)) then
            ind_ket(i)=ibset(ind_ket(i),4+j-1)
        endif
        if (btest(det_ket_dn,4*(i-1)+j-1)) then
            ind_ket(i)=ibset(ind_ket(i),j-1)
        endif
    enddo
    ket_nup_block(i)=nup_dm(ind_ket(i))
    ket_ndn_block(i)=ndn_dm(ind_ket(i))
  enddo

  difference=0

  do i=1,num_blocks
    ind_bra(i)=0_ik
    do j=1,4
        if (btest(det_bra_up,4*(i-1)+j-1)) then
            ind_bra(i)=ibset(ind_bra(i),4+j-1)
        endif
        if (btest(det_bra_dn,4*(i-1)+j-1)) then
            ind_bra(i)=ibset(ind_bra(i),j-1)
        endif
    enddo
    bra_nup_block(i)=nup_dm(ind_bra(i))
    bra_ndn_block(i)=ndn_dm(ind_bra(i))
    if (ind_bra(i) .ne. ind_ket(i)) then
        difference=difference+1
        if (difference .eq. 1) then
            which_diff_1=i
        else
            which_diff_2=i
        endif

        if (difference >2) then
            matrix_element=0._rk
            return
        endif
    endif
  enddo

  if (spr) then
      write (6,*) "nup on blocks of bra", bra_nup_block
      write (6,*) "ndn on blocks of bra", bra_ndn_block
      write (6,*) "ind on blocks of bra", ind_bra
      write (6,*)
      call flush(6)
  endif

  ! If difference is exactly zero
  if (difference .eq. 0) then
    do i=1,num_blocks
        matrix_element=matrix_element+ham_on_2x2(ind_ket(i),ind_ket(i))     ! Add up diagonal terms on each block
    enddo
    if (spr) write (6,*) "Diagonal term being evaluated"
    if (spr) write (6,*) "matrix_element = ",matrix_element
    return
  endif

  ! If difference is exactly one
  if (difference .eq. 1) then
    matrix_element=ham_on_2x2(ind_bra(which_diff_1),ind_ket(which_diff_1))  ! Off diagonal on patch term
    if (spr) write (6,*) "Off diagonal on patch term being evaluated"
    if (spr) write (6,*) "matrix_element = ",matrix_element
    return
  endif

  if (spr) write (6,*) "Off patch term being evaluated"

 ! If difference is exactly 2 then check it connected or not in real space. If not, return 0. If connected, get the orientation of blocks.
 call check_nbr(l_x/2,l_y/2,pbc,which_diff_1,which_diff_2,nbr_types,num_types)       ! In more tools.f90
 if (spr) write (6,*) "Number of ways in which the blocks are connected",num_types
 if (num_types .eq.0) then
       matrix_element=0._rk
       return
 endif

 if (spr) then
     write (6,*) "which_diff_1",which_diff_1
     write (6,*) "which_diff_2",which_diff_2
     call flush(6)
 endif

 up=.false.
 electron_lost=.true.
 ! If connected in real space, check if the number of electrons differs exactly by 1. If not, matrix element =0
 ! If yes, determine what kind of electron it is (up vs down)
 up_diff_1=bra_nup_block(which_diff_1)-ket_nup_block(which_diff_1)
 up_diff_2=bra_nup_block(which_diff_2)-ket_nup_block(which_diff_2)
 dn_diff_1=bra_ndn_block(which_diff_1)-ket_ndn_block(which_diff_1)
 dn_diff_2=bra_ndn_block(which_diff_2)-ket_ndn_block(which_diff_2)

 if (spr) then
     write (6,*) "up_diff_1",up_diff_1
     write (6,*) "up_diff_2",up_diff_2
     write (6,*) "dn_diff_1",dn_diff_1
     write (6,*) "dn_diff_2",dn_diff_2
     call flush(6)
 endif

 if ((dn_diff_1 .eq. 0) .and. (dn_diff_2 .eq. 0)) then
     if ((up_diff_1 .eq. -1) .and. (up_diff_2 .eq. 1)) then
        up=.true.
        electron_lost=.true.
     elseif ((up_diff_1 .eq. 1) .and. (up_diff_2 .eq. -1)) then
        up=.true.
        electron_lost=.false.
     else
        matrix_element=0._rk
        return
     endif
 elseif ((up_diff_1 .eq. 0) .and. (up_diff_2 .eq. 0)) then
     if ((dn_diff_1 .eq. -1) .and. (dn_diff_2 .eq. 1)) then
        up=.false.
        electron_lost=.true.
     elseif ((dn_diff_1 .eq. 1) .and. (dn_diff_2 .eq. -1)) then
        up=.false.
        electron_lost=.false.
     else
        matrix_element=0._rk
        return
     endif
 else
    matrix_element=0._rk
    return
 endif

 if (spr) then
     write (6,*) "up=",up
     call flush(6)
 endif

 matrix_element=0._rk
 do n=1,num_types
    call get_matrix_element_from_table(ind_ket(which_diff_1),ind_ket(which_diff_2),ind_bra(which_diff_1),ind_bra(which_diff_2), nbr_types(n), electron_lost, up, matrix_element_tmp)
    matrix_element=matrix_element+matrix_element_tmp
 enddo

 ! Multiply matrix element by a sign which takes into account other block occupations
 n_tot=0
 if (up) then
    do i=which_diff_1+1,which_diff_2-1
        n_tot=n_tot+ket_nup_block(i)
    enddo
 else
    do i=which_diff_1+1,which_diff_2-1
        n_tot=n_tot+ket_ndn_block(i)
    enddo
 endif

 if (mod(n_tot,2) .ne. 0) matrix_element=-matrix_element

 end subroutine hamiltonian_hubbard_dm
!=========================================================================

! =====================================================================================
  subroutine generate_k_vectors()
! =====================================================================================
  ! ---------------------------------------------------------------------------
  ! Description   : Fills the array k_vectors with columns representing
  !                 possible k-vectors for 2D lattice in k-space. These are
  !                 -Lx+1,-Lx+3,...Lx-1 (Lx odd) or -Lx+2,-Lx+4,...Lx (Lx even)
  !
  ! Author        : A. Holmes, 20 Jun 2011 (similar to Frank's code in heg.f90)
  ! ---------------------------------------------------------------------------

  use types, only : num_words
  implicit none

  ! local variables
  integer                               :: ctr,i,i2,j,k,m,n,ktmp(2)
  integer                               :: tmpkx,tmpky,upleft,dnleft,e_deg,deg_sites
  integer,dimension(nsites)             :: tmp_degeneracy_list
  real(rk)                              :: tmp_energies(nsites)
  integer                               :: tmp_vectors(2,nsites)
  real(rk)                              :: ubyn_tmp,tmp_energy
  real(rk)                              :: upenergies(nsites),dnenergies(nsites)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable              :: conn_dets_up(:),conn_dets_dn(:)
  type(ik_vec)                           :: k_hf_up_disp(1),k_hf_dn_disp(1),rep_up,rep_dn
  type(ik_vec),allocatable               :: tmp_deg_up(:),tmp_deg_dn(:)
#else
  integer(ik), allocatable              :: conn_dets_up(:),conn_dets_dn(:)
  integer(ik)                           :: k_hf_up_disp(1),k_hf_dn_disp(1),rep_up,rep_dn
  integer(ik),allocatable               :: tmp_deg_up(:),tmp_deg_dn(:)
#endif
  integer                               :: nconn_hfdets
  integer                               :: up_site,dn_site,target_up_site
  real(rk),allocatable                  :: matrix_elements(:)
  integer,allocatable                   :: k_hf_con_up_mom(:,:),k_hf_con_dn_mom(:,:),k_hf_con_tot_mom(:,:)
  integer,allocatable                   :: k_energies_rounded(:)
  integer                               :: sortorder(nsites)
  real(rk),allocatable                  :: ham(:,:),eigenvalues(:)
  real(rk)                              :: tol ! tolerance for numerical energy difference between degenerate states.
  integer                               :: nnzero
  character(len=2)                      :: fmt,fmt2

  tol=1.d-15
  write (6,'(''Tolerance for numerical energy difference between degenerate states:'',es10.2)') tol
  call flush(6)

  ubyn=U/real(nsites)

  if (allocated(k_vectors)) then
    deallocate(k_vectors)
  endif

  if (allocated(k_energies)) then
    deallocate(k_energies)
  endif

  allocate(k_vectors(2,nsites))
  allocate(k_energies(nsites))

  do i=1,l_x
    do j=1,l_y
      k_vectors(1,l_y*(i-1)+j) = -l_x+2*i
      k_vectors(2,l_y*(i-1)+j) = -l_y+2*j
    enddo
  enddo

  if (mod(l_x,2)==1) k_vectors(1,:)=k_vectors(1,:)-1
  if (mod(l_y,2)==1) k_vectors(2,:)=k_vectors(2,:)-1

  !write (6,*) "K-vectors:" !call print_int_matrix(2,nsites,k_vectors)

  do i=1,nsites
    if (l_y==1) then! assume that l_x and l_y cannot both be 1
      k_energies(i) = -2._rk*t*(cos(pi*k_vectors(1,i)/real(l_x)))
    elseif (l_x==1) then
      k_energies(i) = -2._rk*t*(cos(pi*k_vectors(2,i)/real(l_y)))
    else
      k_energies(i) = -2._rk*t*(cos(pi*k_vectors(1,i)/real(l_x))+cos(pi*k_vectors(2,i)/real(l_y)))
    endif
  enddo

  !write (6,*) "Rounded energies:" !call print_real_matrix(1,nsites,k_energies)!_rounded)! sort k_vectors,k_energies by energy

  ! Following is order N^2 sort. Could be sped up, but N=nsites is relatively small.
  tmp_energies(:) = k_energies(:)
  tmp_vectors(:,:) = k_vectors(:,:)
  do i=1,nsites
    do j=1,nsites
      if (tmp_energies(j)==minval(tmp_energies)) then
        tmp_energies(j) = maxval(tmp_energies) + 1.0_rk
        sortorder(i)=j
        !k_energies(i) = k_energies(j)
        !k_vectors(:,i) = tmp_vectors(:,j)
        exit
      endif
    enddo
  enddo

  tmp_energies(:) = k_energies(:)
  tmp_vectors(:,:) = k_vectors(:,:)
  do i=1,nsites
    k_energies(i) = tmp_energies(sortorder(i))
    k_vectors(:,i) = tmp_vectors(:,sortorder(i))
  enddo

  allocate(k_energies_rounded(nsites))
  do i=1,nsites
    k_energies_rounded(i) = nint(k_energies(i))
  enddo

  write (6,*) "Sorted energies:"
  call print_real_matrix(1,nsites,k_energies)

  write (6,*) "Sorted, rounded energies:"
  call print_int_matrix(1,nsites,k_energies_rounded)

  write (6,*) "Sorted k-vectors:"
  call print_int_matrix(2,nsites,k_vectors)

  ! create 'map' that takes you from one site to another given a change in momentum:
  ! newsite = kmap(oldsite,delta_kx,delta_ky)
  if ((l_x>1 .and. l_y>1) .or. .true.) then
    if (allocated(kmap)) then
      deallocate(kmap)
    endif
    allocate(kmap(nsites,2*l_x-1,2*l_y-1))
    kmap(:,:,:)=0
    do i=1,nsites
      do j=1,2*l_x-1
        tmpkx = k_vectors(1,i)-2*l_x+2*j
        do k=1,2*l_y-1
          tmpky = k_vectors(2,i)-2*l_y+2*k
          do m=1,nsites
            ! allow move to -l_x
            !if (tmpkx==-l_x)  tmpkx=l_x
            !if (tmpkx==-l_y)  tmpkx=l_y
            if (mod(tmpkx-k_vectors(1,m),2*l_x)==0 .and. mod(tmpky-k_vectors(2,m),2*l_y)==0) then
              kmap(i,j,k) = m
              !write (6,*) "kmap(i,j,k)=m; i,j,k,m=",i,j,k,m
              !call flush(6)
              exit
            endif
          enddo
        enddo
      enddo
    enddo
  endif

  if (space_sym) then
        allocate(c4_map(3,nsites))
        allocate(reflection_map(nsites))

        !Real space - four fold + reflection about y axis
        write (6,*) "Creating real space symmetry maps (not used in the code right now)"
        call flush(6)
        call create_rspace_sym_maps(l_x,c4_map,reflection_map)
        write (6,*) "Finished Creating real space symmetry maps (not used in the code right now)"
        call flush(6)

        do i=1,3
            write (6,*) "C^",i,"(real space)"
            call print_map_on_square(l_x,c4_map(i,:))
        enddo
        write (6,*) "Reflection map (real space)"
        call print_map_on_square(l_x,reflection_map)

        !Momentum space - four fold + reflection about y axis
        write (6,*) "Creating k space symmetry maps "
        call flush(6)
        call create_kspace_sym_maps(k_vectors,c4_map,reflection_map)
        write (6,*) "Finished Creating real space symmetry maps (not used in the code right now)"
        call flush(6)

        do i=1,3
            write (6,*) "C^",i,"(k space)"
            call print_map_on_square(l_x,c4_map(i,:))
        enddo
        write (6,*) "Reflection map (k space)"
        call print_map_on_square(l_x,reflection_map)

   endif

  ! find k_hf_up, k_hf_dn, which are bit-packed representations of lowest energy hf states

  upleft=nup
  dnleft=ndn

  k_hf_up=0_ik
  k_hf_dn=0_ik

  upenergies=k_energies
  dnenergies=k_energies

  ktot(:)=0


  tmp_degeneracy_list(:)=0

  !find which empty spot has lowest energy, subject to momentum conservation.
  !assume total momentum of lowest energy state is what we expect it to be for now.
  !also, assume ground state is nondegenerate
  !and set ktot = the total momentum of the hf state
  tmp_energy = 0._rk
  n=1
  do while (upleft>0)
    do i=1,nsites
      if (abs(upenergies(i)-minval(upenergies))<tol) then
        k_hf_up = ibset(k_hf_up,i-1)
        upenergies(i) = maxval(upenergies) + 1.0 ! excludes site i from being considered again
        if (wf_type .eq. 'sym' .or. wf_type .eq. 'gutz_multi' .or. wf_type .eq. 'cgutz_multi') then
          n=n+1
          if (k_energies(i) /= tmp_energy)  n=1
          tmp_energy = k_energies(i)
          m=i
        endif
        ktot(1) = ktot(1) + k_vectors(1,i)
        ktot(2) = ktot(2) + k_vectors(2,i)
        exit
      endif
    enddo
    upleft=upleft-1
  enddo

  if (nup .eq. ndn) then  ! A hack for open shell

    allocate(k_hf_deg_up(1))
    allocate(k_hf_deg_dn(1))
    allocate(degeneracy_list(1))

    if ((wf_type .eq. 'sym') .or. (wf_type .eq. 'gutz_multi') .or. (wf_type .eq. 'cgutz_multi')) then

      do i=1,nsites ! this makes the total momentum 0
        if (btest(k_hf_up,i-1)) then
          do j=1,nsites
            if ((k_vectors(1,i)+k_vectors(1,j))==0 .and. (k_vectors(2,i)+k_vectors(2,j))==0) then
              k_hf_dn=ibset(k_hf_dn,j-1)
              ktot(1) = ktot(1) + k_vectors(1,j)
              ktot(2) = ktot(2) + k_vectors(2,j)
              dnleft=dnleft-1
              exit
            endif
          enddo
        endif
      enddo

      if (l_x==2 .and. l_y==2 .and. nup==2 .and. ndn==2)  ktot(:)=2
      if (l_x==4 .and. l_y==2 .and. nup==5 .and. ndn==5)  then
            ktot(1)=0
            ktot(2)=4
      endif
      if (l_x==2 .and. l_y==4 .and. nup==5 .and. ndn==5)  then
            ktot(1)=4
            ktot(2)=0
      endif

      ! create list of degenerate lattice sites, for now assuming that nup=ndn
      e_deg = n ! the number of electrons that can be placed in the degenerate lattice states

      j=0
      do k=1,nsites
        if (abs(k_energies(k)-k_energies(m))<tol) then
          j=j+1
          tmp_degeneracy_list(j)=k
        endif
      enddo
      deg_sites = j

      deallocate(degeneracy_list)
      allocate(degeneracy_list(deg_sites))
      do i=1,deg_sites
        degeneracy_list(i) = tmp_degeneracy_list(i)
      enddo

      if (e_deg.ne.deg_sites) then
        write (6,*) "Number of electrons to be placed in degenerate lattice sites:",e_deg
        write (6,*) "Number of degenerate lattice sites:",deg_sites
        write(6,*) "Degenerate lattice sites:",degeneracy_list
        call flush(6)
      endif

      allocate(tmp_deg_up((n_choose_k(deg_sites,e_deg))**2))
      allocate(tmp_deg_dn((n_choose_k(deg_sites,e_deg))**2))
      tmp_deg_up(:)=0
      tmp_deg_dn(:)=0
      m=0
      do i=1,2**deg_sites-1
        n=0
        do j=1,deg_sites
          if (btest(i,j-1))  n=n+1
        enddo
        if (n==e_deg) then
          do i2=1,2**deg_sites-1
            n=0
            do j=1,deg_sites
              if (btest(i2,j-1))  n=n+1
            enddo
            if (n==e_deg) then
              ktmp(:)=0
              do j=1,deg_sites
                if (btest(i,j-1)) then
                  ktmp(1) = ktmp(1) + k_vectors(1,degeneracy_list(j))
                  ktmp(2) = ktmp(2) + k_vectors(2,degeneracy_list(j))
                endif
                if (btest(i2,j-1)) then
                  ktmp(1) = ktmp(1) + k_vectors(1,degeneracy_list(j))
                  ktmp(2) = ktmp(2) + k_vectors(2,degeneracy_list(j))
                endif
              enddo
              if (mod(ktmp(1)-ktot(1),2*l_x)==0 .and. mod(ktmp(2)-ktot(2),2*l_y)==0) then
                m=m+1
                tmp_deg_up(m) = k_hf_up
                tmp_deg_dn(m) = k_hf_dn
                do j=1,deg_sites
                  tmp_deg_up(m) = ibclr(tmp_deg_up(m),degeneracy_list(j)-1)
                  if (btest(i,j-1))  tmp_deg_up(m)=ibset(tmp_deg_up(m),degeneracy_list(j)-1)
                  tmp_deg_dn(m) = ibclr(tmp_deg_dn(m),degeneracy_list(j)-1)
                  if (btest(i2,j-1))  tmp_deg_dn(m)=ibset(tmp_deg_dn(m),degeneracy_list(j)-1)
                enddo
              endif
            endif
          enddo
        endif
      enddo

      ! then allocate k_hf_deg_up
      deallocate(k_hf_deg_up)
      deallocate(k_hf_deg_dn)
      allocate(k_hf_deg_up(m))
      allocate(k_hf_deg_dn(m))
      do i=1,m
        k_hf_deg_up(i) = tmp_deg_up(i)
        k_hf_deg_dn(i) = tmp_deg_dn(i)
      enddo
      deallocate(tmp_deg_up)
      deallocate(tmp_deg_dn)
      ndeg = m

      write (6,*) "Degenerate ground states:"
      write (6,*) "Up:"
      call print_bin(ndeg,nsites,k_hf_deg_up)
      write (6,*) "Dn:"
      call print_bin(ndeg,nsites,k_hf_deg_dn)
      call flush(6)

      write (6,*) "Number of degenerate states :",ndeg
      call flush(6)

      ! now diagonalize hamiltonian constructed between the states with same energy and momentum as the HF ground state to find coefficient list
      allocate(ham(ndeg,ndeg))
      do i=1,ndeg
        do j=1,ndeg
          call hamiltonian_hubbard_k(k_hf_deg_up(i),k_hf_deg_dn(i),k_hf_deg_up(j),k_hf_deg_dn(j),ham(i,j))
        enddo
      enddo

      if (ndeg .lt. 100) then
          write (6,'(/,''Hamiltonian matrix between degenerate HF ground states:'')')
          call print_real_matrix(ndeg,ndeg,ham)
          call flush(6)
      else
          write (6,'(/,''Not printing Hamiltonian matrix between degenerate HF ground states because of big size'')')
      endif

      allocate (eigenvalues(ndeg))
      call real_symmetric_diagonalize_ow_ham(ndeg,ham,eigenvalues)

      write(6,'(''Size of degenerate space is '',i4)') ndeg

      write (6,'(/,''Eigenvalues:'')')
      call print_real_matrix(1,ndeg,eigenvalues)
      call flush(6)

      write (6,'(/,''Eigenvectors:'')')
      call print_real_matrix(ndeg,ndeg,ham)
      call flush(6)

      allocate(c_sym_psi_t(ndeg))
      do i=1,ndeg
        c_sym_psi_t(i) = ham(i,1)
      enddo

      write (6,'(/,''Ground state eigenvector coefficients by order in H matrix:'')')
      do i=1,size(k_hf_deg_up,1)
          write(6,'(2i10,6g14.6)') k_hf_deg_up(i), k_hf_deg_dn(i), c_sym_psi_t(i)
      enddo

     if (l_x .eq. 4 .and. l_y .eq. 4 .and. nup .eq. 4 .and. ndn .eq. 4) then! .and. nshells.eq. 0) then
       c_sym_psi_t(1)=1._rk
       c_sym_psi_t(2)=-1._rk
       c_sym_psi_t(3)=-1._rk
       c_sym_psi_t(4)=1._rk
     endif

     if (l_x .eq. 4 .and. l_y .eq. 4 .and. nup .eq. 2 .and. ndn .eq. 2) then! .and. nshells.eq. 0) then
       c_sym_psi_t(1)=-1._rk
       c_sym_psi_t(2)=1._rk
       c_sym_psi_t(3)=1._rk
       c_sym_psi_t(4)=-1._rk
     endif

      write (6,'(/,''Ground state eigenvector coefficients by weights:'')')
      call sort(c_sym_psi_t, k_hf_deg_up, k_hf_deg_dn)

      do i=1,size(k_hf_deg_up,1)
          write(6,'(2i10,6g14.6)') k_hf_deg_up(i), k_hf_deg_dn(i), c_sym_psi_t(i)
      enddo
      !call print_real_matrix(1,ndeg,c_sym_psi_t)
      call flush(6)
    endif
    if (wf_type .eq. 'gutz')  then ! put dn electrons in lowest energy positions (if nup=ndn, this is the same configuration as the up electons)

    !  do while (dnleft>0)
    !    do i=1,nsites
    !      if (dnenergies(i)==minval(dnenergies)) then
    !        k_hf_dn = ibset(k_hf_dn,i-1)
    !        dnenergies(i) = maxval(dnenergies) + 1.0
    !        ktot(1) = ktot(1) + k_vectors(1,i)
    !        ktot(2) = ktot(2) + k_vectors(2,i)
    !        exit
    !      endif
    !    enddo
    !    dnleft=dnleft-1
    !  enddo

      if (l_x .eq. l_y) then
          do i=1,nsites ! this makes the total momentum 0
            if (btest(k_hf_up,i-1)) then
              do j=1,nsites
                if ((k_vectors(1,i)+k_vectors(1,j))==0 .and. (k_vectors(2,i)+k_vectors(2,j))==0) then
                  k_hf_dn=ibset(k_hf_dn,j-1)
                  ktot(1) = ktot(1) + k_vectors(1,j)
                  ktot(2) = ktot(2) + k_vectors(2,j)
                  dnleft=dnleft-1
                  exit
                endif
              enddo
            endif
          enddo
      else
         k_hf_dn=k_hf_up
         do j=1,nsites
           if (btest(k_hf_dn,j-1)) then
                  ktot(1) = ktot(1) + k_vectors(1,j)
                  ktot(2) = ktot(2) + k_vectors(2,j)
           endif
         enddo
      endif
    endif

  !elseif (.true.) then ! put dn electrons in lowest energy positions (if nup=ndn, this is the same configuration as the up electons)
!  elseif (wf_type .eq. 'gutz')  then ! put dn electrons in lowest energy positions (if nup=ndn, this is the same configuration as the up electons)
!    do while (dnleft>0)
!      do i=1,nsites
!        if (dnenergies(i)==minval(dnenergies)) then
!          k_hf_dn = ibset(k_hf_dn,i-1)
!          dnenergies(i) = maxval(dnenergies) + 1.0
!          ktot(1) = ktot(1) + k_vectors(1,i)
!          ktot(2) = ktot(2) + k_vectors(2,i)
!          exit
!        endif
!      enddo
!      dnleft=dnleft-1
!    enddo

  else ! this is for the 4477 case with ground state momentum = (0,pi)
    do while (dnleft>2)
      do i=1,nsites
        if (dnenergies(i)==minval(dnenergies)) then
          k_hf_dn = ibset(k_hf_dn,i-1)
          dnenergies(i) = maxval(dnenergies) + 1.0
          ktot(1) = ktot(1) + k_vectors(1,i)
          ktot(2) = ktot(2) + k_vectors(2,i)
          exit
        endif
      enddo
      dnleft=dnleft-1
    enddo
    j=0
    do while (dnleft>0)
      do i=1,nsites
        if (dnenergies(i)==minval(dnenergies)) then
          j=j+1
          if (j>3) then
            k_hf_dn = ibset(k_hf_dn,i-1)
            dnenergies(i) = maxval(dnenergies) + 1.0
            ktot(1) = ktot(1) + k_vectors(1,i)
            ktot(2) = ktot(2) + k_vectors(2,i)
            exit
          endif
        endif
      enddo
      dnleft=dnleft-1
    enddo

  endif

  write (6,'(/,''Total momentum of HF state (after gutz):'',2i5)') ktot
! call print_int_matrix(2,1,ktot)
  call flush(6)

  ! now find k_hf_con_up,k_hf_con_dn, which are lists of dets connected to k_hf_up,k_hf_dn
  ! (this is for the first order k-space gutzwiller)

  ! k-space gutz = k_hf + ln(g) * sum(i) [k_hf_matrix _elements*k_hf_con]

  ubyn_tmp=ubyn
  ubyn=1._rk/real(nsites)

  if (allocated(k_hf_con_up)) deallocate(k_hf_con_up)
  if (allocated(k_hf_con_dn)) deallocate(k_hf_con_dn)
  if (allocated(k_hf_matrix_elements)) deallocate(k_hf_matrix_elements)

  allocate(k_hf_con_up((nup*ndn*(nsites-nup))+1))
  allocate(k_hf_con_dn((nup*ndn*(nsites-nup))+1))
  allocate(k_hf_matrix_elements((nup*ndn*(nsites-nup))+1))

  write(6,*) k_hf_up
  write(6,*) k_hf_dn
  call flush(6)

  call find_connected_dets_hubbard_k(k_hf_up,k_hf_dn,n_k_hf_con,k_hf_con_up,k_hf_con_dn,nsites,k_hf_matrix_elements)

  k_hf_matrix_elements(1) = ubyn * real(nup) * real(ndn) ! why this???? Commented by Hitesh - I guess for perturbation code???

  ubyn=ubyn_tmp

  if (space_sym) then
    call get_rep_only(c4_map,reflection_map,z,p,k_hf_up,k_hf_dn,rep_up,rep_dn)
    write (6,*) "rep_up=",rep_up,"rep_dn=",rep_dn
    call flush(6)
    k_hf_up=rep_up
    k_hf_dn=rep_dn
    call symmetry_reduce_hubbardk(n_k_hf_con,k_hf_con_up,k_hf_con_dn,rep_k_hf_con_up,rep_k_hf_con_dn)
    ! Sort the list of reps
    call sort(rep_k_hf_con_up,rep_k_hf_con_dn)
    write (6,*) "Finished sorting reps"
    call flush(6)
    if (allocated (k_hf_matrix_elements)) deallocate (k_hf_matrix_elements)
    n_k_hf_con_sym=size(rep_k_hf_con_up,1)
    allocate (k_hf_matrix_elements(n_k_hf_con))
    write (6,*) "Finished allocating k_hf_matrix_elements"
    call flush(6)
    write (6,*)
    write (6,*) "Representatives connected to the Hartree Fock state"
    call flush(6)
    write (6,*) "n_k_hf_con_sym=",n_k_hf_con_sym
    write (fmt, '(i2)') 9*num_words ! Why 9?
    write (fmt2, '(i2)') 2*num_words
    do i=1,n_k_hf_con_sym
        write(6,'(' // trim(fmt) // 'i6)') k_hf_up,k_hf_dn,rep_k_hf_con_up(i),rep_k_hf_con_dn(i)
        call hamiltonian_hubbard_k_space_sym(k_hf_up,k_hf_dn,rep_k_hf_con_up(i),rep_k_hf_con_dn(i),k_hf_matrix_elements(i),nnzero)
        write(6,'(' // trim(fmt2) // 'i6,f13.8)') rep_k_hf_con_up(i),rep_k_hf_con_dn(i),k_hf_matrix_elements(i)
        call flush(6)
    enddo
  endif

  k_hf_up_mom(:)=0
  k_hf_dn_mom(:)=0

  do k=1,nsites
    if (btest(k_hf_up,k-1)) then
      k_hf_up_mom(1) = k_hf_up_mom(1) + k_vectors(1,k)
      k_hf_up_mom(2) = k_hf_up_mom(2) + k_vectors(2,k)
    endif
  enddo

  do k=1,nsites
    if (btest(k_hf_dn,k-1)) then
      k_hf_dn_mom(1) = k_hf_dn_mom(1) + k_vectors(1,k)
      k_hf_dn_mom(2) = k_hf_dn_mom(2) + k_vectors(2,k)
    endif
  enddo

  k_hf_up_disp(1)=k_hf_up
  k_hf_dn_disp(1)=k_hf_dn

  write (6,*) "HF Up Det:"
  call print_bin(1,nsites,k_hf_up_disp)
  write (6,*) "HF Dn Det:"
  call print_bin(1,nsites,k_hf_dn_disp)
  write (6,*) "HF Up Momentum (multiply by pi/L_x and pi/L_y):"
  call print_int_matrix(1,2,k_hf_up_mom)
  write (6,*) "HF Dn Momentum (multiply by pi/L_x and pi/L_y):"
  call print_int_matrix(1,2,k_hf_dn_mom)
  write (6,*) "HF Total Momentum (multiply by pi/L_x and pi/L_y):"
  call print_int_matrix(1,2,k_hf_up_mom+k_hf_dn_mom)
  call flush(6)

  allocate(k_hf_con_up_mom(2,n_k_hf_con))
  allocate(k_hf_con_dn_mom(2,n_k_hf_con))
  allocate(k_hf_con_tot_mom(2,n_k_hf_con))

  k_hf_con_up_mom(:,:)=0
  k_hf_con_dn_mom(:,:)=0
  k_hf_con_tot_mom(:,:)=0

  do i=1,n_k_hf_con
    do j=1,nsites
      if (btest(k_hf_con_up(i),j-1)) then
        k_hf_con_up_mom(1,i) = k_hf_con_up_mom(1,i) + k_vectors(1,j)
        k_hf_con_up_mom(2,i) = k_hf_con_up_mom(2,i) + k_vectors(2,j)
      endif
      if (btest(k_hf_con_dn(i),j-1)) then
        k_hf_con_dn_mom(1,i) = k_hf_con_dn_mom(1,i) + k_vectors(1,j)
        k_hf_con_dn_mom(2,i) = k_hf_con_dn_mom(2,i) + k_vectors(2,j)
      endif
      k_hf_con_tot_mom(:,i) = k_hf_con_up_mom(:,i) + k_hf_con_dn_mom(:,i)
    enddo
  enddo

  k_hf_up_mom(:)=0
  k_hf_dn_mom(:)=0

  do k=1,nsites
    if (btest(k_hf_up,k-1)) then
      k_hf_up_mom(1) = k_hf_up_mom(1) + k_vectors(1,k)
      k_hf_up_mom(2) = k_hf_up_mom(2) + k_vectors(2,k)
    endif
    if (btest(k_hf_dn,k-1)) then
      k_hf_dn_mom(1) = k_hf_dn_mom(1) + k_vectors(1,k)
      k_hf_dn_mom(2) = k_hf_dn_mom(2) + k_vectors(2,k)
    endif
  enddo

  allocate(up_hf_location(nup))
  allocate(up_hf_empty_location(nsites-nup))
  allocate(dn_hf_location(ndn))

  ctr=1
  do i=1,nsites
    if (btest(k_hf_up,i-1)) then
        up_hf_location(ctr)=i
        ctr=ctr+1
    endif
  enddo

  ctr=1
  do i=1,nsites
    if (btest(k_hf_up,i-1) .eqv. .false.) then
        up_hf_empty_location(ctr)=i
        ctr=ctr+1
    endif
  enddo

  ctr=1
  do i=1,nsites
    if (btest(k_hf_dn,i-1)) then
        dn_hf_location(ctr)=i
        ctr=ctr+1
    endif
  enddo

  allocate(conn_dets_up((nup*ndn*(nsites-nup))+1))
  allocate(conn_dets_dn((nup*ndn*(nsites-nup))+1))
  allocate(matrix_elements((nup*ndn*(nsites-nup))+1))
  call find_connected_dets_hubbard_k(k_hf_up,k_hf_dn,nconn_hfdets,conn_dets_up,conn_dets_dn,nsites,matrix_elements)
  allocate(hf_matrix_elements(nup,ndn,nsites-nup))
  hf_matrix_elements(:,:,:)=0._rk

  do i=1,nconn_hfdets
       if ((conn_dets_up(i) .eq. k_hf_up) .and. (conn_dets_dn(i) .eq. k_hf_dn) ) then
           k_hf_energy=matrix_elements(i)
       else
           do j=1,nup
              if (btest(conn_dets_up(i),up_hf_location(j)-1) .eqv. .false.) up_site=j
           enddo

           do j=1,nsites-nup
              if (btest(conn_dets_up(i),up_hf_empty_location(j)-1)) target_up_site=j
           enddo

           do j=1,ndn
              if (btest(conn_dets_dn(i),dn_hf_location(j)-1) .eqv. .false.) dn_site=j
           enddo

           hf_matrix_elements(up_site,dn_site,target_up_site)=matrix_elements(i)
       endif
  enddo

  write (6,*) "End of generate_k_vectors"
  call flush(6)

  end subroutine generate_k_vectors
!=========================================================================

  ! =====================================================================================

  ! =====================================================================================
  subroutine hamiltonian_hubbard_k(up_bra,dn_bra,up_ket,dn_ket, matrix_element,is_connected)
  ! =====================================================================================
  ! ---------------------------------------------------------------------------
  ! Description   : Generate Hamiltonian matrix element for 2D lattice in k-space
  ! Note          : If is_connected is present, don't check to see if dets are connected.
  ! Author        : A Holmes, 19 Jun 2011 (k-space version of Hitesh's code above)
  ! Modified      : A Holmes, 23 May 2012. Simplified and sped up.
  ! ---------------------------------------------------------------------------
    use tools, only : permutation_factor,count_excitations

    implicit none

    real(rk),   intent(out)     :: matrix_element
    logical,optional,intent(in) :: is_connected
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in)      :: up_bra,dn_bra,up_ket,dn_ket
    type(ik_vec)                 :: det
#else
    integer(ik),intent(in)      :: up_bra,dn_bra,up_ket,dn_ket
    integer(ik)                 :: det
#endif

    integer                     :: i

    !=======================================================
    ! Diagonal term
    !=======================================================
    if (up_bra==up_ket .and. dn_bra==dn_ket) then
      matrix_element = ubyn*nup*ndn
      det = up_bra
      do while (det .ne. 0)
        i = trailz(det)+1
        matrix_element = matrix_element + k_energies(i)
        det = ibclr(det,i-1)
      enddo
      det = dn_bra
      do while (det .ne. 0)
        i = trailz(det)+1
        matrix_element = matrix_element + k_energies(i)
        det = ibclr(det,i-1)
      enddo
      return
    endif

    !=====================================================
    ! Terms separated by 1 up and 1 down electron hop
    !=====================================================

    if (present(is_connected)) then
      matrix_element = ubyn*permutation_factor(up_bra,up_ket)*permutation_factor(dn_bra,dn_ket)
    else
      if (is_connected_hubbard_fast(up_bra,dn_bra,up_ket,dn_ket)) then
        matrix_element = ubyn*permutation_factor(up_bra,up_ket)*permutation_factor(dn_bra,dn_ket)
      else
        matrix_element = 0._rk
      endif
    endif

  end subroutine hamiltonian_hubbard_k

  ! ============================================================================================================
  subroutine hamiltonian_hubbard_k_space_sym(up_bra,dn_bra,up_ket,dn_ket,combined_matrix_element,num_non_zero)
  ! ============================================================================================================
  ! -----------------------------------------------------------------------------------------------------------
  ! Description   : Bit packed version to generate Hamiltonian matrix element for 2D lattice in k-space
  !                 Let us assume that the ket is a representative and the bra is not
  !                 Send out the representative bra
  ! Author        : H.J.Changlani, April 9 2012
  ! -----------------------------------------------------------------------------------------------------------

  implicit none

  ! Dummy
  real(rk),intent(out)        :: combined_matrix_element
  integer,intent(out)         :: num_non_zero
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in)      :: up_ket,dn_ket
  type(ik_vec),intent(inout)   :: up_bra,dn_bra
  type(ik_vec)                 :: new_dets_up(16),new_dets_dn(16),rep_up,rep_dn
  type(ik_vec)                 :: tmp_distinct_up(16),tmp_distinct_dn(16)
#else
  integer(ik),intent(in)      :: up_ket,dn_ket
  integer(ik),intent(inout)   :: up_bra,dn_bra
  integer(ik)                 :: new_dets_up(16),new_dets_dn(16),rep_up,rep_dn
  integer(ik)                 :: tmp_distinct_up(16),tmp_distinct_dn(16)
#endif

  ! Local
  real(rk)                    :: norm_r,norm_s,matrix_element,phase_w_rep,phases(16)
  integer                     :: i,num_distinct,location,counter
  logical                     :: entry_present

  call generate_fourfold_k_configs_efficient(c4_map,reflection_map,z,p,nup,ndn,up_ket,dn_ket,new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm_r,num_distinct)
  call generate_fourfold_k_configs_efficient(c4_map,reflection_map,z,p,nup,ndn,up_bra,dn_bra,new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm_s,num_distinct)

  up_bra=rep_up
  dn_bra=rep_dn

  num_non_zero=0
  counter=0
  combined_matrix_element=0._rk

  tmp_distinct_up(:)=0_ik
  tmp_distinct_dn(:)=0_ik

  if (abs(norm_r) .gt. 1.0e-10 .and. abs(norm_s) .gt. 1.0e-10) then
    do i=1,16
        if (counter .eq. num_distinct) exit
        call is_in_list_mod(new_dets_up(i),new_dets_dn(i),tmp_distinct_up,tmp_distinct_dn,counter,entry_present,location)
        if (entry_present .eqv. .false.) then
            counter=counter+1
            tmp_distinct_up(counter)=new_dets_up(i)
            tmp_distinct_dn(counter)=new_dets_dn(i)
            call hamiltonian_hubbard_k(new_dets_up(i),new_dets_dn(i),up_ket,dn_ket,matrix_element)
            if (abs(matrix_element) .gt. 1.0e-10) then
                combined_matrix_element=combined_matrix_element+(matrix_element*phases(i))                      ! Refer to Lauechli/HJC notes,SD-Frick paper
                num_non_zero=num_non_zero+1
            endif
         endif
    enddo
    combined_matrix_element=combined_matrix_element*phase_w_rep*(norm_s/norm_r)                                 ! Refer to Lauechli/HJC notes,SD-Frick paper
  endif

 end subroutine hamiltonian_hubbard_k_space_sym

!======================================================================================================================
  subroutine off_diagonal_move_hubbard(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j, diagonal_matrix_element, vmc_on)
!======================================================================================================================
    !------------------------------------------------------------------------------------------------------------------
    ! Description : Perform an offdiagonal uniform move in SQMC for Hubbard model
    !               starting at det_i which is split into det_i_up and det_i_dn.
    !		        det_i has weight = +-1. If move fails, for any reason,
    !               weight_j is returned as 0. Otherwise, return new determinant,
    !               det_j, and weight of incoming det, H_ii (diagonal element of Hamiltonian).
    !               If run in the VMC mode certain checking conditions become different.
    !               The weight is not needed in VMC mode
    !               The diagonal_matrix element assumes the role of wf_ratio in VMC mode.
    ! Created     : H.J.Changlani, 28 Nov 2010
    ! Edited last : H.J.Changlani  May 31 2011
    !------------------------------------------------------------------------------------------------------------------
    implicit none

    !dummy variables
    real(rk), intent(out)   :: weight_j
    real(rk), optional,intent(out)   :: diagonal_matrix_element
    logical, optional       :: vmc_on
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up,det_i_dn
    type(ik_vec),intent (out):: det_j_up,det_j_dn
    type(ik_vec)              ::  det
#else
    integer(ik), intent(in) :: det_i_up,det_i_dn
    integer(ik),intent (out):: det_j_up,det_j_dn
    integer(ik)              ::  det
#endif

    !local variables
    real(rk)                 ::  proposal_prob_inv, acceptance_prob
    real(rk)                 ::  matrix_element
    integer                  ::  k,ctr,chosen_site,chosen_spin,target_site,target_k
    integer                  ::  nbr
    integer                  ::  allowed_nbrs(4)
    logical                  ::  allowed,success
    logical                  ::  vmc
    real(rk)                 ::  wf_1,wf_2

    vmc=.false.
    if (present(vmc_on)) vmc=vmc_on

    if (present(diagonal_matrix_element)) then
      diagonal_matrix_element=0._rk
      if (vmc .eqv. .false.) call hamiltonian_hubbard(det_i_up,det_i_dn,det_i_up,det_i_dn,diagonal_matrix_element)
    endif

    weight_j=0

    !initialize det_j to det_i
    det_j_up = det_i_up
    det_j_dn = det_i_dn

    success=.false.

    call choose_random_electron(det_i_up,det_i_dn,chosen_site,chosen_spin)

    if (chosen_spin .eq. 0) det=det_i_dn
    if (chosen_spin .eq. 1) det=det_i_up

    ctr=0
    do k=1,4
        allowed=.false.
        call get_nbr(l_x,l_y,pbc,chosen_site,k,nbr,allowed)
        if (ctr .eq. 0) then
            if ((allowed .eqv. .true.)) then
              if (btest(det,nbr-1) .eqv. .false.) then
                 ctr=ctr+1
                 allowed_nbrs(ctr)=nbr
              endif
            endif
        else
            !if (((allowed .eqv. .true.) .and. (nbr .ne. allowed_nbrs(ctr)))) then        ! wrapping around twice allowed
            if (allowed .eqv. .true.) then
              if (btest(det,nbr-1) .eqv. .false.) then
                 ctr=ctr+1
                 allowed_nbrs(ctr)=nbr
              endif
            endif
        endif
    enddo

    if (ctr .eq. 0) then
       success=.false.
       return
    else
       if (vmc) then
          target_k=random_int(4)
          if (target_k .gt. ctr) then
            success=.false.
            return
          endif
       else
        target_k=random_int(ctr)
       endif

       target_site=allowed_nbrs(target_k)
       if (btest(det,target_site-1) .eqv. .false.) then
           success=.true.
           det=ibset(det,target_site-1)
           det=ibclr(det,chosen_site-1)
       else
           success=.false.
           return
       endif
    endif

    if (success .eqv. .true. ) then
        if (chosen_spin .eq. 0) then
             det_j_dn=det
         else
             det_j_up=det
         endif

        if (vmc.and.present(diagonal_matrix_element)) then
            call wf_calc(det_i_up,det_i_dn,wf_1)
            call wf_calc(det_j_up,det_j_dn,wf_2)
            diagonal_matrix_element=wf_2/wf_1
        else
            call hamiltonian_hubbard(det_i_up,det_i_dn,det_j_up,det_j_dn,matrix_element)
            proposal_prob_inv=(real((nup+ndn)*ctr))
            acceptance_prob=abs(matrix_element)*tau*proposal_prob_inv
!           weight_j = int((acceptance_prob+rannyu())*sign(1._rk,-matrix_element))
            weight_j = acceptance_prob*sign(1._rk,-matrix_element)
        endif
   endif

  end subroutine off_diagonal_move_hubbard
  !===========================================================================================================

!======================================================================================================================
  subroutine off_diagonal_move_hubbard_dm(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j, diagonal_matrix_element)
!======================================================================================================================
!    !------------------------------------------------------------------------------------------------------------------
!    ! Description : Perform an offdiagonal uniform move in SQMC for Hubbard model in a unitary transformed basis
!    !               starting at det_i which is split into det_i_up and det_i_dn.
!    !		        det_i has weight = +-1. If move fails, for any reason,
!    !               weight_j is returned as 0. Otherwise, return new determinant,
!    !               det_j, and weight of incoming det, H_ii (diagonal element of Hamiltonian).
!    ! Created     : H.J.Changlani, Dec 06 2011
!    !------------------------------------------------------------------------------------------------------------------
    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in)  :: det_i_up,det_i_dn
    type(ik_vec),intent (out) :: det_j_up,det_j_dn
#else
    integer(ik), intent(in)  :: det_i_up,det_i_dn
    integer(ik),intent (out) :: det_j_up,det_j_dn
#endif
    real(rk), intent(out)    :: weight_j
    real(rk), optional, intent(out)    :: diagonal_matrix_element

    !local variables
    integer                  ::  i
    real(rk)                 ::  proposal_prob_inv, acceptance_prob,matrix_element
    real(rk)                 ::  rannyu
    integer                  ::  chosen_block,nbr_block,chosen_move,num_blocks
    integer                  ::  nup_block,ndn_block,nup_adj_block,ndn_adj_block
    integer                  ::  nup_new_block,ndn_new_block,nup_new_adj_block,ndn_new_adj_block
    integer                  ::  num_states,num_adj_states,num_connections
    integer (ik)             ::  chosen_vector,chosen_vector_adj,new_vector,new_adj_vector,shift
    integer                  ::  lost
    integer                  ::  nbr_types(4),num_nbr_types
    real(rk)                 ::  extra_factor
    logical                  ::  e_lost,up,success,allowed
    logical                  ::  spr

!   spr=.true.
    spr=.false.

    success=.false.
    up=.false.
    e_lost=.false.
    weight_j=0._rk

    !initialize det_j to det_i
    det_j_up = det_i_up
    det_j_dn = det_i_dn
    if (present(diagonal_matrix_element)) then
      diagonal_matrix_element=0._rk
      call hamiltonian_hubbard_dm(det_i_up,det_i_dn,det_j_up,det_j_dn,diagonal_matrix_element)
    endif
    call flush(6)
    ! Choose a block at random
    num_blocks=nsites/4
    chosen_block=random_int(num_blocks)

    ! Once block is chosen choose a number from 1 to 5
    chosen_move=random_int(5)

   ! 1 = On patch move
   ! 2 = Left move
   ! 3 = Right move
   ! 4 = Up move
   ! 5 = DOwn move

    chosen_vector=0_ik
    do i=1,4
        if (btest(det_i_up,4*(chosen_block-1)+i-1)) chosen_vector=ibset(chosen_vector,4+i-1)
        if (btest(det_i_dn,4*(chosen_block-1)+i-1)) chosen_vector=ibset(chosen_vector,i-1)
    enddo

    nup_block=nup_dm(chosen_vector)
    ndn_block=ndn_dm(chosen_vector)

    num_states=ndm_start_end(nup_block,ndn_block,2)-ndm_start_end(nup_block,ndn_block,1)+1
    shift=ndm_start_end(nup_block,ndn_block,1)
    num_connections=num_states-1 ! Number of off diagonal connections on patch

    if (chosen_move .eq. 1) then                                        ! On patch move
        do                                                              ! Select an off diagonal element on the patch
           if (num_states .gt. 1) then
                new_vector=random_int(num_states)+shift-1
                if (new_vector .ne. chosen_vector) exit
           else
                success=.false.
                return
           endif
        enddo
        ! Change det_j_up and det_j_dn
        do i=1,4
            if (btest(new_vector,i-1)) then
                det_j_dn=ibset(det_j_dn,4*(chosen_block-1)+i-1)
            else
                det_j_dn=ibclr(det_j_dn,4*(chosen_block-1)+i-1)
            endif
        enddo
        do i=5,8
            if (btest(new_vector,i-1)) then
                det_j_up=ibset(det_j_up,4*(chosen_block-1)+i-4-1)
            else
                det_j_up=ibclr(det_j_up,4*(chosen_block-1)+i-4-1)
            endif
        enddo
        success=.true.
        extra_factor=1._rk
    else                                                                ! Off patch moves in 4 directions
        if (rannyu()>0.5_rk) then                                      ! Electron lost or gained
            e_lost=.true.
            lost=-1
        else
            e_lost=.false.
            lost=1
        endif

        if (rannyu()>0.5_rk) then                                     ! Spin type
            up=.true.
        else
            up=.false.
        endif

        call get_nbr(l_x/2,l_y/2,pbc,chosen_block,chosen_move-1,nbr_block,allowed)

        if (allowed) then
            chosen_vector_adj=0_ik
            do i=1,4
                if (btest(det_i_up,4*(nbr_block-1)+i-1)) chosen_vector_adj=ibset(chosen_vector_adj,4+i-1)
                if (btest(det_i_dn,4*(nbr_block-1)+i-1)) chosen_vector_adj=ibset(chosen_vector_adj,i-1)
            enddo

            nup_adj_block=nup_dm(chosen_vector_adj)
            ndn_adj_block=ndn_dm(chosen_vector_adj)

            if (up) then
                nup_new_block=nup_block+lost
                ndn_new_block=ndn_block
                nup_new_adj_block=nup_adj_block-lost
                ndn_new_adj_block=ndn_adj_block
            else
                nup_new_block=nup_block
                nup_new_adj_block=nup_adj_block
                ndn_new_block=ndn_block+lost
                ndn_new_adj_block=ndn_adj_block-lost
            endif

            if ((nup_new_block .lt. 0) .or. (nup_new_block .gt. 4) .or. (ndn_new_block .lt. 0) .or. (ndn_new_block .gt. 4) .or. (nup_new_adj_block .lt. 0) .or. (nup_new_adj_block .gt. 4) .or. (ndn_new_adj_block .lt. 0) .or. (ndn_new_adj_block .gt. 4))then
                success=.false.
                return
            endif

            num_states=ndm_start_end(nup_new_block,ndn_new_block,2)-ndm_start_end(nup_new_block,ndn_new_block,1)+1
            ! Choose a state randomly from these
            shift=ndm_start_end(nup_new_block,ndn_new_block,1)
            new_vector=random_int(num_states)+shift-1

            num_adj_states=ndm_start_end(nup_new_adj_block,ndn_new_adj_block,2)-ndm_start_end(nup_new_adj_block,ndn_new_adj_block,1)+1
            ! Choose a state randomly from these
            shift=ndm_start_end(nup_new_adj_block,ndn_new_adj_block,1)
            new_adj_vector=random_int(num_adj_states)+shift-1

            num_connections=num_states*num_adj_states
            ! Change det_j_up and det_j_dn at chosen block location
            do i=1,4
                if (btest(new_vector,i-1)) then
                    det_j_dn=ibset(det_j_dn,4*(chosen_block-1)+i-1)
                else
                    det_j_dn=ibclr(det_j_dn,4*(chosen_block-1)+i-1)
                endif
            enddo
            do i=5,8
                if (btest(new_vector,i-1)) then
                    det_j_up=ibset(det_j_up,4*(chosen_block-1)+i-4-1)
                else
                    det_j_up=ibclr(det_j_up,4*(chosen_block-1)+i-4-1)
                endif
            enddo

            ! Change det_j_up and det_j_dn at nbr block location
            do i=1,4
                if (btest(new_adj_vector,i-1)) then
                    det_j_dn=ibset(det_j_dn,4*(nbr_block-1)+i-1)
                else
                    det_j_dn=ibclr(det_j_dn,4*(nbr_block-1)+i-1)
                endif
            enddo
            do i=5,8
                if (btest(new_adj_vector,i-1)) then
                    det_j_up=ibset(det_j_up,4*(nbr_block-1)+i-4-1)
                else
                    det_j_up=ibclr(det_j_up,4*(nbr_block-1)+i-4-1)
                endif
            enddo
            success=.true.
            call check_nbr(l_x/2,l_y/2,pbc,chosen_block,nbr_block,nbr_types,num_nbr_types)
            extra_factor=2._rk/real(num_nbr_types) ! 2 for spin x 2 for electron lost or gained, 2 for the number of ways of getting to this state
        else
            success=.false.
        endif

    endif

    if (success .eqv. .true.) then
         call hamiltonian_hubbard_dm(det_i_up,det_i_dn,det_j_up,det_j_dn,matrix_element)
         if (spr) write (6,*) "matrix_element (off diag move) =",matrix_element
         if (spr) call flush(6)
         proposal_prob_inv=(real(num_blocks*5*num_connections*extra_factor))
         if (spr) write(6,*) "proposal_prob_inv",proposal_prob_inv
         acceptance_prob=abs(matrix_element)*tau*proposal_prob_inv
         if (spr) write(6,*) "acceptance_prob",acceptance_prob
         call flush(6)
         weight_j = acceptance_prob*sign(1._rk,-matrix_element)
   endif

  !if (spr) stop
  end subroutine off_diagonal_move_hubbard_dm
  !===========================================================================================================

!======================================================================================================================
  subroutine off_diagonal_move_hubbard_dm2(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j, diagonal_matrix_element)
!======================================================================================================================
!    !------------------------------------------------------------------------------------------------------------------
!    ! Description : Perform an offdiagonal uniform move in SQMC for Hubbard model in a unitary transformaed basis
!    !               starting at det_i which is split into det_i_up and det_i_dn.
!    !		        det_i has weight = +-1. If move fails, for any reason,
!    !               weight_j is returned as 0. Otherwise, return new determinant,
!    !               det_j, and weight of incoming det, H_ii (diagonal element of Hamiltonian).
!    ! Created     : H.J.Changlani, Dec 06 2011
!    !------------------------------------------------------------------------------------------------------------------
    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up,det_i_dn
    type(ik_vec),intent (out):: det_j_up,det_j_dn
#else
    integer(ik), intent(in) :: det_i_up,det_i_dn
    integer(ik),intent (out):: det_j_up,det_j_dn
#endif
!   integer, intent(out)    :: weight_j
    real(rk), intent(out)   :: weight_j
    real(rk), optional,intent(out)   :: diagonal_matrix_element

    !local variables
!    integer                  ::  i,j,counter
    integer                  ::  i,counter
    real(rk)                 ::  out_weight!, out_weight_adj
    real(rk)                 ::  proposal_prob_inv, acceptance_prob,matrix_element
    real(rk)                 ::  rannyu
    integer                  ::  chosen_block,nbr_block,tmp_block,chosen_move
!    integer                  ::  chosen_electron,chosen_block,nbr_block,chosen_move
    integer                  ::  nup_block,ndn_block
    integer                  ::  nup_adj_block,ndn_adj_block,tmp_n_block
    integer                  ::  num_states,num_adj_states,num_connections
    integer (ik)             ::  chosen_vector,chosen_vector_adj,new_vector,new_adj_vector,tmp_vector,shift,shift_adj
!    integer                  ::  n_tot
    integer                  ::  lost,nup_new_block,ndn_new_block,nup_new_adj_block,ndn_new_adj_block
    real(rk)                 ::  extra_factor
    logical                  ::  e_lost,up,success,allowed
    integer                  ::  nbr_types(4),num_nbr_types
    integer(ik)              ::  list_vector(36)
    real(rk)                 ::  weight(36),weight_table(36,36)
    logical                  ::  spr
!   real(rk)                 ::  min_dm_wt=0.01,beta=0.01
    integer                  ::  e_lost_int,updn_int

!   spr=.true.
    spr=.false.

    success=.false.
    up=.false.
    e_lost=.false.
    weight_j=0._rk

    !initialize det_j to det_i
    det_j_up = det_i_up
    det_j_dn = det_i_dn
    if (present(diagonal_matrix_element)) then
      diagonal_matrix_element=0._rk
      call hamiltonian_hubbard_dm(det_i_up,det_i_dn,det_j_up,det_j_dn,diagonal_matrix_element)
    endif
    !write (6,*) "diagonal element=",diagonal_matrix_element
    call flush(6)
    ! Choose electron at random and gets its block and type (up vs down)
    !chosen_electron=random_int(nup+ndn)
    chosen_block=random_int(nsites/4)

   ! Once block is chosen choose a number from 1 to 5
   if (use_eig) then
    chosen_move=random_int(4)+1          !guarantees number from 2-5 as no need for on patch move
   else
    chosen_move=random_int(5)
   endif

   ! 1 = On patch move
   ! 2 = Left move
   ! 3 = Right move
   ! 4 = Up move
   ! 5 = Down move

    chosen_vector=0_ik
    do i=1,4
        if (btest(det_i_up,4*(chosen_block-1)+i-1)) chosen_vector=ibset(chosen_vector,4+i-1)
        if (btest(det_i_dn,4*(chosen_block-1)+i-1)) chosen_vector=ibset(chosen_vector,i-1)
    enddo

    nup_block=nup_dm(chosen_vector)
    ndn_block=ndn_dm(chosen_vector)

    if (chosen_move .eq. 1) then                                        ! On patch move
        weight(:)=0._rk
        num_states=ndm_start_end(nup_block,ndn_block,2)-ndm_start_end(nup_block,ndn_block,1)+1
        shift=ndm_start_end(nup_block,ndn_block,1)
        num_connections=num_states-1 ! Number of off diagonal connections on patch
        if (num_states .eq. 1) return
        counter=0
        do i=1,num_states
            if (i+shift-1.ne. chosen_vector) then
                counter=counter+1
                list_vector(counter)=i+shift-1
                !weight(counter)=max(min_dm_wt,dm_eigenvalues(i+shift-1))
                weight(counter)=abs(ham_on_2x2(list_vector(counter),chosen_vector))
            endif
        enddo
        ! Choose the new state in same number,spin sector weighed by the dm eigenvalues. Do not choose itself!
        call choose_entry_by_its_weight(list_vector,weight(1:counter),counter,new_vector,out_weight,success)
        if (success .eqv. .false.) return
        matrix_element=ham_on_2x2(new_vector,chosen_vector)

        ! Change det_j_up and det_j_dn
        do i=1,4
            if (btest(new_vector,i-1)) then
                det_j_dn=ibset(det_j_dn,4*(chosen_block-1)+i-1)
            else
                det_j_dn=ibclr(det_j_dn,4*(chosen_block-1)+i-1)
            endif
        enddo
        do i=5,8
            if (btest(new_vector,i-1)) then
                det_j_up=ibset(det_j_up,4*(chosen_block-1)+i-4-1)
            else
                det_j_up=ibclr(det_j_up,4*(chosen_block-1)+i-4-1)
            endif
        enddo
        success=.true.
        !write (6,*) "out_weight (on patch)=",out_weight
        extra_factor=5._rk/out_weight
    else                                                                ! Off patch moves in 4 directions
        call get_nbr(l_x/2,l_y/2,pbc,chosen_block,chosen_move-1,nbr_block,allowed)

        if (allowed) then
            chosen_vector_adj=0_ik
            do i=1,4
                if (btest(det_i_up,4*(nbr_block-1)+i-1)) chosen_vector_adj=ibset(chosen_vector_adj,4+i-1)
                if (btest(det_i_dn,4*(nbr_block-1)+i-1)) chosen_vector_adj=ibset(chosen_vector_adj,i-1)
            enddo

            nup_adj_block=nup_dm(chosen_vector_adj)
            ndn_adj_block=ndn_dm(chosen_vector_adj)

            ! Because of the way the matrix tables are built we require that block 1 has a lower numbering than its neighbor
            ! Hence if this is the case one needs to swap chosen and nbr_block and associated quantities
            if (chosen_block>nbr_block) then
                tmp_block=chosen_block
                chosen_block=nbr_block
                nbr_block=tmp_block
                tmp_vector=chosen_vector
                chosen_vector=chosen_vector_adj
                chosen_vector_adj=tmp_vector
                tmp_n_block=nup_block
                nup_block=nup_adj_block
                nup_adj_block=tmp_n_block
                tmp_n_block=ndn_block
                ndn_block=ndn_adj_block
                ndn_adj_block=tmp_n_block
            endif

            if (rannyu() >0.5_rk) then
                e_lost=.true.
                lost=-1
                e_lost_int=1
            else
                e_lost=.false.
                lost=1
                e_lost_int=2
            endif

            if (rannyu() >0.5_rk) then
                up=.true.
                updn_int=1
            else
                up=.false.
                updn_int=2
            endif

            if (up) then
                nup_new_block=nup_block+lost
                ndn_new_block=ndn_block
                nup_new_adj_block=nup_adj_block-lost
                ndn_new_adj_block=ndn_adj_block
            else
                nup_new_block=nup_block
                nup_new_adj_block=nup_adj_block
                ndn_new_block=ndn_block+lost
                ndn_new_adj_block=ndn_adj_block-lost
            endif

            if ((nup_new_block .lt. 0) .or. (nup_new_block .gt. 4) .or. (ndn_new_block .lt. 0) .or. (ndn_new_block .gt. 4) .or. (nup_new_adj_block .lt. 0) .or. (nup_new_adj_block .gt. 4) .or. (ndn_new_adj_block .lt. 0) .or. (ndn_new_adj_block .gt. 4))then
                success=.false.
                return
            endif

            num_states=ndm_start_end(nup_new_block,ndn_new_block,2)-ndm_start_end(nup_new_block,ndn_new_block,1)+1
            num_adj_states=ndm_start_end(nup_new_adj_block,ndn_new_adj_block,2)-ndm_start_end(nup_new_adj_block,ndn_new_adj_block,1)+1
            ! Choose a state randomly from these
            shift=ndm_start_end(nup_new_block,ndn_new_block,1)
            shift_adj=ndm_start_end(nup_new_adj_block,ndn_new_adj_block,1)
            call check_nbr(l_x/2,l_y/2,pbc,chosen_block,nbr_block,nbr_types,num_nbr_types)
            !counter=0
            !do i=1,num_states
            !        counter=counter+1
            !        list_vector(counter)=i+shift-1
            !        if (use_dm)  weight(counter)=max(min_dm_wt,dm_eigenvalues(i+shift-1))
            !        if (use_eig) weight(counter)=exp(-beta*dm_eigenvalues(i+shift-1))
            !enddo
            weight_table(1:num_states,1:num_adj_states)=0._rk
            do i=1,num_nbr_types
                weight_table(1:num_states,1:num_adj_states)=weight_table(1:num_states,1:num_adj_states)+table_block_matel(updn_int,e_lost_int,nbr_types(i),chosen_vector,chosen_vector_adj,1:num_states,1:num_adj_states)
            enddo
            call choose_entry_by_its_weight_rank2(weight_table(1:num_states,1:num_adj_states),num_states,num_adj_states,new_vector,new_adj_vector,out_weight,success)
            if (success .eqv. .false.) return

            !matrix_element=weight_table(new_vector,new_adj_vector)
            new_vector=new_vector+shift-1
            new_adj_vector=new_adj_vector+shift_adj-1
            !n_tot=0
            !do i=chosen_block+1,nbr_block-1
            !    ket=0_ik
            !    do j=1,4
            !        if (btest(det_i_up,4*(i-1)+j-1)) then
            !            ket=ibset(ket,4+j-1)
            !        endif
            !        if (btest(det_i_dn,4*(i-1)+j-1)) then
            !            ket=ibset(ket,j-1)
            !        endif
            !    enddo
            !    if (up) then
            !        n_tot=n_tot+nup_dm(ket)
            !    else
            !        n_tot=n_tot+ndn_dm(ket)
            !    endif
            !enddo

            !if (mod(n_tot,2) .ne. 0) matrix_element=-matrix_element
            !if (spr) write(6,*) "matrix_element =",matrix_element
            ! Choose a state randomly from these
            !shift=ndm_start_end(nup_new_adj_block,ndn_new_adj_block,1)
            !counter=0
            !do i=1,num_adj_states
            !        counter=counter+1
            !        list_vector(counter)=i+shift-1
            !        if (use_dm)  weight(counter)=max(min_dm_wt,dm_eigenvalues(i+shift-1))
            !        if (use_eig) weight(counter)=exp(-beta*dm_eigenvalues(i+shift-1))
            !enddo
            !call choose_entry_by_its_weight(list_vector,weight,counter,new_adj_vector,out_weight_adj)

            ! Change det_j_up and det_j_dn at chosen block location
            do i=1,4
                if (btest(new_vector,i-1)) then
                    det_j_dn=ibset(det_j_dn,4*(chosen_block-1)+i-1)
                else
                    det_j_dn=ibclr(det_j_dn,4*(chosen_block-1)+i-1)
                endif
            enddo
            do i=5,8
                if (btest(new_vector,i-1)) then
                    det_j_up=ibset(det_j_up,4*(chosen_block-1)+i-4-1)
                else
                    det_j_up=ibclr(det_j_up,4*(chosen_block-1)+i-4-1)
                endif
            enddo

            ! Change det_j_up and det_j_dn at nbr block location
            do i=1,4
                if (btest(new_adj_vector,i-1)) then
                    det_j_dn=ibset(det_j_dn,4*(nbr_block-1)+i-1)
                else
                    det_j_dn=ibclr(det_j_dn,4*(nbr_block-1)+i-1)
                endif
            enddo
            do i=5,8
                if (btest(new_adj_vector,i-1)) then
                    det_j_up=ibset(det_j_up,4*(nbr_block-1)+i-4-1)
                else
                    det_j_up=ibclr(det_j_up,4*(nbr_block-1)+i-4-1)
                endif
            enddo
            success=.true.
            if (use_dm)  extra_factor=10._rk/(out_weight*real(num_nbr_types))
            if (use_eig) extra_factor=8._rk/(out_weight*real(num_nbr_types))
        else
            success=.false.
            return
        endif

    endif

    if (success .eqv. .true.) then
         call hamiltonian_hubbard_dm(det_i_up,det_i_dn,det_j_up,det_j_dn,matrix_element)
         if (spr) write (6,*) "matrix_element (off diag move) =",matrix_element
         if (spr) call flush(6)
         proposal_prob_inv=(real((nsites/4)*extra_factor))
         if (spr) write(6,*) "proposal_prob_inv",proposal_prob_inv
         acceptance_prob=abs(matrix_element)*tau*proposal_prob_inv
         if (spr) write(6,*) "acceptance_prob",acceptance_prob
         call flush(6)
!        weight_j = int((acceptance_prob+rannyu())*sign(1._rk,-matrix_element)) ! For integer weights
         weight_j = acceptance_prob*sign(1._rk,-matrix_element)
   endif
  !if (spr) stop

  end subroutine off_diagonal_move_hubbard_dm2
  !===========================================================================================================

!======================================================================================================================================
  subroutine off_diagonal_move_hubbard_k(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j, importance_sampling,proposal_prob_inverse)
!======================================================================================================================================
    !------------------------------------------------------------------------------------------------------------------
    ! Description : Perform an offdiagonal uniform move in SQMC for Hubbard model in k space
    !               starting at det_i which is split into det_i_up and det_i_dn which are k space occupation numbers
    !		        det_i has weight = +-1. If move fails, for any reason,
    !               weight_j is returned as 0. Otherwise, return new determinant,
    !               det_j, and weight of incoming det, H_ii (diagonal element of Hamiltonian).

    ! Note        : This is the preferred move for k space
    !               (over heat bath algorithm, because the number of connected dets is larger in k space)
    ! Created     : H.J.Changlani, July 10 2011
    ! Modified    : A Holmes, 31 May 2012. Sped up.
    !------------------------------------------------------------------------------------------------------------------
    use tools, only : permutation_factor
    use more_tools, only : binary_search
    use common_psi_t, only : psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_g,psi_g_epsilon,psi_g_epsilon_inv

    implicit none

    !dummy variables
    real(rk), intent(out)   :: weight_j
    integer, optional       :: importance_sampling
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up,det_i_dn
    type(ik_vec),intent (out):: det_j_up,det_j_dn
    type(ik_vec)              ::  det
    type(ik_vec)              ::  new_dets_i_up(16),new_dets_i_dn(16),new_dets_j_up(16),new_dets_j_dn(16)
    type(ik_vec)              ::  rep_i_up,rep_i_dn,rep_j_up,rep_j_dn
    type(ik_vec)              ::  tmp_distinct_up(16),tmp_distinct_dn(16)
#else
    integer(ik), intent(in) :: det_i_up,det_i_dn
    integer(ik),intent (out):: det_j_up,det_j_dn
    integer(ik)              ::  det
    integer(ik)              ::  new_dets_i_up(16),new_dets_i_dn(16),new_dets_j_up(16),new_dets_j_dn(16)
    integer(ik)              ::  rep_i_up,rep_i_dn,rep_j_up,rep_j_dn
    integer(ik)              ::  tmp_distinct_up(16),tmp_distinct_dn(16)
#endif
    real(rk),optional,intent(out) :: proposal_prob_inverse

    !local variables
    real(rk)                 ::  proposal_prob_inv, acceptance_prob,guiding_wf_ratio
    real(rk)                 ::  matrix_element,tmp_matrix_element
!   real(rk)                 ::  rannyu
    integer                  ::  i,counter,location
    integer                  ::  chosen_up,chosen_dn,target_up,target_dn
    integer                  ::  q_x,q_y
    integer                  ::  counter_up,counter_dn!,counter_empty_up
    integer                  ::  up_locations(nup),dn_locations(ndn)!,up_empty_locations(nsites-nup)
    logical                  ::  imp,entry_present
    integer                  ::  nnzero
    real(rk)                 ::  norm_i,norm_j,phase_w_rep_i,phase_w_rep_j,phases_i(16),phases_j(16)!,tmp_phases(16)
    integer                  ::  num_distinct_i,num_distinct_j
    integer                  :: n

    imp=.false.
    if (present(importance_sampling)) then
        if (importance_sampling .eq. 1) imp=.true.
    endif

    counter_up=1
    counter_dn=1

    !initialize det_j to det_i
    det_j_up = det_i_up
    det_j_dn = det_i_dn

    ! get index numbers of occupied orbitals
    det = det_i_up
    do while (det.ne.0)
      i = trailz(det)+1
      up_locations(counter_up)=i
      counter_up=counter_up+1
      det = ibclr(det,i-1)
    enddo
    det = det_i_dn
    do while (det.ne.0)
      i = trailz(det)+1
      dn_locations(counter_dn)=i
      counter_dn=counter_dn+1
      det = ibclr(det,i-1)
    enddo

    if (space_sym) then
        call generate_fourfold_k_configs_efficient_given_locations(c4_map,reflection_map,z,p,nup,ndn,det_i_up,det_i_dn,up_locations,dn_locations,new_dets_i_up,new_dets_i_dn,phases_i,rep_i_up,rep_i_dn,phase_w_rep_i,norm_i,num_distinct_i)
    endif

    weight_j=0._rk

    chosen_up = random_int(nup)
    chosen_up = up_locations(chosen_up)

    chosen_dn = random_int(ndn)
    chosen_dn = dn_locations(chosen_dn)

    target_up=random_int(nsites-nup)
    n=0
    do i=1,nsites
        if (btest(det_i_up,i-1).eqv..false.) then
            n=n+1
            if (n==target_up) then
                target_up=i
                exit
            endif
        endif
    enddo

    q_x= -k_vectors(1,chosen_up) + k_vectors(1,target_up) ! q = change in momentum of up det
    q_y= -k_vectors(2,chosen_up) + k_vectors(2,target_up)

    target_dn=kmap(chosen_dn,int(-q_x/2)+l_x,int(-q_y/2)+l_y)

    if (btest(det_i_dn,target_dn-1) .eqv. .false.) then
         det_j_up=ibclr(det_i_up,chosen_up-1)
         det_j_up=ibset(det_j_up,target_up-1)
         det_j_dn=ibclr(det_i_dn,chosen_dn-1)
         det_j_dn=ibset(det_j_dn,target_dn-1)

         proposal_prob_inv=(real((nup)*(ndn)*(nsites-nup)))

         if (present(proposal_prob_inverse))  proposal_prob_inverse = proposal_prob_inv

         if (space_sym) then
            matrix_element=0._rk
            counter=0
            nnzero=0
            tmp_distinct_up(:)=0_ik
            tmp_distinct_dn(:)=0_ik
            !call hamiltonian_hubbard_k_space_sym(det_j_up,det_j_dn,det_i_up,det_i_dn,diagonal_matrix_element,nnzero)
            call generate_fourfold_k_configs_efficient(c4_map,reflection_map,z,p,nup,ndn,det_j_up,det_j_dn,new_dets_j_up,new_dets_j_dn,phases_j,rep_j_up,rep_j_dn,phase_w_rep_j,norm_j,num_distinct_j)
            if (abs(norm_j) .gt. 1.0d-10) then
                if (((rep_j_up .eq. det_i_up) .and. (rep_j_dn .eq. det_i_dn))) return
                !do i=1,num_distinct_i
                do i=1,16
                    if (counter .eq. num_distinct_j) exit
                    call is_in_list_mod(new_dets_j_up(i),new_dets_j_dn(i),tmp_distinct_up,tmp_distinct_dn,counter,entry_present,location)
                    if (entry_present .eqv. .false.) then
                        counter=counter+1
                        tmp_distinct_up(counter)=new_dets_j_up(i)
                        tmp_distinct_dn(counter)=new_dets_j_dn(i)
                        call hamiltonian_hubbard_k(new_dets_j_up(i),new_dets_j_dn(i),det_i_up,det_i_dn,tmp_matrix_element)
                        if (abs(tmp_matrix_element) .gt. 1.0d-10) then
                            matrix_element=matrix_element+(tmp_matrix_element*phases_j(i))  ! Refer to Lauechli/HJC notes,SD-Frick paper
                            nnzero=nnzero+1
                        endif
                    endif
                    !call hamiltonian_hubbard_k(tmp_distinct_up(i),tmp_distinct_dn(i),det_j_up,det_j_dn,tmp_matrix_element)
                    !if (abs(tmp_matrix_element) .gt. 1.0d-10) then
                    !    matrix_element=matrix_element+(tmp_matrix_element*tmp_phases(i))  ! Refer to Lauechli/HJC notes,SD-Frick paper
                    !    nnzero=nnzero+1
                    !endif
                enddo
                det_j_up=rep_j_up
                det_j_dn=rep_j_dn
                if (nnzero .gt. 0) then
                    !matrix_element=matrix_element*phase_w_rep_j*(norm_i/norm_j)   ! Refer to Lauechli/HJC notes,SD-Frick paper
                    matrix_element=matrix_element*phase_w_rep_j*(norm_j/norm_i)  ! Order of i and j depends on which is bra and which is ket
                    !call hamiltonian_hubbard_k_space_sym(det_j_up,det_j_dn,det_i_up,det_i_dn,matrix_element,nnzero) ! det_j_up and det_j_dn will be replaced by its representative
                    !if (((det_j_up .eq. det_i_up) .and. (det_j_dn .eq. det_i_dn)) .or. (nnzero .eq. 0)) return
                    proposal_prob_inv=proposal_prob_inv/nnzero
                    if (present(proposal_prob_inverse))  proposal_prob_inverse = proposal_prob_inv
                else
                    return
                endif
            else
                det_j_up=rep_j_up
                det_j_dn=rep_j_dn
                return
            endif
         else
            ! Edited by AAH on 30 May 2012
            matrix_element = permutation_factor(det_i_up,det_j_up)*permutation_factor(det_i_dn,det_j_dn)*ubyn
            ! End of edit

            if (imp) then

              call binary_search(det_i_up,det_i_dn,psi_t_connected_dets_up,psi_t_connected_dets_dn,location)
              if (location>0) then
                guiding_wf_ratio=1._rk/psi_g(location)
              else
                guiding_wf_ratio=psi_g_epsilon_inv
              endif

              call binary_search(det_j_up,det_j_dn,psi_t_connected_dets_up,psi_t_connected_dets_dn,location)
              if (location>0) then
                matrix_element=matrix_element*guiding_wf_ratio*psi_g(location)
              else
                matrix_element=matrix_element*guiding_wf_ratio*psi_g_epsilon
              endif

            endif

         endif

         acceptance_prob=abs(matrix_element)*tau*proposal_prob_inv
!        weight_j = int((acceptance_prob+rannyu())*sign(1._rk,-matrix_element))
         weight_j = acceptance_prob*sign(1._rk,-matrix_element)
    endif

  end subroutine off_diagonal_move_hubbard_k
  !===========================================================================================================

  !===========================================================================================================
   subroutine wf_calc (det_up,det_dn,wf_value,full_wf)
  ! H Changlani, 2010
  !===========================================================================================================

#ifdef NUM_ORBITALS_GT_127
   type(ik_vec),intent(in)              :: det_up,det_dn
#else
   integer(ik),intent(in)              :: det_up,det_dn
#endif
   real(rk),intent(out)                :: wf_value
   integer                             :: i,k,n
   integer                             :: double_occ
   real(rk)                            :: up_detval,dn_detval
   complex*16                          :: c_up_detval,c_dn_detval
   real(rk)                            :: up_inv(nup,nup),dn_inv(ndn,ndn)
   complex*16                          :: c_up_inv(nup,nup),c_dn_inv(ndn,ndn)
   logical,optional                    :: full_wf
   logical                             :: full
   complex*16                          :: c_wf_value

   if (present(full_wf)) then
    full=full_wf
   else
    full=.true.
   endif

   k=0
   double_occ=0
   wf_value=0._rk
   c_wf_value=dcmplx(0._rk,0._rk)

   if (wf_type .eq. 'gutz_rhf') then
       do i=1,nsites
            if (btest(det_up,i-1)) then
                 k=k+1
                 up_inv(k,1:nup)=rhf_up_orbitals(i,1:nup) ! not yet the inverse until determinant function is called
            endif
       enddo
       k=0
       do i=1,nsites
        if (btest(det_dn,i-1)) then
             k=k+1
             dn_inv(k,1:ndn)=rhf_dn_orbitals(i,1:ndn)
        endif
       enddo
       call determinant(up_inv,nup,up_detval)
       call determinant(dn_inv,ndn,dn_detval)
       wf_value=up_detval*dn_detval
   endif

   if (wf_type .eq. 'gutz_uhf') then
       do i=1,nsites
            if (btest(det_up,i-1)) then
                 k=k+1
                 up_inv(k,1:nup)=uhf_up_orbitals(i,1:nup) ! not yet the inverse until determinant function is called
            endif
       enddo
       k=0
       do i=1,nsites
        if (btest(det_dn,i-1)) then
             k=k+1
             dn_inv(k,1:ndn)=uhf_dn_orbitals(i,1:ndn)
        endif
       enddo
       call determinant(up_inv,nup,up_detval)
       call determinant(dn_inv,ndn,dn_detval)
       wf_value=up_detval*dn_detval
   endif

   if (wf_type .eq. 'gutz_multi') then
       do n=1,nmultidet                                     ! Loop over number of determinants in multi determinant wavefunction
           if (coeffs(n) .ne. 0._rk) then
               k=0
               do i=1,nsites
                    if (btest(det_up,i-1)) then
                         k=k+1
                         up_inv(k,1:nup)=up_orbitals(n,i,1:nup) ! not yet the inverse until determinant function is called
                    endif
               enddo
               k=0
               do i=1,nsites
                if (btest(det_dn,i-1)) then
                     k=k+1
                     dn_inv(k,1:ndn)=dn_orbitals(n,i,1:ndn)
                endif
               enddo
               call determinant(up_inv,nup,up_detval)
               call determinant(dn_inv,ndn,dn_detval)
               wf_value=wf_value+(coeffs(n)*up_detval*dn_detval)
           endif
       enddo
   endif

   if (wf_type .eq. 'cgutz_multi') then
       do n=1,nmultidet                                     ! Loop over number of determinants in multi determinant wavefunction
           if (coeffs(n) .ne. 0._rk) then
               k=0
               do i=1,nsites
                    if (btest(det_up,i-1)) then
                         k=k+1
                         c_up_inv(k,1:nup)=c_up_orbitals(n,i,1:nup) ! not yet the inverse until determinant function is called
                    endif
               enddo
               k=0
               do i=1,nsites
                if (btest(det_dn,i-1)) then
                     k=k+1
                     c_dn_inv(k,1:ndn)=c_dn_orbitals(n,i,1:ndn)
                endif
               enddo
               !call print_complex_matrix(nup,nup,c_up_inv)
               !call print_complex_matrix(ndn,ndn,c_dn_inv)
               call cdeterminant(c_up_inv,nup,c_up_detval)
               call cdeterminant(c_dn_inv,ndn,c_dn_detval)
               c_wf_value=c_wf_value+(coeffs(n)*c_up_detval*c_dn_detval)
           endif
       enddo
       wf_value=real(c_wf_value)
       !write(6,*) "wf_value=", wf_value
       !write(6,*) "c_wf_value=", c_wf_value
       !call flush(6)
   endif


   do i=1,l_x*l_y
        if (btest(det_up,i-1) .and. btest(det_dn,i-1)) double_occ=double_occ+1
   enddo

   if (full) then
    wf_value=(g**real(double_occ))*wf_value
   endif

  end subroutine wf_calc

  !===========================================================================================================
  subroutine get_connected_dets_and_info_hubbard(det_up,det_dn,chosen_e_row,chosen_site,target_site,    &
                                                 & e_types,up_inv,dn_inv,up_detval,dn_detval,wf_detval, &
                                                 & double_occ,diff_double_occs)
  !-----------------------------------------------------------------------------------------------------------
  ! Description : Gets connected determinants and associated information
  ! Created     : Hitesh J. Changlani, Dec. 2010
  !-----------------------------------------------------------------------------------------------------------

   implicit none

   integer                             :: i,j,k,n
#ifdef NUM_ORBITALS_GT_127
   type(ik_vec),intent(in)              :: det_up,det_dn
#else
   integer(ik),intent(in)              :: det_up,det_dn
#endif
   integer,intent(out)                 :: chosen_e_row(4*(nup+ndn)),chosen_site(4*(nup+ndn)),target_site(4*(nup+ndn))
   integer,intent(out)                 :: e_types(4*(nup+ndn)),diff_double_occs(4*(nup+ndn))
   integer,intent(out)                 :: double_occ
   real(rk),intent(out)                :: up_inv(nup,nup),dn_inv(ndn,ndn)
   real(rk),intent(out)                :: up_detval,dn_detval,wf_detval

   ! local variables
   integer                            ::  ind
   integer                            ::  num_e,diff_double_occs_tmp
   integer                            ::  nbr_sites(nup+ndn,1+4)
   logical                            ::  allowed
   complex*16                         ::  c_up_detval,c_dn_detval,c_wf_detval
   complex*16                         ::  c_up_inv(nup,nup),c_dn_inv(ndn,ndn)

   wf_detval=0._rk

   k=0
   if (wf_type .eq. 'gutz_rhf') then
       do i=1,nsites
            if (btest(det_up,i-1)) then
                 k=k+1
                 up_inv(k,1:nup)=rhf_up_orbitals(i,1:nup)
            endif
       enddo

       k=0
       do i=1,nsites
        if (btest(det_dn,i-1)) then
             k=k+1
             dn_inv(k,1:ndn)=rhf_dn_orbitals(i,1:ndn)
        endif
       enddo
       call determinant(up_inv,nup,up_detval)
       call determinant(dn_inv,ndn,dn_detval)
       wf_detval=up_detval*dn_detval
   endif

   if (wf_type .eq. 'gutz_uhf') then
       do i=1,nsites
            if (btest(det_up,i-1)) then
                 k=k+1
                 up_inv(k,1:nup)=uhf_up_orbitals(i,1:nup)
            endif
       enddo
       k=0
       do i=1,nsites
        if (btest(det_dn,i-1)) then
             k=k+1
             dn_inv(k,1:ndn)=uhf_dn_orbitals(i,1:ndn)
        endif
       enddo
       call determinant(up_inv,nup,up_detval)
       call determinant(dn_inv,ndn,dn_detval)
       wf_detval=up_detval*dn_detval
   endif

   if (wf_type .eq. 'gutz_multi') then
       do n=1,nmultidet                                     ! Loop over number of determinants in multi determinant wavefunction
           k=0
           do i=1,nsites
                if (btest(det_up,i-1)) then
                     k=k+1
                     up_inv(k,1:nup)=up_orbitals(n,i,1:nup) ! not yet the inverse until determinant function is called
                endif
           enddo
           k=0
           do i=1,nsites
            if (btest(det_dn,i-1)) then
                 k=k+1
                 dn_inv(k,1:ndn)=dn_orbitals(n,i,1:ndn)
            endif
           enddo
           call determinant(up_inv,nup,up_detval)
           call determinant(dn_inv,ndn,dn_detval)
           wf_detval=wf_detval+(coeffs(n)*up_detval*dn_detval)
       enddo
   endif

  if (wf_type .eq. 'cgutz_multi') then
       do n=1,nmultidet                                     ! Loop over number of determinants in multi determinant wavefunction
           k=0
           do i=1,nsites
                if (btest(det_up,i-1)) then
                     k=k+1
                     c_up_inv(k,1:nup)=c_up_orbitals(n,i,1:nup) ! not yet the inverse until determinant function is called
                endif
           enddo
           k=0
           do i=1,nsites
            if (btest(det_dn,i-1)) then
                 k=k+1
                 c_dn_inv(k,1:ndn)=c_dn_orbitals(n,i,1:ndn)
            endif
           enddo
           call cdeterminant(c_up_inv,nup,c_up_detval)
           call cdeterminant(c_dn_inv,ndn,c_dn_detval)
           c_wf_detval=c_wf_detval+(coeffs(n)*c_up_detval*c_dn_detval)
       enddo
       wf_detval=real(c_wf_detval)
   endif

   num_e=0
   double_occ=0

   do i=1,l_x*l_y
        if (num_e==nup) exit

        if (btest(det_up,i-1)) then

            num_e=num_e+1
            nbr_sites(num_e,1)=i

            if (btest(det_dn,i-1)) then
              double_occ=double_occ+1
              diff_double_occs_tmp=1
            else
              diff_double_occs_tmp=0
            endif

            do j=1,4                    ! Left Right Up Down
                call get_nbr(l_x,l_y,pbc,i,j,nbr_sites(num_e,j+1),allowed)
                ind=(4*(num_e-1))+j
                diff_double_occs(ind)=diff_double_occs_tmp

                !if ((allowed .eqv. .true.) .and. ((nbr_sites(num_e,j+1) .ne. nbr_sites(num_e,j)))) then      ! I HAD put this for length 2 systems but its not consistent with kspace
                if ((allowed .eqv. .true.) ) then                                                             ! Wrapping twice allowed in real space
                        if (btest(det_up,nbr_sites(num_e,j+1)-1) .eqv. .false.) then
                                chosen_site(ind)=i
                                chosen_e_row(ind)=num_e
                                target_site(ind)=nbr_sites(num_e,j+1)
                                e_types(ind)=1
                                if (btest(det_dn,target_site(ind)-1) .eqv. .true. ) then
                                    diff_double_occs(ind)=1-diff_double_occs(ind)
                                else
                                    diff_double_occs(ind)=0-diff_double_occs(ind)
                                endif
                        else
                                e_types(ind)=0
                        endif
                 else
                        e_types(ind)=0
                 endif
            enddo
        endif
    enddo

    num_e=0
    do i=1,l_x*l_y
        if (num_e==ndn) exit

        if (btest(det_dn,i-1)) then
            num_e=num_e+1
            nbr_sites(num_e+nup,1)=i

            if (btest(det_up,i-1)) then
              diff_double_occs_tmp=1
            else
              diff_double_occs_tmp=0
            endif

            do j=1,4                    ! Left Right Up Down
                call get_nbr(l_x,l_y,pbc,i,j,nbr_sites(num_e+nup,j+1),allowed)
                ind=(4*(num_e+nup-1))+j
                diff_double_occs(ind)=diff_double_occs_tmp

                !if ((allowed .eqv. .true.) .and. ((nbr_sites(num_e,j+1) .ne. nbr_sites(num_e,j)))) then      ! I HAD put this for length 2 systems but its not consistent with kspace
                if (allowed .eqv. .true.) then
                        if (btest(det_dn,nbr_sites(num_e+nup,j+1)-1) .eqv. .false.) then
                                chosen_site(ind)=i
                                chosen_e_row(ind)=num_e
                                target_site(ind)=nbr_sites(num_e+nup,j+1)
                                e_types(ind)=-1
                                if (btest(det_up,target_site(ind)-1) .eqv. .true.) then
                                    diff_double_occs(ind)=1-diff_double_occs(ind)
                                else
                                    diff_double_occs(ind)=0-diff_double_occs(ind)
                                endif
                        else
                                e_types(ind)=0
                        endif
                 else
                        e_types(ind)=0
                 endif
            enddo
        endif
    enddo

    end subroutine get_connected_dets_and_info_hubbard
  !=================================================================================================


  subroutine explicit_ratio_hubbard(e_type, n_e,det_in,chosen_e_row,target_site,ratio)
  !-----------------------------------------------------------------------------------------------------------
  ! Description : Computes the ratio of determinants explicitly. This mostly serves as a check for SMW
  ! Created     : Hitesh J. Changlani, Dec. 2010
  !-----------------------------------------------------------------------------------------------------------

   implicit none

   ! Dummy
   integer, intent(in)   :: e_type
   integer, intent(in)   :: n_e
   integer,intent(in)    :: chosen_e_row,target_site
   real(rk),intent(out)  :: ratio

   ! Local
   integer               :: i,k
   integer (ik), intent(in) :: det_in
   real(rk)              :: mat1(n_e,n_e),mat2(n_e,n_e)
   real(rk)              :: det1,det2

   k=0
   if (wf_type .eq. 'gutz_rhf') then
       if (e_type .eq. 1) then
          do i=1,nsites
                if (btest(det_in,i-1)) then
                     k=k+1
                     mat1(k,1:nup)=rhf_up_orbitals(i,1:nup)
                endif
          enddo
       endif

       if (e_type .eq. -1) then
          do i=1,nsites
                if (btest(det_in,i-1)) then
                     k=k+1
                     mat1(k,1:ndn)=rhf_dn_orbitals(i,1:ndn)
                endif
          enddo
       endif

      mat2=mat1
      if (e_type .eq.  1) mat2(chosen_e_row,1:nup)=rhf_up_orbitals(target_site,1:nup)
      if (e_type .eq. -1) mat2(chosen_e_row,1:ndn)=rhf_dn_orbitals(target_site,1:ndn)
  endif

  if (wf_type .eq. 'gutz_uhf') then
       if (e_type .eq. 1) then
          do i=1,nsites
                if (btest(det_in,i-1)) then
                     k=k+1
                     mat1(k,1:nup)=uhf_up_orbitals(i,1:nup)
                endif
          enddo
       endif

       if (e_type .eq. -1) then
          do i=1,nsites
                if (btest(det_in,i-1)) then
                     k=k+1
                     mat1(k,1:ndn)=uhf_dn_orbitals(i,1:ndn)
                endif
          enddo
       endif

      mat2=mat1
      if (e_type .eq.  1) mat2(chosen_e_row,1:nup)=uhf_up_orbitals(target_site,1:nup)
      if (e_type .eq. -1) mat2(chosen_e_row,1:ndn)=uhf_dn_orbitals(target_site,1:ndn)
  endif

  call determinant(mat1,n_e,det1)
  call determinant(mat2,n_e,det2)

  ratio=det2/det1

  end subroutine explicit_ratio_hubbard

  !=================================================================================================
  subroutine explicit_det_hubbard(e_type,n_e,det_in,chosen_e_row,target_site,detval)
  !=================================================================================================
    implicit none

    integer                  :: i,k
    integer, intent(in)      :: e_type,n_e
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_in
#else
    integer(ik), intent(in) :: det_in
#endif
    real(rk)                 :: mat1(n_e,n_e),mat2(n_e,n_e)
    integer,intent(in)       :: chosen_e_row,target_site
    real(rk),intent(out)     :: detval

   k=0

   if (wf_type .eq. 'gutz_rhf') then
       if (e_type .eq. 1) then
           do i=1,nsites
                if (btest(det_in,i-1)) then
                     k=k+1
                     mat1(k,1:nup)=rhf_up_orbitals(i,1:nup)
                endif
           enddo
       endif

       if (e_type .eq. -1) then
           do i=1,nsites
                if (btest(det_in,i-1)) then
                     k=k+1
                     mat1(k,1:ndn)=rhf_dn_orbitals(i,1:ndn)
                endif
           enddo
       endif
   endif

   if (wf_type .eq. 'gutz_uhf') then
       if (e_type .eq. 1) then
           do i=1,nsites
                if (btest(det_in,i-1)) then
                     k=k+1
                     mat1(k,1:nup)=uhf_up_orbitals(i,1:nup)
                endif
           enddo
       endif

       if (e_type .eq. -1) then
           do i=1,nsites
                if (btest(det_in,i-1)) then
                     k=k+1
                     mat1(k,1:ndn)=uhf_dn_orbitals(i,1:ndn)
                endif
           enddo
       endif
   endif

   mat2=mat1

   if (wf_type .eq. 'gutz_rhf') then
        if (e_type .eq. 1) mat2(chosen_e_row,1:nup)=rhf_up_orbitals(target_site,1:nup)
        if (e_type .eq. -1)mat2(chosen_e_row,1:ndn)=rhf_dn_orbitals(target_site,1:ndn)
   endif

   if (wf_type .eq. 'gutz_uhf') then
        if (e_type .eq. 1) mat2(chosen_e_row,1:nup)=uhf_up_orbitals(target_site,1:nup)
        if (e_type .eq. -1)mat2(chosen_e_row,1:ndn)=uhf_dn_orbitals(target_site,1:ndn)
   endif

  call determinant(mat2,n_e,detval)

  end subroutine explicit_det_hubbard


  !=================================================================================================
  subroutine det_ratio_hubbard(e_type,n_e,mat_inv,chosen_e_row,target_site,ratio)
  !=================================================================================================
  ! Calculate ratio using Sherman-Morrison-Woodbury (need only inverse and one row)
  ! Right now no Sherman Morrison for multi det wavefunction
  ! Hitesh Changlani

    implicit none

    integer               :: j
    integer,intent(in)    :: e_type
    integer, intent(in)   :: n_e
    real (rk), intent(in) :: mat_inv(n_e,n_e)
    integer,intent(in)    :: chosen_e_row,target_site
    real(rk),intent(out)  :: ratio

    ratio=0._rk

    if (wf_type .eq. 'gutz_rhf') then
        if (e_type .eq. 1) then
            do j=1,n_e
                ratio = ratio + (mat_inv(j,chosen_e_row)*rhf_up_orbitals(target_site,j))
            enddo
        endif

        if (e_type .eq. -1) then
            do j=1,n_e
                ratio = ratio + (mat_inv(j,chosen_e_row)*rhf_dn_orbitals(target_site,j))
            enddo
        endif
    endif

    if (wf_type .eq. 'gutz_uhf') then
        if (e_type .eq. 1) then
            do j=1,n_e
                ratio = ratio + (mat_inv(j,chosen_e_row)*uhf_up_orbitals(target_site,j))
            enddo
        endif

        if (e_type .eq. -1) then
            do j=1,n_e
                ratio = ratio + (mat_inv(j,chosen_e_row)*uhf_dn_orbitals(target_site,j))
            enddo
        endif
    endif

  end subroutine det_ratio_hubbard

  !=================================================================================================
  subroutine energy_pieces_hubbard(det_i_up, det_i_dn, e_mix_numerator, e_mix_denominator, imp_sampling)
    !-----------------------------------------------------------------------------------------------
    ! Description : Calculate pieces of the local energy for det_i
    !               Numerator   is sum_j=1^{ndet_psi_t} H_ij * trial_wf_weight on det_j
    !               Denominator is trial_wf_weight on det_i
    !               Currently only being used at start of run (local energies are being stored instead)
    !
    ! Created     : H.J.Changlani, 21 Nov 2010 (based on FP's code....(can be made faster)
    !-----------------------------------------------------------------------------------------------
    implicit none

    !dummy variables
    real(rk), intent(out)   :: e_mix_numerator,e_mix_denominator
    integer,optional        :: imp_sampling
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up,det_i_dn
    type(ik_vec)             :: det_i_up_new,det_i_dn_new
    type(ik_vec),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#else
    integer(ik), intent(in) :: det_i_up,det_i_dn
    integer(ik)             :: det_i_up_new,det_i_dn_new
    integer(ik),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#endif

    !local variables
    integer                 :: j,double_occ,location,n_connected_dets
    integer                 :: chosen_e_row(4*(nup+ndn)),chosen_site(4*(nup+ndn)),target_site(4*(nup+ndn)),e_types(4*(nup+ndn))
    integer                 :: diff_double_occs(4*(nup+ndn))
    logical                 :: imp
    real(rk)                :: up_inv(nup,nup),dn_inv(ndn,ndn)
    real(rk)                :: ratio,up_detval,dn_detval,wf_detval,wf_detval_new,gtmp
    real(rk)                :: matrix_element
    real(rk),allocatable    :: matrix_elements(:)
    logical                 :: entry_present

    if (present(imp_sampling)) then
        if (imp_sampling .eq. 0 ) imp=.false.
        if (imp_sampling .eq. 1 ) imp=.true.
    else
        imp=.false.
    endif

    e_mix_numerator   = 0._rk
    e_mix_denominator = 0._rk

    if ((main_wf_type .eq. 'gutz') .and. (imp.eqv. .false.)) then
        !print *,"Here (no imp)"
        call get_connected_dets_and_info_hubbard(  det_i_up,det_i_dn,chosen_e_row,chosen_site,target_site, &
                                                 & e_types,up_inv,dn_inv,up_detval,dn_detval,wf_detval,    &
                                                 & double_occ,diff_double_occs)
        gtmp=g**real(double_occ)
        e_mix_numerator = e_mix_numerator +(U*real(double_occ)*wf_detval)
        do j=1,4*(nup+ndn)
             if ((wf_type .ne. 'gutz_multi') .and. (wf_type .ne. 'cgutz_multi')) then
                 if (e_types(j) .eq. 1) then
                     call explicit_det_hubbard( 1,nup,det_i_up,chosen_e_row(j),target_site(j),ratio)
                     e_mix_numerator = e_mix_numerator -t*ratio*dn_detval*(g**real(diff_double_occs(j)))
                elseif (e_types(j) .eq. -1) then
                     call explicit_det_hubbard(-1,ndn,det_i_dn,chosen_e_row(j),target_site(j),ratio)
                     e_mix_numerator = e_mix_numerator -t*ratio*up_detval*(g**real(diff_double_occs(j)))
                endif
            else
                    det_i_up_new=det_i_up
                    det_i_dn_new=det_i_dn
                    if (e_types(j) .eq. 1) then
                        det_i_up_new=ibset(det_i_up_new,target_site(j)-1)
                        det_i_up_new=ibclr(det_i_up_new,chosen_site(j)-1)
                        call wf_calc(det_i_up_new,det_i_dn_new,wf_detval_new,.false.)
                        call hamiltonian_hubbard(det_i_up,det_i_dn,det_i_up_new,det_i_dn_new,matrix_element)
                        e_mix_numerator = e_mix_numerator+(matrix_element*(wf_detval_new)*(g**real(diff_double_occs(j))))
                    elseif (e_types(j) .eq. -1) then
                        det_i_dn_new=ibset(det_i_dn_new,target_site(j)-1)
                        det_i_dn_new=ibclr(det_i_dn_new,chosen_site(j)-1)
                        call wf_calc(det_i_up_new,det_i_dn_new,wf_detval_new,.false.)
                        call hamiltonian_hubbard(det_i_up,det_i_dn,det_i_up_new,det_i_dn_new,matrix_element)
                        e_mix_numerator = e_mix_numerator+(matrix_element*(wf_detval_new)*(g**real(diff_double_occs(j))))
                    endif
            endif
        enddo
        e_mix_numerator=e_mix_numerator*gtmp
        e_mix_denominator=wf_detval*gtmp

    elseif ((main_wf_type .eq. 'gutz') .and. (imp.eqv. .true.)) then                ! Trial = Guiding !If wavefunction is zero where it should not be there is a bias
        call get_connected_dets_and_info_hubbard(det_i_up,det_i_dn,chosen_e_row,chosen_site,target_site, &
                                                 & e_types,up_inv,dn_inv,up_detval,dn_detval,wf_detval,  &
                                                 & double_occ,diff_double_occs)
        e_mix_numerator = e_mix_numerator +(U*real(double_occ))
        do j=1,4*(nup+ndn)
             if ((wf_type .ne. 'gutz_multi') .and. (wf_type .ne. 'cgutz_multi')) then
                    if (e_types(j) .eq. 1) then
                         call det_ratio_hubbard(1,nup,up_inv,chosen_e_row(j),target_site(j),ratio)
                         e_mix_numerator = e_mix_numerator -t*ratio*(g**real(diff_double_occs(j)))
                    elseif (e_types(j) .eq. -1) then
                         call det_ratio_hubbard(-1,ndn,dn_inv,chosen_e_row(j),target_site(j),ratio)
                         e_mix_numerator = e_mix_numerator -t*ratio*(g**real(diff_double_occs(j)))
                    endif
             else
                    det_i_up_new=det_i_up
                    det_i_dn_new=det_i_dn
                    if (e_types(j) .eq. 1) then
                        det_i_up_new=ibset(det_i_up_new,target_site(j)-1)
                        det_i_up_new=ibclr(det_i_up_new,chosen_site(j)-1)
                        call wf_calc(det_i_up_new,det_i_dn_new,wf_detval_new,.false.)  ! .false. just gives the determinantal part of the wavefunction
                        call hamiltonian_hubbard(det_i_up,det_i_dn,det_i_up_new,det_i_dn_new,matrix_element)
                        e_mix_numerator = e_mix_numerator + (matrix_element*(wf_detval_new/wf_detval)*(g**real(diff_double_occs(j))))
                    elseif (e_types(j) .eq. -1) then
                        det_i_dn_new=ibset(det_i_dn_new,target_site(j)-1)
                        det_i_dn_new=ibclr(det_i_dn_new,chosen_site(j)-1)
                        call wf_calc(det_i_up_new,det_i_dn_new,wf_detval_new,.false.)  ! .false. just gives the determinantal part of the wavefunction
                        call hamiltonian_hubbard(det_i_up,det_i_dn,det_i_up_new,det_i_dn_new,matrix_element)
                        e_mix_numerator = e_mix_numerator + (matrix_element*(wf_detval_new/wf_detval)*(g**real(diff_double_occs(j))))
                    endif
            endif
        enddo
        !if (wf_type .eq. 'gutz_multi') e_mix_numerator=e_mix_numerator/wf_detval
        e_mix_denominator=1._rk

    else
        call find_connected_dets_hubbard(det_i_up,det_i_dn,n_connected_dets,connected_dets_up,connected_dets_dn,connected_matrix_elements)
        ! Very inefficient method of checking. Possibly use properties of Hubbard model! Comment by HJC on May 9 2011
        do j = 1, ndet_psi_t
              call is_in_list_mod(dets_up_psi_t(j),dets_dn_psi_t(j),connected_dets_up,connected_dets_dn, &
                                 & n_connected_dets,entry_present,location)
              if (entry_present) then
                 e_mix_numerator = e_mix_numerator + (matrix_elements(location) * cdet_psi_t(j))
              endif
              if (dets_up_psi_t(j) == det_i_up .and. dets_dn_psi_t(j) == det_i_dn) then
                 !det_i is in the trial wavefunction
                 e_mix_denominator = cdet_psi_t(j)
              endif
        enddo
        deallocate(connected_dets_up)
        deallocate(connected_dets_dn)
        deallocate(matrix_elements)
    endif

    !write (6,*) "e_mix_numerator=",e_mix_numerator
    !write (6,*) "e_mix_denominator=",e_mix_denominator
    !call flush(6)

  end subroutine energy_pieces_hubbard
  !==========================================================================================================================


!============================================================================================================================
  subroutine energy_pieces_hubbard_k(up_ket,dn_ket, e_mix_numerator, e_mix_denominator,importance_sampling)
    !------------------------------------------------------------------------------------------------------------------------
    ! Description : Calculate pieces of the local energy for det_i
    !               Numerator   is sum_j=1^{ndet_psi_t} H_ij * trial_wf_weight on det_j
    !               Denominator is trial_wf_weight on det_i
    ! Created     : A. Holmes, 29 Jun 2011 (based on Hitesh's and Frank's code)
    !------------------------------------------------------------------------------------------------------------------------
    use common_psi_t, only : psi_g_epsilon

    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: up_ket,dn_ket
#else
    integer(ik), intent(in) :: up_ket,dn_ket
#endif
    real(rk), intent(out) :: e_mix_numerator,e_mix_denominator
    integer,optional      :: importance_sampling

    !local variables
    integer            ::i,j,ind,up_site,target_up_site,dn_site
    integer            ::diff_up,diff_dn
    integer            ::updiff,dndiff
    real(rk)           ::e_mix_tmp,d_value
    logical            ::imp
    integer            ::nnzero

    imp=.false.
    if (present(importance_sampling)) then
        if (importance_sampling .eq. 1) imp=.true.
    endif

    e_mix_denominator = 0._rk
    e_mix_numerator   = 0._rk

    if (wf_type .eq. 'sym') then

      if (imp)  d_value=psi_g_epsilon  ! Why this??

      e_mix_denominator = 0._rk
      do i=1,ndet_psi_t ! How about a binary search? - AAH, 31 Oct 2012
        if (up_ket==dets_up_psi_t(i) .and. dn_ket==dets_dn_psi_t(i)) then
          if (imp) then
            e_mix_denominator=1._rk
            d_value=cdet_psi_t(i)  ! Currently using psi_guiding = psi_trial! Need to change this... - AAH, 31 Oct 2012.
          else
            e_mix_denominator=cdet_psi_t(i)
          endif
          exit
        endif
      enddo

      e_mix_numerator = 0._rk
      do i=1,ndet_psi_t
        if (space_sym) then
            call hamiltonian_hubbard_k_space_sym(dets_up_psi_t(i),dets_dn_psi_t(i),up_ket,dn_ket,e_mix_tmp,nnzero)
        else
            call hamiltonian_hubbard_k(dets_up_psi_t(i),dets_dn_psi_t(i),up_ket,dn_ket,e_mix_tmp)
        endif
        e_mix_numerator = e_mix_numerator + cdet_psi_t(i) * e_mix_tmp
      enddo

      if (imp)  e_mix_numerator = e_mix_numerator/d_value

    else                      ! Hartree Fock wavefunction or some modification

        diff_up=0
        diff_dn=0

         if (g .ne. 1.0) then ! we can make this into a recursive function call for higher order gutz later.

           updiff=0 ! number of sites that are different in up det
           dndiff=0 ! number of sites that are different in dn det
           do i=1,nsites
             if (btest(up_ket,i-1) .neqv. btest(k_hf_up,i-1))  updiff=updiff+1
             if (btest(dn_ket,i-1) .neqv. btest(k_hf_dn,i-1))  dndiff=dndiff+1
           enddo

           if (updiff==0 .and. dndiff==0) then
             e_mix_denominator=1._rk+ln_g*k_hf_matrix_elements(1)
           elseif (updiff==2 .and. dndiff==2) then
             do i=2,n_k_hf_con
               if (up_ket==k_hf_con_up(i) .and. dn_ket==k_hf_con_dn(i)) then
                 e_mix_denominator=ln_g*k_hf_matrix_elements(i)
                 exit
               endif
             enddo
           else
             e_mix_denominator = 0._rk
           endif

           call hamiltonian_hubbard_k(k_hf_up,k_hf_dn,up_ket,dn_ket,e_mix_numerator)
           e_mix_numerator = e_mix_numerator * (1._rk + ln_g*k_hf_matrix_elements(1))

           do i=2,n_k_hf_con
             call hamiltonian_hubbard_k(k_hf_con_up(i),k_hf_con_dn(i),up_ket,dn_ket,e_mix_tmp)
             e_mix_numerator = e_mix_numerator + e_mix_tmp*ln_g*k_hf_matrix_elements(i)
           enddo

        else  ! Assume closed shell Hartree Fock wavefunction
            diff_up=0
            diff_dn=0
            if (up_ket==k_hf_up .and. dn_ket==k_hf_dn)  then
                e_mix_denominator=1._rk
                e_mix_numerator=k_hf_energy
                return
            endif

            if (space_sym) then
                call binary_search(up_ket,dn_ket,rep_k_hf_con_up,rep_k_hf_con_dn,ind)
                e_mix_numerator=0._rk
                if (ind .ne. 0 ) e_mix_numerator=k_hf_matrix_elements(ind)
                e_mix_denominator=0._rk
                return
            else
                do j=1,nup
                          if (btest(up_ket,up_hf_location(j)-1) .eqv. .false.) then
                            up_site=j
                            diff_up=diff_up+1
                            if (diff_up .gt. 1) return
                          endif
               enddo

               diff_up=0
               do j=1,nsites-nup
                          if (btest(up_ket,up_hf_empty_location(j)-1)) then
                            target_up_site=j
                            diff_up=diff_up+1
                            if (diff_up .gt. 1) return
                          endif
               enddo

               do j=1,ndn
                          if (btest(dn_ket,dn_hf_location(j)-1) .eqv. .false.) then
                            dn_site=j
                            diff_dn=diff_dn+1
                            if (diff_dn .gt. 1) return
                          endif
               enddo
               e_mix_numerator=hf_matrix_elements(up_site,dn_site,target_up_site)
            endif

       endif

    endif

  end subroutine energy_pieces_hubbard_k

!============================================================================================================================
  subroutine energy_pieces_hubbard_dm(up_ket, dn_ket, e_mix_numerator, e_mix_denominator)
    !------------------------------------------------------------------------------------------------------------------------
    ! Description : Calculate pieces of the local energy for det given by up_ket,dn_ket
    ! Created     : H.J.Changlani Jan 11 2012
    !------------------------------------------------------------------------------------------------------------------------
    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: up_ket,dn_ket
#else
    integer(ik), intent(in) :: up_ket,dn_ket
#endif
    real(rk), intent(out)   :: e_mix_numerator,e_mix_denominator

    ! Local
    integer               :: i
    real(rk)              :: e_mix_num_tmp

    e_mix_denominator = 0._rk
    e_mix_numerator   = 0._rk

    do i=1,ndet_psi_t
        if (up_ket==dets_up_psi_t(i) .and. dn_ket==dets_dn_psi_t(i)) then
            e_mix_denominator=cdet_psi_t(i)
            exit
        endif
    enddo

    do i=1,ndet_psi_t
        call hamiltonian_hubbard_dm(dets_up_psi_t(i),dets_dn_psi_t(i),up_ket,dn_ket,e_mix_num_tmp)
        e_mix_numerator=e_mix_numerator+(cdet_psi_t(i)*e_mix_num_tmp)
    enddo
  end subroutine energy_pieces_hubbard_dm
  !=================================================================================================================================


  !================================================================================================================================
  subroutine generate_trial_wf_half_filled_hubbard()
    !---------------------------------------------------------------------------
    ! Description : Create trial wavefunction.
    ! Created     : H.J. Changlani 21 Nov 2010
    !---------------------------------------------------------------------------
    use common_run, only: tau_multiplier
    implicit none

    integer                 ::      num_configs
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), allocatable::      init_dets_up(:), init_dets_dn(:)
    type(ik_vec)             ::      neel_u_alpha_config,neel_u_beta_config
    type(ik_vec)             ::      up_config,dn_config
    type(ik_vec)             ::      new_config_1,new_config_2
    type(ik_vec)             ::      config,other_config
    type(ik_vec)                     full_config
#else
    integer(ik), allocatable::      init_dets_up(:), init_dets_dn(:)
    integer(ik)             ::      neel_u_alpha_config,neel_u_beta_config
    integer(ik)             ::      up_config,dn_config
    integer(ik)             ::      new_config_1,new_config_2
    integer(ik)             ::      config,other_config
    integer(ik)                     full_config
#endif
    real(rk), allocatable   ::      ham(:,:)
    integer                         nbr(4),i,j,j_ctr,k,location,pos,start,istat
    logical                         allowed,next_allowed,flg,up
    real(rk)                        lowest_eigenvalue
    integer                         start_max

    next_allowed=.true.
    allowed=.true.
    neel_u_alpha_config=0_ik
    neel_u_beta_config=0_ik

    allocate(init_dets_up(12*l_x*l_y+2),stat=istat)
    allocate(init_dets_dn(12*l_x*l_y+2),stat=istat)

    do i=1,l_y
            do j=1,l_x
                pos=(i-1)*(l_x)+j
                if (mod(i-1+j-1,2) .eq. 0) then
                    neel_u_alpha_config=ibset(neel_u_alpha_config,pos-1)
                else
                    neel_u_beta_config=ibset(neel_u_beta_config,pos-1)
                endif
            enddo
    enddo

    if (wf_type .eq. 'defect') then
        ! If only Neel up state is present then we simply consider defects
        ! around the Neel up state (i.e. we perform a symmetry breaking)
        if (neel_up_only) then
            start_max=1
        else
        ! If defects in both Neel states allowed
            start_max=2
        endif

        num_configs=1
        do start=1,start_max

            if (start .eq. 1) then
            up_config=neel_u_alpha_config
            dn_config=neel_u_beta_config
            init_dets_up(num_configs)=neel_u_alpha_config
            init_dets_dn(num_configs)=neel_u_beta_config
            write (6,*) "==================================="
            write (6,*) "         Neel - Up                 "
            write (6,*) "==================================="
            call print_walker(l_x,l_y,init_dets_up(num_configs),init_dets_dn(num_configs))
            print *,init_dets_up(num_configs),init_dets_dn(num_configs)
            num_configs=num_configs+1
            else
            dn_config=neel_u_alpha_config
            up_config=neel_u_beta_config
            init_dets_up(num_configs)=neel_u_beta_config
            init_dets_dn(num_configs)=neel_u_alpha_config
            write (6,*) "==================================="
            write (6,*) "         Neel - Down               "
            write (6,*) "==================================="
            call print_walker(l_x,l_y,init_dets_up(num_configs),init_dets_dn(num_configs))
            num_configs=num_configs+1
            endif

            do j_ctr=1,2*l_x*l_y
            j=(j_ctr+1)/2
            if (mod(j_ctr,2) .ne. 0 ) then
                config=up_config
                other_config=dn_config
                up=.true.
            else
                config=dn_config
                other_config=up_config
                up=.false.
            endif

            if (btest(config,j-1)) then
                !==================================================
                ! States obtained by 1 swap (1 defect)
                !==================================================
                do k=1,4
                call get_nbr(l_x,l_y,pbc,j,k,nbr(k),allowed)
                !print *,"j=",j
                !print *,"k=",k
                !print *,"nbr(k)=",nbr(k)
                next_allowed=.true.

                if (btest(config,j-1) .eqv. btest(config,nbr(k)-1)) then
                   allowed=.false.
                endif

                if (allowed .eqv. .true.) then

                    ! Swap only Right and down to avoid double counting
                    !---------------------------------------------------
                    if (mod(k,2) .ne. 0 ) then
                    new_config_1=ibclr(config,j-1)
                    new_config_1=ibset(new_config_1,nbr(k)-1)
                    new_config_2=ibset(other_config,j-1)
                    new_config_2=ibclr(new_config_2,nbr(k)-1)
                    flg=.false.

                    !if (l_x*l_y .le. 4) then
                        if (up .eqv. .true.) then
                        call is_in_list_mod(new_config_1,new_config_2,init_dets_up,init_dets_dn,num_configs,flg,location)
                        else
                        call is_in_list_mod(new_config_2,new_config_1,init_dets_up,init_dets_dn,num_configs,flg,location)
                        endif
                    !endif

                    if (flg .eqv. .false. ) then
                        if (up .eqv. .true.) then
                        init_dets_up(num_configs)=new_config_1
                        init_dets_dn(num_configs)=new_config_2
                        else
                        init_dets_up(num_configs)=new_config_2
                        init_dets_dn(num_configs)=new_config_1
                        endif
                        !print *,"init_dets_up_(swap)",init_dets_up(num_configs)
                        !print *,"init_dets_dn_(swap)",init_dets_dn(num_configs)
                        num_configs=num_configs+1

                    endif

                    next_allowed=.true.
                    else
                    if (nbr(k) .eq. nbr(k-1) ) then
                        next_allowed =.false.
                    endif
                    endif

                    !===============================================================
                    ! States obtained by 1 excitation (1 holon 1 doublon)
                    !===============================================================
                    if (next_allowed .eqv. .true.) then
                       new_config_1=ibclr(config,j-1)
                       new_config_1=ibset(new_config_1,nbr(k)-1)
                       flg=.false.


                    !if (l_x*l_y .le. 4) then
                        !print *,"Checking walker in list"
                        if (up .eqv. .true.) then
                        call print_walker(l_x,l_y,new_config_1,dn_config)
                        call is_in_list_mod(new_config_1,dn_config,init_dets_up,init_dets_dn,num_configs,flg,location)
                        else
                        call print_walker(l_x,l_y,up_config,new_config_1)
                        call is_in_list_mod(up_config,new_config_1,init_dets_up,init_dets_dn,num_configs,flg,location)
                        endif
                        !print *,"Flg=",flg
                    !endif

                    if (flg .eqv. .false. ) then
                        if (up .eqv. .true.) then
                        init_dets_up(num_configs)=new_config_1
                        init_dets_dn(num_configs)=dn_config
                        else
                        init_dets_up(num_configs)=up_config
                        init_dets_dn(num_configs)=new_config_1
                        endif
                        !print *,"init_dets_up_(excite)",init_dets_up(num_configs)
                        !print *,"init_dets_dn_(excite)",init_dets_dn(num_configs)
                        num_configs=num_configs+1

                    endif
                    endif
                endif
            enddo
            endif
        enddo

        enddo

        num_configs=num_configs-1
        print *,"Num_configs=",num_configs

        ndet_psi_t=num_configs

        allocate(dets_up_psi_t(ndet_psi_t),stat=istat)
        allocate(dets_dn_psi_t(ndet_psi_t),stat=istat)
        allocate(ham(ndet_psi_t,ndet_psi_t),stat=istat)
        allocate(cdet_psi_t(ndet_psi_t),stat=istat)

        dets_up_psi_t=init_dets_up(1:ndet_psi_t)
        dets_dn_psi_t=init_dets_dn(1:ndet_psi_t)

        write(6,'(2i22)') (dets_up_psi_t(i), dets_dn_psi_t(i), i=1,ndet_psi_t)
        !write (6,*) "Call diagonalize_hamiltonian_hubbard 1"
        !call flush(6)
        call diagonalize_hamiltonian_hubbard(dets_up_psi_t, dets_dn_psi_t, ham, cdet_psi_t, lowest_eigenvalue, tau_defect)
        if (e_trial .eq. 0.0_rk) e_trial=lowest_eigenvalue
        write(6,'(''From trial wavefn. e_trial='',f10.6)') e_trial

        deallocate(init_dets_up)
        deallocate(init_dets_dn)
        deallocate(ham)

    endif


    if ((wf_type .eq. 'neel') .and. (neel_up_only .eqv. .false.)) then
        ndet_psi_t=2
        allocate(cdet_psi_t(ndet_psi_t),stat=istat)
        allocate(dets_up_psi_t(ndet_psi_t),stat=istat)
        allocate(dets_dn_psi_t(ndet_psi_t),stat=istat)

        dets_up_psi_t(1)=neel_u_alpha_config
        dets_dn_psi_t(1)=neel_u_beta_config
        dets_up_psi_t(2)=neel_u_beta_config
        dets_dn_psi_t(2)=neel_u_alpha_config
        e_trial=0._rk

        cdet_psi_t(1)=0.707106781_rk
        cdet_psi_t(2)=0.707106781_rk
        print *,"Neel down configuration (in old notation) is",full_config
        call print_walker(l_x,l_y,dets_up_psi_t(2),dets_dn_psi_t(2))
        print *,"Neel up configuration (in old notation) is",full_config
        call print_walker(l_x,l_y,dets_up_psi_t(1),dets_dn_psi_t(1))

    endif

    if ((wf_type .eq. 'neel') .and. (neel_up_only .eqv. .true.)) then
        ndet_psi_t=1
        allocate(cdet_psi_t(ndet_psi_t),stat=istat)
        allocate(dets_up_psi_t(ndet_psi_t),stat=istat)
        allocate(dets_dn_psi_t(ndet_psi_t),stat=istat)

        e_trial=0._rk
        dets_up_psi_t(1)=neel_u_alpha_config
        dets_dn_psi_t(1)=neel_u_beta_config
        cdet_psi_t(1)=1._rk
        print *,"Neel up configuration (in old notation) is",full_config
        call print_walker(l_x,l_y,dets_up_psi_t(1),dets_dn_psi_t(1))

    endif

    end subroutine generate_trial_wf_half_filled_hubbard
    !=============================================================

  !==============================================================================================================================================
  subroutine diagonalize_hamiltonian_hubbard(dets_up, dets_dn, ham, lowest_eigenvector, lowest_eigenvalue, tau_out, print_or_not, sort_or_not)
  !==============================================================================================================================================
  !-----------------------------------------------------------------------------------------------------------------------------------------------
  ! Description : Construct and diagonalize H in a space of dets (broken up into dets_up and dets_dn)
  !               Sort ground state eigenvector by decreasing magnitude (HF first for example).
  !               Put dets_up and dets_dn in the same order.
  !
  ! Created     : H. J. Changlani, 28th Jan 2011 (borrowed from F. Petruzielo, 6 Nov 2010)
  !----------------------------------------------------------------------------------------------------------------------------------------------
    use common_run, only: tau_multiplier
    implicit none

    !dummy arguments
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(inout) :: dets_up(:),dets_dn(:)
#else
    integer(ik), intent(inout) :: dets_up(:),dets_dn(:)
#endif
    real(rk), intent(out) :: ham(:,:),lowest_eigenvector(:),lowest_eigenvalue
    real(rk),intent(out)  :: tau_out
    logical,optional      :: print_or_not
    logical,optional      :: sort_or_not

    !local variables
    real(rk),allocatable  :: eigenvectors(:,:)
    logical               :: ipr,srt
    integer               :: n_det
    integer               :: i, j
    real(rk)              :: tau_optimal_deterministic, tau_optimal_stochastic
    real(rk), parameter   :: epsilon = 1.e-8_rk
    real(rk), allocatable :: eigenvalues(:)

    if (present(print_or_not)) then
        if (print_or_not .eqv. .true.)   ipr=.true.
        if (print_or_not .eqv. .false.)  ipr=.false.
    else
        ipr=.true.
    endif

    if (present(sort_or_not)) then
        if (sort_or_not .eqv. .true.)   srt=.true.
        if (sort_or_not .eqv. .false. ) srt=.false.
    else
        srt=.true.
    endif

    n_det = size(dets_up)
    allocate(eigenvalues(n_det))
    allocate(eigenvectors(n_det,n_det))

    !construct hamiltonian
    ham = 0._rk
    do i = 1, n_det
       if ((hamiltonian_type.eq.'hubbard') .or. (hamiltonian_type .eq. 'hubbard2')) then
         call hamiltonian_hubbard(dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), ham(i,i))  ! diagonal element
       elseif (hamiltonian_type .eq. 'hubbardk') then
         call hamiltonian_hubbard_k(dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), ham(i,i))  ! diagonal element
       endif
       do j = i+1, n_det
         if ((hamiltonian_type.eq.'hubbard') .or. (hamiltonian_type .eq. 'hubbard2')) then
           if (is_connected_hubbard(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j))) then
             call hamiltonian_hubbard(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), ham(i,j))  ! off diagonal element
           else
             ham(i,j) = 0._rk
           endif
         elseif (hamiltonian_type .eq. 'hubbardk') then
           call hamiltonian_hubbard_k(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), ham(i,j))  ! off diagonal element
         endif
         !symmetric matrix
         ham(j,i) = ham(i,j)
       enddo
    enddo

    call real_symmetric_diagonalize(n_det,ham,eigenvectors,eigenvalues)

    if (size(eigenvalues).le.300) then
      write (6,*) "Eigenvalues"
      call print_real_matrix(size(eigenvalues),1,eigenvalues)
      call flush(6)
      write (6,*)
      write (6,*)
    endif

    ! Form symmetric linear combination of eigenvectors if degenerate
    if (size(eigenvalues)>1) then
        if (abs(eigenvalues(1)-eigenvalues(2)) .le. 1.0e-10) then
           lowest_eigenvector = 0.7071*(eigenvectors(:,1)+eigenvectors(:,2))
           lowest_eigenvalue = eigenvalues(1)
        else
        ! Retain lowest eigenvector
           lowest_eigenvector = eigenvectors(:,1)
           lowest_eigenvalue = eigenvalues(1)
       endif
   else
       lowest_eigenvector = eigenvectors(:,1)
       lowest_eigenvalue = eigenvalues(1)
   endif

   if (srt) call sort(lowest_eigenvector, dets_up, dets_dn)

   !if (ipr) then
   ! write(6,'(/,''Wavefunction has'',i9,'' determinants:'',/,'' det_up   det_dn,  cdet'')')  n_det
   ! write(6,'(2i22,f13.9)') (dets_up(i), dets_dn(i), lowest_eigenvector(i), i=1,n_det)
   !endif

! Calculate what tau_optimal_deterministic and tau_optimal_stochastic would be if the wavefn were this trial wavefn.
    tau_optimal_stochastic=1/(eigenvalues(n_det)-eigenvalues(1))
    tau_out=tau_multiplier*tau_optimal_stochastic

    if(n_det.ge.2) then
      tau_optimal_deterministic=2/(eigenvalues(n_det)+eigenvalues(2)-2*eigenvalues(1))
      !tau=tau_optimal_deterministic
      if (ipr) then
            write(6,'(/,''If the wavefn were this trial wf, tau_optimal_deterministic, tau_optimal_stochastic='',9f10.6)') &
           &  tau_optimal_deterministic, tau_optimal_stochastic
      endif
    else
      if (ipr) then
        write(6,'(/,''If the wavefn were this trial wf, tau_optimal_stochastic='',9f10.6)') tau_optimal_stochastic
      endif
    endif

!   tau=tau_optimal_stochastic*tau_multiplier  ! subspace range of eigenvalues is small compared to those for entire space so they cannot be used.

  end subroutine diagonalize_hamiltonian_hubbard
  !================================================================================================================


  !=====================================================================================================================================
  subroutine diagonalize_sparse_hamiltonian_hubbard(dets_up, dets_dn, lowest_eigenvector, lowest_eigenvalue, tau_out, print_or_not, sort_or_not,initial_vector)
  !====================================================================================================================================
  !------------------------------------------------------------------------------------------------------------------
  ! Description : Diagonalize H in a space of dets (broken up into dets_up and dets_dn)
  !               Sort ground state eigenvector by decreasing magnitude (HF first for example).
  !               Put dets_up and dets_dn in the same order.
  !
  ! Created     : A Holmes, 7 Apr 2012 (same as diagonalize_hamiltonian_hubbard, but stores H sparsely)
  ! Modified    : A Holmes, 14 Feb 2013. Optional input initial_vector is the starting vector for Lanczos.
  !------------------------------------------------------------------------------------------------------------------

    use common_run, only: tau_multiplier
    use tools, only : print_excitation_levels_and_wts
    use tools, only : merge_sort2_up_dn
    implicit none

    !dummy arguments
    real(rk), intent(out) :: lowest_eigenvector(:)
    real(rk), intent(out) :: lowest_eigenvalue
    real(rk),intent(out)  :: tau_out
    logical,optional      :: print_or_not
    logical,optional      :: sort_or_not
    real(rk),optional,intent(in) :: initial_vector(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(inout) :: dets_up(:)
    type(ik_vec), intent(inout) :: dets_dn(:)
    type(ik_vec) :: tmp_up,tmp_dn
#else
    integer(ik), intent(inout) :: dets_up(:)
    integer(ik), intent(inout) :: dets_dn(:)
    integer(ik) :: tmp_up,tmp_dn
#endif

    !local variables
    logical               :: ipr
    logical               :: srt
    integer :: n_connected_dets
    integer :: n_det
    integer :: i,j
    real(rk), parameter :: epsilon = 1.e-8_rk
    real(rk), allocatable :: H_values(:)
    integer(i8b), allocatable :: H_indices(:),H_nonzero_elements(:)
    real(rk) :: highest_eigenvalue,second_lowest_eigenvalue
    real(rk) tau_optimal_deterministic, tau_optimal_stochastic
    real(rk),allocatable :: ham(:,:),eigenvalues(:)
    real(rk) :: matrix_element
    integer  :: nnzero
    integer(i8b),parameter :: max_nonzero_elements = int(2e9,i8b)
    real(rk) :: log_num_nonzero_elements
   !integer(i8b) :: num_nonzero_elements
   !type(ik_vec),allocatable    :: temp_i16_up(:),temp_i16_dn(:)!,temp_i_2(:)
   !integer,allocatable        :: iorder(:),temp_i_2(:)

    if (present(print_or_not)) then
        if (print_or_not .eqv. .true.)   ipr=.true.
        if (print_or_not .eqv. .false. ) ipr=.false.
    else
        ipr=.true.
    endif

    if (present(sort_or_not)) then
        if (sort_or_not .eqv. .true.)   srt=.true.
        if (sort_or_not .eqv. .false. ) srt=.false.
    else
        srt=.true.
    endif

    n_det = size(dets_up)

    if (n_det.ge.2000) then
      write (6,'(''Performing matrix_lanczos because n_det>2000, n_det='',i8)') n_det

      ! sort by label so that a binary search can be performed
     !allocate(iorder(n_det))
     !allocate(temp_i16_up((n_det+1)/2))
     !allocate(temp_i16_dn((n_det+1)/2))
     !allocate(temp_i_2((n_det+1)/2))
     !do j=1,n_det
     !  iorder(j)=j
     !enddo

     !call merge_sort2_up_dn(dets_up,dets_dn, iorder, n_det, temp_i16_up, temp_i16_dn, temp_i_2)

     !deallocate(iorder)
     !deallocate(temp_i16_up)
     !deallocate(temp_i16_dn)
     !deallocate(temp_i_2)

      !construct hamiltonian
      call find_connected_dets_hubbard_k(dets_up(1),dets_dn(1),n_connected_dets,connected_dets_up,connected_dets_dn)
      write (6,'(''Number of nonzero elements in sparse H='',i10,'' ='',es9.2,'' Max number of nonzero elements='',i12,'' ='',es9.2,'' Number of connections to HF='',i12,'' ='',es8.1)') int(exp(log_num_nonzero_elements),i8b), real(exp(log_num_nonzero_elements)), max_nonzero_elements, real(max_nonzero_elements), n_connected_dets, real(n_connected_dets); call flush(6)
      ! Estimate whether we can store sparse matrix.
      log_num_nonzero_elements=log(real(n_det))+log(real(min(n_connected_dets,n_det)))-log(2.0)
     !call generate_sparse_ham_hubbardk_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,H_values,hf_to_psit=.false.,max_nonzero_elems=max_nonzero_elements,num_nonzero_elems=num_nonzero_elements)
      write (6,'(''Number of nonzero elements in sparse H='',i10,'' ='',es9.2,'' Max number of nonzero elements='',i12,'' ='',es9.2,'' Number of connections to HF='',i12,'' ='',es8.1)') int(exp(log_num_nonzero_elements)), real(exp(log_num_nonzero_elements)), max_nonzero_elements, real(max_nonzero_elements), n_connected_dets, real(n_connected_dets); call flush(6)
      if (exp(log_num_nonzero_elements)>max_nonzero_elements) then
        write (6,'(''Cannot store Hamiltonian sparsely for Lanczos, so use matrix_lanczos_on_the_fly'')') ; call flush(6)
       !n=size(dets_up)
       !allocate(temp_i16_up((n+1)/2))
       !allocate(temp_i16_dn((n+1)/2))
       !allocate(temp_i_2((n+1)/2))
       !allocate(iorder(n) )
       !call merge_sort2_up_dn(dets_up,dets_dn, iorder, n, temp_i16_up, temp_i16_dn, temp_i_2)
       !deallocate(iorder)
       !deallocate(temp_i16_up)
       !deallocate(temp_i16_dn)
       !deallocate(temp_i_2)
       !write (6,*) size(dets_up),size(dets_dn),size(lowest_eigenvector),lowest_eigenvalue
        if (present(initial_vector)) then
          call matrix_lanczos_on_the_fly(dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,highest_eigenvalue,second_lowest_eigenvalue,initial_vector=initial_vector)
        else
          call matrix_lanczos_on_the_fly(dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,highest_eigenvalue,second_lowest_eigenvalue,short_run_size=int(exp(log(real(max_nonzero_elements))-log_num_nonzero_elements+log(real(n_det)))))
        endif
       !call matrix_lanczos_on_the_fly(dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,highest_eigenvalue,second_lowest_eigenvalue)
      else
        write (6,*) "Using sparsely stored Hamiltonian for Lanczos"
        call generate_sparse_ham_hubbardk_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,H_values,hf_to_psit=.false.) ! When diagonalizing, we don't want to transform the projector first (hence hf_to_psit=false)
        call my_second(2, 'generate_sparse_ham_hubbardk')
        write (6,'(''Performing Lanczos using matrix_lanczos in more_tools, n_det='',i10)') n_det
        call flush(6)
        if (present(initial_vector)) then
          call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
        else
          call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue)
        endif
       !call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue)
        deallocate(H_values)
        deallocate(H_indices)
        deallocate(H_nonzero_elements)
      endif
    else
      write (6,'(''Performing Exact matrix diagonalization (instead of matrix_lanczos in more_tools) because n_det<2000, n_det='',i8)') n_det
      !construct H
      allocate(eigenvalues(n_det))
      allocate(ham(n_det,n_det))
      ham(:,:)=0._rk
      do i=1,n_det
        do j=1,i
          if (space_sym) then
            tmp_up = dets_up(i)
            tmp_dn = dets_dn(i)
            call hamiltonian_hubbard_k_space_sym(tmp_up,tmp_dn,dets_up(j),dets_dn(j),matrix_element,nnzero)
          else
            call hamiltonian_hubbard_k(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j),matrix_element)
          endif
          ham(i,j)=matrix_element
          if (i.ne.j)  ham(j,i)=matrix_element
        enddo
      enddo
      call my_second(2, 'generating full H')
      call flush(6)
      !diagonalize with lapack
      call real_symmetric_diagonalize_ow_ham(n_det,ham,eigenvalues)
      lowest_eigenvalue=eigenvalues(1)
      if (n_det>1)  second_lowest_eigenvalue=eigenvalues(2)
      highest_eigenvalue=eigenvalues(n_det)
      lowest_eigenvector=ham(:,1)
      deallocate(ham,eigenvalues)
    endif

    if (srt) call sort(lowest_eigenvector, dets_up, dets_dn)

!   write(6,'(/,''Wavefunction has'',i9,'' determinants. Dets. with largest'',i3,'' abs coefs are:'',/,'' det_up   det_dn    cdet'')')  n_det,min(n_det,20)
!   write(6,'(2b20,2i15,f13.9)') (dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), lowest_eigenvector(i),i=1,min(n_det,20))
!   write(6,'(/,''Lowest, highest eigenvalues='',9f13.6)') lowest_eigenvalue, highest_eigenvalue
!   call flush(6)

! Calculate number of dets with various excitation levels and the sum of their squared wts.
    if (space_sym .eqv. .false.) then
        call print_excitation_levels_and_wts(n_det,dets_up,dets_dn,lowest_eigenvector,nsites)
    else
        call print_excitation_levels_and_wts_hubbard_sym(n_det,dets_up,dets_dn,lowest_eigenvector,nsites)
    endif

! Calculate what tau_optimal_deterministic and tau_optimal_stochastic would be if the wavefn were this trial wavefn.
    tau_optimal_stochastic=1/(highest_eigenvalue-lowest_eigenvalue)
    tau_out=tau_multiplier*tau_optimal_stochastic
    if (n_det>1) then
      tau_optimal_deterministic=2/(highest_eigenvalue+second_lowest_eigenvalue-2*lowest_eigenvalue)
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_deterministic, tau_optimal_stochastic='',9f10.6)') &
   &  tau_optimal_deterministic, tau_optimal_stochastic
     else
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_stochastic='',9f10.6)') tau_optimal_stochastic
    endif

  end subroutine diagonalize_sparse_hamiltonian_hubbard
  !================================================================================================================

!==========================================================================================================================
  subroutine find_connected_dets_hubbard(det_up, det_dn, n_connected_dets, connected_dets_up, connected_dets_dn, &
                                         & matrix_elements)
!=========================================================================================================================
    !-------------------------------------------------------------------------------------------------------------
    ! Decription :  Return all determinants connected to det (broken into det_up and det_dn)
    !               Additionally, return all of the corresponding matrix elements.
    !
    ! Created     : H. J. Changlani, 21 Nov 2010 (stick to FP's code but do not need matrix element routine)
    !-------------------------------------------------------------------------------------------------------------
    implicit none

    !dummy argument
    real(rk), intent(out),allocatable    ::   matrix_elements(:)
    integer,intent(out)                  ::   n_connected_dets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in)              ::   det_up,det_dn
    type(ik_vec), intent(out),allocatable ::   connected_dets_up(:),connected_dets_dn(:)
    type(ik_vec)             ::   config,other_config
    type(ik_vec)             ::   tmp_connected_dets_up(4*(nup+ndn)+1),tmp_connected_dets_dn(4*(nup+ndn)+1)
#else
    integer(ik), intent(in)              ::   det_up,det_dn
    integer(ik), intent(out),allocatable ::   connected_dets_up(:),connected_dets_dn(:)
    integer(ik)             ::   config,other_config
    integer(ik)             ::   tmp_connected_dets_up(4*(nup+ndn)+1),tmp_connected_dets_dn(4*(nup+ndn)+1)
#endif

    ! local variables
    integer                 ::   ind,i_ctr,i,j,num_e,doubly_occupied
    integer                 ::   nbr_sites(nup+ndn,1+4)
    real(rk)                ::   f_s
    logical                 ::   up,allowed
    real(rk)                ::   tmp_matrix_elements(4*(nup+ndn)+1)

    if(det_up < 0) then
      write(6,'(''det_up is < 0, det_up (incoming)='',i20)') det_up
      stop 'det_up is < 0'
    endif
    if(det_dn < 0) then
      write(6,'(''det_dn is < 0, det_dn (incoming)='',i20)') det_dn
      stop 'det_dn is < 0'
    endif

    if (allocated(connected_dets_up)) then
       deallocate(connected_dets_up)
    endif
    if (allocated(connected_dets_dn)) then
       deallocate(connected_dets_dn)
    endif
    if (allocated(matrix_elements)) then
       deallocate(matrix_elements)
    endif

    num_e=0
    doubly_occupied=0
    up=.false.
    allowed=.false.
    n_connected_dets=0

    do i_ctr=1,2*nsites
        if (num_e==nup+ndn) exit
        i=(i_ctr+1)/2

        if (mod(i_ctr,2) .eq. 0) then
            config=det_dn
            other_config=det_up
            up=.false.
        else
            config=det_up
            other_config=det_dn
            up=.true.
        endif

        if (btest(config,i-1)) then

            num_e=num_e+1
            nbr_sites(num_e,1)=i

            if (up .and. btest(other_config,i-1)) then
              doubly_occupied=doubly_occupied+1
            endif

            do j=1,4                    ! Left Right Up Down
                if (ham_cross) then
                    call get_nbr_from_adj_list(adj_list_12_cross,i,j,nbr_sites(num_e,j+1),allowed)
                elseif (ham_square) then
                    call get_nbr_from_adj_list(adj_list_44_square,i,j,nbr_sites(num_e,j+1),allowed)
                else
                    call get_nbr(l_x,l_y,pbc,i,j,nbr_sites(num_e,j+1),allowed)
                endif

                ind=(4*(num_e-1))+j

                if (allowed .eqv. .true.) then
                        if (btest(config,nbr_sites(num_e,j+1)-1) .eqv. .false.) then
                                f_s=real(fermionic_phase(config,i,nbr_sites(num_e,j+1)))
                                tmp_matrix_elements(ind)=-t*f_s
                                if (up) then
                                    tmp_connected_dets_up(ind)=ibclr(config,i-1)
                                    tmp_connected_dets_up(ind)=ibset(tmp_connected_dets_up(ind),nbr_sites(num_e,j+1)-1)
                                    tmp_connected_dets_dn(ind)=det_dn
                                else
                                    tmp_connected_dets_dn(ind)=ibclr(config,i-1)
                                    tmp_connected_dets_dn(ind)=ibset(tmp_connected_dets_dn(ind),nbr_sites(num_e,j+1)-1)
                                    tmp_connected_dets_up(ind)=det_up
                                endif
                                !print *,"tmp_connected_dets_up (accepted)=",tmp_connected_dets_up(ind)
                                !print *,"tmp_connected_dets_dn (accepted)=",tmp_connected_dets_dn(ind)
                                !print *,"tmp_matrix element   (accepted) =",tmp_matrix_elements(ind)
                                n_connected_dets=n_connected_dets+1
                         else
                                tmp_connected_dets_up(ind)=-1_ik
                                tmp_connected_dets_dn(ind)=-1_ik
                                tmp_matrix_elements(ind)=0._rk
                        endif
                 else
                        tmp_connected_dets_up(ind)=-1_ik
                        tmp_connected_dets_dn(ind)=-1_ik
                        tmp_matrix_elements(ind)=0._rk
                 endif
            enddo
        endif
    enddo

    tmp_matrix_elements(4*(nup+ndn)+1)=U*doubly_occupied
    tmp_connected_dets_up(4*(nup+ndn)+1)=det_up
    tmp_connected_dets_dn(4*(nup+ndn)+1)=det_dn
    n_connected_dets=n_connected_dets+1

    allocate(connected_dets_up(n_connected_dets))
    allocate(connected_dets_dn(n_connected_dets))
    allocate(matrix_elements(n_connected_dets))

    !write(6,'(''n_connected_dets='',i10)') n_connected_dets
    ind=1
    do j=1,(4*(nup+ndn))+1
      if ((tmp_connected_dets_up(j) .ge. 0_ik) .and. (tmp_connected_dets_dn(j) .ge. 0_ik) ) then
         connected_dets_up(ind)=tmp_connected_dets_up(j)
         connected_dets_dn(ind)=tmp_connected_dets_dn(j)
         matrix_elements(ind)=tmp_matrix_elements(j)
    if(connected_dets_up(ind) < 0) then
      write(6,'(''connected_dets_up(ind) is < 0, det_up='',i20)') connected_dets_up(ind)
      stop 'connected_dets_up(ind) is < 0'
    endif
    if(connected_dets_dn(ind) < 0) then
      write(6,'(''connected_dets_dn(ind) is < 0, det_dn='',i20)') connected_dets_dn(ind)
      stop 'connected_dets_dn(ind) is < 0'
    endif
         ind=ind+1
      endif
    enddo

  end subroutine find_connected_dets_hubbard
  !=============================================================


!==========================================================================================================================
  subroutine find_connected_dets_hubbard_k(det_up, det_dn, n_connected_dets, connected_dets_up, connected_dets_dn, n_orb_include ,&
                                         & matrix_elements)!,imp_sampling)
!=========================================================================================================================
    !-------------------------------------------------------------------------------------------------------------
    ! Decription :  Return all determinants connected to det (broken into det_up and det_dn)
    !               Additionally, return all of the corresponding matrix elements.
    !
    ! Created     : A Holmes, 24 Jun 2011 (similar to Hitesh's code above)
    ! Modified    : A Holmes, 17 Jul 2012. Don't allocate connected_dets every time. Before calling this function, connected_dets_up,dn
    !               and matrix_elements (if used) should be allocated to nup*ndn*(nsites-nup)
    !-------------------------------------------------------------------------------------------------------------

    use more_tools, only : binary_search
    use types, only : i16b

    implicit none

    !dummy argument
    integer,intent(out)                           ::   n_connected_dets
    integer, intent(in),optional                  ::   n_orb_include
    real(rk), intent(out)            ,optional    ::   matrix_elements(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in)                       ::   det_up,det_dn
    type(ik_vec), intent(out)                      ::   connected_dets_up(:),connected_dets_dn(:)
    type(ik_vec)                                   ::   allowed_dets_up(nup*(nsites-nup)),tmp
    type(ik_vec)                                   ::   connected_dets_up_sym(nup*(nsites-nup)*ndn),connected_dets_dn_sym(nup*(nsites-nup)*ndn)
    type(ik_vec)                                   ::   rep_up,rep_dn
#else
    integer(ik), intent(in)                       ::   det_up,det_dn
    integer(ik), intent(out)                      ::   connected_dets_up(:),connected_dets_dn(:)
    integer(ik)                                   ::   allowed_dets_up(nup*(nsites-nup)),tmp
    integer(ik)                                   ::   connected_dets_up_sym(nup*(nsites-nup)*ndn),connected_dets_dn_sym(nup*(nsites-nup)*ndn)
    integer(ik)                                   ::   rep_up,rep_dn
#endif

    ! local variables
    integer                                       ::   i,j,k,dimup,q(nup*(nsites-nup),2),countup(nsites),countdn(nsites),countbkup(nsites),countbkdn(nsites),phaseup(nup*(nsites-nup))
    real(rk)                                      ::   matrix_elements_sym(nup*(nsites-nup)*ndn*(nsites-ndn))
    integer                                       ::   norb_include
    integer                                       ::   ctr,nnzero
    logical                                       ::   entry_present
    integer                                       ::   location

    if (present(n_orb_include)) then
      norb_include=n_orb_include
    else
      norb_include=nsites
    endif

    connected_dets_up(:) = 0
    connected_dets_dn(:) = 0

    ! Diagonal term
    connected_dets_up(1) = det_up
    connected_dets_dn(1) = det_dn
    if (present(matrix_elements))  matrix_elements(1)=0._rk
    countbkup(:)=0
    countbkdn(:)=0
    do i=norb_include,1,-1
      if (present(matrix_elements)) then
        if (btest(det_up,i-1))  matrix_elements(1) = matrix_elements(1) + k_energies(i)
        if (btest(det_dn,i-1))  matrix_elements(1) = matrix_elements(1) + k_energies(i)
      endif
      if (i<norb_include) then
        countbkup(i) = countbkup(i+1)
        countbkdn(i) = countbkdn(i+1)
      endif
      if (btest(det_dn,i-1))  countbkdn(i) = countbkdn(i) + 1
      if (btest(det_up,i-1))  countbkup(i) = countbkup(i) + 1
    enddo
    if (present(matrix_elements))  matrix_elements(1) = matrix_elements(1) + ubyn*real(nup)*real(ndn)

    ! Create list of up determinants that are one electron hop away from det_up,
    ! and their corresponding changes in momentum
    countup(:)=0
    countdn(:)=0
    do i=1,norb_include
      if (i>1) then
        countup(i) = countup(i-1)
        countdn(i) = countdn(i-1)
      endif
      if (btest(det_dn,i-1))  countdn(i) = countdn(i) + 1
      if (btest(det_up,i-1))  countup(i) = countup(i) + 1
    enddo

    dimup=0
    q(:,:)=0
    phaseup(:)=0
    do i=1,norb_include
      if (btest(det_up,i-1)) then
        do j=1,norb_include
          if (btest(det_up,j-1) .eqv. .false.) then
            tmp=ibclr(det_up,i-1)
            tmp=ibset(tmp,j-1)
#ifdef NUM_ORBITALS_GT_127
            if (tmp .gt. maskr_vec(norb_include))  cycle
#else
            if (ik==i16b .and. tmp .ge. 2_ik**norb_include)  cycle
#endif
            dimup=dimup+1
            allowed_dets_up(dimup)=tmp
            if (i>j) then
              phaseup(dimup) = (-1)**(countup(j)-countup(i))
            else
              phaseup(dimup) = (-1)**(countbkup(j)-countbkup(i))
            endif
            q(dimup,1) = -k_vectors(1,j) + k_vectors(1,i) ! q = change in momentum of down det
            q(dimup,2) = -k_vectors(2,j) + k_vectors(2,i)
          endif
        enddo
      endif
    enddo

    n_connected_dets=1
    do i=1,norb_include
      if (btest(det_dn,i-1)) then
        do j=1,dimup
          if ((q(j,1)>=-2*l_x+2 .and. q(j,1)<=2*l_x-2) .and. (q(j,2)>=-2*l_y+2 .and. q(j,2)<=2*l_y-2)) then
            if (kmap(i,int(q(j,1)/2)+l_x,int(q(j,2)/2)+l_y) /= 0) then ! if this q vector takes you to another valid site
              k=kmap(i,int(q(j,1)/2)+l_x,int(q(j,2)/2)+l_y)
              if (btest(det_dn,k-1) .eqv. .false.) then ! make sure the "to" site is empty
                !add this pair of determinants to the list and calculate corresponding energy.
                tmp=ibclr(det_dn,i-1)
                tmp=ibset(tmp,k-1)
#ifdef NUM_ORBITALS_GT_127
                if (tmp .gt. maskr_vec(norb_include))  cycle
#else
                if (ik==i16b .and. tmp .ge. 2_ik**norb_include)  cycle
#endif
                n_connected_dets = n_connected_dets + 1
                connected_dets_up(n_connected_dets) = allowed_dets_up(j)
                connected_dets_dn(n_connected_dets)=tmp
                !write(6,*) "n_connected dets = ",n_connected_dets
                call flush(6)
                if (present(matrix_elements)) then
                  if (i>k) then
                    matrix_elements(n_connected_dets) = ubyn * phaseup(j) * (-1)**(countdn(k)-countdn(i))
                  else
                    matrix_elements(n_connected_dets) = ubyn * phaseup(j) * (-1)**(countbkdn(k)-countbkdn(i))
                  endif
                endif
              endif
            endif
          endif
        enddo
      endif
    enddo

    !!! If space-time symmetry then only representatives needed
    if (space_sym) then
        ctr=0
        do i=1,n_connected_dets
            call get_rep_only(c4_map,reflection_map,z,p,connected_dets_up(i),connected_dets_dn(i),rep_up,rep_dn)
            call is_in_list_mod(rep_up,rep_dn,connected_dets_up_sym,connected_dets_dn_sym,ctr,entry_present,location)
            if (entry_present .eqv. .false.) then     ! If determinant is its own representative then compute matrix element
                ctr=ctr+1
                connected_dets_up_sym(ctr) = rep_up
                connected_dets_dn_sym(ctr) = rep_dn
                if (present(matrix_elements)) call hamiltonian_hubbard_k_space_sym(connected_dets_up_sym(ctr),connected_dets_dn_sym(ctr),det_up,det_dn,matrix_elements_sym(ctr),nnzero)
            endif
        enddo
        n_connected_dets=ctr
        connected_dets_up(1:n_connected_dets)=connected_dets_up_sym(1:n_connected_dets)
        connected_dets_dn(1:n_connected_dets)=connected_dets_dn_sym(1:n_connected_dets)
        if (present(matrix_elements)) matrix_elements(1:n_connected_dets)=matrix_elements_sym(1:n_connected_dets)
    endif

  end subroutine find_connected_dets_hubbard_k
  !=============================================================

!===========================================================================================================================
  subroutine find_connected_dets_imp_gutz_hubbard_sr(det_up, det_dn, n_connected_dets, connected_dets_up, connected_dets_dn, &
                                         & matrix_elements, matrix_elements_fn)
    !----------------------------------------------------------------------------------------------------------------------
    ! Decription :  Return all determinants connected to det (broken into det_up and det_dn)
    !               Additionally, return all of the corresponding matrix elements weighted by Gutzwiller function.
    !
    ! Created     : H. J. Changlani, 8 th April 2011
    !----------------------------------------------------------------------------------------------------------------------
    implicit none

    !dummy argument
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in)              ::   det_up,det_dn
    type(ik_vec), intent(out),allocatable ::   connected_dets_up(:),connected_dets_dn(:)
#else
    integer(ik), intent(in)              ::   det_up,det_dn
    integer(ik), intent(out),allocatable ::   connected_dets_up(:),connected_dets_dn(:)
#endif
    real(rk), intent(out),allocatable    ::   matrix_elements(:),matrix_elements_fn(:)
    integer,intent(out)                  ::   n_connected_dets

    !Local
    integer                              ::  chosen_e_row(4*(nup+ndn)),chosen_site(4*(nup+ndn)),target_site(4*(nup+ndn)),e_types(4*(nup+ndn))
    integer                              ::  diff_double_occs(4*(nup+ndn))
    integer                              ::  double_occ
    real(rk)                             ::  ratio
    real(rk)                             ::  up_inv(nup,nup),dn_inv(ndn,ndn)
    real(rk)                             ::  up_detval,dn_detval,wf_detval
    integer                              ::  j,ind,make_abs

    make_abs=0

    if (allocated(connected_dets_up))    deallocate(connected_dets_up)
    if (allocated(connected_dets_dn))    deallocate(connected_dets_dn)
    if (allocated(matrix_elements))      deallocate(matrix_elements)
    if (allocated(matrix_elements_fn))   deallocate(matrix_elements_fn)

    n_connected_dets=0
    call get_connected_dets_and_info_hubbard(det_up,det_dn,chosen_e_row,chosen_site,target_site,          &
                                                 & e_types,up_inv,dn_inv,up_detval,dn_detval,wf_detval,   &
                                                 & double_occ,diff_double_occs)
    do j=1,4*(nup+ndn)
             if (e_types(j) .eq. 1 .or. e_types(j).eq. -1) then
                 n_connected_dets=n_connected_dets+1
            endif
    enddo
    n_connected_dets=n_connected_dets+1

    allocate(connected_dets_up(n_connected_dets))
    allocate(connected_dets_dn(n_connected_dets))
    allocate(matrix_elements(n_connected_dets))
    allocate(matrix_elements_fn(n_connected_dets))

    matrix_elements_fn(1:n_connected_dets)=0._rk                                                          ! Initialize the diagonal element
    matrix_elements(1:n_connected_dets)=0._rk                                                             ! Initialize the diagonal element
    !write(6,'(''n_connected_dets='',i10)') n_connected_dets

    ind=1
    do j=1,(4*(nup+ndn))
         if (e_types(j) .eq. 1) then
             connected_dets_up(ind)=det_up
             connected_dets_up(ind)=ibclr(connected_dets_up(ind),chosen_site(j)-1)
             connected_dets_up(ind)=ibset(connected_dets_up(ind),target_site(j)-1)
             connected_dets_dn(ind)=det_dn
             call det_ratio_hubbard(1,nup,up_inv,chosen_e_row(j),target_site(j),ratio)
         elseif (e_types(j) .eq. -1) then
             connected_dets_dn(ind)=det_dn
             connected_dets_dn(ind)=ibclr(connected_dets_dn(ind),chosen_site(j)-1)
             connected_dets_dn(ind)=ibset(connected_dets_dn(ind),target_site(j)-1)
             connected_dets_up(ind)=det_up
             call det_ratio_hubbard(-1,ndn,dn_inv,chosen_e_row(j),target_site(j),ratio)
             !call explicit_det_hubbard(-1,ndn,det_dn,chosen_e_row(j),target_site(j),ratio)
         endif

         if (e_types(j) .eq. 1 .or. e_types(j) .eq. -1) then
            ratio=ratio*(g**real(diff_double_occs(j)))
            matrix_elements(ind)=-t*ratio
            if (matrix_elements(ind) .gt. 0.0_rk) then   ! Fixed node 1 - put extra weight in diagonal
                matrix_elements_fn(n_connected_dets)=matrix_elements_fn(n_connected_dets)+matrix_elements(ind)
                matrix_elements_fn(ind)=0._rk
            else
                matrix_elements_fn(ind)=matrix_elements(ind)
            endif
            ind=ind+1
         endif
    enddo

    connected_dets_up(n_connected_dets)=det_up
    connected_dets_dn(n_connected_dets)=det_dn
    matrix_elements(n_connected_dets)=matrix_elements(n_connected_dets)+U*real(double_occ)
    matrix_elements_fn(n_connected_dets)=matrix_elements_fn(n_connected_dets)+U*real(double_occ)

  end subroutine find_connected_dets_imp_gutz_hubbard_sr
  !====================================================================================================================


!===========================================================================================================================
  subroutine find_connected_dets_imp_gutz_hubbard(det_up, det_dn, n_connected_dets, connected_dets_up, connected_dets_dn, &
                                         & matrix_elements, fixed_node_option)
    !----------------------------------------------------------------------------------------------------------------------
    ! Decription :  Return all determinants connected to det (broken into det_up and det_dn)
    !               Additionally, return all of the corresponding matrix elements weighted by Gutzwiller function.
    !
    ! Created     : H. J. Changlani, 8 th April 2011
    !----------------------------------------------------------------------------------------------------------------------
    implicit none

    !dummy argument
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in)              ::   det_up,det_dn
    type(ik_vec), intent(out),allocatable ::   connected_dets_up(:),connected_dets_dn(:)
#else
    integer(ik), intent(in)              ::   det_up,det_dn
    integer(ik), intent(out),allocatable ::   connected_dets_up(:),connected_dets_dn(:)
#endif
    real(rk), intent(out),allocatable    ::   matrix_elements(:)
    integer,intent(out)                  ::   n_connected_dets
    integer,optional                     ::   fixed_node_option

    !Local
    integer                              ::  chosen_e_row(4*(nup+ndn)),chosen_site(4*(nup+ndn)),target_site(4*(nup+ndn)),e_types(4*(nup+ndn))
    integer                              ::  diff_double_occs(4*(nup+ndn))
    integer                              ::  double_occ
    real(rk)                             ::  ratio
    real(rk)                             ::  up_inv(nup,nup),dn_inv(ndn,ndn)
    real(rk)                             ::  up_detval,dn_detval
    real(rk)                             ::  wf_detval,wf_detval_new
    integer                              ::  j,ind
    integer                              ::  fn_option
    integer, save                        ::  count_negative=0, count_all=0

    fn_option=0
    if (present(fixed_node_option)) fn_option=fixed_node_option

    if (allocated(connected_dets_up)) deallocate(connected_dets_up)
    if (allocated(connected_dets_dn)) deallocate(connected_dets_dn)
    if (allocated(matrix_elements))   deallocate(matrix_elements)

    n_connected_dets=0
    call get_connected_dets_and_info_hubbard(det_up,det_dn,chosen_e_row,chosen_site,target_site,          &
                                                 & e_types,up_inv,dn_inv,up_detval,dn_detval,wf_detval,   &
                                                 & double_occ,diff_double_occs)
    do j=1,4*(nup+ndn)
             if (e_types(j) .eq. 1 .or. e_types(j).eq. -1) then
                 n_connected_dets=n_connected_dets+1
            endif
    enddo
    n_connected_dets=n_connected_dets+1

    allocate(connected_dets_up(n_connected_dets))
    allocate(connected_dets_dn(n_connected_dets))
    allocate(matrix_elements(n_connected_dets))

    matrix_elements(n_connected_dets)=0._rk                                                             ! Initialize the diagonal element
    !write(6,'(''n_connected_dets='',i10)') n_connected_dets

    ind=1
    do j=1,(4*(nup+ndn))
         if (e_types(j) .eq. 1) then
             connected_dets_up(ind)=det_up
             connected_dets_up(ind)=ibclr(connected_dets_up(ind),chosen_site(j)-1)
             connected_dets_up(ind)=ibset(connected_dets_up(ind),target_site(j)-1)
             connected_dets_dn(ind)=det_dn
             if ((wf_type .ne. 'gutz_multi') .and. (wf_type .ne. 'cgutz_multi')) then
                call det_ratio_hubbard(1,nup,up_inv,chosen_e_row(j),target_site(j),ratio)
                ratio=ratio*(g**real(diff_double_occs(j)))
                matrix_elements(ind)=-t*ratio
             else
                call wf_calc(connected_dets_up(ind),connected_dets_dn(ind),wf_detval_new,.false.)
                call hamiltonian_hubbard(det_up,det_dn,connected_dets_up(ind),connected_dets_dn(ind),matrix_elements(ind))
                matrix_elements(ind)=matrix_elements(ind)*(wf_detval_new/wf_detval)*(g**real(diff_double_occs(j)))
             endif
         elseif (e_types(j) .eq. -1) then
             connected_dets_dn(ind)=det_dn
             connected_dets_dn(ind)=ibclr(connected_dets_dn(ind),chosen_site(j)-1)
             connected_dets_dn(ind)=ibset(connected_dets_dn(ind),target_site(j)-1)
             connected_dets_up(ind)=det_up
             if ((wf_type .ne. 'gutz_multi') .and. (wf_type .ne. 'cgutz_multi')) then
                call det_ratio_hubbard(-1,ndn,dn_inv,chosen_e_row(j),target_site(j),ratio)
                ratio=ratio*(g**real(diff_double_occs(j)))
                matrix_elements(ind)=-t*ratio
             else
                call wf_calc(connected_dets_up(ind),connected_dets_dn(ind),wf_detval_new,.false.)
                call hamiltonian_hubbard(det_up,det_dn,connected_dets_up(ind),connected_dets_dn(ind),matrix_elements(ind))
                matrix_elements(ind)=matrix_elements(ind)*(wf_detval_new/wf_detval)*(g**real(diff_double_occs(j)))
             endif
         endif

         if (e_types(j) .eq. 1 .or. e_types(j) .eq. -1) then
            if ((fn_option==1) .and. (matrix_elements(ind) .gt. 0.0_rk)) then   ! Fixed node 1 - put extra weight in diagonal
                                                                                ! set off diagonal to 0
                matrix_elements(n_connected_dets)=matrix_elements(n_connected_dets)+matrix_elements(ind)
                matrix_elements(ind)=0._rk
            endif
            if (fn_option==2) then                                              ! Fixed node 2 - Make H_ij always negative
                matrix_elements(ind)=-abs(matrix_elements(ind))
            endif
            if ((fn_option==3) .and. (matrix_elements(ind) .gt. 0.0_rk)) then   ! Partial node - Kolodrubretz and Clark
                                                                                ! Need to check this
                matrix_elements(n_connected_dets)=matrix_elements(n_connected_dets)+(partial_node_eps*matrix_elements(ind))
                matrix_elements(ind)=(1._rk-partial_node_eps)*matrix_elements(ind)
            endif
            if (fn_option==4) then                                              ! Fixed node 4 - Make H_ij 0 is positive
                !write(6,*) "FN option 4"
                !call flush(6)
                if (matrix_elements(ind) .gt. 0._rk) then                       ! G. An, van Leeuwen 1991
                    matrix_elements(ind)=0._rk                                  ! off diagonal negative or zero
                endif
            endif
            ind=ind+1
         endif
    enddo

    connected_dets_up(n_connected_dets)=det_up
    connected_dets_dn(n_connected_dets)=det_dn
    matrix_elements(n_connected_dets)=matrix_elements(n_connected_dets)+U*real(double_occ)

    !write(6,'(''diagonal, wf_detval, double_occ='',2es10.2,i3)') matrix_elements(n_connected_dets), wf_detval, double_occ

! For linear projector, if elements of psi_g are arbitrarily small, then there is no finite value of tau that guarantees that
! the diagonal of the projector will never become negative.  Hence, instead of reducing tau, as we used to do, we keep track
! of how many negatives we encounter, and stop the run if it is too large.  It seems we can have many and still have a negligible bias.
    if (fn_option==1) then
        count_all=count_all+1
        if (1+tau*(e_trial-matrix_elements(n_connected_dets)) .lt. 0) then
          count_negative=count_negative+1
          if(count_negative.lt.10000) then
            write(6,'(''Warning: tau too large, diagonal becomes negative: count_negative, negative_frac, diagonal, wf_detval, double_occ='',i9,3es10.2,i3)') count_negative, real(count_negative)/count_all, matrix_elements(n_connected_dets), wf_detval, double_occ
          elseif(count_negative.eq.10000) then
            write(6,'(''Warning: tau too large, diagonal becomes negative: count_negative, negative_frac, diagonal, wf_detval, double_occ='',i9,3es10.2,i3)') count_negative, real(count_negative)/count_all, matrix_elements(n_connected_dets), wf_detval, double_occ
            write(6,'(''No more warning msgs. for tau too large will be written'')')
          endif
          call flush(6)
          if(count_negative.gt.0.001_rk*count_all+10000) then ! Add in 10000 because starting MC configs may have tiny Psi_G and tau is large so there are many negatives
            write(6,'(''Stopped because tau too large, diagonal becomes negative: count_negative, negative_frac, diagonal, wf_detval, double_occ='',i9,3es10.2,i3)') count_negative, real(count_negative)/count_all, matrix_elements(n_connected_dets), wf_detval, double_occ
            stop 'tau too large and maximum number of negative diagonal moves exceeded'
          endif
          matrix_elements(n_connected_dets)=e_trial+1/tau
        endif
    endif
  end subroutine find_connected_dets_imp_gutz_hubbard
  !====================================================================================================================

  function is_connected_hubbard(det_i_up, det_i_dn, det_j_up, det_j_dn)
    !-------------------------------------------------------------------------------------------------------------------
    ! Description : Determine if det_i is connected to det_j by the hamiltonian.
    !               In other words, is det_i=det_j (U) or are they related by a hop (t)
    !               Not presently being used (use is_connected_hubbard_fast instead)
    ! Note        : Very similar to hamiltonian_hubbard, BUT does nbr checking
    ! Created     : H. J. Changlani, 21 Nov 2010
    !-------------------------------------------------------------------------------------------------------------------

    !dummy arguments
    logical :: is_connected_hubbard
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up,det_i_dn,det_j_up,det_j_dn
    type(ik_vec) :: tmp_config_1,tmp_config_2,tmp_config
#else
    integer(ik), intent(in) :: det_i_up,det_i_dn,det_j_up,det_j_dn
    integer(ik) :: tmp_config_1,tmp_config_2,tmp_config
#endif

    !local variables
    integer :: k,n,pos_1,pos_2,nbr,num_ones
    logical :: allowed

    is_connected_hubbard=.false.

    !are determinants equal
    if (det_i_up == det_j_up .and. det_i_dn == det_j_dn) then
       is_connected_hubbard = .true.
       return
    endif

    ! Neither determinant equal - H cannot connect
    if ( (det_i_up .ne. det_j_up) .and. (det_i_dn .ne. det_j_dn) ) then
       is_connected_hubbard= .false.
       return
    endif

   !=====================================================
   ! Terms separated by 1 down electron hop
   !=====================================================

    if (det_i_up==det_j_up) then
       tmp_config_1=iand(det_i_dn,not(det_j_dn))
       tmp_config_2=iand(not(det_i_dn),det_j_dn)
       tmp_config=det_i_dn
    endif

   !=====================================================
   ! Terms separated by 1 up electron hop
   !=====================================================

   if (det_i_dn==det_j_dn) then
     tmp_config_1=iand(det_i_up,not(det_j_up))
     tmp_config_2=iand(not(det_i_up),det_j_up)
     tmp_config=det_i_up
   endif

   num_ones=0
   ! Check 1st tmp config for position 1
   do k=1,l_x*l_y
    if (btest(tmp_config_1,k-1)) then
        num_ones=num_ones+1
        pos_1=k
        if (num_ones .gt. 1) then
            is_connected_hubbard=.false.
            return
        endif
    endif
   enddo

   if (num_ones .eq. 0 ) then
     is_connected_hubbard=.false.
     return
   endif

   num_ones=0
   ! Check 2nd tmp config for position 2
   do k=1,l_x*l_y
    if (btest(tmp_config_2,k-1)) then
        num_ones=num_ones+1
        pos_2=k
        if (num_ones .gt. 1) then
            is_connected_hubbard=.false.
            return
        endif
    endif
   enddo

   if (num_ones .eq. 0 ) then
     is_connected_hubbard=.false.
     return
   endif

   do n=1,4
      call get_nbr(l_x,l_y,pbc,pos_1,n,nbr,allowed)
      if (allowed .and. pos_2 .eq. nbr) then
        is_connected_hubbard = .true.
        return
      endif
   enddo

  end function is_connected_hubbard
  !===========================================================================

  subroutine make_hubbard_matrix_2d(ham,lowest_eigenvec,lowest_eigenvalue,all_dets_up,all_dets_dn,tau_option)
   !--------------------------------------------------------------------------------------------------------
   ! Description   : Generate full Hubbard Hamiltonian matrix for 2D lattice given filling
   !                 Specify lattice dimensions l_x,l_y and boundary condition pbc
   !                 Also specify number f alpha and beta electrons n_alpha and n_beta
   !                 Hubbard parameters t (hopping),U (Coloumb repulsion)needed
   !                 Give size of hamiltonian = total_configs
   !                 Ham matrix is output
   ! Author        : H.J. Changlani
   ! ------------------------------------------------------------------------------------------------------
        implicit none

        real(rk),allocatable,intent(out)     ::  ham(:,:)
        real(rk),allocatable,intent(out)     ::  lowest_eigenvec(:)
        real(rk),intent(out)                 ::  lowest_eigenvalue
        logical,optional                     ::   tau_option
#ifdef NUM_ORBITALS_GT_127
        type(ik_vec), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        type(ik_vec), allocatable             ::  dets_up(:),dets_dn(:)
#else
        integer(ik), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        integer(ik), allocatable             ::  dets_up(:),dets_dn(:)
#endif

        !Local variables
        integer                              ::  i,j,istat
        logical                              ::  get_tau
        real(rk)                             ::  tau_out

        write (6,*) "nsites = ",nsites
        write (6,*) "nup    = ",nup
        write (6,*) "ndn    = ",ndn

        n_det_up=int(n_choose_k(nsites,nup),i4b)
        n_det_dn=int(n_choose_k(nsites,ndn),i4b)

        if (present(tau_option)) then
            get_tau=tau_option
        else
            get_tau=.false.
        endif

        allocate(dets_up(n_det_up),stat=istat)
        allocate(dets_dn(n_det_dn),stat=istat)
        allocate(all_dets_up(n_det_up*n_det_dn), stat=istat)
        allocate(all_dets_dn(n_det_up*n_det_dn), stat=istat)
        allocate(ham(n_det_up*n_det_dn,n_det_up*n_det_dn),stat=istat)
        allocate(lowest_eigenvec(n_det_up*n_det_dn),stat=istat)

        call constrained_dets(nup,n_det_up,dets_up)
        call constrained_dets(ndn,n_det_dn,dets_dn)

        do i=1,n_det_up
           do j=1,n_det_dn
             all_dets_up((i-1)*n_det_dn+j)=dets_up(i)
             all_dets_dn((i-1)*n_det_dn+j)=dets_dn(j)
           enddo
        enddo
        call diagonalize_hamiltonian_hubbard(all_dets_up,all_dets_dn,ham,lowest_eigenvec,lowest_eigenvalue,tau_out,.false.,.false.)
        if (get_tau) then
            tau_diag=tau_out
        endif

        deallocate(dets_up)
        deallocate(dets_dn)

   end subroutine make_hubbard_matrix_2d
  !===========================================================================


  subroutine make_hubbard_hamiltonian_k(ktot,ham,all_dets_up,all_dets_dn,lowest_eigenvec,lowest_eigenvalue)
   !-----------------------------------------------------------------------------------
   ! Description   : Generates the block of the k-space Hamiltonian matrix
   !                 with total momentum "ktot"; stores this block as "ham"
   !                 Also generates the lists of determinants with this total
   !                 momentum, "all_dets_up","all_dets_dn"
   !                 (Specify lattice dimensions l_x,l_y and boundary condition pbc
   !                 Also specify number of alpha and beta electrons n_alpha and n_beta
   !                 Hubbard parameters t (hopping),U (Coloumb repulsion)needed)
   !
   ! Author        : A. Holmes (inspired by Hitesh's subroutine above)
   ! ------------------------------------------------------------------------------------
        implicit none

        integer,intent(in)                   ::  ktot(2)
        real(rk),allocatable,intent(out)     ::  ham(:,:)
        real(rk),allocatable,intent(out)     ::  lowest_eigenvec(:)
        real(rk),intent(out)                 ::  lowest_eigenvalue
#ifdef NUM_ORBITALS_GT_127
        type(ik_vec), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        type(ik_vec), allocatable             ::  dets_up(:),dets_dn(:),all_dets_up_tmp(:),all_dets_dn_tmp(:)
        type(ik_vec),allocatable              ::  rep_all_dets_up(:),rep_all_dets_dn(:)
#else
        integer(ik), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        integer(ik), allocatable             ::  dets_up(:),dets_dn(:),all_dets_up_tmp(:),all_dets_dn_tmp(:)
        integer(ik),allocatable              ::  rep_all_dets_up(:),rep_all_dets_dn(:)
#endif

        !Local variables
        integer                              ::  i,j,istat,m,ms
        integer,allocatable                  ::  order(:)
        real(rk)                             ::  mat_elem
        integer                              ::  ktmp(2)
        real(rk),allocatable                     ::  eigenvalues(:),eigenvectors(:,:)
        integer                              ::  nnzero
        real(rk),allocatable                 ::  ham_s(:,:)
        real(rk)                             ::  second_lowest_eigenvalue,highest_eigenvalue
        mat_elem=0._rk

        n_det_up=int(n_choose_k(nsites,nup),i4b)
        n_det_dn=int(n_choose_k(nsites,ndn),i4b)

        allocate(order(n_det_up*n_det_dn))
        allocate(dets_up(n_det_up),stat=istat)
        allocate(dets_dn(n_det_dn),stat=istat)
        allocate(all_dets_up_tmp(n_det_up*n_det_dn), stat=istat)
        allocate(all_dets_dn_tmp(n_det_up*n_det_dn), stat=istat)

        call constrained_dets(nup,n_det_up,dets_up)
        call constrained_dets(ndn,n_det_dn,dets_dn)

        do i=1,n_det_up
          do j=1,n_det_dn
            all_dets_up_tmp((i-1)*n_det_dn+j)=dets_up(i)
            all_dets_dn_tmp((i-1)*n_det_dn+j)=dets_dn(j)
          enddo
        enddo

        ! m = number of combinations of up,dn dets with total momentum equal to ktot
        m=0
        do j=1,n_det_up*n_det_dn
          ktmp(:) = 0
          do i=1,nsites
            if (btest(all_dets_up_tmp(j),i-1)) ktmp = ktmp + k_vectors(:,i)
            if (btest(all_dets_dn_tmp(j),i-1)) ktmp = ktmp + k_vectors(:,i)
          enddo
          if ((mod(ktmp(1)-ktot(1),2*l_x)==0) .and. (mod(ktmp(2)-ktot(2),2*l_y)==0)) then
            m = m + 1
            order(m) = j
          endif
        enddo

        allocate(all_dets_up(m))
        allocate(all_dets_dn(m))

        do i=1,m
          all_dets_up(i)=all_dets_up_tmp(order(i))
          all_dets_dn(i)=all_dets_dn_tmp(order(i))
        enddo

        if (space_sym) then
            call symmetry_reduce_hubbardk(m,all_dets_up,all_dets_dn,rep_all_dets_up,rep_all_dets_dn)
            write (6,*) "Unsymmetrized list"
            do i=1,m
                write(6,'(9i6)') all_dets_up(i),all_dets_dn(i)
            enddo
            write (6,*)
            write (6,*) "Representative list"
            ms=size(rep_all_dets_up,1)
            do i=1,ms
                write(6,'(9i6)') rep_all_dets_up(i),rep_all_dets_dn(i)
            enddo
        endif
        ! now, calculate the Hamiltonian matrix "ham"
        allocate(ham(m,m))
        if (space_sym) allocate (ham_s(ms,ms))

        do i=1,m
          do j=1,m
               call hamiltonian_hubbard_k(all_dets_up(i),all_dets_dn(i),all_dets_up(j),all_dets_dn(j),ham(i,j))
          enddo
        enddo

        if (space_sym) then
            do i=1,ms
              do j=1,ms
                 call hamiltonian_hubbard_k_space_sym(rep_all_dets_up(i),rep_all_dets_dn(i),rep_all_dets_up(j),rep_all_dets_dn(j),ham_s(i,j),nnzero)
              enddo
            enddo
        endif


        if (space_sym) then
            !diagonalize with lapack routine
            allocate(eigenvalues(ms))
            allocate(eigenvectors(ms,ms))

            write (6,*) "Dimension of symmetrized Hamiltonian:",ms
            flush(6)
            if (ms<2000) then
                call real_symmetric_diagonalize(ms,ham_s,eigenvectors,eigenvalues)
            else
                if (allocated(lowest_eigenvec)) deallocate(lowest_eigenvec)
                allocate(lowest_eigenvec(ms))
                call matrix_lanczos(ms,lowest_eigenvec,lowest_eigenvalue,ham_s,highest_eigenvalue,second_lowest_eigenvalue)
            endif

            if (ms<200) then
              write(6,*)
              write(6,*) "Symmetrized Hamiltonian"
              call print_real_matrix(ms,ms,ham_s)
              write (6,*) "Eigenvalues of symmetrized Hamiltonian"
              call print_real_matrix(size(eigenvalues),1,eigenvalues)
              write (6,*)
              write (6,*)
              write(6,*) "Symmetrized Ground state eigenvector"
              call print_real_matrix(size(eigenvectors(:,1)),1,eigenvectors(:,1))
              write (6,*)
            else
              write (6,'(''Exact Ground State Energy:'',f10.6)') eigenvalues(1)
              write (6,*) "Eigenvalues of symmetrized Hamiltonian"
              call print_real_matrix(size(eigenvalues),1,eigenvalues)
            endif
            deallocate(eigenvalues)
            deallocate(eigenvectors)
        endif

        !diagonalize with lapack routine
        allocate(eigenvalues(m))
        allocate(eigenvectors(m,m))

        write (6,*)
        write (6,*) "Dimension of Hamiltonian:",m
        flush(6)
        if (m<2000) then
                call real_symmetric_diagonalize(m,ham,eigenvectors,eigenvalues)
                if (allocated(lowest_eigenvec)) deallocate(lowest_eigenvec)
                allocate(lowest_eigenvec(m),stat=istat)
                lowest_eigenvec(:)=eigenvectors(:,1)
                lowest_eigenvalue=eigenvalues(1)
                if (m<500) then
                  write(6,*)  "Unsymmetrized Hamiltonian"
                  call print_real_matrix(m,m,ham)
                  write (6,*) "Unsymmetrized Eigenvalues"
                  call print_real_matrix(size(eigenvalues),1,eigenvalues)
                  write (6,*)
                  write(6,*) "Unsymmetrized Ground state eigenvector"
                  call print_real_matrix(size(eigenvectors(:,1)),1,eigenvectors(:,1))
                  write (6,*)
                  write (6,*)
                else
                  write (6,'(''Exact Ground State Energy:'',f10.6)') eigenvalues(1)
                  write (6,*) "Eigenvalues of Hamiltonian"
                  call print_real_matrix(size(eigenvalues),1,eigenvalues)
                endif
        else
                if (allocated(lowest_eigenvec)) deallocate(lowest_eigenvec)
                allocate(lowest_eigenvec(m))
                call matrix_lanczos(m,lowest_eigenvec,lowest_eigenvalue,ham,highest_eigenvalue,second_lowest_eigenvalue)
        endif

        deallocate(eigenvalues)
        deallocate(eigenvectors)
        deallocate(order)
        deallocate(dets_up)
        deallocate(dets_dn)
        deallocate(all_dets_up_tmp)
        deallocate(all_dets_dn_tmp)

   end subroutine make_hubbard_hamiltonian_k
  !=====================================================================================
   subroutine lanczos_hubbard(lowest_eigenvec,lowest_eigenvalue,all_dets_up,all_dets_dn)
   !------------------------------------------------------------------------------------
   ! Description   : Perform Lanczos for Hubbard Hamiltonian 2D lattice given filling
   !                 in real or momentum space
   !                 Specify lattice dimensions l_x,l_y and boundary condition pbc
   !                 Also specify number f alpha and beta electrons n_alpha and n_beta
   !                 Hubbard parameters t (hopping),U (Coloumb repulsion)needed
   ! Author        : H.J. Changlani
   ! ------------------------------------------------------------------------------------
        implicit none

        real(rk),allocatable,intent(out)     ::  lowest_eigenvec(:)
        real(rk),intent(out)                 ::  lowest_eigenvalue
#ifdef NUM_ORBITALS_GT_127
        type(ik_vec), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        type(ik_vec), allocatable             ::  tmp_all_dets_up(:),tmp_all_dets_dn(:)
        type(ik_vec)                          ::  rep_up,rep_dn
        type(ik_vec)                          ::  det_i_up,det_i_dn
        type(ik_vec)                          ::  tmp_up,tmp_dn
#else
        integer(ik), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        integer(ik), allocatable             ::  tmp_all_dets_up(:),tmp_all_dets_dn(:)
        integer(ik)                          ::  rep_up,rep_dn
        integer(ik)                          ::  det_i_up,det_i_dn
        integer(ik)                          ::  tmp_up,tmp_dn
#endif
        !Local variables
        integer(ik),allocatable              ::  inverse_k_map(:),order(:)
        integer,allocatable                  ::  inverse_map_up(:),inverse_map_dn(:)
        integer                              ::  n_connected_dets
        integer                              ::  i,j,k,n,ncycles,ind_1,ind_2,it,vnum,istat,degen
        integer(ik)                           ::  ind
        integer                              ::  ktmp(2)
        integer                              ::  num_degen,iterations=10,num_cycles=25,npower=0
        real(rk)                             ::  norm,norm_inv,shift,rannyu
        real(rk),allocatable                 ::  w(:),alphas(:),betas(:),v(:,:)
        real(rk),allocatable                 ::  tri(:,:),eigenvalues(:)
        real(rk)                             ::  energy,lowest_eigenvalue_prev,lowest_eigenvalue_1,lowest_eigenvalue_2
        real(rk),allocatable                 ::  lowest_eigenvec_1(:),lowest_eigenvec_2(:)
        logical,allocatable                  ::  allowed_det(:)
        logical                              ::  converged=.false.
        integer                              ::  nnzero
        real(rk)                             ::  previous_e_t,e_t,matrix_element

        if (npower>iterations) npower=iterations -1              !A check

        !shift=U*real(min(nup,ndn))
        shift=0._rk
        e_t=e_trial
        previous_e_t=e_t
        write(6,*) "tau = ",tau
        call flush(6)

        if (ham_cross) then
            write (6,*) "Considering 12 site Hamiltonian"
            write (6,*) nup,ndn
            n_det_up=int(n_choose_k(12,nup),i4b)
            n_det_dn=int(n_choose_k(12,ndn),i4b)
            num_degen=1
        elseif (ham_square) then
            write (6,*) "Considering 16 site square Hamiltonian"
            write (6,*) nup,ndn
            n_det_up=int(n_choose_k(16,nup),i4b)
            n_det_dn=int(n_choose_k(16,ndn),i4b)
            num_degen=1
        else
            n_det_up=int(n_choose_k(l_x*l_y,nup),i4b)
            n_det_dn=int(n_choose_k(l_x*l_y,ndn),i4b)
            num_degen=1
        endif
        !ndet_hubbard=n_det_up*n_det_dn
        call flush(6)
        iterations=min(iterations,n_det_up*n_det_dn)

        allocate(tmp_all_dets_up(n_det_up),stat=istat)
        allocate(tmp_all_dets_dn(n_det_dn),stat=istat)
        allocate(alphas(iterations+1-npower),stat=istat)
        allocate(betas(iterations+1-npower),stat=istat)
        allocate(allowed_det(n_det_up*n_det_dn),stat=istat)
        allocate(inverse_k_map(n_det_up*n_det_dn),stat=istat)
        allocate(order(n_det_up*n_det_dn),stat=istat)
        allocate(tri(1,1),stat=istat)                             ! To avoid warning messages since allocated under if condition
        allocate(eigenvalues(1),stat=istat)                       ! To avoid warning messages since allocated under if condition

        if (allocated(stored_guiding_wf)) deallocate (stored_guiding_wf)
        allocate (stored_guiding_wf(n_det_up*n_det_dn))

        call constrained_dets(nup,n_det_up,tmp_all_dets_up)
        call constrained_dets(ndn,n_det_dn,tmp_all_dets_dn)

        if ((hamiltonian_type .eq. 'hubbard2').or. (ham_cross .eqv. .true.) .or. (ham_square .eqv. .true.)) then
            n=n_det_up*n_det_dn
            write (6,*) "In Lanczos with indexing: The total number of determinants in real space is",n
        elseif (hamiltonian_type .eq. 'hubbardk') then
            !n = number of combinations of up,dn dets with total momentum equal to ktot
            n=0
            do j=1,n_det_up*n_det_dn
              ind_1=((j-1)/n_det_dn)+1
              ind_2=j-(n_det_dn*(ind_1-1))
              ktmp(:) = 0
              do i=1,nsites
                if (btest(tmp_all_dets_up(ind_1),i-1)) ktmp = ktmp + k_vectors(:,i)
                if (btest(tmp_all_dets_dn(ind_2),i-1)) ktmp = ktmp + k_vectors(:,i)
              enddo
              if ((mod(ktmp(1)-ktot(1),2*l_x)==0) .and. (mod(ktmp(2)-ktot(2),2*l_y)==0)) then
                if (space_sym) then
                    call get_rep_only(c4_map,reflection_map,z,p,tmp_all_dets_up(ind_1),tmp_all_dets_dn(ind_2),rep_up,rep_dn)
                    if ((tmp_all_dets_up(ind_1) .eq. rep_up) .and. (tmp_all_dets_dn(ind_2) .eq. rep_dn)) then
                        n = n + 1
                        order(n) = j
                        allowed_det(j)=.true.
                    else
                        allowed_det(j)=.false.
                    endif
                else
                    n = n + 1
                    order(n) = j
                    allowed_det(j)=.true.
                endif
              else
                allowed_det(j)=.false.
              endif
            enddo
            write (6,*) "In Lanczos with indexing: The total number of determinants in momentum space is",n
            call flush(6)
        endif

        allocate(lowest_eigenvec(n),stat=istat)
        allocate(lowest_eigenvec_1(n),stat=istat)
        allocate(lowest_eigenvec_2(n),stat=istat)
        allocate(w(n),stat=istat)
        allocate(v(n,iterations+1-npower),stat=istat)
        allocate(inverse_map_up(2**nsites),stat=istat)
        allocate(inverse_map_dn(2**nsites),stat=istat)
        allocate(all_dets_up(n),stat=istat)
        allocate(all_dets_dn(n),stat=istat)

        v(:,:)=0._rk
        w(:)=0._rk

        if (hamiltonian_type .eq. 'hubbardk') then
            k=0
            do i=1,n_det_up*n_det_dn
                if (allowed_det(i) .eqv. .true.) then
                    k=k+1
                    ind_1=((i-1)/n_det_dn)+1
                    ind_2=i-(n_det_dn*(ind_1-1))
                    all_dets_up(k)=tmp_all_dets_up(ind_1)
                    all_dets_dn(k)=tmp_all_dets_dn(ind_2)
                    inverse_k_map(i)=k
                endif
            enddo
        endif

        if ((hamiltonian_type .eq. 'hubbard2') .or. (ham_cross .eqv. .true.) .or. (ham_square .eqv. .true.)) then
            do i=1,n_det_up*n_det_dn
                    ind_1=((i-1)/n_det_dn)+1
                    ind_2=i-(n_det_dn*(ind_1-1))
                    all_dets_up(i)=tmp_all_dets_up(ind_1)
                    all_dets_dn(i)=tmp_all_dets_dn(ind_2)
            enddo
        endif

        do i=1,n_det_up
#ifdef NUM_ORBITALS_GT_127
           inverse_map_up(tmp_all_dets_up(i)%v(1))=i
#else
           inverse_map_up(tmp_all_dets_up(i))=i
#endif
        enddo

        do i=1,n_det_dn
#ifdef NUM_ORBITALS_GT_127
           inverse_map_dn(tmp_all_dets_dn(i)%v(1))=i
#else
           inverse_map_dn(tmp_all_dets_dn(i))=i
#endif
        enddo

        if (ham_cross .or. ham_square) then
            allocate(inverse_map_up12(2**nsites),stat=istat)
            allocate(inverse_map_dn12(2**nsites),stat=istat)
            inverse_map_up12(:)=inverse_map_up(:)
            inverse_map_dn12(:)=inverse_map_dn(:)
        endif

        !print *,"K_HF_energy",k_hf_energy !print *,"nup=",nup !print *,"ndn=",ndn
        norm=0._rk
        do i=1,n
               det_i_up=all_dets_up(i)
               det_i_dn=all_dets_dn(i)
               if ((ham_cross .eqv. .true.) .or. (ham_square .eqv. .true.)) then
                    !write(6,*) "Ham cross",ham_cross !write(6,*) "Ham square",ham_square
                    v(i,1)=rannyu()
               else
                   call flush(6)
                   if ((hamiltonian_type .eq. 'hubbard2') .and. (main_wf_type .eq. 'gutz')) then ! variational wavefunction as start guess
                       call wf_calc(det_i_up,det_i_dn,v(i,1))
                       stored_guiding_wf(i)=v(i,1)
                   elseif (hamiltonian_type .eq. 'hubbardk') then
                     if (space_sym) then
                        tmp_up = det_i_up
                        tmp_dn = det_i_dn
                        call hamiltonian_hubbard_k_space_sym(tmp_up,tmp_dn,det_i_up,det_i_dn,energy,nnzero)
                     else
                        call hamiltonian_hubbard_k(det_i_up,det_i_dn,det_i_up,det_i_dn,energy)
                    endif
                     if (abs(energy-k_hf_energy)<1.0e-10) then                      ! If the determinant is another equal energy state
                         write(6,'(g14.6)') energy
                         v(i,1)=rannyu()
                     endif
                   else
                       v(i,1)=rannyu()                                              ! Else choose random vector as starting guess
                   endif
               endif
               norm=norm+v(i,1)*v(i,1)
        enddo

        if (abs(norm) .le. 1.0e-10) then
            write (6,*) "Norm of start vector should never be so small. Most likely your k space code has an inconsistency"
            write (6,*) "I am still proceeding, but dont expect the correct result"
            do i=1,n
                v(i,1)=rannyu()
            enddo
            norm=dot_product(v(:,1),v(:,1))
        endif

        norm_inv=1._rk/(norm**(0.5_rk))                                              ! Normalize v_1
        v(:,1)=v(:,1)*norm_inv

        betas(1)=0._rk

        write(6,'(/,''Diagonalizing in Lanczos_hubbard (with indexing)'')')

        do degen=1,num_degen
            write(6,'(/,''Computing state number '',i2)') degen
            do ncycles=1,num_cycles
                write (6,*)
                write(6,'(/,''In Cycle '',i2)') ncycles
                do it=1,iterations
                    call flush(6)
                    if (it .gt. npower) then
                        vnum=it-npower
                        !shift=0._rk
                    else
                        vnum=1
                    endif
                    call flush(6)
                    !Computing H*v_p - This is the bulk of the operation
                    if ((hamiltonian_type .eq. 'hubbard2') .or. (ham_cross .eqv. .true.) .or. (ham_square .eqv. .true.)) then
                         !write (6,*) "Ham_cross",ham_cross
                         do i=1,n
                            ind_1=((i-1)/n_det_dn)+1
                            ind_2=i-(n_det_dn*(ind_1-1))
                            det_i_up=tmp_all_dets_up(ind_1)
                            det_i_dn=tmp_all_dets_dn(ind_2)
                            if (abs(v(i,vnum)) .ge. w_min) then
                                call find_connected_dets_hubbard(det_i_up,det_i_dn,n_connected_dets,connected_dets_up,connected_dets_dn,connected_matrix_elements)
                                do k=1,n_connected_dets
#ifdef NUM_ORBITALS_GT_127
                                    ind_1=inverse_map_up(connected_dets_up(k)%v(1))
                                    ind_2=inverse_map_dn(connected_dets_dn(k)%v(1))
#else
                                    ind_1=inverse_map_up(connected_dets_up(k))
                                    ind_2=inverse_map_dn(connected_dets_dn(k))
#endif
                                    ind=(n_det_dn*(ind_1-1))+ind_2
                                    if (use_projector .eqv. .true.) then                            ! Use projector = 1 + tau (E_T-H)
                                        w(ind)=w(ind)-(tau*connected_matrix_elements(k)*v(i,vnum))  ! -tau H
                                        if (ind==i) then
                                            w(ind)=w(ind)+v(i,vnum)                                 ! Identity
                                            w(ind)=w(ind)+(tau*e_t*v(i,vnum))                       ! +tau E_T*idenity
                                        endif
                                    else
                                        w(ind)=w(ind)+(connected_matrix_elements(k)*v(i,vnum))      ! H
                                        if (ind==i) w(ind)=w(ind)-(shift*v(i,vnum))                 ! -S (shift)
                                    endif
                                enddo
                            endif
                        enddo
                    elseif (hamiltonian_type .eq. 'hubbardk') then
                        !write (6,*) "Action of Hamiltonian"
                        call flush(6)
                        do i=1,n
                            det_i_up=all_dets_up(i)
                            det_i_dn=all_dets_dn(i)
                            if (abs(v(i,vnum)) .ge. w_min) then                                     ! If an initiator
                                call find_connected_dets_hubbard_k(det_i_up,det_i_dn,n_connected_dets,connected_dets_up,connected_dets_dn,nsites,connected_matrix_elements)
                                do k=1,n_connected_dets
#ifdef NUM_ORBITALS_GT_127
                                    ind_1=inverse_map_up(connected_dets_up(k)%v(1))
                                    ind_2=inverse_map_dn(connected_dets_dn(k)%v(1))
#else
                                    ind_1=inverse_map_up(connected_dets_up(k))
                                    ind_2=inverse_map_dn(connected_dets_dn(k))
#endif
                                    ind=(n_det_dn*(ind_1-1))+ind_2
                                    ind=inverse_k_map(ind)
                                    if (use_projector .eqv. .true.) then                            ! Use projector = 1 + tau (E_T-H)
                                        w(ind)=w(ind)-(tau*connected_matrix_elements(k)*v(i,vnum))  ! -tau H
                                        if (ind==i) then
                                            w(ind)=w(ind)+v(i,vnum)                                 ! Identity
                                            w(ind)=w(ind)+(tau*e_t*v(i,vnum))                       ! +tau E_T*identity
                                        endif
                                    else
                                        w(ind)=w(ind)+(connected_matrix_elements(k)*v(i,vnum))
                                        if (ind==i) w(ind)=w(ind)-(shift*v(i,vnum))
                                    endif
                                enddo
                            else                                                                    ! If not an initiator
                                call hamiltonian_hubbard_k(det_i_up,det_i_dn,det_i_up,det_i_dn,matrix_element)
                                if (use_projector .eqv. .true.) then                                ! Use projector = 1 + tau (E_T-H)
                                        w(i)=w(i)-(tau*matrix_element*v(i,vnum))                    ! -tau H_ii
                                        w(i)=w(i)+v(i,vnum)                                         ! Identity
                                        w(i)=w(i)+(tau*e_t*v(i,vnum))                               ! +tau E_T*identity
                                endif
                            endif
                        enddo
                    endif

                    if (it .gt. npower) then
                        if (it .gt. npower+1) w(:)=w(:)-betas(vnum)*v(:,vnum-1)
                        alphas(vnum)=dot_product(w,v(:,vnum))
                        w(:)=w(:)-alphas(vnum)*v(:,vnum)
                        norm=dot_product(w,w)
                        betas(vnum+1)=norm**(0.5_rk)
                        norm_inv=1._rk/betas(vnum+1)
                        v(:,vnum+1)=w(:)*norm_inv
                        w(:)=v(:,vnum+1)
                        do i=1,it-npower
                            norm=dot_product(v(:,vnum+1),v(:,i))
                            !write (6,*) "q=",norm
                            call flush(6)
                            w(:)=w(:)-norm*v(:,i)
                        enddo
                        if (degen .eq. 2) then
                            norm=dot_product(v(:,vnum+1),lowest_eigenvec_1(:))
                            w(:)=w(:)-norm*lowest_eigenvec_1(:)
                        endif
                        v(:,vnum+1)=w(:)
                        w(:)=0._rk
                        norm=dot_product(v(:,vnum+1),v(:,vnum+1))
                        !write (6,*) "Norm after reorth = ",norm
                        norm_inv=1._rk/(norm**(0.5_rk))
                        v(:,vnum+1)=v(:,vnum+1)*norm_inv

                        if (allocated(tri)) deallocate(tri)
                        allocate(tri(it-npower,it-npower))
                        tri(:,:)=0._rk

                        if (allocated(eigenvalues)) deallocate(eigenvalues)
                        allocate(eigenvalues(it-npower))
                        eigenvalues(:)=0._rk

                        do i=1,it-npower
                           do j=1,it-npower
                            tri(i,j)=0._rk
                            enddo
                        enddo

                        do i=1,it-npower
                            tri(i,i)=alphas(i)
                            if (i<it-npower) then
                                tri(i,i+1)=betas(i+1)
                                tri(i+1,i)=betas(i+1)
                            endif
                        enddo

                        call real_symmetric_diagonalize_ow_ham(it-npower,tri,eigenvalues)

                        lowest_eigenvalue=eigenvalues(1)
                        !call print_real_matrix(size(eigenvalues),1,eigenvalues)
                        if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<1.0e-10) then
                            converged=.true.
                            exit
                        else
                            lowest_eigenvalue_prev=lowest_eigenvalue
                            write(6,'(''Iteration, Eigenvalue='',i3,f15.9)') it, lowest_eigenvalue+shift
                            call flush(6)

                            !if (it-npower>2) shift=-(eigenvalues(2)+eigenvalues(it-npower))/2._rk
                            if (it<iterations) then
                                deallocate(eigenvalues)
                                deallocate(tri)
                            endif
                         endif
                 else   ! Power method
                        e_t=dot_product(v(:,1),w(:)) ! Energy= <v|w>/<v|v>= <v|1+tau(E_t-H)|v>/<v|v> but v is already normalized
                        e_t=e_t-1._rk
                        e_t=e_t/tau
                        e_t=previous_e_t-e_t
                        previous_e_t=e_t
                        write(6,*) "E_T = ",e_t      ! Try power method for hubbardk to compare with Cyrus' observations
                        v(:,1)=w(:)                  ! for 8x8 13,13
                        w(:)=0._rk
                        norm=dot_product(v(:,1),v(:,1))
                        norm_inv=1._rk/(norm**(0.5_rk))
                        v(:,1)=v(:,1)*norm_inv
                 endif ! it > npower
               enddo ! it=1,iterations

               it=min(it,iterations)
               lowest_eigenvalue=lowest_eigenvalue+shift
               write(6,'(''n, Lowest eigenvalue ='',i10,f16.10)') n, lowest_eigenvalue
               call flush(6)

               !write(6,*) "iterations",it
               !write (6,*) "tridiagonal matrix rows",size(tri,1)
               !write (6,*) "tridiagonal matrix cols",size(tri,2)
               !write (6,*) v
               call flush(6)

               if (degen .eq. 1) lowest_eigenvec_1(:)=matmul(v(:,1:it-npower),tri(:,1))
               if (degen .eq. 2) lowest_eigenvec_2(:)=matmul(v(:,1:it-npower),tri(:,1))

               if (converged) then
                  exit
               else
                  if (degen .eq. 1) v(:,1)=lowest_eigenvec_1(:)
                  if (degen .eq. 2) v(:,1)=lowest_eigenvec_2(:)
               endif
               deallocate(eigenvalues)
               deallocate(tri)

             enddo ! cycles

             if (degen .eq. 1) then
                 lowest_eigenvalue_1=lowest_eigenvalue
                 norm=0._rk
                 do i=1,n
                    v(i,1)=rannyu()
                 enddo
                 norm=dot_product(v(:,1),lowest_eigenvec_1)
                 v(:,1)=v(:,1)-norm*lowest_eigenvec_1(:)
                 norm=dot_product(v(:,1),v(:,1))
                 norm_inv=1._rk/sqrt(norm)
                 v(:,1)=v(:,1)*norm_inv
                 converged=.false.                               ! Reset convergence
            endif

            if (degen .eq. 2) lowest_eigenvalue_2=lowest_eigenvalue
         enddo ! degeneracy

         lowest_eigenvec(:)=lowest_eigenvec_1(:)
         lowest_eigenvalue=lowest_eigenvalue_1

         deallocate(inverse_map_up)
         deallocate(inverse_map_dn)
         deallocate(w)
         deallocate(v)
         deallocate(alphas)
         deallocate(betas)

         call my_second(2,'lanczos_hubbard') ; call flush(6)

   end subroutine lanczos_hubbard

!-----------------------------------------------------------------------------------------

   subroutine arnoldi_hubbard_binary_search(lowest_eigenvec,lowest_eigenvalue,all_dets_up,all_dets_dn)
   !------------------------------------------------------------------------------------
   ! Description   : Perform Lanczos for Hubbard Hamiltonian 2D lattice given filling
   !                 in real or momentum space
   !                 Specify lattice dimensions l_x,l_y and boundary condition pbc
   !                 Also specify number f alpha and beta electrons n_alpha and n_beta
   !                 Hubbard parameters t (hopping),U (Coloumb repulsion)needed
   ! Author        : H.J. Changlani
   ! ------------------------------------------------------------------------------------
        implicit none

        real(rk),allocatable,intent(out)     ::  lowest_eigenvec(:)
        real(rk),intent(out)                 ::  lowest_eigenvalue
#ifdef NUM_ORBITALS_GT_127
        type(ik_vec), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        type(ik_vec), allocatable             ::  tmp_all_dets_up(:),tmp_all_dets_dn(:)
        type(ik_vec)                          ::  rep_up,rep_dn
        type(ik_vec)                          ::  det_i_up,det_i_dn
        type(ik_vec)                          ::  tmp_up,tmp_dn
#else
        integer(ik), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        integer(ik), allocatable             ::  tmp_all_dets_up(:),tmp_all_dets_dn(:)
        integer(ik)                          ::  rep_up,rep_dn
        integer(ik)                          ::  det_i_up,det_i_dn
        integer(ik)                          ::  tmp_up,tmp_dn
#endif

        !Local variables
        integer(ik),allocatable              ::  inverse_k_map(:),order(:)
        integer                              ::  n_connected_dets
        integer                              ::  c,i,j,k,n,ncycles,ind_1,ind_2,it,vnum,istat,degen
        integer(ik)                           ::  ind
       !type(ik_vec)                          ::  ind,rep_up,rep_dn
        integer                              ::  ktmp(2)
        integer                              ::  num_degen,iterations=50,num_cycles=25,npower=0
        real(rk)                             ::  norm,norm_inv,shift,rannyu
        real(rk),allocatable                 ::  w(:),v(:,:)
        real(rk)                             ::  full_hessen(1000,1000)
        real(rk),allocatable                 ::  hessen(:,:),eigenvalues(:),im_eigenvalues(:)
        real(rk)                             ::  energy,lowest_eigenvalue_prev,lowest_eigenvalue_1!,lowest_eigenvalue_2
        real(rk),allocatable                 ::  lowest_eigenvec_1(:),lowest_eigenvec_2(:)
        logical,allocatable                  ::  allowed_det(:)
        logical                              ::  converged=.false.
        integer                              ::  nnzero
        real(rk)                             ::  e_mix,e_mix_numerator,e_mix_denominator

        if (npower>iterations) npower=iterations -1              !A check

        !shift=U*real(min(nup,ndn))
        shift=0._rk

        n_det_up=int(n_choose_k(l_x*l_y,nup),i4b)
        n_det_dn=int(n_choose_k(l_x*l_y,ndn),i4b)
        num_degen=1
        call flush(6)
        iterations=min(iterations,n_det_up*n_det_dn)

        allocate(tmp_all_dets_up(n_det_up),stat=istat)
        allocate(tmp_all_dets_dn(n_det_dn),stat=istat)
        allocate(allowed_det(n_det_up*n_det_dn),stat=istat)
        allocate(inverse_k_map(n_det_up*n_det_dn),stat=istat)
        allocate(order(n_det_up*n_det_dn),stat=istat)
        allocate(hessen(1,1),stat=istat)                             ! To avoid warning messages since allocated under if condition
        allocate(eigenvalues(1),stat=istat)                          ! To avoid warning messages since allocated under if condition
        allocate(im_eigenvalues(1),stat=istat)                       ! To avoid warning messages since allocated under if condition

        call constrained_dets(nup,n_det_up,tmp_all_dets_up)
        call constrained_dets(ndn,n_det_dn,tmp_all_dets_dn)
        if (allocated(stored_guiding_wf)) deallocate (stored_guiding_wf)
        allocate (stored_guiding_wf(n_det_up*n_det_dn))

        if (hamiltonian_type .eq. 'hubbard2') then
            n=n_det_up*n_det_dn
            write (6,*) "In Binary search Arnoldi: The total number of determinants in real space is: ",n
        elseif (hamiltonian_type .eq. 'hubbardk') then
            !n = number of combinations of up,dn dets with total momentum equal to ktot
            n=0
            do j=1,n_det_up*n_det_dn
              ind_1=((j-1)/n_det_dn)+1
              ind_2=j-(n_det_dn*(ind_1-1))
              ktmp(:) = 0
              do i=1,nsites
                if (btest(tmp_all_dets_up(ind_1),i-1)) ktmp = ktmp + k_vectors(:,i)
                if (btest(tmp_all_dets_dn(ind_2),i-1)) ktmp = ktmp + k_vectors(:,i)
              enddo
              if ((mod(ktmp(1)-ktot(1),2*l_x)==0) .and. (mod(ktmp(2)-ktot(2),2*l_y)==0)) then
                if (space_sym) then
                    call get_rep_only(c4_map,reflection_map,z,p,tmp_all_dets_up(ind_1),tmp_all_dets_dn(ind_2),rep_up,rep_dn)
                    if ((tmp_all_dets_up(ind_1) .eq. rep_up) .and. (tmp_all_dets_dn(ind_2) .eq. rep_dn)) then
                        n = n + 1
                        order(n) = j
                        allowed_det(j)=.true.
                    else
                        allowed_det(j)=.false.
                    endif
                else
                    n = n + 1
                    order(n) = j
                    allowed_det(j)=.true.
                endif
              else
                allowed_det(j)=.false.
              endif
            enddo
            write (6,*) "In Binary search Arnoldi: The total number of determinants in momentum space is",n
            call flush(6)
        endif

        allocate(lowest_eigenvec(n),stat=istat)
        allocate(lowest_eigenvec_1(n),stat=istat)
        allocate(lowest_eigenvec_2(n),stat=istat)
        allocate(w(n),stat=istat)
        allocate(v(n,iterations+1-npower),stat=istat)
        allocate(all_dets_up(n),stat=istat)
        allocate(all_dets_dn(n),stat=istat)

        v(:,:)=0._rk
        w(:)=0._rk

        if (hamiltonian_type .eq. 'hubbardk') then
            k=0
            do i=1,n_det_up*n_det_dn
                if (allowed_det(i) .eqv. .true.) then
                    k=k+1
                    ind_1=((i-1)/n_det_dn)+1
                    ind_2=i-(n_det_dn*(ind_1-1))
                    all_dets_up(k)=tmp_all_dets_up(ind_1)
                    all_dets_dn(k)=tmp_all_dets_dn(ind_2)
                    inverse_k_map(i)=k
                endif
            enddo
        endif

        if (hamiltonian_type .eq. 'hubbard2') then
            do i=1,n_det_up*n_det_dn
                    ind_1=((i-1)/n_det_dn)+1
                    ind_2=i-(n_det_dn*(ind_1-1))
                    all_dets_up(i)=tmp_all_dets_up(ind_1)
                    all_dets_dn(i)=tmp_all_dets_dn(ind_2)
            enddo
        endif

        norm=0._rk
        do i=1,n
               det_i_up=all_dets_up(i)
               det_i_dn=all_dets_dn(i)
                   call flush(6)
                   if ((hamiltonian_type .eq. 'hubbard2') .and. (main_wf_type .eq. 'gutz')) then ! variational wavefunction as start guess
                     !if (global_fn_option .eq. -1) then
                     call wf_calc(det_i_up,det_i_dn,v(i,1))
                     stored_guiding_wf(i)=v(i,1)
                     v(i,1)=v(i,1)*v(i,1)
                     !else
                     !   v(i,1)=1._rk             ! Eigenvector is approximately psi**2 if importance sampling used
                     !endif
                   elseif (hamiltonian_type .eq. 'hubbardk') then
                     if (space_sym) then
                        tmp_up = det_i_up
                        tmp_dn = det_i_dn
                        call hamiltonian_hubbard_k_space_sym(tmp_up,tmp_dn,det_i_up,det_i_dn,energy,nnzero)
                     else
                        call hamiltonian_hubbard_k(det_i_up,det_i_dn,det_i_up,det_i_dn,energy)
                    endif
                     if (abs(energy-k_hf_energy)<1.0e-10) then                      ! If the determinant is another equal energy state
                         write(6,'(g14.6)') energy
                         v(i,1)=rannyu()
                     endif
                   else
                       v(i,1)=rannyu()                                              ! Else choose random vector as starting guess
                   endif
               norm=norm+v(i,1)*v(i,1)
        enddo

        if (abs(norm) .le. 1.0e-10) then
            write (6,*) "Norm of start vector should never be so small. Most likely your k space code has an inconsistency"
            write (6,*) "I am still proceeding, but dont expect the correct result"
            do i=1,n
                v(i,1)=rannyu()
            enddo
            norm=dot_product(v(:,1),v(:,1))
        endif

        norm_inv=1._rk/(norm**(0.5_rk))                                              ! Normalize v_1
        v(:,1)=v(:,1)*norm_inv

        write(6,'(/,''Diagonalizing in Arnoldi_hubbard (with binary search)'')')
        do degen=1,num_degen
            write(6,'(/,''Computing state number '',i2)') degen
            do ncycles=1,num_cycles
                write (6,*)
                write(6,'(/,''In Cycle '',i2)') ncycles
                do it=1,iterations
                    call flush(6)
                    if (it .gt. npower) then
                        vnum=it-npower
                        !shift=0._rk
                    else
                        vnum=1
                    endif
                    call flush(6)
                    !Computing H*v_p - This is the bulk of the operation
                    if (hamiltonian_type .eq. 'hubbard2') then
                         do i=1,n
                            ind_1=((i-1)/n_det_dn)+1
                            ind_2=i-(n_det_dn*(ind_1-1))
                            det_i_up=tmp_all_dets_up(ind_1)
                            det_i_dn=tmp_all_dets_dn(ind_2)
                            if (abs(v(i,vnum)) .gt. 1.0e-15) then
                                if (global_fn_option .ne. -1) then
                                    call find_connected_dets_imp_gutz_hubbard(det_i_up, det_i_dn, n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, global_fn_option)
                                else
                                    call find_connected_dets_hubbard(det_i_up, det_i_dn, n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements)
                                endif
                                do k=1,n_connected_dets
                                    call binary_search_single(connected_dets_up(k),tmp_all_dets_up,ind_1)
                                    call binary_search_single(connected_dets_dn(k),tmp_all_dets_dn,ind_2)
                                    ind=(n_det_dn*(ind_1-1))+ind_2
                                    w(ind)=w(ind)+(connected_matrix_elements(k)*v(i,vnum))
                                    if (ind==i) w(ind)=w(ind)-(shift*v(i,vnum))
                                enddo
                            endif
                        enddo
                    elseif (hamiltonian_type .eq. 'hubbardk') then
                        do i=1,n
                            det_i_up=all_dets_up(i)
                            det_i_dn=all_dets_dn(i)
                            if (abs(v(i,vnum)) .gt. 1.0e-15) then
                                call find_connected_dets_hubbard_k(det_i_up,det_i_dn,n_connected_dets,connected_dets_up,connected_dets_dn,nsites,connected_matrix_elements)
                                do k=1,n_connected_dets
                                    call binary_search_single(connected_dets_up(k),tmp_all_dets_up,ind_1)
                                    call binary_search_single(connected_dets_dn(k),tmp_all_dets_dn,ind_2)
                                    ind=(n_det_dn*(ind_1-1))+ind_2
                                    ind=inverse_k_map(ind)
                                    w(ind)=w(ind)+(connected_matrix_elements(k)*v(i,vnum))
                                    if (ind==i) w(ind)=w(ind)-(shift*v(i,vnum))
                                enddo
                            endif
                        enddo
                    endif

                    if (it .gt. npower) then
                        write (6,*) "vnum=",vnum
                        call flush(6)
                        do i=1,vnum
                            full_hessen(i,vnum)=dot_product(v(:,i),w(:))
                            call flush(6)
                            w(:)=w(:)-full_hessen(i,vnum)*v(:,i)
                        enddo
                        v(:,vnum+1)=w(:)
                        w(:)=0._rk
                        norm=dot_product(v(:,vnum+1),v(:,vnum+1))
                        full_hessen(vnum+1,vnum)=norm**(0.5_rk)
                        norm_inv=1._rk/(norm**(0.5_rk))
                        v(:,vnum+1)=v(:,vnum+1)*norm_inv

                        if (allocated(hessen)) deallocate(hessen)
                        allocate(hessen(it-npower,it-npower))
                        hessen(:,:)=0._rk

                        if (allocated(eigenvalues))    deallocate(eigenvalues)
                        if (allocated(im_eigenvalues)) deallocate(im_eigenvalues)
                        allocate(eigenvalues(it-npower))
                        allocate(im_eigenvalues(it-npower))
                        eigenvalues(:)=0._rk
                        im_eigenvalues(:)=0._rk

                        hessen(1:it-npower,1:it-npower)=full_hessen(1:it-npower,1:it-npower)

                        write (6,*) "Hessenberg matrix"
                        call print_real_matrix(it-npower,it-npower,hessen)
                        call real_general_diagonalize_ow_ham(it-npower,hessen,eigenvalues,im_eigenvalues)
                        write (6,*) "Eigenvalues"
                        do c=1,it-npower
                            write(6,*) eigenvalues(c),im_eigenvalues(c)
                        enddo

                        lowest_eigenvalue=eigenvalues(1)
                        !call print_real_matrix(size(eigenvalues),1,eigenvalues)
                        if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<1.0e-10) then
                            converged=.true.
                            exit
                        else
                            lowest_eigenvalue_prev=lowest_eigenvalue
                            write(6,'(''Iteration, Eigenvalue='',i3,f15.9)') it, lowest_eigenvalue+shift
                            call flush(6)

                            !if (it-npower>2) shift=-(eigenvalues(2)+eigenvalues(it-npower))/2._rk
                            if (it<iterations) then
                                deallocate(eigenvalues)
                                deallocate(im_eigenvalues)
                                deallocate(hessen)
                            endif
                         endif
                 else   ! Power method
                        v(:,1)=w(:)
                        w(:)=0._rk
                        norm=dot_product(v(:,1),v(:,1))
                        norm_inv=1._rk/(norm**(0.5_rk))
                        v(:,1)=v(:,1)*norm_inv
                 endif
               enddo

                 it=min(it,iterations)
                 lowest_eigenvalue=lowest_eigenvalue+shift
                 write(6,'(''n, Lowest eigenvalue ='',i10, f16.10)') n, lowest_eigenvalue
                 call flush(6)

                 if (degen .eq. 1) lowest_eigenvec_1(:)=matmul(v(:,1:it-npower),hessen(:,1))
                 if (converged) then
                    exit
                 else
                    if (degen .eq. 1) v(:,1)=lowest_eigenvec_1(:)
                 endif
                 deallocate(eigenvalues)
                 deallocate(im_eigenvalues)
                 deallocate(hessen)

             enddo ! cycles

             if (degen .eq. 1) then
                 lowest_eigenvalue_1=lowest_eigenvalue
                 norm=0._rk
                 do i=1,n
                    v(i,1)=rannyu()
                 enddo
                 norm=dot_product(v(:,1),lowest_eigenvec_1)
                 v(:,1)=v(:,1)-norm*lowest_eigenvec_1(:)
                 norm=dot_product(v(:,1),v(:,1))
                 norm_inv=1._rk/sqrt(norm)
                 v(:,1)=v(:,1)*norm_inv
                 converged=.false.                               ! Reset convergence
            endif
         enddo ! degeneracy

         lowest_eigenvec(:)=lowest_eigenvec_1(:)
         lowest_eigenvalue=lowest_eigenvalue_1

         ! Get the mixed energy for hubbard2
         if (hamiltonian_type .eq. 'hubbard2') then
             w(:)=0._rk
             do i=1,n
                if (abs(stored_guiding_wf(i)) .gt. 1.0e-15) lowest_eigenvec(i)=lowest_eigenvec(i)/stored_guiding_wf(i)                           ! Scale eigenvector back by guiding wavefunction
                ind_1=((i-1)/n_det_dn)+1
                ind_2=i-(n_det_dn*(ind_1-1))
                det_i_up=tmp_all_dets_up(ind_1)
                det_i_dn=tmp_all_dets_dn(ind_2)
                if (abs(stored_guiding_wf(i)) .gt. 1.0e-15) then
                    call find_connected_dets_hubbard(det_i_up, det_i_dn, n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements)
                    do k=1,n_connected_dets
                        call binary_search_single(connected_dets_up(k),tmp_all_dets_up,ind_1)
                        call binary_search_single(connected_dets_dn(k),tmp_all_dets_dn,ind_2)
                        ind=(n_det_dn*(ind_1-1))+ind_2
                        w(ind)=w(ind)+(connected_matrix_elements(k)*stored_guiding_wf(i))           ! H acting on guiding = trial wavefunction
                    enddo
                endif
            enddo

            ! W= H acting on stored_wf
            e_mix_numerator=dot_product(w,lowest_eigenvec)
            e_mix_denominator=dot_product(stored_guiding_wf,lowest_eigenvec)
            e_mix=e_mix_numerator/e_mix_denominator
            write(6,*) "Mixed energy numerator",e_mix_numerator
            write(6,*) "Mixed energy denominator",e_mix_denominator
            write(6,*) "Mixed energy estimator",e_mix
         endif

         deallocate(w)
         deallocate(v)

         call my_second(2,'arnoldi_hubbard_binary_search') ; call flush(6)

   end subroutine arnoldi_hubbard_binary_search

!--------------------------------------------------------------------------------------------------

   subroutine lanczos_hubbard_binary_search(lowest_eigenvec,lowest_eigenvalue,all_dets_up,all_dets_dn)
   !------------------------------------------------------------------------------------
   ! Description   : Perform Lanczos for Hubbard Hamiltonian 2D lattice given filling
   !                 in real or momentum space
   !                 Specify lattice dimensions l_x,l_y and boundary condition pbc
   !                 Also specify number f alpha and beta electrons n_alpha and n_beta
   !                 Hubbard parameters t (hopping),U (Coloumb repulsion)needed
   ! Author        : H.J. Changlani
   ! ------------------------------------------------------------------------------------
        implicit none

        real(rk),allocatable,intent(out)     ::  lowest_eigenvec(:)
        real(rk),intent(out)                 ::  lowest_eigenvalue
#ifdef NUM_ORBITALS_GT_127
        type(ik_vec), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        type(ik_vec), allocatable             ::  tmp_all_dets_up(:),tmp_all_dets_dn(:)
        type(ik_vec)                          ::  rep_up,rep_dn
        type(ik_vec)                          ::  det_i_up,det_i_dn
        type(ik_vec)                          ::  tmp_up,tmp_dn
#else
        integer(ik), allocatable,intent(out) ::  all_dets_up(:),all_dets_dn(:)
        integer(ik), allocatable             ::  tmp_all_dets_up(:),tmp_all_dets_dn(:)
        integer(ik)                          ::  rep_up,rep_dn
        integer(ik)                          ::  det_i_up,det_i_dn
        integer(ik)                          ::  tmp_up,tmp_dn
#endif

        !Local variables
        integer(ik),allocatable              ::  inverse_k_map(:),order(:)
        integer                              ::  n_connected_dets
        integer                              ::  i,j,k,n,ncycles,ind_1,ind_2,it,vnum,istat,degen
        integer(ik)                           ::  ind
        integer                              ::  ktmp(2)
        integer                              ::  num_degen,iterations=50,num_cycles=25,npower=0
        real(rk)                             ::  norm,norm_inv,shift,rannyu
        real(rk),allocatable                 ::  w(:),alphas(:),betas(:),v(:,:)
        real(rk),allocatable                 ::  tri(:,:),eigenvalues(:)
        real(rk)                             ::  energy,lowest_eigenvalue_prev,lowest_eigenvalue_1,lowest_eigenvalue_2
        real(rk),allocatable                 ::  lowest_eigenvec_1(:),lowest_eigenvec_2(:)
        logical,allocatable                  ::  allowed_det(:)
        logical                              ::  converged=.false.
        integer                              ::  nnzero

        if (npower>iterations) npower=iterations -1              !A check

        shift=U*real(min(nup,ndn))

        n_det_up=int(n_choose_k(l_x*l_y,nup),i4b)
        n_det_dn=int(n_choose_k(l_x*l_y,ndn),i4b)
        num_degen=1
        call flush(6)
        iterations=min(iterations,n_det_up*n_det_dn)

        allocate(tmp_all_dets_up(n_det_up),stat=istat)
        allocate(tmp_all_dets_dn(n_det_dn),stat=istat)
        allocate(alphas(iterations+1),stat=istat)
        allocate(betas(iterations+1),stat=istat)
        allocate(allowed_det(n_det_up*n_det_dn),stat=istat)
        allocate(inverse_k_map(n_det_up*n_det_dn),stat=istat)
        allocate(order(n_det_up*n_det_dn),stat=istat)
        allocate(tri(1,1),stat=istat)                             ! To avoid warning messages since allocated under if condition
        allocate(eigenvalues(1),stat=istat)                       ! To avoid warning messages since allocated under if condition

        call constrained_dets(nup,n_det_up,tmp_all_dets_up)
        call constrained_dets(ndn,n_det_dn,tmp_all_dets_dn)

        if (hamiltonian_type .eq. 'hubbard2') then
            n=n_det_up*n_det_dn
            write (6,*) "In Binary search Lanczos: The total number of determinants in real space is",n
        elseif (hamiltonian_type .eq. 'hubbardk') then
            !n = number of combinations of up,dn dets with total momentum equal to ktot
            n=0
            do j=1,n_det_up*n_det_dn
              ind_1=((j-1)/n_det_dn)+1
              ind_2=j-(n_det_dn*(ind_1-1))
              ktmp(:) = 0
              do i=1,nsites
                if (btest(tmp_all_dets_up(ind_1),i-1)) ktmp = ktmp + k_vectors(:,i)
                if (btest(tmp_all_dets_dn(ind_2),i-1)) ktmp = ktmp + k_vectors(:,i)
              enddo
              if ((mod(ktmp(1)-ktot(1),2*l_x)==0) .and. (mod(ktmp(2)-ktot(2),2*l_y)==0)) then
                if (space_sym) then
                    call get_rep_only(c4_map,reflection_map,z,p,tmp_all_dets_up(ind_1),tmp_all_dets_dn(ind_2),rep_up,rep_dn)
                    if ((tmp_all_dets_up(ind_1) .eq. rep_up) .and. (tmp_all_dets_dn(ind_2) .eq. rep_dn)) then
                        n = n + 1
                        order(n) = j
                        allowed_det(j)=.true.
                    else
                        allowed_det(j)=.false.
                    endif
                else
                    n = n + 1
                    order(n) = j
                    allowed_det(j)=.true.
                endif
              else
                allowed_det(j)=.false.
              endif
            enddo
            write (6,*) "In Binary search Lanczos: The total number of determinants in momentum space is",n
            call flush(6)
        endif

        allocate(lowest_eigenvec(n),stat=istat)
        allocate(lowest_eigenvec_1(n),stat=istat)
        allocate(lowest_eigenvec_2(n),stat=istat)
        allocate(w(n),stat=istat)
        allocate(v(n,iterations+1-npower),stat=istat)
        allocate(all_dets_up(n),stat=istat)
        allocate(all_dets_dn(n),stat=istat)

        v(:,:)=0._rk
        w(:)=0._rk

        if (hamiltonian_type .eq. 'hubbardk') then
            k=0
            do i=1,n_det_up*n_det_dn
                if (allowed_det(i) .eqv. .true.) then
                    k=k+1
                    ind_1=((i-1)/n_det_dn)+1
                    ind_2=i-(n_det_dn*(ind_1-1))
                    all_dets_up(k)=tmp_all_dets_up(ind_1)
                    all_dets_dn(k)=tmp_all_dets_dn(ind_2)
                    inverse_k_map(i)=k
                endif
            enddo
        endif

        if (hamiltonian_type .eq. 'hubbard2') then
            do i=1,n_det_up*n_det_dn
                    ind_1=((i-1)/n_det_dn)+1
                    ind_2=i-(n_det_dn*(ind_1-1))
                    all_dets_up(i)=tmp_all_dets_up(ind_1)
                    all_dets_dn(i)=tmp_all_dets_dn(ind_2)
            enddo
        endif

        norm=0._rk
        do i=1,n
               det_i_up=all_dets_up(i)
               det_i_dn=all_dets_dn(i)
                   call flush(6)
                   if ((hamiltonian_type .eq. 'hubbard2') .and. (main_wf_type .eq. 'gutz')) then ! variational wavefunction as start guess
                       call wf_calc(det_i_up,det_i_dn,v(i,1))
                   elseif (hamiltonian_type .eq. 'hubbardk') then
                     if (space_sym) then
                        tmp_up = det_i_up
                        tmp_dn = det_i_dn
                        call hamiltonian_hubbard_k_space_sym(tmp_up,tmp_dn,det_i_up,det_i_dn,energy,nnzero)
                     else
                        call hamiltonian_hubbard_k(det_i_up,det_i_dn,det_i_up,det_i_dn,energy)
                    endif
                     if (abs(energy-k_hf_energy)<1.0e-10) then                      ! If the determinant is another equal energy state
                         write(6,'(g14.6)') energy
                         v(i,1)=rannyu()
                     endif
                   else
                       v(i,1)=rannyu()                                              ! Else choose random vector as starting guess
                   endif
               norm=norm+v(i,1)*v(i,1)
        enddo

        if (abs(norm) .le. 1.0e-10) then
            write (6,*) "Norm of start vector should never be so small. Most likely your k space code has an inconsistency"
            write (6,*) "I am still proceeding, but dont expect the correct result"
            do i=1,n
                v(i,1)=rannyu()
            enddo
            norm=dot_product(v(:,1),v(:,1))
        endif

        norm_inv=1._rk/(norm**(0.5_rk))                                              ! Normalize v_1
        v(:,1)=v(:,1)*norm_inv

        betas(1)=0._rk

        write(6,'(/,''Diagonalizing in Lanczos_hubbard (with binary search)'')')
        do degen=1,num_degen
            write(6,'(/,''Computing state number '',i2)') degen
            do ncycles=1,num_cycles
                write (6,*)
                write(6,'(/,''In Cycle '',i2)') ncycles
                do it=1,iterations
                    call flush(6)
                    if (it .gt. npower) then
                        vnum=it-npower
                        !shift=0._rk
                    else
                        vnum=1
                    endif
                    call flush(6)
                    !Computing H*v_p - This is the bulk of the operation
                    if ((hamiltonian_type .eq. 'hubbard2') .or. (ham_cross .eqv. .true.) .or. (ham_square .eqv. .true.)) then
                         !write (6,*) "Ham_cross",ham_cross
                         do i=1,n
                            ind_1=((i-1)/n_det_dn)+1
                            ind_2=i-(n_det_dn*(ind_1-1))
                            det_i_up=tmp_all_dets_up(ind_1)
                            det_i_dn=tmp_all_dets_dn(ind_2)
                            if (abs(v(i,vnum)) .gt. 1.0e-15) then
                                call find_connected_dets_hubbard(det_i_up,det_i_dn,n_connected_dets,connected_dets_up,connected_dets_dn,connected_matrix_elements)
                                do k=1,n_connected_dets
                                    call binary_search_single(connected_dets_up(k),tmp_all_dets_up,ind_1)
                                    call binary_search_single(connected_dets_dn(k),tmp_all_dets_dn,ind_2)
                                    ind=(n_det_dn*(ind_1-1))+ind_2
                                    w(ind)=w(ind)+(connected_matrix_elements(k)*v(i,vnum))
                                    if (ind==i) w(ind)=w(ind)-(shift*v(i,vnum))
                                enddo
                            endif
                        enddo
                    elseif (hamiltonian_type .eq. 'hubbardk') then
                        do i=1,n
                            det_i_up=all_dets_up(i)
                            det_i_dn=all_dets_dn(i)
                            if (abs(v(i,vnum)) .gt. 1.0e-15) then
                                call find_connected_dets_hubbard_k(det_i_up,det_i_dn,n_connected_dets,connected_dets_up,connected_dets_dn,nsites,connected_matrix_elements)
                                do k=1,n_connected_dets
                                    call binary_search_single(connected_dets_up(k),tmp_all_dets_up,ind_1)
                                    call binary_search_single(connected_dets_dn(k),tmp_all_dets_dn,ind_2)
                                    ind=(n_det_dn*(ind_1-1))+ind_2
                                    ind=inverse_k_map(ind)
                                    w(ind)=w(ind)+(connected_matrix_elements(k)*v(i,vnum))
                                    if (ind==i) w(ind)=w(ind)-(shift*v(i,vnum))
                                enddo
                            endif
                        enddo
                    endif

                    if (it .gt. npower) then
                        if (it .gt. npower+1) w(:)=w(:)-betas(vnum)*v(:,vnum-1)
                        alphas(vnum)=dot_product(w,v(:,vnum))
                        w(:)=w(:)-alphas(vnum)*v(:,vnum)
                        norm=dot_product(w,w)
                        betas(vnum+1)=norm**(0.5_rk)
                        norm_inv=1._rk/betas(vnum+1)
                        v(:,vnum+1)=w(:)*norm_inv
                        w(:)=v(:,vnum+1)
                        do i=1,it-npower
                            norm=dot_product(v(:,vnum+1),v(:,i))
                            !write (6,*) "q=",norm
                            call flush(6)
                            w(:)=w(:)-norm*v(:,i)
                        enddo
                        if (degen .eq. 2) then
                            norm=dot_product(v(:,vnum+1),lowest_eigenvec_1(:))
                            w(:)=w(:)-norm*lowest_eigenvec_1(:)
                        endif
                        v(:,vnum+1)=w(:)
                        w(:)=0._rk
                        norm=dot_product(v(:,vnum+1),v(:,vnum+1))
                        !write (6,*) "Norm after reorth = ",norm
                        norm_inv=1._rk/(norm**(0.5_rk))
                        v(:,vnum+1)=v(:,vnum+1)*norm_inv

                        if (allocated(tri)) deallocate(tri)
                        allocate(tri(it-npower,it-npower))
                        tri(:,:)=0._rk

                        if (allocated(eigenvalues)) deallocate(eigenvalues)
                        allocate(eigenvalues(it-npower))
                        eigenvalues(:)=0._rk

                        do i=1,it-npower
                           do j=1,it-npower
                            tri(i,j)=0._rk
                            enddo
                        enddo

                        do i=1,it-npower
                            tri(i,i)=alphas(i)
                            if (i<it-npower) then
                                tri(i,i+1)=betas(i+1)
                                tri(i+1,i)=betas(i+1)
                            endif
                        enddo

                        call real_symmetric_diagonalize_ow_ham(it-npower,tri,eigenvalues)

                        lowest_eigenvalue=eigenvalues(1)
                        !call print_real_matrix(size(eigenvalues),1,eigenvalues)
                        if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<1.0e-10) then
                            converged=.true.
                            exit
                        else
                            lowest_eigenvalue_prev=lowest_eigenvalue
                            write(6,'(''Iteration, Eigenvalue='',i3,f15.9)') it, lowest_eigenvalue+shift
                            call flush(6)

                            !if (it-npower>2) shift=-(eigenvalues(2)+eigenvalues(it-npower))/2._rk
                            if (it<iterations) then
                                deallocate(eigenvalues)
                                deallocate(tri)
                            endif
                         endif
                 else
                        v(:,1)=w(:)
                        w(:)=0._rk
                        norm=dot_product(v(:,1),v(:,1))
                        norm_inv=1._rk/(norm**(0.5_rk))
                        v(:,1)=v(:,1)*norm_inv
                 endif
               enddo

                 it=min(it,iterations)
                 lowest_eigenvalue=lowest_eigenvalue+shift
                 write(6,'(''n, Lowest eigenvalue ='',i10, f16.10)') n, lowest_eigenvalue
                 call flush(6)

                 !write(6,*) "iterations",it
                 !write (6,*) "tridiagonal matrix rows",size(tri,1)
                 !write (6,*) "tridiagonal matrix cols",size(tri,2)
                 !write (6,*) v
                 call flush(6)

                 if (degen .eq. 1) lowest_eigenvec_1(:)=matmul(v(:,1:it-npower),tri(:,1))
                 if (degen .eq. 2) lowest_eigenvec_2(:)=matmul(v(:,1:it-npower),tri(:,1))

                 if (converged) then
                    exit
                 else
                    if (degen .eq. 1) v(:,1)=lowest_eigenvec_1(:)
                    if (degen .eq. 2) v(:,1)=lowest_eigenvec_2(:)
                 endif
                 deallocate(eigenvalues)
                 deallocate(tri)

             enddo ! cycles

             if (degen .eq. 1) then
                 lowest_eigenvalue_1=lowest_eigenvalue
                 norm=0._rk
                 do i=1,n
                    v(i,1)=rannyu()
                 enddo
                 norm=dot_product(v(:,1),lowest_eigenvec_1)
                 v(:,1)=v(:,1)-norm*lowest_eigenvec_1(:)
                 norm=dot_product(v(:,1),v(:,1))
                 norm_inv=1._rk/sqrt(norm)
                 v(:,1)=v(:,1)*norm_inv
                 converged=.false.                               ! Reset convergence
            endif

            if (degen .eq. 2) lowest_eigenvalue_2=lowest_eigenvalue
         enddo ! degeneracy

         lowest_eigenvec(:)=lowest_eigenvec_1(:)
         lowest_eigenvalue=lowest_eigenvalue_1

         deallocate(w)
         deallocate(v)
         deallocate(alphas)
         deallocate(betas)

         call my_second(2,'lanczos_hubbard_binary_search') ; call flush(6)

   end subroutine lanczos_hubbard_binary_search


!===========================================================================
  subroutine test_k_space_code()
!===========================================================================

  integer                     :: i,j,k,n,m
  integer                     :: n_connected_dets
  real(rk),allocatable        :: ham(:,:)!,matrix_elements(:)
  real(rk),allocatable        :: ham2(:,:)!,matrix_elements2(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable     :: dets_up(:),dets_dn(:)
  type(ik_vec)                 :: det_up(1),det_dn(1)
  type(ik_vec),allocatable     :: con_up(:),con_dn(:)
#else
  integer(ik),allocatable     :: dets_up(:),dets_dn(:)
  integer(ik)                 :: det_up(1),det_dn(1)
  integer(ik),allocatable     :: con_up(:),con_dn(:)
#endif
  real(rk),allocatable        :: mat(:)
  real(rk),allocatable        :: eigenvectors(:,:),eigenvalues(:)
  real(rk), allocatable       :: gs_eigenvec(:)
  real(rk)                    :: gs_energy

  call generate_k_vectors()

  write (6,*) "K map:"
  call print_int_matrix(nsites,2*l_x-1,kmap(:,:,1))
  call flush(6)

  call make_hubbard_hamiltonian_k(ktot,ham,dets_up,dets_dn,gs_eigenvec,gs_energy)

  m=size(ham,1)
  allocate(ham2(m,m))

  ham2(:,:)=0

  allocate(con_up(nup*ndn*(nsites-nup)))
  allocate(con_dn(nup*ndn*(nsites-nup)))
  allocate(mat(nup*ndn*(nsites-nup)))

  do i=1,m
    call find_connected_dets_hubbard_k(dets_up(i),dets_dn(i),n,con_up,con_dn,nsites,mat)
    do j=1,n
      do k=1,m
        if (con_up(j)==dets_up(k) .and. con_dn(j)==dets_dn(k))  ham2(i,k)=mat(j)
      enddo
    enddo
  enddo
  allocate(eigenvalues(m))
  allocate(eigenvectors(m,m))

  write (6,*) "Dimension of Hamiltonian:",m
  flush(6)

  call real_symmetric_diagonalize(m,ham2,eigenvectors,eigenvalues)

  if (m<21) then
    write(6,*) "Hamiltonian"
    call print_real_matrix(m,m,ham2)
    write (6,*) "Eigenvalues"
    call print_real_matrix(size(eigenvalues),1,eigenvalues)
    write (6,*)
    write (6,*)
  else
    write (6,'(''Exact Ground State Energy:'',f10.6)') eigenvalues(1)
  endif

  deallocate(eigenvalues)
  deallocate(eigenvectors)

  if (size(dets_up)<10) then
    write (6,*) "Basis states:"
    write (6,*) "Up:"
    call print_bin(size(dets_up),nsites,dets_up)
    write (6,*) "Down:"
    call print_bin(size(dets_dn),nsites,dets_dn)
  endif

  if (nsites<20) then
    write (6,*) "k-vectors:"
    call print_int_matrix(2,nsites,k_vectors)
    write (6,*) "k-energies:"
    call print_real_matrix(1,nsites,k_energies)
  endif

  det_up(1) = dets_up(1)
  det_dn(1) = dets_dn(1)

 !allocate(connected_dets_up(nup*ndn*(nsites-nup)))
 !allocate(connected_dets_dn(nup*ndn*(nsites-nup)))
 !allocate(matrix_elements(nup*ndn*(nsites-nup)))

  call find_connected_dets_hubbard_k(det_up(1),det_dn(1),n_connected_dets,connected_dets_up,connected_dets_dn,nsites,connected_matrix_elements)

  write (6,*) "Up Det:"
  call print_bin(1,nsites,det_up)

  write (6,*) "Down Det:"
  call print_bin(1,nsites,det_dn)

  write (6,*) "Number of Connected Dets:",n_connected_dets

  write (6,*) "Connected Up Dets:"
  call print_bin(n_connected_dets,nsites,connected_dets_up)

  write (6,*) "Connected Down Dets:"
  call print_bin(n_connected_dets,nsites,connected_dets_dn)

  write (6,*) "Matrix Elements:"
  call print_real_matrix(1,n_connected_dets,connected_matrix_elements)

  call flush(6)

  endsubroutine test_k_space_code
  !===========================================================================

  !===========================================================================
  subroutine test_hamiltonian_hubbard_dm()
    !---------------------------------------------------------------------------
    ! Description : Test hamiltonian construction of 4,4 1,1 in DM code
    ! Created     : H. J. Changlani Jan 18 2012
    !---------------------------------------------------------------------------
    implicit none

    integer                    :: i,i1,j,j1,k,l,n,num,num_blocks,total_nup_poss,total_ndn_poss
    real(rk),allocatable       :: ham(:,:),eigenvecs(:,:),eigenvalues(:),lowest_eigenvector(:)
    integer                    :: n_det_up_sav,n_det_dn_sav,nup_sav,ndn_sav,nsites_sav
    integer                    :: l_x_sav,l_y_sav
    logical                    :: pbc_sav
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),allocatable    :: all_dm_dets_up(:),all_dm_dets_dn(:)
#else
    integer(ik),allocatable    :: all_dm_dets_up(:),all_dm_dets_dn(:)
#endif
    integer                    :: nup_poss(625,4),ndn_poss(625,4)
    integer                    :: counter
    integer,allocatable        :: combined_nup_ndn_poss(:,:)
    integer                    :: shift,total,repeat_frequency
    integer                    :: combined_ind,how_many_times
    integer                    :: beta_basis_start_end(625*625,2)
    real(rk)                   :: lowest_eigenvalue

    n_det_up_sav=n_det_up
    n_det_dn_sav=n_det_dn
    nsites_sav=nsites
    nup_sav=nup
    ndn_sav=ndn
    l_x_sav=l_x
    l_y_sav=l_y
    pbc_sav=pbc

    l_x=4
    l_y=4
    pbc=.true.
    nsites=16
    nup=1
    ndn=1

     num_blocks=4
     counter=1
     ! Make nup possibility list on blocks 1 2 3 4
     do i=0,624
          num=i
          do k=1,4
            nup_poss(counter,k)=mod(num,5)
            num=num/5
          enddo
          if (sum(nup_poss(counter,:)) .eq. nup) then
            counter=counter+1
          else
            nup_poss(counter,:)=0
          endif
     enddo
     total_nup_poss=counter-1

     ! Make ndn possibility list on blocks 1 2 3 4
     counter=1
     do i=0,624
          num=i
          do k=1,4
            ndn_poss(counter,k)=mod(num,5)
            num=num/5
          enddo
          if (sum(ndn_poss(counter,:)) .eq. ndn) then
            counter=counter+1
          else
            ndn_poss(counter,:)=0
          endif
     enddo
     total_ndn_poss=counter-1

     write (6,*) "Nup possibilities on the blocks"
     do i=1,total_nup_poss
                 write(6,'(4i5)') (nup_poss(i,j),j=1,4)
     enddo
     call flush(6)

     write (6,*) "Ndn possibilities on the blocks"
     do i=1,total_ndn_poss
             write(6,'(4i5)') (ndn_poss(i,j),j=1,4)
     enddo
     call flush(6)

     allocate(combined_nup_ndn_poss(total_nup_poss*total_ndn_poss,8))

     write(6,'(/,'' Generating all nup,ndn combinations on blocks.......'')')
     call flush(6)

     ! Combine all possibilities on blocks 1 up 2up 3up 4up 1dn 2dn 3dn 4dn
     counter=1
     do i=1,total_nup_poss
                do j=1,total_ndn_poss
                    do i1=1,4
                        combined_nup_ndn_poss(counter,i1)=nup_poss(i,i1)
                        combined_nup_ndn_poss(counter,i1+4)=ndn_poss(j,i1)
                    enddo
                    counter=counter+1
                enddo
     enddo

     ! Expand maps for all_dets_dm
     counter=1
     shift=0

     write(6,'(/,'' Generating full DM basis.......'')')
     call flush(6)

     n=int(n_choose_k(nsites,nup)*n_choose_k(nsites,ndn),i4b)
     allocate(ham(n,n))
     allocate(eigenvecs(n,n))
     allocate(eigenvalues(n))
     allocate(all_dm_dets_up(n))
     allocate(all_dm_dets_dn(n))
     all_dm_dets_up(:)=0_ik
     all_dm_dets_dn(:)=0_ik

     do i=1,total_ndn_poss*total_nup_poss
        !write (6,*) "i (combined_nup_ndn)=",i
        call flush(6)
        total=1
        do j=1,4 ! which block
            total=total*(-ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),1)+ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),2)+1)
        enddo
        !write (6,*) "total=", total
        call flush(6)
        ! total possibilities with these numbers on all blocks
        combined_ind=0
        do j=1,4
            combined_ind=combined_ind+(5**(j-1))*(combined_nup_ndn_poss(i,4+j))+(5**(4+j-1))*combined_nup_ndn_poss(i,j)
        enddo
        beta_basis_start_end(combined_ind,1)=shift+1
        beta_basis_start_end(combined_ind,2)=shift+total

        how_many_times=1
        do j=1,4    ! over blocks
            repeat_frequency=total/(ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),2)-ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),1)+1)
            repeat_frequency=repeat_frequency/how_many_times
            counter=1
            !write (6,*) "repeat frequency=",repeat_frequency
            !write (6,*) "how_many_times=",how_many_times
            call flush(6)
            do j1=1,repeat_frequency
                do k=ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),1),ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),2)
                  ! k is a number between 1 and 256 indicating which dm eigenvector we are looking at
                  do l=1,how_many_times
                      do i1=1,8
                        if (i1 .gt. 4 .and. btest(k,i1-1)) all_dm_dets_up(counter+shift)=ibset(all_dm_dets_up(counter+shift),4*(j-1)+i1-5)
                        if (i1 .le. 4 .and. btest(k,i1-1)) all_dm_dets_dn(counter+shift)=ibset(all_dm_dets_dn(counter+shift),4*(j-1)+i1-1)
                      enddo
                      counter=counter+1
                  enddo
                enddo
             enddo
             how_many_times=how_many_times*(ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),2)-ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),1)+1)
        enddo
        !write (6,*) "shift=",shift
        shift=shift+total
     enddo

     write(6,'(/,'' Finished full DM basis.......'')')
     call flush(6)

    write(6,'(/,'' Computing full Hamiltonian matrix.......'')')
    call flush(6)
    ham(:,:)=0._rk
    do i=1,n
       do j=1,n
         call hamiltonian_hubbard_dm(all_dm_dets_up(i),all_dm_dets_dn(i),all_dm_dets_up(j),all_dm_dets_dn(j),ham(i,j))
       enddo
    enddo
    write(6,'(/,'' Finished computing full Hamiltonian matrix.......'')')
    call flush(6)

    write(6,'(/,'' Verifying symmetry of Hamiltonian.......'')')
    call flush(6)

    ! Check if hamiltonian is symmetric
    do i=1,n
        do j=1,n
            if (abs(ham(i,j)-ham(j,i))>1.0e-10) then
                write (6,*) "There is a bug because Hamiltonian is not symmetric"
                call print_real_matrix(n,n,ham)
                call flush(6)
                stop
            endif
        enddo
    enddo

    write(6,'(/,'' Hamiltonian passed symmetry test.......'')')
    write (6,*)
    call flush(6)

    if (n .le. 2000) then
        write (6,*) "Hamiltonian"
        call print_real_matrix(n,n,ham,20,20)
        write (6,*) "Testing",l_x,",",l_y,",",nup,",",ndn," Hamiltonian in DM basis with full diagonalization"
        call real_symmetric_diagonalize(n,ham,eigenvecs,eigenvalues)
        write (6,*) "All Eigenvalues"
        call print_real_matrix(size(eigenvalues),1,eigenvalues)
        call flush(6)
    endif

    allocate(lowest_eigenvector(n))
    call matrix_lanczos(n,lowest_eigenvector,lowest_eigenvalue,ham)

    write (6,*) "Testing",l_x,",",l_y,",",nup,",",ndn," Hamiltonian in DM basis with matrix Lanczos"
    write (6,*) ""
    call flush(6)
    write (6,*) "Final eigenvector in DM basis"
    call sort(lowest_eigenvector,all_dm_dets_up,all_dm_dets_dn)
    call print_real_matrix(size(lowest_eigenvector),1,lowest_eigenvector)
    call flush(6)

    n_det_up=n_det_up_sav
    n_det_dn=n_det_dn_sav
    nsites=nsites_sav
    nup=nup_sav
    ndn=ndn_sav
    l_x=l_x_sav
    l_y=l_y_sav
    pbc=pbc_sav

    !stop

  end subroutine test_hamiltonian_hubbard_dm
  !===========================================================================


!===========================================================================
  subroutine sort_electrons(n,det_old,det_new,inv)
    !---------------------------------------------------------------------------
    ! Description : Given bitpacked strings det_old,det_new, shuffles columns
    !               of inv so that all n electrons are sorted accoring to
    !               site number.
    !
    ! Created     : A. Holmes, 2 Jun 2011
    !---------------------------------------------------------------------------

  implicit none

  integer,intent(in) :: n
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(in) :: det_old,det_new
#else
  integer(ik), intent(in) :: det_old,det_new
#endif
  real(rk),dimension(n,n),intent(out) :: inv

  real(rk),dimension(n,n) :: tmp
  integer,dimension(n)    :: order
  integer :: i,j,k,l

  l=0
  k=1
  do i=1,nsites
    if ((btest(det_old,i-1) .eqv. .true.) .and. (btest(det_new,i-1) .eqv. .true.)) then ! electron didn't move
      l = l+1
      order(l) = k
      k = k+1
    endif
    if ((btest(det_old,i-1) .eqv. .true.) .and. (btest(det_new,i-1) .eqv. .false.)) then ! electron moved from this site
      l = l+1
      order(l) = -1
    endif
    if ((btest(det_old,i-1) .eqv. .false.) .and. (btest(det_new,i-1) .eqv. .true.)) then ! electron moved to this site
      j=k
      k=k+1
    endif
  enddo
  do i=1,n
    if (order(i) .eq. -1) then
      order(i) = j
    endif
  enddo

  do i=1,n
    tmp(:,order(i)) = inv(:,i)
  enddo

  inv = tmp

  end subroutine sort_electrons

!===============================================================================================
  subroutine density_matrix_2by2(cluster)
    !-------------------------------------------------------------------------------------------
    ! Description : Compute the density matrix of a 2 by 2 square embedded in a 12 site lattice
    ! Purpose     : At the moment, to get some intuition on how the DM looks!
    ! Created     : H.J. Changlani Oct 6 2011
    !-------------------------------------------------------------------------------------------

!                              6----7
!                              |    |
!                          5---1----2---8
!                          |   |    |   |
!                          9---3----4---10
!                              |    |
!                              11---12

!                                or
!
!                          13---6----7---14
!                           |   |    |   |
!                           5---1----2---8
!                           |   |    |   |
!                           9---3----4---10
!                           |   |    |   |
!                           15--11---12--16


!  where the DM computed is that of the cluster 1234 in the ground state of the cluster
!  We have tried to be intelligent about the numbering so that we dont have to worry about fermion signs when computing the Density Matrix
   character*16,intent(in)     :: cluster
   integer                     :: i,j,k,l,shift,i_up,i_dn,j_up,j_dn,num
   integer(ik)                  :: i16,j16,nup_sys,ndn_sys,nup_sys2,ndn_sys2,num_int16
   integer(ik)                 :: env_up,env_dn
   integer                     :: nup_env,ndn_env
   integer                     :: ind1_int4,ind2_int4
   integer(ik)                 :: ind1,ind2
   real(rk)                    :: density_matrix(256,256)
   real(rk)                    :: energy_cross
   real(rk),allocatable        :: lowest_eigenvec_cross(:)
   real(rk),allocatable        :: lowest_evec_2x2(:),tmp_ham_2x2(:,:)
   real(rk)                    :: lowest_eval_2x2
#ifdef NUM_ORBITALS_GT_127
   type(ik_vec),allocatable     :: all_dets_up_cross(:),all_dets_dn_cross(:)
   type(ik_vec),allocatable     :: all_dets_up_2x2(:),all_dets_dn_2x2(:)
   type(ik_vec),allocatable    :: dets_up(:),dets_dn(:)
#else
   integer(ik),allocatable     :: all_dets_up_cross(:),all_dets_dn_cross(:)
   integer(ik),allocatable     :: all_dets_up_2x2(:),all_dets_dn_2x2(:)
   integer(ik),allocatable    :: dets_up(:),dets_dn(:)
#endif
   integer                     :: nup_back,ndn_back,nsites_back
   integer                     :: n_det_up_12,n_det_dn_12,nup_12,ndn_12
   logical                     :: chosen,marked(256)
   integer                     :: chosen_j,counter
   integer                     :: nup_dm_tmp(256),ndn_dm_tmp(256)
   real(rk)                    :: dm_eigenvalues_tmp(256),dm_eigenvecs_tmp(256,256)
   real(rk)                    :: dm_eig_present
   integer                     :: nup_dm_present,ndn_dm_present,old_labels(256)
   integer(ik)                 :: env_max
   integer                     :: sites_in_env
   integer                     :: n_det_up_sav,n_det_dn_sav,nup_sav,ndn_sav,nsites_sav
   integer                     :: l_x_sav,l_y_sav
   logical                     :: pbc_sav

!   if (use_dm .or. use_eig) then
   if (use_dm) then
           if (cluster .eq. "cross") then
               ham_cross=.true.
               nup_back=nup
               ndn_back=ndn
               nsites_back=nsites
               nup=max(int(real(nup)*12.0/real(l_x*l_y)),1)  ! Closest up filling on the 12 site lattice - must be at least 1
               ndn=max(int(real(ndn)*12.0/real(l_x*l_y)),1)  ! Closest dn filling on the 12 site lattice - must be at least 1
               nsites=12
               nup_12=nup
               ndn_12=ndn
               n_det_up_12=int(n_choose_k(nsites,nup),i4b)
               n_det_dn_12=int(n_choose_k(nsites,ndn),i4b)
               call lanczos_hubbard(lowest_eigenvec_cross,energy_cross,all_dets_up_cross,all_dets_dn_cross)   ! Lanczos diagonalization in real space on the cross of squares
               nup=nup_back                           ! Restore nup
               ndn=ndn_back                           ! Restore ndn
               nsites=nsites_back
               ham_cross=.false.                      ! Do not do "cross lattice" calculations from here on
          else                                       ! Default is square
               ham_square=.true.
               nup_back=nup
               ndn_back=ndn
               nsites_back=nsites
               nup=max(int(real(nup)*16.0/real(l_x*l_y)),1)  ! Closest up filling on the 16 site lattice - must be at least 1
               ndn=max(int(real(ndn)*16.0/real(l_x*l_y)),1)  ! Closest dn filling on the 16 site lattice - must be at least 1
               nsites=16
               nup_12=nup
               ndn_12=ndn
               n_det_up_12=int(n_choose_k(nsites,nup),i4b)
               n_det_dn_12=int(n_choose_k(nsites,ndn),i4b)
               !print *,"Ham square before Lanczos",ham_square
               call lanczos_hubbard(lowest_eigenvec_cross,energy_cross,all_dets_up_cross,all_dets_dn_cross)   ! Lanczos diagonalization in real space on the cross of squares
               !print *,"Ham square after Lanczos",ham_square
               nup=nup_back                           ! Restore nup
               ndn=ndn_back                           ! Restore ndn
               nsites=nsites_back
               ham_square=.false.                      ! Do not do "square lattice" calculations from here on
               diag_44_done=.true.
          endif
   endif

   ! Now the lowest eigenvector has been computed and all its components are real
   ! The DM for this case has linear dimension 4^4 i.e. it is a 256x256 matrix. However it is block diagonal since the Hubbard model conserves number of particles and total Sz

   density_matrix(:,:)=0._rk

   ordered_dets=0_ik
   counter=1

   do i=0,4
        allocate (dets_up(n_choose_k(4,i)))
        dets_up=0_ik
        call constrained_dets(i,size(dets_up),dets_up)
     do j=0,4
        allocate (dets_dn(n_choose_k(4,j)))
        call constrained_dets(j,size(dets_dn),dets_dn)
        do k=1,size(dets_up)
            do l=1,size(dets_dn)
#ifdef NUM_ORBITALS_GT_127
                num_int16=16*dets_up(k)+dets_dn(l)%v(1)
#else
                num_int16=16*dets_up(k)+dets_dn(l)
#endif
                ordered_dets(counter)=num_int16
                inverse_dm_map(num_int16+1)=counter
                counter=counter+1
            enddo
        enddo
        deallocate(dets_dn)
     enddo
     deallocate(dets_up)
   enddo

   write (6,*) "Ordered dets"
   write (6,*) ordered_dets
   call flush(6)

   write (6,*) "Inverse map"
   write (6,*)
   do i=1,256
        write(6,'(i5,i5)') i-1,inverse_dm_map(i)
   enddo
   write (6,*)
   call flush(6)

!   if (use_eig) then            !!! Hack
!        if (allocated(lowest_k_eigenvec)) deallocate(lowest_k_eigenvec)
!        allocate (lowest_k_eigenvec(size(lowest_eigenvec_cross,1)))
!        lowest_k_eigenvec=lowest_eigenvec_cross
!        if (allocated(all_dets_up)) deallocate(all_dets_up)
!        if (allocated(all_dets_dn)) deallocate(all_dets_dn)
!        allocate(all_dets_up(size(all_dets_up_cross,1)))
!        allocate(all_dets_dn(size(all_dets_dn_cross,1)))
!        all_dets_up=all_dets_up_cross
!        all_dets_dn=all_dets_dn_cross
!   endif

   if (use_dm) then               ! Begin if use_dm
       ! The process is inefficient - but it probably doesnt matter for such small lattices!
       write (6,*) "Computing the density matrix....."
       call flush(6)

       if (cluster .eq. "cross") then           ! Cross
        env_max=255_ik
        sites_in_env=8
       else                                     ! Default is square
        print *,"Square cluster"
        env_max=4095_ik
        sites_in_env=12
        if (allocated(lowest_k_eigenvec)) deallocate(lowest_k_eigenvec)
        allocate (lowest_k_eigenvec(size(lowest_eigenvec_cross,1)))
        lowest_k_eigenvec=lowest_eigenvec_cross
        if (allocated(all_dets_up)) deallocate(all_dets_up)
        if (allocated(all_dets_dn)) deallocate(all_dets_dn)
        allocate(all_dets_up(size(all_dets_up_cross,1)))
        allocate(all_dets_dn(size(all_dets_dn_cross,1)))
        all_dets_up=all_dets_up_cross
        all_dets_dn=all_dets_dn_cross
       endif

       do env_up=0_ik,env_max
           nup_env=0
           do k=1,sites_in_env
               if (btest(env_up,k-1)) nup_env=nup_env+1
           enddo
           !write (6,*) "nup_env=",nup_env

           do env_dn=0_ik,env_max
               ndn_env=0
               do k=1,sites_in_env
                    if (btest(env_dn,k-1)) ndn_env=ndn_env+1
               enddo
               !write (6,*) "ndn_env=",ndn_env

               if (nup_env .le. nup_12 .and. ndn_env .le. ndn_12) then
                   do i_up=0,15_ik
                        nup_sys=0
                        do k=1,4
                            if (btest(i_up,k-1)) nup_sys=nup_sys+1
                        enddo
                        !write (6,*) "nup_sys=",nup_sys

                        if (nup_sys+nup_env .eq. nup_12) then
                            !write (6,*) "Condition 1 satisfied ","nup_sys=",nup_sys,"nup_env=",nup_env
                            do i_dn=0,15_ik
                                ndn_sys=0
                                do k=1,4
                                    if (btest(i_dn,k-1)) ndn_sys=ndn_sys+1
                                enddo
                                if (ndn_sys+ndn_env .eq. ndn_12) then
                                    !write (6,*) "Condition 2 satisfied ","ndn_sys=",ndn_sys,"ndn_env=",ndn_env
                                    i16=(16_ik*i_up)+i_dn+1
                                    ind1=(inverse_map_up12((16_ik*env_up+i_up))-1)*n_det_dn_12 + inverse_map_dn12(((16_ik*env_dn)+i_dn))
                                    !write (6,*) "i_up=",i_up,"i_dn=",i_dn,"i=",i,"ind1=",ind1
                                    call flush(6)
                                    if (abs(lowest_eigenvec_cross(ind1)).gt. 1.0e-15) then
                                        do j_up=0,15_ik
                                            nup_sys2=0
                                            do k=1,4
                                                if (btest(j_up,k-1)) nup_sys2=nup_sys2+1
                                            enddo
                                            if (nup_sys2 .eq. nup_sys) then
                                                !write (6,*) "Condition 3 satisfied"
                                                do j_dn=0,15_ik
                                                    ndn_sys2=0
                                                    do k=1,4
                                                        if (btest(j_dn,k-1)) ndn_sys2=ndn_sys2+1
                                                    enddo
                                                    if (ndn_sys .eq. ndn_sys2) then
                                                        !write (6,*) "Condition 4 satisfied"
                                                        ind2=(inverse_map_up12((16_ik*env_up+j_up))-1)*n_det_dn_12 + inverse_map_dn12((16_ik*env_dn)+j_dn)
                                                        j16=(16_ik*j_up)+j_dn+1
                                                        density_matrix(inverse_dm_map(i16),inverse_dm_map(j16))=density_matrix(inverse_dm_map(i16),inverse_dm_map(j16))+lowest_eigenvec_cross(ind1)*lowest_eigenvec_cross(ind2)
                                                    endif
                                                enddo
                                            endif
                                        enddo
                                    endif
                                endif
                            enddo
                        endif
                   enddo
                endif
           enddo
       enddo

       write (6,*) "TRACE: Finished Computing the density matrix....."
       call flush(6)
  endif  ! If use_dm

   if (allocated(dm_eigenvecs)) deallocate(dm_eigenvecs)
   allocate(dm_eigenvecs(256,256))
   if (allocated(dm_eigenvalues)) deallocate(dm_eigenvalues)
   allocate(dm_eigenvalues(256))

  if (use_dm) dm_eigenvecs=density_matrix

  if (allocated(nup_dm)) deallocate(nup_dm)
  if (allocated(ndn_dm)) deallocate(ndn_dm)

  allocate(nup_dm(256))
  allocate(ndn_dm(256))

  shift=0
  counter=1
  do i=0,4
       do j=0,4
         num=int(n_choose_k(4,i)*n_choose_k(4,j),i4b)
         ndm_start_end(i,j,1)=shift+1
         ndm_start_end(i,j,2)=shift+num
         if (use_dm) call real_symmetric_diagonalize(num,density_matrix(shift+1:shift+num,shift+1:shift+num),dm_eigenvecs(shift+1:shift+num,shift+1:shift+num),dm_eigenvalues(shift+1:shift+num))
         shift=shift+num
         do k=1,num
            nup_dm(counter)=i
            ndn_dm(counter)=j
            counter=counter+1
         enddo
        enddo
   enddo

   shift=0

   if (use_dm) then
       write (6,*) "DM eigenvalues with nup and ndn"
    !   call print_real_array_with_index(256,dm_eigenvalues)
        do i=1,256
                 write(6,'(i5,f10.5,i5,i5)') i,dm_eigenvalues(i),nup_dm(i),ndn_dm(i)
        enddo
       call flush(6)

      ! Sort by increasing nup,ndn then dm eigenvalue
      marked(:)=.false.
      do i=1,256
        chosen=.false.
        do j=1,256
          if (marked(j) .eqv. .false.) then
              if (chosen .eqv. .false.) then
                nup_dm_present=nup_dm(j)
                ndn_dm_present=ndn_dm(j)
                dm_eig_present=dm_eigenvalues(j)
                chosen=.true.
                chosen_j=j
              else
                if (nup_dm(j) .lt. nup_dm_present) then
                    nup_dm_present=nup_dm(j)
                    ndn_dm_present=ndn_dm(j)
                    dm_eig_present=dm_eigenvalues_tmp(j)
                    chosen_j=j
                elseif (nup_dm(j) .eq. nup_dm_present .and. ndn_dm(j) .lt. ndn_dm_present) then
                    nup_dm_present=nup_dm(j)
                    ndn_dm_present=ndn_dm(j)
                    dm_eig_present=dm_eigenvalues(j)
                    chosen_j=j
                elseif (nup_dm(j) .eq. nup_dm_present .and. ndn_dm(j) .eq. ndn_dm_present .and. dm_eig_present .lt. dm_eigenvalues(j)) then
                    nup_dm_present=nup_dm(j)
                    ndn_dm_present=ndn_dm(j)
                    dm_eig_present=dm_eigenvalues(j)
                    chosen_j=j
                endif
              endif
          endif
         enddo
         old_labels(i)=chosen_j
         marked(chosen_j)=.true.
      enddo

      do i=1,256
        dm_eigenvecs_tmp(:,i)=dm_eigenvecs(:,old_labels(i))
        nup_dm_tmp(i)=nup_dm(old_labels(i))
        ndn_dm_tmp(i)=ndn_dm(old_labels(i))
        dm_eigenvalues_tmp(i)=dm_eigenvalues(old_labels(i))
      enddo

      nup_dm=nup_dm_tmp
      ndn_dm=ndn_dm_tmp
      dm_eigenvalues=dm_eigenvalues_tmp
      dm_eigenvecs=dm_eigenvecs_tmp

      write (6,*)
      write (6,*) "DM eigenvalues with nup and ndn after sorting"
      nup_dm_present=0
      ndn_dm_present=0
      i=1
      write(6,'(i5,f10.5,i5,i5)') i,dm_eigenvalues(i),nup_dm(i),ndn_dm(i)
      do i=2,256
                 if (nup_dm_present .ne. nup_dm(i) .or. ndn_dm_present .ne. ndn_dm(i)) then
                    nup_dm_present=nup_dm(i)
                    ndn_dm_present=ndn_dm(i)
                    write (6,*)
                 endif
                 write(6,'(i5,f10.5,i5,i5)') i,dm_eigenvalues(i),nup_dm(i),ndn_dm(i)
      enddo
      write (6,*)
  endif        ! If use_dm

  if (use_eye) dm_eigenvecs=eye_256 ! A hack for debugging code - must give real space result

  write (6,*) "Starting and ending points and size of nup,ndn sectors"
  do i=0,4
    do j=0,4
             write (6,*)
             write(6,'(/,''================================================================================================================================================================================================'')')
             write(6,'(5i5)') i,j,ndm_start_end(i,j,1),ndm_start_end(i,j,2),ndm_start_end(i,j,2)-ndm_start_end(i,j,1)+1
             write (6,*)
             if (use_dm) then
                write(6,'(/,''            Density Matrix               '')')
                write (6,*)
                call print_real_matrix(ndm_start_end(i,j,2)-ndm_start_end(i,j,1)+1,ndm_start_end(i,j,2)-ndm_start_end(i,j,1)+1,density_matrix(ndm_start_end(i,j,1):ndm_start_end(i,j,2),ndm_start_end(i,j,1):ndm_start_end(i,j,2)))
                write (6,*)
             endif
             write(6,'(/,''            DM/Energy Eigenvectors              '')')
             write (6,*)
             call print_real_matrix(ndm_start_end(i,j,2)-ndm_start_end(i,j,1)+1,ndm_start_end(i,j,2)-ndm_start_end(i,j,1)+1,dm_eigenvecs(ndm_start_end(i,j,1):ndm_start_end(i,j,2),ndm_start_end(i,j,1):ndm_start_end(i,j,2)))
             write(6,'(/,''================================================================================================================================================================================================'')')
             write (6,*)
    enddo
  enddo
  write (6,*)

  if (allocated(ham_on_2x2)) deallocate(ham_on_2x2)
  allocate(ham_on_2x2(256,256))
  ham_on_2x2(:,:)=0._rk

  n_det_up_sav=n_det_up
  n_det_dn_sav=n_det_dn
  nsites_sav=nsites
  nup_sav=nup
  ndn_sav=ndn
  l_x_sav=l_x
  l_y_sav=l_y
  pbc_sav=pbc

  l_x=2
  l_y=2
  pbc=.false.
  nsites=4
  do i=0,4
      nup=i
      do j=0,4
        ndn=j
        call make_hubbard_matrix_2d(tmp_ham_2x2,lowest_evec_2x2,lowest_eval_2x2,all_dets_up_2x2,all_dets_dn_2x2)
        do k=1,size(lowest_evec_2x2)
#ifdef NUM_ORBITALS_GT_127
            ind1=(16_ik*all_dets_up_2x2(k))+all_dets_dn_2x2(k)%v(1)+1
#else
            ind1=(16_ik*all_dets_up_2x2(k))+all_dets_dn_2x2(k)+1
#endif
            ind1_int4=inverse_dm_map(ind1)
            do l=1,size(lowest_evec_2x2)
#ifdef NUM_ORBITALS_GT_127
                ind2=(16_ik*all_dets_up_2x2(l))+all_dets_dn_2x2(l)%v(1)+1
#else
                ind2=(16_ik*all_dets_up_2x2(l))+all_dets_dn_2x2(l)+1
#endif
                ind2_int4=inverse_dm_map(ind2)
                ham_on_2x2(ind1_int4,ind2_int4)=tmp_ham_2x2(k,l)
            enddo
        enddo
      enddo
  enddo

  n_det_up=n_det_up_sav
  n_det_dn=n_det_dn_sav
  nsites=nsites_sav
  nup=nup_sav
  ndn=ndn_sav
  l_x=l_x_sav
  l_y=l_y_sav
  pbc=pbc_sav
  call flush(6)

  if (use_eig) dm_eigenvecs=ham_on_2x2

  if (use_eig) then
      shift=0
      dm_eigenvalues(:)=0._rk
      do i=0,4
           do j=0,4
             num=int(n_choose_k(4,i)*n_choose_k(4,j),i4b)
             call real_symmetric_diagonalize(num, ham_on_2x2(shift+1:shift+num,shift+1:shift+num),dm_eigenvecs(shift+1:shift+num,shift+1:shift+num),dm_eigenvalues(shift+1:shift+num))
             shift=shift+num
            enddo
       enddo
       shift=0
  endif

  ham_on_2x2=matmul(ham_on_2x2,dm_eigenvecs)
  ham_on_2x2=matmul(transpose(dm_eigenvecs),ham_on_2x2)

  write (6,*) "Similarity Transformed Hamiltonian on 2x2 patch organized by sector"
  do i=0,4
    do j=0,4
             write (6,*)
             write(6,'(/,''================================================================================================================================================================================================'')')
             write(6,'(5i5)') i,j,ndm_start_end(i,j,1),ndm_start_end(i,j,2),ndm_start_end(i,j,2)-ndm_start_end(i,j,1)+1
             write (6,*)
             write(6,'(/,''            Hamiltonian on patch               '')')
             write (6,*)
             call print_real_matrix(ndm_start_end(i,j,2)-ndm_start_end(i,j,1)+1,ndm_start_end(i,j,2)-ndm_start_end(i,j,1)+1,ham_on_2x2(ndm_start_end(i,j,1):ndm_start_end(i,j,2),ndm_start_end(i,j,1):ndm_start_end(i,j,2)))
             write (6,*)
    enddo
  enddo
  write (6,*)

  end subroutine density_matrix_2by2

!=!=====================================================================================================================================================================
  subroutine project_44wf_dmbasis()

  integer                 :: i,i1,j,j1,k,l,num,ind(4),num_blocks,total_nup_poss,total_ndn_poss
  real(rk)                :: prod
  integer                 :: nup_poss(625,4),ndn_poss(625,4)
  integer                 :: counter
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable :: all_dets_up_dm(:),all_dets_dn_dm(:)
#else
  integer(ik),allocatable :: all_dets_up_dm(:),all_dets_dn_dm(:)
#endif
  integer                 :: nup_block(4),ndn_block(4)
  integer,allocatable     :: combined_nup_ndn_poss(:,:)
  integer                 :: shift,total,repeat_frequency
  real(rk)                :: norm
  integer                 :: combined_ind,how_many_times
  integer                 :: beta_basis_start_end(625*625,2)
  integer(ik)             :: which_vector

  if (l_x .ne. 4 .or. l_y .ne. 4) return  ! Valid only for 4x4

  ! We have the lowest eigenvector of the 4x4 Hubbard model and we want to find its expansion in the density matrix basis
  ! Load wavefunction from disk - we dont want to be doing Lanczos everytime we check the projection
  allocate(lowest_vec_dmbasis(size(lowest_k_eigenvec,1)))
  lowest_vec_dmbasis(:)=0._rk

  num_blocks=4

  counter=1
  ! Make nup possibility list on blocks 1 2 3 4
  do i=0,624
      num=i
      do k=1,4
        nup_poss(counter,k)=mod(num,5)
        num=num/5
      enddo
      if (sum(nup_poss(counter,:)) .eq. nup) then
        counter=counter+1
      else
        nup_poss(counter,:)=0
      endif
  enddo
  total_nup_poss=counter-1

  ! Make ndn possibility list on blocks 1 2 3 4
  counter=1
  do i=0,624
      num=i
      do k=1,4
        ndn_poss(counter,k)=mod(num,5)
        num=num/5
      enddo
      if (sum(ndn_poss(counter,:)) .eq. ndn) then
        counter=counter+1
      else
        ndn_poss(counter,:)=0
      endif
  enddo
  total_ndn_poss=counter-1

  write (6,*) "Nup possibilities on the blocks"
  do i=1,total_nup_poss
             write(6,'(4i5)') (nup_poss(i,j),j=1,4)
  enddo
  call flush(6)

  write (6,*) "Ndn possibilities on the blocks"
  do i=1,total_ndn_poss
             write(6,'(4i5)') (ndn_poss(i,j),j=1,4)
  enddo
  call flush(6)

  ! Expand maps
  allocate(all_dets_up_dm(size(lowest_k_eigenvec)))
  allocate(all_dets_dn_dm(size(lowest_k_eigenvec)))
  allocate(combined_nup_ndn_poss(total_nup_poss*total_ndn_poss,8))

  call sort(lowest_k_eigenvec,all_dets_up,all_dets_dn)
  write(6,'(/,''==============================================================='')')
  write(6,'(/,''                     Eigenvector from ED                       '')')
  write(6,'(/,''==============================================================='')')
  do i=1,size(lowest_k_eigenvec)
    write(6,'(2i10,10g20.5)') all_dets_up(i), all_dets_dn(i), lowest_k_eigenvec(i)
    call flush(6)
  enddo
  write (6,*)
  call flush(6)

  write(6,'(/,'' Generating all nup,ndn combinations on blocks.......'')')
  call flush(6)

  ! Combine all possibilities on blocks 1 up 2up 3up 4up 1dn 2dn 3dn 4dn
  counter=1
  do i=1,total_nup_poss
            do j=1,total_ndn_poss
                do i1=1,4
                    combined_nup_ndn_poss(counter,i1)=nup_poss(i,i1)
                    combined_nup_ndn_poss(counter,i1+4)=ndn_poss(j,i1)
                enddo
                counter=counter+1
            enddo
  enddo

!  write(6,'(/,''==============================================================='')')
!  write (6,*) " all nup,ndn possibilities on the blocks"
!  write(6,'(/,''  1u    2u    3u   4u   1d   2d   3d   4d                      '')')
!  write(6,'(/,''==============================================================='')')
!  do i=1,total_ndn_poss*total_nup_poss
!             write(6,'(8i5)') (combined_nup_ndn_poss(i,j),j=1,8)
!             call flush(6)
!  enddo
!  write (6,*)
!  call flush(6)

  all_dets_up_dm=0_ik
  all_dets_dn_dm=0_ik

  ! Expand maps for all_dets_dm
  counter=1
  shift=0

  if (use_dm)  write(6,'(/,'' Generating full  2x2 DM     basis.......'')')
  if (use_eig) write(6,'(/,'' Generating full  2x2 Eig    basis.......'')')
  call flush(6)

  do i=1,total_ndn_poss*total_nup_poss
    call flush(6)
    total=1
    do j=1,4 ! which block
        total=total*(-ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),1)+ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),2)+1)
    enddo
    call flush(6)
    ! total possibilities with these numbers on all blocks
    combined_ind=0
    do j=1,4
        combined_ind=combined_ind+(5**(j-1))*(combined_nup_ndn_poss(i,4+j))+(5**(4+j-1))*combined_nup_ndn_poss(i,j)
    enddo
    beta_basis_start_end(combined_ind,1)=shift+1
    beta_basis_start_end(combined_ind,2)=shift+total

    how_many_times=1
    do j=1,4    ! over blocks
        repeat_frequency=total/(ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),2)-ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),1)+1)
        repeat_frequency=repeat_frequency/how_many_times
        counter=1
        call flush(6)
        do j1=1,repeat_frequency
            do k=ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),1),ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),2)
              ! k is a number between 1 and 256 indicating which dm eigenvector we are looking at
              do l=1,how_many_times
                  do i1=1,8
                    if (i1 .gt. 4 .and. btest(k,i1-1)) all_dets_up_dm(counter+shift)=ibset(all_dets_up_dm(counter+shift),4*(j-1)+i1-5)
                    if (i1 .le. 4 .and. btest(k,i1-1)) all_dets_dn_dm(counter+shift)=ibset(all_dets_dn_dm(counter+shift),4*(j-1)+i1-1)
                  enddo
                  counter=counter+1
              enddo
            enddo
         enddo
         how_many_times=how_many_times*(ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),2)-ndm_start_end(combined_nup_ndn_poss(i,j),combined_nup_ndn_poss(i,4+j),1)+1)
    enddo
    shift=shift+total
  enddo

! Print beta basis
!  write (6,*) "Beta basis start and end"
!  do i=1,size(beta_basis_start_end,1)
!        if (beta_basis_start_end(i,1) .ne. 0 .and. beta_basis_start_end(i,2) .ne. 0) write(6,'(3i10)') i,beta_basis_start_end(i,1),beta_basis_start_end(i,2)
!  enddo

! Print expanded basis for checking
!  write (6,*) "Expanded map"
!do i=1,size(all_dets_up_dm)
!    do j=1,4 ! over blocks
!        which_vector=0_ik
!        do i1=1,4     ! scan bits
!            if (btest(all_dets_up_dm(i),4*(j-1)+i1-1)) which_vector=ibset(which_vector,4+i1-1)
!            if (btest(all_dets_dn_dm(i),4*(j-1)+i1-1)) which_vector=ibset(which_vector,i1-1)
!        enddo
!        ind(j)=which_vector
!    enddo
!    write(6,'(4i5)') (ind(j),j=1,4)
!enddo
!  write (6,*)

  if (use_dm)  write(6,'(/,'' Performing projection from real space to 2x2 DM basis.......'')')
  if (use_eig) write(6,'(/,'' Performing projection from real space to 2x2 Eig basis.......'')')
  call flush(6)

! Perform projection
  do j=1,size(lowest_k_eigenvec,1)
      if (abs(lowest_k_eigenvec(j)) .gt. 1.0e-10) then
          do i=1,4  ! over blocks
             nup_block(i)=0
             ndn_block(i)=0
             ind(i)=0_ik
             do k=0,3  ! 4 bits represent up, 4 bits represent down configuration on block
                if (btest(all_dets_up(j),4*(i-1)+k)) then
                    nup_block(i)=nup_block(i)+1
                    ind(i)=ibset(ind(i),4+k)
                endif
                if (btest(all_dets_dn(j),4*(i-1)+k)) then
                    ndn_block(i)=ndn_block(i)+1
                    ind(i)=ibset(ind(i),k)
                endif
             enddo
          enddo
          combined_ind=0
          do i=1,4
            combined_ind=combined_ind+((5**(i-1))*(ndn_block(i)))+((5**(4+i-1))*nup_block(i))
          enddo
          do k=beta_basis_start_end(combined_ind,1),beta_basis_start_end(combined_ind,2)
            !write (6,*) "========================================================================"
            prod=1._rk
            do i=1,4           !over blocks
                which_vector=0_ik
                do i1=1,4
                    if (btest(all_dets_up_dm(k),4*(i-1)+i1-1)) which_vector=ibset(which_vector,4+i1-1)
                    if (btest(all_dets_dn_dm(k),4*(i-1)+i1-1)) which_vector=ibset(which_vector,i1-1)
                enddo
                !write (6,*) "block        =",i
                !write (6,*) "which_vector =",which_vector
                if (which_vector .gt. 0) then
                    prod=prod*dm_eigenvecs(inverse_dm_map(ind(i)+1),which_vector)
                    if (abs(prod) .le. 1.0e-10) exit
                else
                    prod=0._rk
                endif
                !write (6,*) "prod         =",prod
            enddo
            lowest_vec_dmbasis(k)=lowest_vec_dmbasis(k)+(prod*lowest_k_eigenvec(j))
          enddo
       endif
  enddo

  call sort(lowest_vec_dmbasis,all_dets_up_dm,all_dets_dn_dm)

  write (6,*)
  if (use_dm)  write (6,*) "Sorted eigenvector in 2x2 DM  basis obtained from projection of exact vector from ED"
  if (use_eig) write (6,*) "Sorted eigenvector in 2x2 Eig basis obtained from projection of exact vector from ED"
  do i=1,size(lowest_vec_dmbasis)
    write(6,'(2i10,10g20.5)') all_dets_up_dm(i), all_dets_dn_dm(i), lowest_vec_dmbasis(i)
    call flush(6)
  enddo
  write (6,*)
  call flush(6)

  norm=dot_product(lowest_vec_dmbasis,lowest_vec_dmbasis)

  write (6,*) "Norm of the projected vector is",norm
  write(6,*)
  write (6,*) "If norm is not 1 your code has a BUG!"
  write(6,*)
  stop
  end subroutine project_44wf_dmbasis

!====================================================================================================================

  subroutine ham_most_probable_basis()
!====================================================================================================================
 ! Description : Retain only few states per block (nsb=10 means 10 states per block) based on an energy or DM criterion
 !               Given this local truncated basis construct the Hamiltonian in the tensor product basis.
 !               Diagonalize this Hamiltonian and look at the lowest eigenvalue and eigen energy
 ! Created by  : Hitesh J. Changlani
 ! Date        : Feb 1, 2012
!====================================================================================================================
  implicit none

  ! Local variables
  integer                   :: i,j,i1,counter
  integer                   :: nsb=4                                                   ! Number of states per block
  integer(ik),allocatable   :: list_of_states(:)
  integer                   :: nup_tot,ndn_tot,ni_states,retained_states                ! Total number of "important states" on all blocks combined
  real(rk),allocatable      :: ham(:,:),small_ham(:,:),eigenvectors_ham(:,:),eigenvalues_ham(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable   :: dets_up(:),dets_dn(:)
#else
  integer(ik),allocatable   :: dets_up(:),dets_dn(:)
#endif
  integer(ik)               :: which_state(nsites/4)
  real(rk)                  :: dmevals_copy(256)
  integer(ik)               :: sorted_states(256)
  real(rk),allocatable      :: lowest_eigenvector(:)
  real(rk),allocatable      :: trunc_lowest_eigenvector(:)
  real(rk)                  :: lowest_eigenvalue,maxv,trunc_energy
  real(rk)                  :: navg

  nsb=min(nsb,256)
  write  (6,*) "Constructing Hamiltonian in the most probable/ lowest energy basis..............."
  write  (6,*) "In case of energy this is like one conventional RG transformation with truncation"
  call flush(6)

  ! Get the most important list of states based on DM eigenvalues or energy criterion. There are nsb such states
  allocate(list_of_states(nsb))
  allocate(dets_up(nsb**(nsites/4)))
  allocate(dets_dn(nsb**(nsites/4)))

  dmevals_copy(:)=dm_eigenvalues(:)

  do i=1,256
    sorted_states(i)=i
  enddo

  if (use_dm)  then
    call sort(dmevals_copy,sorted_states,sorted_states,.false.)
  endif

  navg=real(nup)/real(nsites)
  if (use_eig) then
    do i=1,256
        dmevals_copy(i)=dmevals_copy(i)-(navg)*U*real(nup_dm(i)+ndn_dm(i))+0.1*(real(nup_dm(i)-4*navg)**2)+0.1*(real(ndn_dm(i)-4*navg)**2)  ! Add chemical potential+ cheat potential
    enddo
    dmevals_copy=-dmevals_copy
    call sort(dmevals_copy,sorted_states,sorted_states,.true.)  ! I have flipped the energy - so that lowest energy is at top of list
    dmevals_copy=-dmevals_copy
  endif

  write (6,*)
  call flush(6)

  do i=1,256
    if (use_eig) write(6,'(f10.5,9i6)') dmevals_copy(i),sorted_states(i),nup_dm(sorted_states(i)),ndn_dm(sorted_states(i))
    if (use_dm) then
        write(6,'(f10.5,f10.5,9i6)') dmevals_copy(i),ham_on_2x2(sorted_states(i),sorted_states(i)),sorted_states(i),nup_dm(sorted_states(i)),ndn_dm(sorted_states(i))
    endif
  enddo

  call flush(6)

  counter=0
  do i=1,nsb
     list_of_states(i)=sorted_states(i)
  enddo

  write (6,*)
  write (6,*) "Important states on block"
  call flush(6)
  do i=1,nsb
    write(6,'(9i6)') list_of_states(i),nup_dm(list_of_states(i)),ndn_dm(list_of_states(i))
    call flush(6)
  enddo
  call flush(6)

  write (6,*) "Enumerating states on all blocks combined...."
  call flush(6)
  counter=0
  write (6,*) "nsb**(nsites/4)",nsb**(nsites/4)

  do i=0,(nsb**(nsites/4))-1                                 ! Number of states per block ^ Number of blocks
     ! Get string of states - Convert integer i to base nsb
     !write (6,*) "i=",i
     call flush(6)
     call convert_integer_to_nary(int(i,ik),nsb,nsites/4,which_state)
     !write (6,*) which_state
     call flush(6)
     nup_tot=0
     ndn_tot=0
     do j=1,nsites/4
        which_state(j)=list_of_states(which_state(j)+1)    ! Convert it to the index in the full list
        nup_tot=nup_tot+nup_dm(which_state(j))
        ndn_tot=ndn_tot+ndn_dm(which_state(j))
     enddo
     !write (6,*) "After substitution", which_state
     call flush(6)

     if (nup_tot .eq. nup .and. ndn_tot .eq. ndn) then   ! allowed state since it has the right particle numbers
        counter=counter+1
        dets_up(counter)=0_ik
        dets_dn(counter)=0_ik
        do j=1,4                    !over blocks
                do i1=1,4
                    if (btest(which_state(j),i1-1)) dets_dn(counter)=ibset(dets_dn(counter),4*(j-1)+i1-1)
                enddo
                do i1=5,8
                    if (btest(which_state(j),i1-1)) dets_up(counter)=ibset(dets_up(counter),4*(j-1)+i1-4-1)
                enddo
        enddo
     endif
  enddo

  ni_states=counter
  write (6,*) "Number of important states is ",ni_states

  allocate (ham(ni_states,ni_states))

  write (6,*) "Constructing Hamiltonian in important basis"
  ham(:,:)=0._rk
  do i=1,ni_states
    do j=1,ni_states
        call hamiltonian_hubbard_dm(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j),ham(i,j))
    enddo
  enddo
  write (6,*) "Finished constructing Hamiltonian in important basis"

  write (6,*) " Now diagonalizing... "
  call flush(6)
  allocate(lowest_eigenvector(ni_states))

  ! Diagonalize this Hamiltonian
  if (ni_states .le. 2000) then
      allocate (eigenvalues_ham(ni_states))
      allocate (eigenvectors_ham(ni_states,ni_states))
      call real_symmetric_diagonalize(ni_states,ham,eigenvectors_ham,eigenvalues_ham)
      call matrix_lanczos(ni_states,lowest_eigenvector,lowest_eigenvalue,ham)
      write (6,*) "Finished diagonalizing"
      call flush(6)

      write (6,*) "All Eigenvalues of Hamiltonian in important basis"
      call print_real_matrix(size(eigenvalues_ham),1,eigenvalues_ham)
      call flush(6)

      write (6,*) "Lowest eigenvector of Hamiltonian in important basis"
      call sort(lowest_eigenvector,dets_up,dets_dn)
      call print_real_matrix(size(lowest_eigenvector),1,lowest_eigenvector)
      call flush(6)
  else
      call matrix_lanczos(ni_states,lowest_eigenvector,lowest_eigenvalue,ham)
      write (6,*) "Lowest eigenvector of Hamiltonian in important basis"
      call sort(lowest_eigenvector,dets_up,dets_dn)
      call print_real_matrix(size(lowest_eigenvector),1,lowest_eigenvector)
  endif

  maxv=abs(lowest_eigenvector(1))
  counter=0
  do i=1,ni_states
    if (abs(lowest_eigenvector(i))>maxv/100.0) then
        counter=counter+1
    else
        exit
    endif
  enddo
  retained_states=min(counter,10)
  write (6,*) "Number of retained states is ",retained_states
  allocate (trunc_lowest_eigenvector(retained_states))
  trunc_lowest_eigenvector(:)=0._rk
  trunc_lowest_eigenvector(:)=lowest_eigenvector(1:retained_states)
  ! After sorting of eigenvector, reconstruct Hamiltonian or be smart about it and extract from prev calculation
  allocate (small_ham(retained_states,retained_states))

  write (6,*) "Constructing Hamiltonian in retained basis... Inefficient!!!"
  small_ham(:,:)=0._rk
  do i=1,retained_states
    do j=1,retained_states
        call hamiltonian_hubbard_dm(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j),small_ham(i,j))
    enddo
  enddo
  write (6,*) "Finished constructing Hamiltonian in retained basis"

  trunc_energy=dot_product(trunc_lowest_eigenvector,matmul(small_ham,trunc_lowest_eigenvector))/dot_product(trunc_lowest_eigenvector,trunc_lowest_eigenvector)
  write (6,*) "Energy of truncated lowest eigenvector is ",trunc_energy

  ndet_psi_t=retained_states
  allocate(cdet_psi_t(ndet_psi_t))
  allocate(dets_up_psi_t(ndet_psi_t))
  allocate(dets_dn_psi_t(ndet_psi_t))
  dets_up_psi_t(:)=dets_up(1:retained_states)
  dets_dn_psi_t(:)=dets_dn(1:retained_states)
  cdet_psi_t(:)=trunc_lowest_eigenvector(:)

  end subroutine ham_most_probable_basis

!==================================================================================================
   subroutine cpmc

   implicit none

   integer              :: nbr,i,j
   logical              :: allowed

   ! Initialization
   allocate(h_k(nsites,nsites))
   allocate(exp_h_k(nsites,nsites))
   allocate(cpmc_orbitals(nsites,nsites))
   allocate(phi_t_up(nsites,nup))
   allocate(phi_t_dn(nsites,ndn))
   allocate(k_eigenvalues(nsites))

   h_k(:,:)=0._rk

   do i=1,nsites
        do j=1,4
            call get_nbr(l_x,l_y,pbc,i,j,nbr,allowed)
            if ((allowed) .and. (nbr>i)) then
                h_k(i,nbr)=-t
                h_k(nbr,i)=-t
            endif
         enddo
   enddo

   ! RHF orbitals chosen as CPMC trial wavefunction
   ! In general one could choose UHF or multideterminant wavefunction
   call real_symmetric_diagonalize(nsites,h_k,cpmc_orbitals,k_eigenvalues)

   ! Half K projector
   exp_h_k(:,:)=0._rk
   do i=1,nsites
    exp_h_k(i,i)=exp(-tau*0.5*k_eigenvalues(i))
   enddo

   exp_h_k=matmul(cpmc_orbitals,exp_h_k)
   exp_h_k=matmul(exp_h_k,transpose(cpmc_orbitals))

   phi_t_up(:,1:nup)=cpmc_orbitals(:,1:nup)
   phi_t_dn(:,1:nup)=cpmc_orbitals(:,1:ndn)

   ! Trial Kinetic energy
   e_k_cpmc=sum(k_eigenvalues(1:nup))+sum(k_eigenvalues(1:ndn))

   allocate(n_r_up(nsites))
   allocate(n_r_dn(nsites))

   n_r_up(:)=0
   n_r_dn(:)=0

   ! Diagonal entries of one body density matrix
   do i=1,nsites
    do j=1,nup
        n_r_up(i)=n_r_up(i)+(phi_t_up(i,j)*phi_t_up(i,j))
    enddo
   enddo

   do i=1,nsites
    do j=1,ndn
        n_r_dn(i)=n_r_dn(i)+(phi_t_dn(i,j)*phi_t_dn(i,j))
    enddo
   enddo

   e_v_cpmc=U*dot_product(n_r_up,n_r_dn)         ! Trial Potential energy
   e_t_cpmc=e_k_cpmc+e_v_cpmc                    ! Total energy

   allocate(phi_up(2*nwalk_cpmc,nsites,nup))
   allocate(phi_dn(2*nwalk_cpmc,nsites,ndn))
   allocate(walk_wt_cpmc(2*nwalk_cpmc))
   allocate(overlaps_cpmc(2*nwalk_cpmc))

   phi_up(:,:,:)=0._rk
   phi_dn(:,:,:)=0._rk
   walk_wt_cpmc(:)=0._rk
   overlaps_cpmc(:)=0._rk

   do i=1,nwalk_cpmc
        phi_up(i,:,:)=phi_t_up(:,:)
        phi_dn(i,:,:)=phi_t_dn(:,:)
        walk_wt_cpmc(i)=1._rk
        overlaps_cpmc(i)=1._rk
   enddo

   ! Data needed for Auxiliary field
   fac_norm=(e_t_cpmc-0.5_rk*U*nelec)*tau
   gamma_=cosh(exp(0.5*tau*U))

   do i=1,2
       do j=1,2
            aux_fld(i,j)=exp(gamma_*(-1)**(i+j))
       enddo
   enddo

   aux_fld_even=exp(gamma_)
   aux_fld_odd=exp(-gamma_)

   end subroutine cpmc
!===================================================================================================

   subroutine act_halfK(phi_up,phi_dn,walk_wt,overlap)

   implicit none

   real(rk),intent(inout):: phi_up(nsites,nup),phi_dn(nsites,ndn)
   real(rk),intent(inout):: walk_wt
   real(rk),intent(inout):: overlap

   real(rk) :: overlap_ratio,overlap_new,detup,detdn
   real(rk) :: inv0_matrix_up(nup,nup),inv0_matrix_dn(ndn,ndn)

   ! Act e^-tau*K/2 on the walker
   phi_up(:,:)=matmul(exp_h_K,phi_up(:,:))
   phi_dn(:,:)=matmul(exp_h_K,phi_dn(:,:))

   ! Update overlap
    inv0_matrix_up=matmul(transpose(phi_t_up),phi_up)
    call matinv(inv0_matrix_up,nup,detup)
    inv0_matrix_dn=matmul(transpose(phi_t_dn),phi_dn)
    call matinv(inv0_matrix_dn,ndn,detdn)
    overlap_new=1._rk/(detup*detdn);
    overlap_ratio=overlap_new/overlap;

   ! Apply constrained path
   if (overlap_ratio .gt. 0._rk) then
     overlap=overlap_new
     walk_wt=walk_wt*overlap_ratio
   else
     walk_wt=0._rk
   endif

   end subroutine act_halfK

!===================================================================================================
!   subroutine act_V(phi_up,phi_dn,walk_wt,overlap)

!   real(rk),intent(inout) ::phi_up(nsites,nup),phi_dn(nsites,ndn)
!   real(rk),intent(inout) :: walk_wt,overlap
!   end subroutine act_V

!===================================================================================================
  subroutine check_momentum_conservation(det_up,det_dn,valid_combination)
!===================================================================================================

    ! AAH, 7 Mar 2012
    ! Not presently being used.

    implicit none

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
#else
    integer(ik),intent(in) :: det_up,det_dn
#endif
    logical,intent(out) :: valid_combination

    integer :: i,tmp_momentum(2)

    tmp_momentum(:) = 0

    do i=1,nsites
      if (btest(det_up,i-1)) then
        tmp_momentum(1) = tmp_momentum(1) + k_vectors(1,i)
        tmp_momentum(2) = tmp_momentum(2) + k_vectors(2,i)
      endif
      if (btest(det_dn,i-1)) then
        tmp_momentum(1) = tmp_momentum(1) + k_vectors(1,i)
        tmp_momentum(2) = tmp_momentum(2) + k_vectors(2,i)
      endif
    enddo

    valid_combination=.false.
    if (mod(tmp_momentum(1)-ktot(1),2*l_x)==0.and.mod(tmp_momentum(2)-ktot(2),2*l_y)==0)  valid_combination=.true.

  end subroutine check_momentum_conservation

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
  integer,     dimension((nwalk+1)/2), intent (out) :: temp_i_2                 ! temporary array actually neither in nor out
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: key_up(nwalk), key_dn(nwalk)
  type(ik_vec), dimension((nwalk+1)/2), intent (out) :: temp_i16_up, temp_i16_dn ! temporary array actually neither in nor out
  type(ik_vec) t0
#else
  integer(ik), intent(inout) :: key_up(nwalk), key_dn(nwalk)
  integer(ik), dimension((nwalk+1)/2), intent (out) :: temp_i16_up, temp_i16_dn ! temporary array actually neither in nor out
  integer(ik) t0
#endif
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
#else
  integer(ik), intent(inout) :: a_up(na), a_dn(na), c_up(nc), c_dn(nc) ! b_up overlays c_up(na+1:nc), and similarly for b_dn, c_dn
  integer(ik), intent(in)    :: b_up(nb), b_dn(nb)
#endif
  integer    , intent(inout) :: a2(na), c2(nc)                         ! b2   overlays c2(na+1:nc)
  integer    , intent(in)    :: b2(nb)
  integer :: i,j,k

  i = 1; j = 1; k = 1;
  do while(i <= abs(na) .and. j <= abs(nb))
    if (abs(a_up(i)) < abs(b_up(j)) .or. (a_up(i) == b_up(j) .and. abs(a_dn(i)) <= abs(b_dn(j)))) then
      c_up(k) = a_up(i)
      c_dn(k) = a_dn(i)
      c2(k) = a2(i)
      i = i+1
    else
      c_up(k) = b_up(j)
      c_dn(k) = b_dn(j)
      c2(k) = b2(j)
      j = j+1
    endif
    k = k + 1
  enddo

  do while (i <= abs(na))
    c_up(k) = a_up(i)
    c_dn(k) = a_dn(i)
    c2(k) = a2(i)
    i = i + 1
    k = k + 1
  enddo

  return
end subroutine merge2_up_dn
!-------------------------------------------------------------------------------

  !===========================================================================
  subroutine symmetry_reduce_hubbardk(nimp,imp_up,imp_dn,rep_imp_up,rep_imp_dn)
  ! Purpose : Reduce list by accounting for symmetrized determinants
  ! Created : Hitesh Changlani, April 13 2012

  !===========================================================================

  implicit none
  integer                            :: nimp
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in)             :: imp_up(nimp),imp_dn(nimp)
  type(ik_vec),allocatable,intent(out):: rep_imp_up(:),rep_imp_dn(:)
  type(ik_vec)                        :: tmp_rep_up(nimp),tmp_rep_dn(nimp)
  type(ik_vec)                        :: new_dets_up(16),new_dets_dn(16)
  type(ik_vec)                        :: rep_up,rep_dn
#else
  integer(ik),intent(in)             :: imp_up(nimp),imp_dn(nimp)
  integer(ik),allocatable,intent(out):: rep_imp_up(:),rep_imp_dn(:)
  integer(ik)                        :: tmp_rep_up(nimp),tmp_rep_dn(nimp)
  integer(ik)                        :: new_dets_up(16),new_dets_dn(16)
  integer(ik)                        :: rep_up,rep_dn
#endif
  integer                            :: i,nsym_imp
  real(rk)                           :: phases(16)
  real(rk)                           :: norm,phase_w_rep
  integer                            :: num_distinct

  nsym_imp=0

  do i=1,nimp
    call generate_fourfold_k_configs(c4_map,reflection_map,z,p,imp_up(i),imp_dn(i),new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)
    if (hubbard_ipr.ge.1) then
      write (6,*) "c4_map=",c4_map
      write (6,*) "reflection_map=",reflection_map
      write (6,*) "z,p=",z,p
      write (6,*) "imp_up,dn=",imp_up(i),imp_dn(i)
      write (6,*) "new_dets_up,dn=",new_dets_up,new_dets_dn
      write (6,*) "phases=",phases
      write (6,*) "rep_up,dn=",rep_up,rep_dn
      write (6,*) "phase_w_rep=",phase_w_rep
      write (6,*) "norm=",norm
      write (6,*) "num_distinct=",num_distinct
    endif
    if ((imp_up(i) .eq. rep_up) .and. (imp_dn(i) .eq. rep_dn)) then
        nsym_imp=nsym_imp+1
        tmp_rep_up(nsym_imp)=rep_up
        tmp_rep_dn(nsym_imp)=rep_dn
    endif
  enddo

  write (6,*) "nsym_imp=",nsym_imp

  allocate(rep_imp_up(nsym_imp))
  allocate(rep_imp_dn(nsym_imp))

  do i=1,nsym_imp
    rep_imp_up(i)=tmp_rep_up(i)
    rep_imp_dn(i)=tmp_rep_dn(i)
  enddo

  end subroutine symmetry_reduce_hubbardk

  !===========================================================================
  subroutine symmetry_reduce_and_replace_hubbardk(imp_up,imp_dn)
  ! Purpose : Reduce and replace determinant list by accounting for symmetrized determinants
  ! Created : Hitesh Changlani, April 13 2012

  implicit none
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable,intent(inout) :: imp_up(:),imp_dn(:)
  type(ik_vec),allocatable               :: rep_imp_up(:),rep_imp_dn(:)
#else
  integer(ik),allocatable,intent(inout) :: imp_up(:),imp_dn(:)
  integer(ik),allocatable               :: rep_imp_up(:),rep_imp_dn(:)
#endif

  ! Local
  integer                               :: nimp

  nimp=size(imp_up)

  call symmetry_reduce_hubbardk(nimp,imp_up,imp_dn,rep_imp_up,rep_imp_dn)

  nimp=size(rep_imp_up)

  deallocate(imp_up)
  deallocate(imp_dn)
  allocate(imp_up(nimp))
  allocate(imp_dn(nimp))

  imp_up(:)=rep_imp_up(:)
  imp_dn(:)=rep_imp_dn(:)

  end subroutine symmetry_reduce_and_replace_hubbardk

  !===========================================================================
  subroutine symmetry_reduce_hubbardk_without_reallocating(nimp,imp_up,imp_dn)
  ! Purpose : Reduce list by accounting for symmetrized determinants
  ! Created : A Holmes, 10 Sep 2012

  !===========================================================================

  implicit none
  integer,intent(inout)              :: nimp
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(inout)             :: imp_up(nimp),imp_dn(nimp)
  type(ik_vec)                        :: tmp_rep_up(nimp),tmp_rep_dn(nimp)
  type(ik_vec)                        :: new_dets_up(16),new_dets_dn(16)
  type(ik_vec)                        :: rep_up,rep_dn
#else
  integer(ik),intent(inout)             :: imp_up(nimp),imp_dn(nimp)
  integer(ik)                        :: tmp_rep_up(nimp),tmp_rep_dn(nimp)
  integer(ik)                        :: new_dets_up(16),new_dets_dn(16)
  integer(ik)                        :: rep_up,rep_dn
#endif
  integer                            :: i,nsym_imp
  real(rk)                           :: phases(16)
  real(rk)                           :: norm,phase_w_rep
  integer                            :: num_distinct

  nsym_imp=0

  do i=1,nimp
    call generate_fourfold_k_configs(c4_map,reflection_map,z,p,imp_up(i),imp_dn(i),new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)
    if ((imp_up(i) .eq. rep_up) .and. (imp_dn(i) .eq. rep_dn)) then
      nsym_imp=nsym_imp+1
      tmp_rep_up(nsym_imp)=rep_up
      tmp_rep_dn(nsym_imp)=rep_dn
    endif
  enddo

  do i=1,nsym_imp
    imp_up(i)=tmp_rep_up(i)
    imp_dn(i)=tmp_rep_dn(i)
  enddo

  nimp = nsym_imp

  end subroutine symmetry_reduce_hubbardk_without_reallocating


  !=============================================================
  subroutine generate_sparse_ham_hubbardk(dets_up,dets_dn,H_indices,H_nonzero_elements,H_values,max_nonzero_elems,num_nonzero_elems,importance_sampling)
    !---------------------------------------------------------------------------
    ! Description : Compute the Hamiltonian in a given space and store it sparsely.
    !
    ! Created     : A Holmes, 7 Apr 2012
    ! Modified    : A Holmes, 15 Aug 2012. If max_nonzero_elems and num_nonzero_elems are present,
    !               then count the number of nonzero elements in H and store it as num_nonzero_elems.
    !               If that number exceeds max_nonzero_elems, return before constructing H.
    !---------------------------------------------------------------------------

  use types, only : i8b
  use common_psi_t, only : psi_g,psi_g_epsilon,psi_g_epsilon_inv,psi_t_connected_dets_up,psi_t_connected_dets_dn
  use mpi_routines, only : ncores

  implicit none

  integer(i8b),allocatable,intent(out) :: H_indices(:),H_nonzero_elements(:)
  real(rk),allocatable,intent(out) :: H_values(:)
  integer(i8b),optional,intent(in) :: max_nonzero_elems
  integer(i8b),optional,intent(out) :: num_nonzero_elems
  logical,optional,intent(in) :: importance_sampling
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
  type(ik_vec) :: tmp_up,tmp_dn
#else
  integer(ik),intent(in) :: dets_up(:),dets_dn(:)
  integer(ik) :: tmp_up,tmp_dn
#endif
  integer :: i,j,k,nnzero,n_det
  integer(i8b) :: nonzero_elements,isparse,itotal
  real(rk) :: matrix_element
  integer,allocatable :: H_rows(:),row_count(:)
  integer :: n_connected_dets
  integer :: connected_indices(max_connected_dets)
  integer :: numops
  logical :: brute_force = .false. ! whether or not to use brute force method (as opposed to binary search). For Hubbard, typical matrices are sparse enough that binary search usually wins.
  logical :: imp
  real(rk) :: guiding_wf_ratio

  call my_second(1,"generate_sparse_ham_hubbardk")
  call flush(6)

  imp = .false.
  if (present(importance_sampling))  imp=importance_sampling
  if (imp)  brute_force=.true. ! FIX THIS LATER; for now, only brute force method works with importance sampling.

  n_det=size(dets_up)

  allocate(H_nonzero_elements(n_det))
  H_nonzero_elements(:)=0

  allocate(row_count(n_det))

  ! Count number of nonzero elements of H
  nonzero_elements=0

  call find_connected_dets_hubbard_k(dets_up(1),dets_dn(1),n_connected_dets,connected_dets_up,connected_dets_dn)

 !call my_second(1,"generate_sparse_ham_hubbard")

  write (6,*) "N_det=",n_det,"N_con=",n_connected_dets
  call flush(6)

  if (brute_force) then

    do i = 1, n_det
      if (n_det>10) then
        if (mod(i,n_det/10)==0.or.i==n_det/100.or.i==(3*n_det)/100.or.i==(5*n_det)/100) then
          write (6,*) (100*i)/n_det,"% done"
        endif
      endif
      do j = 1, i
        if (is_connected_hubbard_fast(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j))) then
          nonzero_elements = nonzero_elements+2
          H_nonzero_elements(i)=H_nonzero_elements(i)+1
          H_nonzero_elements(j)=H_nonzero_elements(j)+1
          if (i==j) then
            nonzero_elements=nonzero_elements-1
            H_nonzero_elements(i)=H_nonzero_elements(i)-1
          endif
        endif
      enddo
    enddo

  else

   ! find connected dets and do binary search
   do i = 1, n_det
      if (n_det>10) then
        if (mod(i,n_det/10)==0.or.i==n_det/100.or.i==(3*n_det)/100.or.i==(5*n_det)/100) then
          write (6,*) (100*i)/n_det,"% done"
        endif
      endif
      call find_connected_dets_hubbard_k(dets_up(i),dets_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn)
      do k=1,n_connected_dets
        call binary_search(connected_dets_up(k),connected_dets_dn(k),dets_up(1:n_det),dets_dn(1:n_det),j)
        if (j>0.and.i>=j) then
          nonzero_elements = nonzero_elements+2
          H_nonzero_elements(i)=H_nonzero_elements(i)+1
          H_nonzero_elements(j)=H_nonzero_elements(j)+1
          if (i==j) then
            nonzero_elements=nonzero_elements-1
            H_nonzero_elements(i)=H_nonzero_elements(i)-1
          endif
        endif
      enddo
    enddo

  endif

  if (present(max_nonzero_elems).and.present(num_nonzero_elems)) then
   !call find_connected_dets_hubbard_k(dets_up(1),dets_dn(1),n_connected_dets,connected_dets_up,connected_dets_dn)
   !num_nonzero_elems=n_connected_dets*n_det
    num_nonzero_elems=nonzero_elements
    if (num_nonzero_elems>max_nonzero_elems)  return
  endif

  ! Generate H

  write(6,'(''allocating H_indices etc. arrays of size nonzero_elements='',i10,'' ='',es11.4,'' in generate_sparse_ham_hubbard'')') nonzero_elements, real(nonzero_elements)
  call flush(6)
  if (ncores>1.and.real(nonzero_elements)>7.e6_rk) then
    write(6,'(''Error: attempted to allocate too much memory ('',i10,'' elements) on each of too many cores ('',i4,''). Try generating the local energies and deterministic projector with a single-core run, then read in the elements during the multi-core run.'')') nonzero_elements, ncores
    stop "Attempted to allocate too much memory on each of too many cores!"
  endif

  allocate(H_indices(nonzero_elements))
  allocate(H_rows(nonzero_elements))
  allocate(H_values(nonzero_elements))
  H_indices(:)=0
  H_values(:)=0._rk

  !call my_second(2, 'allocating H_indices etc.')

  isparse = 0
  itotal = 0

  if (brute_force) then

    do i = 1, n_det
       if (imp) then
         call binary_search(dets_up(i),dets_dn(i),psi_t_connected_dets_up,psi_t_connected_dets_dn,j)
         if (j>0) then
           guiding_wf_ratio=psi_g(j)
         else
           guiding_wf_ratio=psi_g_epsilon
         endif
       endif
       if (i>1) then
         itotal=itotal+H_nonzero_elements(i-1)
         H_rows(isparse+1:itotal)=H_rows(isparse)
         isparse = itotal
       endif
       do j = 1, i
          if (is_connected_hubbard_fast(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j))) then
              if (space_sym) then
                tmp_up=dets_up(i)
                tmp_dn=dets_dn(i)
                call hamiltonian_hubbard_k_space_sym(tmp_up,tmp_dn, dets_up(j), dets_dn(j), matrix_element,nnzero)
              else
                call hamiltonian_hubbard_k(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element,.true.)
              endif
              isparse=isparse+1
              H_indices(isparse)=j
              H_rows(isparse)=i
              if (imp) then
                call binary_search(dets_up(j),dets_dn(j),psi_t_connected_dets_up,psi_t_connected_dets_dn,k)
                if (k>0) then
                  H_values(isparse)=matrix_element*guiding_wf_ratio/psi_g(k)
                else
                  H_values(isparse)=matrix_element*guiding_wf_ratio*psi_g_epsilon_inv
                endif
              else
                H_values(isparse)=matrix_element
              endif
              if (i.ne.j) then
                H_values(row_count(H_indices(isparse))) = H_values(isparse)
                H_indices(row_count(H_indices(isparse))) = H_rows(isparse)
                row_count(H_indices(isparse)) = row_count(H_indices(isparse)) + 1
              endif
              if (i==j)  row_count(i)=int(isparse)+1 ! First location of a nonzero element in the i'th row
          endif
       enddo
    enddo

  else

    ! find connected dets and do binary search
    do i = 1, n_det
       if (i>1) then
         itotal=itotal+H_nonzero_elements(i-1)
         H_rows(isparse+1:itotal)=H_rows(isparse)
         isparse = itotal
       endif
       call find_connected_dets_hubbard_k(dets_up(i),dets_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn)
       do k=1,n_connected_dets
          call binary_search(connected_dets_up(k),connected_dets_dn(k),dets_up(1:n_det),dets_dn(1:n_det),j)
          connected_indices(k) = j
       enddo
       ! sort connected_indices
       call sort(n_connected_dets,connected_indices(1:n_connected_dets),numops)
       do k=1,n_connected_dets
          j = connected_indices(n_connected_dets-k+1)
          if (j>0.and.i>=j) then
             if (space_sym) then
               tmp_up = dets_up(i)
               tmp_dn = dets_dn(i)
               call hamiltonian_hubbard_k_space_sym(tmp_up,tmp_dn, dets_up(j), dets_dn(j), matrix_element,nnzero)
             else
               call hamiltonian_hubbard_k(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element,.true.)
             endif
             isparse=isparse+1
             H_indices(isparse)=j
             H_rows(isparse)=i
             H_values(isparse)=matrix_element
             if (i.ne.j) then
               H_values(row_count(H_indices(isparse))) = H_values(isparse)
               H_indices(row_count(H_indices(isparse))) = H_rows(isparse)
               row_count(H_indices(isparse)) = row_count(H_indices(isparse)) + 1
             endif
             if (i==j)  row_count(i)=int(isparse)+1 ! First location of a nonzero element in the i'th row
          endif
       enddo
    enddo

  endif

  call my_second(2, 'generating sparse Hubbardk Hamiltonian')

  end subroutine generate_sparse_ham_hubbardk


  !=============================================================
  subroutine generate_sparse_ham_hubbardk_upper_triangular(dets_up,dets_dn,H_indices,H_nonzero_elements,H_values,hf_to_psit,max_nonzero_elems,num_nonzero_elems,importance_sampling)
    !---------------------------------------------------------------------------
    ! Description : Compute the Hamiltonian in a given space and store it sparsely, as an upper triangular matrix.
    !               Not only does this save us a factor of 2 in storage, but it also makes it easier to work with
    !               non-square deterministic projectors.
    !
    ! Created     : A Holmes, 14 Dec 2012
    ! Modified    : A Holmes, 17 Dec 2012. Input variable hf_to_psit replaces the first state (HF) with psi_trial.
    !---------------------------------------------------------------------------

  use types, only : i8b
  use common_psi_t, only : psi_g,psi_g_epsilon,psi_g_epsilon_inv,psi_t_connected_dets_up,psi_t_connected_dets_dn
  use mpi_routines, only : ncores

  implicit none

  integer(i8b),allocatable,intent(out) :: H_indices(:),H_nonzero_elements(:)
  real(rk),allocatable,intent(out) :: H_values(:)
  integer(i8b),optional,intent(in) :: max_nonzero_elems
  integer(i8b),optional,intent(out) :: num_nonzero_elems
  logical,optional,intent(in) :: importance_sampling
  logical,intent(in) :: hf_to_psit
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
  type(ik_vec) :: tmp_up,tmp_dn
#else
  integer(ik),intent(in) :: dets_up(:),dets_dn(:)
  integer(ik) :: tmp_up,tmp_dn
#endif
  integer :: i,j,k,nnzero,n_det
  integer(i8b) :: nonzero_elements,isparse,itotal
  real(rk) :: matrix_element
  integer :: n_connected_dets
  integer :: connected_indices(max_connected_dets)
  integer :: numops
  logical :: brute_force = .false. ! whether or not to use brute force method (as opposed to binary search). For Hubbard, typical matrices are sparse enough that binary search usually wins.
  logical :: imp
  real(rk) :: guiding_wf_ratio

  call my_second(1,"generate_sparse_ham_hubbardk_upper_triangular")
  call flush(6)

  imp = .false.
  if (present(importance_sampling))  imp=importance_sampling
  if (imp)  brute_force=.true. ! FIX THIS LATER; for now, only brute force method works with importance sampling.

  if (brute_force) then
    write (6,*) "Generating matrix using brute force method"
  else
    write (6,*) "Generating matrix using binary search method"
  endif

  n_det=size(dets_up)

  allocate(H_nonzero_elements(n_det))
  H_nonzero_elements(:)=0

  !allocate(row_count(n_det))

  ! Count number of nonzero elements of H
  nonzero_elements=0

  call find_connected_dets_hubbard_k(dets_up(1),dets_dn(1),n_connected_dets,connected_dets_up,connected_dets_dn)

 !call my_second(1,"generate_sparse_ham_hubbard")

  write (6,*) "N_det=",n_det,"N_con=",n_connected_dets
  call flush(6)

  if (brute_force) then

    do i = 1, n_det
      if (n_det>10) then
        if (mod(i,n_det/10)==0.or.i==n_det/100.or.i==(3*n_det)/100.or.i==(5*n_det)/100) then
          write (6,*) (100*i)/n_det,"% done"
        endif
      endif
      do j = i, n_det
        if (.not.hf_to_psit.or.(i.ne.1.and.j.ne.1).or.(i==1.and.j==1)) then
          if (space_sym) then
            if (is_connected_hubbard_fast_sym(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j))) then
              nonzero_elements = nonzero_elements+1
              H_nonzero_elements(i)=H_nonzero_elements(i)+1
            endif
          else
            if (is_connected_hubbard_fast(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j))) then
              nonzero_elements = nonzero_elements+1
              H_nonzero_elements(i)=H_nonzero_elements(i)+1
            endif
          endif
        endif
      enddo
    enddo

  else

   ! find connected dets and do binary search
   do i = 1, n_det
      if (n_det>10) then
        if (mod(i,n_det/10)==0.or.i==n_det/100.or.i==(3*n_det)/100.or.i==(5*n_det)/100) then
          write (6,*) (100*i)/n_det,"% done"
        endif
      endif
      if (i==1) then
        if (hf_to_psit) then
          nonzero_elements = nonzero_elements + 1
          H_nonzero_elements(1) = H_nonzero_elements(1) + 1
          cycle
        endif
      endif
      call find_connected_dets_hubbard_k(dets_up(i),dets_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn)
      do k=1,n_connected_dets
        call binary_search(connected_dets_up(k),connected_dets_dn(k),dets_up(1:n_det),dets_dn(1:n_det),j)
        if (j>0.and.i<=j) then
          nonzero_elements = nonzero_elements+1
          H_nonzero_elements(i)=H_nonzero_elements(i)+1
        endif
      enddo
    enddo

  endif

  if (present(max_nonzero_elems).and.present(num_nonzero_elems)) then
   !call find_connected_dets_hubbard_k(dets_up(1),dets_dn(1),n_connected_dets,connected_dets_up,connected_dets_dn)
   !num_nonzero_elems=n_connected_dets*n_det
    num_nonzero_elems=nonzero_elements
    if (num_nonzero_elems>max_nonzero_elems)  return
  endif

  ! Generate H

  write(6,'(''allocating H_indices etc. arrays of size nonzero_elements='',i10,'' ='',es11.4,'' in generate_sparse_ham_hubbard'')') nonzero_elements, real(nonzero_elements)
  call flush(6)
  if (ncores>1.and.real(nonzero_elements)>7.e6_rk) then
    write(6,'(''Error: attempted to allocate too much memory ('',i10,'' elements) on each of too many cores ('',i4,''). Try generating the local energies and deterministic projector with a single-core run, then read in the elements during the multi-core run.'')') nonzero_elements, ncores
    stop "Attempted to allocate too much memory on each of too many cores!"
  endif

  allocate(H_indices(nonzero_elements))
  allocate(H_values(nonzero_elements))
  H_indices(:)=0
  H_values(:)=0._rk

  !call my_second(2, 'allocating H_indices etc.')

  isparse = 0
  itotal = 0

  if (brute_force) then

    do i = 1, n_det
       if (imp) then
         call binary_search(dets_up(i),dets_dn(i),psi_t_connected_dets_up,psi_t_connected_dets_dn,j)
         if (j>0) then
           guiding_wf_ratio=psi_g(j)
         else
           guiding_wf_ratio=psi_g_epsilon
         endif
       endif
       do j = i, n_det
         if (.not.hf_to_psit.or.(i.ne.1.and.j.ne.1)) then
           if (space_sym) then
             if (is_connected_hubbard_fast_sym(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j))) then
               tmp_up=dets_up(i)
               tmp_dn=dets_dn(i)
               call hamiltonian_hubbard_k_space_sym(tmp_up,tmp_dn, dets_up(j), dets_dn(j), matrix_element,nnzero)
             else
               cycle
             endif
           else
             if (is_connected_hubbard_fast(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j))) then
               call hamiltonian_hubbard_k(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element,.true.)
             else
               cycle
             endif
             isparse=isparse+1
             H_indices(isparse)=j
             if (imp) then
               call binary_search(dets_up(j),dets_dn(j),psi_t_connected_dets_up,psi_t_connected_dets_dn,k)
               if (k>0) then
                 H_values(isparse)=matrix_element*guiding_wf_ratio/psi_g(k)
               else
                 H_values(isparse)=matrix_element*guiding_wf_ratio*psi_g_epsilon_inv
               endif
             else
               H_values(isparse)=matrix_element
             endif
           endif
         endif
         if (hf_to_psit.and.i==1.and.j==1) then
           isparse=isparse+1
           H_indices(isparse)=j
           H_values(isparse)=0._rk
         endif
       enddo
    enddo

  else

    ! find connected dets and do binary search
    do i = 1, n_det
      if (i==1) then
        if (hf_to_psit) then
          H_indices(1) = 1
          H_values(1) = 0._rk
          isparse = 1
          cycle
        endif
      endif
      call find_connected_dets_hubbard_k(dets_up(i),dets_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn)
      do k=1,n_connected_dets
         call binary_search(connected_dets_up(k),connected_dets_dn(k),dets_up(1:n_det),dets_dn(1:n_det),j)
         connected_indices(k) = j
      enddo
      ! sort connected_indices
      call sort(n_connected_dets,connected_indices(1:n_connected_dets),numops)
      do k=1,n_connected_dets
         j = connected_indices(n_connected_dets-k+1)
         if (j>0.and.i<=j) then
            if (space_sym) then
              tmp_up = dets_up(i)
              tmp_dn = dets_dn(i)
              call hamiltonian_hubbard_k_space_sym(tmp_up,tmp_dn, dets_up(j), dets_dn(j), matrix_element,nnzero)
            else
              call hamiltonian_hubbard_k(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element,.true.)
            endif
            isparse=isparse+1
            H_indices(isparse)=j
            H_values(isparse)=matrix_element
         endif
      enddo
    enddo

  endif

  call my_second(2, 'generating sparse Hubbardk Hamiltonian')

  end subroutine generate_sparse_ham_hubbardk_upper_triangular


! !===========================================================================
  function is_connected_hubbard_fast(up1,dn1,up2,dn2)
    ! Description:  Checks whether dets 1 and 2 are connected using intelligent bit counting implemented by Frank in tools.f90. Assumes det 1 and det 2 are already in same symmetry sector (e.g., same total momentum)
    ! Created    :  A Holmes, 24 April 2012

    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: up1,dn1,up2,dn2
    type(ik_vec) :: diff
#else
    integer(ik), intent(in) :: up1,dn1,up2,dn2
    integer(ik) :: diff
#endif
    logical :: is_connected_hubbard_fast

    integer :: n_set

    if (up1==up2.and.dn1==dn2) then
      is_connected_hubbard_fast=.true.
      return
    endif

    n_set = 0
    diff = iand(up1,not(up2))
    do while (diff .ne. 0.and.n_set<2)
      n_set = n_set + 1
      diff = iand(diff, diff - 1)
    enddo

    if (n_set.ne.1) then
      is_connected_hubbard_fast=.false.
      return
    endif

    n_set = 0
    diff = iand(dn1,not(dn2))
    do while (diff .ne. 0.and.n_set<2)
      n_set = n_set + 1
      diff = iand(diff, diff - 1)
    enddo

    if (n_set==1) then
      is_connected_hubbard_fast=.true.
    else
      is_connected_hubbard_fast=.false.
    endif

  end function is_connected_hubbard_fast

  !===========================================================================
  function is_connected_hubbard_fast_sym(up1,dn1,up2,dn2)
    ! Description:  Checks whether dets 1 and symmetry related configurations of 2 are connected using intelligent bit counting implemented by Frank in tools.f90.
    !               Assumes det 1 and det 2 are already in same symmetry sector (e.g., same total momentum)
    ! Created    :  A Holmes, 24 April 2012

    implicit none
    logical :: is_connected_hubbard_fast_sym
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: up1,dn1,up2,dn2
    type(ik_vec)              :: up3,dn3,diff
    type(ik_vec)              :: new_dets_up(16),new_dets_dn(16)            ! Symmetry related determinants
    type(ik_vec)              :: rep_up,rep_dn                              ! Representatives
#else
    integer(ik), intent(in) :: up1,dn1,up2,dn2
    integer(ik)              :: up3,dn3,diff
    integer(ik)              :: new_dets_up(16),new_dets_dn(16)            ! Symmetry related determinants
    integer(ik)              :: rep_up,rep_dn                              ! Representatives
#endif
    integer                  :: n_set
    real(rk)                 :: phases(16)                                 ! Phase factors in the linear combination
    real(rk)                 :: phase_w_rep,norm                           ! Norm
    integer                  :: num_distinct                               ! Number of distinct configurations
    integer                  :: i

    call generate_fourfold_k_configs(c4_map,reflection_map,z,p,up2,dn2,new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)

    do i=1,16
        up3=new_dets_up(i)
        dn3=new_dets_dn(i)
        if (up1==up3.and.dn1==dn3) then
            is_connected_hubbard_fast_sym=.true.
            return
        endif
        n_set = 0
        diff = iand(up1,not(up3))
        do while (diff .ne. 0.and.n_set<2)
          n_set = n_set + 1
          diff = iand(diff, diff - 1)
        enddo

        if (n_set.ne.1) then
          is_connected_hubbard_fast_sym=.false.
          cycle
        endif

        n_set = 0
        diff = iand(dn1,not(dn3))
        do while (diff .ne. 0.and.n_set<2)
          n_set = n_set + 1
          diff = iand(diff, diff - 1)
        enddo

        if (n_set==1) then
          is_connected_hubbard_fast_sym=.true.
          return
        else
          is_connected_hubbard_fast_sym=.false.
        endif
     enddo

  end function is_connected_hubbard_fast_sym


!=====================================================================================================================
 subroutine matrix_lanczos_on_the_fly(dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,highest_eigenvalue,second_lowest_eigenvalue,short_run_size,initial_vector)
 ! Created by  : Hitesh Changlani, March 2012
 ! Modified by : A Holmes, 15 Aug 2012. Modified input to allow for computing H on the fly.
 !               A Holmes, 12 Feb 2013. When doing matrix_lanczos_on_the_fly, first seed with a shorter Lanczos run in a truncated space.
!=====================================================================================================================
 implicit none

 ! Dummy
 real(rk),intent(out) :: lowest_eigenvector(:)
 real(rk),intent(out) :: lowest_eigenvalue
 real(rk),intent(out),optional :: highest_eigenvalue,second_lowest_eigenvalue
 integer,intent(in),optional :: short_run_size ! if present, seed long Lanczos run with a short Lanczos run on this many of the first dets.
 real(rk),intent(in),optional :: initial_vector(:)
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
 type(ik_vec),allocatable    :: temp_i16_up(:),temp_i16_dn(:)!,temp_i_2(:)
#else
 integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
 integer(ik),allocatable    :: temp_i16_up(:),temp_i16_dn(:)!,temp_i_2(:)
#endif

 ! Local
 integer                    :: i,it
!real(rk)                   :: rannyu
 real(rk)                   :: energy_shift
 real(rk)                   :: norm,norm_inv
 real(rk),allocatable       :: w(:),v(:,:)
 real(rk),allocatable       :: alphas(:),betas(:)
 integer                    :: iterations
 real(rk)                   :: lowest_eigenvalue_prev
 logical                    :: converged=.false.
 integer                    :: len_work,info
 real(rk),allocatable       :: work(:),eigenvalues(:),tridiag(:,:)
 real(rk)                   :: epsilon=1.0e-6  ! for Lanczos
 integer                    :: n
 real(rk), allocatable :: H_values(:)
 integer(i8b), allocatable :: H_indices(:),H_nonzero_elements(:)
 real(rk),allocatable :: tmp_lowest_eigenvector(:)
 integer,allocatable        :: iorder(:),temp_i_2(:)

  n=size(dets_up)

  iterations=30          ! User option
  iterations=min(n,iterations)
  allocate (v(n,iterations+1))
  allocate (w(n))
  allocate(alphas(iterations+1))
  allocate(betas(iterations+1))
  w(:)=0._rk

  if (present(short_run_size)) then
    write (6,*) "In matrix_lanczos_on_the_fly, performing short run on the first",short_run_size,"dets"

    ! sort by label so that a binary search can be performed
    n=short_run_size
    allocate(temp_i16_up((n+1)/2))
    allocate(temp_i16_dn((n+1)/2))
    allocate(temp_i_2((n+1)/2))
    allocate(iorder(n))
    call merge_sort2_up_dn(dets_up(1:n),dets_dn(1:n), iorder, n, temp_i16_up, temp_i16_dn, temp_i_2)
    deallocate(iorder)
    deallocate(temp_i16_up)
    deallocate(temp_i16_dn)
    deallocate(temp_i_2)

    call generate_sparse_ham_hubbardk_upper_triangular(dets_up(1:short_run_size),dets_dn(1:short_run_size),H_indices,H_nonzero_elements,H_values,hf_to_psit=.false.) ! When diagonalizing, we don't want to transform the projector first (hence hf_to_psit=false)
    allocate(tmp_lowest_eigenvector(size(dets_up)))
    tmp_lowest_eigenvector(:)=0._rk
    call matrix_lanczos(short_run_size,tmp_lowest_eigenvector(1:short_run_size),lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue)
    call my_second(2,'short Lanczos run')

    ! sort by label so that a binary search can be performed
    n=size(dets_up)
    allocate(temp_i16_up((n+1)/2))
    allocate(temp_i16_dn((n+1)/2))
    allocate(temp_i_2((n+1)/2))
    allocate(iorder(n) )
    do i=1,n
      iorder(i)=i
    enddo
    call merge_sort2_up_dn(dets_up,dets_dn, iorder, n, temp_i16_up, temp_i16_dn, temp_i_2)
    deallocate(temp_i16_up)
    deallocate(temp_i16_dn)
    deallocate(temp_i_2)

    v(1:n,1)=tmp_lowest_eigenvector(iorder(1:n))
    deallocate(iorder)
    deallocate(tmp_lowest_eigenvector)

  elseif (present(initial_vector)) then
    v(:,1)=initial_vector(:)
  else
    v(:,1)=0._rk
    v(1,1)=1._rk
  endif

  energy_shift=0._rk
  betas(1)=0._rk
  allocate (tridiag(iterations,iterations))

  allocate(eigenvalues(iterations))
  len_work = 3*iterations-1
  allocate(work(len_work))

  converged=.false.

  if (n>1) then
      do it=1,iterations
         call apply_H_on_the_fly(dets_up,dets_dn,v(:,it),w(:))
         if (it .gt. 1) w(:)=w(:)-betas(it)*v(:,it-1)
         alphas(it)=dot_product(w,v(:,it))
         w(:)=w(:)-alphas(it)*v(:,it)
         norm=dot_product(w,w)
         if (norm<(1.e-12_rk))  converged=.true.
         betas(it+1)=norm**(0.5_rk)
         norm_inv=1._rk/betas(it+1)
         v(:,it+1)=w(:)*norm_inv
         w(:)=v(:,it+1)
         do i=1,it        ! Reorthogonalization
             norm=dot_product(v(:,it+1),v(:,i))
             !write (6,*) "q=",norm
             call flush(6)
             w(:)=w(:)-norm*v(:,i)
         enddo
         v(:,it+1)=w(:)
         w(:)=0._rk
         norm=dot_product(v(:,it+1),v(:,it+1))
         norm_inv=1._rk/(norm**(0.5_rk))
         v(:,it+1)=v(:,it+1)*norm_inv
         tridiag(:,:)=0._rk

         eigenvalues(:)=0._rk

         do i=1,it
             tridiag(i,i)=alphas(i)
             if (i<it) then
                 tridiag(i,i+1)=betas(i+1)
                 tridiag(i+1,i)=betas(i+1)
             endif
         enddo

         !diagonalize with lapack routine
         len_work = 3*it-1

         call dsyev('V', 'U', it, tridiag(1:it,1:it), it, eigenvalues, work, len_work, info)

         lowest_eigenvalue=eigenvalues(1)
         if (present(highest_eigenvalue))  highest_eigenvalue=eigenvalues(it)
         if (present(second_lowest_eigenvalue).and.it>1)  second_lowest_eigenvalue=eigenvalues(2)
         if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<epsilon) then
             converged=.true.
             exit
         else
             lowest_eigenvalue_prev=lowest_eigenvalue
             write(6,'(''Iteration, Eigenvalue='',i3,f15.9)') it, lowest_eigenvalue
             call flush(6)
         endif
         if (converged)  exit

      enddo

      it=min(it,iterations)
      write(6,'(''n, Lowest eigenvalue ='',i10, f16.10)') n, lowest_eigenvalue
      call flush(6)

      v(:,1)=matmul(v(:,1:it),tridiag(1:it,1))

      if (allocated(eigenvalues)) deallocate(eigenvalues)
      if (allocated(tridiag))     deallocate(tridiag)
      if (allocated(work))        deallocate(work)

  else
      lowest_eigenvalue = eigenvalues(1)
  endif

   lowest_eigenvector(:)=v(:,1)

   call my_second(2,'matrix_lanczos_on_the_fly') ; call flush(6)

  end subroutine matrix_lanczos_on_the_fly

!=====================================================================================================================
  subroutine apply_H_on_the_fly(dets_up,dets_dn,v1,v2)
    ! A Holmes, 15 Aug 2012
    ! v2 = H * v1
    ! Assumes dets_up,dn already sorted by label
!=====================================================================================================================

    implicit none

    real(rk),intent(in) :: v1(:)
    real(rk),intent(out) :: v2(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec) :: tmp_up,tmp_dn
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik) :: tmp_up,tmp_dn
#endif

    integer :: i,j,k,n_connected_dets
    real(rk) :: matrix_element
    integer :: nnzero

    v2(:)=0._rk
    do i=1,size(dets_up)
      ! generate connections = connected_dets_up,dn
      call find_connected_dets_hubbard_k(dets_up(i),dets_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn,nsites)
      ! loop over connections: binary search the dets list for each connection.
      do j=1,n_connected_dets
        call binary_search(connected_dets_up(j),connected_dets_dn(j),dets_up,dets_dn,k)
        if (k>0) then
          if (space_sym) then
            tmp_up=dets_up(i)
            tmp_dn=dets_dn(i)
            call hamiltonian_hubbard_k_space_sym(tmp_up,tmp_dn,dets_up(k),dets_dn(k),matrix_element,nnzero)
          else
            call hamiltonian_hubbard_k(dets_up(i),dets_dn(i),dets_up(k),dets_dn(k),matrix_element,.true.)
          endif
          v2(k) = v2(k) + matrix_element * v1(i)
        endif
      enddo
    enddo

  end subroutine apply_H_on_the_fly


!=====================================================================================================================
  subroutine generate_hfb(n_max,dets_up,dets_dn)
    ! Generates HFB spectrum using HF as reference, i.e., generates all HFB states (where every state is composed of pairs of
    ! up and dn electrons with 0 total momentum per pair).
    ! Not presently being used.
    ! A Holmes, 29 Nov 2012.
!=====================================================================================================================
   implicit none

   integer,optional,intent(in) :: n_max
#ifdef NUM_ORBITALS_GT_127
   type(ik_vec),allocatable,intent(out) :: dets_up(:),dets_dn(:)
   type(ik_vec) :: tmp_det
#else
   integer(ik),allocatable,intent(out) :: dets_up(:),dets_dn(:)
   integer(ik) :: tmp_det
#endif
   integer :: n,i,j,old,new

   if(n_choose_k(nsites,nup).le.2147483647) then !  2^31-1 = 2147483647
     n = min(n_max,int(n_choose_k(nsites,nup),i4b))
   else
     stop 'n_choose_k(nsites,nup) > 2^31-1'
   endif
   allocate(dets_up(n))
   allocate(dets_dn(n))

   dets_dn(:) = 0_ik
   dets_up(1) = ibset(0_ik,nup)-1_ik
   dets_dn(1) = dets_dn(1)

   do i=2,n
     dets_up(i) = dets_up(i-1)
     call twiddle_determinant(dets_up(i))
     tmp_det = dets_up(i)
     do j=1,nup
       old = trailz(tmp_det)+1
       new = kmap(old,mod(-2*k_vectors(1,old),2*l_x),mod(-2*k_vectors(2,old),2*l_y))
       dets_dn(i) = ibset(dets_dn(i),new-1)
       tmp_det = ibclr(tmp_det,old-1)
     enddo
   enddo

  end subroutine generate_hfb

  !=============================================================
  subroutine twiddle_determinant(det_spin)
    !---------------------------------------------------------------------------
    ! Description : Return the next lexographically ordered bit pattern after det_spin.
    !               Note det_spin is a determinant of either up or dn spin.
    !               Algorithm from Bit Twiddling Hacks website.
    !
    ! Created     : F. Petruzielo, 14 Mar 2011
    !---------------------------------------------------------------------------
    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(inout) :: det_spin
    type(ik_vec) :: temp
#else
    integer(ik), intent(inout) :: det_spin
    integer(ik) :: temp
#endif

    temp = ior(det_spin, det_spin - 1) + 1
    det_spin = ior(temp,ishft(iand(temp, -temp) / iand(det_spin,-det_spin),-1) - 1)

  end subroutine twiddle_determinant

    subroutine print_excitation_levels_and_wts_hubbard_sym(n_det,dets_up,dets_dn,lowest_eigenvector,norb,orbital_symmetries)
    !=====================================================================================================================
    ! Description   : Calculate number of dets with various excitation levels and the sum of their squared wts.
    ! Author        : Cyrus Umrigar, 3 Jul 2012. (moved to this subroutine from chemistry.f90 by A Holmes)
    ! ---------------------------------------------------------------------------------------------------------------------
        use more_tools,only  : generate_fourfold_k_configs
        implicit none

        integer,intent(in) :: n_det
        real(rk),intent(in) :: lowest_eigenvector(:)
        integer,intent(in) :: norb
        integer,optional,intent(in) :: orbital_symmetries(:)
#ifdef NUM_ORBITALS_GT_127
        type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
        type(ik_vec)              :: new_dets_up(16),new_dets_dn(16)   ! Symmetry related determinants
        type(ik_vec)              :: rep_up,rep_dn
        type(ik_vec)              :: tmp_distinct_up(16),tmp_distinct_dn(16)
#else
        integer(ik),intent(in) :: dets_up(:),dets_dn(:)
        integer(ik)              :: new_dets_up(16),new_dets_dn(16)   ! Symmetry related determinants
        integer(ik)              :: rep_up,rep_dn
        integer(ik)              :: tmp_distinct_up(16),tmp_distinct_dn(16)
#endif

        integer :: i,j,iorb
        real(rk), allocatable :: sum_w2_orb(:)
        integer :: excite_level
        integer, parameter :: max_excite_level=20 ! Used for printing statistics only
        integer :: popul_excit_level(0:max_excite_level)
        real(rk) ::sum_w2_excit_level(0:max_excite_level)
        real(rk) :: scal

        ! Local
        real(rk)                 :: phases(16)                        ! Phase factors in the linear combination
        real(rk)                 :: phase_w_rep,norm                  ! Norm
        integer                  :: num_distinct

        ! Local
        integer                  :: location,counter
        logical                  :: entry_present

        counter=0
        do i=1,16
          tmp_distinct_up(i)=0_ik
          tmp_distinct_dn(i)=0_ik
        enddo

        popul_excit_level(1:max_excite_level)=0
        sum_w2_excit_level(1:max_excite_level)=0
        do i=1,n_det
          excite_level=min(popcnt(iand(dets_up(i),not(dets_up(1))))+popcnt(iand(dets_dn(i),not(dets_dn(1)))),max_excite_level)
          popul_excit_level(excite_level)=popul_excit_level(excite_level)+1
          sum_w2_excit_level(excite_level)=sum_w2_excit_level(excite_level)+lowest_eigenvector(i)**2
        enddo
        write(6,'(/,''Excit_level  # dets  sum_w2_excit_level'')')
        do excite_level=0,max_excite_level
          if(sum_w2_excit_level(excite_level).ne.0_rk) write(6,'(i6,i10,f14.8)') excite_level, popul_excit_level(excite_level), sum_w2_excit_level(excite_level)
        enddo
    ! Calculate sum of squared wts on each orbital
        allocate(sum_w2_orb(norb))
        sum_w2_orb(1:norb)=0
        !write (6,*) "n_det=",n_det
        !do i = 1, min(n_det,100)
        !    write(6,'(2i22,f18.14)') dets_up(i), dets_dn(i), lowest_eigenvector(i)
        !enddo
        !call flush(6)
        do i=1,n_det
            ! Total occupation being measured - instead of up and down separately
                call generate_fourfold_k_configs(c4_map,reflection_map,z,p,dets_up(i),dets_dn(i),new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)
                scal=norm/sqrt(16._rk)
                do j=1,16
                  tmp_distinct_up(j)=0_ik
                  tmp_distinct_dn(j)=0_ik
                enddo
                counter=0
                if (abs(norm) .gt. 1.0e-10 ) then
                    do j=1,16
                        if (counter .eq. num_distinct) exit
                            call is_in_list_mod(new_dets_up(j),new_dets_dn(j),tmp_distinct_up,tmp_distinct_dn,counter,entry_present,location)
                            if (entry_present .eqv. .false.) then
                                counter=counter+1
                                tmp_distinct_up(counter)=new_dets_up(j)
                                tmp_distinct_dn(counter)=new_dets_dn(j)
                                !write(6,'(2i22,f18.14)') new_dets_up(j), new_dets_dn(j), lowest_eigenvector(j)*scal
                                do iorb=1,norb
                                    if(btest(new_dets_up(j),iorb-1)) sum_w2_orb(iorb)=sum_w2_orb(iorb)+(scal*lowest_eigenvector(i))**2
                                    if(btest(new_dets_dn(j),iorb-1)) sum_w2_orb(iorb)=sum_w2_orb(iorb)+(scal*lowest_eigenvector(i))**2
                                enddo
                            endif
                    enddo
                endif
        enddo
        if (present(orbital_symmetries)) then
          write(6,'(/,''Orbital orb_sym sum_w2_orb'')')
          do iorb=1,norb
            if(sum_w2_orb(iorb).ne.0_rk) write(6,'(2i6,f12.6)') iorb, orbital_symmetries(iorb), sum_w2_orb(iorb)
          enddo
        else
          write(6,'(/,''Orbital sum_w2_orb'')')
          do iorb=1,norb
            if(sum_w2_orb(iorb).ne.0_rk) write(6,'(i6,f12.6)') iorb, sum_w2_orb(iorb)
          enddo
        endif
        deallocate(sum_w2_orb)
        call flush(6)

    end subroutine print_excitation_levels_and_wts_hubbard_sym


end module hubbard
