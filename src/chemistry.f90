module chemistry
  use types, only: rk, ik, ik_vec, i8b, rs_absH, int_vec
#ifdef NUM_ORBITALS_GT_127
  use overload
#endif
  use generic_sort, only : sort
  use tools, only : random_int,merge_sort2_up_dn,merge_original_with_spawned3
  use more_tools, only : print_real_matrix,matrix_lanczos
  use common_run, only : ipr, proposal_method, tau, run_type, max_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, connected_matrix_elements_fn, connected_diag_elems_info
  use common_ham, only : nelec, nup, ndn, diagonalize_ham,norb, n_core_orb,max_energy
  use common_psi_t, only : ndet_psi_t, dets_up_psi_t, dets_dn_psi_t, cdet_psi_t, e_trial,trial_wf_iters, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf
  use common_selected_ci, only : det_sel_iters, norb_det_sel, n_sym_uniq_det_det_sel,lanczos_iters,lanczos_initiators,lanczos_truncate
  use common_imp, only: imp_iters,norb_imp, n_imp,imp_up,imp_dn

  implicit none
  save
  public ::  read_chem, system_setup_chem, off_diagonal_move_chem, find_connected_dets_chem, find_important_connected_dets_chem, energy_pieces_chem, generate_sparse_ham_chem, &  !subroutines
           & generate_sparse_ham_chem_upper_triangular, is_connected_chem,symmetry_reduce_and_replace_chem, &
           & off_diagonal_move_chem_cauchySchwarz,setup_orb_by_symm, off_diagonal_move_chem_efficient_heatbath, &
           & get_new_diag_elem, is_occupied, occ_up, occ_dn, &
           & diagonalize_sparse_hamiltonian_chem_excited, integral_value, rs_absH, orb_order, orb_order_inv

  !set up pointer for the integrals array
  real(rk), dimension(:),pointer,contiguous      :: integrals
  real(rk)                                       :: nuclear_nuclear_energy,mp2_energy_correction
  real(rk), allocatable                          :: orbital_energies(:)
  integer, allocatable                           :: orbital_symmetries(:)
  integer, allocatable                           :: product_table_elems(:, :)
!*** AR: for now use larger length, make sure trim is always used
! character(len=3)                               :: point_group
  character(16)                                  :: point_group

  integer                                        :: spatial_symmetry_wf
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec)                                   :: last_det_up, last_det_dn
#else
  integer(ik)                                    :: last_det_up, last_det_dn
#endif
  real(rk)                                       :: last_energy
  real(rk)                                       :: orb_hf_energy ! used to store hf_energy appximately calculated
  integer                                        :: counter_two_body  ! used for improved two-body energy routine
  integer                                        :: z
  logical                                        :: time_sym
  real(rk)                                       :: sqrt2,sqrt2inv
  integer                                        :: tmp_ctr
  integer, allocatable                           :: filled_up(:), empty_up(:)
  integer, allocatable                           :: filled_dn(:), empty_dn(:)
  integer                                        :: init_multiplier,trunc_multiplier
  integer,allocatable                            :: orb_order(:)
  integer,allocatable                            :: orb_order_inv(:)
  integer                                        :: d_infinity_h ! 1 if using d_infinity_h; 0 otherwise
  integer,allocatable                            :: hf_up_orbs(:),hf_dn_orbs(:)

  real(rk),allocatable :: one_orbital_probabilities(:)
  real(rk),allocatable :: two_orbital_probabilities(:,:)

  real(rk),allocatable :: three_orbital_probabilities_same_spin(:,:,:),three_orbital_probabilities_opposite_spin(:,:,:)
  real(rk) :: c_three_orbital_probabilities_same_spin,c_three_orbital_probabilities_opposite_spin
  integer,allocatable :: J_three_orbital_probabilities_same_spin(:,:,:),J_three_orbital_probabilities_opposite_spin(:,:,:)
  real(rk),allocatable :: q_three_orbital_probabilities_same_spin(:,:,:),q_three_orbital_probabilities_opposite_spin(:,:,:)

  real,pointer,dimension(:),contiguous :: four_orbital_probabilities_same_spin
  real(rk) :: c_four_orbital_probabilities_same_spin !Should we make this real to be consistent?

  integer,pointer,dimension(:),contiguous :: J_four_orbital_probabilities_same_spin
  real,pointer,dimension(:),contiguous :: q_four_orbital_probabilities_same_spin

  real,pointer,dimension(:),contiguous :: four_orbital_probabilities_opposite_spin
  real(rk) :: c_four_orbital_probabilities_opposite_spin !Should we make this real to be consistent?

  integer,pointer,dimension(:),contiguous :: J_four_orbital_probabilities_opposite_spin
  real,pointer,dimension(:),contiguous :: q_four_orbital_probabilities_opposite_spin

  real(rk),allocatable :: Htot_same(:,:)
  real(rk),allocatable :: Htot_opposite(:,:,:)

  integer,allocatable :: occ_up(:),occ_dn(:)
  integer,allocatable :: elecs(:)
  real(rk),allocatable :: e1_prob(:),c_e1_prob(:)
  real(rk),allocatable :: e2_prob(:),c_e2_prob(:)
  integer,allocatable :: pairs_e1(:),pairs_e2(:)
  real(rk),allocatable :: pairs_prob(:),c_pairs_prob(:)
  integer,allocatable :: J_pairs_prob(:)
  real(rk),allocatable :: q_pairs_prob(:)
  real(rk) :: proposal_prob_single
  real(rk) :: proposal_prob_double
  integer :: n_pairs ! Number of pairs of electrons = nelec*(nelec-1)/2
  integer,allocatable :: tmp1(:),tmp2(:)

  logical                                        :: uniform_sampling
  integer                                        :: n_group_elements
  integer,allocatable                            :: which_orb_by_sym(:,:)
  real(rk),allocatable                           :: sym_sum_cs_sqrt(:,:)
  real(rk),allocatable                           :: cs_sqrt_orb(:)
  real(rk),allocatable                           :: cs_sqrt_prime(:),cs_sqrt_prime_spin(:)
  real(rk),allocatable                           :: sqrt_integrals(:,:)
  real(rk)                                       :: sum_cs_sqrt_prime
  integer,allocatable                            :: num_orb_by_sym(:)
  integer,allocatable                            :: occ_orb(:)
  integer,allocatable                            :: occ_orb_up(:),occ_orb_dn(:)
  integer,allocatable                            :: num_occ_up_by_sym(:),num_occ_dn_by_sym(:)
  integer,allocatable                            :: occ_orb_up_by_sym(:,:),occ_orb_dn_by_sym(:,:)
  integer,allocatable                            :: occ_orb_up_by_sym_copy(:,:),occ_orb_dn_by_sym_copy(:,:)
  integer,allocatable :: combine_2(:,:)

! integer, save :: icount=0

  ! For deterministic heat-bath routines. Combine target orbitals {r,s} with the absolute value of the matrix element that takes you there, absH


  type(rs_absH),pointer,dimension(:),contiguous :: dtm_hb
  integer,allocatable :: pq_ind(:) ! index of start of pq terms in dtm_hb
  integer,allocatable :: pq_count(:) ! count of start of pq terms in dtm_hb
  real(rk) :: max_double

  integer,allocatable :: orbs_by_sym(:,:),norbs_by_sym(:)

contains

  !===========================================================================
  subroutine read_chem
    !---------------------------------------------------------------------------
    ! Description : Read input relevant to chemistry.
    !---------------------------------------------------------------------------
!*** Edited by AR [7/24/12]: changed I/O to handle parallelism
    use mpi_routines, only:master_core,mpi_bsend,mpi_bsend_diffroot,mpi_barr,mpi_stop
    use common_selected_ci, only : g,u,lz,up,dn
!***
    implicit none
    integer :: i

    tmp_ctr=0
    counter_two_body = -10  ! set negative here so that we use regular two_body routine initially


    if(master_core) read(5,*) nelec, nup
    call mpi_bsend(nelec)
    call mpi_bsend(nup)

    ndn=nelec-nup
    allocate(hf_up_orbs(nup))
    allocate(hf_dn_orbs(ndn))
    hf_up_orbs(:) = 0
    hf_dn_orbs(:) = 0
    if (up(1).ne.0) then
      do i=1,nup
        hf_up_orbs(i) = up(i)
      enddo
      do i=1,ndn
        hf_dn_orbs(i) = dn(i)
      enddo
    endif

    write(6,'(''nelec, nup, ndn='',9i5)') nelec, nup, ndn
!    write(6,'(''nelec, nup, ndn='',9i5)') nelec, nup, ndn


    if (nup < ndn) stop 'nup < ndn'

    if(master_core) then
        read(5,*) point_group
        write(6,'(''point_group='',a4)') point_group
        read(5,*) time_sym
    endif
    call mpi_bsend(point_group)
    call mpi_bsend(time_sym)
    if (point_group.ne.'dih') then
      if (g.or.u.or.lz.ne.0) then
        write (6,*) "point_group=",point_group,"lz=",lz,"g=",g,"u=",u
        write (6,*) "Lz and g/u specification only works with Lz integrals! Specify target irrep explicitly using Molpro convention"; call flush(6)
        call mpi_barr()
        call mpi_stop("Lz and g/u specification only works with Lz integrals! Specify target irrep explicitly using Molpro convention")
      endif
    endif


    if (time_sym) then

       if(master_core) then
         write(6,*) 'Time reversal symmetry of chem being used'
         read(5,*) z            ! Time reversal and reflection (parity)
         write(6,'(''z(time reversal)'', 9i4)') z
       endif
       call mpi_bsend(z)

       if (nup .ne. ndn) then
         write (6,*) "Time symmetry for chemistry will not work when nup and ndn are not equal! Stopping the program...."
         stop 'Time symmetry for chemistry will not work when nup and ndn are not equal! Stopping the program....'
       endif
       if(abs(z).ne.1) then
         write(6,'(''z(time reversal) must be either 1 or -1'')')
         stop 'z(time reversal) must be either 1 or -1'
       endif
    else
       write(6,*) 'Time reversal of chemistry not being used'
    endif

    if(master_core) then
        read(5,*) norb
     endif
     call mpi_bsend(norb)


    if(run_type.ne.'hci') then
      if(master_core) then
          read(5,*) n_core_orb
          read(5,*) trial_wf_iters
       endif
       call mpi_bsend(n_core_orb)
       call mpi_bsend(trial_wf_iters)

       allocate(norb_trial_wf(max(0,trial_wf_iters)))
       allocate(n_initiators_trial_wf(max(0,trial_wf_iters)))
       allocate(n_truncate_trial_wf(max(0,trial_wf_iters)))
!
       if(master_core) then
           read(5,*) norb_trial_wf
           read(5,*) n_initiators_trial_wf
           read(5,*) n_truncate_trial_wf
           write (6,*) "trial_wf_iters",trial_wf_iters
           write (6,*) "norb_trial_wf",norb_trial_wf
           write (6,*) "n_initiators_trial_wf",n_initiators_trial_wf
           write (6,*) "n_truncate_trial_wf",n_truncate_trial_wf
           write(6,'(''n_core_orb='',i5)') n_core_orb
       endif
       call mpi_bsend_diffroot(norb_trial_wf,0)
       call mpi_bsend_diffroot(n_initiators_trial_wf,0)
       call mpi_bsend_diffroot(n_truncate_trial_wf,0)

       if(maxval(norb_trial_wf).gt.norb) then
          write(6,'(''maxval(norb_trial_wf) must be <= norb; norb_trial_wf, norb='',9i5)') maxval(norb_trial_wf), norb
         stop 'maxval(norb_imp) must be <= norb'
       endif
    endif ! run_type.ne.'hci'

    allocate(orbital_symmetries(norb))
    allocate(orbital_energies(norb))

    if(master_core) then
        read(5,*) orbital_symmetries
        write(6,'(''orbital_symmetries='',1000i3)') orbital_symmetries
        read(5,*) spatial_symmetry_wf
        write(6,'(''spatial_symmetry_wf='',1000i3)') spatial_symmetry_wf
        call flush(6)
        if(run_type.ne.'hci') then
          read(5,*) diagonalize_ham
          write(6,'(''diagonalize_ham (chemistry) ='',i4)') diagonalize_ham
        endif
    endif
    call mpi_bsend_diffroot(orbital_symmetries,0)
    call init_point_group()
    call mpi_bsend(spatial_symmetry_wf)
    call mpi_bsend(diagonalize_ham)

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
      call mpi_bsend_diffroot(norb_det_sel,0)
      call mpi_bsend_diffroot(n_sym_uniq_det_det_sel,0)
      call mpi_bsend(init_multiplier)
      call mpi_bsend(trunc_multiplier)

      norb_det_sel(:) = init_multiplier*norb_det_sel(:)
      n_sym_uniq_det_det_sel(:) = trunc_multiplier*n_sym_uniq_det_det_sel(:)

      write (6,*) "det_sel_iters",det_sel_iters
      write (6,*) "norb_det_sel",norb_det_sel
      write (6,*) "n_sym_uniq_det_det_sel",n_sym_uniq_det_det_sel
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
    call flush(6)
  end subroutine read_chem
  !===========================================================================

  subroutine system_setup_chem(hf_up,hf_dn)
    !---------------------------------------------------------------------------
    ! Description : Define system dependent parameters.
    !               ndn :: Number of down electrons
    !               nuclear_nuclear_energy :: constant nuclear-nuclear term in H
    !
    !               Read in one and two body integrals
    !               Initialize point group information
    !               Generate trial wavefunction
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    ! Modified    : A Holmes, 5 July 2016. Output the hf determinant
    !---------------------------------------------------------------------------
    use common_run, only: tau_multiplier,use_efficient_heatbath
    use common_ham, only: diagonal_ham_lowest, diagonal_ham_highest
    use common_run, only : connected_dets_up,connected_dets_dn,connected_diag_elems_info,n_connected_dets_hf
    use common_selected_ci, only : hf_det, get_auto_hf, hf_symmetry, up, dn
    use mpi_routines, only:mpi_bsend,mpi_stop,mpi_barr,mpi_barr_in_node,master_core,mpi_allred
    use more_tools, only : get_occ_orbs

    implicit none

    !local variables
    real(rk) :: one_body_energy, two_body_energy
!   real(rk) :: diagonal_ham_lowest, diagonal_ham_highest
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(out) :: hf_up,hf_dn
    type(ik_vec) :: max_det_up, max_det_dn
    type(ik_vec) :: current_det_up, current_det_dn
    type(ik_vec),allocatable :: all_dets_up(:),all_dets_dn(:)
#else
    integer(ik),intent(out) :: hf_up,hf_dn
    integer(ik) :: max_det_up, max_det_dn
    integer(ik) :: current_det_up, current_det_dn
    integer(ik),allocatable :: all_dets_up(:),all_dets_dn(:)
#endif
    integer :: n_excite_part, n_excite
    logical :: valid_combination
    real(rk),allocatable :: exact_lowest_eigenvector(:)
    real(rk) :: exact_lowest_eigenvalue,exact_tau
    integer  :: n_det_exact
    !type(ik_vec),allocatable    :: temp_i16_up(:),temp_i16_dn(:)!,temp_i_2(:)
    !integer,allocatable        :: iorder(:),temp_i_2(:)
    logical :: use_fcidump_instead ! whether to use FCIDUMP instead of integrals.dat to get the integrals
    integer :: i,j
    integer :: n_single
    integer :: n_double_up
    integer :: n_double_dn
    integer :: n_double_both
    integer :: n_double
    integer :: n_total

    call my_second(1,'system_setup_chem')

    ndn = nelec - nup
    allocate(occ_up(nup))
    allocate(occ_dn(ndn))

    n_pairs = (nelec*(nelec-1))/2
    allocate(tmp1(n_pairs))
    allocate(tmp2(n_pairs))
    sqrt2=sqrt(2._rk)
    sqrt2inv=1._rk/sqrt2

    ! allocate the connected dets vectors
    max_connected_dets = upper_bound_connections()
    allocate(connected_dets_up(max_connected_dets))
    allocate(connected_dets_dn(max_connected_dets))
    allocate(connected_matrix_elements(max_connected_dets))
    allocate(connected_diag_elems_info(max_connected_dets))

    ! also allocate these vectors, which are used to find connected dets
    allocate(filled_up(nup))
    allocate(filled_dn(ndn))
    allocate(empty_up(norb-nup))
    allocate(empty_dn(norb-ndn))

    allocate(orb_order(norb+1))
    allocate(orb_order_inv(norb+1))
    do i=1,norb+1
      orb_order(i) = i
      orb_order_inv(i) = i
    enddo

    ! Initialize combine_2(i,j)
    allocate(combine_2(norb+1,norb+1))
    do i=1,norb
      do j=1,norb
        if (i.gt.j) then
          combine_2(i,j)=(i*(i-1))/2+j
        else
          combine_2(i,j)=(j*(j-1))/2+i
        endif
      enddo ! j
    enddo ! i
    combine_2(norb+1,norb+1) = ((norb+1)*norb)/2+norb+1

    call read_integrals(hf_up,hf_dn,use_fcidump_instead)

    nuclear_nuclear_energy = integrals(integral_index(norb+1,norb+1,norb+1,norb+1))


    call setup_orb_by_symm()
    last_det_up = 0
    last_det_dn = 0
    last_energy = 0._rk


    !write pieces of HF energy
    call one_body(hf_up, hf_dn, one_body_energy)
    call two_body(hf_up, hf_dn, two_body_energy)
    call mp2_deterministic(hf_up, hf_dn, mp2_energy_correction)
    if (time_sym) then
      call hamiltonian_chem_time_sym(hf_up, hf_dn, hf_up, hf_dn, diagonal_ham_lowest)
    else
      call hamiltonian_chem(hf_up, hf_dn, hf_up, hf_dn, 0, diagonal_ham_lowest)
    endif
    write(6,*)
    write(6,'(''One-Body Energy        ='',f14.8)') one_body_energy
    write(6,'(''Two-Body Energy        ='',f14.8)') two_body_energy
    write(6,'(''Nuclear-Nuclear Energy ='',f14.8)') nuclear_nuclear_energy
    write(6,'(''HF Energy              ='',f14.8)') diagonal_ham_lowest !one_body_energy+two_body_energy+nuclear_nuclear_energy
    write(6,'(''MP2 Energy             ='',f14.8)') one_body_energy+two_body_energy+nuclear_nuclear_energy+mp2_energy_correction
    write(6,*)
    call my_second(2,'HF pieces')
    call flush(6)
    n_excite = nelec - 2 * n_core_orb
    n_excite_part = n_excite / 2  !want integer division here
    if (n_excite_part > ndn - n_core_orb) then
       n_excite_part = ndn - n_core_orb
    endif

    ! Estimate spectrum as the difference between highest and lowest diagonal elements
    ! Estimate the largest diagonal element by calculating the diagonal element of the determinant formed by occupying the highest-energy orbitals.
    ! Warning: may have wrong symmetry!

#ifdef NUM_ORBITALS_GT_127
    max_det_up = maskr_vec(norb)-maskr_vec(norb+n_core_orb-nup)+maskr_vec(n_core_orb)
    max_det_dn = maskr_vec(norb)-maskr_vec(norb+n_core_orb-ndn)+maskr_vec(n_core_orb)
#else
    max_det_up = maskr(norb,ik)-maskr(norb+n_core_orb-nup,ik)+maskr(n_core_orb,ik)
    max_det_dn = maskr(norb,ik)-maskr(norb+n_core_orb-ndn,ik)+maskr(n_core_orb,ik)
#endif

    call hamiltonian_chem(max_det_up, max_det_dn, max_det_up, max_det_dn, 0, diagonal_ham_highest)
    max_energy = diagonal_ham_highest

    if(run_type.eq.'hci') then
       write(6,'(/,''diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau='',2f13.6,9f10.6)') diagonal_ham_lowest, diagonal_ham_highest
    else
      write(6,'(/,''before setting tau='',d12.4)') tau
      if(tau.eq.0.d0) then
        tau=tau_multiplier/(diagonal_ham_highest-diagonal_ham_lowest)
         write(6,'(/,''diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau='',2f13.6,9f10.6)') diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau
       else
         write(6,'(/,''diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, deduced tau='',2f13.6,9f10.6)') diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau_multiplier/(diagonal_ham_highest-diagonal_ham_lowest)
         write(6,'(''Input value (actually used) of tau='',9f12.6)') tau
      endif
    endif

    if (diagonalize_ham .eq. 1) then
       write (6,*) "Performing diagonalization of full hamiltonian"
      n_det_exact = 0
      current_det_up = 2_ik ** nup - 1
      do while (current_det_up <= max_det_up)
        current_det_dn = 2_ik ** ndn - 1
        do while (current_det_dn <= max_det_dn)
           !check for correct symmetry and make sure that determinant is not too highly excited
           call filter_dets(current_det_up, current_det_dn, valid_combination)
           if (time_sym) then
             if (current_det_up .gt. current_det_dn) valid_combination=.false.
           endif
           if (valid_combination) then
              n_det_exact = n_det_exact + 1
              if(n_det_exact.gt.2000000) exit
           endif
           call twiddle_determinant(current_det_dn)
        enddo
        call twiddle_determinant(current_det_up)
      enddo
      write (6,*) "n_det_exact",n_det_exact
      if (n_det_exact .le. 2000000) then
          allocate (all_dets_up(n_det_exact))
          allocate (all_dets_dn(n_det_exact))
          allocate (exact_lowest_eigenvector(n_det_exact))
          n_det_exact = 0
          current_det_up = 2_ik ** nup - 1
          do while (current_det_up <= max_det_up)
            current_det_dn = 2_ik ** ndn - 1
            do while (current_det_dn <= max_det_dn)
               !check for correct symmetry and make sure that determinant is not too highly excited
               call filter_dets(current_det_up, current_det_dn, valid_combination)
               if (time_sym) then
                 if (current_det_up .gt. current_det_dn) valid_combination=.false.
               endif
               if (valid_combination) then
                  n_det_exact = n_det_exact + 1
                  all_dets_up(n_det_exact)=current_det_up
                  all_dets_dn(n_det_exact)=current_det_dn
               endif
               call twiddle_determinant(current_det_dn)
            enddo
            call twiddle_determinant(current_det_up)
          enddo
          call sort(all_dets_up,all_dets_dn)
          if (n_det_exact .le. 2000) then
              call diagonalize_hamiltonian(all_dets_up, all_dets_dn, exact_lowest_eigenvector, exact_lowest_eigenvalue)
          else
              call diagonalize_sparse_hamiltonian_chem(all_dets_up, all_dets_dn, n_det_exact, exact_lowest_eigenvector, exact_lowest_eigenvalue, exact_tau, time_sym)
          endif
          write(6,'(/,''n_det, exact energy ='',i8,f14.8)') n_det_exact, exact_lowest_eigenvalue
      else
           write (6,*) "Skipping diagonalization of full Hamiltonian because space size exceeded inbuilt threshold of 2*10^6"
      endif
    endif

    if (.not.use_efficient_heatbath) then

      n_single = (nup-n_core_orb)*(norb - nup) + (ndn-n_core_orb)*(norb - ndn)
      n_double_up = (nup-n_core_orb)*(nup -n_core_orb-1)*(norb-nup)*(norb-nup-1)/4
      n_double_dn = (ndn-n_core_orb)*(ndn -n_core_orb-1)*(norb-ndn)*(norb-ndn-1)/4
      n_double_both = (nup-n_core_orb)*(norb-nup)*(ndn-n_core_orb)*(norb-ndn)

      n_double = n_double_up + n_double_dn + n_double_both
      n_total = n_single + n_double

      proposal_prob_single=(n_single)/real(n_total,rk)
      proposal_prob_double=(n_double)/real(n_total,rk)        ! or 1 - proposal_prob_single

    endif

    call find_connected_dets_chem(hf_up, hf_dn, n_connected_dets_hf, connected_dets_up, connected_dets_dn, norb)

    call my_second(2,'system_setup_chem')

  end subroutine system_setup_chem
  !=============================================================

  !=============================================================================
  subroutine read_integrals(hf_up,hf_dn,use_fcidump_instead)
    !-------------------------------------------------------------------------
    ! Description : Reads and stores the 1-body and 2-body integrals from the file
    !               This can be done in a more memory efficient way.
    !
    ! Created     : Mihir Khadilkar,   12 Nov 2010 (modified by F. Petruzielo on Nov 22 2010)
    !-------------------------------------------------------------------------
    use types, only : num_words,i8b
    use mpi_routines, only:master_core,master_core_node,mpi_bsend,shmem_allocate,mpi_barr,mpi_stop,&
      mpi_barr_in_node,mpi_bsend_between_nodes
    use common_selected_ci, only : get_auto_hf, up, dn, hf_symmetry, irreps, n_irrep, irrep_occs_up, irrep_occs_dn
    use more_tools, only : get_occ_orbs

    !local variables
    logical,intent(out) :: use_fcidump_instead ! Returns false if integrals.dat file is present; else true
    integer:: stat, p, q, r, s
    real(rk) :: integral_value
    integer :: i
    integer(i8b) :: iSize, n
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(out) :: hf_up,hf_dn
    type(ik_vec) :: tmp
#else
    integer(ik),intent(out) :: hf_up,hf_dn
    integer(ik) :: tmp
#endif
    integer :: a,b,j
    logical :: singlet_triplet_error
    integer :: tmp_sym
    character(1600) :: string


! If integrals.dat exists, read integrals from it.  Otherwise read from FCIDUMP.
    inquire(file="integrals.dat", exist=use_fcidump_instead)
    use_fcidump_instead=.not.(use_fcidump_instead)

    if (use_fcidump_instead) then
      write (6,'(/,''Reading integrals from FCIDUMP file since integrals.dat does not exist, and then sorting'',/)')
    else
      write (6,'(/,''Reading integrals from integrals.dat file'',/)')
    endif
    call flush(6)

    if (norb > num_words*(bit_size(1_ik) - 1)) then
       write(6,'(''Integers used in the program cannot represent the requested number of orbitals'')')
       write(6,'(i5,a,i3)') norb, " > ", num_words*(bit_size(1_ik) - 1)
       stop 'Consider compiling with higher precision integers'
    endif

    iSize = integral_index(norb+1,norb+1,norb+1,norb+1)
    write (6,'(''About to allocate integrals to size'',i10)') iSize; call flush(6)
    !allocate(integrals(iSize))

    !Allocate shared memory using MPI
    call shmem_allocate(integrals,iSize)
    write (6,'(''Finished allocating integrals'',/)') ; call flush(6)

! Edited by AAH on 7 Jan 2014. read from FCIDUMP instead of integrals.dat file.
! Note: add lines listing the occupied orbitals for up and down electrons in HF before the "END" in the FCIDUMP header, e.g.:
    if (use_fcidump_instead) then
      call my_second(1,'Read from FCIDUMP')
      !each master_core_node reads the file
      !if(master_core) open(1, file= 'FCIDUMP', status='old')
      if (master_core) then
        open(1, file= 'FCIDUMP', status='old')
        do i=1,3
          read(1,'(A)') string
          write(6,*) trim(string) ; call flush(6)
        enddo
!       If HF occupancies have not been read in from the input, they can be read from FCIDUMP as before.
!       Putting HF det in input file instead is recommended, so that distinct FCIDUMP files
!       are not needed for runs for different symmetry sectors
        write (6,'(/,''Reading HF occupations from FCIDUMP file (if they exist)'')') ; call flush(6)
        do
          read(1,'(A)') string
          write(6,*) trim(string) ; call flush(6)
          if (up(1)==0) then
            if(index(string,'hf_up').ne.0) then
              read(string,*) hf_up_orbs(1:nup)
              write (6,'(''HF up orbitals occupied:'',999i4)') hf_up_orbs(1:nup)
            endif
            if(index(string,'hf_dn').ne.0) then
              read(string,*) hf_dn_orbs(1:ndn)
              write (6,'(''HF dn orbitals occupied:'',999i4)') hf_dn_orbs(1:ndn)
            endif
          endif
          if(index(string,'&').ne.0) exit
        enddo
      endif
      call mpi_bsend(hf_up_orbs)
      call mpi_bsend(hf_dn_orbs)


      ! If time-reversal symmetry and hf_det specified, make sure they're
      ! consistent with each other!
      if (time_sym.and.hf_symmetry==999) then
        if (z<0) then
          singlet_triplet_error = .true.
          do i=1,nup
            write (6,*) "hf_up_orbs=",hf_up_orbs
            write (6,*) "hf_dn_orbs=",hf_dn_orbs
            if (hf_up_orbs(i).ne.hf_dn_orbs(i)) then
              singlet_triplet_error = .false.
              exit
            endif
          enddo
          if (singlet_triplet_error) then
            close(1)
            write (6,*) ""
            write (6,*) "Error: Input asked for a triplet state starting from a singlet det"
            write (6,*) "z=",z
            write (6,*) "hf_up=",hf_up_orbs(1:nup)
            write (6,*) "hf_dn=",hf_dn_orbs(1:ndn)
            write (6,*) "Consider specifying symmetry sector in input file to have program choose starting det automatically, e.g.:"
            write (6,*) "'&hf_det hf_symmetry=6'"
            write (6,*) "Stopping program..."; call flush(6)
            call mpi_barr()
            call mpi_stop("Singlet/triplet error!")
          endif
        endif
      endif


      ! Read FCIDUMP
      if (master_core) then
        do
          !only master_core reads
          if(master_core) read(1,*,iostat=stat) integral_value,p,q,r,s
          if (p==0) p=norb+1
          if (q==0) q=norb+1
          if (r==0) r=norb+1
          if (s==0) s=norb+1

          if(stat < 0) then ! end of file, so exit infinite do loop
            close(1)
            exit
          elseif (stat > 0) then
            write(6,'(''Error reading FCIDUMP file'')')
            stop 'Error reading FCIDUMP file'
          endif

          !only master_core_node assigns integrals a value
          if (abs(integral_value)>1.e-9_rk) then ! Avoids zeroing out integrals for which some permutations of the indices are zero in D_inf_h
            integrals(integral_index(p,q,r,s)) = integral_value
          endif
        enddo

      endif !if master core
      !All master cores need to get the integrals file
      if (master_core_node) then
        call mpi_bsend_between_nodes(integrals)
      endif
      call mpi_barr_in_node()
      call my_second(2,'Read from FCIDUMP')

     !if (get_auto_hf) then
        if (up(1)==0.and.hf_up_orbs(1)==0) then
       !if (up(1)==0.and.hf_up_orbs(1)==0.and.n_irrep==0) then
          ! Assign first few orbitals to be occupied in HF det
          do i=1,nup
            hf_up_orbs(i) = i
          enddo
          do i=1,ndn
            hf_dn_orbs(i) = i
          enddo
        else
          ! Use the ones in the input file
          if (up(1)>0) then
            hf_up_orbs(1:nup) = up(1:nup)
            hf_dn_orbs(1:ndn) = dn(1:ndn)
          endif
        endif
     !endif
     !write (6,*) "nup,ndn=",nup,ndn; call flush(6)
     !write (6,*) "hf_up_orbs=",hf_up_orbs(1:nup)
     !write (6,*) "hf_dn_orbs=",hf_dn_orbs(1:ndn)

#ifdef NUM_ORBITALS_GT_127
      hf_up = 0
      hf_dn = 0
      do i=1,nup
        hf_up%v(1) = ibset(hf_up%v(1),hf_up_orbs(i)-1)
      enddo
      do i=1,ndn
        hf_dn%v(1) = ibset(hf_dn%v(1),hf_dn_orbs(i)-1)
      enddo
#else
      hf_up = 0_ik
      hf_dn = 0_ik
      do i=1,nup
        hf_up = ibset(hf_up,hf_up_orbs(i)-1)
      enddo
      do i=1,ndn
        hf_dn = ibset(hf_dn,hf_dn_orbs(i)-1)
      enddo
#endif
      if (master_core) then
        if (n_irrep>0) then
          ! Use the numbers of occupations per irrep
          write (6,*) ""
          write (6,*) "Assigning HCI starting det orbs by irrep occupancies:"
          write (6,*) "Irreps=",irreps(1:n_irrep)
          write (6,*) "Up electron occupancies=",irrep_occs_up(1:n_irrep)
          write (6,*) "Dn electron occupancies=",irrep_occs_dn(1:n_irrep)
          call assign_hf_occs_by_irrep(n_irrep,irreps,irrep_occs_up,irrep_occs_dn,orbital_symmetries,hf_up,hf_dn)
          write (6,*) "HF det from input irrep occupancies=",hf_up,hf_dn; call flush(6)
          call get_occ_orbs(hf_up,hf_up_orbs)
          call get_occ_orbs(hf_dn,hf_dn_orbs)
        endif
      endif
      call mpi_bsend(hf_up_orbs)
      call mpi_bsend(hf_dn_orbs)
      hf_up = 0_ik
      hf_dn = 0_ik
      do i=1,nup
        hf_up = ibset(hf_up,hf_up_orbs(i)-1)
      enddo
      do i=1,ndn
        hf_dn = ibset(hf_dn,hf_dn_orbs(i)-1)
      enddo
      write(6,*)
      write (6,*) "Occupied up orbitals=", hf_up_orbs(1:nup),"Symmetries=",orbital_symmetries(hf_up_orbs(1:nup))
      write (6,*) "Occupied dn orbitals=", hf_dn_orbs(1:ndn),"Symmetries=",orbital_symmetries(hf_dn_orbs(1:ndn))

      ! Reassign HF automatically if requested
      if (get_auto_hf.and.hf_symmetry.ne.999) then
        write (6,*) ""
        write (6,*) "Reassigning HF occupancies to make HF have target symmetry"; call flush(6)
        if (master_core) then
          call auto_assign_hci0_occs(nup, ndn, norb, hf_symmetry, hf_up, hf_dn)
          write (6,*) "Automatically chosen HF det=",hf_up,hf_dn; call flush(6)
          call get_occ_orbs(hf_up,hf_up_orbs)
          call get_occ_orbs(hf_dn,hf_dn_orbs)
        endif
        call mpi_bsend(hf_up_orbs)
        call mpi_bsend(hf_dn_orbs)
        write (6,*) "Occupied up orbitals=", hf_up_orbs(1:nup),"Symmetries=",orbital_symmetries(hf_up_orbs(1:nup))
        write (6,*) "Occupied dn orbitals=", hf_dn_orbs(1:ndn),"Symmetries=",orbital_symmetries(hf_dn_orbs(1:ndn))
        tmp_sym = 1
        do i=1,nup
          tmp_sym = product_table(tmp_sym, orbital_symmetries(hf_up_orbs(i)))
        enddo
        do i=1,ndn
          tmp_sym = product_table(tmp_sym, orbital_symmetries(hf_dn_orbs(i)))
        enddo
        write (6,*) "HF symmetry=",tmp_sym; call flush(6)
        hf_up = 0_ik
        hf_dn = 0_ik
        do i=1,nup
          hf_up = ibset(hf_up,hf_up_orbs(i)-1)
        enddo
        do i=1,ndn
          hf_dn = ibset(hf_dn,hf_dn_orbs(i)-1)
        enddo
      endif ! get_auto_hf

      write(6,*)
      write (6,*) "Before reordering, hf_up,hf_dn=",hf_up,hf_dn; call flush(6)
      call sort_integrals(hf_up,hf_dn,orb_order,orb_order_inv)

      ! If time-reversal symmetry, return the correct representative
      if (time_sym) then
        if (hf_dn<hf_up) then
          tmp = hf_up
          hf_up = hf_dn
          hf_dn = tmp
        endif
      endif

      write (6,*) "After reordering, hf_up,hf_dn=",hf_up,hf_dn
      call mpi_barr()

    else ! read from integrals.dat

      call my_second(1,'Read from integrals.dat')
      !only master_core_node reads
      if(master_core_node) open(1, file= 'integrals.dat', status='old')
      if (master_core_node) then
        n = 0
        do
           n = n+1
           read(1,*,iostat = stat) p, q, r, s, integral_value

           if(stat < 0) then
             !done reading
             close(1)
             exit
           elseif (stat > 0) then
               if(master_core) write(6,'(''Error reading integrals.dat file'')')
               stop 'Error reading integrals.dat file'
           endif

           !only master_core_node assigns value
           integrals(integral_index(p,q,r,s)) = integral_value
           if (mod(n,1000000)==0) then
             write (6,*) "Finished reading",n/1000000,"million integrals"; call flush(6)
           endif

        enddo
      endif ! master_core_node

      call mpi_barr()
      call my_second(2,'Read from integrals.dat')

#ifdef NUM_ORBITALS_GT_127
      hf_up = maskr_vec(nup)
      hf_dn = maskr_vec(ndn)
#else
      hf_up = maskr(nup,ik)
      hf_dn = maskr(ndn,ik)
#endif
      write (6,*) "7hf_up,hf_dn=",hf_up,hf_dn

      call compute_orbital_energies(hf_up,hf_dn,orbital_energies)

    endif ! integrals.dat

    ! Now, make combine_2(i,j)
    do i=1,norb
      a = orb_order(i)
      do j=1,norb
        b = orb_order(j)
        if (a.gt.b) then
          combine_2(i,j)=(a*(a-1))/2+b
        else
          combine_2(i,j)=(b*(b-1))/2+a
        endif
      enddo ! j
    enddo ! i
    call mpi_barr()

  end subroutine read_integrals
  !=============================================================

  subroutine setup_efficient_heatbath()
    ! Setup efficient heatbath probabilities
    ! A Holmes, Mar 2015 (moved to its own subroutine on 16 Oct 2015)

      use types, only : i8b
      use more_tools, only : setup_alias
      use mpi_routines, only:whoami,master_core,master_core_node,mpi_bsend,shmem_allocate,mpi_barr,mpi_stop,shmem_reallocate

      integer :: i,j,k,l
      real(rk) :: H_ijkl
      integer(i8b) :: size_same,size_opposite
      integer,allocatable :: orb_tmp1(:),orb_tmp2(:) ! For setup_alias routine
      logical :: is_heatbath_unbiased
#ifdef NUM_ORBITALS_GT_127
      type(ik_vec) :: zero_det
#else
      integer(ik) :: zero_det
#endif
      integer :: p,q,r,s
      integer(i8b) :: occ_pair_ind,unocc_pair_ind
      integer :: sym_p,sym_q,sym_r,s_ind
      integer(i8b) :: max_unocc,max_occ
      integer(i8b) :: dtm_hb_size

      if (run_type.ne.'hci' .and. proposal_method.ne.'fast_heatbath')  return

      call my_second(1,'setup_efficient_heatbath')

      if (run_type.eq.'hci') then

        allocate(norbs_by_sym(n_group_elements))
        norbs_by_sym(:) = 0
        allocate(orbs_by_sym(n_group_elements,norb))
        do i=1,norb
          norbs_by_sym(orbital_symmetries(i)) = norbs_by_sym(orbital_symmetries(i)) + 1
          orbs_by_sym(orbital_symmetries(i),norbs_by_sym(orbital_symmetries(i))) = i
        enddo

        dtm_hb_size=combine_2_indices(norb,2*norb)*combine_2_indices(norb,norb)
        write (6,*) "About to call shmem_allocate. size=",dtm_hb_size; call flush(6)
        call shmem_allocate(dtm_hb,dtm_hb_size) ! assumes p<q, since combine_2_indices is invariant under interchange of p and q
        allocate(pq_ind(int(combine_2_indices(norb,2*norb))))
        allocate(pq_count(int(combine_2_indices(norb,2*norb))))
        !MJO since we are using shared memory with dtm_hb, we only want master_core_node to do the setup

        if (master_core_node) then

          dtm_hb(:)%absH = 0._rk
          max_double = 0._rk
          unocc_pair_ind = 0

          ! Same spin
          do p=1,norb
            sym_p = orbital_symmetries(p)
            do q=p+1,norb
              sym_q = product_table(sym_p,orbital_symmetries(q))
              occ_pair_ind = int(combine_2_indices(p,q))
              pq_ind(occ_pair_ind) = unocc_pair_ind + 1
              pq_count(occ_pair_ind) = 0
              do r=1,norb
                sym_r = product_table(sym_q,get_inverse_dih(orbital_symmetries(r)))
                do s_ind=1,norbs_by_sym(sym_r)
                  s = orbs_by_sym(sym_r,s_ind)
                  if (s<r)  cycle
                  H_ijkl = abs(double_excitation_matrix_element_no_ref(p,q,r,s))
                  if (H_ijkl.ne.0._rk) then
                    unocc_pair_ind = unocc_pair_ind + 1
                    pq_count(occ_pair_ind) = pq_count(occ_pair_ind) + 1
                    dtm_hb(unocc_pair_ind)%r = r
                    dtm_hb(unocc_pair_ind)%s = s
                    dtm_hb(unocc_pair_ind)%absH = H_ijkl
                  endif
                enddo ! s
              enddo ! r
              ! Now, sort dtm_hb(occ_pair_ind,:) in decreasing order by absH
              if (pq_count(occ_pair_ind)>1)  call sort_rs_absH(dtm_hb(pq_ind(occ_pair_ind):pq_ind(occ_pair_ind)+pq_count(occ_pair_ind)-1))
              if (dtm_hb(pq_ind(occ_pair_ind))%absH>max_double)  max_double = dtm_hb(pq_ind(occ_pair_ind))%absH
            enddo ! q
          enddo ! p

          ! Opposite spin
          do p=1,norb
            sym_p = orbital_symmetries(p)
            do q=norb+p,2*norb
              sym_q = product_table(sym_p,orbital_symmetries(q-norb))
              occ_pair_ind = int(combine_2_indices(p,q))
              pq_ind(occ_pair_ind) = unocc_pair_ind + 1
              pq_count(occ_pair_ind) = 0
              do r=1,norb
                sym_r = product_table(sym_q,get_inverse_dih(orbital_symmetries(r)))
                do s_ind=1,norbs_by_sym(sym_r)
                  s = orbs_by_sym(sym_r,s_ind)+norb
                  H_ijkl = abs(double_excitation_matrix_element_no_ref(p,q,r,s))
                  if (H_ijkl.ne.0._rk) then
                    unocc_pair_ind = unocc_pair_ind + 1
                    pq_count(occ_pair_ind) = pq_count(occ_pair_ind) + 1
                    dtm_hb(unocc_pair_ind)%r = r
                    dtm_hb(unocc_pair_ind)%s = s
                    dtm_hb(unocc_pair_ind)%absH = H_ijkl
                  endif
                enddo ! s
              enddo ! r
              ! Now, sort dtm_hb(occ_pair_ind,:) in decreasing order by absH
              if (pq_count(occ_pair_ind)>1)  call sort_rs_absH(dtm_hb(pq_ind(occ_pair_ind):pq_ind(occ_pair_ind)+pq_count(occ_pair_ind)-1))
              if (dtm_hb(pq_ind(occ_pair_ind))%absH>max_double)  max_double = dtm_hb(pq_ind(occ_pair_ind))%absH
            enddo ! q
          enddo ! p

        endif ! master_core_node

        !MJO master_core (not master_core_node), true master, sends out his max_double value
        !    this could be mpi_bsend_in_node, but it is unimportant, since everyone just the needs
        !    the one equivalent value
        call mpi_bsend(max_double)
        call mpi_bsend(pq_ind)
        call mpi_bsend(pq_count)
        call mpi_bsend(unocc_pair_ind)
        write (6,*) "Max double excitation magnitude=",max_double
        write (6,*) "Number of entries in dtm_hb=",unocc_pair_ind; call flush(6)

        ! Finally, allocate dtm_hb to correct size
        call shmem_reallocate(dtm_hb,unocc_pair_ind)

        allocate(pairs_e1((nelec*(nelec-1))/2))
        allocate(pairs_e2((nelec*(nelec-1))/2))

      else ! run_type .ne. 'hci'

        ! Now, allocate all quantities needed for run
        allocate(orb_tmp1(2*norb))
        allocate(orb_tmp2(2*norb))

        ! Double excitation PDFs
        ! 2*norb because these are in terms of spin-orbitals. First norb are nup, last norb are ndn.
        call my_second(1,"Efficient heat-bath startup routines"); call flush(6)

        write (6,*) "Allocating partial sums for efficient heat-bath sampling"; call flush(6)
        allocate(one_orbital_probabilities(norb))
        allocate(two_orbital_probabilities(2*norb,2*norb))

        write (6,*) "Finished allocating 2-index partial sums for efficient heat-bath sampling"; call flush(6)

        allocate(three_orbital_probabilities_same_spin(norb,norb,norb))
        allocate(three_orbital_probabilities_opposite_spin(norb,norb,norb))

        allocate(J_three_orbital_probabilities_same_spin(norb,norb,norb))
        allocate(q_three_orbital_probabilities_same_spin(norb,norb,norb))
        allocate(J_three_orbital_probabilities_opposite_spin(norb,norb,norb))
        allocate(q_three_orbital_probabilities_opposite_spin(norb,norb,norb))

        allocate(Htot_same(combine_2_indices(norb,norb),norb))
        allocate(Htot_opposite(norb,norb,norb))

        write (6,*) "Finished allocating 3-index partial sums for efficient heat-bath sampling"; call flush(6)

        size_same = same_index(norb,norb,norb,norb)
        call shmem_allocate(four_orbital_probabilities_same_spin,size_same)
        call shmem_allocate(J_four_orbital_probabilities_same_spin,size_same)
        call shmem_allocate(q_four_orbital_probabilities_same_spin,size_same)

        size_opposite = opposite_index(norb,norb,norb,norb)
        call shmem_allocate(four_orbital_probabilities_opposite_spin,size_opposite)
        call shmem_allocate(J_four_orbital_probabilities_opposite_spin,size_opposite)
        call shmem_allocate(q_four_orbital_probabilities_opposite_spin,size_opposite)

        write (6,*) "Finished allocating 4-index partial sums for efficient heat-bath sampling"; call flush(6)

        one_orbital_probabilities(:) = 0._rk

        two_orbital_probabilities(:,:) = 0._rk

        three_orbital_probabilities_same_spin(:,:,:) = 0._rk
        three_orbital_probabilities_opposite_spin(:,:,:) = 0._rk

        if (master_core_node) then
          four_orbital_probabilities_same_spin(:) = 0._rk
          four_orbital_probabilities_opposite_spin(:) = 0._rk
        endif

        Htot_same(:,:) = 0._rk
        Htot_opposite(:,:,:) = 0._rk

        write (6,*) "Generating efficient heatbath probabilities tensors"; call flush(6)
        zero_det = 0_ik

        ! Same spin
        do i=n_core_orb+1,norb
          if (mod(i,2)==0)  write (6,*) nint(100*float(i-1)/float(norb)),"% done"; call flush(6)
          do j=n_core_orb+1,norb
            do k=n_core_orb+1,norb
              do l=n_core_orb+1,norb

                if (i==j .or. k==l .or. i==k .or. i==l .or. j==k .or. j==l)  cycle
                ! No time-reversal symmetry here because there are no determinants! (only orbitals ijkl)
                call hamiltonian_chem(ibset(ibset(zero_det,(i-1)),(j-1)),zero_det,ibset(ibset(zero_det,(k-1)),(l-1)),zero_det,2,H_ijkl,nosign=.true.)

                one_orbital_probabilities(i) = one_orbital_probabilities(i) + abs(H_ijkl)

                two_orbital_probabilities(i,j) = two_orbital_probabilities(i,j) + abs(H_ijkl)
                two_orbital_probabilities(norb+i,norb+j) = two_orbital_probabilities(norb+i,norb+j) + abs(H_ijkl)

                three_orbital_probabilities_same_spin(i,j,k) = three_orbital_probabilities_same_spin(i,j,k) + abs(H_ijkl)

                if (master_core_node)  four_orbital_probabilities_same_spin(same_index(i,j,k,l)) = real(four_orbital_probabilities_same_spin(same_index(i,j,k,l)) + abs(H_ijkl))

              enddo
            enddo
          enddo
        enddo

        ! Opposite spin
        do i=n_core_orb+1,norb
          do j=n_core_orb+1,norb
            do k=n_core_orb+1,norb
              do l=n_core_orb+1,norb
                ! i,k have same spin, j,l have same spin
                if (.not.(i==k.or.j==l)) then

                  ! No time-reversal symmetry here because there are no determinants! (only orbitals ijkl)
                  call hamiltonian_chem(ibset(zero_det,i-1),ibset(zero_det,j-1),ibset(zero_det,k-1),ibset(zero_det,l-1),2,H_ijkl,nosign=.true.)
                  one_orbital_probabilities(i) = one_orbital_probabilities(i) + abs(H_ijkl)

                  two_orbital_probabilities(i,norb+j) = two_orbital_probabilities(i,norb+j) + abs(H_ijkl)
                  two_orbital_probabilities(norb+i,j) = two_orbital_probabilities(norb+i,j) + abs(H_ijkl)

                  three_orbital_probabilities_opposite_spin(i,j,k) = three_orbital_probabilities_opposite_spin(i,j,k) + abs(H_ijkl)
                  if (master_core_node)  four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,l)) = real(four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,l)) + abs(H_ijkl))
                endif
                ! i,l have same spin, j,k have same spin
                if (.not.(i==l.or.j==k)) then
                  ! No time-reversal symmetry here because there are no determinants! (only orbitals ijkl)
                  call hamiltonian_chem(ibset(zero_det,i-1),ibset(zero_det,j-1),ibset(zero_det,l-1),ibset(zero_det,k-1),2,H_ijkl,nosign=.true.)
                  one_orbital_probabilities(i) = one_orbital_probabilities(i) + abs(H_ijkl)
                  two_orbital_probabilities(i,norb+j) = two_orbital_probabilities(i,norb+j) + abs(H_ijkl)
                  two_orbital_probabilities(norb+i,j) = two_orbital_probabilities(norb+i,j) + abs(H_ijkl)
                endif
              enddo
            enddo
          enddo
        enddo

        call mpi_barr()
        write (6,*) "Done generating efficient heatbath probabilites tensors"; call flush(6)

        ! Compute cumulative three orbital probabilities
        do i=n_core_orb+1,norb
          do j=n_core_orb+1,norb
            do k=n_core_orb+1,norb
              do l=n_core_orb+1,norb
                if (i>j) then
                  Htot_same(combine_2_indices(i,j),k) = Htot_same(combine_2_indices(i,j),k) + four_orbital_probabilities_same_spin(same_index(i,j,k,l)) * 2._rk
                endif
                Htot_opposite(i,j,k) = Htot_opposite(i,j,k) + four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,l)) * 2._rk
              enddo
              ! Same spin
              if (k==n_core_orb+1) then
                c_three_orbital_probabilities_same_spin = three_orbital_probabilities_same_spin(i,j,k)
              else
                c_three_orbital_probabilities_same_spin = three_orbital_probabilities_same_spin(i,j,k) + c_three_orbital_probabilities_same_spin
              endif
              ! Opposite spin
              if (k==n_core_orb+1) then
                c_three_orbital_probabilities_opposite_spin = three_orbital_probabilities_opposite_spin(i,j,k)
              else
                c_three_orbital_probabilities_opposite_spin = three_orbital_probabilities_opposite_spin(i,j,k) + c_three_orbital_probabilities_opposite_spin
              endif
            enddo ! k
            if (i<=norb .and. j<=norb) then
              ! Same spin
              if (c_three_orbital_probabilities_same_spin > 0._rk) then
                three_orbital_probabilities_same_spin(i,j,:) = three_orbital_probabilities_same_spin(i,j,:)/c_three_orbital_probabilities_same_spin
              endif
              ! Opposite spin
              if (c_three_orbital_probabilities_opposite_spin > 0._rk) then
                three_orbital_probabilities_opposite_spin(i,j,:) = three_orbital_probabilities_opposite_spin(i,j,:)/c_three_orbital_probabilities_opposite_spin
              endif
            endif
          enddo
        enddo
        ! Now setup three,four_orb_probabilities for fast sampling
        do i=n_core_orb+1,norb
          do j=n_core_orb+1,norb
            call setup_alias(norb,three_orbital_probabilities_same_spin(i,j,:),orb_tmp1,orb_tmp2,J_three_orbital_probabilities_same_spin(i,j,:),q_three_orbital_probabilities_same_spin(i,j,:))
            call setup_alias(norb,three_orbital_probabilities_opposite_spin(i,j,:),orb_tmp1,orb_tmp2,J_three_orbital_probabilities_opposite_spin(i,j,:),q_three_orbital_probabilities_opposite_spin(i,j,:))
          enddo
        enddo

        write (6,*) "Normalizing 4-index tensors"; call flush(6)

        ! 4 index tensors
        if (master_core_node) then
          do i=n_core_orb+1,norb
           !write (6,*) "i=",i; call flush(6)
            do j=n_core_orb+1,norb
              do k=n_core_orb+1,norb
                do l=n_core_orb+1,norb
                  if (l==n_core_orb+1) then
                    c_four_orbital_probabilities_same_spin = four_orbital_probabilities_same_spin(same_index(i,j,k,l))
                  else
                    c_four_orbital_probabilities_same_spin = c_four_orbital_probabilities_same_spin + four_orbital_probabilities_same_spin(same_index(i,j,k,l))
                  endif
                enddo
                if (c_four_orbital_probabilities_same_spin > 0._rk) then
                  four_orbital_probabilities_same_spin(same_index(i,j,k,1):same_index(i,j,k,norb)) = real(four_orbital_probabilities_same_spin(same_index(i,j,k,1):same_index(i,j,k,norb))/c_four_orbital_probabilities_same_spin)
                  call setup_alias(norb,four_orbital_probabilities_same_spin(same_index(i,j,k,1):same_index(i,j,k,norb)),orb_tmp1,orb_tmp2,J_four_orbital_probabilities_same_spin(same_index(i,j,k,1):same_index(i,j,k,norb)),q_four_orbital_probabilities_same_spin(same_index(i,j,k,1):same_index(i,j,k,norb)))
                endif
                do l=n_core_orb+1,norb
                  if (l==n_core_orb+1) then
                    c_four_orbital_probabilities_opposite_spin = four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,l))
                  else
                    c_four_orbital_probabilities_opposite_spin = c_four_orbital_probabilities_opposite_spin + four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,l))
                  endif
                enddo
                if (c_four_orbital_probabilities_opposite_spin > 0._rk) then
                  four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,1):opposite_index(i,j,k,norb)) = real(four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,1):opposite_index(i,j,k,norb))/c_four_orbital_probabilities_opposite_spin)
                  call setup_alias(norb,four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,1):opposite_index(i,j,k,norb)),orb_tmp1,orb_tmp2,J_four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,1):opposite_index(i,j,k,norb)),q_four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,1):opposite_index(i,j,k,norb)))
                endif
              enddo! k
            enddo! j
          enddo ! i
        endif ! master_core_node

        call mpi_barr()
        write (6,*) "Done normalizing 4-index tensors"; call flush(6)

        call my_second(2,"Efficient heat-bath startup routines"); call flush(6)
        call mpi_barr()

        ! Allocate stuff

        allocate(elecs(nelec))
        allocate(e1_prob(nelec))
        allocate(e2_prob(nelec))
        allocate(c_e1_prob(nelec))
        allocate(c_e2_prob(nelec))

        allocate(pairs_e1((nelec*(nelec-1))/2))
        allocate(pairs_e2((nelec*(nelec-1))/2))
        allocate(pairs_prob((nelec*(nelec-1))/2))
        allocate(c_pairs_prob((nelec*(nelec-1))/2))
        allocate(J_pairs_prob((nelec*(nelec-1))/2))
        allocate(q_pairs_prob((nelec*(nelec-1))/2))

        ! Check to make sure that heatbath is unbiased for this system
        write (6,*) "Checking whether there is any potential bias for using heatbath on this system..." ; call flush(6)
        call check_heatbath_unbiased(is_heatbath_unbiased)
        if (is_heatbath_unbiased) then
          write (6,*) "Check passed! Heatbath is unbiased for this system!" ; call flush(6)
        else
          write (6,*) "Check failed! Heatbath is (potentially) biased for this system. Aborting run..."; call flush(6)
          if (run_type.ne.'hci')  call mpi_stop("Heatbath may be biased for this system!")
        endif

      endif ! run_type


      call my_second(2,'setup_efficient_heatbath')

  end subroutine setup_efficient_heatbath


  !=============================================================
  function integral_value(p, q, r, s)
    !---------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    ! Description : Returns 1-body or 2-body integral with indices pqrs.
    !               Note, if two of the indices are norb+1 then
    !               you are getting one-body integrals. If all of the indices in the
    !               array are norb+1 then you are getting the nuclear-nuclear energy
    !               This can be done in a more memory efficient way. (To be done)
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    !
    ! Orbital reordering and compact integral indexing scheme implemented by A Holmes
    ! D_infinity_h conservation of L_z check implemented by A Holmes on 16 May 2014
    !---------------------------------------------------------------------------
    implicit none

    !dummy arguments
    integer, intent(in) :: p,q,r,s
    real(rk) :: integral_value

    integral_value = integrals(integral_index(p,q,r,s))

  end function integral_value
  !=============================================================

  !===========================================================================
  subroutine hamiltonian_chem(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level, matrix_element,nosign)
    !---------------------------------------------------------------------------
    ! Description : Calculate matrix element of hamiltonian for det_i and det_j (split into up and down spin)
    !               Note det_i_up, det_i_dn, det_j_up, det_j_dn are integers with their binary
    !               representation giving the occupancies of the spatial orbitals.
    !               excite_level specifies the level of excitation
    !               ***Assume**** that det_i and det_j are equal or related by
    !               either a single or double excitation. This is denoted by excite_level.
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(in) :: det_j_up
    type(ik_vec), intent(in) :: det_j_dn
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(in) :: det_j_up
    integer(ik), intent(in) :: det_j_dn
#endif
    integer, intent(in) :: excite_level
    real(rk), intent(out) :: matrix_element
    logical, optional, intent(in) :: nosign ! if no sign, then don't compute phase factor, so could return a right-sign or wrong-sign element of either sign

    !local
    real(rk) :: one_body_energy, two_body_energy

    if (excite_level /= 2 .and. present(nosign)) then
      write (6,*) "STOP: No-sign Hamiltonian elements only works for double excitations"; call flush(6)
      stop "No-sign Hamiltonian elements only works for double excitations"
    endif

    if (excite_level == 0) then
       !diagonal element of hamiltonian
       call one_body(det_i_up, det_i_dn, one_body_energy)
       call two_body(det_i_up, det_i_dn, two_body_energy)
       matrix_element = one_body_energy + two_body_energy + nuclear_nuclear_energy
    elseif (excite_level == 1) then
       !single excitation element of hamiltonian
       call one_body_single(det_i_up, det_i_dn, det_j_up, det_j_dn, one_body_energy)
       call two_body_single(det_i_up, det_i_dn, det_j_up, det_j_dn, two_body_energy)
       matrix_element = one_body_energy + two_body_energy
    elseif (excite_level == 2) then
       !double excitation element of hamiltonian
       if (present(nosign)) then
         call two_body_double(det_i_up, det_i_dn, det_j_up, det_j_dn, two_body_energy,nosign)
       else
         call two_body_double(det_i_up, det_i_dn, det_j_up, det_j_dn, two_body_energy)
       endif
       matrix_element = two_body_energy
   !else
   !   write(6,*) "Determinants ", det_i_up, det_i_dn, det_j_up, det_j_dn, " are related by an excitation of level ", excite_level
   !   stop "Determinants are not connected by the Hamiltonian!"
    endif

  end subroutine hamiltonian_chem
  !===========================================================================

  subroutine hamiltonian_chem_time_sym(det_i_up, det_i_dn, det_j_up, det_j_dn, matrix_element)
    !---------------------------------------------------------------------------
    ! Description : Time symmetrized matrix element
    ! Created     : H. J. Changlani, April 24 2012
    !---------------------------------------------------------------------------
    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up,det_i_dn,det_j_up,det_j_dn
#else
    integer(ik), intent(in) :: det_i_up,det_i_dn,det_j_up,det_j_dn
#endif
    real(rk), intent(out)   :: matrix_element

    real(rk)                :: matrix_element1,matrix_element2,norm_bra,norm_ketinv
    logical                 :: check
    integer                 :: excite_level

    matrix_element1=0._rk
    matrix_element2=0._rk
    norm_ketinv=1._rk
    norm_bra=1._rk
    check=.true.
    if (det_j_up.eq. det_j_dn) norm_ketinv=sqrt2inv
    if (det_i_up .eq. det_i_dn) then
        norm_bra=sqrt2
        check=.false.
    endif

    if ((det_i_up .eq. det_j_up) .and. (det_i_dn .eq. det_j_dn)) then
        excite_level=0
    else
        call excitation_level(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level)
    endif

    if (excite_level >= 0) then
          !det_i is connected to a det in the trial wavefunction
          call hamiltonian_chem(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level, matrix_element1)
    endif

    if (check) then
      if (det_j_up .ne. det_j_dn) then
        call excitation_level(det_i_dn, det_i_up, det_j_up, det_j_dn, excite_level)
        if (excite_level >= 0) then
          !det_i is connected to a det in the trial wavefunction
          call hamiltonian_chem(det_i_dn, det_i_up, det_j_up, det_j_dn, excite_level, matrix_element2)
        endif
      else
        matrix_element2=matrix_element1
      endif
    endif
    matrix_element=(norm_bra*norm_ketinv)*(matrix_element1+(z*matrix_element2))

  end subroutine hamiltonian_chem_time_sym
  !===========================================================================


  !===========================================================================
  subroutine one_body(det_i_up, det_i_dn, one_body_energy)
    !---------------------------------------------------------------------------
    ! Description : Calculate one_body energy for det_i (split into up and down spin)
    !               This includes both kinetic energy and electron-nuclear-potential
    !               Note det_i_up and det_i_dn are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    ! Modified    : A Holmes, 24 May 2012. Sped up.
    !---------------------------------------------------------------------------
    implicit none

    !dummy variables
    real(rk), intent(out) :: one_body_energy
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec) :: det
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik) :: det
#endif

    !local variables
    integer :: i
    real(rk) :: energy

    energy = 0._rk

    !handle up spins
    det = det_i_up
    do while (det .ne. 0)
      i = trailz(det)+1
      energy = energy + integral_value(i, i, norb + 1, norb + 1)
      det = ibclr(det,i-1)
    enddo

    !handle down spins
    if (det_i_dn == det_i_up) then
       !double one_body energy
       energy = energy * 2._rk
    elseif (det_i_dn /= 0) then
      det = det_i_dn
      do while (det .ne. 0)
        i = trailz(det)+1
        energy = energy + integral_value(i, i, norb + 1, norb + 1)
        det = ibclr(det,i-1)
      enddo
    endif

    one_body_energy = energy

  end subroutine one_body
  !===========================================================================

  !===========================================================================
  subroutine one_body_single(det_i_up, det_i_dn, det_j_up, det_j_dn, one_body_single_energy)
    !---------------------------------------------------------------------------
    ! Description : Calculate one_body term of hamiltonian for det_i and det_j
    !               where det_i and det_j are related by a single excitation.
    !               Note det_i_up and det_i_dn are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    ! Modified    : A Holmes, 11 Jun 2012. Sped up.
    !---------------------------------------------------------------------------
    use tools, only : permutation_factor

    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(in) :: det_j_up
    type(ik_vec), intent(in) :: det_j_dn
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(in) :: det_j_up
    integer(ik), intent(in) :: det_j_dn
#endif
    real(rk), intent(out) :: one_body_single_energy

    !local variables
    integer :: i_bit, j_bit

    if (det_i_up /= det_j_up) then
      i_bit = trailz(iand(det_i_up,not(det_j_up)))
      j_bit = trailz(iand(det_j_up,not(det_i_up)))
      one_body_single_energy = permutation_factor(det_i_up,det_j_up) * integral_value(i_bit + 1, j_bit + 1, norb + 1, norb + 1)
    else
      i_bit = trailz(iand(det_i_dn,not(det_j_dn)))
      j_bit = trailz(iand(det_j_dn,not(det_i_dn)))
      one_body_single_energy = permutation_factor(det_i_dn,det_j_dn) * integral_value(i_bit + 1, j_bit + 1, norb + 1, norb + 1)
    endif

  end subroutine one_body_single
  !===========================================================================

  !===========================================================================
  subroutine exchange_correction(com_up, com_dn, k, exchange_correction_energy)
  !===========================================================================
    !---------------------------------------------------------------------------
    ! Description : Calculate exchange correction to two_body energy for the
    !               orbital k . Used for improved two-body energy routine.
    !                     abs(k) is the spatial index of the orbital,while sign
    !               of k being positive (negative) denotes that it is a
    !               spin-up (down) orbital. Used in speeding up two_body
    !               routine. Here com_up and com_dn are common up and down
    !               parts of the current and last determinant.
    !
    ! Created     : Mihir Khadilkar, 8 Apr 2011
    ! Modified    : A Holmes, 24 May 2012. Sped up.
    !---------------------------------------------------------------------------
    implicit none
    !dummy variables
    integer, intent(in) :: k
    real(rk), intent(out) :: exchange_correction_energy
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: com_up, com_dn
    type(ik_vec) :: det
#else
    integer(ik), intent(in) :: com_up, com_dn
    integer(ik) :: det
#endif

    !local variables
    integer :: i,l
    real(rk) :: energy

    energy = 0._rk

    if(k.gt.0) then
      ! k is an up spin orbital
      det = com_up
      do while (det .ne. 0)
        i = trailz(det)+1
        energy = energy - integral_value(i,k,k,i)
        det = ibclr(det,i-1)
      enddo
    else
      ! k is a down spin orbital
      l = abs(k)
      det = com_dn
      do while (det .ne. 0)
        i = trailz(det)+1
        energy = energy - integral_value(i,l,l,i)
        det = ibclr(det,i-1)
      enddo
    endif

    exchange_correction_energy = energy

  end subroutine exchange_correction
  !===========================================================================

  !===========================================================================
  subroutine direct_correction(com_up, com_dn, k, direct_correction_energy)
  !===========================================================================
    !---------------------------------------------------------------------------
    ! Description : Calculate direct correction to two_body energy for the
    !               orbital k . abs(k) is the spatial index of the orbital
    !               while sign of k being positive (negative) denotes that it
    !               is a spin-up (down) orbital. Used in speeding up two_body
    !               routine. Here com_up and com_dn are common up and down
    !               parts of the current and last determinant.
    !
    ! Created     : Mihir Khadilkar, 8 Apr 2011
    ! Modified    : A Holmes, 24 May 2012. Sped up.
    !---------------------------------------------------------------------------
    implicit none
    !dummy variables
    integer, intent(in) :: k
    real(rk), intent(out) :: direct_correction_energy
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: com_up, com_dn
    type(ik_vec) :: det
#else
    integer(ik), intent(in) :: com_up, com_dn
    integer(ik) :: det
#endif

    !local variables
    integer :: i, l
    real(rk) :: energy

    energy = 0._rk

    if(k.gt.0) then
      ! k is an up spin orbital
      det = com_up
      do while (det .ne. 0)
        i = trailz(det)+1
        if (i.ne.k)  energy = energy + integral_value(i,i,k,k)
        det = ibclr(det,i-1)
      enddo
      det = com_dn
      do while (det .ne. 0)
        i = trailz(det)+1
        energy = energy + integral_value(i,i,k,k)
        det = ibclr(det,i-1)
      enddo
    else
      ! k is an dn spin orbital
      l = abs(k)
      det = com_up
      do while (det .ne. 0)
        i = trailz(det)+1
        energy = energy + integral_value(i,i,l,l)
        det = ibclr(det,i-1)
      enddo
      det = com_dn
      do while (det .ne. 0)
        i = trailz(det)+1
        if (l.ne.i)  energy = energy + integral_value(i,i,l,l)
        det = ibclr(det,i-1)
      enddo
    endif

    direct_correction_energy = energy

  end subroutine direct_correction
 !===========================================================================

  !===========================================================================
  subroutine two_body(det_i_up, det_i_dn, two_body_energy)
    !---------------------------------------------------------------------------
    ! Description : Calculate two_body energy for det_i (split into up and down spin)
    !               This includes both direct and exchange terms.
    !               Note det_i_up and det_i_dn are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    ! Modified    : Mihir Khadilkar, 15 May 2011
    !               The routine now just calculates the difference between the two
    !                body energy of the last det it (routine) was called with
    !                and the current det. Uses exchange and direct correction
    !                routines to calculate the same.
    !---------------------------------------------------------------------------
    implicit none

    !dummy variables
    real(rk), intent(out) :: two_body_energy
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec) :: unique_det_i_up, unique_det_i_dn, unique_last_det_up, unique_last_det_dn, com_up, com_dn, temp_com_up, temp_com_dn
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik) :: unique_det_i_up, unique_det_i_dn, unique_last_det_up, unique_last_det_dn, com_up, com_dn, temp_com_up, temp_com_dn
#endif

    !local variables
    real(rk) :: direct_energy, exchange_energy, direct_correction_energy, exchange_correction_energy
    integer :: i, j, i_i, i_j, i_a, i_b, speedup_level
    logical :: use_speedup

    if (counter_two_body.ge.0) then
       counter_two_body = counter_two_body + 1
    endif
    if (counter_two_body.lt.0) then ! currently this just prevents modified version being used for the initial part, for accuracy purposes.
       use_speedup = .false.
    else
       use_speedup = .true.
    endif
    speedup_level = 0
    if ((last_det_up== det_i_up).and.(last_det_dn==det_i_dn)) then
       two_body_energy = last_energy
       use_speedup =.false.
       speedup_level = -1
    endif
    if(use_speedup) then
       unique_det_i_up    = iand(det_i_up, not(last_det_up))
       unique_det_i_dn    = iand(det_i_dn, not(last_det_dn))
       unique_last_det_up = iand(not(det_i_up), last_det_up)
       unique_last_det_dn = iand(not(det_i_dn), last_det_dn)
       com_up = iand(det_i_up, last_det_up)
       com_dn = iand(det_i_dn, last_det_dn)
       i_i = 0
       i_j = 0
       i_a = 0
       i_b = 0
       do i = 1, norb
          if(btest(unique_last_det_up,i-1)) then
             if(i_i==0) then
                i_i = i
                speedup_level = 1
             elseif(i_j==0) then
                i_j = i
                speedup_level = 2
             else
                speedup_level = 0
                exit  ! too excited to use speedup, need to use the normal metho
             endif
          endif
       enddo
       do i = 1, norb
          if(btest(unique_last_det_dn,i-1)) then
             if(i_i==0) then
                i_i = -i
                speedup_level = 1
             elseif(i_j==0) then
                i_j = -i
                speedup_level = 2
             else
                speedup_level = 0
                exit  ! too excited to use speedup, need to use the normal metho
             endif
          endif
       enddo
       do i = 1, norb
          if(btest(unique_det_i_up,i-1)) then
             if(i_a==0) then
                i_a = i
                speedup_level = 1
             elseif(i_b==0) then
                i_b = i
                speedup_level = 2
             else
                speedup_level =0
                exit  ! too excited to use speedup, need to use the normal metho
             endif
          endif
          if(btest(unique_det_i_dn,i-1)) then
             if(i_a==0) then
                i_a = -i
                speedup_level = 1
             elseif(i_b==0) then
                i_b = -i
                speedup_level = 2
             else
                speedup_level = 0
                exit  ! too excited to use speedup, need to use the normal metho
             endif
          endif
       enddo
       if (speedup_level==1) then
          two_body_energy = last_energy
          call direct_correction(com_up, com_dn, i_i, direct_correction_energy)
          two_body_energy = two_body_energy - direct_correction_energy
          call direct_correction(com_up, com_dn, i_a, direct_correction_energy)
          two_body_energy = two_body_energy + direct_correction_energy
          call exchange_correction(com_up, com_dn, i_i, exchange_correction_energy)
          two_body_energy = two_body_energy - exchange_correction_energy
          call exchange_correction(com_up, com_dn, i_a, exchange_correction_energy)
          two_body_energy = two_body_energy + exchange_correction_energy
          last_energy = two_body_energy
          last_det_up = det_i_up
          last_det_dn = det_i_dn

       elseif(speedup_level==2) then
          if ((i_i*i_j.lt.0).and.(i_i*i_a.lt.0)) then
             i_a = i_a + i_b
             i_b = i_a - i_b
             i_a = i_a - i_b
          endif
          two_body_energy = last_energy
          call direct_correction(com_up, com_dn, i_i, direct_correction_energy)
          two_body_energy = two_body_energy - direct_correction_energy
          call direct_correction(com_up, com_dn, i_a, direct_correction_energy)
          two_body_energy = two_body_energy + direct_correction_energy
          call exchange_correction(com_up, com_dn, i_i, exchange_correction_energy)
          two_body_energy = two_body_energy - exchange_correction_energy
          call exchange_correction(com_up, com_dn, i_a, exchange_correction_energy)
          two_body_energy = two_body_energy + exchange_correction_energy

          temp_com_up = com_up
          temp_com_dn = com_dn
          call direct_correction(temp_com_up, temp_com_dn, i_j, direct_correction_energy)
          two_body_energy = two_body_energy - direct_correction_energy
          call direct_correction(temp_com_up, temp_com_dn, i_b, direct_correction_energy)
          two_body_energy = two_body_energy + direct_correction_energy
          call exchange_correction(temp_com_up, temp_com_dn, i_j, exchange_correction_energy)
          two_body_energy = two_body_energy - exchange_correction_energy
          call exchange_correction(temp_com_up, temp_com_dn, i_b, exchange_correction_energy)
          two_body_energy = two_body_energy + exchange_correction_energy

          two_body_energy = two_body_energy + integral_value(abs(i_a),abs(i_a),abs(i_b),abs(i_b))
          two_body_energy = two_body_energy - integral_value(abs(i_i),abs(i_i),abs(i_j),abs(i_j))
          if ((i_a*i_b).gt.0) then
             two_body_energy = two_body_energy - integral_value(abs(i_b),abs(i_a),abs(i_a),abs(i_b))
          endif
          if ((i_i*i_j).gt.0) then
             two_body_energy = two_body_energy + integral_value(abs(i_j),abs(i_i),abs(i_i),abs(i_j))
          endif
          last_energy = two_body_energy
          last_det_up = det_i_up
          last_det_dn = det_i_dn
       endif
    endif
    if (((.not.use_speedup).or.speedup_level==0).and.(speedup_level.ne.-1)) then
       ! no speedup used. Do it the usual way
       exchange_energy = 0._rk
       direct_energy = 0._rk

       !exchange
       !handle up spins
       do i = 1, norb
          if (btest(det_i_up, i - 1)) then
             do j = i + 1, norb
                if (btest(det_i_up, j - 1)) then
                   exchange_energy = exchange_energy - integral_value(i, j, j, i)
                endif
             enddo
          endif
       enddo

       !handle down spins
       if (det_i_dn == det_i_up) then
          !double two_body energy
          exchange_energy = exchange_energy * 2._rk
       elseif (det_i_dn /= 0) then
          do i = 1, norb
             if (btest(det_i_dn, i - 1)) then
                do j = i + 1, norb
                   if (btest(det_i_dn, j - 1)) then
                      exchange_energy = exchange_energy - integral_value(i, j, j, i)
                   endif
                enddo
             endif
          enddo
          !else system is completely spin polarized
       endif

       !direct
       do i = 1, norb
          if (btest(det_i_up, i - 1)) then
             do j = i + 1, norb
                !up interacting with up
                if (btest(det_i_up, j - 1)) then
                   direct_energy = direct_energy + integral_value(i, i, j, j)
                endif
             enddo
             do j = 1, norb
                !up interacting with dn
                if (btest(det_i_dn, j - 1)) then
                   direct_energy = direct_energy + integral_value(i, i, j, j)
                endif
             enddo
          endif
          if (btest(det_i_dn, i - 1)) then
             do j = i + 1, norb
                !dn interacting with dn
                if (btest(det_i_dn, j - 1)) then
                   direct_energy = direct_energy + integral_value(i, i, j, j)
                endif
             enddo
          endif
       enddo

       two_body_energy = exchange_energy + direct_energy
       last_energy = two_body_energy
       last_det_up = det_i_up
       last_det_dn = det_i_dn
    endif
    use_speedup = .true.
  end subroutine two_body
  !===========================================================================

  !===========================================================================
  subroutine two_body_single(det_i_up, det_i_dn, det_j_up, det_j_dn, two_body_single_energy)
    !---------------------------------------------------------------------------
    ! Description : Calculate two_body term of hamiltonian for det_i and det_j
    !               where det_i and det_j are related by a single excitation.
    !               Note det_i_up and det_i_dn are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    ! Modified    : A Holmes, 12 Jun 2012. Sped up, but new order of operations
    !               gives different roundoff error (not statistically signficant).
    !---------------------------------------------------------------------------
    use tools, only : permutation_factor

    implicit none

    !dummy variables
    real(rk), intent(out) :: two_body_single_energy
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(in) :: det_j_up
    type(ik_vec), intent(in) :: det_j_dn
    type(ik_vec) :: det
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(in) :: det_j_up
    integer(ik), intent(in) :: det_j_dn
    integer(ik) :: det
#endif

    !local variables
    integer :: i, i_bit, j_bit
    real(rk) :: energy

    energy = 0._rk

    if (det_i_up /= det_j_up) then

      i_bit = trailz(iand(det_i_up,not(det_j_up)))+1
      j_bit = trailz(iand(det_j_up,not(det_i_up)))+1

      det = det_i_up
      do while (det .ne. 0)
        i = trailz(det)+1
        if (i /= i_bit .and. i /= j_bit) then
          energy = energy - integral_value(i_bit, i, i, j_bit) + integral_value(i_bit, j_bit, i, i) ! exchange + direct
        endif
        det = ibclr(det,i-1)
      enddo

      det = det_i_dn
      do while (det .ne. 0)
        i = trailz(det)+1
        energy = energy + integral_value(i_bit, j_bit, i, i) ! direct
        det = ibclr(det,i-1)
      enddo

      two_body_single_energy = permutation_factor(det_i_up,det_j_up) * energy

    else

      i_bit = trailz(iand(det_i_dn,not(det_j_dn)))+1
      j_bit = trailz(iand(det_j_dn,not(det_i_dn)))+1

      det = det_i_dn
      do while (det .ne. 0)
        i = trailz(det)+1
        if (i /= i_bit .and. i /= j_bit) then
          energy = energy - integral_value(i_bit, i, i, j_bit) + integral_value(i_bit, j_bit, i, i) ! exchange + direct
        endif
        det = ibclr(det,i-1)
      enddo

      det = det_i_up
      do while (det .ne. 0)
        i = trailz(det)+1
        energy = energy + integral_value(i_bit, j_bit, i, i) ! direct
        det = ibclr(det,i-1)
      enddo

      two_body_single_energy = permutation_factor(det_i_dn,det_j_dn) * energy

    endif

  end subroutine two_body_single
  !===========================================================================

  !===========================================================================
  subroutine two_body_double(det_i_up, det_i_dn, det_j_up, det_j_dn, two_body_double_energy,nosign)
    !---------------------------------------------------------------------------
    ! Description : Calculate two_body_double energy for det_i and det_j(split into up and down spin)
    !               Note det_i_up, det_i_dn, det_j_up, det_j_dn are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !               ***Assume**** that det_i and det_j are related by a double excitation
    !
    ! Created     : F. Petruzielo, 28 Oct 2010
    ! Modified    : A Holmes, 13 Jun 2012. Simplified and sped up.
    !---------------------------------------------------------------------------
    use tools, only : permutation_factor,permutation_factor2

    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(in) :: det_j_up
    type(ik_vec), intent(in) :: det_j_dn
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(in) :: det_j_up
    integer(ik), intent(in) :: det_j_dn
#endif
    real(rk), intent(out) :: two_body_double_energy
    logical, optional, intent(in) :: nosign

    !local variables
    integer :: first_i_bit, second_i_bit, first_j_bit, second_j_bit, gamma

    if (det_i_up==det_j_up) then

      if (present(nosign)) then
        call permutation_factor2(det_i_dn,det_j_dn,gamma,first_i_bit,second_i_bit,first_j_bit,second_j_bit,nosign)
        two_body_double_energy = (integral_value(first_i_bit + 1, first_j_bit + 1, second_i_bit + 1, second_j_bit + 1) - integral_value(first_i_bit + 1, second_j_bit + 1, second_i_bit + 1, first_j_bit + 1))
      else
        call permutation_factor2(det_i_dn,det_j_dn,gamma,first_i_bit,second_i_bit,first_j_bit,second_j_bit)
        two_body_double_energy = gamma * (integral_value(first_i_bit + 1, first_j_bit + 1, second_i_bit + 1, second_j_bit + 1) - integral_value(first_i_bit + 1, second_j_bit + 1, second_i_bit + 1, first_j_bit + 1))
      endif

    elseif (det_i_dn==det_j_dn) then

      if (present(nosign)) then
        call permutation_factor2(det_i_up,det_j_up,gamma,first_i_bit,second_i_bit,first_j_bit,second_j_bit,nosign)
        two_body_double_energy = (integral_value(first_i_bit + 1, first_j_bit + 1, second_i_bit + 1, second_j_bit + 1) - integral_value(first_i_bit + 1, second_j_bit + 1, second_i_bit + 1, first_j_bit + 1))
      else
        call permutation_factor2(det_i_up,det_j_up,gamma,first_i_bit,second_i_bit,first_j_bit,second_j_bit)
        two_body_double_energy = gamma * (integral_value(first_i_bit + 1, first_j_bit + 1, second_i_bit + 1, second_j_bit + 1) - integral_value(first_i_bit + 1, second_j_bit + 1, second_i_bit + 1, first_j_bit + 1))
      endif

    else

      first_i_bit = trailz(iand(det_i_up,not(det_j_up)))
      first_j_bit = trailz(iand(det_j_up,not(det_i_up)))
      second_i_bit = trailz(iand(det_i_dn,not(det_j_dn)))
      second_j_bit = trailz(iand(det_j_dn,not(det_i_dn)))

      if (present(nosign)) then
        two_body_double_energy = integral_value(first_i_bit + 1, first_j_bit + 1, second_i_bit + 1, second_j_bit + 1)
      else
        two_body_double_energy = permutation_factor(det_i_up, det_j_up) * permutation_factor(det_i_dn,det_j_dn) * integral_value(first_i_bit + 1, first_j_bit + 1, second_i_bit + 1, second_j_bit + 1)
      endif

    endif

  end subroutine two_body_double
  !=====================================================================================================

  !=====================================================================================================
  subroutine is_connected_chem(det_i_up,det_i_dn,det_j_up,det_j_dn,connected,excite_level,proposal_prob)
  !-------------------------------------------------------------------------------------------------------
  ! Created by : Hitesh Changlani
  ! Purpose    : Checks connectedness of two states. If they are possibly connected, check if they have right symmetry
  !              If they have the right symmetry then compute the proposal probability of getting det_j from det_i
  !              This proposal probability must be consistent with that in off_diagonal_move_chem
  ! Date       : April 30 2012
  !-------------------------------------------------------------------------------------------------------
  ! Updated by : Matt Otten
  ! Reason     : Sped up by utilizing popcnt and trailz instead of do loops. See revision 812 for previous version
  ! Date       : September 11 2015
  !-------------------------------------------------------------------------------------------------------

  use tools, only : count_bits_imp1

  implicit none
  logical,intent(out)    :: connected
  integer,intent(out)    :: excite_level
  real(rk),intent(out),optional :: proposal_prob
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: det_i_up,det_i_dn,det_j_up,det_j_dn
  type(ik_vec)            :: unique_i_up,unique_j_up,unique_i_dn,unique_j_dn
  type(ik_vec)            :: det,det1,det2
#else
  integer(ik),intent(in) :: det_i_up,det_i_dn,det_j_up,det_j_dn
  integer(ik)            :: unique_i_up,unique_j_up,unique_i_dn,unique_j_dn
  integer(ik)            :: det,det1,det2
#endif

  ! Local
  logical                :: check_up,check_dn
  integer                :: up_count1,up_count2,dn_count1,dn_count2
  integer                :: excite_from_1_i,excite_from_2_i,excite_to_1,excite_to_2
  integer                :: diff_up_1,diff_up_2,diff_up_3,diff_up_4,diff_dn_1,diff_dn_2,diff_dn_3,diff_dn_4
  integer                :: i,i_open,i_open2,sym1,diff_1,diff_2,diff_3,diff_4
  real(rk)               :: tmp_prob,tmp_prob2,denominator

  !write (6,*) "norb =",norb
  connected=.true.
  excite_level=-1
  if(present(proposal_prob)) proposal_prob=0._rk
  i_open=0
  i_open2=0
  tmp_prob=0._rk
  tmp_prob2=0._rk

  up_count1=0
  up_count2=0
  dn_count1=0
  dn_count2=0

  check_up=.true.
  check_dn=.true.

  if ((det_i_up .eq. det_j_up)) check_up=.false.                                                       ! Do I need to check up to compute excitation level?
  if ((det_i_dn .eq. det_j_dn)) check_dn=.false.                                                       ! Do I need to check downs to compute excitation level ?

  if (check_up) then
        unique_i_up = iand(det_i_up, not(det_j_up))                                            ! Count in how many places up determinants differ and in what location(s)
        unique_j_up = iand(det_j_up, not(det_i_up))
        !Count the differing places and return if more than 2 are different.
        up_count1 = popcnt(unique_i_up)
        up_count2 = popcnt(unique_j_up)
        if((up_count1>2).or.(up_count1.ne.up_count2).or.(up_count2>2)) then
           connected = .false.
           excite_level = -1
           return
        endif
        !find where the set bits are
        diff_up_1 = trailz(unique_i_up)
        unique_i_up = ibclr(unique_i_up,diff_up_1)
        diff_up_3 = trailz(unique_j_up)
        unique_j_up = ibclr(unique_j_up,diff_up_3)

        !Counting from 1, not from 0
        diff_up_1 = diff_up_1 + 1
        diff_up_3 = diff_up_3 + 1

        if(up_count1.eq.2) then
           diff_up_2 = trailz(unique_i_up)
           diff_up_4 = trailz(unique_j_up)
           !Counting from 1, not from 0
           diff_up_2 = diff_up_2 + 1
           diff_up_4 = diff_up_4 + 1
        endif
  endif

  if (check_dn) then
        unique_i_dn = iand(det_i_dn, not(det_j_dn))                                            ! Count in how many places dn determinants differ and in what location(s)
        unique_j_dn = iand(det_j_dn, not(det_i_dn))
        !Count the differing places and return if more than 2 are different.
        dn_count1 = popcnt(unique_i_dn)
        dn_count2 = popcnt(unique_j_dn)
        if((dn_count1>2).or.(dn_count1.ne.dn_count2).or.(dn_count2>2)) then
           connected = .false.
           excite_level = -1
           return
        endif
        !find Where the set bits are
        diff_dn_1 = trailz(unique_i_dn)
        unique_i_dn = ibclr(unique_i_dn,diff_dn_1)
        diff_dn_3 = trailz(unique_j_dn)
        unique_j_dn = ibclr(unique_j_dn,diff_dn_3)

        !Counting from 1, not from 0
        diff_dn_1 = diff_dn_1 + 1
        diff_dn_3 = diff_dn_3 + 1

        if(dn_count1.eq.2) then
           diff_dn_2 = trailz(unique_i_dn)
           diff_dn_4 = trailz(unique_j_dn)
           !Counting from 1, not from 0
           diff_dn_2 = diff_dn_2 + 1
           diff_dn_4 = diff_dn_4 + 1
        endif
  endif

  excite_level=up_count1+dn_count1
  if (excite_level .gt. 2) then
     excite_level=-1
     connected=.false.
     return
  endif

  if (present(proposal_prob)) then
     ! Proposal probability has to be consistent with off_diagonal_move_chem
     ! The multiplication by (n_single/n_total) done in off_diagonal_move_chem
     if (excite_level .eq. 1) then
        det=det_i_dn
        diff_1=diff_dn_1
        diff_2=diff_dn_3
        if (up_count1 .eq. 1) then
           det=det_i_up
           diff_1=diff_up_1
           diff_2=diff_up_3
        endif
        sym1=orbital_symmetries(diff_1)
        if (sym1.ne. orbital_symmetries(diff_2)) then
           connected=.false.
           return
        endif
        i_open=0
        do i=1,norb
           if (btest(det,i-1) .eqv. .false.) then
              if (orbital_symmetries(i) .eq. sym1) i_open=i_open+1                        ! The open orbitals must have the same spatial symmetry as the chosen electron
           endif
        enddo
        proposal_prob=1._rk/((nelec-2*n_core_orb)*(i_open))
        return
     endif

     ! Proposal probability has to be consistent with off_diagonal_move_chem
     ! The multiplication by (n_double/n_total) done in off_diagonal_move_chem

     if (excite_level .eq. 2) then
        if (up_count1 .eq. 2) then          ! Both electrons up
           !write (6,*) "A case of up up Exc 2"                               ! One electron up and one electron down
           diff_1=diff_up_1
           diff_2=diff_up_2
           diff_3=diff_up_3
           diff_4=diff_up_4
           det1=det_i_up
           det1=ibset(det1,diff_3-1)
           det2=det_i_up
           det2=ibset(det2,diff_4-1)
           tmp_prob=2._rk/(norb-nup)
           tmp_prob2=2._rk/(norb-nup)
        elseif (dn_count1 .eq. 2) then     ! Both electrons down
           !write (6,*) "A case of dn dn Exc 2"                               ! One electron up and one electron down
           diff_1=diff_dn_1
           diff_2=diff_dn_2
           diff_3=diff_dn_3
           diff_4=diff_dn_4
           det1=det_i_dn
           det1=ibset(det1,diff_3-1)
           det2=det_i_dn
           det2=ibset(det2,diff_4-1)
           tmp_prob=2._rk/(norb-ndn)
           tmp_prob2=2._rk/(norb-ndn)
        else
           !write (6,*) "A case of up dn Exc 2"                               ! One electron up and one electron down
           diff_1=diff_up_1
           diff_2=diff_dn_1
           diff_3=diff_up_3
           diff_4=diff_dn_3
           det2=det_i_up
           det1=det_i_dn
           !det1=ibset(det1,diff_3-1)
           !det2=ibset(det2,diff_4-1)
           tmp_prob=1._rk/(norb-nup)
           tmp_prob2=1._rk/(norb-ndn)
        endif
        sym1=product_table(orbital_symmetries(diff_1),orbital_symmetries(diff_2))
        if (sym1.ne. product_table(orbital_symmetries(diff_3),orbital_symmetries(diff_4))) then
           connected=.false.
           return
        endif
        !if(uniform_sampling) then
        if(proposal_method.ne.'CauchySchwarz') then
           do i=1,norb
              if (btest(det1,i-1) .eqv. .false.) then
                 if (product_table(orbital_symmetries(diff_3),orbital_symmetries(i))==sym1) i_open=i_open+1
              endif
           enddo
           do i=1,norb
              if (btest(det2,i-1) .eqv. .false.) then
                 if (product_table(orbital_symmetries(i),orbital_symmetries(diff_4))==sym1) i_open2=i_open2+1
              endif
           enddo
           !if ((i_open .eq. 0)  .and.   (i_open2 .ne. 0))   proposal_prob=(2._rk/((nelec-2*n_core_orb)*(nelec-2*n_core_orb-1)))*((tmp_prob2/(i_open2)))                               ! Consider both pathways of exciting
           !if ((i_open2 .eq. 0) .and.   (i_open  .ne. 0))   proposal_prob=(2._rk/((nelec-2*n_core_orb)*(nelec-2*n_core_orb-1)))*((tmp_prob/(i_open)))                                 ! Consider both pathways of exciting
           !if ((i_open2 .ne. 0) .and.   (i_open  .ne. 0))   proposal_prob=(2._rk/((nelec-2*n_core_orb)*(nelec-2*n_core_orb-1)))*((tmp_prob/(i_open))+(tmp_prob2/(i_open2)))           ! Consider both pathways of exciting
           if ((i_open .eq. 0)  .and.   (i_open2 .ne. 0))   proposal_prob=(1._rk/((nelec-2*n_core_orb)*(nelec-2*n_core_orb-1)))*((tmp_prob2/(i_open2)))                               ! Consider both pathways of exciting
           if ((i_open2 .eq. 0) .and.   (i_open  .ne. 0))   proposal_prob=(1._rk/((nelec-2*n_core_orb)*(nelec-2*n_core_orb-1)))*((tmp_prob/(i_open)))                                 ! Consider both pathways of exciting
           if ((i_open2 .ne. 0) .and.   (i_open  .ne. 0))   proposal_prob=(1._rk/((nelec-2*n_core_orb)*(nelec-2*n_core_orb-1)))*((tmp_prob/(i_open))+(tmp_prob2/(i_open2)))           ! Consider both pathways of exciting
           return
        else !CS sampling

           if (up_count1 .eq. 2) then          ! Both electrons up
              !Note:
              !diff_up_1 = occ_orb(excite_from_1_i) = i
              !diff_up_2 = occ_orb(excite_from_2_i) = j
              !diff_up_3 = excite_to_1 = k
              !diff_up_4 = excite_to_2 = l


              !See if there are open moves
              sym1 = orbital_symmetries(diff_up_4)
              i_open = num_orb_by_sym(sym1) - num_occ_up_by_sym(sym1)
              if(sym1.eq.orbital_symmetries(diff_up_3)) i_open = i_open-1
              if(i_open.eq.0) return !No moves possible

              !Translate diff_up_1 into an index, similar to excite_from_1_i
              do i = 1,nup
                 if(occ_orb_up(i).eq.diff_up_1) excite_from_1_i = i
                 if(occ_orb_up(i).eq.diff_up_2) excite_from_2_i = i
              enddo

              !Calculate P(i,j)
              proposal_prob = (cs_sqrt_prime_spin(excite_from_1_i)/sum_cs_sqrt_prime*&
                   &cs_sqrt_prime_spin(excite_from_2_i)/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_1_i)) +&
                   &cs_sqrt_prime_spin(excite_from_2_i)/sum_cs_sqrt_prime*&
                   &cs_sqrt_prime_spin(excite_from_1_i)/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_2_i)))

              excite_to_1 = diff_up_3
              excite_to_2 = diff_up_4

              !Calculate the denominator for P(j->l|i->k)
              denominator = 0
              do i = 1,num_occ_up_by_sym(sym1)
                 denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_up_by_sym(sym1,i))&
                      &+sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym1,i))
              enddo
              denominator = sym_sum_cs_sqrt(sym1,occ_orb(excite_from_2_i)) +&
                   &sym_sum_cs_sqrt(sym1,occ_orb(excite_from_1_i))- denominator
              if(sym1.eq.orbital_symmetries(diff_up_3)) denominator = denominator - &
                   &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1) -&
                   &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)

              !Calculate P(i->k)P(j->l|i->k)
              tmp_prob = (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)+&
                   &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1))&
                   &/(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))*&
                   &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)+&
                   &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2))&
                   &/denominator

              sym1 = orbital_symmetries(diff_up_3)
              i_open = num_orb_by_sym(sym1) - num_occ_up_by_sym(sym1)
              if(sym1.eq.orbital_symmetries(diff_up_4)) i_open = i_open-1
              if(i_open.ne.0) then
                 !Calculate P(i->l)P(j->k|i->l)
                 !If the the symmetries are the same, the denominator only changes a little

                 if(sym1.eq.orbital_symmetries(diff_up_4)) then
                    denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
                         &-sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
                         &+sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)&
                         &-sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)
                 else
                    denominator = 0
                    do i = 1,num_occ_up_by_sym(sym1)
                       denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_up_by_sym(sym1,i))&
                            &-sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym1,i))
                    enddo
                    denominator = sym_sum_cs_sqrt(sym1,occ_orb(excite_from_2_i)) +&
                         &sym_sum_cs_sqrt(sym1,occ_orb(excite_from_1_i)) - denominator
                 endif

                 tmp_prob = tmp_prob + (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)+&
                      &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2))&
                      &/(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))*&
                      &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
                      &+sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1))&
                      &/denominator

              endif
           elseif (dn_count1 .eq. 2) then     ! Both electrons down
              !Note:
              !diff_dn_1 = occ_orb(excite_from_1_i) = i
              !diff_dn_2 = occ_orb(excite_from_2_i) = j
              !diff_dn_3 = excite_to_1 = k
              !diff_dn_4 = excite_to_2 = l


              !See if there are open moves
              sym1 = orbital_symmetries(diff_dn_4)
              i_open = num_orb_by_sym(sym1) - num_occ_dn_by_sym(sym1)
              if(sym1.eq.orbital_symmetries(diff_dn_3)) i_open = i_open-1
              if(i_open.eq.0) return !No moves possible

              !Translate diff_dn_1 into an index, similar to excite_from_1_i
              do i = 1,ndn
                 if(occ_orb_dn(i).eq.diff_dn_1) excite_from_1_i = i + nup
                 if(occ_orb_dn(i).eq.diff_dn_2) excite_from_2_i = i + nup
              enddo

              proposal_prob = (cs_sqrt_prime_spin(excite_from_1_i)/sum_cs_sqrt_prime*&
                   &cs_sqrt_prime_spin(excite_from_2_i)/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_1_i)) +&
                   &cs_sqrt_prime_spin(excite_from_2_i)/sum_cs_sqrt_prime*&
                   &cs_sqrt_prime_spin(excite_from_1_i)/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_2_i)))

              excite_to_1 = diff_dn_3
              excite_to_2 = diff_dn_4

              !Calculate the denominator for P(j->l|i->k)
              denominator = 0
              do i = 1,num_occ_dn_by_sym(sym1)
                 denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_dn_by_sym(sym1,i))&
                      &+sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym(sym1,i))
              enddo

              denominator = sym_sum_cs_sqrt(sym1,occ_orb(excite_from_2_i)) &
                   &+sym_sum_cs_sqrt(sym1,occ_orb(excite_from_1_i)) - denominator

              if(sym1.eq.orbital_symmetries(diff_dn_3)) denominator = denominator - &
                   &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1) -&
                   &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)

              !Calculate P(i->k)P(j->l|i->k)
              tmp_prob = (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)+&
                   &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1))&
                   &/(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))*&
                   &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
                   &+sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2))&
                   &/denominator

              sym1 = orbital_symmetries(diff_dn_3)
              i_open = num_orb_by_sym(sym1) - num_occ_dn_by_sym(sym1)
              if(sym1.eq.orbital_symmetries(diff_dn_4)) i_open = i_open-1
              if(i_open.ne.0) then
                 !Calculate P(i->l)P(j->k|i->l)
                 !If the the symmetries are the same, the denominator only changes a little

                 if(sym1.eq.orbital_symmetries(diff_dn_4)) then
                    denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
                         &-sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
                         &+sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)&
                         &-sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)
                 else
                    denominator = 0 !Just a temporary variable
                    do i = 1,num_occ_dn_by_sym(sym1)
                       denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_dn_by_sym(sym1,i))&
                            &+sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym(sym1,i))
                    enddo
                    denominator = sym_sum_cs_sqrt(sym1,occ_orb(excite_from_2_i)) &
                         &+sym_sum_cs_sqrt(sym1,occ_orb(excite_from_1_i)) - denominator
                 endif

                 tmp_prob = tmp_prob + (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)&
                      &+sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2))&
                      &/(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))*&
                      &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
                      &+sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1))&
                      &/denominator

              endif
           else !Anti parallel excitation
              !Note:
              !diff_up_1 = occ_orb(excite_from_1_i) = i
              !diff_dn_1 = occ_orb(excite_from_2_i) = j
              !diff_up_3 = excite_to_1 = k
              !diff_dn_3 = excite_to_2 = l

              !See if there are open moves
              sym1 = orbital_symmetries(diff_dn_3)
              i_open = num_orb_by_sym(sym1) - num_occ_dn_by_sym(sym1)
              if(i_open.eq.0) return !No moves possible

              !Translate diff_dn_1 and diff_up_1into an index, similar to excite_from_1_i
              do i = 1,ndn
                 if(occ_orb_dn(i).eq.diff_dn_1) excite_from_2_i = i+nup
              enddo

              do i = 1,nup
                 if(occ_orb_up(i).eq.diff_up_1) excite_from_1_i = i
              enddo

              proposal_prob = (cs_sqrt_prime_spin(excite_from_1_i)/sum_cs_sqrt_prime*&
                   &cs_sqrt_prime_spin(excite_from_2_i)/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_1_i)) +&
                   &cs_sqrt_prime_spin(excite_from_2_i)/sum_cs_sqrt_prime*&
                   &cs_sqrt_prime_spin(excite_from_1_i)/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_2_i)))

              excite_to_1 = diff_up_3
              excite_to_2 = diff_dn_3

              !Calculate the denominator for P(j->l|i->k)
              denominator = 0
              do i = 1,num_occ_dn_by_sym(sym1)
                 denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_dn_by_sym(sym1,i))&
                      &+ sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym(sym1,i))
              enddo

              denominator = sym_sum_cs_sqrt(sym1,occ_orb(excite_from_2_i)) &
                   &+ sym_sum_cs_sqrt(sym1,occ_orb(excite_from_1_i))- denominator

              tmp_prob = (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)+&
                   &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1))&
                   &/(cs_sqrt_prime_spin(excite_from_1_i)+cs_sqrt_prime_spin(excite_from_2_i))*&
                   &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
                   &+sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2))&
                   &/denominator

              denominator = 0 !Just a temporary variable
              sym1 = orbital_symmetries(diff_up_3)
              i_open = num_orb_by_sym(sym1) - num_occ_up_by_sym(sym1)
              if(i_open.ne.0) then
                 do i = 1,num_occ_up_by_sym(sym1)
                    denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_up_by_sym(sym1,i))&
                         &+sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym1,i))
                 enddo

                 denominator = sym_sum_cs_sqrt(sym1,occ_orb(excite_from_2_i)) &
                      &+ sym_sum_cs_sqrt(sym1,occ_orb(excite_from_1_i)) - denominator

                 tmp_prob = tmp_prob + (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)&
                      &+sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2))&
                      &/(cs_sqrt_prime_spin(excite_from_1_i)+cs_sqrt_prime_spin(excite_from_2_i))*&
                      &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
                      &+sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1))&
                      &/denominator
              endif
           endif
           proposal_prob = proposal_prob*tmp_prob

        endif ! CauchySchwarz
     endif ! excite_level .eq. 2
  endif ! present(proposal_prob)
  !write (6,*) "Proposal prob", proposal_prob
  !call flush(6)

  return

  end subroutine is_connected_chem

!===================================================================================================
  subroutine setup_orb_by_symm()
    !----------------------------------------------------------------------------------------------------------
    ! Description : Sets up a table so that orbitals by a specific symmetry can be found
    !               swiftly.
    !
    ! Created     : Matt Otten, 19 Jan 2015
    ! last edited : Matt Otten, 19 Jan 2015
    !-----------------------------------------------------------------------------------------------------------
    implicit none
    integer i,sym1,j
    integer(i8b) :: int_index

    !Allocate all arrays
    allocate(num_orb_by_sym(n_group_elements))
    allocate(which_orb_by_sym(n_group_elements,norb))
    allocate(occ_orb(nup+ndn))
    allocate(occ_orb_up(nup))
    allocate(occ_orb_dn(ndn))
    allocate(num_occ_up_by_sym(n_group_elements))
    allocate(num_occ_dn_by_sym(n_group_elements))
    allocate(occ_orb_up_by_sym(n_group_elements,(nup+1)))
    allocate(occ_orb_dn_by_sym(n_group_elements,(ndn+1)))
    allocate(occ_orb_up_by_sym_copy(n_group_elements,(nup+1)))
    allocate(occ_orb_dn_by_sym_copy(n_group_elements,(ndn+1)))
    allocate(sym_sum_cs_sqrt(n_group_elements,norb))
    allocate(cs_sqrt_orb(norb))
    allocate(cs_sqrt_prime(nup+ndn))
    allocate(cs_sqrt_prime_spin(nup+ndn))
    allocate(sqrt_integrals(norb,norb))

    num_orb_by_sym       = 0
    which_orb_by_sym     = 0
    sym_sum_cs_sqrt      = 0
    cs_sqrt_orb          = 0
    cs_sqrt_prime        = 0
    cs_sqrt_prime_spin   = 0
    sqrt_integrals       = 0

    !Set up table of orbitals by symmetry
    do i = 1, norb
       sym1                                        = orbital_symmetries(i)
       num_orb_by_sym(sym1)                        = num_orb_by_sym(sym1) + 1
       which_orb_by_sym(sym1,num_orb_by_sym(sym1)) = i
    enddo


    !Set up sums of sqrt of integrals for Cauchy Schwarz sampling
    if(proposal_method.eq.'CauchySchwarz') then
      do i=1,norb
         do j=1,norb
            sym1                    = orbital_symmetries(j)
            !Check for negative integrals
            int_index               = integral_index(i,j,i,j)
            if(integrals(int_index).lt.-1e-6) then
               write(6,*) 'Negative integrals!'
               stop
            endif
            if(integrals(int_index).lt.0) integrals(int_index) = 0
            sqrt_integrals(i,j)     = sqrt(integrals(int_index))
            sym_sum_cs_sqrt(sym1,i) = sym_sum_cs_sqrt(sym1,i) + sqrt_integrals(i,j)
            cs_sqrt_orb(i)          = cs_sqrt_orb(i) + sqrt_integrals(i,j)
         enddo
      enddo
    endif

    uniform_sampling = .false.
  end subroutine setup_orb_by_symm

!===============================================================================================================
  subroutine off_diagonal_move_chem_cauchySchwarz(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j, importance_sampling, iwalk_child)
    !----------------------------------------------------------------------------------------------------------
    ! Description : Perform an offdiagonal move in SQMC for a slater determinant starting at det_i
    !               which is split into det_i_up and det_i_dn. det_i has weight = +-1.
    !               Probability is calculated using Cauchy-Schwarz from Adam's notes.
    !
    ! Created     : Matthew Otten, 21 January 2015
    !
    !-----------------------------------------------------------------------------------------------------------
    use more_tools, only : binary_search
    use common_psi_t, only : psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_g,psi_g_epsilon,psi_g_epsilon_inv
    implicit none

    !dummy variables
!   integer, intent(out) :: weight_j
    integer, intent(in)  :: iwalk_child
    real(rk), intent(out) :: weight_j
   !real(rk),optional, intent(out) :: diagonal_matrix_element
    integer,optional,intent(in) :: importance_sampling
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(out) :: det_j_up
    type(ik_vec), intent(out) :: det_j_dn
    type(ik_vec) :: tmp_det
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(out) :: det_j_up
    integer(ik), intent(out) :: det_j_dn
    integer(ik) :: tmp_det
#endif

    !local variables
    !excite_from_1_i and excite_from_2_i are both indices, not abolute orbital numbers
    integer :: excite_from_1_i, excite_from_2_i, excite_to_1,excite_to_2, excite_level, excite_level_sym, spatial_index1, spatial_index2
    integer :: n_single, n_double, n_total, sym1, sym2, sym3
    integer :: n_double_up, n_double_dn, n_double_both
    integer :: tot_spin_from !in units of 1/2
    integer :: i, i_open, j, tmp_orb
    integer,save :: n_occ_up,n_occ_dn
    ! 100 -> number of symmetries, norb -> maximum of norb of a specific symmetry
    real(rk), parameter :: epsilon = 1.e-10_rk
!   real(rk) :: ran, rannyu, proposal_prob, acceptance_prob, matrix_element, temp1,temp2
    real(rk) :: proposal_prob, proposal_prob_sym, matrix_element,matrix_element1,matrix_element2, temp1
    real(rk) :: proposal_prob_single, proposal_prob_double
    real(rk) :: norm_i,norm_j
    real(rk) :: temp_prob,electron_prob,denominator
    real(rk) :: rannyu,rand_number
!   logical ::  too_excited
    integer ::  nnzero
    logical ::  are_they_connected
    logical ::  spr,debug,foundOrb
    logical :: imp
    integer :: location
    real(rk) :: guiding_wf_ratio

    imp=.false.
    if (present(importance_sampling)) then
        if (importance_sampling .eq. 1) imp=.true.
    endif

    !occ_orb_up_by_sym_copy = 0
    !occ_orb_dn_by_sym_copy = 0
!    if(iwalk_child.eq.1) then
       !Setup matrices / vectors
       !We only have to do this once per determinant; so, only the first
       !time it is called (iwalk_child=1) per determinant
       n_occ_up          = 0
       n_occ_dn          = 0
       occ_orb_up        = 0
       occ_orb_dn        = 0
       num_occ_up_by_sym = 0
       num_occ_dn_by_sym = 0

       tmp_det           = det_i_up

       do while (tmp_det > 0)
          i                                               = trailz(tmp_det)+1
          sym1                                            = orbital_symmetries(i)
          num_occ_up_by_sym(sym1)                         = num_occ_up_by_sym(sym1) + 1
          n_occ_up                                        = n_occ_up+1
          occ_orb_up_by_sym(sym1,num_occ_up_by_sym(sym1)) = i
          occ_orb_up(n_occ_up)                            = i
          tmp_det                                         = ibclr(tmp_det,i-1)
       enddo

       tmp_det = det_i_dn
       do while (tmp_det > 0)
          i                                               = trailz(tmp_det)+1
          sym1                                            = orbital_symmetries(i)
          num_occ_dn_by_sym(sym1)                         = num_occ_dn_by_sym(sym1) + 1
          n_occ_dn                                        = n_occ_dn+1
          occ_orb_dn_by_sym(sym1,num_occ_dn_by_sym(sym1)) = i
          occ_orb_dn(n_occ_dn)                            = i
          tmp_det                                         = ibclr(tmp_det,i-1)
       enddo

       !Store all occupied orbitals in one vector
       do i =1,n_occ_up
          occ_orb(i) = occ_orb_up(i)
       enddo

       do i=1,n_occ_dn
          occ_orb(i+n_occ_up) = occ_orb_dn(i)
       enddo

       n_occ_up = n_occ_up - n_core_orb
       n_occ_dn = n_occ_dn - n_core_orb

!    endif


    !debug=.true.
    debug=.false.
    !spr=.true.
    spr=.false.

    excite_level = 0
    excite_from_1_i = 0
    excite_from_2_i = 0
    excite_to_1 = 0
    excite_to_2 = 0
    proposal_prob = 1._rk
    nnzero=0
    !initialize det_j to det_i and weight_j to 0

    det_j_up = det_i_up
    det_j_dn = det_i_dn
    weight_j = 0
    if(spr) write (6,*) "Incoming determinant",det_i_up,det_i_dn
    if(debug) write (6,*) "Incoming determinant",det_i_up,det_i_dn

    if (time_sym) then
        if (det_i_up .eq. det_i_dn) then
             if (z .eq. 1) then
                norm_i=sqrt2
             else
                norm_i=0._rk
                write (6,*) "This state is inconsistent with the time symmetry you put in and should never enter this routine as an incoming state"
                stop
             endif
        else
             norm_i=1._rk
        endif
    endif

    n_single = (nup-n_core_orb)*(norb - nup) + (ndn-n_core_orb)*(norb - ndn)
    n_double_up = (nup-n_core_orb)*(nup -n_core_orb-1)*(norb-nup)*(norb-nup-1)/4
    n_double_dn = (ndn-n_core_orb)*(ndn -n_core_orb-1)*(norb-ndn)*(norb-ndn-1)/4
    n_double_both = (nup-n_core_orb)*(norb-nup)*(ndn-n_core_orb)*(norb-ndn)

    n_double = n_double_up + n_double_dn + n_double_both
    n_total = n_single + n_double

    proposal_prob_single=(n_single)/real(n_total,rk)
    proposal_prob_double=(n_double)/real(n_total,rk)        ! or 1 - proposal_prob_single

    if (random_int(n_total)> n_single) then
       excite_level = 2
       proposal_prob = n_double /real(n_total,rk)

       !Calculate S'_2 = S_2(i) - sum_{j.ne.i,j occ} sqrt(<ij|ij>)
       !S'_2 is stored as all of the occupied up orbitals first,
       !followed by the down orbitals.
       !i.e., when looping through up orbitals, use i
       !when looping through down orbitals, use i+n_occ_up

       !First do up electrons
       do i = 1,nup
          cs_sqrt_prime(i) = cs_sqrt_orb(occ_orb_up(i))
          do j = 1,nup
             cs_sqrt_prime(i) = cs_sqrt_prime(i) - &
                     &sqrt_integrals(occ_orb_up(i),occ_orb_up(j))
          enddo
       enddo

       !Now do down electrons
       do i = 1,ndn
          cs_sqrt_prime(i+nup) = cs_sqrt_orb(occ_orb_dn(i))
          do j = 1,ndn
             cs_sqrt_prime(i+nup) = cs_sqrt_prime(i+nup) - &
                  &sqrt_integrals(occ_orb_dn(i),occ_orb_dn(j))
          enddo
       enddo

       !Calculate S'_2 for spin flip
       !S'_2spin = 2*S_2(i) + - sum_{js.ne.is,j occ} sqrt(<ij|ij)
       !but, now j is looped over BOTH spins!
       do i = 1,nup+ndn
          cs_sqrt_prime_spin(i) = 2*cs_sqrt_orb(occ_orb(i))
          do j = 1,nup+ndn
             cs_sqrt_prime_spin(i) = cs_sqrt_prime_spin(i) - &
                  &sqrt_integrals(occ_orb(i),occ_orb(j))
          enddo
       enddo

       !Find the sum for normalization of the probability

       sum_cs_sqrt_prime = 0
       do i=1+n_core_orb,nup
          sum_cs_sqrt_prime = sum_cs_sqrt_prime + cs_sqrt_prime_spin(i)
       enddo

       do i=1+nup+n_core_orb,nup+ndn
          sum_cs_sqrt_prime = sum_cs_sqrt_prime + cs_sqrt_prime_spin(i)
       enddo


       !Pick electron i with probability
       ! P_i = S'_2(i)/sum_k S'_2(k)

       !Draw a random number
       rand_number = rannyu()
      !call random_number(rand_number)

       !Search for which electron was selected
       electron_prob = 0
       foundOrb      = .false.
       do i = 1+n_core_orb,nup

          if(uniform_sampling) then
             electron_prob = electron_prob + 1/(real(n_occ_up+n_occ_dn,rk))!Temporary test of uniform prob
          else
             electron_prob = electron_prob + cs_sqrt_prime_spin(i)/sum_cs_sqrt_prime
          endif

          !Is it < or <=?
          if(rand_number.le.electron_prob) then
             excite_from_1_i = i
             foundOrb = .true.
             exit
          endif
       enddo

       if(foundOrb.eqv..false.) then
          do i = 1+nup+n_core_orb,nup+ndn

             if(uniform_sampling) then
                electron_prob = electron_prob + 1/(real(n_occ_up+n_occ_dn,rk))!Temporary test of uniform prob
             else
                electron_prob = electron_prob + cs_sqrt_prime_spin(i)/sum_cs_sqrt_prime
             endif

             !Is it < or <=?
             if(rand_number.le.electron_prob) then
                excite_from_1_i = i
                exit
             endif
          enddo
       endif

       !Translate excite_from_1_i and store spin
       tot_spin_from = 0
       if (excite_from_1_i > nup) then
          !spin down
          tot_spin_from = tot_spin_from - 1
!          excite_from_1_i = excite_from_1_i + 2*n_core_orb
       else
          !spin up
          tot_spin_from = tot_spin_from + 1
!          excite_from_1_i = excite_from_1_i + n_core_orb
       endif

       !Now, pick 2nd electron j with probability
       ! P(j|i) = S'_2(j) / (sum_k S'_2(k) - S'_2(i))

       !Draw a random number
       call random_number(rand_number)

       !Search for which electron was selected
       electron_prob = 0
       if(uniform_sampling) then
          do i = 1,n_occ_up+n_occ_dn-1!-1 might be temporary?
             electron_prob = electron_prob + 1/(real(n_occ_up+n_occ_dn-1,rk))!Temporary test of uniform prob

             !Is it < or <=?
             if(rand_number.le.electron_prob) then
                excite_from_2_i = i
                exit
             endif
          enddo
          if(tot_spin_from.eq.1) then
             if(excite_from_2_i.eq.excite_from_1_i-n_core_orb) excite_from_2_i = nelec-2*n_core_orb
          else
             if(excite_from_2_i.eq.excite_from_1_i-2*n_core_orb) excite_from_2_i = nelec-2*n_core_orb
          endif
          if (excite_from_2_i > nup-n_core_orb) then
             !spin down
             excite_from_2_i = excite_from_2_i + 2*n_core_orb
          else
             !spin up
             excite_from_2_i = excite_from_2_i + n_core_orb
          endif

       else

          foundOrb = .false.
          do i = 1+n_core_orb,nup
             !excite_from_1_i is still an index at this point; so, compare indices
             if(i.ne.excite_from_1_i) then
                electron_prob = electron_prob + cs_sqrt_prime_spin(i)&
                     &/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_1_i))
                !Is it < or <=?
                if(rand_number.le.electron_prob) then
                   excite_from_2_i = i
                   foundOrb = .true.
                   exit
                endif
             endif
          enddo

          if(foundOrb.eqv..false.) then
             do i = 1+nup+n_core_orb,nup+ndn
                !excite_from_1_i is still an index at this point; so, compare indices
                if(i.ne.excite_from_1_i) then
                   electron_prob = electron_prob + cs_sqrt_prime_spin(i)&
                        &/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_1_i))
                   !Is it < or <=?
                   if(rand_number.le.electron_prob) then
                      excite_from_2_i = i
                      exit
                   endif
                endif
             enddo
          endif
       endif

       !Translate excite_from_2_i and store spin
       if (excite_from_2_i > nup) then
          !spin down
          tot_spin_from = tot_spin_from - 1
!          excite_from_2_i = excite_from_2_i + 2*n_core_orb
       else
          !spin up
          tot_spin_from = tot_spin_from + 1
!          excite_from_2_i = excite_from_2_i + n_core_orb
       endif

       !Am I consistently doing this, with n_core_orb.ne.0?

       !We could have chosen them the reverse way, so calculate probability as
       ! P(i,j) = P(i)P(j|i) + P(j)P(i|j)

       if(uniform_sampling) then
       !Temporary test of uniform probability
          proposal_prob = proposal_prob*2.0 /(1.0_rk*(nelec-2*n_core_orb)*(nelec-2*n_core_orb-1))
       else
          proposal_prob = proposal_prob*(cs_sqrt_prime_spin(excite_from_1_i)/sum_cs_sqrt_prime*&
               &cs_sqrt_prime_spin(excite_from_2_i)/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_1_i)) +&
               &cs_sqrt_prime_spin(excite_from_2_i)/sum_cs_sqrt_prime*&
               &cs_sqrt_prime_spin(excite_from_1_i)/(sum_cs_sqrt_prime-cs_sqrt_prime_spin(excite_from_2_i)))

       endif

    else

       excite_level = 1
       proposal_prob = proposal_prob * n_single /(n_total*1.0_rk)
       ! excitation now is a single excitation
       ! choose electron to excite from at random
       excite_from_1_i = random_int(nelec-2*n_core_orb)
       ! calculate if the spin is up or down
       !tot_spin_from = 0

       if (excite_from_1_i > nup-n_core_orb) then
          !spin down
          tot_spin_from = -1
          excite_from_1_i = excite_from_1_i + 2*n_core_orb
       else
          !spin up
          tot_spin_from = 1
          excite_from_1_i = excite_from_1_i + n_core_orb
       endif

    endif

   ! remove excited electrons in det_j
    ! set excite_from_2_i > excite_from_1_i
    i =  excite_from_1_i + excite_from_2_i
    excite_from_1_i = i - max(excite_from_1_i,excite_from_2_i)
    excite_from_2_i = i - excite_from_1_i

    sym1 = 1 ! this variable will store total symmetry of the excitation.
    ! for the single case, just the symmetry of excite_from_1_i
    ! for the double case, product of symmetries of excite_from_1_i and excite_from_2_i
    !calculate total spin of the excitation
    !pick single or double excitation

    ! Edited by MJO on 16 Jan 2015. Sped up by using new data structures.
    if(excite_level==1) then
       !Because the excite_from* were reorded, and are initially set to 0,
       !if we have a single excitation, then excite_from_1_i is 0! So,
       !use excite_from_2_i, instead
       i = occ_orb(excite_from_2_i)

       if(excite_from_2_i.le.nup) then
          det_j_up=ibclr(det_j_up, i-1)
       else
          det_j_dn=ibclr(det_j_dn, i-1)
       endif
       sym1 = orbital_symmetries(i)

    else !excite_level==2
       i = occ_orb(excite_from_1_i)

       if(excite_from_1_i.le.nup) then
          det_j_up=ibclr(det_j_up, i-1)
       else
          det_j_dn=ibclr(det_j_dn, i-1)
       endif

       sym1 = orbital_symmetries(i)

       i = occ_orb(excite_from_2_i)
       if(excite_from_2_i.le.nup) then
          det_j_up=ibclr(det_j_up, i-1)
       else
          det_j_dn=ibclr(det_j_dn, i-1)
       endif

       sym1 = product_table(sym1,orbital_symmetries(i))
    endif

    if(debug.eqv..true.)print*,'before excite',excite_from_1_i,excite_from_2_i,proposal_prob

! At this point, we have selected a single or a double excitation
! and removed appropriate electrons from det_j_up and/or det_j_dn
!----------------------------------------------------------------------------
    i_open = 0
    if (excite_level==1) then

       proposal_prob = proposal_prob /  (nelec-2*n_core_orb)
       if (tot_spin_from ==1) then
          if(debug.eqv..true.) write(6,*)'single,totspin1'
          if(spr.eqv..true.) write(6,*)'single,totspin1'
          i_open = num_orb_by_sym(sym1) - num_occ_up_by_sym(sym1)
          ! now select one of them for excite_to_1
          if (i_open.ne.0) then
             excite_to_1 = random_int(i_open)
             proposal_prob = proposal_prob / i_open
          else
             return  ! no valid move possible, return with zero weight_j
          endif

          do i=1,num_occ_up_by_sym(sym1)
             if(occ_orb_up_by_sym(sym1,i).le.which_orb_by_sym(sym1,excite_to_1)) then
                excite_to_1=excite_to_1+1
             else
                exit
             endif
          enddo

          i = which_orb_by_sym(sym1,excite_to_1)
          !add electron to this orbital
          det_j_up = ibset(det_j_up,i-1)
       else ! total spin -1
          if(debug.eqv..true.) write(6,*)'single,totspin-1'
          if(spr.eqv..true.) write(6,*)'single,totspin-1'
          i_open = num_orb_by_sym(sym1) - num_occ_dn_by_sym(sym1)
          ! now select one of them for excite_to_1
          if (i_open.ne.0) then
             excite_to_1 = random_int(i_open)
             proposal_prob = proposal_prob / i_open
          else
             return  ! no valid move possible, return with zero weight_j
          endif

          do i=1,num_occ_dn_by_sym(sym1)
             if(occ_orb_dn_by_sym(sym1,i).le.which_orb_by_sym(sym1,excite_to_1)) then
                excite_to_1=excite_to_1+1
             else
                exit
             endif
          enddo

          i = which_orb_by_sym(sym1,excite_to_1)
          !add electron to this orbital
          det_j_dn = ibset(det_j_dn,i-1)
       endif

       ! end of single excitation part
    else
       ! now filling up double excitation: excite_to_1 and excite_to_2

       if (tot_spin_from == 2) then

          !have to excite to spinup orbital
          if(debug.eqv..true.) write(6,*)'double,totspin2'
          if(spr.eqv..true.) write(6,*)'double,totspin2'

          !Choose unoccupied orbital k with same spin as electron i with probability
          !P(i->k || j->k) = (sqrt( <ik|ik>)+sqrt(<jk|jk>))/(S'_2(i)+S'_2(j))

          !Draw a random number
          call random_number(rand_number)

          electron_prob = 0
          j             = 1
          do i=1,norb
             !We have to skip occupied orbitals
             if(j.le.nup) then
                tmp_orb = occ_orb_up(j)!Just a temporary variable
             else
                tmp_orb = 0
             endif

             if(i.eq.tmp_orb) then
                   !skip this orbital, as it is occupied
                   j=j+1
             else
                if(uniform_sampling) then
                   !Temporary test of uniform probability
                   electron_prob = electron_prob + 1/real(norb-nup,rk)
                else
                   electron_prob = electron_prob + (sqrt_integrals(occ_orb(excite_from_1_i),i) +&
                        sqrt_integrals(occ_orb(excite_from_2_i),i))/&
                        &(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))
                endif
                !Is it < or <=?
                if(rand_number.le.electron_prob) then
                   excite_to_1 = i
                   exit
                endif
             endif
          enddo

          spatial_index1 = excite_to_1
          det_j_up       = ibset(det_j_up, excite_to_1-1)

          !Choose unoccupied orbital l with correct symmetry and same spin as electron j with prob
          !P(j->l|i->k || i->l|j->k) = (sqrt( <jl | jl>)+sqrt( <il|il>)) /
          !(G_s(j) - sum_{(m.ne.j),sym(m)=s}sqrt( <jm|jm>) + G_s(i) - sum_{(m.ne.i),sum(m)=s}sqrt(<im|im>))

          !Calculate symmetry
          !When d_inf_h is used, the irreducible representation is no longer its own inverse
          if(d_infinity_h==1) then
             sym2 = product_table(get_inverse_dih(orbital_symmetries(spatial_index1)),sym1)
          else
             sym2 = product_table(orbital_symmetries(spatial_index1),sym1)
          endif


          i_open = num_orb_by_sym(sym2) - num_occ_up_by_sym(sym2)
          if(sym2.eq.orbital_symmetries(spatial_index1)) i_open = i_open-1
          if(i_open.eq.0) return ! no valid move possible
          !Temporary test of uniform probability
          if(uniform_sampling) temp1 = 1.0_rk / i_open ! to be later used for proposal probability

          !Draw a random number
          call random_number(rand_number)
          !If the two new orbitals are of the same symmetry, we have to do a little extra
          !Because the number of occupied orbitals is now one higher.
          if(sym2.eq.orbital_symmetries(spatial_index1)) then
             do i=1,num_occ_up_by_sym(sym2)
                occ_orb_up_by_sym_copy(sym2,i) = occ_orb_up_by_sym(sym2,i)
             enddo

             foundOrb=.false.
             do i=1,num_occ_up_by_sym(sym2)
                if(occ_orb_up_by_sym_copy(sym2,i).ge.spatial_index1) then
                   do j=num_occ_up_by_sym(sym2),i,-1
                      occ_orb_up_by_sym_copy(sym2,j+1) = occ_orb_up_by_sym_copy(sym2,j)
                   enddo
                   occ_orb_up_by_sym_copy(sym2,i)   = spatial_index1
                   foundOrb=.true.
                   exit
                endif
             enddo
             if(.not.foundOrb) occ_orb_up_by_sym_copy(sym2,num_occ_up_by_sym(sym2)+1) = spatial_index1
             !Calculate denominator first
             !Calculate sum
             electron_prob = 0 !Just a temporary variable
             do i = 1,num_occ_up_by_sym(sym2)+1
                electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_up_by_sym_copy(sym2,i))&
                     &+sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym_copy(sym2,i))
             enddo

             denominator = sym_sum_cs_sqrt(sym2,occ_orb(excite_from_2_i)) +&
                  &sym_sum_cs_sqrt(sym2,occ_orb(excite_from_1_i)) - electron_prob

             j             = 1
             excite_to_2   = 0
             electron_prob = 0
             do i=1,num_orb_by_sym(sym2)
                !This can probably be done more efficiently
                if(j.le.num_occ_up_by_sym(sym2)+1) then
                   tmp_orb = occ_orb_up_by_sym_copy(sym2,j)
                else
                   tmp_orb = 0
                endif
                if(which_orb_by_sym(sym2,i).eq.tmp_orb) then
                   j      = j+1
                else
                   if(uniform_sampling) then
                      !Temporary test of uniform probability
                      electron_prob = electron_prob + 1/(real(num_orb_by_sym(sym2) - num_occ_up_by_sym(sym2)-1,rk))
                   else
                      electron_prob = electron_prob + &
                           &(sqrt_integrals(occ_orb(excite_from_2_i),which_orb_by_sym(sym2,i))+&
                           &sqrt_integrals(occ_orb(excite_from_1_i),which_orb_by_sym(sym2,i)))&
                           &/denominator
                   endif
                   !Is it < or <=?
                   if(rand_number.le.electron_prob) then
                      excite_to_2 = which_orb_by_sym(sym2,i)
                      exit
                   endif
                endif
             enddo
          else
             !Calculate denominator first
             !Calculate sum
             electron_prob = 0 !Just a temporary variable
             do i = 1,num_occ_up_by_sym(sym2)
                electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_up_by_sym(sym2,i))&
                     & + sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym2,i))
             enddo

             denominator = sym_sum_cs_sqrt(sym2,occ_orb(excite_from_2_i)) +&
                  &sym_sum_cs_sqrt(sym2,occ_orb(excite_from_1_i))- electron_prob

             !Now find our new orbital
             j             = 1
             excite_to_2   = 0
             electron_prob = 0
             do i=1,num_orb_by_sym(sym2)

                !This can probably be done more efficiently
                if(j.le.num_occ_up_by_sym(sym2)) then
                   tmp_orb = occ_orb_up_by_sym(sym2,j)
                else
                   tmp_orb = 0
                endif

                if(which_orb_by_sym(sym2,i).eq.tmp_orb) then
                   j   = j+1
                else
                   !Temporary test of uniform probability
                   if(uniform_sampling) then
                      electron_prob = electron_prob + 1/(real(num_orb_by_sym(sym2)-num_occ_up_by_sym(sym2),rk))
                   else
                      electron_prob = electron_prob + &
                           &(sqrt_integrals(occ_orb(excite_from_2_i),which_orb_by_sym(sym2,i))&
                           &+sqrt_integrals(occ_orb(excite_from_1_i),which_orb_by_sym(sym2,i)))&
                           &/denominator
                   endif
                   !Is it < or <=?
                   if(rand_number.le.electron_prob) then
                      excite_to_2 = which_orb_by_sym(sym2,i)
                      exit
                   endif
                endif
             enddo
          endif

          spatial_index2 = excite_to_2
          det_j_up = ibset(det_j_up, excite_to_2-1)

          !Since these two unoccupied orbitals could have been selected in either order
          !We must calculate the probability for both pathways
          !P((i,j)->(k,l)|i,j) = 1/2 *[ P(i->k)P(j->l|i->k)
          ! + P(i->l)P(j->k|i->l)
          ! We don't need the 3rd and 4th terms right now, because we have a definite
          ! spawning order
          ! + P(j->k)P(i->l|j->k)
          ! + P(j->l)P(i->k|j->l)
          ! Note:
          ! i = occ_orb(excite_from_1_i)
          ! j = occ_orb(excite_from_2_i)
          ! k = excite_to_1
          ! l = excite_to_2

          !Calculate P(i->k)P(j->l|i->k)
          !Denominator was already calculated before
          temp_prob = (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)+&
               &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1))&
               &/(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))*&
               &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)+&
               &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2))&
               &/denominator


          !Calculate P(i->l)P(j->k|i->l)
          !If the the symmetries are the same, the denominator only changes a little

          if(sym2.eq.orbital_symmetries(spatial_index1)) then
             denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
                  &-sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
                  &+sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)&
                  &-sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)
          else
             sym3 = orbital_symmetries(spatial_index1)
             electron_prob = 0 !Just a temporary variable
             do i = 1,num_occ_up_by_sym(sym3)
                electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_up_by_sym(sym3,i))&
                     &+sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym3,i))
             enddo
             denominator = sym_sum_cs_sqrt(sym3,occ_orb(excite_from_2_i))&
                  &+sym_sum_cs_sqrt(sym3,occ_orb(excite_from_1_i))- electron_prob
          endif

          temp_prob = temp_prob + (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)+&
               &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2))&
               &/(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))*&
               &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)+&
               &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1))&
               &/denominator

          !Calculate P(j->k)P(i->l|j->k)
          !Calculate new sum for denominator; now, we need to use the symmetry of spatial_index1
          !and the G sum of J
          ! if(d_infinity_h==1) then
          !    i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_up_by_sym(product_table(get_inverse_dih(sym2),sym1))
          ! else
          !    i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_up_by_sym(product_table(sym2,sym1))
          ! endif

          ! if(i_open.ne.0) then
          !    electron_prob = 0 !Just a temporary variable
          !    do i = 1,num_occ_up_by_sym(orbital_symmetries(spatial_index1))
          !       electron_prob = electron_prob + &
          !            &sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(orbital_symmetries(spatial_index1),i))
          !    enddo

          !    denominator = sym_sum_cs_sqrt(orbital_symmetries(spatial_index1),occ_orb(excite_from_1_i)) - electron_prob
          ! !If the symmetries are the same, we need to subtract off one more integral
          ! ! if(sym2.eq.orbital_symmetries(spatial_index1)) then
          ! !    !excite_to_1 should never equal excite_from_1_i, so I don't need to check
          ! !    !for equality
          ! !       denominator = denominator - sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)
          ! ! endif

          ! ! temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
          ! !      &/cs_sqrt_prime(excite_from_2_i)*&
          ! !      &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)&
          ! !      &/denominator

          ! ! !Calculate P(j->l)P(i->k|j->l)
          ! ! !Calculate sum for denominator
          ! ! !Denominator was already calculated before; if the symmetry of k and
          ! ! !l are the same, then we need to add back in the contribution of k and
          ! ! !subtract out the contribution of l, since it was a sum over occupied orbitals
          ! ! !If the the symmetries are not the same, the denominator does not change
          !    if(sym2.eq.orbital_symmetries(spatial_index1)) then
          !       !excite_from_1_i should never equal excite_to_1 or excite_to_2
          !       !So, I don't need to check for equality
          !       denominator = denominator - sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)
          !    endif

          !    temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
          !         &/cs_sqrt_prime(excite_from_2_i)*&
          !         &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)&
          !         &/denominator


          if(uniform_sampling) then
             !Temporary test of uniform probability
             proposal_prob = proposal_prob / real(norb - nup,rk)
             if(d_infinity_h==1) then
                i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_up_by_sym(product_table(get_inverse_dih(sym2),sym1))
             else
                i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_up_by_sym(product_table(sym2,sym1))
             endif
             if(sym2.eq.orbital_symmetries(spatial_index1)) i_open = i_open - 1
             if (i_open==0) then
                ! not possible to excite otherwise
                proposal_prob = proposal_prob*temp1 ! correcting for 2 pathways to get this excitation
             else
                proposal_prob = proposal_prob*(temp1+(1.0_rk/i_open)) ! correcting for 2 pathways to get this excitation
             endif
          else
             proposal_prob = proposal_prob*temp_prob
          endif

          if(debug) write(6,*)'prop',proposal_prob,sym1,spatial_index2,spatial_index1,product_table(orbital_symmetries(spatial_index2),sym1),i_open
          if(debug) write(6,'(4B64)')det_i_up,det_i_dn,det_j_up,det_j_dn
          !write(6,*)'symmstuff',sym1,spatial_index1,spatial_index2,orbital_symmetries(spatial_index1),orbital_symmetries(spatial_index2),product_table(orbital_symmetries(spatial_index2),sym1)
       elseif (tot_spin_from == -2) then
          if(debug.eqv..true.) write(6,*)'double,totspin-2'
          if(spr.eqv..true.) write(6,*)'double,totspin-2'
          !Choose unoccupied orbital k with same spin as electron i with probability
          !P(i->k) = sqrt( <ik|ik>)/S'_2(i)

          !Draw a random number
          call random_number(rand_number)

          electron_prob = 0
          j             = 1
          do i=1,norb
             !We have to skip occupied orbitals
             if(j.le.ndn) then
                tmp_orb = occ_orb_dn(j)!Just a temporary variable
             else
                tmp_orb = 0
             endif
             if(i.eq.tmp_orb) then
                   !skip this orbital, as it is occupied
                   j=j+1
             else
                !Temporary test of uniform probability
                if(uniform_sampling) then
                   electron_prob = electron_prob + 1/real(norb-ndn,rk)
                else
                   electron_prob = electron_prob + (sqrt_integrals(occ_orb(excite_from_1_i),i)+&
                        &sqrt_integrals(occ_orb(excite_from_2_i),i))/&
                        &(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))
                endif
                !is it < or <=?
                if(rand_number.le.electron_prob) then
                   excite_to_1 = i
                   exit
                endif
             endif
          enddo

          spatial_index1 = excite_to_1
          det_j_dn       = ibset(det_j_dn, i-1)

          !Choose unoccupied orbital l with correct symmetry and same spin as electron j with prob
          !P(j->l|i->k || i->l|j->k) = (sqrt( <jl | jl>)+sqrt( <il|il>)) /
          !(G_s(j) - sum_{(m.ne.j),sym(m)=s}sqrt( <jm|jm>) + G_s(i) - sum_{(m.ne.i),sum(m)=s}sqrt(<im|im>))


          !Calculate symmetry
          !When d_inf_h is used, the irreducible representation is no longer its own inverse
          if(d_infinity_h==1) then
             sym2 = product_table(get_inverse_dih(orbital_symmetries(spatial_index1)),sym1)
          else
             sym2 = product_table(orbital_symmetries(spatial_index1),sym1)
          endif

          i_open = num_orb_by_sym(sym2) - num_occ_dn_by_sym(sym2)
          if(sym2.eq.orbital_symmetries(spatial_index1)) i_open = i_open-1
          if(i_open.eq.0) return !No moves possible

          if(uniform_sampling) temp1 = 1.0_rk / i_open ! Temporary test of uniform

          !Draw a random number
          call random_number(rand_number)
          !If the two new orbitals are of the same symmetry, we have to do a little extra
          !Because the number of occupied orbitals is now one higher.
          if(sym2.eq.orbital_symmetries(spatial_index1)) then
             do i=1,num_occ_dn_by_sym(sym2)
                occ_orb_dn_by_sym_copy(sym2,i) = occ_orb_dn_by_sym(sym2,i)
             enddo

             foundOrb=.false.
             do i=1,num_occ_dn_by_sym(sym2)
                if(occ_orb_dn_by_sym_copy(sym2,i).ge.spatial_index1) then
                   do j=num_occ_dn_by_sym(sym2),i,-1
                      occ_orb_dn_by_sym_copy(sym2,j+1) = occ_orb_dn_by_sym_copy(sym2,j)
                   enddo
                   occ_orb_dn_by_sym_copy(sym2,i)   = spatial_index1
                   foundOrb=.true.
                   exit
                endif
             enddo
             if(.not.foundOrb) occ_orb_dn_by_sym_copy(sym2,num_occ_dn_by_sym(sym2)+1) = spatial_index1
             !Calculate denominator first
             !Calculate sum
             electron_prob = 0 !Just a temporary variable
             do i = 1,num_occ_dn_by_sym(sym2)+1
                electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_dn_by_sym_copy(sym2,i))&
                     &+ sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym_copy(sym2,i))
             enddo

             denominator = sym_sum_cs_sqrt(sym2,occ_orb(excite_from_2_i))&
                  &+sym_sum_cs_sqrt(sym2,occ_orb(excite_from_1_i)) - electron_prob

             j             = 1
             excite_to_2   = 0
             electron_prob = 0
             do i=1,num_orb_by_sym(sym2)

                !This can probably be done more efficiently
                if(j.le.num_occ_dn_by_sym(sym2)+1) then
                   tmp_orb = occ_orb_dn_by_sym_copy(sym2,j)
                else
                   tmp_orb = 0
                endif
                if(which_orb_by_sym(sym2,i).eq.tmp_orb) then
                   j      = j+1
                else
                   !Temporary test of uniform probability
                   if(uniform_sampling) then
                      electron_prob = electron_prob + 1/(real(num_orb_by_sym(sym2)-num_occ_dn_by_sym(sym2)-1,rk))
                   else
                      electron_prob = electron_prob + &
                           &(sqrt_integrals(occ_orb(excite_from_2_i),which_orb_by_sym(sym2,i))+&
                           &sqrt_integrals(occ_orb(excite_from_1_i),which_orb_by_sym(sym2,i)))&
                           &/denominator
                   endif
                   !Is it < or <=?
                   if(rand_number.le.electron_prob) then
                      excite_to_2 = which_orb_by_sym(sym2,i)
                      exit
                   endif
                endif
             enddo
          else
             !Calculate denominator first
             !Calculate sum
             electron_prob = 0 !Just a temporary variable
             do i = 1,num_occ_dn_by_sym(sym2)
                electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_dn_by_sym(sym2,i))+&
                     &sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym(sym2,i))
             enddo

             denominator = sym_sum_cs_sqrt(sym2,occ_orb(excite_from_2_i))&
                  &+sym_sum_cs_sqrt(sym2,occ_orb(excite_from_1_i)) - electron_prob

             !Now find our new orbital
             j             = 1
             excite_to_2   = 0
             electron_prob = 0
             do i=1,num_orb_by_sym(sym2)

                !This can probably be done more efficiently
                if(j.le.num_occ_dn_by_sym(sym2)) then
                   tmp_orb = occ_orb_dn_by_sym(sym2,j)
                else
                   tmp_orb = 0
                endif

                if(which_orb_by_sym(sym2,i).eq.tmp_orb) then
                   j   = j+1
                else
                   if(uniform_sampling) then
                      !Temporary test of uniform probability
                      electron_prob = electron_prob + 1/(real(num_orb_by_sym(sym2)-num_occ_dn_by_sym(sym2),rk))
                   else
                      electron_prob = electron_prob + &
                           &(sqrt_integrals(occ_orb(excite_from_2_i),which_orb_by_sym(sym2,i))&
                           &+sqrt_integrals(occ_orb(excite_from_1_i),which_orb_by_sym(sym2,i)))&
                           &/denominator
                   endif
                   !is it < or <=?
                   if(rand_number.le.electron_prob) then
                      excite_to_2 = which_orb_by_sym(sym2,i)
                      exit
                   endif
                endif
             enddo
          endif


          spatial_index2 = excite_to_2
          det_j_dn = ibset(det_j_dn, excite_to_2-1)

          !Since these two unoccupied orbitals could have been selected in either order
          !We must calculate the probability for both pathways
          !P((i,j)->(k,l)|i,j) = 1/2 *[ P(i->k)P(j->l|i->k)
          ! + P(i->l)P(j->k|i->l)
          ! We do not need the 3rd and 4th terms right now
          ! + P(j->k)P(i->l|j->k)
          ! + P(j->l)P(i->k|j->l)
          ! Note:
          ! i = occ_orb(excite_from_1_i)
          ! j = occ_orb(excite_from_2_i)
          ! k = excite_to_1
          ! l = excite_to_2

          !Calculate P(i->k||j->l)P(j->l|i->k||i->l|j->k)
          !Denominator was already calculated before
          temp_prob = (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)+&
               &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1))&
               &/(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))*&
               &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)+&
               &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2))&
               &/denominator


          !Calculate P(i->l||j->l)P(j->k|i->l||i->k|j->l)
          !If the the symmetries are the same, the denominator only changes a little

          if(sym2.eq.orbital_symmetries(spatial_index1)) then
             denominator = denominator + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
                  &-sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
                  &+sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)&
                  &-sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)
          else
             sym3 = orbital_symmetries(spatial_index1)
             electron_prob = 0 !Just a temporary variable
             do i = 1,num_occ_dn_by_sym(sym3)
                electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_dn_by_sym(sym3,i))&
                     &+ sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym(sym3,i))
             enddo
             denominator = sym_sum_cs_sqrt(sym3,occ_orb(excite_from_2_i))+&
                  &sym_sum_cs_sqrt(sym3,occ_orb(excite_from_1_i)) - electron_prob
          endif

          temp_prob = temp_prob + (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)&
               &+sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2))&
               &/(cs_sqrt_prime(excite_from_1_i)+cs_sqrt_prime(excite_from_2_i))*&
               &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)+&
               &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1))&
               &/denominator

          !Calculate P(j->k)P(i->l|j->k)
          !Calculate new sum for denominator; now, we need to use the symmetry of i, not j
          !and the G sum of J
          ! if(d_infinity_h==1) then
          !    i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_dn_by_sym(product_table(get_inverse_dih(sym2),sym1))
          ! else
          !    i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_dn_by_sym(product_table(sym2,sym1))
          ! endif

          ! if(i_open.ne.0) then
          !    electron_prob = 0 !Just a temporary variable
          !    do i = 1,num_occ_dn_by_sym(orbital_symmetries(spatial_index1))
          !       electron_prob = electron_prob + &
          !            &sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym(orbital_symmetries(spatial_index1),i))
          !    enddo

          !    denominator = sym_sum_cs_sqrt(orbital_symmetries(spatial_index1),occ_orb(excite_from_1_i)) - electron_prob
          !If the symmetries are the same, we need to subtract off one more integral
          ! if(sym2.eq.orbital_symmetries(spatial_index1)) then
          !    !excite_to_1 should never equal excite_from_1_i, so I don't need to check
          !    !for equality
          !       denominator = denominator - sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)
          ! endif

          ! temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
          !      &/cs_sqrt_prime(excite_from_2_i)*&
          !      &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)&
          !      &/denominator

          ! !Calculate P(j->l)P(i->k|j->l)
          ! !Calculate sum for denominator
          ! !Denominator was already calculated before; if the symmetry of k and
          ! !l are the same, then we need to add back in the contribution of k and
          ! !subtract out the contribution of l, since it was a sum over occupied orbitals
          ! !If the the symmetries are not the same, the denominator does not change
             ! if(sym2.eq.orbital_symmetries(spatial_index1)) then
             !    !excite_from_1_i should never equal excite_to_1 or excite_to_2
             !    !So, I don't need to check for equality
             !    ! denominator = denominator + sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)&
             !    !      &-sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)
             !    denominator = denominator-sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)
             ! endif

             ! temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
             !      &/cs_sqrt_prime(excite_from_2_i)*&
             !      &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)&
             !      &/denominator


          !Temporary test of uniform probability
          if(uniform_sampling) then
             proposal_prob = proposal_prob /real(norb - ndn,rk)
             if(d_infinity_h==1) then
                i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_dn_by_sym(product_table(get_inverse_dih(sym2),sym1))
             else
                i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_dn_by_sym(product_table(sym2,sym1))
             endif
             if(sym2.eq.orbital_symmetries(spatial_index1)) i_open = i_open - 1

             !          if(debug) write(6,'(4B64)')det_i_up,det_i_dn,det_j_up,det_j_dn
             if (i_open==0) then
                ! not possible to excite otherwise
                proposal_prob = proposal_prob*temp1 ! correcting for 2 pathways to get this excitation
             else
                proposal_prob = proposal_prob*(temp1+(1.0_rk/i_open)) ! correcting for 2 pathways to get this excitation
             endif

          else
             proposal_prob = proposal_prob*temp_prob
          endif

          if(debug) write(6,*)'prop',proposal_prob,sym1,spatial_index2,spatial_index1,product_table(orbital_symmetries(spatial_index2),sym1),i_open
          if(debug) write(6,'(4B64)')det_i_up,det_i_dn,det_j_up,det_j_dn
       else
          !up-down case of a double excitation
          !Temporary test of uniform probability
          if(uniform_sampling) proposal_prob = proposal_prob*  1.0_rk/(2*norb - ndn-nup)
          sym2 =1

          !Draw a random number
          call random_number(rand_number)

          electron_prob = 0
          j             = 1
          foundOrb      = .false.
          !First, loop over up orbitals
          do i=1,norb
             !We have to skip occupied orbitals
             if(j.le.nup) then
                tmp_orb = occ_orb_up(j)!Just a temporary variable
             else
                tmp_orb = 0
             endif
             if(i.eq.tmp_orb) then
                   !skip this orbital, as it is occupied
                   j=j+1
             else
                !Temporary test of uniform probability
                if(uniform_sampling) then
                   electron_prob = electron_prob + 1/real(2*norb - nup - ndn,rk)
                else
                   electron_prob = electron_prob + (sqrt_integrals(occ_orb(excite_from_1_i),i)+&
                        &sqrt_integrals(occ_orb(excite_from_2_i),i))/&
                        &(cs_sqrt_prime_spin(excite_from_1_i)+cs_sqrt_prime_spin(excite_from_2_i))
                endif
                !is it < or <=?
                if(rand_number.le.electron_prob) then
                   excite_to_1 = i
                   foundOrb    = .true.
                   exit
                endif
             endif
          enddo

          j = 1
          !First, now, loop over down orbitals, if we didn't find it before
          if(foundOrb.eqv..false.) then
             do i=1,norb
                !We have to skip occupied orbitals
                if(j.le.ndn) then
                   tmp_orb = occ_orb_dn(j)!Just a temporary variable
                else
                   tmp_orb = 0
                endif
                if(i.eq.tmp_orb) then
                   !skip this orbital, as it is occupied
                   j=j+1
                else
                   !Temporary test of uniform probability
                   if(uniform_sampling) then
                      electron_prob = electron_prob + 1/real(2*norb - nup - ndn,rk)
                   else
                      electron_prob = electron_prob + (sqrt_integrals(occ_orb(excite_from_1_i),i)+&
                           &sqrt_integrals(occ_orb(excite_from_2_i),i))/&
                           &(cs_sqrt_prime_spin(excite_from_1_i)+cs_sqrt_prime_spin(excite_from_2_i))
                   endif
                   !Is it < or <=?
                   if(rand_number.le.electron_prob) then
                      excite_to_1 = i + norb
                      exit
                   endif
                endif
             enddo
          endif
          if (excite_to_1.le.norb) then
             !Excite to up electron first
             if(debug.eqv..true.) write(6,*)'double,totspin1'
             if(spr.eqv..true.) write(6,*)'double,totspin1'

             ! proposal_prob = proposal_prob * sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)/&
             !      &(cs_sqrt_prime_spin(excite_from_1_i))

             det_j_up = ibset(det_j_up, excite_to_1-1)
             spatial_index1 = excite_to_1
             sym2 = orbital_symmetries(excite_to_1)

             !try to find a spin_down orbital that conserves symmetry
             !When d_inf_h is used, the irreducible representation is no longer its own inverse
             if(d_infinity_h==1) then
                sym2 = product_table(sym1,get_inverse_dih(sym2))
             else
                sym2 = product_table(sym1,sym2)
             endif

             i_open = num_orb_by_sym(sym2) - num_occ_dn_by_sym(sym2)
             if(i_open.eq.0) return !No moves possible
             temp1 = 1.0_rk / i_open

             !Draw a random number
             call random_number(rand_number)
             !We DON'T have to check for same symmetry, because the up hole
             !will not make less available moves for the down hole
             !Calculate denominator first
             !Calculate sum
             electron_prob = 0 !Just a temporary variable

             do i = 1,num_occ_dn_by_sym(sym2)
                electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_dn_by_sym(sym2,i))&
                     &+ sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym(sym2,i))
             enddo

             denominator = sym_sum_cs_sqrt(sym2,occ_orb(excite_from_2_i))&
                  &+sym_sum_cs_sqrt(sym2,occ_orb(excite_from_1_i)) - electron_prob

             !Now find our new orbital
             j             = 1
             excite_to_2   = 0
             electron_prob = 0

             do i=1,num_orb_by_sym(sym2)
                !This can probably be done more efficiently
                if(j.le.num_occ_dn_by_sym(sym2)) then
                   tmp_orb = occ_orb_dn_by_sym(sym2,j)
                else
                   tmp_orb = 0
                endif

                if(which_orb_by_sym(sym2,i).eq.tmp_orb) then
                   j   = j+1
                else
                   !Temporary test of uniform probability
                   if(uniform_sampling) then
                      electron_prob = electron_prob + 1/(real(num_orb_by_sym(sym2)-num_occ_dn_by_sym(sym2),rk))
                   else
                      electron_prob = electron_prob + &
                           &(sqrt_integrals(occ_orb(excite_from_2_i),which_orb_by_sym(sym2,i))+&
                           &sqrt_integrals(occ_orb(excite_from_1_i),which_orb_by_sym(sym2,i)))&
                           &/denominator
                   endif
                   !is it < or <=?
                   if(rand_number.le.electron_prob) then
                      excite_to_2 = which_orb_by_sym(sym2,i)
                      exit
                   endif
                endif
             enddo


             spatial_index2 = excite_to_2
             det_j_dn = ibset(det_j_dn, excite_to_2-1)

             ! Need to calculate proposal probaility for this move

             !Since these two unoccupied orbitals could have been selected in either order
             !We must calculate the probability for both pathways
             !P((i,j)->(k,l)|i,j) = 1/2 *[ P(i->k)P(j->l|i->k)
             ! We do not need the 2nd and 3rd terms right now
             ! + P(i->l)P(j->k|i->l)
             ! + P(j->k)P(i->l|j->k)
             ! + P(j->l)P(i->k|j->l)
             ! Note:
             ! i = occ_orb(excite_from_1_i)
             ! j = occ_orb(excite_from_2_i)
             ! k = excite_to_1 = up
             ! l = excite_to_2 = down


          !Calculate P(i->k||j->l)P(j->l|i->k||i->l|j->k)
          !Denominator was already calculated before
             temp_prob = (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)+&
                  &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1))&
                  &/(cs_sqrt_prime_spin(excite_from_1_i)+cs_sqrt_prime_spin(excite_from_2_i))*&
                  &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)+&
                  &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2))&
                  &/denominator

             ! temp_prob = sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
             !      &/denominator


             !Calculate P(i->l)P(j->k|i->l)

             !Calculate new denominator first - because we had opposite spins!
             !Calculate sum
             electron_prob = 0 !Just a temporary variable
             sym3 = orbital_symmetries(spatial_index1)
             i_open = num_orb_by_sym(sym3) - num_occ_up_by_sym(sym3)
             if(i_open.ne.0) then
                do i = 1,num_occ_up_by_sym(sym3)
                   electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_up_by_sym(sym3,i))&
                        &+sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym3,i))
                enddo

                denominator = sym_sum_cs_sqrt(sym3,occ_orb(excite_from_2_i))&
                     &+sym_sum_cs_sqrt(sym3,occ_orb(excite_from_1_i)) - electron_prob

                temp_prob = temp_prob + (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)+&
                     &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2))&
                     &/(cs_sqrt_prime_spin(excite_from_1_i)+cs_sqrt_prime_spin(excite_from_2_i))*&
                     &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)+&
                     &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1))&
                     &/denominator
             endif

             ! temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
             !      &/denominator

             !calculate P(j->k)P(i->l|j->k)
             !Calculate new sum for denominator; now, we need to use the symmetry of spatial_index1
             !and the G sum of j
             ! electron_prob = 0 !Just a temporary variable
             ! do i = 1,num_occ_up_by_sym(orbital_symmetries(spatial_index1))
             !    if(occ_orb_up_by_sym(orbital_symmetries(spatial_index1),i).ne.occ_orb(excite_from_1_i)) then
             !       electron_prob = electron_prob + &
             !            &sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(orbital_symmetries(spatial_index1),i))
             !    endif
             ! enddo

             ! denominator = sym_sum_cs_sqrt(orbital_symmetries(spatial_index1),occ_orb(excite_from_1_i)) - electron_prob

             ! temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
             !      &/cs_sqrt_prime_spin(excite_from_2_i)*&
             !      &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)&
             !      &/denominator

             !Calculate P(j->l)P(i->k|j->l)
             ! if(d_infinity_h==1) then
             !    i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_up_by_sym(product_table(get_inverse_dih(sym2),sym1))
             ! else
             !    i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_up_by_sym(product_table(sym2,sym1))
             ! endif

             ! if(i_open.ne.0) then
             !    sym3 = orbital_symmetries(spatial_index1)
             !    electron_prob = 0 !Just a temporary variable
             !    do i = 1,num_occ_up_by_sym(sym3)
             !       electron_prob = electron_prob + &
             !            &sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym3,i))
             !    enddo

             !    denominator = sym_sum_cs_sqrt(sym3,occ_orb(excite_from_1_i)) - electron_prob

             !    temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
             !         &/cs_sqrt_prime_spin(excite_from_2_i)*&
             !         &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)&
             !         &/denominator
             ! endif


             !Temporary test of uniform probability
            ! Need to modify proposal probability for this move. For this, we calculate the probability of opposite excitation
            ! i.e. the excitation in which the electron at spatial_index2 was filled first
            ! sym2   = product_table(orbital_symmetries(spatial_index2),sym1)
             if(uniform_sampling) then
                if(d_infinity_h==1) then
                   i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_up_by_sym(product_table(get_inverse_dih(sym2),sym1))
                else
                   i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_up_by_sym(product_table(sym2,sym1))
                endif

                if (i_open.ne.0) then
                   proposal_prob = proposal_prob*(temp1 + (1.0_rk/i_open))! correcting for 2 pathways to get this excita
                else
                   proposal_prob = proposal_prob*temp1! correcting for 2 pathways to get this excita
                endif
             else
                proposal_prob = proposal_prob*temp_prob
             endif

             if(debug.eqv..true.) write(6,*)'prop',proposal_prob,sym1,spatial_index2,spatial_index1,product_table(orbital_symmetries(spatial_index2),sym1),i_open
             if(debug) write(6,'(4B64)')det_i_up,det_i_dn,det_j_up,det_j_dn
          else
             !Down orbital first
             if(debug.eqv..true.) write(6,*)'double,totspin-1'
             if(spr.eqv..true.) write(6,*)'double,totspin-1'

             excite_to_1 = excite_to_1 - norb

             ! proposal_prob = proposal_prob * sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)/&
             !      &(cs_sqrt_prime_spin(excite_from_1_i))

             !excite to spin-down orbital
             det_j_dn = ibset(det_j_dn, excite_to_1-1)
             spatial_index1 = excite_to_1
             sym2 = orbital_symmetries(excite_to_1)


             !try to find a spin_up orbital that conserves symmetry
             !When d_inf_h is used, the irreducible representation is no longer its own inverse
             if(d_infinity_h==1) then
                sym2 = product_table(sym1,get_inverse_dih(sym2))
             else
                sym2 = product_table(sym1,sym2)
             endif

             i_open = num_orb_by_sym(sym2) - num_occ_up_by_sym(sym2)
             if(i_open.eq.0) return !No moves possible
             temp1 = 1.0_rk/i_open

             !Draw a random number
             call random_number(rand_number)
             !We DON'T have to check for same symmetry, because the up hole
             !will not make less available moves for the down hole
             !Calculate denominator first
             !Calculate sum
             electron_prob = 0 !Just a temporary variable

             do i = 1,num_occ_up_by_sym(sym2)
                electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_up_by_sym(sym2,i))&
                     &+sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym2,i))
             enddo

             denominator = sym_sum_cs_sqrt(sym2,occ_orb(excite_from_2_i)) +&
                  &sym_sum_cs_sqrt(sym2,occ_orb(excite_from_1_i)) - electron_prob

             !Now find our new orbital
             j             = 1
             excite_to_2   = 0
             electron_prob = 0

             do i=1,num_orb_by_sym(sym2)
                !This can probably be done more efficiently
                if(j.le.num_occ_up_by_sym(sym2)) then
                   tmp_orb = occ_orb_up_by_sym(sym2,j)
                else
                   tmp_orb = 0
                endif

                if(which_orb_by_sym(sym2,i).eq.tmp_orb) then
                   j   = j+1
                else
                   !Temporary test of uniform probability
                   if(uniform_sampling) then
                      electron_prob = electron_prob + 1/real(num_orb_by_sym(sym2)-num_occ_up_by_sym(sym2),rk)
                   else
                      electron_prob = electron_prob + &
                           &(sqrt_integrals(occ_orb(excite_from_2_i),which_orb_by_sym(sym2,i))+&
                           &sqrt_integrals(occ_orb(excite_from_1_i),which_orb_by_sym(sym2,i)))&
                           &/denominator
                   endif
                   !is it < or <=?
                   if(rand_number.le.electron_prob) then
                      excite_to_2 = which_orb_by_sym(sym2,i)
                      exit
                   endif
                endif
             enddo

             spatial_index2 = excite_to_2

             det_j_up = ibset(det_j_up, excite_to_2-1)



             ! Need to calculate proposal probaility for this move

             !Since these two unoccupied orbitals could have been selected in either order
             !We must calculate the probability for both pathways
             !P((i,j)->(k,l)|i,j) = 1/2 *[ P(i->k)P(j->l|i->k)
             ! + P(i->l)P(j->k|i->l)
             ! We do not need the 3rd and 4th terms right now
             ! + P(j->k)P(i->l|j->k)
             ! + P(j->l)P(i->k|j->l)
             ! Note:
             ! i = occ_orb(excite_from_1_i)
             ! j = occ_orb(excite_from_2_i)
             ! k = excite_to_1 = down
             ! l = excite_to_2 = up

             !Calculate P(i->k||j->l)P(j->l|i->k||i->l|j->k)
             !Denominator was already calculated before
             temp_prob = (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)+&
                  &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1))&
                  &/(cs_sqrt_prime_spin(excite_from_1_i)+cs_sqrt_prime_spin(excite_from_2_i))*&
                  &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)+&
                  &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2))&
                  &/denominator

             ! temp_prob = sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
             !      &/denominator

             !Calculate new denominator - have to use different spin!
             !Calculate sum
             sym3 = orbital_symmetries(spatial_index1)
             i_open = num_orb_by_sym(sym3) - num_occ_dn_by_sym(sym3)
             if(i_open.ne.0) then
                electron_prob = 0 !Just a temporary variable
                do i = 1,num_occ_dn_by_sym(sym3)
                   electron_prob = electron_prob + sqrt_integrals(occ_orb(excite_from_2_i),occ_orb_dn_by_sym(sym3,i))&
                        &+sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym(sym3,i))
                enddo

                denominator = sym_sum_cs_sqrt(sym3,occ_orb(excite_from_2_i))&
                     &+sym_sum_cs_sqrt(sym3,occ_orb(excite_from_1_i))- electron_prob

                !Calculate P(i->l)P(j->k|i->l)
                temp_prob = temp_prob + (sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)+&
                     &sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2))&
                     &/(cs_sqrt_prime_spin(excite_from_1_i)+cs_sqrt_prime_spin(excite_from_2_i))*&
                     &(sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)+&
                     &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1))&
                     &/denominator
             endif

             ! temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
             !      &/denominator

             !calculate P(j->k)P(i->l|j->k)
             !Calculate new sum for denominator; now, we need to use the symmetry of i, not j
             !and the G sum of j
             ! sym1 = orbital_symmetries(spatial_index1)
             ! electron_prob = 0 !Just a temporary variable
             ! do i = 1,num_occ_up_by_sym(sym1)
             !    if(excite_from_1_i.le.n_occ_up) then
             !       if(occ_orb_up_by_sym(sym1,i).ne.occ_orb(excite_from_1_i)) then
             !          electron_prob = electron_prob + &
             !               &sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym1,i))
             !       endif
             !    else
             !          electron_prob = electron_prob + &
             !               &sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_up_by_sym(sym1,i))
             !    endif
             ! enddo

             ! if(excite_from_1_i.le.n_occ_up) then
             !    denominator = sym_sum_cs_sqrt(sym1,occ_orb(excite_from_1_i)) - electron_prob
             ! else
             !    if(orbital_symmetries(occ_orb(excite_from_1_i)).eq.sym1) then
             !       denominator = sym_sum_cs_sqrt(sym1,occ_orb(excite_from_1_i)) &
             !            &+ sqrt_integrals(occ_orb(excite_from_1_i),occ_orb(excite_from_1_i))&
             !            &- electron_prob
             !    else
             !       denominator = sym_sum_cs_sqrt(sym1,occ_orb(excite_from_1_i)) - electron_prob
             !    endif
             ! endif

             ! temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_1)&
             !      &/cs_sqrt_prime_spin(excite_from_2_i)*&
             !      &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_2)&
             !      &/denominator

             ! !Calculate P(j->l)P(i->k|j->l)
             ! if(d_infinity_h==1) then
             !    i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_dn_by_sym(product_table(get_inverse_dih(sym2),sym1))
             ! else
             !    i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_dn_by_sym(product_table(sym2,sym1))
             ! endif

             ! if(i_open.ne.0) then

             !    sym3 = orbital_symmetries(spatial_index1)
             !    electron_prob = 0 !Just a temporary variable
             !    do i = 1,num_occ_dn_by_sym(sym3)
             !       electron_prob = electron_prob + &
             !            &sqrt_integrals(occ_orb(excite_from_1_i),occ_orb_dn_by_sym(sym3,i))
             !    enddo

             !    denominator = sym_sum_cs_sqrt(sym3,occ_orb(excite_from_1_i)) - electron_prob

             !    temp_prob = temp_prob + sqrt_integrals(occ_orb(excite_from_2_i),excite_to_2)&
             !         &/cs_sqrt_prime_spin(excite_from_2_i)*&
             !         &sqrt_integrals(occ_orb(excite_from_1_i),excite_to_1)&
             !         &/denominator

             ! endif

             !temporary test of uniform probabilty
            ! Need to modify proposal probability for this move. For this, we calculate the probability of opposite excitation
             ! i.e. the excitation in which the electron at spatial_index2 was filled first
             if(uniform_sampling) then

                if(d_infinity_h==1) then
                   i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_dn_by_sym(product_table(get_inverse_dih(sym2),sym1))
                else
                   i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_dn_by_sym(product_table(sym2,sym1))
                endif


                if (i_open.ne.0) then
                   proposal_prob = proposal_prob*(temp1 + (1.0_rk/i_open))! correcting for 2 pathways to get this excita
                else
                   proposal_prob = proposal_prob*temp1! correcting for 2 pathways to get this excita
                endif
             else
                proposal_prob = proposal_prob*temp_prob
             endif
             if(debug) write(6,*)'prop',proposal_prob,sym1,spatial_index2,product_table(orbital_symmetries(spatial_index2),sym1),i_open
             if(debug) write(6,'(4B16)') det_i_up,det_i_dn,det_j_up,det_j_dn
          endif
       endif
    endif

   !call restrict_excitation(det_i_up,det_i_dn,too_excited)
   !if (too_excited) then
   !   return
   !else
       if (time_sym) then
           if ((det_j_up .eq. det_i_up) .and. (det_j_dn .eq. det_i_dn)) return  ! Already accounted for in diagonal
           if ((det_j_dn .eq. det_i_up) .and. (det_j_up .eq. det_i_dn)) return  ! Already accounted for in diagonal
           norm_j=1._rk
           if (det_j_up .eq. det_j_dn) then
              if (z .eq. 1) then
                 call is_connected_chem(det_i_up,det_i_dn,det_j_dn,det_j_up,are_they_connected,excite_level_sym,proposal_prob_sym)
                 norm_j = 1/proposal_prob
                 if(excite_level_sym.eq.1) then
                    proposal_prob = (proposal_prob + proposal_prob_sym*proposal_prob_single)/2
                 else
                    proposal_prob = (proposal_prob + proposal_prob_sym*proposal_prob_double)/2
                 endif
                 norm_j = 2*norm_j/sqrt2 * (proposal_prob)
                !norm_j=sqrt2 !Should be 1/sqrt(2) * (proposal_prob + proposal_prob_sym*proposal_prob_double/single) : switch off of excite_level_sym
              else
                norm_j=0._rk
                return
              endif
              call hamiltonian_chem(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level, matrix_element)               ! Proposed state is special, dont need to time reverse
              matrix_element=(norm_j/norm_i)*matrix_element
              !write (6,*) "Matrix element (off_diagonal move chem)",matrix_element
           else
               norm_j=1._rk
               call hamiltonian_chem(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level, matrix_element1)
! Definitely connected, because it was proposed
               if (spr) then
                    write (6,*) "Checking proposal probability of is_connected"
                    call is_connected_chem(det_i_up,det_i_dn,det_j_up,det_j_dn,are_they_connected,excite_level_sym,proposal_prob_sym)
                    write(6,*) "Proposal prob (Mihir)", proposal_prob
                    if (excite_level_sym .eq. 1) write(6,*) "Proposal prob (HJC) Exc 1",proposal_prob_sym*proposal_prob_single
                    if (excite_level_sym .eq. 2) write(6,*) "Proposal prob (HJC) Exc 2",proposal_prob_sym*proposal_prob_double
                    are_they_connected=.false.
                    excite_level_sym=-1
                    proposal_prob_sym=0._rk
                    call flush(6)
               endif
               if (abs(matrix_element1) .gt. 1.0e-10) nnzero=nnzero+1
               call is_connected_chem(det_i_up,det_i_dn,det_j_dn,det_j_up,are_they_connected,excite_level_sym,proposal_prob_sym)
               if (spr) then
                    write (6,*) "Are they connected?", det_i_up,det_i_dn,det_j_dn,det_j_up,are_they_connected
                    write (6,*) "Excite level sym (Method 1)",excite_level_sym
                    call excitation_level(det_i_up,det_i_dn,det_j_dn,det_j_up,excite_level_sym)
                    write (6,*) "Excite level sym (Method 2)",excite_level_sym
               endif
               if (are_they_connected) then
                  call hamiltonian_chem(det_i_up, det_i_dn, det_j_dn, det_j_up, excite_level_sym, matrix_element2)
                  !if (abs(matrix_element2) .gt. 1.0e-10) then
                  if (spr) write (6,*) "Proposal prob before",proposal_prob
                  nnzero=nnzero+1
                  if (excite_level_sym .eq. 1) proposal_prob=proposal_prob+(proposal_prob_sym*proposal_prob_single) ! Proposal probability of the second route
                  if (excite_level_sym .eq. 2) proposal_prob=proposal_prob+(proposal_prob_sym*proposal_prob_double) ! Proposal probability of the second route
                  matrix_element=(norm_j/norm_i)*(matrix_element1+z*matrix_element2)
                  if (spr) write (6,*) "Proposal prob after",proposal_prob
                  !else
                  !     matrix_element=(norm_j/norm_i)*(matrix_element1)
                  !     matrix_element2=0._rk
                  !endif
              else
                  matrix_element=(norm_j/norm_i)*(matrix_element1)
              endif
          endif
          if (det_j_up .gt. det_j_dn) then                 ! Convert to its representative
             tmp_det=det_j_up
             det_j_up=det_j_dn
             det_j_dn=tmp_det
             matrix_element=matrix_element*z
          endif
       else ! no time-reversal symmetry
          !write(6,'(B64)') det_i_up
          !write(6,'(B64)') det_i_dn
          !write(6,'(B64)') det_j_up
          !write(6,'(B64)') det_j_dn
          !write(6,*) det_i_up,det_i_dn,det_j_up,det_j_dn,excite_level,matrix_element
            call hamiltonian_chem(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level, matrix_element)

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
          if(debug) write (6,*) "Outgoing determinant",det_j_up,det_j_dn
       if (spr) then
          write (6,*) "Outgoing determinant",det_j_up,det_j_dn
         write (6,*) "Tau=",tau
         write (6,*) "Abs (matrix element) by Method 1 ",abs(matrix_element)
         if (time_sym) then
           write (6,*) "Confirming matrix element"
           call hamiltonian_chem_time_sym(det_i_up,det_i_dn,det_j_up,det_j_dn,matrix_element2)
           write (6,*) "Abs (matrix element) by Method 2",abs(matrix_element2)
           if (abs(matrix_element-matrix_element2) .gt. 1.0e-10) write (6,*) "Ham INCONSISTENT!"
         endif
         write (6,*) "Proposal probability",proposal_prob
       endif
!      weight_j = int((acceptance_prob+rannyu())*sign(1._rk,-matrix_element))
       weight_j = - tau *  matrix_element / proposal_prob
       if (spr) write (6,*) "weight_j",weight_j
       if (spr) call flush(6)
   !endif


       if(debug) then
          tmp_ctr = tmp_ctr+1
          if(tmp_ctr.gt.1000) stop
       endif
   if (spr) then
     tmp_ctr=tmp_ctr+1
     if (tmp_ctr .gt. 1000) stop
   endif

!  if(ipr.ge.1) then
!    icount=icount+1
!    if(mod(icount,1000).eq.0) then
!      if(excite_level.eq.1) then
!         write(11,'(''excite_level, det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j/tau, proposal_prob'',i2,4i20,f9.1,es9.2)') excite_level, det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j/tau, proposal_prob
!      endif
!      if(excite_level.eq.2) then
!         write(12,'(''excite_level, det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j/tau, proposal_prob'',i2,4i20,f9.1,es9.2)') excite_level, det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j/tau, proposal_prob
!      endif
!      write(13,'(''excite_level, det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j/tau, proposal_prob'',i2,4i20,f9.1,es9.2)') excite_level, det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j/tau, proposal_prob
!    endif
!    if(weight_j.gt.1000*tau) write(14,'(''excite_level, det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j/tau, proposal_prob'',i2,4i20,es9.1,es9.2,9f5.2)') excite_level, det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j/tau, proposal_prob, proposal_prob_single, proposal_prob_double
!  endif

 end subroutine off_diagonal_move_chem_cauchySchwarz


!===============================================================================================================
  subroutine off_diagonal_move_chem(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j, importance_sampling, iwalk_child)
    !----------------------------------------------------------------------------------------------------------
    ! Description : Perform an offdiagonal move in SQMC for a slater determinant starting at det_i
    !               which is split into det_i_up and det_i_dn. det_i has weight = +-1.
    !               since the looping over walkers on same determinant is done
    !               in general part of the code. If move fails, for any reason,
    !               weight_j is returned as 0. Otherwise, return new determinant,
    !               det_j, and its weight, H_jj (diagonal element of Hamiltonian).
    !
    ! Created     : Mihir Khadilkar, 1 Dec 2010
    ! last edited : A Holmes, 21 May 2012
    !-----------------------------------------------------------------------------------------------------------
    use more_tools, only : binary_search
    use common_psi_t, only : psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_g,psi_g_epsilon,psi_g_epsilon_inv
    implicit none

    !dummy variables
!   integer, intent(out) :: weight_j
    integer, intent(in)  :: iwalk_child
    real(rk), intent(out) :: weight_j
   !real(rk),optional, intent(out) :: diagonal_matrix_element
    integer,optional,intent(in) :: importance_sampling
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(out) :: det_j_up
    type(ik_vec), intent(out) :: det_j_dn
    type(ik_vec) :: tmp_det
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(out) :: det_j_up
    integer(ik), intent(out) :: det_j_dn
    integer(ik) :: tmp_det
#endif

    !local variables
    !excite_from_1_i and excite_from_2_i are both indices, not abolute orbital numbers
    integer :: excite_from_1_i, excite_from_2_i, excite_to_1,excite_to_2, excite_level, excite_level_sym, spatial_index1, spatial_index2
    integer :: n_single, n_double, n_total, sym1, sym2
    integer :: n_double_up, n_double_dn, n_double_both
    integer :: tot_spin_from !in units of 1/2
    integer :: i, i_open, j
    integer,save :: n_occ_up,n_occ_dn
    ! 100 -> number of symmetries, norb -> maximum of norb of a specific symmetry
    real(rk), parameter :: epsilon = 1.e-10_rk
!   real(rk) :: ran, rannyu, proposal_prob, acceptance_prob, matrix_element, temp1,temp2
    real(rk) :: proposal_prob, proposal_prob_sym, matrix_element,matrix_element1,matrix_element2, temp1
    real(rk) :: proposal_prob_single, proposal_prob_double
    real(rk) :: norm_i,norm_j
!   logical ::  too_excited
    integer ::  nnzero
    logical ::  are_they_connected
    logical ::  spr,debug,foundOrb
    logical :: imp
    integer :: location
    real(rk) :: guiding_wf_ratio

    uniform_sampling=.true.
    imp=.false.
    if (present(importance_sampling)) then
        if (importance_sampling .eq. 1) imp=.true.
    endif

    !occ_orb_up_by_sym_copy = 0
    !occ_orb_dn_by_sym_copy = 0
    if(iwalk_child.eq.1) then
       !Setup matrices / vectors
       !We only have to do this once per determinant; so, only the first
       !time it is called (iwalk_child=1) per determinant
       n_occ_up          = 0
       n_occ_dn          = 0
       occ_orb_up        = 0
       occ_orb_dn        = 0
       num_occ_up_by_sym = 0
       num_occ_dn_by_sym = 0

       tmp_det           = det_i_up

       do while (tmp_det > 0)
          i                                               = trailz(tmp_det)+1
          sym1                                            = orbital_symmetries(i)
          num_occ_up_by_sym(sym1)                    = num_occ_up_by_sym(sym1) + 1
          n_occ_up                                        = n_occ_up+1
          occ_orb_up_by_sym(sym1,num_occ_up_by_sym(sym1)) = i
          occ_orb_up(n_occ_up)                            = i
          tmp_det                                         = ibclr(tmp_det,i-1)
       enddo

       tmp_det = det_i_dn
       do while (tmp_det > 0)
          i                                               = trailz(tmp_det)+1
          sym1                                            = orbital_symmetries(i)
          num_occ_dn_by_sym(sym1)                 = num_occ_dn_by_sym(sym1) + 1
          occ_orb_dn_by_sym(sym1,num_occ_dn_by_sym(sym1)) = i
          n_occ_dn                                        = n_occ_dn+1
          occ_orb_dn(n_occ_dn)                            = i
          tmp_det                                         = ibclr(tmp_det,i-1)
       enddo

       !Store all occupied orbitals in one vector
       do i =1,n_occ_up
          occ_orb(i) = occ_orb_up(i)
       enddo

       do i=1,n_occ_dn
          occ_orb(i+n_occ_up) = occ_orb_dn(i)
       enddo
    endif

    debug=.true.
    debug=.false.
    !spr=.true.
    spr=.false.
    excite_level = 0
    excite_from_1_i = 0
    excite_from_2_i = 0
    excite_to_1 = 0
    excite_to_2 = 0
    proposal_prob = 1.0_rk
    nnzero=0
    !initialize det_j to det_i and weight_j to 0

    det_j_up = det_i_up
    det_j_dn = det_i_dn
    weight_j = 0
    if(spr) write (6,*) "Incoming determinant",det_i_up,det_i_dn
    if(debug) write (6,*) "Incoming determinant",det_i_up,det_i_dn

    if (time_sym) then
        if (det_i_up .eq. det_i_dn) then
             if (z .eq. 1) then
                norm_i=sqrt2
             else
                norm_i=0._rk
                write (6,*) "This state is inconsistent with the time symmetry you put in and should never enter this routine as an incoming state"
                stop
             endif
        else
               norm_i=1._rk
        endif
    endif

    n_single = (nup-n_core_orb)*(norb - nup) + (ndn-n_core_orb)*(norb - ndn)
    n_double_up = (nup-n_core_orb)*(nup -n_core_orb-1)*(norb-nup)*(norb-nup-1)/4
    n_double_dn = (ndn-n_core_orb)*(ndn -n_core_orb-1)*(norb-ndn)*(norb-ndn-1)/4
    n_double_both = (nup-n_core_orb)*(norb-nup)*(ndn-n_core_orb)*(norb-ndn)

    n_double = n_double_up + n_double_dn + n_double_both
    n_total = n_single + n_double

    proposal_prob_single=(n_single)/real(n_total,rk)
    proposal_prob_double=(n_double)/real(n_total,rk)        ! or 1 - proposal_prob_single

    if (random_int(n_total)> n_single) then
       excite_level = 2
       proposal_prob = n_double /real(n_total,rk)
       ! the excitation is now a double excitation
       ! choose 2 electrons to excite from at random

       ! Edited by AAH on 21 May 2012
       ! The number of valence electrons nval=2*n_core_orb
       ! If excite_from_1_i=nval then excite_from_2_i != excite_from_1_i and we have picked it with prob 1/(nval-1)
       ! If excite_from_1_i!=nval then excite_from_2_i = excite_from_1_i with prob 1/(nval-1) and we change excite_from_2_i to nval
       ! Since excite_from_2_i = excite_from_1_i with prob 1/(nval-1) was picked with prob 1/(nval-1), so is the new location, nval.
       excite_from_1_i = random_int(nelec-2*n_core_orb)
       excite_from_2_i = random_int(nelec-2*n_core_orb-1)
       if (excite_from_2_i==excite_from_1_i)  excite_from_2_i=nelec-2*n_core_orb
      !excite_from_1_i = random_int(nelec-2*n_core_orb)
      !do
      !   excite_from_2_i = random_int(nelec-2*n_core_orb)
      !   if (excite_from_1_i /= excite_from_2_i) exit
      !enddo
       ! End of edit

       tot_spin_from = 0
       if (excite_from_1_i > nup-n_core_orb) then
          !spin down
          tot_spin_from = tot_spin_from - 1
          excite_from_1_i = excite_from_1_i + 2*n_core_orb
       else
          !spin up
          tot_spin_from = tot_spin_from + 1
          excite_from_1_i = excite_from_1_i + n_core_orb
       endif

       if (excite_from_2_i > nup-n_core_orb) then
          !spin down
          tot_spin_from = tot_spin_from - 1
          excite_from_2_i = excite_from_2_i + 2*n_core_orb
       else
          !spin up
          tot_spin_from = tot_spin_from + 1
          excite_from_2_i = excite_from_2_i + n_core_orb
       endif

    else

       excite_level = 1
       proposal_prob = proposal_prob * n_single /(n_total*1.0_rk)
       ! excitation now is a single excitation
       ! choose electron to excite from at random
       excite_from_1_i = random_int(nelec-2*n_core_orb)
       ! calculate if the spin is up or down
       !tot_spin_from = 0

       if (excite_from_1_i > nup-n_core_orb) then
          !spin down
          tot_spin_from = -1
          excite_from_1_i = excite_from_1_i + 2*n_core_orb
       else
          !spin up
          tot_spin_from = 1
          excite_from_1_i = excite_from_1_i + n_core_orb
       endif

    endif

   ! remove excited electrons in det_j
    ! set excite_from_2_i > excite_from_1_i
    i =  excite_from_1_i + excite_from_2_i
    excite_from_1_i = i - max(excite_from_1_i,excite_from_2_i)
    excite_from_2_i = i - excite_from_1_i

    sym1 = 1 ! this variable will store total symmetry of the excitation.
    ! for the single case, just the symmetry of excite_from_1_i
    ! for the double case, product of symmetries of excite_from_1_i and excite_from_2_i
    !calculate total spin of the excitation
    !pick single or double excitation

    ! Edited by MJO on 16 Jan 2015. Sped up by using new data structures.
    if(excite_level==1) then
       !Because the excite_from* were reorded, and are initially set to 0,
       !if we have a single excitation, then excite_from_1_i is 0! So,
       !use excite_from_2_i, instead
       i = occ_orb(excite_from_2_i)

       if(excite_from_2_i.le.n_occ_up) then
          det_j_up=ibclr(det_j_up, i-1)
       else
          det_j_dn=ibclr(det_j_dn, i-1)
       endif
       sym1 = orbital_symmetries(i)

    else !excite_level==2
       i = occ_orb(excite_from_1_i)

       if(excite_from_1_i.le.n_occ_up) then
          det_j_up=ibclr(det_j_up, i-1)
       else
          det_j_dn=ibclr(det_j_dn, i-1)
       endif

       sym1 = orbital_symmetries(i)

       i = occ_orb(excite_from_2_i)
       if(excite_from_2_i.le.n_occ_up) then
          det_j_up=ibclr(det_j_up, i-1)
       else
          det_j_dn=ibclr(det_j_dn, i-1)
       endif

       sym1 = product_table(sym1,orbital_symmetries(i))
    endif

    !write(6,*)'sym1',sym1
! At this point, we have selected a single or a double excitation
! and removed appropriate electrons from det_j_up and/or det_j_dn
!----------------------------------------------------------------------------
    i_open = 0
    if (excite_level==1) then
       proposal_prob = proposal_prob /  (nelec-2*n_core_orb)
       if (tot_spin_from ==1) then
          if(debug.eqv..true.) write(6,*)'single,totspin1'
          if(spr.eqv..true.) write(6,*)'single,totspin1'

          i_open = num_orb_by_sym(sym1) - num_occ_up_by_sym(sym1)
          ! now select one of them for excite_to_1
          if (i_open.ne.0) then
             excite_to_1 = random_int(i_open)
             proposal_prob = proposal_prob / i_open
          else
             return  ! no valid move possible, return with zero weight_j
          endif

          do i=1,num_occ_up_by_sym(sym1)
             if(occ_orb_up_by_sym(sym1,i).le.which_orb_by_sym(sym1,excite_to_1)) then
                excite_to_1=excite_to_1+1
             else
                exit
             endif
          enddo

          i = which_orb_by_sym(sym1,excite_to_1)
          !add electron to this orbital
          det_j_up = ibset(det_j_up,i-1)

       else ! total spin -1
          if(debug.eqv..true.) write(6,*)'single,totspin-1'
          if(spr.eqv..true.) write(6,*)'single,totspin-1'

          i_open = num_orb_by_sym(sym1) - num_occ_dn_by_sym(sym1)
          ! now select one of them for excite_to_1
          if (i_open.ne.0) then
             excite_to_1 = random_int(i_open)
             proposal_prob = proposal_prob / i_open
          else
             return  ! no valid move possible, return with zero weight_j
          endif

          do i=1,num_occ_dn_by_sym(sym1)
             if(occ_orb_dn_by_sym(sym1,i).le.which_orb_by_sym(sym1,excite_to_1)) then
                excite_to_1=excite_to_1+1
             else
                exit
             endif
          enddo

          i = which_orb_by_sym(sym1,excite_to_1)
          !add electron to this orbital
          det_j_dn = ibset(det_j_dn,i-1)
       endif

       ! end of single excitation part
    else
       ! now filling up double excitation: excite_to_1 and excite_to_2

       proposal_prob = proposal_prob*2.0 /(1.0_rk*(nelec-2*n_core_orb)*(nelec-2*n_core_orb-1))
       if(debug.eqv..true.)print*,'before excite',excite_from_1_i,excite_from_2_i,proposal_prob
       !find a single spin orbital to excite to
       if (tot_spin_from == 2) then
          !have to excite to spinup orbital
          if(debug.eqv..true.) write(6,*)'double,totspin2'
          if(spr.eqv..true.) write(6,*)'double,totspin2'
          excite_to_1 = random_int(norb - nup)
          proposal_prob = proposal_prob / (norb - nup)

          ! immediately fill excite_to_1
          do i=1,n_occ_up
             if(occ_orb_up(i).le.excite_to_1) then
                excite_to_1=excite_to_1+1
             else
                exit
             endif
          enddo

          i = excite_to_1
          spatial_index1 = i
          det_j_up = ibset(det_j_up, i-1)

          !When d_inf_h is used, the irreducible representation is no longer its own inverse
          if(d_infinity_h==1) then
             sym2 = product_table(get_inverse_dih(orbital_symmetries(spatial_index1)),sym1)
          else
             sym2 = product_table(orbital_symmetries(spatial_index1),sym1)
          endif

          i_open = num_orb_by_sym(sym2) - num_occ_up_by_sym(sym2)
          ! If the two new orbitals have the same symmetry, we have to subtract by one, since we
          ! just filled one of the orbitals
          if(sym2.eq.orbital_symmetries(spatial_index1)) i_open = i_open - 1

          if (i_open.ne.0) then
             ! now among these, choose the 2nd orbital to excite to
             excite_to_2 = random_int(i_open)
             temp1 = 1.0_rk / i_open ! to be later used for proposal probability
          else
             return ! no valid move possible, return with zero weight_j
          endif

          if(sym2.eq.orbital_symmetries(spatial_index1)) then
             do i=1,num_occ_up_by_sym(sym2)
                occ_orb_up_by_sym_copy(sym2,i) = occ_orb_up_by_sym(sym2,i)
             enddo

             foundOrb=.false.
             do i=1,num_occ_up_by_sym(sym2)
                if(occ_orb_up_by_sym_copy(sym2,i).ge.spatial_index1) then
                   do j=num_occ_up_by_sym(sym2),i,-1
                      occ_orb_up_by_sym_copy(sym2,j+1) = occ_orb_up_by_sym_copy(sym2,j)
                   enddo
                   occ_orb_up_by_sym_copy(sym2,i)   = spatial_index1
                   foundOrb=.true.
                   exit
                endif
             enddo
             if(.not.foundOrb) occ_orb_up_by_sym_copy(sym2,num_occ_up_by_sym(sym2)+1) = spatial_index1

!             num_occ_up_by_sym(sym2) = num_occ_up_by_sym(sym2) + 1
             do i=1,num_occ_up_by_sym(sym2) + 1
                if(occ_orb_up_by_sym_copy(sym2,i).le.which_orb_by_sym(sym2,excite_to_2)) then
                   excite_to_2=excite_to_2+1
                else
                   exit
                endif
             enddo
          else
             do i=1,num_occ_up_by_sym(sym2)
                if(occ_orb_up_by_sym(sym2,i).le.which_orb_by_sym(sym2,excite_to_2)) then
                   excite_to_2=excite_to_2+1
                else
                   exit
                endif
             enddo
          endif

          ! !Fix because we already filled one of the orbitals
          ! if(sym2.eq.orbital_symmetries(spatial_index1)) then
          !    foundOrb=.false.
          !    do i=1,num_occ_up_by_sym(sym2)
          !       if(occ_orb_up_by_sym(sym2,i).ge.spatial_index1) then
          !          do j=num_occ_up_by_sym(sym2),i,-1
          !             occ_orb_up_by_sym(sym2,j+1) = occ_orb_up_by_sym(sym2,j)
          !          enddo
          !          occ_orb_up_by_sym(sym2,i)   = spatial_index1
          !          foundOrb=.true.
          !          exit
          !       endif
          !    enddo
          !    if(.not.foundOrb) occ_orb_up_by_sym(sym2,num_occ_up_by_sym(sym2)+1) = spatial_index1
          !    num_occ_up_by_sym(sym2) = num_occ_up_by_sym(sym2) + 1
          ! endif

          ! do i=1,num_occ_up_by_sym(sym2)
          !    if(occ_orb_up_by_sym(sym2,i).le.which_orb_by_sym(sym2,excite_to_2)) then
          !       excite_to_2=excite_to_2+1
          !    else
          !       exit
          !    endif
          ! enddo

          i = which_orb_by_sym(sym2,excite_to_2)
          spatial_index2 = i
          det_j_up = ibset(det_j_up, i-1)


          ! Need to modify proposal probability for this move. For this, we calculate the probability of opposite excitation
          ! i.e. the excitation in which the electron at spatial_index2 was filled first

          if(d_infinity_h==1) then
             i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_up_by_sym(product_table(get_inverse_dih(sym2),sym1))
          else
             i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_up_by_sym(product_table(sym2,sym1))
          endif

          if(sym2.eq.orbital_symmetries(spatial_index1)) i_open = i_open - 1
          !write(6,*)'symmstuff',sym1,spatial_index1,spatial_index2,orbital_symmetries(spatial_index1),orbital_symmetries(spatial_index2),product_table(orbital_symmetries(spatial_index2),sym1)
          if (i_open==0) then
             ! not possible to excite otherwise
             proposal_prob = proposal_prob*temp1 ! correcting for 2 pathways to get this excitation
          else
             proposal_prob = proposal_prob*(temp1+(1.0_rk/i_open)) ! correcting for 2 pathways to get this excitation
          endif
          if(debug) write(6,*)'prop',proposal_prob,sym1,spatial_index2,spatial_index1,product_table(orbital_symmetries(spatial_index2),sym1),i_open
          if(debug) write(6,'(4B64)')det_i_up,det_i_dn,det_j_up,det_j_dn

          !write(6,*)'prop',proposal_prob
       elseif (tot_spin_from == -2) then
          if(debug.eqv..true.) write(6,*)'double,totspin-2'
          if(spr.eqv..true.) write(6,*)'double,totspin-2'
          excite_to_1 = random_int(norb - ndn)
          proposal_prob = proposal_prob /(norb - ndn)

          ! immediately fill excite_to_1
          do i=1,n_occ_dn
             if(occ_orb_dn(i).le.excite_to_1) then
                excite_to_1=excite_to_1+1
             else
                exit
             endif
          enddo

          i = excite_to_1
          spatial_index1 = i
          det_j_dn = ibset(det_j_dn, i-1)

          !When d_inf_h is used, the irreducible representation is no longer its own inverse
          if(d_infinity_h==1) then
             sym2 = product_table(sym1,get_inverse_dih(orbital_symmetries(spatial_index1)))
          else
             sym2 = product_table(orbital_symmetries(spatial_index1),sym1)
          endif


          i_open = num_orb_by_sym(sym2) - num_occ_dn_by_sym(sym2)
          ! If the two new orbitals have the same symmetry, we have to subtract by one, since we
          ! just filled one of the orbital
          if(sym2.eq.orbital_symmetries(spatial_index1)) i_open = i_open - 1

          if (i_open.ne.0) then
             ! now among these, choose the 2nd orbital to excite to
             excite_to_2 = random_int(i_open)
             temp1 = 1.0_rk / i_open ! to be later used for proposal probability
          else
             return ! no valid move possible, return with zero weight_j
          endif

          !Fix for same symmetry
          if(sym2.eq.orbital_symmetries(spatial_index1)) then
             do i=1,num_occ_dn_by_sym(sym2)
                occ_orb_dn_by_sym_copy(sym2,i) = occ_orb_dn_by_sym(sym2,i)
             enddo

             foundOrb=.false.
             do i=1,num_occ_dn_by_sym(sym2)
                if(occ_orb_dn_by_sym_copy(sym2,i).ge.spatial_index1) then
                   do j=num_occ_dn_by_sym(sym2),i,-1
                      occ_orb_dn_by_sym_copy(sym2,j+1) = occ_orb_dn_by_sym_copy(sym2,j)
                   enddo
                   occ_orb_dn_by_sym_copy(sym2,i)   = spatial_index1
                   foundOrb=.true.
                   exit
                endif
             enddo
             if(.not.foundOrb) occ_orb_dn_by_sym_copy(sym2,num_occ_dn_by_sym(sym2)+1) = spatial_index1

!             num_occ_dn_by_sym(sym2) = num_occ_dn_by_sym(sym2) + 1
             do i=1,num_occ_dn_by_sym(sym2) + 1
                if(occ_orb_dn_by_sym_copy(sym2,i).le.which_orb_by_sym(sym2,excite_to_2)) then
                   excite_to_2=excite_to_2+1
                else
                   exit
                endif
             enddo
          else
             do i=1,num_occ_dn_by_sym(sym2)
                if(occ_orb_dn_by_sym(sym2,i).le.which_orb_by_sym(sym2,excite_to_2)) then
                   excite_to_2=excite_to_2+1
                else
                   exit
                endif
             enddo
          endif




          !fix because we already filled one of the orbitals
          ! if(sym2.eq.orbital_symmetries(spatial_index1)) then
          !    foundOrb=.false.

          !    do i=1,num_occ_dn_by_sym(sym2)
          !       if(occ_orb_dn_by_sym(sym2,i).ge.spatial_index1) then
          !          do j=num_occ_dn_by_sym(sym2),i,-1
          !             occ_orb_dn_by_sym(sym2,j+1) = occ_orb_dn_by_sym(sym2,j)
          !          enddo
          !          occ_orb_dn_by_sym(sym2,i)   = spatial_index1
          !          foundOrb=.true.
          !          exit
          !       endif
          !    enddo
          !    if(.not.foundOrb) occ_orb_dn_by_sym(sym2,num_occ_dn_by_sym(sym2)+1) = spatial_index1
          !    num_occ_dn_by_sym(sym2) = num_occ_dn_by_sym(sym2) + 1
          ! endif

          ! do i=1,num_occ_dn_by_sym(sym2)
          !    if(occ_orb_dn_by_sym(sym2,i).le.which_orb_by_sym(sym2,excite_to_2)) then
          !       excite_to_2=excite_to_2+1
          !    else
          !       exit
          !    endif
          ! enddo

          i = which_orb_by_sym(sym2,excite_to_2)
          spatial_index2 = i
          det_j_dn = ibset(det_j_dn, i-1)


          ! Need to modify proposal probability for this move. For this, we calculate the probability of opposite excitation
          ! i.e. the excitation in which the electron at spatial_index2 was filled first
          ! It seems to take longer if I do this assignment
          !sym2  = product_table(orbital_symmetries(spatial_index2),sym1)
          if(d_infinity_h==1) then
             i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_dn_by_sym(product_table(get_inverse_dih(sym2),sym1))
          else
             i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_dn_by_sym(product_table(sym2,sym1))
          endif


          if(sym2.eq.orbital_symmetries(spatial_index1)) i_open = i_open - 1

          if (i_open==0) then
             ! not possible to excite otherwise
             proposal_prob = proposal_prob*temp1 ! correcting for 2 pathways to get this excitation
          else
             proposal_prob = proposal_prob*(temp1+(1.0_rk/i_open)) ! correcting for 2 pathways to get this excitation
          endif
          if(debug) write(6,*)'prop',proposal_prob,sym1,spatial_index2,spatial_index1,product_table(orbital_symmetries(spatial_index2),sym1),i_open
          if(debug) write(6,'(4B64)')det_i_up,det_i_dn,det_j_up,det_j_dn
       else
          !up-down case of a double excitation
          proposal_prob = proposal_prob*  1.0_rk/(2*norb - ndn-nup)
          sym2 =1
          !can excite to any open spin orbital
          excite_to_1 = random_int(2*norb - nup - ndn)

          if (excite_to_1 <= (norb - nup)) then
             if(debug.eqv..true.) write(6,*)'double,totspin1'
             if(spr.eqv..true.) write(6,*)'double,totspin1'

             !excite to spin-up orbital
             do i=1,n_occ_up
                if(occ_orb_up(i).le.excite_to_1) then
                   excite_to_1=excite_to_1+1
                else
                   exit
                endif
             enddo
             i = excite_to_1
             det_j_up = ibset(det_j_up, i-1)
             spatial_index1 = i
             sym2 = orbital_symmetries(i)


             !try to find a spin_down orbital that conserves symmetry
             !When d_inf_h is used, the irreducible representation is no longer its own inverse
             if(d_infinity_h==1) then
                sym2 = product_table(sym1,get_inverse_dih(sym2))
             else
                sym2 = product_table(sym1,sym2)
             endif

             i_open = num_orb_by_sym(sym2) - num_occ_dn_by_sym(sym2)
             if (i_open.ne.0) then
                ! total possibilities counted
                excite_to_2 = random_int(i_open)
                temp1 = 1.0_rk /  i_open
             else
                return ! no valid move possible, return with zero weight_j
             endif

             do i=1,num_occ_dn_by_sym(sym2)
                if(occ_orb_dn_by_sym(sym2,i).le.which_orb_by_sym(sym2,excite_to_2)) then
                   excite_to_2=excite_to_2+1
                else
                   exit
                endif
             enddo

             i = which_orb_by_sym(sym2,excite_to_2)
             spatial_index2 = i
             det_j_dn = ibset(det_j_dn, i-1)

            ! Need to modify proposal probability for this move. For this, we calculate the probability of opposite excitation
            ! i.e. the excitation in which the electron at spatial_index2 was filled first
            ! sym2   = product_table(orbital_symmetries(spatial_index2),sym1)
             if(d_infinity_h==1) then
                i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_up_by_sym(product_table(get_inverse_dih(sym2),sym1))
             else
                i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_up_by_sym(product_table(sym2,sym1))
             endif

            if (i_open.ne.0) then
               proposal_prob = proposal_prob*(temp1 + (1.0_rk/i_open))! correcting for 2 pathways to get this excita
            else
               proposal_prob = proposal_prob*temp1! correcting for 2 pathways to get this excita
            endif
             if(debug.eqv..true.) write(6,*)'prop',proposal_prob,sym1,spatial_index2,spatial_index1,product_table(orbital_symmetries(spatial_index2),sym1),i_open
             if(debug) write(6,'(4B64)')det_i_up,det_i_dn,det_j_up,det_j_dn
          else
             if(debug.eqv..true.) write(6,*)'double,totspin-1'
             if(spr.eqv..true.) write(6,*)'double,totspin-1'

             !excite to spin-down orbital

             do i=1,n_occ_dn
                if(occ_orb_dn(i).le.excite_to_1-(norb-nup)) then
                   excite_to_1=excite_to_1+1
                else
                   exit
                endif
             enddo
             i = excite_to_1-(norb-nup)
             det_j_dn = ibset(det_j_dn, i-1)
             spatial_index1 = i
             sym2 = orbital_symmetries(i)


             !When d_inf_h is used, the irreducible representation is no longer its own inverse
             if(d_infinity_h==1) then
                sym2 = product_table(sym1,get_inverse_dih(sym2))
             else
                sym2 = product_table(sym1,sym2)
             endif

             i_open = num_orb_by_sym(sym2) - num_occ_up_by_sym(sym2)
             if (i_open.ne.0) then
                ! total possibilities counted
                excite_to_2 = random_int(i_open)
                temp1 = 1.0_rk /  i_open
             else
                return ! no valid move possible, return with zero weight_j
             endif

             do i=1,num_occ_up_by_sym(sym2)
                if(occ_orb_up_by_sym(sym2,i).le.which_orb_by_sym(sym2,excite_to_2)) then
                   excite_to_2=excite_to_2+1
                else
                   exit
                endif
             enddo

             i = which_orb_by_sym(sym2,excite_to_2)
             spatial_index2 = i
             det_j_up = ibset(det_j_up, i-1)

            ! Need to modify proposal probability for this move. For this, we calculate the probability of opposite excitation
            ! i.e. the excitation in which the electron at spatial_index2 was filled first

             if(d_infinity_h==1) then
                i_open = num_orb_by_sym(product_table(get_inverse_dih(sym2),sym1)) - num_occ_dn_by_sym(product_table(get_inverse_dih(sym2),sym1))
             else
                i_open = num_orb_by_sym(product_table(sym2,sym1)) - num_occ_dn_by_sym(product_table(sym2,sym1))
             endif

            if (i_open.ne.0) then
               proposal_prob = proposal_prob*(temp1 + (1.0_rk/i_open))! correcting for 2 pathways to get this excita
            else
               proposal_prob = proposal_prob*temp1! correcting for 2 pathways to get this excita
            endif
             if(debug) write(6,*)'prop',proposal_prob,sym1,spatial_index2,product_table(orbital_symmetries(spatial_index2),sym1),i_open
             if(debug) write(6,'(4B16)') det_i_up,det_i_dn,det_j_up,det_j_dn

          endif
       endif
    endif

   !call restrict_excitation(det_i_up,det_i_dn,too_excited)
   !if (too_excited) then
   !   return
   !else
       if (time_sym) then
           if ((det_j_up .eq. det_i_up) .and. (det_j_dn .eq. det_i_dn)) return  ! Already accounted for in diagonal
           if ((det_j_dn .eq. det_i_up) .and. (det_j_up .eq. det_i_dn)) return  ! Already accounted for in diagonal
           norm_j=1._rk
           if (det_j_up .eq. det_j_dn) then
              if (z .eq. 1) then
                norm_j=sqrt2
              else
                norm_j=0._rk
                return
              endif
              call hamiltonian_chem(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level, matrix_element)               ! Proposed state is special, dont need to time reverse
              matrix_element=(norm_j/norm_i)*matrix_element
              !write (6,*) "Matrix element (off_diagonal move chem)",matrix_element
           else
               norm_j=1._rk
               call hamiltonian_chem(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level, matrix_element1)
! Definitely connected, because it was proposed
               if (spr) then
                    write (6,*) "Checking proposal probability of is_connected"
                    call is_connected_chem(det_i_up,det_i_dn,det_j_up,det_j_dn,are_they_connected,excite_level_sym,proposal_prob_sym)
                    write(6,*) "Proposal prob (Mihir)", proposal_prob
                    if (excite_level_sym .eq. 1) write(6,*) "Proposal prob (HJC) Exc 1",proposal_prob_sym*proposal_prob_single
                    if (excite_level_sym .eq. 2) write(6,*) "Proposal prob (HJC) Exc 2",proposal_prob_sym*proposal_prob_double
                    are_they_connected=.false.
                    excite_level_sym=-1
                    proposal_prob_sym=0._rk
                    call flush(6)
               endif
               if (abs(matrix_element1) .gt. 1.0e-10) nnzero=nnzero+1
               call is_connected_chem(det_i_up,det_i_dn,det_j_dn,det_j_up,are_they_connected,excite_level_sym,proposal_prob_sym)
               if (spr) then
                    write (6,*) "Are they connected?", det_i_up,det_i_dn,det_j_dn,det_j_up,are_they_connected
                    write (6,*) "Excite level sym (Method 1)",excite_level_sym
                    call excitation_level(det_i_up,det_i_dn,det_j_dn,det_j_up,excite_level_sym)
                    write (6,*) "Excite level sym (Method 2)",excite_level_sym
               endif
               if (are_they_connected) then
                  call hamiltonian_chem(det_i_up, det_i_dn, det_j_dn, det_j_up, excite_level_sym, matrix_element2)
                  !if (abs(matrix_element2) .gt. 1.0e-10) then
                  if (spr) write (6,*) "Proposal prob before",proposal_prob
                  nnzero=nnzero+1
                  if (excite_level_sym .eq. 1) proposal_prob=proposal_prob+(proposal_prob_sym*proposal_prob_single) ! Proposal probability of the second route
                  if (excite_level_sym .eq. 2) proposal_prob=proposal_prob+(proposal_prob_sym*proposal_prob_double) ! Proposal probability of the second route
                  matrix_element=(norm_j/norm_i)*(matrix_element1+z*matrix_element2)
                  if (spr) write (6,*) "Proposal prob after",proposal_prob
                  !else
                  !     matrix_element=(norm_j/norm_i)*(matrix_element1)
                  !     matrix_element2=0._rk
                  !endif
              else
                  matrix_element=(norm_j/norm_i)*(matrix_element1)
              endif
          endif
          if (det_j_up .gt. det_j_dn) then                 ! Convert to its representative
                   tmp_det=det_j_up
                   det_j_up=det_j_dn
                   det_j_dn=tmp_det
                   matrix_element=matrix_element*z
          endif
       else ! no time-reversal symmetry
          !write(6,'(B64)') det_i_up
          !write(6,'(B64)') det_i_dn
          !write(6,'(B64)') det_j_up
          !write(6,'(B64)') det_j_dn
          !write(6,*) det_i_up,det_i_dn,det_j_up,det_j_dn,excite_level,matrix_element
            call hamiltonian_chem(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level, matrix_element)

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
          if(debug) write (6,*) "Outgoing determinant",det_j_up,det_j_dn
       if (spr) then
          write (6,*) "Outgoing determinant",det_j_up,det_j_dn
         write (6,*) "Tau=",tau
         write (6,*) "Abs (matrix element) by Method 1 ",abs(matrix_element)
         if (time_sym) then
           write (6,*) "Confirming matrix element"
           call hamiltonian_chem_time_sym(det_i_up,det_i_dn,det_j_up,det_j_dn,matrix_element2)
           write (6,*) "Abs (matrix element) by Method 2",abs(matrix_element2)
           if (abs(matrix_element-matrix_element2) .gt. 1.0e-10) write (6,*) "Ham INCONSISTENT!"
         endif
         write (6,*) "Proposal probability",proposal_prob
       endif
!      weight_j = int((acceptance_prob+rannyu())*sign(1._rk,-matrix_element))
       weight_j = - tau *  matrix_element / proposal_prob
       if (spr) write (6,*) "weight_j",weight_j
       if (spr) call flush(6)
   !endif
       if(debug) then
          tmp_ctr = tmp_ctr+1
          if(tmp_ctr.gt.100000) stop
       endif
   if (spr) then
     tmp_ctr=tmp_ctr+1
     if (tmp_ctr .gt. 1000) stop
   endif

  end subroutine off_diagonal_move_chem
!=============================================================

!===============================================================================================================
  subroutine off_diagonal_move_chem_efficient_heatbath(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j, iwalk_child, importance_sampling, n_new_dets, excite_level)
    !----------------------------------------------------------------------------------------------------------
    ! Description : Perform an offdiagonal move in SQMC for a slater determinant starting at det_i
    !               which is split into det_i_up and det_i_dn. det_i has weight = +-1.
    !               since the looping over walkers on same determinant is done
    !               in general part of the code. If move fails, for any reason,
    !               weight_j is returned as 0. Otherwise, return new determinant,
    !               det_j, and its weight, H_jj (diagonal element of Hamiltonian).
    !
    ! Created     : A Holmes, 18 Mar 2015
    !-----------------------------------------------------------------------------------------------------------
    use more_tools, only : sample_discrete_distribution,binary_search
    use more_tools, only : sample_alias,setup_alias
    use mpi_routines, only : mpi_stop
    implicit none

    !dummy variables
    real(rk), intent(out) :: weight_j(2)
    integer,intent(in) :: iwalk_child
    integer,optional,intent(in) :: importance_sampling
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(out) :: det_j_up(2)
    type(ik_vec), intent(out) :: det_j_dn(2)
    type(ik_vec) :: tmp_det
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(out) :: det_j_up(2)
    integer(ik), intent(out) :: det_j_dn(2)
    integer(ik) :: tmp_det
#endif
    integer, intent(out) :: n_new_dets ! = 1 or 2. We have the option to return both a single and double excitation if single matrix element is large
    integer,intent(out) :: excite_level(2)

    !local variables
    integer :: excite_from_1, excite_from_2, excite_to_1,excite_to_2
    integer :: i, i_elec
    real(rk), parameter :: epsilon = 1.e-10_rk
    real(rk) :: proposal_prob, matrix_element
    real(rk) :: rannyu,sing_num,sing_den
    integer :: excite_spin
    real(rk) :: e1_prob_sav
    real(rk) :: c_e2_prob_sav
    real(rk) :: normalization
    logical :: same_spin
    real(rk) :: prob_1_then_2

    excite_from_1 = 0
    excite_from_2 = 0
    excite_to_1 = 0
    excite_to_2 = 0

    det_j_up = det_i_up
    det_j_dn = det_i_dn
    weight_j = 0

    excite_level(:) = -1

    if (iwalk_child==1) then
      tmp_det = det_i_up
      i_elec = 0
      do while (tmp_det .ne. 0)
        i = trailz(tmp_det)+1
        i_elec = i_elec + 1
        occ_up(i_elec) = i
        tmp_det = ibclr(tmp_det,i-1)
      enddo

      tmp_det = det_i_dn
      i_elec = 0
      do while (tmp_det .ne. 0)
        i = trailz(tmp_det)+1
        i_elec = i_elec + 1
        occ_dn(i_elec) = i
        tmp_det = ibclr(tmp_det,i-1)
      enddo

      ! set up first electron probabilities

      do i=1,nup
        elecs(i) = occ_up(i)
        e1_prob(i) = one_orbital_probabilities(elecs(i))
        c_e1_prob(i) = one_orbital_probabilities(elecs(i))
        if (i>1)  c_e1_prob(i) = c_e1_prob(i) + c_e1_prob(i-1)
      enddo

      do i=1,ndn
        elecs(i+nup) = occ_dn(i)+norb
        e1_prob(i+nup) = one_orbital_probabilities(occ_dn(i))
        c_e1_prob(i+nup) = one_orbital_probabilities(occ_dn(i))
        c_e1_prob(i+nup) = c_e1_prob(i+nup) + c_e1_prob(i+nup-1)
      enddo

      e1_prob(:) = e1_prob(:)/c_e1_prob(nup+ndn)
      c_e1_prob(:) = c_e1_prob(:)/c_e1_prob(nup+ndn)

    endif ! iwalk_child==1


    ! Choose electron 1 first, then choose electron2

    i = sample_discrete_distribution(c_e1_prob)
    excite_from_1 = elecs(i)
    proposal_prob = e1_prob(i)
    e1_prob_sav = e1_prob(i)

    ! Now for electron 2

    do i=1,nup
      if (elecs(i)==excite_from_1) then
        e2_prob(i) = 0._rk
        c_e2_prob(i) = 0._rk
      else
        e2_prob(i) = two_orbital_probabilities(excite_from_1,occ_up(i))
        c_e2_prob(i) = two_orbital_probabilities(excite_from_1,occ_up(i))
      endif
      if (i>1)  c_e2_prob(i) = c_e2_prob(i) + c_e2_prob(i-1)
    enddo

    do i=1,ndn
      if (elecs(i+nup)==excite_from_1) then
        e2_prob(i+nup) = 0._rk
        c_e2_prob(i+nup) = 0._rk
      else
        e2_prob(i+nup) = two_orbital_probabilities(excite_from_1,occ_dn(i)+norb)
        c_e2_prob(i+nup) = two_orbital_probabilities(excite_from_1,occ_dn(i)+norb)
      endif
      c_e2_prob(i+nup) = c_e2_prob(i+nup) + c_e2_prob(i+nup-1)
    enddo

    c_e2_prob_sav = c_e2_prob(nelec)
    e2_prob(:) = e2_prob(:)/c_e2_prob(nelec)
    c_e2_prob(:) = c_e2_prob(:)/c_e2_prob(nelec)

    i = sample_discrete_distribution(c_e2_prob)
    excite_from_2 = elecs(i)
    proposal_prob = proposal_prob*e2_prob(i)


    ! Now we have the two orbitals to excite from, excite_from_1,2
    ! Next, get first hole

    call choose_first_hole(excite_from_1,excite_from_2,excite_to_1,excite_spin)

    ! If first "hole" is actually occupied, no move is proposed
    if (excite_to_1<=norb) then
      if (btest(det_i_up,excite_to_1-1)) then
        weight_j = 0._rk
        return
      endif
    else
      if (btest(det_i_dn,excite_to_1-norb-1)) then
        weight_j = 0._rk
        return
      endif
    endif

    ! Now we have selected two electrons and the first hole
    ! Finally, select the second hole

    ! Choose hole 2 from the PDF given by four_orbital_probabilities(excite_from_1,excite_from_2,excite_from_1,:)
    ! Automatically takes into account electron spin and orbital symmetry, but does not take into account whether hole 2 is already occupied by current det_i

    ! Incorporate single excitations:
    ! Assume that excite_from_1 and excite_to_1 are of the same spin!

    if (excite_to_1<=norb) then
      call hamiltonian_chem(det_i_up,det_i_dn,ibset(ibclr(det_i_up,excite_from_1-1),excite_to_1-1),det_i_dn,1,matrix_element)
    else
      call hamiltonian_chem(det_i_up,det_i_dn,det_i_up,ibset(ibclr(det_i_dn,excite_from_1-norb-1),excite_to_1-norb-1),1,matrix_element)
    endif

    same_spin = (excite_spin.ne.0)

    sing_num = abs(matrix_element)

    sing_den = sing_num + Htot(excite_from_1,excite_from_2,excite_to_1) ! Single excitation probability denominator

    if (sing_num>(sing_den-sing_num)) then ! Single excitation matrix element greater than sum of double excitation matrix elements, so make both a single and double excitation.
      prob_1_then_2 = proposal_prob
      n_new_dets = 2
      normalization = e1_prob_sav / c_e2_prob_sav
      call generate_heatbath_single(det_i_up,det_i_dn,det_j_up(1),det_j_dn(1),excite_from_1,excite_to_1,matrix_element,proposal_prob,normalization) ! Does not change matrix_element!
      if (time_sym)  call apply_time_reversal_symmetry(det_i_up,det_i_dn,det_j_up(1),det_j_dn(1),occ_up,occ_dn,matrix_element,proposal_prob)
      weight_j(1) = -tau * matrix_element / proposal_prob
      excite_level(1) = 1

      ! The following should only be a problem if single move e.v. is added in
      if (Htot(excite_from_1,excite_from_2,excite_to_1)==0._rk) then
        weight_j(2) = 0._rk
        n_new_dets = 1
        return
      endif

      excite_to_2 = choose_second_hole(excite_from_1,excite_from_2,excite_to_1,same_spin)
      ! Remove moves to already-occupied orbitals
      if (excite_to_2<=norb) then
        if (btest(det_i_up,excite_to_2-1)) then
          n_new_dets = 1
          weight_j(2) = 0._rk
          return
        endif
      else
        if (btest(det_i_dn,excite_to_2-norb-1)) then
          n_new_dets = 1
          weight_j(2) = 0._rk
          return
        endif
      endif

      call generate_heatbath_double(det_i_up,det_i_dn,det_j_up(2),det_j_dn(2),excite_from_1,excite_from_2,excite_to_1,excite_to_2,occ_up,occ_dn,matrix_element,prob_1_then_2,proposal_prob)
      if (time_sym)  call apply_time_reversal_symmetry(det_i_up,det_i_dn,det_j_up(2),det_j_dn(2),occ_up,occ_dn,matrix_element,proposal_prob)
      weight_j(2) = -tau * matrix_element / proposal_prob
      excite_level(2) = 2
      return

    else ! Choose only one excitation (either single or double)
      n_new_dets = 1

      proposal_prob_single = p_single_excit(sing_num,sing_den-sing_num)
      proposal_prob_double = 1._rk - proposal_prob_single

      if (rannyu()<proposal_prob_single) then ! single excitation

        excite_level(1) = 1

        normalization = e1_prob_sav / c_e2_prob_sav
        call generate_heatbath_single(det_i_up,det_i_dn,det_j_up(1),det_j_dn(1),excite_from_1,excite_to_1,matrix_element,proposal_prob,normalization)

      else ! double excitation

        excite_level(1) = 2

        excite_to_2 = choose_second_hole(excite_from_1,excite_from_2,excite_to_1,same_spin)
        ! Remove moves to already-occupied orbitals
        if (excite_to_2<=norb) then
          if (btest(det_i_up,excite_to_2-1)) then
            n_new_dets = 0
            weight_j(1) = 0._rk
            return
          endif
        else
          if (btest(det_i_dn,excite_to_2-norb-1)) then
            n_new_dets = 0
            weight_j(1) = 0._rk
            return
          endif
        endif

        prob_1_then_2 = proposal_prob
        call generate_heatbath_double(det_i_up,det_i_dn,det_j_up(1),det_j_dn(1),excite_from_1,excite_from_2,excite_to_1,excite_to_2,occ_up,occ_dn,matrix_element,prob_1_then_2,proposal_prob)

      endif ! single or double excitation

      if (time_sym)  call apply_time_reversal_symmetry(det_i_up,det_i_dn,det_j_up(1),det_j_dn(1),occ_up,occ_dn,matrix_element,proposal_prob)
      weight_j(1) = -tau * matrix_element / proposal_prob

    endif ! n_new_dets

  end subroutine off_diagonal_move_chem_efficient_heatbath
!=============================================================

  subroutine apply_time_reversal_symmetry(det_i_up,det_i_dn,det_j_up,det_j_dn,occ_up,occ_dn,matrix_element,proposal_prob)
  ! Checks to see if the move from i_up,i_dn to j_dn,j_up is also possible,
  ! and if so, updates the matrix element and proposal probability to include this possibility.
  ! Also, replaces det_j with its representative (i.e., flips det_j_up and det_j_dn if need be)
  ! Original symmetry code written by Hitesh
  ! Moved here by A Holmes, 23 Sep 2015

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_i_up,det_i_dn
    type(ik_vec),intent(inout) :: det_j_up,det_j_dn
    type(ik_vec) :: tmp_det
#else
    integer(ik),intent(in) :: det_i_up,det_i_dn
    integer(ik),intent(inout) :: det_j_up,det_j_dn
    integer(ik) :: tmp_det
#endif
    integer,intent(in) :: occ_up(:),occ_dn(:)
    real(rk),intent(inout) :: matrix_element,proposal_prob

    real(rk) :: matrix_element2,proposal_prob_sym
    integer :: excite_level_sym
    real(rk) :: norm_i,norm_j


      if ((det_j_up .eq. det_i_up) .and. (det_j_dn .eq. det_i_dn)) then  ! Already accounted for in diagonal
        matrix_element = 0._rk
        return
      endif
      if ((det_j_dn .eq. det_i_up) .and. (det_j_up .eq. det_i_dn)) then  ! Already accounted for in diagonal
        matrix_element = 0._rk
        return
      endif

      if (det_i_up .eq. det_i_dn) then
        if (z .eq. 1) then
          norm_i=sqrt2
        else
          norm_i=0._rk
          write (6,*) "This state is inconsistent with the time symmetry you put in and should never enter this routine as an incoming state"
          stop
        endif
      else
        norm_i=1._rk
      endif

      norm_j=1._rk
      if (det_j_up .eq. det_j_dn) then
         if (z .eq. 1) then
           norm_j=sqrt2
         else
           norm_j=0._rk
           return
         endif
         matrix_element=(norm_j/norm_i)*matrix_element
      else
         norm_j=1._rk
         call excitation_level(det_i_up,det_i_dn,det_j_dn,det_j_up,excite_level_sym)
          if (excite_level_sym>=0) then
            call hamiltonian_chem(det_i_up, det_i_dn, det_j_dn, det_j_up, excite_level_sym, matrix_element2)
            if (abs(matrix_element2) .gt. 1.0e-10) then
              call proposal_prob_efficient_heatbath(det_i_up,det_i_dn,det_j_dn,det_j_up,occ_up,occ_dn,excite_level_sym,matrix_element2,proposal_prob_sym) ! matrix_element2 used only if single excitation
              proposal_prob = proposal_prob + proposal_prob_sym
              matrix_element=(norm_j/norm_i)*(matrix_element+z*matrix_element2)
            else
              matrix_element=(norm_j/norm_i)*matrix_element
            endif
         else
            matrix_element=(norm_j/norm_i)*(matrix_element)
         endif
      endif
      if (det_j_up .gt. det_j_dn) then                 ! Convert to its representative
         tmp_det=det_j_up
         det_j_up=det_j_dn
         det_j_dn=tmp_det
         matrix_element=matrix_element*z
      endif

  end subroutine apply_time_reversal_symmetry


  !=============================================================
  subroutine proposal_prob_efficient_heatbath(det_i_up,det_i_dn,det_j_up,det_j_dn,occ_up,occ_dn,excite_level,off_diag_elem,proposal_prob)
    ! Returns the efficient heatbath proposal probability of proposing this off-diagonal move.
    ! Note: Should only be used in cases where this move has not already been proposed, (e.g.,
    !       time-reversal symmetry), since it re-computes some of the same sums that are needed
    !       to propose the move.
    ! A Holmes, 16 Apr 2015

    use mpi_routines, only : mpi_stop

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_i_up,det_i_dn,det_j_up,det_j_dn
    type(ik_vec) :: tmp_det
#else
    integer(ik),intent(in) :: det_i_up,det_i_dn,det_j_up,det_j_dn
    integer(ik) :: tmp_det
#endif
    integer,intent(in) :: occ_up(:)
    integer,intent(in) :: occ_dn(:)
    integer,intent(in) :: excite_level
    real(rk),intent(in) :: off_diag_elem
    real(rk),intent(out) :: proposal_prob

    integer :: excite_from_1,excite_from_2,excite_to_1,excite_to_2
    integer :: excite_spin,i

    real(rk) :: prob_1_then_2
    real(rk) :: one_elec_prob_1
    real(rk) :: prob_of_choosing_1_and_2
    real(rk) :: sum_of_one_elec_probs
    real(rk) :: sum_of_probs_of_choosing_1_and_other

    real(rk) :: normalization,single_elem
    logical :: same_spin

    ! get excite_from_1,excite_from_2,excite_to_1,excite_to_2
    if (excite_level==1) then
      if (det_i_up==det_j_up) then
        excite_spin = -1
        excite_from_1 = trailz(iand(det_i_dn,not(det_j_dn))) + norb + 1
        excite_to_1 = trailz(iand(det_j_dn,not(det_i_dn))) + norb + 1
      else
        excite_spin = 1
        excite_from_1 = trailz(iand(det_i_up,not(det_j_up))) + 1
        excite_to_1 = trailz(iand(det_j_up,not(det_i_up))) + 1
      endif
    else
      if (det_i_up==det_j_up) then
        excite_spin = -1
        tmp_det = iand(det_i_dn,not(det_j_dn))
        excite_from_1 = trailz(tmp_det) + 1
        excite_from_2 = trailz(ibclr(tmp_det,excite_from_1-1)) + 1
        tmp_det = iand(det_j_dn,not(det_i_dn))
        excite_to_1 = trailz(tmp_det) + 1
        excite_to_2 = trailz(ibclr(tmp_det,excite_to_1-1)) + 1
        excite_from_1 = excite_from_1 + norb
        excite_from_2 = excite_from_2 + norb
        excite_to_1 = excite_to_1 + norb
        excite_to_2 = excite_to_2 + norb
      elseif (det_i_dn==det_j_dn) then
        excite_spin = 1
        tmp_det = iand(det_i_up,not(det_j_up))
        excite_from_1 = trailz(tmp_det) + 1
        excite_from_2 = trailz(ibclr(tmp_det,excite_from_1-1)) + 1
        tmp_det = iand(det_j_up,not(det_i_up))
        excite_to_1 = trailz(tmp_det) + 1
        excite_to_2 = trailz(ibclr(tmp_det,excite_to_1-1)) + 1
      else
        excite_spin = 0
        excite_from_1 = trailz(iand(det_i_up,not(det_j_up))) + 1
        excite_from_2 = trailz(iand(det_i_dn,not(det_j_dn))) + norb + 1
        excite_to_1 = trailz(iand(det_j_up,not(det_i_up))) + 1
        excite_to_2 = trailz(iand(det_j_dn,not(det_i_dn))) + norb + 1
      endif
    endif

    proposal_prob = 0._rk
    sum_of_one_elec_probs = 0._rk
    sum_of_probs_of_choosing_1_and_other = 0._rk
    if (excite_from_1<=norb) then
      do i=1,nup
        sum_of_one_elec_probs = sum_of_one_elec_probs + one_orbital_probabilities(occ_up(i))
        if (excite_from_1==occ_up(i)) then
          one_elec_prob_1 = one_orbital_probabilities(occ_up(i))
          cycle
        endif
        sum_of_probs_of_choosing_1_and_other = sum_of_probs_of_choosing_1_and_other + two_orbital_probabilities(excite_from_1,occ_up(i))
      enddo
      do i=1,ndn
        sum_of_one_elec_probs = sum_of_one_elec_probs + one_orbital_probabilities(occ_dn(i))
        sum_of_probs_of_choosing_1_and_other = sum_of_probs_of_choosing_1_and_other + two_orbital_probabilities(excite_from_1,occ_dn(i)+norb)
      enddo
    else
      do i=1,nup
        sum_of_one_elec_probs = sum_of_one_elec_probs + one_orbital_probabilities(occ_up(i))
        sum_of_probs_of_choosing_1_and_other = sum_of_probs_of_choosing_1_and_other + two_orbital_probabilities(excite_from_1,occ_up(i))
      enddo
      do i=1,ndn
        sum_of_one_elec_probs = sum_of_one_elec_probs + one_orbital_probabilities(occ_dn(i))
        if (excite_from_1==occ_dn(i)+norb) then
          one_elec_prob_1 = one_orbital_probabilities(occ_dn(i))
          cycle
        endif
        sum_of_probs_of_choosing_1_and_other = sum_of_probs_of_choosing_1_and_other + two_orbital_probabilities(excite_from_1,occ_dn(i)+norb)
      enddo
    endif

    if (excite_level==1) then
      normalization = one_elec_prob_1 / sum_of_one_elec_probs / sum_of_probs_of_choosing_1_and_other ! This is the probability of selecting the first electron divided by the sum of probabilities of selecting any second electron
      call prob_heatbath_single(excite_from_1,excite_to_1,off_diag_elem,proposal_prob,normalization)
    else ! double excitation
      single_elem = compute_single_elem(det_i_up,det_i_dn,excite_from_1,excite_to_1)
      prob_of_choosing_1_and_2 = two_orbital_probabilities(excite_from_1,excite_from_2)
      prob_1_then_2 = one_elec_prob_1 / sum_of_one_elec_probs * prob_of_choosing_1_and_2 / sum_of_probs_of_choosing_1_and_other ! This is the probability of selecting the ordered pair of electrons 1,2
      same_spin = (excite_spin.ne.0)
      call prob_heatbath_double(det_i_up,det_i_dn,excite_from_1,excite_from_2,excite_to_1,excite_to_2,occ_up,occ_dn,single_elem,prob_1_then_2,proposal_prob,same_spin)
    endif


  end subroutine proposal_prob_efficient_heatbath
!===========================================================================


  subroutine generate_heatbath_single(det_i_up,det_i_dn,det_j_up,det_j_dn,excite_from_1,excite_to_1,matrix_element,proposal_prob,normalization)

   ! To be called by off_diagonal_move_chem_efficient_heatbath
   ! A Holmes, 14 Sep 2015

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_i_up,det_i_dn
    type(ik_vec),intent(out) :: det_j_up,det_j_dn
#else
    integer(ik),intent(in) :: det_i_up,det_i_dn
    integer(ik),intent(out) :: det_j_up,det_j_dn
#endif
    integer,intent(in) :: excite_from_1,excite_to_1
    real(rk),intent(in) :: matrix_element
    real(rk),intent(out) :: proposal_prob
    real(rk),intent(in) :: normalization ! If choose_pairs, then normalization = 1 / c_pairs_prob(n_pairs) / 2._rk; else, normalization = e1_prob_sav / c_e2_prob_sav

      ! Generate new determinants
      if (excite_from_1<=norb) then ! up excitation
        det_j_up = ibclr(det_i_up,excite_from_1-1)
        det_j_up = ibset(det_j_up,excite_to_1-1)
        det_j_dn = det_i_dn
      else ! dn excitation
        det_j_dn = ibclr(det_i_dn,excite_from_1-norb-1)
        det_j_dn = ibset(det_j_dn,excite_to_1-norb-1)
        det_j_up = det_i_up
      endif

      ! Compute proposal probability
      call prob_heatbath_single(excite_from_1,excite_to_1,matrix_element,proposal_prob,normalization)

  end subroutine generate_heatbath_single


  subroutine prob_heatbath_single(excite_from_1,excite_to_1,matrix_element,proposal_prob,normalization)

   ! To be called by off_diagonal_move_chem_efficient_heatbath
   ! A Holmes, 23 Sep 2015

    integer,intent(in) :: excite_from_1,excite_to_1
    real(rk),intent(in) :: matrix_element
    real(rk),intent(out) :: proposal_prob
    real(rk),intent(in) :: normalization ! If choose_pairs, then normalization = 1 / c_pairs_prob(n_pairs) / 2._rk; else, normalization = e1_prob_sav / c_e2_prob_sav
    real(rk) :: sing_num
    integer :: i

      ! Compute proposal probability
      proposal_prob = 0._rk
      sing_num = abs(matrix_element) ! Single excitation numerator for calculating probability of choosing single rather than double

      if (excite_from_1<=norb) then
        do i=1,nup
          if (excite_from_1==occ_up(i))  cycle
          proposal_prob = proposal_prob + two_orbital_probabilities(excite_from_1,occ_up(i))*three_orbital_probabilities_same_spin(excite_from_1,occ_up(i),excite_to_1)*p_single_excit(sing_num,Htot_same(combine_2_indices(excite_from_1,occ_up(i)),excite_to_1))
        enddo
        do i=1,ndn
          proposal_prob = proposal_prob + two_orbital_probabilities(excite_from_1,occ_dn(i)+norb)*three_orbital_probabilities_opposite_spin(excite_from_1,occ_dn(i),excite_to_1)*p_single_excit(sing_num,Htot_opposite(excite_from_1,occ_dn(i),excite_to_1))
        enddo
      else
        do i=1,nup
          proposal_prob = proposal_prob + two_orbital_probabilities(excite_from_1,occ_up(i))*three_orbital_probabilities_opposite_spin(excite_from_1-norb,occ_up(i),excite_to_1-norb)*p_single_excit(sing_num,Htot_opposite(excite_from_1-norb,occ_up(i),excite_to_1-norb))
        enddo
        do i=1,ndn
          if (excite_from_1==occ_dn(i)+norb)  cycle
          proposal_prob = proposal_prob + two_orbital_probabilities(excite_from_1,occ_dn(i)+norb)*three_orbital_probabilities_same_spin(excite_from_1-norb,occ_dn(i),excite_to_1-norb)*p_single_excit(sing_num,Htot_same(combine_2_indices(excite_from_1-norb,occ_dn(i)),excite_to_1-norb))
        enddo
      endif

      proposal_prob = proposal_prob * normalization


  end subroutine prob_heatbath_single


  subroutine generate_heatbath_double(det_i_up,det_i_dn,det_j_up,det_j_dn,excite_from_1,excite_from_2,excite_to_1,excite_to_2,occ_up,occ_dn,matrix_element,prob_1_then_2,proposal_prob)

   ! To be called by off_diagonal_move_chem_efficient_heatbath.
   ! proposal_prob:  input is probability of choosing ordered pair of electrons
   !                 output is probability of choosing this double excitation
   ! matrix_element: input is matrix element of single excitation excite_from_1 to excite_to_1
   !                 output is matrix element of double excitation
   ! A Holmes, 14 Sep 2015

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_i_up,det_i_dn
    type(ik_vec),intent(out) :: det_j_up,det_j_dn
#else
    integer(ik),intent(in) :: det_i_up,det_i_dn
    integer(ik),intent(out) :: det_j_up,det_j_dn
#endif
    integer,intent(in) :: excite_from_1,excite_to_1
    integer,intent(in) :: excite_from_2,excite_to_2
    integer,intent(in) :: occ_up(:),occ_dn(:)
    real(rk),intent(inout) :: matrix_element
    real(rk),intent(in) :: prob_1_then_2
    real(rk),intent(out) :: proposal_prob

    integer :: excite_spin ! number of up excitations - number of dn excitations
    logical :: same_spin

      ! Generate new determinants
      det_j_up = det_i_up
      det_j_dn = det_i_dn

      excite_spin = 0
      if (excite_from_1<=norb) then
        det_j_up = ibclr(det_j_up,excite_from_1-1)
        excite_spin = excite_spin + 1
      else
        det_j_dn = ibclr(det_j_dn,excite_from_1-norb-1)
        excite_spin = excite_spin - 1
      endif
      if (excite_from_2<=norb) then
        det_j_up = ibclr(det_j_up,excite_from_2-1)
        excite_spin = excite_spin + 1
      else
        det_j_dn = ibclr(det_j_dn,excite_from_2-norb-1)
        excite_spin = excite_spin - 1
      endif

      if (excite_to_1<=norb) then
        det_j_up = ibset(det_j_up,excite_to_1-1)
      else
        det_j_dn = ibset(det_j_dn,excite_to_1-norb-1)
      endif
      if (excite_to_2<=norb) then
        det_j_up = ibset(det_j_up,excite_to_2-1)
      else
        det_j_dn = ibset(det_j_dn,excite_to_2-norb-1)
      endif

      ! Compute proposal probability
      same_spin = (excite_spin.ne.0)
      call prob_heatbath_double(det_i_up,det_i_dn,excite_from_1,excite_from_2,excite_to_1,excite_to_2,occ_up,occ_dn,matrix_element,prob_1_then_2,proposal_prob,same_spin)
      ! Compute matrix element
      call hamiltonian_chem(det_i_up,det_i_dn,det_j_up,det_j_dn,2,matrix_element)


  end subroutine generate_heatbath_double


  subroutine prob_heatbath_double(det_i_up,det_i_dn,excite_from_1,excite_from_2,excite_to_1,excite_to_2,occ_up,occ_dn,matrix_element,prob_1_then_2,proposal_prob,same_spin)

   ! To be called by off_diagonal_move_chem_efficient_heatbath.
   ! proposal_prob:  input is probability of choosing ordered pair of electrons
   !                 output is probability of choosing this double excitation
   ! matrix_element: input is matrix element of single excitation excite_from_1 to excite_to_1
   ! A Holmes, 23 Sep 2015

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_i_up,det_i_dn
#else
    integer(ik),intent(in) :: det_i_up,det_i_dn
#endif
    integer,intent(in) :: excite_from_1,excite_to_1
    integer,intent(in) :: excite_from_2,excite_to_2
    integer,intent(in) :: occ_up(:),occ_dn(:)
    real(rk),intent(in) :: matrix_element
    real(rk),intent(in) :: prob_1_then_2 ! P(choose 1 then 2); prob_2_then_1 = P(choose 2 then 1)
    real(rk),intent(out) :: proposal_prob
    logical,intent(in) :: same_spin

    real(rk) :: prob_2_then_1 ! P(choose 2 then 1)
    real(rk) :: sing_num ! single excitation probability numerator
    real(rk) :: term1,term2,term3,term4
    integer :: i
    real(rk) :: one_elec_prob_2
    real(rk) :: prob_of_choosing_1_and_2
    real(rk) :: sum_of_one_elec_probs
    real(rk) :: sum_of_probs_of_choosing_2_and_other

      ! First get probability of selecting electron 2 first
      if (excite_from_2<=norb) then
        sum_of_probs_of_choosing_2_and_other = 0._rk
        do i=1,nup
          sum_of_one_elec_probs = sum_of_one_elec_probs + one_orbital_probabilities(occ_up(i))
          if (excite_from_2==occ_up(i)) then
            one_elec_prob_2 = one_orbital_probabilities(occ_up(i))
            cycle
          endif
          sum_of_probs_of_choosing_2_and_other = sum_of_probs_of_choosing_2_and_other + two_orbital_probabilities(excite_from_2,occ_up(i))
        enddo
        do i=1,ndn
          sum_of_one_elec_probs = sum_of_one_elec_probs + one_orbital_probabilities(occ_dn(i))
          sum_of_probs_of_choosing_2_and_other = sum_of_probs_of_choosing_2_and_other + two_orbital_probabilities(excite_from_2,occ_dn(i)+norb)
        enddo
      else
        do i=1,nup
          sum_of_one_elec_probs = sum_of_one_elec_probs + one_orbital_probabilities(occ_up(i))
          sum_of_probs_of_choosing_2_and_other = sum_of_probs_of_choosing_2_and_other + two_orbital_probabilities(excite_from_2,occ_up(i))
        enddo
        do i=1,ndn
          sum_of_one_elec_probs = sum_of_one_elec_probs + one_orbital_probabilities(occ_dn(i))
          if (excite_from_2==occ_dn(i)+norb) then
            one_elec_prob_2 = one_orbital_probabilities(occ_dn(i))
            cycle
          endif
          sum_of_probs_of_choosing_2_and_other = sum_of_probs_of_choosing_2_and_other + two_orbital_probabilities(excite_from_2,occ_dn(i)+norb)
        enddo
      endif
      prob_of_choosing_1_and_2 = two_orbital_probabilities(excite_from_1,excite_from_2)
      prob_2_then_1 = one_elec_prob_2 / sum_of_one_elec_probs * prob_of_choosing_1_and_2 / sum_of_probs_of_choosing_2_and_other ! This is the probability of selecting the ordered pair of electrons 2,1


      ! Term 1: Choose from_1->to_1, then from_2->to_2
      term1 = p_first_hole(excite_from_1,excite_from_2,excite_to_1)*p_second_hole(excite_from_1,excite_from_2,excite_to_1,excite_to_2)
      sing_num = abs(matrix_element) ! = abs(compute_single_elem(det_i_up,det_i_dn,excite_from_1,excite_to_1)) - checked.
      term1 = term1 * (p_double_excit(sing_num,Htot(excite_from_1,excite_from_2,excite_to_1)))

      ! Term 2: Choose from_1->to_2, then from_2->to_1
      if (same_spin) then
        term2 = p_first_hole(excite_from_1,excite_from_2,excite_to_2)*p_second_hole(excite_from_1,excite_from_2,excite_to_2,excite_to_1)
        sing_num = abs(compute_single_elem(det_i_up,det_i_dn,excite_from_1,excite_to_2))
        term2 = term2 * (p_double_excit(sing_num,Htot(excite_from_1,excite_from_2,excite_to_2)))
      else
        term2 = 0._rk
      endif

      ! Term 3: Choose from_2->to_1, then from_1->to_2
      if (same_spin) then
        term3 = p_first_hole(excite_from_2,excite_from_1,excite_to_1)*p_second_hole(excite_from_2,excite_from_1,excite_to_1,excite_to_2)
        sing_num = abs(compute_single_elem(det_i_up,det_i_dn,excite_from_2,excite_to_1))
        term3 = term3 * (p_double_excit(sing_num,Htot(excite_from_2,excite_from_1,excite_to_1)))
      else
        term3 = 0._rk
      endif

      ! Term 4: Choose from_2->to_2, then from_1->to_1
      term4 = p_first_hole(excite_from_2,excite_from_1,excite_to_2)*p_second_hole(excite_from_2,excite_from_1,excite_to_2,excite_to_1)
      sing_num = abs(compute_single_elem(det_i_up,det_i_dn,excite_from_2,excite_to_2))
      term4 = term4 * (p_double_excit(sing_num,Htot(excite_from_2,excite_from_1,excite_to_2)))

      proposal_prob = prob_1_then_2 * (term1 + term2) + prob_2_then_1 * (term3 + term4)

  end subroutine prob_heatbath_double


  !===========================================================================
  real(rk) function p_single_excit(single_elem, sum_of_double_elems)
    !---------------------------------------------------------------------------
    ! Return the probability of selecting a single excitation once two electrons and one hole have already been selected
    ! A Holmes, 26 Apr 2015
    ! Modified by A Holmes on 9 Oct 2015 to work with proposal of two off-diagonal moves
    !---------------------------------------------------------------------------
    real(rk),intent(in) :: single_elem, sum_of_double_elems

      p_single_excit = abs(single_elem)/(abs(single_elem)+sum_of_double_elems)
      if (p_single_excit>0.5_rk)  p_single_excit = 1._rk ! Large single excitation element, so propose two moves

  end function p_single_excit
!===========================================================================

  !===========================================================================
  real(rk) function p_double_excit(single_elem, sum_of_double_elems)
    !---------------------------------------------------------------------------
    ! Return the probability of selecting a double excitation once two electrons and one hole have already been selected
    ! A Holmes, 9 Oct 2015
    !---------------------------------------------------------------------------
    real(rk),intent(in) :: single_elem, sum_of_double_elems

      p_double_excit = sum_of_double_elems/(abs(single_elem)+sum_of_double_elems)
      if (p_double_excit<0.5_rk)  p_double_excit = 1._rk ! Large single excitation element, so propose two moves

  end function p_double_excit
!===========================================================================


  !===========================================================================
  subroutine energy_pieces_chem(det_i_up, det_i_dn, e_mix_numerator, e_mix_denominator)
    !---------------------------------------------------------------------------
    ! Description : Calculate pieces of the local energy for det_i with weight_i
    !               Numerator is sum_j=1^{ndet_psi_t} H_ij  *  trial_wf_weight on det_j
    !               Denominator is weight_i * trial_wf_weight on det_i
    !               Currently only being used at start of run (local energies are being stored instead)
    !
    ! Created     : F. Petruzielo, 9 Nov 2010
    ! Modified    : for time symmetry, by Hitesh Changlani April 24 2012
    !---------------------------------------------------------------------------
    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
#endif
    real(rk), intent(out) :: e_mix_numerator
    real(rk), intent(out) :: e_mix_denominator

    !local variables

    integer :: j, excite_level
    real(rk) :: matrix_element

    e_mix_numerator = 0._rk
    e_mix_denominator = 0._rk

    do j = 1, ndet_psi_t
       if (time_sym) then                                   ! If time reversal symmetry, use the symmetrized Hamiltonian matrix element
            call hamiltonian_chem_time_sym(det_i_up, det_i_dn, dets_up_psi_t(j), dets_dn_psi_t(j), matrix_element)
            e_mix_numerator = e_mix_numerator + matrix_element * cdet_psi_t(j)
            !write (6,*) "matrix element",matrix_element
            !write (6,*) "cdet_psi_t",cdet_psi_t(j)
            !call flush(6)
            if (dets_up_psi_t(j) == det_i_up .and. dets_dn_psi_t(j) == det_i_dn) then
                     !det_i is in the trial wavefunction
                     e_mix_denominator = cdet_psi_t(j)
                     !write (6,*) "emix_denominator",e_mix_denominator
                     !call flush(6)
            endif
       else                                                 ! Not touching Frank's unsymmetrized code
           call excitation_level(det_i_up, det_i_dn, dets_up_psi_t(j), dets_dn_psi_t(j), excite_level)
           if (excite_level >= 0) then
              !det_i is connected to a det in the trial wavefunction
              call hamiltonian_chem(det_i_up, det_i_dn, dets_up_psi_t(j), dets_dn_psi_t(j), excite_level, matrix_element)
              e_mix_numerator = e_mix_numerator + matrix_element * cdet_psi_t(j)
              if (dets_up_psi_t(j) == det_i_up .and. dets_dn_psi_t(j) == det_i_dn) then
                 !det_i is in the trial wavefunction
                 e_mix_denominator = cdet_psi_t(j)
              endif
           endif
       endif
    enddo

  end subroutine energy_pieces_chem
  !=============================================================

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

    !remove core
    det_spin = ishft(det_spin, -n_core_orb)
    !next lexographically ordered bit pattern from Bit Twiddling Hacks
    temp = ior(det_spin, det_spin - 1) + 1
    det_spin = ior(temp,ishft(iand(temp, -temp) / iand(det_spin,-det_spin),-1) - 1)
    !add back core
    det_spin = ishft(det_spin, n_core_orb) + 2_ik ** n_core_orb - 1

  end subroutine twiddle_determinant

  !=============================================================

  !=============================================================
  subroutine diagonalize_hamiltonian(dets_up, dets_dn, lowest_eigenvector, lowest_eigenvalue)
    !---------------------------------------------------------------------------
    ! Description : Construct and diagonalize H in a space of dets (broken up into dets_up and dets_dn)
    !               Sort ground state eigenvector by decreasing magnitude (HF first for example).
    !               Put dets_up and dets_dn in the same order.
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    !---------------------------------------------------------------------------
    use types, only : num_words
    implicit none

    !dummy arguments
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(inout) :: dets_up(:)
    type(ik_vec), intent(inout) :: dets_dn(:)
#else
    integer(ik), intent(inout) :: dets_up(:)
    integer(ik), intent(inout) :: dets_dn(:)
#endif
    real(rk), intent(out) :: lowest_eigenvector(:)
    real(rk), intent(out) :: lowest_eigenvalue

    !local variables
!   character(len=4) :: fmt !temp
    integer :: n_det, excite_level
    integer :: len_work, info
    integer :: i, j
    real(rk) tau_optimal_deterministic, tau_optimal_stochastic
    real(rk), parameter :: epsilon = 1.e-8_rk
    real(rk), allocatable :: ham(:,:)
    real(rk), allocatable :: eigenvalues(:), work(:)
    character(len=2) :: fmt

    n_det = size(dets_up)
    write (6,*) "Number of dets to diagonalize:",n_det
    call flush(6)
    allocate(ham(n_det, n_det))
    allocate(eigenvalues(n_det))

    !construct hamiltonian
    ham = 0._rk
    do i = 1, n_det
       if (time_sym) then
            call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), ham(i,i))  ! diagonal element
       else
            call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), 0, ham(i,i))  ! diagonal element
       endif
       do j = i+1, n_det
          if (time_sym) then
             call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), ham(i,j))
          else
              call excitation_level(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level)
              if (excite_level >= 0) then
                 !determinants are connected by a single or double excitation so evaluate matrix elementx
                 call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, ham(i,j))
              else
                 ham(i,j) = 0._rk
              endif
          endif
          !symmetric matrix
          ham(j,i) = ham(i,j)
       enddo
!      write(fmt, '(i4)') n_det
!      write(6,'(i8,i8,' // trim(fmt) // 'f16.12)')  dets_up(i), dets_dn(i), ham(i,:)
    enddo

!   call print_real_matrix(n_det,n_det,ham,n_det,n_det)

    write (6,'(/,''Performing Lanczos assuming dense H, using matrix_lanczos in more_tools, n_det='',i10)') n_det
    call flush(6)
    call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,ham)
    call sort(lowest_eigenvector, dets_up, dets_dn) ! Sort in order of decreasing absolute coef
    write(6,'(/,''Wavefunction has'',i9,'' determinants. Dets. with largest'',i3,'' abs coefs are:'',/,''     dets_up             dets_dn        dets_up      dets_dn  excit_level cdet'')')  n_det,min(n_det,20)
    write (fmt, '(i2)') 2*num_words
    write(6,'(' // trim(fmt) // 'b20,' // trim(fmt) // 'i25,i3,f13.9)') (dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), popcnt(iand(dets_up(i),not(dets_up(1))))+popcnt(iand(dets_dn(i),not(dets_dn(1)))), lowest_eigenvector(i),i=1,min(n_det,20))
    write(6,'(/,''n_det, Lowest eigenvalue ='',i10,f12.6)') n_det, lowest_eigenvalue ; call flush(6)

    if (n_det<1000) then
       !construct hamiltonian
       ham = 0._rk
       do i = 1, n_det
          if (time_sym) then
            call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), ham(i,i))  ! diagonal element
          else
            call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), 0, ham(i,i))  ! diagonal element
          endif
          do j = i+1, n_det
             if (time_sym) then
                call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), ham(i,j))  ! diagonal element
             else
                 call excitation_level(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level)
                 if (excite_level >= 0) then
                    !determinants are connected by a single or double excitation so evaluate matrix elementx
                    call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, ham(i,j))
                 else
                    ham(i,j) = 0._rk
                 endif
             endif
             !symmetric matrix
             ham(j,i) = ham(i,j)
          enddo
          !      write(fmt, '(i4)') n_det
          !      write(6,'(i8,i8,' // trim(fmt) // 'f16.12)')  dets_up(i), dets_dn(i), ham(i,:)
       enddo

      !diagonalize with lapack routine
      len_work = 3 * n_det -1
      allocate(work(len_work))

      write(6,'(/,''n_det, len_work='',9i8)') n_det,len_work ; call flush(6)

      write (6,'(''Performing LAPACK diag, n_det='')') n_det
      call dsyev('V', 'U', n_det, ham, n_det, eigenvalues, work, len_work, info)
      if (info /= 0) then
         write(6,*) "info = ", info
         stop 'Diagonalization Failed!'
      endif

      lowest_eigenvector = ham(:,1)
      lowest_eigenvalue = eigenvalues(1)
      call sort(lowest_eigenvector, dets_up, dets_dn) ! Sort in order of decreasing absolute coef

      write(6,'(/,''Wavefunction has'',i9,'' determinants. Dets. with largest'',i3,'' abs coefs are:'',/,''     dets_up             dets_dn        dets_up      dets_dn  excit_level cdet'')')  n_det,min(n_det,20)
      write (fmt, '(i2)') 2*num_words
      write(6,'(' // trim(fmt) // 'b20,' // trim(fmt) // 'i25,i3,f13.9)') (dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), popcnt(iand(dets_up(i),not(dets_up(1))))+popcnt(iand(dets_dn(i),not(dets_dn(1)))), lowest_eigenvector(i),i=1,min(n_det,20))
      write(6,'(/,''Lowest'',i3,'' eigenvalues='',20f13.6)') min(n_det,20),(eigenvalues(i),i=1,min(n_det,20))
      write(6,'(''Lowest, highest eigenvalues='', 9f13.6)') eigenvalues(1),eigenvalues(n_det)

    endif

! Calculate what tau_optimal_deterministic and tau_optimal_stochastic would be if the wavefn were this trial wavefn.
    tau_optimal_stochastic=1/(eigenvalues(n_det)-eigenvalues(1))
    if(n_det.ge.2) then
      tau_optimal_deterministic=2/(eigenvalues(n_det)+eigenvalues(2)-2*eigenvalues(1))
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_deterministic, tau_optimal_stochastic='',9f10.6)') &
   &  tau_optimal_deterministic, tau_optimal_stochastic
     else
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_stochastic='',9f10.6)') tau_optimal_stochastic
    endif

  end subroutine diagonalize_hamiltonian
  !=============================================================


  subroutine diagonalize_sparse_hamiltonian_chem_excited(dets_up, dets_dn, n_det, n_states, lowest_eigenvectors, lowest_eigenvalues, time_sym_on, initial_vectors, sparse_ham, remote_det_map, local_det_map, average_connections,use_hash)

    use types, only : i8b, num_words
    use more_tools, only : real_symmetric_diagonalize_ow_ham,davidson_sparse,davidson_sparse_mpi,parpack_diagonalize,&
         davidson_sparse_mpi2
    use tools, only : print_excitation_levels_and_wts
    use common_selected_ci, only : max_nonzero_elements,sparse_mat
    use mpi_routines, only : det_map,det_map_l
    implicit none

    !dummy arguments
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(inout) :: dets_up(:), dets_dn(:)
#else
    integer(ik), intent(inout) :: dets_up(:), dets_dn(:)
#endif
    integer,intent(in) :: n_det,n_states
    real(rk),optional,intent(inout) :: average_connections
    real(rk), intent(out) :: lowest_eigenvectors(:,:)
    real(rk), intent(out) :: lowest_eigenvalues(:)
    logical,intent(in)   :: time_sym_on
    real(rk),intent(in) :: initial_vectors(:,:)
    type(det_map), optional,intent(inout) :: remote_det_map
    type(det_map_l), optional,intent(inout) :: local_det_map
    logical, optional, intent(in)  :: use_hash
    type(sparse_mat),intent(inout) :: sparse_ham

    !local variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),allocatable :: dets_up_srt(:), dets_dn_srt(:)
#else
    integer(ik),allocatable :: dets_up_srt(:), dets_dn_srt(:)
#endif
    real(rk), allocatable :: lowest_eigenvector_srt(:)
    real(rk), parameter :: epsilon = 1.e-8_rk
    real(rk), allocatable :: H_values(:)
    integer(i8b), allocatable :: H_indices(:),H_nonzero_elements(:)
    real(rk) :: highest_eigenvalue,second_lowest_eigenvalue
    logical  :: sym_on
    real(rk),allocatable :: ham(:,:),eigenvalues(:)
    real(rk) :: matrix_element
    integer :: i,j,excite_level
    integer(i8b) :: ndets_global
    integer :: n_connected_dets
    character(len=2) :: fmt

    call my_second(1,'diagonalize_sparse_hamiltonian_chem')

    sym_on=.false.
    if (time_sym_on .and. time_sym) sym_on=time_sym_on

    if (present(remote_det_map)) then
      if (present(use_hash)) then
        call generate_sparse_ham_chem_upper_triangular_mpi(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,H_values,sym_on,hf_to_psit=.false.,sparse_ham=sparse_ham,remote_det_map=remote_det_map,local_det_map=local_det_map,use_hash=use_hash)
      else
        call generate_sparse_ham_chem_upper_triangular_mpi(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,H_values,sym_on,hf_to_psit=.false.,sparse_ham=sparse_ham,remote_det_map=remote_det_map,local_det_map=local_det_map)
      endif
    else
      call generate_sparse_ham_chem_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,H_values,sym_on,hf_to_psit=.false.,sparse_ham=sparse_ham,average_connections=average_connections) ! Pass in old H, return new H
    endif


    write (6,'(/,''Performing Davidson using sparsely stored H and matrix_lanczos in more_tools, n_det='',i10,es10.2)') n_det, real(n_det)
    call flush(6)
    if (present(remote_det_map)) then
      call davidson_sparse_mpi2(remote_det_map,local_det_map,n_states,lowest_eigenvectors,lowest_eigenvalues,H_indices,H_nonzero_elements,H_values,initial_vectors)
       !call parpack_diagonalize(remote_det_map,local_det_map,n_states,lowest_eigenvectors,lowest_eigenvalues,H_indices,H_nonzero_elements,H_values,initial_vectors)
    else
       call davidson_sparse(n_det,n_states,lowest_eigenvectors,lowest_eigenvalues,H_indices,H_nonzero_elements,H_values,initial_vectors)
    endif
   !call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
    deallocate(H_values)
    deallocate(H_indices)
    deallocate(H_nonzero_elements)

    if(present(remote_det_map)) then
      !MJO Need to work with the global list in parallel, rather than the local lists
      ndets_global = local_det_map%ndets_global
      if(.true.) then !ndets_global.le.100000) then
        allocate(lowest_eigenvector_srt(ndets_global)) ; allocate(dets_up_srt(ndets_global)) ; allocate(dets_dn_srt(ndets_global))
        do j=1,n_states
          lowest_eigenvector_srt(1:ndets_global)=remote_det_map%remote_wts(1:ndets_global,j) !MJO FIXME Make this work
!          lowest_eigenvector_srt(1:ndets_global)=0 !lowest_eigenvectors(1:ndets,j)
          dets_up_srt(1:ndets_global)=remote_det_map%global_dets_up(1:ndets_global)
          dets_dn_srt(1:ndets_global)=remote_det_map%global_dets_dn(1:ndets_global)
          call sort(lowest_eigenvector_srt, dets_up_srt, dets_dn_srt) ! Sort in order of decreasing absolute coef
          write (6,*) "State",j,":"
          write(6,'(/,''Wavefunction has'',i9,'' determinants. Dets. with largest'',i3,'' abs coefs are:'',/,''     dets_up             dets_dn        dets_up      dets_dn  excit_level cdet'')')  ndets_global,min(ndets_global,20)
          write (fmt, '(i2)') 2*num_words
          write(6,'(' // trim(fmt) // 'b20,' // trim(fmt) // 'i25,i3,f13.9)') (dets_up_srt(i), dets_dn_srt(i), dets_up_srt(i), dets_dn_srt(i), popcnt(iand(dets_up_srt(i),not(local_det_map%hf_up)))+popcnt(iand(dets_dn_srt(i),not(local_det_map%hf_dn))), lowest_eigenvector_srt(i),i=1,min(ndets_global,20))
        enddo
        deallocate(lowest_eigenvector_srt) ; deallocate(dets_up_srt) ; deallocate(dets_dn_srt)
      else
        write(6,'(/,''Wavefunction has'',i9,'' determinants'')') ndets_global
      endif
      write(6,'(/,''Lowest eigenvalues='',10f13.6)') lowest_eigenvalues(1:n_states)
      call print_excitation_levels_and_wts(int(ndets_global,4),remote_det_map%global_dets_up,remote_det_map%global_dets_dn,remote_det_map%remote_wts(:,1),norb,orbital_symmetries,local_det_map%hf_up,local_det_map%hf_dn)
    else

      if(n_det.le.100000) then
        allocate(lowest_eigenvector_srt(n_det)) ; allocate(dets_up_srt(n_det)) ; allocate(dets_dn_srt(n_det))
        do j=1,n_states
          lowest_eigenvector_srt(1:n_det)=lowest_eigenvectors(1:n_det,j)
          dets_up_srt(1:n_det)=dets_up(1:n_det)
          dets_dn_srt(1:n_det)=dets_dn(1:n_det)
          call sort(lowest_eigenvector_srt, dets_up_srt, dets_dn_srt) ! Sort in order of decreasing absolute coef
          write (6,*) "State",j,":"
          write(6,'(/,''Wavefunction has'',i9,'' determinants. Dets. with largest'',i3,'' abs coefs are:'',/,''     dets_up             dets_dn        dets_up      dets_dn  excit_level cdet'')')  n_det,min(n_det,20)
          write (fmt, '(i2)') 2*num_words
          write(6,'(' // trim(fmt) // 'b20,' // trim(fmt) // 'i25,i3,f13.9)') (dets_up_srt(i), dets_dn_srt(i), dets_up_srt(i), dets_dn_srt(i), popcnt(iand(dets_up_srt(i),not(dets_up(1))))+popcnt(iand(dets_dn_srt(i),not(dets_dn(1)))), lowest_eigenvector_srt(i),i=1,min(n_det,20))
        enddo
        deallocate(lowest_eigenvector_srt) ; deallocate(dets_up_srt) ; deallocate(dets_dn_srt)
      else
        write(6,'(/,''Wavefunction has'',i9,'' determinants'')') n_det
      endif
      write(6,'(/,''Lowest eigenvalues='',10f13.6)') lowest_eigenvalues(1:n_states)
      call print_excitation_levels_and_wts(n_det,dets_up,dets_dn,lowest_eigenvectors(:,1),norb,orbital_symmetries)
    endif


! Calculate number of dets with various excitation levels and the sum of their squared wts.


    call my_second(2,'diagonalize_sparse_hamiltonian_chem')

  end subroutine diagonalize_sparse_hamiltonian_chem_excited


  !=================================================================================================================================
  subroutine diagonalize_sparse_hamiltonian_chem(dets_up, dets_dn, n_det, lowest_eigenvector, lowest_eigenvalue, tau_out, time_sym_on, sort_or_not, initial_vector, epsilon3, sparse_ham)
    !------------------------------------------------------------------------------------------------------------
    ! Description : Construct and diagonalize (find lowest eigenvector of) H in a space of dets (broken up into dets_up and dets_dn)
    !               Sort ground state eigenvector by decreasing magnitude (HF first for example).
    !               Put dets_up and dets_dn in the same order.
    !
    ! Created     : A Holmes, 7 Apr 2012 (same as diagonalize_hamiltonian, but stores H sparsely)
    ! Modified    : A Holmes, 10 Jan 2013. Added optional input sort_or_not, whether to sort the lowest_eigenvector.
    !             : A Holmes, 14 Feb 2013. Optional input initial_vector is the starting vector for Lanczos.
    !             : A Holmes, 12 Oct 2016. Optional input epsilon3 allows on the fly Lanczos to be fast by skipping over tiny contributions that are smaller than epsilon3
    !             : Warning! Final energy is dependent upon initial vector when epsilon3 is nonzero!
    !----------------------------------------------------------------------------------------------------------

    use types, only : i8b, num_words
    use more_tools, only : real_symmetric_diagonalize_ow_ham,cyrus_diagonalize_sparse,davidson_sparse_single
    use common_run, only: tau_multiplier
    use tools, only : print_excitation_levels_and_wts
    use common_selected_ci, only : max_nonzero_elements,sparse_mat
    implicit none

    !dummy arguments
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(inout) :: dets_up(:), dets_dn(:)
#else
    integer(ik), intent(inout) :: dets_up(:), dets_dn(:)
#endif
    integer,intent(in) :: n_det
    real(rk), intent(out) :: lowest_eigenvector(:)
    real(rk), intent(out) :: lowest_eigenvalue
    real(rk),intent(out) :: tau_out
    logical,intent(in)   :: time_sym_on
    logical,optional,intent(in) :: sort_or_not ! Whether to sort by decreasing order by magnitude
    real(rk),optional,intent(in) :: initial_vector(:)
    real(rk),optional,intent(in) :: epsilon3
    type(sparse_mat),optional,intent(inout) :: sparse_ham
    !logical,optional      :: time_sym_on

    !local variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),allocatable :: dets_up_srt(:), dets_dn_srt(:)
#else
    integer(ik),allocatable :: dets_up_srt(:), dets_dn_srt(:)
#endif
    real(rk), allocatable :: lowest_eigenvector_srt(:)
    real(rk), parameter :: epsilon = 1.e-8_rk
    real(rk), allocatable :: H_values(:)
    integer(i8b), allocatable :: H_indices(:),H_nonzero_elements(:)
    real(rk) :: highest_eigenvalue,second_lowest_eigenvalue
    real(rk) :: tau_optimal_deterministic, tau_optimal_stochastic
    logical  :: sym_on
    real(rk),allocatable :: ham(:,:),eigenvalues(:)
    real(rk) :: matrix_element
    integer :: i,j,excite_level
    real(rk) :: log_num_nonzero_elements
   !type(ik_vec),allocatable    :: temp_i16_up(:),temp_i16_dn(:)!,temp_i_2(:)
   !integer,allocatable        :: iorder(:),temp_i_2(:)
   !integer :: n
    logical               :: srt
    integer :: n_connected_dets
    character(len=2) :: fmt

    call my_second(1,'diagonalize_sparse_hamiltonian_chem')

    sym_on=.false.
    !if (time_sym .and. present(time_sym_on)) sym_on=time_sym_on
    if (time_sym_on .and. time_sym) sym_on=time_sym_on

    if (present(sort_or_not)) then
      srt=sort_or_not
    else
      srt=.true.
    endif

    if (.true.) then

      call find_connected_dets_chem(dets_up(1), dets_dn(1), n_connected_dets, connected_dets_up, connected_dets_dn, norb)

      ! Estimate whether we can store sparse matrix.
      log_num_nonzero_elements=log(real(n_det))+log(real(min(100,min(n_connected_dets,n_det))))-log(2.0)
      write (6,'(''N_con_HF='',i10)') n_connected_dets
      write(6,'(/,''n_det='',i12,es10.2,'' estim # of nonzero elem in H='',i14,es10.2,'', Max # of nonzero elem='',es10.2,'', estim average # of nonzero elem='',f8.1,/)') n_det, real(n_det), int(exp(log_num_nonzero_elements),i8b), exp(log_num_nonzero_elements), real(max_nonzero_elements), exp(log_num_nonzero_elements)/n_det

      if (exp(log_num_nonzero_elements)>max_nonzero_elements) then
        write (6,'(/,''Cannot store Hamiltonian sparsely for Lanczos, so use matrix_lanczos_on_the_fly, n_det='',i10,es10.2)') n_det, real(n_det) ; call flush(6)
        if (present(initial_vector)) then
          call matrix_lanczos_on_the_fly(dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,highest_eigenvalue,second_lowest_eigenvalue,sym_on=sym_on,initial_vector=initial_vector,eps=epsilon3)
        else
          call matrix_lanczos_on_the_fly(dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,highest_eigenvalue,second_lowest_eigenvalue,int(exp(log(real(max_nonzero_elements))-log_num_nonzero_elements+log(real(n_det)))),sym_on)
        endif
      else
!       write (6,*) "Using sparsely stored Hamiltonian for Lanczos"
        if (present(sparse_ham)) then
          call generate_sparse_ham_chem_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,H_values,sym_on,hf_to_psit=.false.,sparse_ham=sparse_ham) ! Pass in old H, return new H
        else
          call generate_sparse_ham_chem_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,H_values,sym_on,hf_to_psit=.false.)
        endif
        write (6,'(/,''Performing Lanczos using sparsely stored H and matrix_lanczos in more_tools, n_det='',i10,es10.2)') n_det, real(n_det)
        call flush(6)
       !call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue)
        if (present(initial_vector)) then
          call davidson_sparse_single(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,initial_vector)
         !call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
         !call cyrus_diagonalize_sparse(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
        else
          call davidson_sparse_single(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue)
         !call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue)
         !call cyrus_diagonalize_sparse(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue)
        endif
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
          if (i.eq.j) then
             if (sym_on) then
                call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
             else
                call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), 0, matrix_element)
             endif
             ham(i,j)=matrix_element
          else
             if (sym_on) then
                call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                ham(i,j)=matrix_element
                ham(j,i)=matrix_element
             else
                call excitation_level(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level)
                if (excite_level >= 0) then
                    call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, matrix_element)
                    ham(i,j)=matrix_element
                    ham(j,i)=matrix_element
                endif
             endif
          endif
        enddo
      enddo
      if (n_det<20) then
        write (6,*) "Hamiltonian in space of determinants I'm diagonalizing:"
        do i=1,n_det
          write (6,*) ham(i,1:n_det)
        enddo
        call flush(6)
      endif
      !diagonalize with lapack
      !write (6,*) "Ham"
      !call print_real_matrix(n_det,n_det,ham)
      call real_symmetric_diagonalize_ow_ham(n_det,ham,eigenvalues)
      lowest_eigenvalue=eigenvalues(1)
      if (n_det>1)  second_lowest_eigenvalue=eigenvalues(2)
      highest_eigenvalue=eigenvalues(n_det)
      lowest_eigenvector=ham(:,1)
      deallocate(ham,eigenvalues)
    endif

    if (srt) then
      call sort(lowest_eigenvector, dets_up, dets_dn) ! Sort in order of decreasing absolute coef
      write(6,'(/,''Wavefunction has'',i9,'' determinants. Dets. with largest'',i3,'' abs coefs are:'',/,''     dets_up             dets_dn        dets_up      dets_dn  excit_level cdet'')')  n_det,min(n_det,20)
      write (fmt, '(i2)') 2*num_words
      write(6,'(' // trim(fmt) // 'b20,' // trim(fmt) // 'i25,i3,f13.9)') (dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), popcnt(iand(dets_up(i),not(dets_up(1))))+popcnt(iand(dets_dn(i),not(dets_dn(1)))), lowest_eigenvector(i),i=1,min(n_det,20))
    elseif(n_det.le.100000) then
      allocate(lowest_eigenvector_srt(n_det))
      allocate(dets_up_srt(n_det))
      allocate(dets_dn_srt(n_det))
      lowest_eigenvector_srt(1:n_det)=lowest_eigenvector(1:n_det)
      dets_up_srt(1:n_det)=dets_up(1:n_det)
      dets_dn_srt(1:n_det)=dets_dn(1:n_det)
      call sort(lowest_eigenvector_srt, dets_up_srt, dets_dn_srt) ! Sort in order of decreasing absolute coef
      write(6,'(/,''Wavefunction has'',i9,'' determinants. Dets. with largest'',i3,'' abs coefs are:'',/,''     dets_up             dets_dn        dets_up      dets_dn  excit_level cdet'')')  n_det,min(n_det,20)
      write (fmt, '(i2)') 2*num_words
      write(6,'(' // trim(fmt) // 'b20,' // trim(fmt) // 'i25,i3,f13.9)') (dets_up_srt(i), dets_dn_srt(i), dets_up_srt(i), dets_dn_srt(i), popcnt(iand(dets_up_srt(i),not(dets_up_srt(1))))+popcnt(iand(dets_dn_srt(i),not(dets_dn_srt(1)))), lowest_eigenvector_srt(i),i=1,min(n_det,20))
      deallocate(lowest_eigenvector_srt) ; deallocate(dets_up_srt) ; deallocate(dets_dn_srt)
    else
      write(6,'(/,''Wavefunction has'',i9,'' determinants'')') n_det
    endif

    write(6,'(/,''Lowest, highest eigenvalue='',9f13.6)') lowest_eigenvalue, highest_eigenvalue

! Calculate number of dets with various excitation levels and the sum of their squared wts.
    call print_excitation_levels_and_wts(n_det,dets_up,dets_dn,lowest_eigenvector,norb,orbital_symmetries)

! Calculate what tau_optimal_deterministic and tau_optimal_stochastic would be if the wavefn were this trial wavefn.
    tau_optimal_stochastic=1/(highest_eigenvalue-lowest_eigenvalue)
    tau_out=tau_multiplier*tau_optimal_stochastic
    if(n_det.ge.2) then
      tau_optimal_deterministic=2/(highest_eigenvalue+second_lowest_eigenvalue-2*lowest_eigenvalue)
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_deterministic, tau_optimal_stochastic='',9f10.6)') &
   &  tau_optimal_deterministic, tau_optimal_stochastic
     else
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_stochastic='',9f10.6)') tau_optimal_stochastic
    endif

    call my_second(2,'diagonalize_sparse_hamiltonian_chem')

  end subroutine diagonalize_sparse_hamiltonian_chem
  !=============================================================

  !=============================================================
  subroutine filter_dets(det_i_up, det_i_dn, valid_combination)
    !---------------------------------------------------------------------------
    ! Description : Determine if det_i has spatial symmetry given by spatial_symmetry_wf
    !               Note det_i_up, det_i_dn, are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !
    ! Created     : F. Petruzielo, 26 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
#endif
    logical , intent(out) :: valid_combination

    !local variables
    integer :: spatial_symmetry_det
    integer :: j

   !write (6,*) trim(point_group); call flush(6); stop "DONE"

!!  if (trim(point_group).eq.'dih') then
!!    stop "Exact diagonalization for d_inf_h not yet working"
!      write (6,*) "In dih"; call flush(6)
!
!      do n=1,10 !j = 1, norb
!        !if (btest(det_i_up, j - 1)) then
!            lz = 0
!            gu = 0 ! 0=g,1=u
!           !n=orbital_symmetries(j)
!            gu_tmp=mod(n-2,4)
!            gu=abs(gu-mod(gu_tmp+1,2))
!            lz=lz+((n-3)/4+1)*(1-2*((gu_tmp-1)/2))
!            write (6,*) "up j=",j,"n=",n,"lz=",lz,"gu=",gu
!        !endif
!        !if (btest(det_i_dn, j - 1)) then
!            lz = 0
!            gu = 0 ! 0=g,1=u
!           !n=orbital_symmetries(j)
!            gu_tmp=mod(n-2,4)
!            gu=abs(gu-mod(gu_tmp+1,2))
!            lz=lz+((n-3)/4+1)*(1-2*((gu_tmp-1)/2))
!            write (6,*) "dn j=",j,"n=",n,"lz=",lz,"gu=",gu
!        !endif
!      enddo
!
!      spatial_symmetry_det = 4*abs(lz)+1+gu
!      if (lz>0)  spatial_symmetry_det = spatial_symmetry_det-2
!      write (6,*) "lz=",lz,"gu=",gu,"spatial_symmetry_det=",spatial_symmetry_det
!      call flush(6)
!      stop "DONE"

!!  else

      valid_combination = .true.
      spatial_symmetry_det = 1

      do j = 1, norb
         if (btest(det_i_up, j - 1)) then
            spatial_symmetry_det = product_table(spatial_symmetry_det, orbital_symmetries(j))
         endif
         if (btest(det_i_dn, j - 1)) then
            spatial_symmetry_det = product_table(spatial_symmetry_det, orbital_symmetries(j))
         endif
      enddo

!!  endif

    if (spatial_symmetry_det /= spatial_symmetry_wf) then
       valid_combination = .false.
    endif

  end subroutine filter_dets
  !====================================================================================

  !====================================================================================
  subroutine find_connected_dets_chem(det_up, det_dn, n_connected_dets, connected_dets_up, connected_dets_dn, norb_include, matrix_elements, sym)
    !---------------------------------------------------------------------------
    ! Decription :  Return all determinants, formed from norb_include orbitals,
    !               connected to det (broken into det_up and det_dn).
    !               det cannot have higher orbitals than norb_include!
    !               Note that if n_core_orb (global) is greater than 0 then the current
    !               det does not connect to determinants with unoccupied core orbitals.
    !               Additionally, return all of the corresponding matrix elements.
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    ! Modified    : A Holmes, 1 Aug 2012. No longer allocate connected_dets, etc.
    !               Allocate these ahead of time with upper_bound_connections().
    ! Modified    : A Holmes, 4 May 2017. Optional input sym returns only
    !               connections of irrep equal to sym
    !---------------------------------------------------------------------------
    implicit none

    !dummy argument
    integer, intent(out) :: n_connected_dets
    integer, intent(in) :: norb_include
    real(rk), intent(out), optional :: matrix_elements(:)
    integer, intent(in), optional :: sym
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_up
    type(ik_vec), intent(in) :: det_dn
    type(ik_vec), intent(out) :: connected_dets_up(:), connected_dets_dn(:)
    type(ik_vec) :: temp_det_up, temp_det_dn
    type(ik_vec), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#else
    integer(ik), intent(in) :: det_up
    integer(ik), intent(in) :: det_dn
    integer(ik), intent(out) :: connected_dets_up(:), connected_dets_dn(:)
    integer(ik) :: temp_det_up, temp_det_dn
    integer(ik), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#endif

    !local variables
    integer :: n_double_excite_det
    integer :: count_empty_dn, count_filled_dn, count_empty_up, count_filled_up
    integer :: i, j, k, l
    real(rk)::  norm_inv
    integer,allocatable      :: iorder(:),temp_i_2(:)
    real(rk)                 :: tmp_element
    logical :: keep_det

    if (norb_include.ne.norb) then
   !if (norb_include.ne.norb-n_core_orb) then
      write (6,*) "find_connected_dets called with norb_include =",norb_include
      call flush(6)
    endif

    !if a element of array isn't set below then -1 will cause errors. This is a good safety check.
    filled_up = -1
    filled_dn = -1
    empty_up = -1
    empty_dn = -1

    !determine where up electrons are and are not
    count_filled_up = 0
    count_empty_up = 0
    do i = 1, norb_include
       if (btest(det_up, i - 1)) then
          count_filled_up = count_filled_up + 1
          filled_up(count_filled_up) = i - 1
       else
          count_empty_up = count_empty_up + 1
          empty_up(count_empty_up) = i - 1  !which bit
       endif
    enddo

    !handle case where up determinant has occupied orbitals greater than norb_include
    do i = norb_include + 1, norb
       if (btest(det_up, i - 1)) then
          count_filled_up = count_filled_up + 1
          filled_up(count_filled_up) = i - 1
       endif
    enddo

    !determine where dn electrons are and are not
    count_filled_dn = 0
    count_empty_dn = 0
    do i = 1, norb_include
       if (btest(det_dn, i - 1)) then
          count_filled_dn = count_filled_dn + 1
          filled_dn(count_filled_dn) = i - 1
       else
          count_empty_dn = count_empty_dn + 1
          empty_dn(count_empty_dn) = i - 1  !which bit
       endif
    enddo

    !handle case where dn determinant has occupied orbitals greater than norb_include
    do i = norb_include + 1, norb
       if (btest(det_dn, i - 1)) then
          count_filled_dn = count_filled_dn + 1
          filled_dn(count_filled_dn) = i - 1
       endif
    enddo

    do i=1,n_connected_dets
      connected_dets_up(i) = 0
      connected_dets_dn(i) = 0
    enddo
    n_connected_dets = 1
    connected_dets_up(n_connected_dets) = det_up
    connected_dets_dn(n_connected_dets) = det_dn

    !generate all spin 1 double excitations
    do i = 1 + n_core_orb, nup - 1
       do j = i + 1, nup
          do k = 1, count_empty_up - 1
             do l = k+1, count_empty_up
                temp_det_up = ibclr(det_up,filled_up(i))
                temp_det_up = ibclr(temp_det_up,filled_up(j))
                temp_det_up = ibset(temp_det_up,empty_up(k))
                temp_det_up = ibset(temp_det_up,empty_up(l))
                keep_det = .false.
                if (present(sym)) then
                  if (det_sym(temp_det_up,det_dn)==sym)  keep_det = .true.
                else
                  if (product_table(orbital_symmetries(filled_up(i)+1), orbital_symmetries(filled_up(j)+1)) == product_table(orbital_symmetries(empty_up(k)+1), orbital_symmetries(empty_up(l)+1)))  keep_det = .true.
                endif
                if (keep_det) then
#ifdef NUM_ORBITALS_GT_127
                   if (temp_det_up .gt. maskr_vec(norb_include) .or. det_dn .gt. maskr_vec(norb_include)) cycle !det cannot have higher orbitals than norb_include!
#else
                   if (temp_det_up .ge. 2_ik**norb_include .or. det_dn .ge. 2_ik**norb_include) cycle !det cannot have higher orbitals than norb_include!
#endif
                   n_connected_dets = n_connected_dets + 1
                   connected_dets_up(n_connected_dets) = temp_det_up
                   connected_dets_dn(n_connected_dets) = det_dn
                endif
             enddo
          enddo
       enddo
    enddo

    !generate all spin -1 double excitations
    do i = 1 + n_core_orb, ndn - 1
       do j = i + 1, ndn
          do k = 1, count_empty_dn - 1
             do l = k+1, count_empty_dn
                temp_det_dn = ibclr(det_dn,filled_dn(i))
                temp_det_dn = ibclr(temp_det_dn,filled_dn(j))
                temp_det_dn = ibset(temp_det_dn,empty_dn(k))
                temp_det_dn = ibset(temp_det_dn,empty_dn(l))
                keep_det = .false.
                if (present(sym)) then
                  if (det_sym(det_up,temp_det_dn)==sym)  keep_det = .true.
                else
                  if (product_table(orbital_symmetries(filled_dn(i)+1), orbital_symmetries(filled_dn(j)+1)) == product_table(orbital_symmetries(empty_dn(k)+1), orbital_symmetries(empty_dn(l)+1)))  keep_det = .true.
                endif
                if (keep_det) then
#ifdef NUM_ORBITALS_GT_127
                   if (temp_det_up .gt. maskr_vec(norb_include) .or. det_dn .gt. maskr_vec(norb_include)) cycle !det cannot have higher orbitals than norb_include!
#else
                   if (det_up .ge. 2_ik**norb_include .or. temp_det_dn .ge. 2_ik**norb_include) cycle
#endif
                   n_connected_dets = n_connected_dets + 1
                   connected_dets_up(n_connected_dets) = det_up
                   connected_dets_dn(n_connected_dets) = temp_det_dn
                endif
             enddo
          enddo
       enddo
    enddo

    !generate all spin 0 double excitations
    do i = 1 + n_core_orb, nup
       do j = 1 + n_core_orb, ndn
          do k = 1, count_empty_up
             do l = 1, count_empty_dn
                temp_det_up = ibclr(det_up,filled_up(i))
                temp_det_up = ibset(temp_det_up,empty_up(k))
                temp_det_dn = ibclr(det_dn,filled_dn(j))
                temp_det_dn = ibset(temp_det_dn,empty_dn(l))
                keep_det = .false.
                if (present(sym)) then
                  if (det_sym(temp_det_up,temp_det_dn)==sym)  keep_det = .true.
                else
                  if (product_table(orbital_symmetries(filled_up(i)+1), orbital_symmetries(filled_dn(j)+1)) == product_table(orbital_symmetries(empty_up(k)+1), orbital_symmetries(empty_dn(l)+1)))  keep_det = .true.
                endif
                if (keep_det) then
#ifdef NUM_ORBITALS_GT_127
                   if (temp_det_up .gt. maskr_vec(norb_include) .or. det_dn .gt. maskr_vec(norb_include)) cycle !det cannot have higher orbitals than norb_include!
#else
                   if (temp_det_up .ge. 2_ik**norb_include .or. temp_det_dn .ge. 2_ik**norb_include) cycle
#endif
                   n_connected_dets = n_connected_dets + 1
                   connected_dets_up(n_connected_dets) = temp_det_up
                   connected_dets_dn(n_connected_dets) = temp_det_dn
                endif
             enddo
          enddo
       enddo
    enddo
    n_double_excite_det = n_connected_dets - 1

    !generate all single excitations
    do i = 1 + n_core_orb, nup
       do k = 1, count_empty_up
          temp_det_up = ibclr(det_up,filled_up(i))
          temp_det_up = ibset(temp_det_up,empty_up(k))
          keep_det = .false.
          if (present(sym)) then
            if (det_sym(temp_det_up,det_dn)==sym)  keep_det = .true.
          else
            if (orbital_symmetries(filled_up(i)+1) == orbital_symmetries(empty_up(k)+1))  keep_det = .true.
          endif
          if (keep_det) then
#ifdef NUM_ORBITALS_GT_127
             if (temp_det_up .gt. maskr_vec(norb_include) .or. det_dn .gt. maskr_vec(norb_include)) cycle !det cannot have higher orbitals than norb_include!
#else
             if (temp_det_up .ge. 2_ik**norb_include .or. det_dn .ge. 2_ik**norb_include) cycle
#endif
             n_connected_dets = n_connected_dets + 1
             connected_dets_up(n_connected_dets) = temp_det_up
             connected_dets_dn(n_connected_dets) = det_dn
          endif
       enddo
    enddo

    do i = 1 + n_core_orb, ndn
       do k = 1, count_empty_dn
          temp_det_dn = ibclr(det_dn,filled_dn(i))
          temp_det_dn = ibset(temp_det_dn,empty_dn(k))
          keep_det = .false.
          if (present(sym)) then
            if (det_sym(det_up,temp_det_dn)==sym)  keep_det = .true.
          else
            if (orbital_symmetries(filled_dn(i)+1) == orbital_symmetries(empty_dn(k)+1))  keep_det = .true.
          endif
          if (keep_det) then
#ifdef NUM_ORBITALS_GT_127
             if (temp_det_up .gt. maskr_vec(norb_include) .or. det_dn .gt. maskr_vec(norb_include)) cycle !det cannot have higher orbitals than norb_include!
#else
             if (det_up .ge. 2_ik**norb_include .or. temp_det_dn .ge. 2_ik**norb_include) cycle
#endif
             n_connected_dets = n_connected_dets + 1
             connected_dets_up(n_connected_dets) = det_up
             connected_dets_dn(n_connected_dets) = temp_det_dn
          endif
       enddo
    enddo

    !cleanup
   !allocate(temp_dets(n_connected_dets))

   !temp_dets = connected_dets_up(1:n_connected_dets)
   !deallocate(connected_dets_up)
   !allocate(connected_dets_up(n_connected_dets))
   !connected_dets_up = temp_dets

   !temp_dets = connected_dets_dn(1:n_connected_dets)
   !deallocate(connected_dets_dn)
   !allocate(connected_dets_dn(n_connected_dets))
   !connected_dets_dn = temp_dets

   !deallocate(temp_dets)

    !calculate matrix elements for the connected determinants
    if (present(matrix_elements)) then
      !allocate(matrix_elements(n_connected_dets))
       matrix_elements = 0._rk
       !diagonal element
       call hamiltonian_chem(det_up, det_dn, connected_dets_up(1), connected_dets_dn(1), 0, matrix_elements(1))
       do i = 2, n_double_excite_det + 1
          !double excitation
          call hamiltonian_chem(det_up, det_dn, connected_dets_up(i), connected_dets_dn(i), 2, matrix_elements(i))
       enddo
       do i = n_double_excite_det + 2, n_connected_dets
          !single excitation
          call hamiltonian_chem(det_up, det_dn, connected_dets_up(i), connected_dets_dn(i), 1, matrix_elements(i))
       enddo
    endif

    ! If time symmetry is turned on.......
    ! HJC comments : We have a list of connected determinants and matrix elements upto this point. We just have to convert all determinants to their representatives and add up or subtact out matrix elements depending on z. We have to approporiately weigh the matrix elements depending on whether det_up=det_dn or not. Convention for representative is det_up<det_dn
    ! As a convention the first entry MUST be the diagonal entry since find doubly excited makes this assumption in calculation of numerator and denominator of energy_pieces to store
    if (time_sym) then
        !write (6,*) "Time symmetry on"
        !call flush(6)
        norm_inv=1._rk
        allocate(iorder(n_connected_dets))
        do j=1,n_connected_dets
           iorder(j)=j
        enddo
        allocate(temp_i16_up((n_connected_dets+1)/2))
        allocate(temp_i16_dn((n_connected_dets+1)/2))
        allocate(temp_i_2((n_connected_dets+1)/2))
        !write (6,*) "n_connected_dets = ",n_connected_dets
        !write (6,*) "Before time reversal symmetry imposed"
        !write (6,'(''n connected dets before ='' i10)') n_connected_dets
        if (connected_dets_up(1) .eq. connected_dets_dn(1)) norm_inv=sqrt2inv
        do i=1,n_connected_dets
            !write(6,'(2i15,f15.8)') connected_dets_up(i),connected_dets_dn(i),matrix_elements(i)
            if (connected_dets_up(i)>connected_dets_dn(i)) then         ! If up is greater than down then swap them
                temp_det_up=connected_dets_up(i)
                connected_dets_up(i)=connected_dets_dn(i)
                connected_dets_dn(i)=temp_det_up
                if (present(matrix_elements)) matrix_elements(i)=z*matrix_elements(i)
            endif
            if (present(matrix_elements)) then
                if (connected_dets_up(i) .eq. connected_dets_dn(i)) matrix_elements(i)=sqrt2*matrix_elements(i)
                matrix_elements(i)=norm_inv*matrix_elements(i)
            endif
        enddo
        call merge_sort2_up_dn(connected_dets_up,connected_dets_dn,iorder,n_connected_dets,temp_i16_up,temp_i16_dn,temp_i_2)
        if (present(matrix_elements)) then
            matrix_elements(1:n_connected_dets)=matrix_elements(iorder(1:n_connected_dets))
            call merge_original_with_spawned3(connected_dets_up,connected_dets_dn,matrix_elements,n_connected_dets) ! Why do we need to merge here? -AAH
        else
            call merge_original_with_spawned3(n_connected_dets,connected_dets_up,connected_dets_dn)
        endif

        deallocate(iorder,temp_i16_up,temp_i16_dn,temp_i_2)
        !Bring diagonal to top
        if ((connected_dets_up(1) .ne.det_up) .or. (connected_dets_dn(1) .ne. det_dn)) then
            do i=2,n_connected_dets
                if ((connected_dets_up(i) .eq. det_up) .and. (connected_dets_dn(i) .eq. det_dn)) then
                     connected_dets_up(i)=connected_dets_up(1)
                     connected_dets_dn(i)=connected_dets_dn(1)
                     connected_dets_up(1)=det_up
                     connected_dets_dn(1)=det_dn
                     if (present(matrix_elements)) then
                        tmp_element=matrix_elements(i)
                        matrix_elements(i)=matrix_elements(1)
                        matrix_elements(1)=tmp_element
                     endif
                     exit
                endif
            enddo
        endif

        !write (6,*)
        !write (6,*) "After time reversal symmetry imposed"
        !write (6,'(''n connected dets after ='' i10)') n_connected_dets
        !do i=1,n_connected_dets
        !    write(6,'(2i15,f15.8)') connected_dets_up(i),connected_dets_dn(i),matrix_elements(i)
        !enddo
        !call flush(6)
        !stop
    endif

  end subroutine find_connected_dets_chem
  !==============================================================================================================


  subroutine find_important_connected_dets_chem(det_up, det_dn, eps, n_connected_dets, connected_dets_up, connected_dets_dn, matrix_elements, ref_diag_elem, diag_elems_info, eps_big, matrix_elements_big, min_H_already_done_elem, core_up, core_dn, virt_up, virt_dn, active_only)
  ! Same as find_connected_dets_chem, but instead of finding *all*
  ! connected dets, it finds all single excitations but only the
  ! important double excitations, i.e., only those with value
  ! greater than eps
  !
  ! Requires dtm_hb, the deterministic heatbath data structure
  !
  ! Result is unsorted
  ! 
  ! Does not contain repeats unless time-reversal symmetry is used
  ! 
  ! If time-sym is used, then apply usual heat-bath algorithm to only the
  ! representatives, so new time-reversed basis states are selected if at least
  ! one matrix element connecting it to the reference exceeds eps. When this
  ! subroutine is used for PT (and matrix elements are returned), each
  ! excitation carries its contribution to the time-symmetrized matrix element,
  ! so when they are merged and sorted, all components of the matrix elements
  ! exceeeding eps will be included
  !
  ! A Holmes, 22 Feb 2016
  ! Modified : A Holmes, 16 May 2017. Added optional arguments core_up/dn,
  !            virt_up/dn, which are masks of core or virtual orbitals (which must be
  !            either occupied or unoccupied in the connected dets list)
  !            Optional variable active_only is used as follows. If core_up,
  !            etc. are set, active_only is assumed to be true, meaning only
  !            dets in the active space are generated. If active_only is false,
  !            then only dets NOT in the active space are generated

    use common_run, only : diag_elem_info,connected_diag_elems_info
    use more_tools, only : get_occ_orbs

    implicit none

    !dummy argument
    real(rk), intent(in) :: eps
    integer, intent(out) :: n_connected_dets
    real(rk), intent(out), optional :: matrix_elements(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_up,det_dn
    type(ik_vec), intent(out) :: connected_dets_up(:), connected_dets_dn(:)
    type(ik_vec) :: new_up, new_dn, tmp_det
    type(ik_vec), optional, intent(in) :: core_up, core_dn, virt_up, virt_dn
#else
    integer(ik), intent(in) :: det_up,det_dn
    integer(ik), intent(out) :: connected_dets_up(:), connected_dets_dn(:)
    integer(ik) :: new_up, new_dn, tmp_det
    integer(ik), optional, intent(in) :: core_up, core_dn, virt_up, virt_dn
#endif
    real(rk),optional,intent(in) :: ref_diag_elem
    type(diag_elem_info),optional,intent(out) :: diag_elems_info(:)
    real(rk),optional,intent(in) :: eps_big ! If present, remove contributions to sum that exceed eps_big
    real(rk),optional,intent(out) :: matrix_elements_big(:)
    real(rk),optional,intent(in) :: min_H_already_done_elem
    logical,optional,intent(in) :: active_only

    integer :: p,q,r,s,p_el
    real(rk) :: matrix_element
    integer :: ipair,npairs
    integer :: i
    integer :: epair,hpair
    integer :: p2,q2,r_tmp
    integer :: excite_level
    logical :: active_only_local=.true.
    logical :: keep

    if (present(min_H_already_done_elem).and.present(matrix_elements)) then
      stop "Find important connected dets called with both min H and matrix elements (should not get here)!"
    endif

    if (present(active_only))  active_only_local = active_only

    ! First, add diagonal element, to be consistent with find_connected_dets_chem
    ! WARNING: Setting matrix_elements(1) to 0 since it is usually not needed for anything
    n_connected_dets = 1
    connected_dets_up(1) = det_up
    connected_dets_dn(1) = det_dn
    if (present(matrix_elements))  matrix_elements(1) = 0._rk
    if (present(matrix_elements_big))  matrix_elements_big(1) = 0._rk
    if (present(diag_elems_info))  diag_elems_info(1)%old_diag_elem = 1.e51_rk

    ! Get occ_up, occ_dn
    call get_occ_orbs(det_up,det_dn,occ_up,occ_dn)

    ! First, add single excitations
    ! This can be sped up using symmetries (TODO)
    ! If time-rev symmetry, only need one of these
    do p=1,nelec
      if (p<=nup) then
        p_el = occ_up(p)
      else
        p_el = occ_dn(p-nup)
      endif
      do r=1,norb
        if (p<=nup) then
          if (btest(det_up,r-1))  cycle
        else
          if (btest(det_dn,r-1))  cycle
        endif
        if (orbital_symmetries(p_el).ne.orbital_symmetries(r))  cycle
        if (p<=nup) then
          new_up = ibset(ibclr(det_up,p_el-1),r-1)
          new_dn = det_dn
        else
          new_up = det_up
          new_dn = ibset(ibclr(det_dn,p_el-1),r-1)
        endif
        if (.not.active_only_local)  keep = .false.
        if (present(core_up)) then
          if (active_only_local) then
            if (iand(core_up,new_up).ne.core_up)  cycle ! core up excitations
            if (iand(core_dn,new_dn).ne.core_dn)  cycle ! core dn excitations
          else
            if (iand(core_up,new_up).ne.core_up)  keep = .true.
            if (iand(core_dn,new_dn).ne.core_dn)  keep = .true.
          endif
        endif
        if (present(virt_up)) then
          if (active_only_local) then
            if (iand(virt_up,new_up).ne.0_ik)  cycle ! virtual up excitations
            if (iand(virt_dn,new_dn).ne.0_ik)  cycle ! virtual dn excitations
          else
            if (iand(virt_up,new_up).ne.0_ik)  keep = .true.
            if (iand(virt_dn,new_dn).ne.0_ik)  keep = .true.
          endif
        endif
        if (.not.active_only_local) then
          if (.not.keep)  cycle
        endif
    
        if (time_sym) then
          if (new_up==new_dn.and.z<0)  cycle
          if (det_up==new_dn.and.det_dn==new_up)  cycle ! time-reversed move is diagonal
        endif

        call hamiltonian_chem(det_up,det_dn,new_up,new_dn,1,matrix_element)

        if (abs(matrix_element)<eps)  cycle
        if (present(min_H_already_done_elem)) then
          if (abs(matrix_element) > min_H_already_done_elem)  cycle
        endif

        if (time_sym.and.present(matrix_elements)) then ! Only compute this contribution to the off-diagonal matrix element
          if (det_up==det_dn.and.new_up.ne.new_dn)  matrix_element = sqrt2inv*matrix_element
          if (new_up==new_dn.and.det_up.ne.det_dn)  matrix_element = sqrt2*matrix_element
        endif

        if (time_sym.and.new_up>new_dn) then
          tmp_det = new_up
          new_up = new_dn
          new_dn = tmp_det
          matrix_element = z*matrix_element
        endif

        n_connected_dets = n_connected_dets + 1
        connected_dets_up(n_connected_dets) = new_up
        connected_dets_dn(n_connected_dets) = new_dn
        if (present(matrix_elements))  matrix_elements(n_connected_dets) = matrix_element
        if (present(eps_big)) then
          if (abs(matrix_element)>eps_big) then
            matrix_elements_big(n_connected_dets) = matrix_element
          else
            matrix_elements_big(n_connected_dets) = 0._rk
          endif
        endif
      enddo
    enddo
    
    ! Initialize diagonal elements corresponding to single excitations to "unknown" valuefor now (even though computing new diagonal element for single excitations is easy; I just haven't done it yet) (TODO)
    if (present(ref_diag_elem)) then
      do i=1,n_connected_dets
        diag_elems_info(i)%old_diag_elem = 1.e51_rk
      enddo
    endif
  
    ! Next, add double excitations
    if (eps>max_double)  return
  
    npairs = 0
    ! Same spin, up
    do p=1,nup
      do q=p+1,nup
        npairs = npairs + 1
        pairs_e1(npairs) = occ_up(p)
        pairs_e2(npairs) = occ_up(q)
      enddo
    enddo
    ! Same spin, dn
    do p=1,ndn
      do q=p+1,ndn
        npairs = npairs + 1
        pairs_e1(npairs) = occ_dn(p)+norb
        pairs_e2(npairs) = occ_dn(q)+norb
      enddo
    enddo
    ! Opposite spin
    do p=1,nup
      do q=1,ndn
        npairs = npairs + 1
        pairs_e1(npairs) = occ_up(p)
        pairs_e2(npairs) = occ_dn(q)+norb
      enddo
    enddo
  
    do ipair=1,npairs
      p = pairs_e1(ipair)
      q = pairs_e2(ipair)
  
      p2 = p
      q2 = q
      if (p>norb.and.q>norb) then
        p2 = p-norb
        q2 = q-norb
      endif
      if (p<=norb.and.q>norb.and.p>q-norb) then
        p2 = q-norb
        q2 = p+norb
      endif
  
      epair = int(combine_2_indices(p2,q2))

      do hpair=1,pq_count(epair)

        if (dtm_hb(pq_ind(epair)+hpair-1)%absH<=eps)  exit

        if (present(min_H_already_done_elem)) then
          if (dtm_hb(pq_ind(epair)+hpair-1)%absH > min_H_already_done_elem)  cycle
        endif

        r = dtm_hb(pq_ind(epair)+hpair-1)%r
        s = dtm_hb(pq_ind(epair)+hpair-1)%s
        if (p>norb.and.q>norb) then
          r = r+norb
          s = s+norb
        endif
        if (p<=norb.and.q>norb.and.p>q-norb) then
          r_tmp = s-norb
          s = r+norb
          r = r_tmp
        endif
        if (is_occupied(det_up,det_dn,r))  cycle
        if (is_occupied(det_up,det_dn,s))  cycle
  
        ! Generate new determinants
        new_up = det_up
        new_dn = det_dn
        if (p<=norb) then
          new_up = ibclr(new_up,p-1)
        else
          new_dn = ibclr(new_dn,p-norb-1)
        endif
        if (q<=norb) then
          new_up = ibclr(new_up,q-1)
        else
          new_dn = ibclr(new_dn,q-norb-1)
        endif
        if (r<=norb) then
          new_up = ibset(new_up,r-1)
        else
          new_dn = ibset(new_dn,r-norb-1)
        endif
        if (s<=norb) then
          new_up = ibset(new_up,s-1)
        else
          new_dn = ibset(new_dn,s-norb-1)
        endif
  
        ! Check whether excitation is in active space
        if (.not.active_only_local)  keep = .false.
        if (present(core_up)) then
          if (active_only_local) then
            if (iand(core_up,new_up).ne.core_up)  cycle ! core up excitations
            if (iand(core_dn,new_dn).ne.core_dn)  cycle ! core dn excitations
          else
            if (iand(core_up,new_up).ne.core_up)  keep = .true.
            if (iand(core_dn,new_dn).ne.core_dn)  keep = .true.
          endif
        endif
        if (present(virt_up)) then
          if (active_only_local) then
            if (iand(virt_up,new_up).ne.0_ik)  cycle ! virtual up excitations
            if (iand(virt_dn,new_dn).ne.0_ik)  cycle ! virtual dn excitations
          else
            if (iand(virt_up,new_up).ne.0_ik)  keep = .true.
            if (iand(virt_dn,new_dn).ne.0_ik)  keep = .true.
          endif
        endif
        if (.not.active_only_local) then
          if (.not.keep)  cycle
        endif

        if (time_sym) then
          if (new_up==new_dn.and.z<0)  cycle
          if (det_up==new_dn.and.det_dn==new_up)  cycle ! time-reversed move is diagonal
        endif
  
        ! Get matrix element
        if (present(matrix_elements)) then
         !matrix_element = dtm_hb(epair,hpair)%absH ! can't just do this, because sign matters
          call hamiltonian_chem(det_up,det_dn,new_up,new_dn,2,matrix_element)
          if (time_sym) then
            if (det_up==det_dn.and.new_up.ne.new_dn)  matrix_element = sqrt2inv*matrix_element
            if (new_up==new_dn.and.det_up.ne.det_dn)  matrix_element = sqrt2*matrix_element
          endif
        endif

        if (time_sym) then
          if (new_up>new_dn) then
            tmp_det = new_up
            new_up = new_dn
            new_dn = tmp_det
            matrix_element = z*matrix_element
          endif
        endif
  
        n_connected_dets = n_connected_dets + 1
        connected_dets_up(n_connected_dets)=new_up
        connected_dets_dn(n_connected_dets)=new_dn
        if (present(matrix_elements))  matrix_elements(n_connected_dets) = matrix_element
        if (present(eps_big)) then
          if (abs(matrix_element)>eps_big) then
            matrix_elements_big(n_connected_dets) = matrix_element
          else
            matrix_elements_big(n_connected_dets) = 0._rk
          endif
        endif
  
        ! Save info for computing new diagonal element in O(N) time
        if (present(ref_diag_elem)) then
          diag_elems_info(n_connected_dets)%old_diag_elem = ref_diag_elem
          diag_elems_info(n_connected_dets)%p = p
          diag_elems_info(n_connected_dets)%q = q
          diag_elems_info(n_connected_dets)%r = r
          diag_elems_info(n_connected_dets)%s = s
        endif

      enddo ! hpair (r,s)

    enddo ! ipair (p,q)

  end subroutine find_important_connected_dets_chem

  !===============================================================================================================
  subroutine excitation_level(det_i_up, det_i_dn, det_j_up, det_j_dn, excite_level)
    !--------------------------------------------------------------------------------------------------------------
    ! Description : If not generating_psi_t  then
    !               determine if det_i is connected to det_j by the hamiltonian.
    !               If det_i is connected to det_j return the excite_level
    !               relating them. Otherwise return -1 if not connected.
    !               Note this is called only with determinants of the correct spatial symmetry.
    !               If generating_psi_t then determine if det_i is maximally excited with respect to det_j.
    !               The routine is called with det_i_up as some excited det, det_i_dn as HF, and det_j as HF.
    !               This is strictly for improving efficiency in generating of psi_t
    !
    ! Created     : F. Petruzielo, 22 Nov 2010
    ! Modified    : A Holmes, 13 Oct 2015. Simplified and sped up, using same method as in count_excitations in tools.f90
    !-------------------------------------------------------------------------------------------------------------

    !dummy arguments
    integer, intent(out) :: excite_level
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(in) :: det_j_up
    type(ik_vec), intent(in) :: det_j_dn
    type(ik_vec) :: tmp
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(in) :: det_j_up
    integer(ik), intent(in) :: det_j_dn
    integer(ik) :: tmp
#endif

    excite_level = 0

    tmp = iand(det_i_up,not(det_j_up))
    do while(tmp.ne.0)
      excite_level = excite_level + 1
      if (excite_level>2) then
        excite_level = -1
        return
      endif
      tmp = iand(tmp,tmp-1)
    enddo

    tmp = iand(det_i_dn,not(det_j_dn))
    do while(tmp.ne.0)
      excite_level = excite_level + 1
      if (excite_level>2) then
        excite_level = -1
        return
      endif
      tmp = iand(tmp,tmp-1)
    enddo

  !excite_level = popcnt(iand(det_i_up,not(det_j_up)))
  !if (excite_level>2) then
  !  excite_level = -1
  !  return
  !endif

  !excite_level = excite_level + popcnt(iand(det_i_dn,not(det_j_dn)))
  !if (excite_level>2) then
  !  excite_level = -1
  !  return
  !endif

  end subroutine excitation_level
  !============================================================================================


  !=============================================================
  subroutine init_point_group()
    !---------------------------------------------------------------------------
    ! Description : Determine if det_i is connected to det_j by the hamiltonian.
    !               If det_i is connected to det_j return the excite_level
    !               relating them. Otherwise return -1 if not connected.
    !               Note this is called only with determinants of the correct spatial symmetry.
    !
    ! Created     : F. Petruzielo, 26 Nov 2010
    !---------------------------------------------------------------------------

    !dummy arguments

    !local variables
    integer :: a,b,i,old,new,lz,gu

    d_infinity_h = 0

    select case (trim(point_group))
    case ('c1')
       n_group_elements = 1
       allocate(product_table_elems(n_group_elements, n_group_elements))
       product_table_elems(1,:) = (/1 /)
    case ('cs')
       n_group_elements = 2
       allocate(product_table_elems(n_group_elements, n_group_elements))
       product_table_elems(1,:) = (/1, 2 /)
       product_table_elems(2,:) = (/2, 1 /)
    case ('c2v')
       n_group_elements = 4
       allocate(product_table_elems(n_group_elements, n_group_elements))
       product_table_elems(1,:) = (/1, 2, 3, 4 /)
       product_table_elems(2,:) = (/2, 1, 4, 3 /)
       product_table_elems(3,:) = (/3, 4, 1, 2 /)
       product_table_elems(4,:) = (/4, 3, 2, 1 /)
    case ('c2h')
       n_group_elements = 4
       allocate(product_table_elems(n_group_elements, n_group_elements))
       product_table_elems(1,:) = (/1, 2, 3, 4 /)
       product_table_elems(2,:) = (/2, 1, 4, 3 /)
       product_table_elems(3,:) = (/3, 4, 1, 2 /)
       product_table_elems(4,:) = (/4, 3, 2, 1 /)
    case ('d2h')
       n_group_elements = 8
       allocate(product_table_elems(n_group_elements, n_group_elements))
       product_table_elems(1,:) = (/1, 2, 3, 4, 5, 6, 7, 8 /)
       product_table_elems(2,:) = (/2, 1, 4, 3, 6, 5, 8, 7 /)
       product_table_elems(3,:) = (/3, 4, 1, 2, 7, 8, 5, 6 /)
       product_table_elems(4,:) = (/4, 3, 2, 1, 8, 7, 6, 5 /)
       product_table_elems(5,:) = (/5, 6, 7, 8, 1, 2, 3, 4 /)
       product_table_elems(6,:) = (/6, 5, 8, 7, 2, 1, 4, 3 /)
       product_table_elems(7,:) = (/7, 8, 5, 6, 3, 4, 1, 2 /)
       product_table_elems(8,:) = (/8, 7, 6, 5, 4, 3, 2, 1 /)
    case ('dih') ! d_infinity_h, with following index correspondence:

       ! Sandeep's indices:  1, 2, 5, 6,-5,-6, 7, 8,-7, -8,  9, 10, -9,-10, 11, 12,-11,-12
       ! Adam's indices:     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
       ! angular momentum:   0, 0, 1, 1,-1,-1, 2, 2,-2, -2,  3,  3, -3, -3,  4,  4, -4, -4
       ! g/u (g=0,u=1):      0, 1, 0, 1, 0, 1, 0, 1, 0,  1,  0,  1,  0,  1,  0,  1,  0,  1

       if (minval(orbital_symmetries)<0) then ! there are negative indices, so assume that the indices given are Sandeep's. This converts them to my indices.
         do i=1,norb
           if (orbital_symmetries(i).ne.1.and.orbital_symmetries(i).ne.2) then
             old = orbital_symmetries(i)
             a = abs(old)/2
             b = (abs(old)+1)/2
             new = a+3*b-8
             if (old<0)  new=new+2
             orbital_symmetries(i) = new
           endif
         enddo
       endif

       call get_lz(maxval(orbital_symmetries),lz,gu)
       !Max number of group elements is based off of angular momentum selection rules
       !We must take 3 * maximum angular momentum, as that is the maximum angular momentum
       !possible, and 4 * that, as there are 4 indices per angular momentum,
       !and an additional 2 for the 0 lz state.

       n_group_elements=12*abs(lz)+2

       d_infinity_h = 1

       ! Following product_table is implemented in the code in a way that allows for arbitrarily
       ! large L_z:

      !n_group_elements = 18
      !allocate(product_table_elems(n_group_elements, n_group_elements))
      !product_table_elems(1,:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18/)
      !product_table_elems(2,:) = (/ 2, 1, 4, 3, 6, 5, 8, 7,10, 9,12,11,14,13,16,15,18,17/)
      !product_table_elems(3,:) = (/ 3, 4, 7, 8, 1, 2,11,12, 5, 6,15,16, 9,10,19,20,13,14/)
      !product_table_elems(4,:) = (/ 4, 3, 8, 7, 2, 1,12,11, 6, 5,16,15,10, 9,20,19,14,13/)
      !product_table_elems(5,:) = (/ 5, 6, 1, 2, 9,10, 3, 4,13,14, 7, 8,17,18,11,12,21,22/)
      !product_table_elems(6,:) = (/ 6, 5, 2, 1,10, 9, 4, 3,14,13, 8, 7,18,17,12,11,22,21/)
      !product_table_elems(7,:) = (/ 7, 8,11,12, 3, 4,15,16, 1, 2,19,20, 5, 6,23,24, 9,10/)
      !product_table_elems(8,:) = (/ 8, 7,12,11, 4, 3,16,15, 2, 1,20,19, 6, 5,24,23,10, 9/)
      !product_table_elems(9,:) = (/ 9,10, 5, 6,13,14, 1, 2,17,18, 3, 4,21,22, 7, 8,25,26/)
      !product_table_elems(10,:)= (/10, 9, 6, 5,14,13, 2, 1,18,17, 4, 3,22,21, 8, 7,26,25/)
      !product_table_elems(11,:)= (/11,12,15,16, 7, 8,19,20, 3, 4,23,24, 1, 2,27,28, 5, 6/)
      !product_table_elems(12,:)= (/12,11,16,15, 8, 7,20,19, 4, 3,24,23, 2, 1,28,27, 6, 5/)
      !product_table_elems(13,:)= (/13,14, 9,10,17,18, 5, 6,21,22, 1, 2,25,26, 3, 4,29,30/)
      !product_table_elems(14,:)= (/14,13,10, 9,18,17, 6, 5,22,21, 2, 1,26,25, 4, 3,30,29/)
      !product_table_elems(15,:)= (/15,16,19,20,11,12,23,24, 7, 8,27,28, 3, 4,31,32, 1, 2/)
      !product_table_elems(16,:)= (/16,15,20,19,12,11,24,23, 8, 7,28,27, 4, 3,32,31, 2, 1/)
      !product_table_elems(17,:)= (/17,18,13,14,21,22, 9,10,25,26, 5, 6,29,30, 1, 2,33,34/)
      !product_table_elems(18,:)= (/18,17,14,13,22,21,10, 9,26,25, 6, 5,30,29, 2, 1,34,33/)

    case default
       write(6,*) "The point group ", point_group, " is not implemented"
       stop "The point group is not implemented"
    end select

  end subroutine init_point_group
  !===========================================================================
  !=========================================================================================
  subroutine symmetry_reduce_and_replace_chem(imp_up,imp_dn)
  !-----------------------------------------------------------------------------------------
  ! Purpose : Reduce and replace determinant list by accounting for symmetrized determinants
  ! Created : Hitesh Changlani, April 13 2012
  !-----------------------------------------------------------------------------------------

  implicit none
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable,intent(inout) :: imp_up(:),imp_dn(:)
  type(ik_vec)                           :: tmp_rep_up(size(imp_up)),tmp_rep_dn(size(imp_up))
#else
  integer(ik),allocatable,intent(inout) :: imp_up(:),imp_dn(:)
  integer(ik)                           :: tmp_rep_up(size(imp_up)),tmp_rep_dn(size(imp_up))
#endif

  ! Local
  integer                               :: nimp,nsym_imp
  integer                               :: i

  nimp=size(imp_up)
  nsym_imp=0

  do i=1,nimp
    if ((imp_up(i) .le. imp_dn(i))) then
        nsym_imp=nsym_imp+1
        tmp_rep_up(nsym_imp)=imp_up(i)
        tmp_rep_dn(nsym_imp)=imp_dn(i)
    endif
  enddo

  deallocate(imp_up)
  deallocate(imp_dn)
  allocate(imp_up(nsym_imp))
  allocate(imp_dn(nsym_imp))

  do i=1,nsym_imp
   imp_up(i)=tmp_rep_up(i)
   imp_dn(i)=tmp_rep_dn(i)
  enddo

  end subroutine symmetry_reduce_and_replace_chem

  !=====================================================================================================
  subroutine generate_sparse_ham_chem(dets_up,dets_dn,H_indices,H_nonzero_elements,H_values,time_sym_on,max_nonzero_elems,num_nonzero_elems)
    !---------------------------------------------------------------------------
    ! Description : Compute the Hamiltonian in a given space and store it sparsely.
    !
    ! Created     : A Holmes, 7 Apr 2012
    ! Modified    : A Holmes, 15 Aug 2012. If max_nonzero_elems and num_nonzero_elems are present,
    !               then count the number of nonzero elements in H and store it as num_nonzero_elems.
    !               If that number exceeds max_nonzero_elems, return before constructing H.
    !---------------------------------------------------------------------------

      use types, only : i8b
      use more_tools, only : binary_search
      use mpi_routines, only : ncores

      implicit none

      ! Dummy
#ifdef NUM_ORBITALS_GT_127
      type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
#else
      integer(ik),intent(in) :: dets_up(:),dets_dn(:)
#endif
      integer(i8b),allocatable,intent(out) :: H_indices(:),H_nonzero_elements(:)
      real(rk),allocatable,intent(out) :: H_values(:)
      logical,intent(in)               :: time_sym_on
      integer(i8b),optional,intent(in) :: max_nonzero_elems
      integer(i8b),optional,intent(out) :: num_nonzero_elems

      ! Local
      integer :: i,j,k
      integer :: excite_level,n_det
      integer(i8b) :: nonzero_elements,isparse,itotal
      real(rk) :: matrix_element
      integer,allocatable :: H_rows(:),row_count(:)
      logical             :: sym_on
      integer :: n_connected_dets
      integer :: connected_indices(max_connected_dets)
      integer :: numops
      logical :: brute_force ! whether or not to use the brute force method (instead of binary search)

      call my_second(1, 'generate_sparse_ham_chem')

      sym_on=.false.
      !if (present(time_sym_on) .and. time_sym) sym_on=time_sym_on
      if (time_sym_on .and. time_sym) sym_on=time_sym_on

      n_det=size(dets_up)

      allocate(H_nonzero_elements(n_det))
      H_nonzero_elements(:)=0

      allocate(row_count(n_det))

      ! Count number of nonzero elements of H
      nonzero_elements=0

      call find_connected_dets_chem(dets_up(1), dets_dn(1), n_connected_dets, connected_dets_up, connected_dets_dn, norb)

      ! Determine which method to use based on n_det and n_connected_dets (in general, larger matrices are faster with binary search)
      if (real(n_det)>2.1*real(n_connected_dets)*log(real(n_det))) then
        brute_force = .false.
      else
        brute_force = .true.
      endif

      if (brute_force) then

        do i = 1, n_det
           if (n_det>10) then
             if (mod(i,n_det/10)==0.or.i==n_det/100.or.i==(3*n_det)/100.or.i==(5*n_det)/100) then
               write (6,*) (100*i)/n_det,"% done"
             endif
           endif
           do j = 1, i
              matrix_element=0._rk
              if (i.eq.j) then
                 if (sym_on) then
                      call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                 else
                      call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), 0, matrix_element)
                 endif
                 if (abs(matrix_element).gt.1.d-16) then
                   nonzero_elements = nonzero_elements+1
                   H_nonzero_elements(i)=H_nonzero_elements(i)+1
                 endif
              else
                 if (sym_on) then
                      call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                 else
                     call excitation_level(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level)
                     if (excite_level >= 0) call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, matrix_element)
                 endif
                 if (abs(matrix_element).gt.1.d-16) then
                      nonzero_elements = nonzero_elements+2
                      H_nonzero_elements(i)=H_nonzero_elements(i)+1
                      H_nonzero_elements(j)=H_nonzero_elements(j)+1
                 endif
              endif
           enddo
        enddo

      else

        do i = 1, n_det
           if (n_det>10) then
             if (mod(i,n_det/10)==0.or.i==n_det/100.or.i==(3*n_det)/100.or.i==(5*n_det)/100) then
               write (6,*) (100*i)/n_det,"% done"
             endif
           endif
           call find_connected_dets_chem(dets_up(i), dets_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn, norb)
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
        num_nonzero_elems=nonzero_elements
        if (num_nonzero_elems>max_nonzero_elems)  return
      endif

      write(6,'(''allocating H_indices etc. arrays of size nonzero_elements='',i10,'' ='',es9.2,'' in generate_sparse_ham_chem'')') nonzero_elements, real(nonzero_elements)
      call flush(6)
      if (ncores>1.and.real(nonzero_elements)>7.e6_rk) then
        write(6,'(''Error: attempted to allocate too much memory ('',i10,'' elements) on each of too many cores ('',i4,''). Try generating the local energies and deterministic projector with a single-core run, then read in the elements during the multi-core run.'')') nonzero_elements, ncores
        stop "Attempted to allocate too much memory on each of too many cores!"
      endif

      call my_second(1, 'allocating H_indices etc.')

      allocate(H_indices(nonzero_elements))
      allocate(H_rows(nonzero_elements))
      allocate(H_values(nonzero_elements))
      H_indices(:)=0
      H_values(:)=0._rk

      call my_second(2, 'allocating H_indices etc.')

      ! Generate H

      isparse = 0
      itotal = 0

      if (brute_force) then

        do i = 1, n_det
           if (i>1) then
             itotal=itotal+H_nonzero_elements(i-1)
             H_rows(isparse+1:itotal)=H_rows(isparse)
             isparse = itotal
           endif
           do j = 1, i
              matrix_element=0._rk
              if (i.eq.j) then
                 if (sym_on) then
                      call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                 else
                      call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), 0, matrix_element)
                 endif
                 if (abs(matrix_element).gt.1.d-16) then
                   isparse=isparse+1
                   H_indices(isparse)=j
                   H_rows(isparse)=i
                   H_values(isparse)=matrix_element
                 endif
                 row_count(i)=int(isparse)+1 ! First location of a zero element in the i'th row
              else
                 if (sym_on) then
                     call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                 else
                     call excitation_level(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level)
                     if (excite_level >= 0) call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, matrix_element)
                 endif
                 if (abs(matrix_element).gt.1.d-16) then
                       isparse=isparse+1
                       H_indices(isparse)=j
                       H_rows(isparse)=i
                       H_values(isparse)=matrix_element
                       H_values(row_count(H_indices(isparse))) = H_values(isparse)
                       H_indices(row_count(H_indices(isparse))) = H_rows(isparse)
                       row_count(H_indices(isparse)) = row_count(H_indices(isparse)) + 1
                 endif
              endif
           enddo
        enddo

      else

        do i = 1, n_det
           if (i>1) then
             itotal=itotal+H_nonzero_elements(i-1)
             H_rows(isparse+1:itotal)=H_rows(isparse)
             isparse = itotal
           endif
           call find_connected_dets_chem(dets_up(i), dets_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn, norb)
           do k=1,n_connected_dets
              call binary_search(connected_dets_up(k),connected_dets_dn(k),dets_up(1:n_det),dets_dn(1:n_det),j)
              connected_indices(k) = j
           enddo
           ! sort connected_indices
           call sort(n_connected_dets,connected_indices(1:n_connected_dets),numops)
           do k=1,n_connected_dets
              j = connected_indices(n_connected_dets-k+1)
              if (j>0.and.i>=j) then
                   matrix_element=0._rk
                   if (sym_on) then
                        call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                   else
                      if (i==j) then
                        call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), 0, matrix_element)
                      else
                        call excitation_level(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j), excite_level)
                        if (excite_level >= 0) then
                          call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, matrix_element)
                        endif
                      endif
                   endif
                   isparse=isparse+1
                   H_indices(isparse)=j
                   H_rows(isparse)=i
                   H_values(isparse)=matrix_element
                   row_count(i)=int(isparse)+1 ! First location of a zero element in the i'th row
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

      call my_second(2, 'generate_sparse_ham_chem')

  end subroutine generate_sparse_ham_chem


  !=====================================================================================================
  subroutine generate_sparse_ham_chem_upper_triangular(dets_up,dets_dn,H_indices,H_nonzero_elements,H_values,time_sym_on,hf_to_psit,max_nonzero_elems,num_nonzero_elems,sparse_ham,average_connections)
    !---------------------------------------------------------------------------
    ! Description : Compute the Hamiltonian in a given space and store it sparsely, as an upper triangular matrix.
    !
    ! Created     : A Holmes, 17 Dec 2012
    ! Modified    : A Holmes, 11 Jan 2013. Input hf_to_psit replaces the first state (HF) with psi_trial.
    !             : A Holmes, Nov 2016. Implemented Sandeep's partial connections
    !             : A Holmes, 27 Feb 2017. Added time-reversal symmetry for partial connections
    !---------------------------------------------------------------------------

      use types, only : i8b
      use more_tools, only : binary_search
      use mpi_routines, only : ncores
      use generic_sort, only : sort_by_first_argument
      use common_selected_ci, only : sparse_mat

      implicit none

      ! Dummy
#ifdef NUM_ORBITALS_GT_127
      type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
      type(ik_vec),allocatable :: beta(:),alpha_m1(:)
      type(ik_vec),allocatable :: dets_up_time_sym(:)
#else
      integer(ik),intent(in) :: dets_up(:),dets_dn(:)
      integer(ik),allocatable :: beta(:),alpha_m1(:)
      integer(ik),allocatable :: dets_up_time_sym(:)
#endif
      integer(i8b),allocatable,intent(out) :: H_indices(:),H_nonzero_elements(:)
      real(rk),allocatable,intent(out) :: H_values(:)
      logical,intent(in)               :: time_sym_on
      integer(i8b),optional,intent(in) :: max_nonzero_elems
      integer(i8b),optional,intent(out) :: num_nonzero_elems
      logical,intent(in) :: hf_to_psit
      type(sparse_mat),optional,intent(inout) :: sparse_ham
      real(rk),optional,intent(inout) :: average_connections
      ! Local
      integer(i8b),allocatable :: tmp_H_indices(:)
      real(rk),allocatable :: tmp_H_values(:)
      integer(i8b) :: i
      integer :: k, n_unique_beta_new
      integer(i8b) :: j
      integer :: excite_level
      integer(i8b) :: n_det
      integer(i8b) :: nonzero_elements,old_nonzero_elements,isparse
      real(rk) :: matrix_element
      logical             :: sym_on
      integer :: n_connected_dets
      integer, allocatable :: connected_indices(:)
      integer :: numops
      logical :: brute_force ! whether or not to use the brute force method (instead of binary search)
      logical :: use_partial_connections ! whether or not to use Sandeep's partial connections idea
      integer,allocatable :: beta_ind(:),alpha_m1_ind(:)
      integer :: isparse_old,n_copy
      integer(i8b) :: ndet_old,ndet_new
      logical, allocatable :: is_included(:)
      integer :: ndet_target

      call my_second(1, 'generate_sparse_ham_chem_upper_triangular_serial')

      sym_on=.false.
      !if (present(time_sym_on) .and. time_sym) sym_on=time_sym_on
      if (time_sym_on .and. time_sym) sym_on=time_sym_on

      n_det=size(dets_up)
      if (time_sym) then
        allocate(connected_indices(2*n_det))
        allocate(H_nonzero_elements(2*n_det))
        H_nonzero_elements(1:2*n_det)=0
      else
        allocate(connected_indices(n_det))
        allocate(H_nonzero_elements(n_det))
        H_nonzero_elements(1:n_det)=0
      endif

     !allocate(row_count(n_det))

      ! Count number of nonzero elements of H
      nonzero_elements=0

      call find_connected_dets_chem(dets_up(1), dets_dn(1), n_connected_dets, connected_dets_up, connected_dets_dn, norb)

      use_partial_connections = .true.

      if (use_partial_connections) then
        write (6,'(''Constructing sparse Hamiltonian using partial connections'')') ; call flush(6)
      else
        ! Determine which method to use based on n_det and n_connected_dets (in general, larger matrices are faster with binary search)
        if (real(n_det)>2.1*real(n_connected_dets)*log(real(n_det))) then
          brute_force = .false.
          write (6,'(''Constructing sparse Hamiltonian by generating connections and performing binary search'')') ; call flush(6)
        else
          brute_force = .true.
          write (6,'(''Constructing sparse Hamiltonian by checking for connectedness between all (row,column) pairs'')') ; call flush(6)
        endif
      endif

      !nonzero_elements = int(n_det*real(min(n_det,n_connected_dets))**0.5,i8b) ! Just a starting guess
      !MJO Overallocate by 5% so that we never reallocate. This is useful for having reproduceable results
      !    The copy can cause small variations
      nonzero_elements  = int(n_det*(average_connections+1)*1.05,i8b) ! Just a starting guess
      if (present(max_nonzero_elems).and.present(num_nonzero_elems)) then
        num_nonzero_elems=nonzero_elements
        if (num_nonzero_elems>max_nonzero_elems)  return
      endif

      write(6,'(''allocating H_indices etc. arrays of size nonzero_elements (estimated)='',i10,'' ='',es9.2,'' in generate_sparse_ham_chem'')') nonzero_elements, real(nonzero_elements)
      call flush(6)
      if (ncores>1.and.real(nonzero_elements)>7.e6_rk) then
        write(6,'(''Error: attempted to allocate too much memory ('',i10,'' elements) on each of too many cores ('',i4,''). Try generating the local energies and deterministic projector with a single-core run, then read in the elements during the multi-core run.'')') nonzero_elements, ncores
        stop "Attempted to allocate too much memory on each of too many cores!"
      endif

      call my_second(1, 'allocating H_indices etc.')

      allocate(H_indices(nonzero_elements))
      allocate(H_values(nonzero_elements))
      H_indices(:)=0
      H_values(:)=0._rk

      call my_second(2, 'allocating H_indices etc.')

      ! Generate H

      isparse = 0

      if (use_partial_connections) then

        call my_second(1,'create beta and alpha-1 helper lists')

        if (present(sparse_ham)) then
          ndet_old = sparse_ham%ndet
          ndet_new = n_det - ndet_old
        else
          ndet_old = 0
          ndet_new = n_det
        endif

        ! First, sort a list of (dets_dn(i),i) by dets_dn
        allocate(beta(ndet_new))
        allocate(beta_ind(ndet_new))

        beta(1:ndet_new) = dets_dn(ndet_old + 1:n_det)
        do i=1,ndet_new
          beta_ind(i) = ndet_old + i
        enddo

        call sort_by_first_argument(int(ndet_new),beta,beta_ind)

        ! Next, generate alpha_m1, a list of all up configurations reachable by removing one up electron
        allocate(alpha_m1(ndet_new*nup))
        allocate(alpha_m1_ind(ndet_new*nup))
        call get_n_minus_1_configs(ndet_new,nup,dets_up(ndet_old+1:n_det),alpha_m1)
        do i=1,ndet_new
          alpha_m1_ind(nup*(i-1)+1:nup*i) = ndet_old + i
        enddo
        call sort_by_first_argument(int(ndet_new*nup),alpha_m1,alpha_m1_ind)

! Write out the number of unique beta strings
        n_unique_beta_new = 0
        do i=2,ndet_new
           if(beta(i-1).ne.beta(i)) n_unique_beta_new = n_unique_beta_new+1
        enddo
        write(6,'(/,''n_unique_beta_new,ndet_new='',i8,i10)') n_unique_beta_new,ndet_new

        allocate(is_included(n_det))
        is_included = .false.

        call my_second(2,'create beta and alpha-1 helper lists')

        ! Finally, use helper lists to construct H

        if (present(sparse_ham)) then
          isparse_old = 0
        endif

        do i=1,n_det

          ! Copy over old matrix elements
          if (present(sparse_ham)) then
            if (i<=ndet_old) then
              ! Copy old dets
              n_copy = sparse_ham%nonzero_elements(i)
              isparse_old=isparse_old+n_copy
              isparse=isparse+n_copy
              if (isparse>nonzero_elements) then ! increment nonzero_elements by 41%
                write (6,'(''chem_copy: Reallocating from'',i12,'' ='',es11.4,'' to'',i12'' ='',es11.4)') nonzero_elements, real(nonzero_elements), nonzero_elements+n_copy, real(nonzero_elements+n_copy)
                allocate(tmp_H_indices(nonzero_elements))
                allocate(tmp_H_values(nonzero_elements))
                tmp_H_indices(1_i8b:nonzero_elements) = H_indices(1_i8b:nonzero_elements)
                tmp_H_values(1_i8b:nonzero_elements) = H_values(1_i8b:nonzero_elements)
                deallocate(H_indices,H_values)
                old_nonzero_elements = nonzero_elements
                nonzero_elements = nonzero_elements + n_copy
                allocate(H_indices(nonzero_elements))
                allocate(H_values(nonzero_elements))
                H_indices(1_i8b:old_nonzero_elements) = tmp_H_indices(1_i8b:old_nonzero_elements)
                H_values(1_i8b:old_nonzero_elements) = tmp_H_values(1_i8b:old_nonzero_elements)
                deallocate(tmp_H_indices,tmp_H_values)
              endif
              H_indices(isparse-n_copy+1:isparse) = sparse_ham%indices(isparse_old-n_copy+1:isparse_old)
              H_values(isparse-n_copy+1:isparse) = sparse_ham%values(isparse_old-n_copy+1:isparse_old)
              H_nonzero_elements(i) = sparse_ham%nonzero_elements(i)
            endif
          endif ! present(sparse_ham)

          ! Following subroutine call doesn't need to be changed because dets_up/dn are only used for computing matrix elements,
          ! and the beta and alpha_m1 lists are already shortened to only the information pertaining to the new dets
          call get_connected_dets_in_list(int(i),dets_up(i),dets_dn(i),dets_up(1:n_det),dets_dn(1:n_det),beta,beta_ind,alpha_m1,alpha_m1_ind,n_connected_dets,connected_indices,connected_matrix_elements,is_included)

          ! Add to Hamiltonian if connected_index >= i
          do j=1,n_connected_dets
            if (connected_indices(j)>i.or.(connected_indices(j)==i.and.i>ndet_old)) then
              ! Add to H
              isparse=isparse+1_i8b
              if (isparse>nonzero_elements) then ! increment nonzero_elements by 41%
                write (6,'(''chem_pc: Reallocating from'',i12,'' ='',es11.4,'' to'',i12,'' ='',es11.4)') nonzero_elements, real(nonzero_elements), max(3,int(1.05*nonzero_elements,i8b)), 1.05*nonzero_elements
                allocate(tmp_H_indices(nonzero_elements))
                allocate(tmp_H_values(nonzero_elements))
                tmp_H_indices(1_i8b:nonzero_elements) = H_indices(1_i8b:nonzero_elements)
                tmp_H_values(1_i8b:nonzero_elements) = H_values(1_i8b:nonzero_elements)
                deallocate(H_indices,H_values)
                old_nonzero_elements = nonzero_elements
                nonzero_elements = max(3,int(1.05*nonzero_elements,i8b))
                allocate(H_indices(nonzero_elements))
                allocate(H_values(nonzero_elements))
                H_indices(1_i8b:old_nonzero_elements) = tmp_H_indices(1_i8b:old_nonzero_elements)
                H_values(1_i8b:old_nonzero_elements) = tmp_H_values(1_i8b:old_nonzero_elements)
                deallocate(tmp_H_indices,tmp_H_values)
              endif
              H_indices(isparse)=connected_indices(j)
              H_values(isparse)=connected_matrix_elements(j)
              H_nonzero_elements(i)=H_nonzero_elements(i)+1
              is_included(connected_indices(j)) = .false.
            endif
          enddo ! connected dets

        enddo ! reference dets

      else

        if (brute_force) then

          do i = 1, n_det
             do j = i, n_det
                matrix_element=0._rk
                if (.not.hf_to_psit.or.(i.ne.1.and.j.ne.1)) then
                  if (sym_on) then
                    call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                  else
                    if (i.eq.j) then
                      call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), 0, matrix_element)
                    else
                      call excitation_level(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level)
                      if (excite_level >= 0) call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, matrix_element)
                    endif
                  endif
                endif
                if (abs(matrix_element).gt.1.d-12.or.(hf_to_psit.and.i==1.and.j==1)) then
                  isparse=isparse+1
                  if (isparse>nonzero_elements) then ! increment nonzero_elements by 5%
                    write (6,'(''chem_bf: Reallocating from size'',i14,'' ='',es11.4,'' to size'',i14,'' ='',es11.4)') nonzero_elements, real(nonzero_elements), int(1.05*nonzero_elements,i8b), 1.05*nonzero_elements
                    !Only reallocate by a small amount, since our estimate is quite close
                    allocate(tmp_H_indices(nonzero_elements))
                    allocate(tmp_H_values(nonzero_elements))
                    tmp_H_indices(1_i8b:nonzero_elements) = H_indices(1_i8b:nonzero_elements)
                    tmp_H_values(1_i8b:nonzero_elements) = H_values(1_i8b:nonzero_elements)
                    deallocate(H_indices,H_values)
                    old_nonzero_elements = nonzero_elements
                    nonzero_elements = int(1.05*nonzero_elements,i8b)
                    allocate(H_indices(nonzero_elements))
                    allocate(H_values(nonzero_elements))
                    H_indices(1_i8b:old_nonzero_elements) = tmp_H_indices(1_i8b:old_nonzero_elements)
                    H_values(1_i8b:old_nonzero_elements) = tmp_H_values(1_i8b:old_nonzero_elements)
                    deallocate(tmp_H_indices,tmp_H_values)
                  endif
                  H_indices(isparse)=j
                  H_values(isparse)=matrix_element
                  H_nonzero_elements(i)=H_nonzero_elements(i)+1
                endif
             enddo
          enddo

        else ! binary search

          do i = 1, n_det
             if (i==1) then
               if (hf_to_psit) then
                 H_indices(1) = 1
                 H_values(1) = 0._rk
                 isparse = 1
                 cycle
               endif
             endif
             call find_connected_dets_chem(dets_up(i), dets_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn, norb)
             do k=1,n_connected_dets
                call binary_search(connected_dets_up(k),connected_dets_dn(k),dets_up(1:n_det),dets_dn(1:n_det),j)
                connected_indices(k) = j
             enddo
             ! sort connected_indices
             call sort(n_connected_dets,connected_indices(1:n_connected_dets),numops)
             do k=1,n_connected_dets
                j = connected_indices(n_connected_dets-k+1)
                if (j>0.and.i<=j) then
                   matrix_element=0._rk
                   if (sym_on) then
                        call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                   else
                      if (i==j) then
                        call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), 0, matrix_element)
                      else
                        call excitation_level(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j), excite_level)
                        if (excite_level >= 0) then
                          call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, matrix_element)
                        endif
                      endif
                   endif
                   if (abs(matrix_element).gt.1.d-12.or.(hf_to_psit.and.i==1.and.j==1)) then
                     isparse=isparse+1_i8b
                     if (isparse>nonzero_elements) then ! increment nonzero_elements by 20%
                       write (6,'(''chem_bs: Reallocating from size'',i14,'' ='',es11.4,'' to size'',i14,'' ='',es11.4)') nonzero_elements, real(nonzero_elements), int(1.05*nonzero_elements,i8b), 1.05*nonzero_elements
                       !Only reallocate by a small amount, since our estimate is quite close
                       allocate(tmp_H_indices(nonzero_elements))
                       allocate(tmp_H_values(nonzero_elements))
                       tmp_H_indices(1_i8b:nonzero_elements) = H_indices(1_i8b:nonzero_elements)
                       tmp_H_values(1_i8b:nonzero_elements) = H_values(1_i8b:nonzero_elements)
                       deallocate(H_indices,H_values)
                       old_nonzero_elements = nonzero_elements
                       nonzero_elements = int(1.05*nonzero_elements,i8b)
                       allocate(H_indices(nonzero_elements))
                       allocate(H_values(nonzero_elements))
                       H_indices(1_i8b:old_nonzero_elements) = tmp_H_indices(1_i8b:old_nonzero_elements)
                       H_values(1_i8b:old_nonzero_elements) = tmp_H_values(1_i8b:old_nonzero_elements)
                       deallocate(tmp_H_indices,tmp_H_values)
                     endif
                     H_indices(isparse)=j
                     H_values(isparse)=matrix_element
                     H_nonzero_elements(i)=H_nonzero_elements(i)+1
                  endif
                endif
             enddo
          enddo

        endif ! brute_force

      endif ! use_partial_connections


      if (present(sparse_ham)) then ! Save new sparse Hamiltonian

        sparse_ham%ndet = n_det

        if (allocated(sparse_ham%nonzero_elements))  deallocate(sparse_ham%nonzero_elements)
        if (allocated(sparse_ham%indices))  deallocate(sparse_ham%indices)
        if (allocated(sparse_ham%values))  deallocate(sparse_ham%values)

        allocate(sparse_ham%nonzero_elements(n_det))
        allocate(sparse_ham%indices(isparse))
        allocate(sparse_ham%values(isparse))

        sparse_ham%nonzero_elements = H_nonzero_elements
        sparse_ham%indices = H_indices
        sparse_ham%values = H_values

      endif ! present(sparse_ham)

      write(6,'(/,''n_det='',i12,es10.2,'', # of nonzero elem in H='',i14,es10.2,'', average # of nonzero elem='',f8.1,/)') n_det, real(n_det), isparse, real(isparse), real(isparse)/n_det
      average_connections = real(isparse)/n_det
      call my_second(2, 'generate_sparse_ham_chem_upper_triangular_serial')

  end subroutine generate_sparse_ham_chem_upper_triangular

  !=====================================================================================================
  subroutine generate_sparse_ham_chem_upper_triangular_mpi(dets_up,dets_dn,H_indices,H_nonzero_elements,H_values,time_sym_on,hf_to_psit,max_nonzero_elems,num_nonzero_elems,sparse_ham,remote_det_map,local_det_map,use_hash)
    !---------------------------------------------------------------------------
    ! Description : Compute the Hamiltonian in a given space and store it sparsely, as an upper triangular matrix.
    !
    ! Created     : A Holmes, 17 Dec 2012
    ! Modified    : A Holmes, 11 Jan 2013. Input hf_to_psit replaces the first state (HF) with psi_trial.
    !---------------------------------------------------------------------------

      use types, only : i8b, int_vec
      use more_tools, only : binary_search, get_occ_orbs
      use mpi_routines, only : ncores,shmem_allocate,mpi_barr_in_node,shmem_deallocate,master_core_node
      use generic_sort, only : sort_by_first_argument
      use common_selected_ci, only : sparse_mat
      use mpi_routines, only : det_map,mpi_barr,det_map_l
      use fhash_module__ik_int_list
      use ik_hash_module
      use tools, only : sort_and_merge
      implicit none

      ! Dummy
#ifdef NUM_ORBITALS_GT_127
      type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
      type(ik_vec),pointer,dimension(:),contiguous :: beta,alpha_m1
      type(ik_vec) :: tmp_alpha,tmp_alpha2
#else
      integer(ik),intent(in) :: dets_up(:),dets_dn(:)
      !integer(ik),allocatable :: beta(:),alpha_m1(:)
      integer(ik),pointer,dimension(:),contiguous :: beta,alpha_m1
      integer(ik) :: tmp_alpha,tmp_alpha2
#endif
      type(det_map), optional,intent(in) :: remote_det_map
      type(det_map_l), optional,intent(inout) :: local_det_map
      integer(i8b),allocatable,intent(out) :: H_indices(:),H_nonzero_elements(:)
      real(rk),allocatable,intent(out) :: H_values(:)
      logical,intent(in)               :: time_sym_on
      integer(i8b),optional,intent(in) :: max_nonzero_elems
      integer(i8b),optional,intent(out) :: num_nonzero_elems
      logical, optional, intent(in)     :: use_hash
      logical,intent(in) :: hf_to_psit
      type(sparse_mat),optional,intent(inout) :: sparse_ham
!      integer(i8b), allocatable :: H_rows(:)
      ! Local
      integer(i8b),allocatable :: tmp_H_indices(:)
      real(rk),allocatable :: tmp_H_values(:)
      integer(i8b) :: i
      integer :: k, l, i_det, n_unique_beta_new, i_el
      integer(i8b) :: j
      integer :: excite_level
      integer(i8b) :: n_det
      integer(i8b) :: nonzero_elements,old_nonzero_elements,isparse
      real(rk) :: matrix_element
      logical             :: sym_on
      integer :: n_connected_dets
      integer, allocatable :: connected_indices(:)
      integer :: numops
      logical :: brute_force ! whether or not to use the brute force method (instead of binary search)
      logical :: use_partial_connections ! whether or not to use Sandeep's partial connections idea
      !integer,allocatable :: beta_ind(:),alpha_m1_ind(:)
      integer,pointer,dimension(:),contiguous :: beta_ind,alpha_m1_ind
      integer :: isparse_old,n_copy
      integer(i8b) :: ndet_old,ndet_new,n_det_global
      logical, allocatable :: is_included(:)
      type(fhash_type__ik_int_list) :: beta_hash, alpha_m1_hash, alpha_hash,alpha_prime_hash
      type(fhash_type_iterator__ik_int_list) :: hash_it
      type(int_vec) :: det_indices, det_indices2, det_indices3
      integer,allocatable :: occ(:),unocc(:)
      logical success
      integer n_connected_dets2
      integer, allocatable :: connected_indices2(:)
      allocate(connected_indices2(400000))
      call my_second(1, 'generate_sparse_ham_chem_upper_triangular_mpi')

      sym_on=.false.
      !if (present(time_sym_on) .and. time_sym) sym_on=time_sym_on
      if (time_sym_on .and. time_sym) sym_on=time_sym_on

      n_det=size(dets_up)

      !MJO - set n_det_global
      n_det_global = local_det_map%ndets_global

      allocate(connected_indices(n_det_global))
      allocate(H_nonzero_elements(n_det))
      H_nonzero_elements(:)=0



     !allocate(row_count(n_det))

      ! Count number of nonzero elements of H
      nonzero_elements=0

      use_partial_connections = .true.

      if (use_partial_connections) then
        write (6,'(''Constructing sparse Hamiltonian using partial connections'')') ; call flush(6)
      else
        ! Determine which method to use based on n_det and n_connected_dets (in general, larger matrices are faster with binary search)
        if (real(n_det)>2.1*real(n_connected_dets)*log(real(n_det))) then
          brute_force = .false.
          write (6,'(''Constructing sparse Hamiltonian by generating connections and performing binary search'')') ; call flush(6)
        else
          brute_force = .true.
          write (6,'(''Constructing sparse Hamiltonian by checking for connectedness between all (row,column) pairs'')') ; call flush(6)
        endif
      endif

!      nonzero_elements = int(n_det*real(min(n_det,n_connected_dets))**0.5,i8b) !Just a starting guess
      nonzero_elements = int(n_det*(local_det_map%average_connections+1),i8b) !Just a starting guess, +1 for a little breathing room

      if (present(max_nonzero_elems).and.present(num_nonzero_elems)) then
        num_nonzero_elems=nonzero_elements
        if (num_nonzero_elems>max_nonzero_elems)  return
      endif
      if (nonzero_elements.lt.250) then
        !MJO - Put a sufficiently large starting guess for the very small determinants case. Avoids excessive malloc
        nonzero_elements = 250
      endif
      write(6,'(''allocating H_indices etc. arrays of size nonzero_elements (estimated)='',i10,'' ='',es9.2,'' in generate_sparse_ham_chem'')') nonzero_elements, real(nonzero_elements)
      call flush(6)
      if (ncores>1.and.real(nonzero_elements)>7.e6_rk) then
         !MJO TEMP Removed this error message, it doesn't seem to be a real error anymore
!        write(6,'(''Error: attempted to allocate too much memory ('',i10,'' elements) on each of too many cores ('',i4,''). Try generating the local energies and deterministic projector with a single-core run, then read in the elements during the multi-core run.'')') nonzero_elements, ncores
        !stop "Attempted to allocate too much memory on each of too many cores!" !MJO TEMP PARALLEL
      endif

      call my_second(1, 'allocating H_indices etc.')

      allocate(H_indices(nonzero_elements))
      allocate(H_values(nonzero_elements))
      H_indices(:)=0
      H_values(:)=0._rk

      call my_second(2, 'allocating H_indices etc.')

      ! Generate H

      isparse = 0

      if (use_partial_connections) then

        call my_second(1,'create beta and alpha-1 helper lists')

        if (present(sparse_ham)) then
          ndet_old = local_det_map%ndets_global_old
          ndet_new = local_det_map%ndets_global - ndet_old
        else
          ndet_old = 0
          ndet_new = n_det_global
        endif


        if (present(use_hash)) then
          call beta_hash%reserve(2*int(ndet_new,4))
          do i=ndet_old+1,n_det_global
            call beta_hash%get(remote_det_map%global_dets_dn(i),det_indices)
            call det_indices%append(int(i,4))
            call beta_hash%set(remote_det_map%global_dets_dn(i),det_indices)
          enddo


          call alpha_m1_hash%reserve(2*int(nup*ndet_new,4))

          allocate(occ(nelec))
          do i = ndet_old+1,n_det_global
            call get_occ_orbs(remote_det_map%global_dets_up(i),occ)
            do i_el = 1,nup
              call alpha_m1_hash%get(ibclr(remote_det_map%global_dets_up(i),occ(i_el)-1),det_indices)
              call det_indices%append(int(i,4))
              call alpha_m1_hash%set(ibclr(remote_det_map%global_dets_up(i),occ(i_el)-1),det_indices)
            enddo ! i_el
          enddo ! i_det

          allocate(is_included(n_det_global))
          is_included = .false.

          ! call alpha_hash%reserve(2*int(ndet_new,4))
          ! do i=ndet_old+1,n_det_global
          !   call alpha_hash%get(remote_det_map%global_dets_up(i),det_indices)
          !   call det_indices%append(int(i,4))
          !   call alpha_hash%set(remote_det_map%global_dets_up(i),det_indices)
          ! enddo

          ! MJO Setup lists of indices singly connected to alpha
          ! call alpha_prime_hash%reserve(2*int(ndet_new,4))
          ! do i = ndet_old+1,n_det_global
          !   det_indices2%n = 0 !MJO Clear list
          !   call get_occ_orbs(remote_det_map%global_dets_up(i),occ)
          !   call alpha_hash%get(remote_det_map%global_dets_up(i),det_indices3)
          !   do k=1,det_indices3%n
          !     is_included(det_indices3%list(k)) = .true.
          !   enddo
          !   do i_el = 1,nup
          !     ! MJO now we get the now filled lists from the tables
          !     tmp_alpha = ibclr(remote_det_map%global_dets_up(i),occ(i_el)-1)
          !     call alpha_m1_hash%get(tmp_alpha,det_indices)
          !     do k=1,det_indices%n
          !       if (.not.is_included(det_indices%list(k) )) then
          !         call det_indices2%append(det_indices%list(k))
          !         is_included(det_indices%list(k)) = .true.
          !       endif
          !     enddo
          !   enddo ! i_el
          !   ! MJO reset is_included
          !   do k=1,det_indices3%n
          !     is_included(det_indices3%list(k)) = .false.
          !   enddo
          !   do k=1,det_indices2%n
          !     is_included(det_indices2%list(k)) = .false.
          !   enddo
          !   ! MJO Store in table
          !   call alpha_prime_hash%set(remote_det_map%global_dets_up(i),det_indices2)
          ! enddo ! i_det

          ! call hash_it%begin(alpha_hash)
          ! call alpha_m1_hash%reserve(2*int(nup*ndet_new,4))

          ! ! MJO create temporary alpha_m1 lists
          ! allocate(occ(nup))
          ! do i = 1,alpha_hash%key_count()
          !   det_indices2%n = 0 ! MJO clear list
          !   call hash_it%next(tmp_alpha,det_indices)
          !   call get_occ_orbs(tmp_alpha,occ)
          !   do i_el = 1,nup
          !     tmp_alpha2 = ibclr(tmp_alpha,occ(i_el)-1)
          !     call alpha_m1_hash%get(tmp_alpha2,det_indices2)
          !     do k=1,det_indices%n
          !       call det_indices2%append(det_indices%list(k))
          !     enddo
          !     call alpha_m1_hash%set(tmp_alpha2,det_indices2)
          !   enddo ! i_el
          ! enddo ! i_det

          ! call alpha_prime_hash%reserve(2*int(ndet_new,4))
          ! call hash_it%begin(alpha_hash)
          ! is_included = .false.

          ! do i = 1,alpha_hash%key_count()
          !   det_indices2%n = 0 ! MJO clear list
          !   call hash_it%next(tmp_alpha,det_indices3)
          !   do k=1,det_indices3%n
          !     is_included(det_indices3%list(k)) = .true.
          !   enddo
          !   call get_occ_orbs(tmp_alpha,occ)
          !   do i_el = 1,nup
          !     ! MJO now we get the filled lists from the table
          !     tmp_alpha2 = ibclr(tmp_alpha,occ(i_el)-1)
          !     call alpha_m1_hash%get(tmp_alpha2,det_indices)
          !     do k=1,det_indices%n
          !       if (.not.is_included(det_indices%list(k) )) then
          !         call det_indices2%append(det_indices%list(k))
          !         is_included(det_indices%list(k)) = .true.
          !       endif
          !     enddo
          !   enddo ! i_el
          !     ! MJO reset is_included
          !     do k=1,det_indices3%n
          !       is_included(det_indices3%list(k)) = .false.
          !     enddo
          !     do k=1,det_indices2%n
          !       is_included(det_indices2%list(k)) = .false.
          !     enddo
          !   ! MJO Store in table
          !   call alpha_prime_hash%set(tmp_alpha,det_indices2)
          ! enddo ! i_det

          write(6,*) 'beta key / collisions: ',beta_hash%key_count(),beta_hash%n_collisions(),' alpham1 key/collisions: ',alpha_m1_hash%key_count(),alpha_m1_hash%n_collisions()
        else
          allocate(is_included(n_det_global))
          is_included = .false.

          call shmem_allocate(alpha_m1,ndet_new*nup)
          call shmem_allocate(alpha_m1_ind,ndet_new*nup)
          call shmem_allocate(beta,ndet_new)
          call shmem_allocate(beta_ind,ndet_new)

          if (master_core_node) then
            beta(:) = remote_det_map%global_dets_dn(ndet_old+1:n_det_global)
            do i=1,ndet_new
              beta_ind(i) = ndet_old + i
            enddo

            call sort_by_first_argument(int(ndet_new),beta,beta_ind)
            ! Next, generate alpha_m1, a list of all up configurations reachable by removing one up electron
            !allocate(alpha_m1(ndet_new*nup))
            !allocate(alpha_m1_ind(ndet_new*nup))
            call get_n_minus_1_configs(ndet_new,nup,remote_det_map%global_dets_up(ndet_old+1:n_det_global),alpha_m1)
            do i=1,ndet_new
              alpha_m1_ind(nup*(i-1)+1:nup*i) = ndet_old + i
            enddo
            call sort_by_first_argument(int(ndet_new*nup),alpha_m1,alpha_m1_ind)

            ! Write out the number of unique beta strings
            n_unique_beta_new = 0
            do i=2,ndet_new
              if(beta(i-1).ne.beta(i)) n_unique_beta_new = n_unique_beta_new+1
            enddo
            write(6,'(/,''n_unique_beta_new, ndet_new='',i8,i10)') n_unique_beta_new, ndet_new

          endif

          call mpi_barr_in_node()
        endif
        call my_second(2,'create beta and alpha-1 helper lists')

        ! Finally, use helper lists to construct H

        if (present(sparse_ham)) then
          isparse_old = 0
        endif
        do i=1,local_det_map%ndets
          i_det = local_det_map%local_indices(i)
          ! Copy over old matrix elements
          if (present(sparse_ham)) then
            if (i_det<=ndet_old) then
              ! Copy old dets
              n_copy = sparse_ham%nonzero_elements(i)
              isparse_old=isparse_old+n_copy
              isparse=isparse+n_copy
              if (isparse>nonzero_elements) then ! increment nonzero_elements by 5%
                write (6,'(''chem_mpi_copy: Reallocating from size'',i14,'' ='',es11.4,'' to size'',i14,'' ='',es11.4)') nonzero_elements, real(nonzero_elements), max(3,int(1.05*nonzero_elements,i8b)), 1.05*nonzero_elements
                allocate(tmp_H_indices(nonzero_elements))
                allocate(tmp_H_values(nonzero_elements))
                tmp_H_indices(1_i8b:nonzero_elements) = H_indices(1_i8b:nonzero_elements)
                tmp_H_values(1_i8b:nonzero_elements) = H_values(1_i8b:nonzero_elements)
                deallocate(H_indices,H_values)
                old_nonzero_elements = nonzero_elements
                nonzero_elements = max(3,int(1.05*nonzero_elements,i8b))
                !nonzero_elements = nonzero_elements + n_copy
                allocate(H_indices(nonzero_elements))
                allocate(H_values(nonzero_elements))
                H_indices(1_i8b:old_nonzero_elements) = tmp_H_indices(1_i8b:old_nonzero_elements)
                H_values(1_i8b:old_nonzero_elements) = tmp_H_values(1_i8b:old_nonzero_elements)
                deallocate(tmp_H_indices,tmp_H_values)
              endif
              H_indices(isparse-n_copy+1:isparse) = sparse_ham%indices(isparse_old-n_copy+1:isparse_old)
              H_values(isparse-n_copy+1:isparse) = sparse_ham%values(isparse_old-n_copy+1:isparse_old)
              H_nonzero_elements(i) = sparse_ham%nonzero_elements(i)
            endif
          endif ! present(sparse_ham)

        ! Following subroutine call doesn't need to be changed because dets_up/dn are only used for computing matrix elements,
        ! and the beta and alpha_m1 lists are already shortened to only the information pertaining to the new dets
          !call get_connected_dets_in_list(int(i_det),dets_up(i),dets_dn(i),remote_det_map%global_dets_up(1:n_det_global),remote_det_map%global_dets_dn(1:n_det_global),beta,beta_ind,alpha_m1,alpha_m1_ind,n_connected_dets,connected_indices,connected_matrix_elements)
          if (present(use_hash)) then
            call get_connected_dets_in_list_hash(int(i_det),dets_up(i),dets_dn(i),remote_det_map%global_dets_up(1:n_det_global),remote_det_map%global_dets_dn(1:n_det_global),beta_hash,alpha_m1_hash,n_connected_dets,connected_indices,connected_matrix_elements,is_included)
          else
            call get_connected_dets_in_list(int(i_det),dets_up(i),dets_dn(i),remote_det_map%global_dets_up(1:n_det_global),remote_det_map%global_dets_dn(1:n_det_global),beta,beta_ind,alpha_m1,alpha_m1_ind,n_connected_dets,connected_indices,connected_matrix_elements,is_included)
          endif
         ! is_included = .false.


         !  is_included = .false.
!          call get_connected_dets_in_list_hash_new(int(i_det),dets_up(i),dets_dn(i),remote_det_map%global_dets_up(1:n_det_global),remote_det_map%global_dets_dn(1:n_det_global),beta_hash,alpha_hash,alpha_prime_hash,n_connected_dets2,connected_indices2,connected_matrix_elements,is_included)


        ! Add to Hamiltonian if connected_index >= i
        do j=1,n_connected_dets
           if (connected_indices(j)>i_det.or.(connected_indices(j)==i_det.and.i_det>ndet_old)) then
              ! Add to H
              isparse=isparse+1_i8b
              if (isparse>nonzero_elements) then ! increment nonzero_elements by 41%
                write (6,'(''chem_mpi_pc: Reallocating from size'',i14,'' ='',es11.4,'' to size'',i14,'' ='',es11.4)') nonzero_elements, real(nonzero_elements), max(3,int(1.05*nonzero_elements,i8b)), 1.05*nonzero_elements
                allocate(tmp_H_indices(nonzero_elements))
                allocate(tmp_H_values(nonzero_elements))
                tmp_H_indices(1_i8b:nonzero_elements) = H_indices(1_i8b:nonzero_elements)
                tmp_H_values(1_i8b:nonzero_elements) = H_values(1_i8b:nonzero_elements)
                deallocate(H_indices,H_values)
                old_nonzero_elements = nonzero_elements
                nonzero_elements = max(3,int(1.05*nonzero_elements,i8b))
                allocate(H_indices(nonzero_elements))
                allocate(H_values(nonzero_elements))
                H_indices(1_i8b:old_nonzero_elements) = tmp_H_indices(1_i8b:old_nonzero_elements)
                H_values(1_i8b:old_nonzero_elements) = tmp_H_values(1_i8b:old_nonzero_elements)
                deallocate(tmp_H_indices,tmp_H_values)
              endif
              H_indices(isparse)=connected_indices(j)
              H_values(isparse)=connected_matrix_elements(j)
              H_nonzero_elements(i)=H_nonzero_elements(i)+1
              is_included(connected_indices(j)) = .false.
           endif
          enddo ! connected dets
        enddo ! reference dets

!      do i=1,isparse
!        write(6,*) H_rows(i),H_indices(i),H_values(i)
!      enddo
!      call flush(6)
!      call mpi_barr
!      stop

      else

        if (brute_force) then

          do i = 1, n_det
             do j = i, n_det
                matrix_element=0._rk
                if (.not.hf_to_psit.or.(i.ne.1.and.j.ne.1)) then
                  if (sym_on) then
                    call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                  else
                    if (i.eq.j) then
                      call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), 0, matrix_element)
                    else
                      call excitation_level(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level)
                      if (excite_level >= 0) call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, matrix_element)
                    endif
                  endif
                endif
                if (abs(matrix_element).gt.1.d-12.or.(hf_to_psit.and.i==1.and.j==1)) then
                  isparse=isparse+1
                  if (isparse>nonzero_elements) then ! increment nonzero_elements by 20%
                    write (6,'(''chem_mpi_bf: Reallocating from size'',i14,'' ='',es11.4,'' to size'',i14,'' ='',es11.4)') nonzero_elements, real(nonzero_elements), int(1.05*nonzero_elements,i8b), 1.05*nonzero_elements
                    allocate(tmp_H_indices(nonzero_elements))
                    allocate(tmp_H_values(nonzero_elements))
                    tmp_H_indices(1_i8b:nonzero_elements) = H_indices(1_i8b:nonzero_elements)
                    tmp_H_values(1_i8b:nonzero_elements) = H_values(1_i8b:nonzero_elements)
                    deallocate(H_indices,H_values)
                    old_nonzero_elements = nonzero_elements
                    nonzero_elements = int(1.05*nonzero_elements,i8b)
                    allocate(H_indices(nonzero_elements))
                    allocate(H_values(nonzero_elements))
                    H_indices(1_i8b:old_nonzero_elements) = tmp_H_indices(1_i8b:old_nonzero_elements)
                    H_values(1_i8b:old_nonzero_elements) = tmp_H_values(1_i8b:old_nonzero_elements)
                    deallocate(tmp_H_indices,tmp_H_values)
                  endif
                  H_indices(isparse)=j
                  H_values(isparse)=matrix_element
                  H_nonzero_elements(i)=H_nonzero_elements(i)+1
                endif
             enddo
          enddo

        else ! binary search

          do i = 1, n_det
             if (i==1) then
               if (hf_to_psit) then
                 H_indices(1) = 1
                 H_values(1) = 0._rk
                 isparse = 1
                 cycle
               endif
             endif
             call find_connected_dets_chem(dets_up(i), dets_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn, norb)
             do k=1,n_connected_dets
                call binary_search(connected_dets_up(k),connected_dets_dn(k),dets_up(1:n_det),dets_dn(1:n_det),j)
                connected_indices(k) = j
             enddo
             ! sort connected_indices
             call sort(n_connected_dets,connected_indices(1:n_connected_dets),numops)
             do k=1,n_connected_dets
                j = connected_indices(n_connected_dets-k+1)
                if (j>0.and.i<=j) then
                   matrix_element=0._rk
                   if (sym_on) then
                        call hamiltonian_chem_time_sym(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
                   else
                      if (i==j) then
                        call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), 0, matrix_element)
                      else
                        call excitation_level(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j), excite_level)
                        if (excite_level >= 0) then
                          call hamiltonian_chem(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), excite_level, matrix_element)
                        endif
                      endif
                   endif
                   if (abs(matrix_element).gt.1.d-12.or.(hf_to_psit.and.i==1.and.j==1)) then
                     isparse=isparse+1_i8b
                     if (isparse>nonzero_elements) then ! increment nonzero_elements by 20%
                       write (6,'(''chem_mpi_bs: Reallocating from size'',i14,'' ='',es11.4'' to size'',i14,'' ='',es11.4)') nonzero_elements, real(nonzero_elements), int(1.05*nonzero_elements,i8b), 1.05*nonzero_elements
                       allocate(tmp_H_indices(nonzero_elements))
                       allocate(tmp_H_values(nonzero_elements))
                       tmp_H_indices(1_i8b:nonzero_elements) = H_indices(1_i8b:nonzero_elements)
                       tmp_H_values(1_i8b:nonzero_elements) = H_values(1_i8b:nonzero_elements)
                       deallocate(H_indices,H_values)
                       old_nonzero_elements = nonzero_elements
                       nonzero_elements = int(1.05*nonzero_elements,i8b)
                       allocate(H_indices(nonzero_elements))
                       allocate(H_values(nonzero_elements))
                       H_indices(1_i8b:old_nonzero_elements) = tmp_H_indices(1_i8b:old_nonzero_elements)
                       H_values(1_i8b:old_nonzero_elements) = tmp_H_values(1_i8b:old_nonzero_elements)
                       deallocate(tmp_H_indices,tmp_H_values)
                     endif
                     H_indices(isparse)=j
                     H_values(isparse)=matrix_element
                     H_nonzero_elements(i)=H_nonzero_elements(i)+1
                  endif
                endif
             enddo
          enddo

        endif ! brute_force

      endif ! use_partial_connections

      if (present(sparse_ham)) then ! Save new sparse Hamiltonian

        sparse_ham%ndet = n_det

        if (allocated(sparse_ham%nonzero_elements))  deallocate(sparse_ham%nonzero_elements)
        if (allocated(sparse_ham%indices))  deallocate(sparse_ham%indices)
        if (allocated(sparse_ham%values))  deallocate(sparse_ham%values)

        allocate(sparse_ham%nonzero_elements(n_det))
        allocate(sparse_ham%indices(isparse))
        allocate(sparse_ham%values(isparse))

        sparse_ham%nonzero_elements = H_nonzero_elements
        sparse_ham%indices = H_indices
        sparse_ham%values = H_values

      endif ! present(sparse_ham)

      write(6,'(/,''n_det='',i12,es10.2,'', # of nonzero elem in local H='',i14,es10.2,'', average # of local nonzero elem='',f8.1,/)') local_det_map%ndets_global, real(local_det_map%ndets_global), isparse, real(isparse), real(isparse)/n_det
      !MJO Store our new average connections
      local_det_map%average_connections = real(isparse)/n_det

      if (present(use_hash)) then
        call beta_hash%clear()
        call alpha_m1_hash%clear()
        !      call alpha_prime_hash%clear()
        !      call alpha_hash%clear()
      else
        call shmem_deallocate(beta_ind)
        call shmem_deallocate(alpha_m1_ind)
        call shmem_deallocate(beta)
        call shmem_deallocate(alpha_m1)
      endif


      call my_second(2, 'generate_sparse_ham_chem_upper_triangular_mpi')

  end subroutine generate_sparse_ham_chem_upper_triangular_mpi


  !=====================================================================================================
  subroutine generate_sparse_ham_chem2(dets_up,dets_dn,H_indices,H_nonzero_elements,H_values,time_sym_on)
    !---------------------------------------------------------------------------
    ! Description :  Compute the Hamiltonian in a given space and store it sparsely.
    ! Created     :  H.J. Changlani May 23,2012 - modifying Adam's code, 7 Apr 2012
    !---------------------------------------------------------------------------

      use more_tools,only :binary_search
      use mpi_routines, only : ncores

      implicit none
      ! Dummy
#ifdef NUM_ORBITALS_GT_127
      type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
#else
      integer(ik),intent(in) :: dets_up(:),dets_dn(:)
#endif
      integer(i8b),allocatable,intent(out) :: H_indices(:),H_nonzero_elements(:)
      real(rk),allocatable,intent(out) :: H_values(:)
      logical,intent(in)               :: time_sym_on

      ! Local
     !type(ik_vec),allocatable          :: conn_dets_up(:),conn_dets_dn(:)
     !real(rk),allocatable             :: mat_elts(:)
      integer                          :: nconnected_dets
      integer                          :: i,j,ind
      integer                          :: nonzero_elements,isparse,itotal,n_det
      integer,allocatable              :: H_rows(:),row_count(:)
      logical                          :: sym_on

      call my_second(1, 'generate_sparse_ham_chem2')

      sym_on=.false.
      if (time_sym_on .and. time_sym) sym_on=time_sym_on

      n_det=size(dets_up)

      allocate(H_nonzero_elements(n_det))
      H_nonzero_elements(:)=0

      allocate(row_count(n_det))

      ! Count number of nonzero elements of H
      nonzero_elements=0
      do i = 1, n_det
         call find_connected_dets_chem(dets_up(i), dets_dn(i), nconnected_dets, connected_dets_up, connected_dets_dn, norb, connected_matrix_elements)
         if (abs(connected_matrix_elements(1)).gt.1.d-10) then                         ! Diagonal element
           nonzero_elements = nonzero_elements+1
           H_nonzero_elements(i)=H_nonzero_elements(i)+1
         endif
         do j = 2, nconnected_dets
           call binary_search(connected_dets_up(j),connected_dets_dn(j),dets_up,dets_dn,ind)
           if (abs(connected_matrix_elements(j)).gt.1.d-10) then
                nonzero_elements = nonzero_elements+2
                H_nonzero_elements(i)=H_nonzero_elements(i)+1
                H_nonzero_elements(ind)=H_nonzero_elements(ind)+1
           endif
         enddo
      enddo

      write(6,'(''allocating H_indices etc. arrays of size nonzero_elements='',i10,'' ='',es9.2,'' in generate_sparse_ham_chem2'')') nonzero_elements, real(nonzero_elements)
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

      isparse = 0
      itotal = 0
      do i = 1, n_det
         call find_connected_dets_chem(dets_up(i),dets_dn(i), nconnected_dets, connected_dets_up,connected_dets_dn, norb, connected_matrix_elements)
         if (i>1) then
           itotal=itotal+H_nonzero_elements(i-1)
           H_rows(isparse+1:itotal)=H_rows(isparse)
           isparse = itotal
         endif
         if (abs(connected_matrix_elements(1)).gt.1.d-10) then                         ! Diagonal element
           isparse=isparse+1
           H_indices(isparse)=j
           H_rows(isparse)=i
           H_values(isparse)=connected_matrix_elements(1)
         endif
         do j = 2, nconnected_dets
           call binary_search(connected_dets_up(j),connected_dets_dn(j),dets_up,dets_dn,ind)
           if (abs(connected_matrix_elements(j)).gt.1.d-10) then
                 isparse=isparse+1
                 H_indices(isparse)=ind
                 H_rows(isparse)=i
                 H_values(isparse)=connected_matrix_elements(j)
                 H_values(row_count(H_indices(isparse))) = H_values(isparse)
                 H_indices(row_count(H_indices(isparse))) = H_rows(isparse)
                 row_count(H_indices(isparse)) = row_count(H_indices(isparse)) + 1
           endif
         enddo
      enddo

      call my_second(1, 'generate_sparse_ham_chem2')

  end subroutine generate_sparse_ham_chem2

  !=====================================================================================================
  function upper_bound_connections()
    !use simple combinatorics argument to get upper bound on number of connected dets (can exclude core excitations with n_core_orb > 0)
    implicit none

    integer :: upper_bound_connections
    integer :: count_empty_up,count_empty_dn

    count_empty_up = norb - nup
    count_empty_dn = norb - ndn

    upper_bound_connections = 1 &   !current det
         &         + (nup - n_core_orb) * (nup - n_core_orb - 1) / 2 * (count_empty_up) * (count_empty_up - 1) / 2 &  !double excitations
         &         + (ndn - n_core_orb) * (ndn - n_core_orb - 1) / 2 * (count_empty_dn) * (count_empty_dn - 1) / 2 &
         &         + (nup - n_core_orb) * (ndn - n_core_orb) * (count_empty_up) * (count_empty_dn) &
         &         + (nup - n_core_orb) * (count_empty_up) + (ndn - n_core_orb) * (count_empty_dn) !single excitations

    if (time_sym)  upper_bound_connections = 2*upper_bound_connections
    write (6,'(/,''nup,ndn,norb,empty_up,empty_dn,max_connections='',5i6,i12)') nup,ndn,norb,count_empty_up,count_empty_dn,upper_bound_connections

  end function upper_bound_connections



!=====================================================================================================================
 subroutine matrix_lanczos_on_the_fly(dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,highest_eigenvalue,second_lowest_eigenvalue,short_run_size,sym_on,initial_vector,eps)
 ! Created by  : Hitesh Changlani, March 2012
 ! Modified by : A Holmes, 15 Aug 2012. Modified input to allow for computing H on the fly.
 !               A Holmes, 12 Feb 2013. When doing matrix_lanczos_on_the_fly, first seed with a shorter Lanczos run in a truncated space.
 !               A Holmes, 12 Oct 2016. Added input eps, which allows v' = H v to be
 !               quickly approximated by skipping over contributions to it that are
 !               smaller in magnitude than eps.
!=====================================================================================================================

 use types, only : i8b
 use more_tools, only : matrix_lanczos

 implicit none

 ! Dummy
 real(rk),intent(out) :: lowest_eigenvector(:)
 real(rk),intent(out) :: lowest_eigenvalue
 real(rk),intent(out),optional :: highest_eigenvalue,second_lowest_eigenvalue
 integer,intent(in),optional :: short_run_size ! if present, seed long Lanczos run with a short Lanczos run on this many of the first dets.
 logical,intent(in),optional :: sym_on
 real(rk),intent(in),optional :: initial_vector(:)
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
 type(ik_vec),allocatable    :: temp_i16_up(:),temp_i16_dn(:)!,temp_i_2(:)
#else
 integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
 integer(ik),allocatable    :: temp_i16_up(:),temp_i16_dn(:)!,temp_i_2(:)
#endif
 real(rk),optional,intent(in)  :: eps

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
 real(rk)                   :: epsilon=1.0e-10 ! for Lanczos
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
    write (6,'(''In matrix_lanczos_on_the_fly, performing short run on the first'',i10,'' of'',i12,'' dets'')') short_run_size, size(dets_up)

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

    call generate_sparse_ham_chem_upper_triangular(dets_up(1:short_run_size),dets_dn(1:short_run_size),H_indices,H_nonzero_elements,H_values,sym_on,hf_to_psit=.false.) ! When diagonalizing, we don't want to transform the projector first (hence hf_to_psit=false)
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
    v(1,1)=1._rk ! start with HF
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
         call apply_H_on_the_fly(dets_up,dets_dn,v(:,it),w(:),eps)
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
         write(6,'(''Iteration, Eigenvalue='',i3,f15.9)') it, lowest_eigenvalue
         call flush(6)
         if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<epsilon) then
             converged=.true.
             exit
         endif
         lowest_eigenvalue_prev = lowest_eigenvalue

      enddo

      it=min(it,iterations)
      write(6,'(''n, Lowest eigenvalue ='',i10, f16.10)') n, lowest_eigenvalue ; call flush(6)

      v(:,1)=matmul(v(:,1:it),tridiag(1:it,1))

      if (allocated(eigenvalues)) deallocate(eigenvalues)
      if (allocated(tridiag))     deallocate(tridiag)
      if (allocated(work))        deallocate(work)

  else
      lowest_eigenvalue = eigenvalues(1)
  endif

   lowest_eigenvector(:)=v(:,1)

  end subroutine matrix_lanczos_on_the_fly

  !=====================================================================================================
  subroutine apply_H_on_the_fly(dets_up,dets_dn,v1,v2,eps)
    ! A Holmes, 15 Aug 2012
    ! v2 = H * v1
    ! Assumes dets_up,dn already sorted by label
    ! Edited by A Holmes, 12 Oct 2016. Added optional input eps. If present,
    ! only include contributions in the sum (v2)_i = sum_j H_{ij} (v1)_j
    ! that exceed eps in magnitude

    use more_tools, only : binary_search

    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
#endif
    real(rk),intent(in) :: v1(:)
    real(rk),intent(out) :: v2(:)
    real(rk),optional,intent(in) :: eps

    integer :: i,j,k,n_connected_dets
    real(rk) :: matrix_element
    integer :: excite_level
    real(rk) :: eps_local

    if (present(eps)) then
      eps_local = eps
    else
      eps_local = 0._rk
    endif

    v2(:)=0._rk
    do i=1,size(dets_up)
      ! generate connections = connected_dets_up,dn
      if (v1(i).ne.0._rk) then
        call find_important_connected_dets_chem(dets_up(i), dets_dn(i), eps_local/abs(v1(i)), n_connected_dets, connected_dets_up, connected_dets_dn)
       !call find_connected_dets_chem(dets_up(i), dets_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn,norb)
        ! loop over connections: binary search the dets list for each connection.
        do j=1,n_connected_dets
          call binary_search(connected_dets_up(j),connected_dets_dn(j),dets_up,dets_dn,k)
          if (k>0) then
              if (time_sym .eqv. .true.) then
                      call hamiltonian_chem_time_sym(dets_up(i),dets_dn(i),dets_up(k),dets_dn(k), matrix_element)
              else
                      call excitation_level(dets_up(i),dets_dn(i),dets_up(k),dets_dn(k), excite_level)
                      if (excite_level >=0) call hamiltonian_chem(dets_up(i),dets_dn(i),dets_up(k),dets_dn(k), excite_level, matrix_element)
              endif
              v2(k) = v2(k) + matrix_element * v1(i)
          endif
        enddo
      endif ! v1(i) .ne. 0
    enddo

  end subroutine apply_H_on_the_fly

    !=============================================================================
  subroutine sort_integrals(hf_up,hf_dn,orb_order,orb_order_inv)
    ! Generates the vector orb_order, which maps the integral labels in the dump file to the indices that are in approximately the right order energetically.
    ! e.g., if orbital label 27 in the dump file is the 3rd lowest energy, then orb_order(3)=27.
    ! The energy of each orbital is estimated by putting one up and one down electron in that orbital.
    ! Update the hf determinant to use the sorted integrals
    ! A Holmes, 15 Jan 2014

      use mpi_routines, only : master_core,mpi_bsend
      use more_tools, only : get_occ_orbs

#ifdef NUM_ORBITALS_GT_127
      type(ik_vec),intent(inout) :: hf_up,hf_dn
#else
      integer(ik),intent(inout) :: hf_up,hf_dn
#endif
      integer,intent(out) :: orb_order(norb+1)
      integer,intent(out) :: orb_order_inv(norb+1)
      integer :: i,j
      real(rk) :: tmp_orbital_energies(norb)
      integer :: tmp_orb_sym(norb)
      integer :: occ_up(nup)
      integer :: occ_dn(ndn)

      ! Save old HF det
      call get_occ_orbs(hf_up,hf_dn,occ_up,occ_dn)

      if (master_core) then

        do i=1,norb+1
          orb_order(i) = i
        enddo

        call compute_orbital_energies(hf_up,hf_dn,orbital_energies)

        tmp_orbital_energies=orbital_energies

! If the orbital is occupied in the HF up or dn determinant, make its temporary energy very low
! so that is gets included among the first nup orbitals at the next step (we assume nup >= ndn)
! Caution: don't set it so low that all digits of its energy are lost
        do i=1,norb
          if (btest(hf_up,i-1)) then
            tmp_orbital_energies(i) = tmp_orbital_energies(i)-1.e9_rk
          endif
          if (btest(hf_dn,i-1)) then
            tmp_orbital_energies(i) = tmp_orbital_energies(i)-1.e9_rk
          endif
        enddo

! Now sort the orbitals.  Once an orbital is used set its energy very high
        do i=1,norb
          do j=1,norb
            if (tmp_orbital_energies(j)==minval(tmp_orbital_energies)) then
              orb_order(i)=j
              orb_order_inv(j)=i
              tmp_orbital_energies(j) = 1.e99_rk
              exit
            endif
          enddo
        enddo

        write (6,'(/,''Orbitals reordered: norb='',i6)') norb
        do j=1,norb
          write (6,'(''orb_order('',i4,'')='',i4,'' energy='',f16.8,'' symmetry='',i2)') j, orb_order(j), orbital_energies(orb_order(j)), orbital_symmetries(orb_order(j))
        enddo

        call flush(6)

        ! Now sort orbital_symmetries to the new order, such that
        ! new_orb_sym(unsorted_indices)=old_orb_sym(sorted_indices)

        !tmp_orb_sym = orbital_symmetries
        orbital_symmetries(1:norb) = orbital_symmetries(orb_order(1:norb))
        orbital_energies(1:norb) = orbital_energies(orb_order(1:norb))
        do i=1,norb
          !orbital_symmetries(i) = tmp_orb_sym(orb_order(i))
          if(i.ge.2) then
            if (orbital_energies(i)+1.e-9.lt.orbital_energies(i-1)) then
              write(6,'(''Warning: orbital_energy'',i4,f10.6,'' <orbital_energy'',i4,f10.6,''. The HF filled orbitals are usually the ones with lowest energy, but not always'')') i, orbital_energies(i), i-1, orbital_energies(i-1)
            endif
          endif
        enddo

      endif ! master_core

      call mpi_bsend(orb_order)
      call mpi_bsend(orb_order_inv)
      call mpi_bsend(orbital_symmetries)

      ! Update HF det
      hf_up = 0_ik
      do i=1,nup
        occ_up(i) = orb_order_inv(occ_up(i))
        hf_up = ibset(hf_up,occ_up(i)-1)
      enddo

      hf_dn = 0_ik
      do i=1,ndn
        occ_dn(i) = orb_order_inv(occ_dn(i))
        hf_dn = ibset(hf_dn,occ_dn(i)-1)
      enddo

  end subroutine sort_integrals

  !=============================================================================
  function get_inverse_dih(i)
  ! Returns the value of the inverse of the dih irreducible representation
  ! M Otten, 20 Jan 2015
    integer, intent(in) :: i
    integer             :: get_inverse_dih,lz,gu

    !I am making explicit use of Adam's indices.
    !Where the inverse is:
    !itself, if lz=0
    !itself+2, if lz is positive
    !itself-2, if lz is negative
    if(d_infinity_h==1) then
       call get_lz(i,lz,gu)
       if(lz.gt.0) then
          get_inverse_dih = i+2
       elseif (lz.lt.0) then
          get_inverse_dih = i-2
       else
          get_inverse_dih = i
       endif
    else
       get_inverse_dih = i
    endif

    return
  end function get_inverse_dih

  !=============================================================================
  function product_table(i,j)
  ! Returns the value of the product table element
  ! A Holmes, 16 Jan 2014
    integer,intent(in) :: i,j
    integer :: product_table
    integer :: lz_i,lz_j,gu_i,gu_j

    if (d_infinity_h==1) then
      call get_lz(i,lz_i,gu_i)
      call get_lz(j,lz_j,gu_j)
      call get_ind(lz_i+lz_j,mod(gu_i+gu_j,2),product_table)
      !write(6,'(''i,j,lz_i,lz_j,gu_i,gu_j,product_table'',9i6)') i,j,lz_i,lz_j,gu_i,gu_j,product_table
    else
      product_table = product_table_elems(i,j)
    endif

  end function product_table

  !=============================================================================
  subroutine get_lz(ind,lz,gu)
  ! Given an index number in d_inf_h, compute the corresponding L_z (lz) and g/u (gu=0 for g, gu=1 for u)
  ! A Holmes, 16 Jan 2014
    integer,intent(in) :: ind
    integer,intent(out) :: lz,gu ! gu = 0g, 1u

    if (ind<=2) then
      lz = 0
      gu = ind-1
    else
      lz = 1+(ind-3)/4
      if (mod((ind-1)/2,2)==0) lz=-lz
      gu = mod(ind+1,2)
    endif

  end subroutine get_lz

  !=============================================================================
  subroutine get_ind(lz,gu,ind)
  ! Given the L_z and whether it is a g or a u, compute the index in d_inf_h
  ! A Holmes, 16 Jan 2014
    integer,intent(in) :: lz,gu ! gu = 0g, 1u
    integer,intent(out) :: ind

    if (lz==0) then
      ind = gu+1
    else
      ind = 4*abs(lz)-1+gu
      if (lz<0)  ind=ind+2
    endif

  end subroutine get_ind

  !=============================================================================
  integer(i8b) function integral_index(i,j,k,l)
     ! Get the index of physical order UMAT element <IJ|KL>.
     ! Indices are internally reordered such that I>K, J>L,(I,K)>(J,L)
     ! Note: (i,k)>(j,l) := (k>l) || ((k==l)&&(i>j))
     ! In:
     !    I,J,K,L: orbital indices. These refer to spin orbitals in
     !      unrestricted calculations and spatial orbitals in restricted
     !      calculations.
     ! Out: a unique index that would be unchanged under the valid permutations of (i,j,k,l)
     ! From George Booth

     implicit none
     integer, intent(in) :: i,j,k,l
     integer(i8b) a,b

     ! Combine indices i and j, ensuring i>j
     a = combine_2(i,j)

     ! Combine indices k and l, ensuring k>l
     b = combine_2(k,l)

     ! Combine (ij) and (kl) in a unique way  (k > l or if k = l then i > j)
     if(a.gt.b) then
       integral_index=(a*(a-1))/2+b
     else
       integral_index=(b*(b-1))/2+a
     endif

  end function integral_index

  !=============================================================================
  integer(i8b) function combine_2_indices(i,j)
     ! Out: a unique index that would be unchanged under the valid permutations of (i,j)
     ! A Holmes, 2 Apr 2015. Simplified 2-index version of George Booth's routine above

     implicit none
     integer, intent(in) :: i,j

     !Combine indices i and j, ensuring i>=j
     if(i.gt.j) then
         combine_2_indices=(i*(i-1))/2+j
     else
         combine_2_indices=(j*(j-1))/2+i
     endif

  end function combine_2_indices

  !=============================================================================
  integer(i8b) function same_index(excite_from_1,excite_from_2,excite_to_1,excite_to_2)
    ! Combines all 4 indices into a single index of four_orbitals_same_spin
    ! A Holmes, 5 Apr 2015

    integer,intent(in) :: excite_from_1,excite_from_2,excite_to_1,excite_to_2
    integer :: i,j,k,l

    i = excite_from_1
    if (i>norb)  i = i-norb
    j = excite_from_2
    if (j>norb)  j = j-norb
    k = excite_to_1
    if (k>norb)  k = k-norb
    l = excite_to_2
    if (l>norb)  l = l-norb

    same_index = (combine_2_indices(i,j)-1)*norb**2+(k-1)*norb+l

  end function same_index

  !=============================================================================
  integer(i8b) function opposite_index(excite_from_1,excite_from_2,excite_to_1,excite_to_2)
    ! Combines all 4 indices into a single index of four_orbitals_opposite_spin
    ! A Holmes, 5 Apr 2015

    integer,intent(in) :: excite_from_1,excite_from_2,excite_to_1,excite_to_2
    integer :: i,j,k,l

    i = excite_from_1
    if (i>norb)  i = i-norb
    j = excite_from_2
    if (j>norb)  j = j-norb
    k = excite_to_1
    if (k>norb)  k = k-norb
    l = excite_to_2
    if (l>norb)  l = l-norb

    opposite_index = (i-1)*int(norb,i8b)**3+(j-1)*norb**2+(k-1)*norb+l

  end function opposite_index

  !=============================================================================
  real(rk) function p_first_hole(excite_from_1,excite_from_2,excite_to_1)
    ! Returns probability of choosing first hole
    ! A Holmes, 14 Sep 2015

    integer,intent(in) :: excite_from_1,excite_from_2,excite_to_1

      if (excite_from_1<=norb) then
        if (excite_from_2<=norb) then ! same spin, up
          p_first_hole = three_orbital_probabilities_same_spin(excite_from_1,excite_from_2,excite_to_1)
        else ! opposite spin, 1=up, 2=dn
          p_first_hole = three_orbital_probabilities_opposite_spin(excite_from_1,excite_from_2-norb,excite_to_1)
        endif
      else
        if (excite_from_2<=norb) then ! opposite spin, 1=dn, 2=up
          p_first_hole = three_orbital_probabilities_opposite_spin(excite_from_1-norb,excite_from_2,excite_to_1-norb)
        else ! same spin, dn
          p_first_hole = three_orbital_probabilities_same_spin(excite_from_1-norb,excite_from_2-norb,excite_to_1-norb)
        endif
      endif

  end function p_first_hole

  !=============================================================================
  real(rk) function p_second_hole(excite_from_1,excite_from_2,excite_to_1,excite_to_2)
    ! Returns probability of choosing second hole
    ! A Holmes, 14 Sep 2015

    integer,intent(in) :: excite_from_1,excite_from_2,excite_to_1,excite_to_2

      if ((excite_from_1<=norb.and.excite_from_2<=norb).or.(excite_from_1>norb.and.excite_from_2>norb)) then
        p_second_hole = four_orbital_probabilities_same_spin(same_index(excite_from_1,excite_from_2,excite_to_1,excite_to_2))
      else
        p_second_hole = four_orbital_probabilities_opposite_spin(opposite_index(excite_from_1,excite_from_2,excite_to_1,excite_to_2))
      endif

  end function p_second_hole

  !=============================================================================
  real(rk) function Htot(excite_from_1,excite_from_2,excite_to_1)
    ! Returns Htot, the sum of magnitudes of double excitation matrix elements including the indices excite_from_1,excite_from_2,excite_to_1
    ! A Holmes, 14 Sep 2015

    integer,intent(in) :: excite_from_1,excite_from_2,excite_to_1

      if (excite_from_1<=norb) then
        if (excite_from_2<=norb) then ! same spin, up
          Htot = Htot_same(combine_2_indices(excite_from_1,excite_from_2),excite_to_1)
        else ! opposite spin, 1=up, 2=dn
          Htot = Htot_opposite(excite_from_1,excite_from_2-norb,excite_to_1)
        endif
      else
        if (excite_from_2<=norb) then ! opposite spin, 1=dn, 2=up
          Htot = Htot_opposite(excite_from_1-norb,excite_from_2,excite_to_1-norb)
        else ! same spin, dn
          Htot = Htot_same(combine_2_indices(excite_from_1-norb,excite_from_2-norb),excite_to_1-norb)
        endif
      endif

  end function Htot

  !=============================================================================
  real(rk) function compute_single_elem(det_up,det_dn,excite_from_1,excite_to_1)
    ! Returns the single excitation matrix element corresponding to the excitation from excite_from_1 to excite_from_2
    ! Note: excite_from/to_1 > norb corresponds to dn excitation
    ! A Holmes, 14 Sep 2015

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
#else
    integer(ik),intent(in) :: det_up,det_dn
#endif
    integer,intent(in) :: excite_from_1,excite_to_1

      if (excite_from_1<=norb) then
        call hamiltonian_chem(det_up,det_dn,ibset(ibclr(det_up,excite_from_1-1),excite_to_1-1),det_dn,1,compute_single_elem)
      else
        call hamiltonian_chem(det_up,det_dn,det_up,ibset(ibclr(det_dn,excite_from_1-norb-1),excite_to_1-norb-1),1,compute_single_elem)
      endif

  end function compute_single_elem

  !=============================================================================
  subroutine choose_first_hole(excite_from_1,excite_from_2,excite_to_1,excite_spin)
    ! Chooses excite_from_2 using Alias method
    ! A Holmes, 24 Feb 2016
    use more_tools, only : sample_alias

    integer,intent(in) :: excite_from_1,excite_from_2
    integer,intent(out) :: excite_to_1
    integer,intent(out) :: excite_spin

    if (excite_from_1<=norb .and. excite_from_2<=norb) then ! same spin, up
      excite_spin = 1
      excite_to_1 = sample_alias(norb,J_three_orbital_probabilities_same_spin(excite_from_1,excite_from_2,:),q_three_orbital_probabilities_same_spin(excite_from_1,excite_from_2,:))
    elseif (excite_from_1>norb .and. excite_from_2>norb) then ! same spin, dn
      excite_spin = -1
      excite_to_1 = sample_alias(norb,J_three_orbital_probabilities_same_spin(excite_from_1-norb,excite_from_2-norb,:),q_three_orbital_probabilities_same_spin(excite_from_1-norb,excite_from_2-norb,:)) + norb
    else ! opposite spin
      excite_spin = 0
      ! Make excite_to_1 be of the same spin as excite_from_1
      if (excite_from_1>norb) then
        excite_to_1 = sample_alias(norb,J_three_orbital_probabilities_opposite_spin(excite_from_1-norb,excite_from_2,:),q_three_orbital_probabilities_opposite_spin(excite_from_1-norb,excite_from_2,:))+norb
      else
        excite_to_1 = sample_alias(norb,J_three_orbital_probabilities_opposite_spin(excite_from_1,excite_from_2-norb,:),q_three_orbital_probabilities_opposite_spin(excite_from_1,excite_from_2-norb,:))
      endif
    endif

  end subroutine choose_first_hole

  !=============================================================================
  integer function choose_second_hole(excite_from_1,excite_from_2,excite_to_1,same_spin)
    ! Chooses excite_from_2 using Alias method
    ! A Holmes, 14 Sep 2015
    use more_tools, only : sample_alias

    integer,intent(in) :: excite_from_1,excite_from_2,excite_to_1
    logical,intent(in) :: same_spin

      if (.not.same_spin) then ! opposite spin excitation

        choose_second_hole = sample_alias(norb,J_four_orbital_probabilities_opposite_spin(opposite_index(excite_from_1,excite_from_2,excite_to_1,1):opposite_index(excite_from_1,excite_from_2,excite_to_1,norb)),q_four_orbital_probabilities_opposite_spin(opposite_index(excite_from_1,excite_from_2,excite_to_1,1):opposite_index(excite_from_1,excite_from_2,excite_to_1,norb)))
        if (excite_from_1<=norb) then ! opposite spin and 1 is up, so make 2 dn
          choose_second_hole = choose_second_hole + norb
        endif
      else ! same spin excitation
        choose_second_hole = sample_alias(norb,J_four_orbital_probabilities_same_spin(same_index(excite_from_1,excite_from_2,excite_to_1,1):same_index(excite_from_1,excite_from_2,excite_to_1,norb)),q_four_orbital_probabilities_same_spin(same_index(excite_from_1,excite_from_2,excite_to_1,1):same_index(excite_from_1,excite_from_2,excite_to_1,norb)))
        if (excite_from_1>norb) then ! same spin and 1 is dn so make 2 dn too
          choose_second_hole = choose_second_hole + norb
        endif
      endif

  end function choose_second_hole

  !=============================================================================
  subroutine check_heatbath_unbiased(is_heatbath_unbiased)
    ! Checks a sufficient (but not necessary) condition for efficient heatbath to work
    ! on the current system: is_heatbath_unbiased = (max(nup,ndn) > N_{orb.uniq.sym.}),
    ! where N_{orb.uniq.sym.} is the number of spatial orbitals of unique symmetry.
    ! A Holmes, 9 Oct 2015

    logical,intent(out) :: is_heatbath_unbiased
    real(rk) :: a
    integer :: p,q,s,n_orb_uniq_sym
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec) :: zero_det
    zero_det = 0_ik
#endif

      n_orb_uniq_sym = 0
      do p=1,norb
        s = 0
        do q=1,norb
          if (p==q)  cycle
#ifdef NUM_ORBITALS_GT_127
          call hamiltonian_chem(ibset(zero_det,p-1),zero_det,ibset(zero_det,q-1),zero_det,1,a)
#else
          call hamiltonian_chem(ibset(0_ik,p-1),0_ik,ibset(0_ik,q-1),0_ik,1,a)
#endif
          if (abs(a)>1.e-10_rk) then
            s = s + 1
            exit
          endif
        enddo
        if (s==0) then
          n_orb_uniq_sym = n_orb_uniq_sym + 1
          write (6,*) "Orbital",p,"has unique symmetry",orbital_symmetries(p); call flush(6)
        endif
      enddo

      write (6,*) "max(nup,ndn)=",max(nup,ndn); call flush(6)
      write (6,*) "number of orbitals of unique symmetry=",n_orb_uniq_sym; call flush(6)

      if (max(nup,ndn)>n_orb_uniq_sym) then
        is_heatbath_unbiased = .true.
      else
        is_heatbath_unbiased = .false.
      endif


  end subroutine check_heatbath_unbiased

  !=============================================================================
  subroutine compute_orbital_energies(hf_up,hf_dn,orbital_energies)
    ! A Holmes, 11 Feb 2016
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: hf_up,hf_dn
#else
    integer(ik),intent(in) :: hf_up,hf_dn
#endif
    real(rk),intent(out) :: orbital_energies(:)
    real(rk) :: direct_energy,exchange_energy
    integer :: i,j

      write (6,*) "computing orbital energies with hf_up,hf_dn=",hf_up,hf_dn; call flush(6)

      do i=1,norb

        orbital_energies(i) = integral_value(i,i,norb+1,norb+1)

        exchange_energy = 0._rk
        direct_energy = 0._rk

        !exchange
        !up electron only exchanges with other up spins
        do j = 1, norb
           if (j.ne.i.and.btest(hf_up, j - 1)) then
              exchange_energy = exchange_energy - integral_value(i, j, j, i)
           endif
           if (j.ne.i.and.btest(hf_dn, j - 1)) then
              exchange_energy = exchange_energy - integral_value(i, j, j, i)
           endif
        enddo

        !direct
        do j = 1, norb
           !up interacting with up
           if (j.ne.i.and.btest(hf_up, j - 1)) then
              direct_energy = direct_energy + integral_value(i, i, j, j)
           endif
        enddo
        do j = 1, norb
           !up interacting with dn
           if (btest(hf_dn, j - 1)) then
              direct_energy = direct_energy + integral_value(i, i, j, j)
           endif
        enddo

        do j = 1, norb
           !dn interacting with dn
           if (j.ne.i.and.btest(hf_dn, j - 1)) then
              direct_energy = direct_energy + integral_value(i, i, j, j)
           endif
        enddo
        do j = 1, norb
           !dn interacting with up
           if (btest(hf_up, j - 1)) then
              direct_energy = direct_energy + integral_value(i, i, j, j)
           endif
        enddo

!       write(6,'(''orbital_energies(i), exchange_energy, direct_energy'',9f12.6)') orbital_energies(i), exchange_energy, direct_energy

        orbital_energies(i) = orbital_energies(i) + .5_rk*(exchange_energy + direct_energy)

      enddo

  end subroutine compute_orbital_energies

  !===========================================================================
  subroutine mp2_deterministic(hf_up,hf_dn,energy_correction)
  ! A Holmes, 10 Feb 2016

  use common_walk, only : n_connected_dets
  use common_run, only : connected_dets_up,connected_dets_dn

#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: hf_up,hf_dn
  type(ik_vec) :: tmp
#else
  integer(ik),intent(in) :: hf_up,hf_dn
  integer(ik) :: tmp
#endif
  real(rk),intent(out) :: energy_correction
  real(rk) :: den
  integer :: i,j

    call find_connected_dets_chem(hf_up, hf_dn, n_connected_dets, connected_dets_up, connected_dets_dn, norb, connected_matrix_elements)
    energy_correction = 0._rk

    do i=2,n_connected_dets ! Start from 2 because 1 is the HF det
      den = 0._rk

      ! up electron(s)
      tmp = iand(hf_up,not(connected_dets_up(i)))
      if (tmp.ne.0_ik) then
        j = trailz(tmp)+1
        den = den + orbital_energies(j)
        tmp = ibclr(tmp,j-1)
      endif
      if (tmp.ne.0_ik) then
        j = trailz(tmp)+1
        den = den + orbital_energies(j)
      endif

      ! dn electron(s)
      tmp = iand(hf_dn,not(connected_dets_dn(i)))
      if (tmp.ne.0_ik) then
        j = trailz(tmp)+1
        den = den + orbital_energies(j)
        tmp = ibclr(tmp,j-1)
      endif
      if (tmp.ne.0_ik) then
        j = trailz(tmp)+1
        den = den + orbital_energies(j)
      endif

      ! up valence orbital(s)
      tmp = iand(not(hf_up),connected_dets_up(i))
      if (tmp.ne.0_ik) then
        j = trailz(tmp)+1
        den = den - orbital_energies(j)
        tmp = ibclr(tmp,j-1)
      endif
      if (tmp.ne.0_ik) then
        j = trailz(tmp)+1
        den = den - orbital_energies(j)
      endif

      ! dn valence orbital(s)
      tmp = iand(not(hf_dn),connected_dets_dn(i))
      if (tmp.ne.0_ik) then
        j = trailz(tmp)+1
        den = den - orbital_energies(j)
        tmp = ibclr(tmp,j-1)
      endif
      if (tmp.ne.0_ik) then
        j = trailz(tmp)+1
        den = den - orbital_energies(j)
      endif

      energy_correction = energy_correction + connected_matrix_elements(i)**2/den

    enddo

  end subroutine mp2_deterministic

  !===========================================================================
  logical function is_occupied(det_up,det_dn,iorb)
  ! Returns true if orbital "iorb" is occupied in det_up/dn
  ! Uses the convention used in the efficient heatbath code,
  ! where iorb<=norb is an up orbital and iorb>norb corresponds
  ! to dn orbital iorb-norb
  ! A Holmes, 24 Feb 2016

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
#else
    integer(ik),intent(in) :: det_up,det_dn
#endif
    integer,intent(in) :: iorb

    if (iorb<=norb) then
      is_occupied = btest(det_up,iorb-1)
    else
      is_occupied = btest(det_dn,iorb-norb-1)
    endif

  end function is_occupied

  !=============================================================================
  real(rk) function double_excitation_matrix_element(det_up,det_dn,p,q,r,s)
  ! Returns the double excitation matrix element corresponding to the
  ! excitation (p,q -> r,s) from the determinant det_up/dn
  ! A Holmes, 24 Feb 2016

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
#else
    integer(ik),intent(in) :: det_up,det_dn
#endif
    integer,intent(in) :: p,q,r,s

    if (p<=norb) then
      if (r>norb) then
        double_excitation_matrix_element = 0._rk
        return
      endif
      if (q<=norb) then
        if (s>norb) then
          double_excitation_matrix_element = 0._rk
          return
        endif
        if (time_sym) then
          call hamiltonian_chem_time_sym(det_up,det_dn,ibset(ibset(ibclr(ibclr(det_up,p-1),q-1),r-1),s-1),det_dn,double_excitation_matrix_element)
        else
          call hamiltonian_chem(det_up,det_dn,ibset(ibset(ibclr(ibclr(det_up,p-1),q-1),r-1),s-1),det_dn,2,double_excitation_matrix_element)
        endif
      else ! q>norb
        if (s<=norb) then
          double_excitation_matrix_element = 0._rk
          return
        endif
        if (time_sym) then
          call hamiltonian_chem_time_sym(det_up,det_dn,ibset(ibclr(det_up,p-1),r-1),ibset(ibclr(det_dn,q-norb-1),s-norb-1),double_excitation_matrix_element)
        else
          call hamiltonian_chem(det_up,det_dn,ibset(ibclr(det_up,p-1),r-1),ibset(ibclr(det_dn,q-norb-1),s-norb-1),2,double_excitation_matrix_element)
        endif
      endif
    else ! p>norb
      if (r<=norb) then
        double_excitation_matrix_element = 0._rk
        return
      endif
      if (q<=norb) then
        if (s>norb) then
          double_excitation_matrix_element = 0._rk
          return
        endif
        if (time_sym) then
          call hamiltonian_chem_time_sym(det_up,det_dn,ibset(ibclr(det_up,q-1),s-1),ibset(ibclr(det_dn,p-norb-1),r-norb-1),double_excitation_matrix_element)
        else
          call hamiltonian_chem(det_up,det_dn,ibset(ibclr(det_up,q-1),s-1),ibset(ibclr(det_dn,p-norb-1),r-norb-1),2,double_excitation_matrix_element)
        endif
      else ! q>norb
        if (s<=norb) then
          double_excitation_matrix_element = 0._rk
          return
        endif
        if (time_sym) then
          call hamiltonian_chem_time_sym(det_up,det_dn,det_up,ibset(ibset(ibclr(ibclr(det_dn,p-norb-1),q-norb-1),r-norb-1),s-norb-1),double_excitation_matrix_element)
        else
          call hamiltonian_chem(det_up,det_dn,det_up,ibset(ibset(ibclr(ibclr(det_dn,p-norb-1),q-norb-1),r-norb-1),s-norb-1),2,double_excitation_matrix_element)
        endif
      endif
    endif

  end function double_excitation_matrix_element

  !=============================================================================
  real(rk) function double_excitation_matrix_element_no_ref(p,q,r,s)
  ! Returns the double excitation matrix element corresponding to the
  ! excitation (p,q -> r,s) from the determinant det_up/dn
  ! A Holmes, 11 Mar 2016

    integer,intent(in) :: p,q,r,s

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec) :: zero_det
#else
    integer(ik) :: zero_det
#endif
    real(rk) :: elem

    if (p==q.or.r==s.or.p==r.or.q==s.or.p==s.or.q==r) then
      double_excitation_matrix_element_no_ref = 0._rk
      return
    endif

    zero_det = 0_ik

    if (p<=norb.and.q<=norb) then ! same spin, up
      call hamiltonian_chem(ibset(ibset(zero_det,p-1),q-1),zero_det,ibset(ibset(zero_det,r-1),s-1),zero_det,2,elem)
    elseif (p>norb.and.q>norb) then ! same spin, dn
      call hamiltonian_chem(zero_det,ibset(ibset(zero_det,p-norb-1),q-norb-1),zero_det,ibset(ibset(zero_det,r-norb-1),s-norb-1),2,elem)
    else ! opposite spin
      call hamiltonian_chem(ibset(zero_det,p-1),ibset(zero_det,q-norb-1),ibset(zero_det,r-1),ibset(zero_det,s-norb-1),2,elem)
    endif

    double_excitation_matrix_element_no_ref = elem

  end function double_excitation_matrix_element_no_ref

  !===========================================================================
  subroutine get_new_diag_elem(info,new_det_up,new_det_dn,new_diag_elem)
  ! Get new diagonal matrix element in O(N) time by using the diagonal matrix element from
  ! a reference determinant that is one double excitation away
  ! Computing the new diagonal matrix element from scratch is an O(N^2) operation
  ! This is most important for quickly computing the 2nd-order PT expression of the energy in HCI,
  ! but it is likely worth doing for SQMC as well
  ! A Holmes, 9 Mar 2016

    use more_tools, only : get_occ_orbs
    use common_run, only : diag_elem_info

    implicit none

    type(diag_elem_info),intent(in) :: info
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: new_det_up,new_det_dn
#else
    integer(ik),intent(in) :: new_det_up,new_det_dn
#endif
    real(rk),intent(out) :: new_diag_elem
    real(rk) :: old_diag_elem
    integer :: p_in,q_in,r_in,s_in
    integer :: p,q,r,s
    logical :: p_up,q_up
    integer :: i,j

    ! Unpack 'info' input
    old_diag_elem = info%old_diag_elem
    p_in = info%p
    q_in = info%q
    r_in = info%r
    s_in = info%s

    p = p_in
    p_up = (p<=norb)
    if (p>norb)  p = p-norb

    q = q_in
    q_up = (q<=norb)
    if (q>norb)  q = q-norb

    r = r_in
    if (r>norb)  r = r-norb

    s = s_in
    if (s>norb)  s = s-norb

    ! One-body part: E += h(r) + h(s) - h(p) - h(q)
    new_diag_elem = old_diag_elem + integral_value(r,r,norb+1,norb+1) + integral_value(s,s,norb+1,norb+1) - integral_value(p,p,norb+1,norb+1) - integral_value(q,q,norb+1,norb+1)

    ! O(1) Two-body direct part: E += direct(r,s) - direct(p,q)
    new_diag_elem = new_diag_elem + integral_value(r,r,s,s) - integral_value(p,p,q,q)

    ! O(1) Two-body exchange part: if (same_spin) then E += -exchange(r,s) + exchange(p,q)
    if (p_up.eqv.q_up) then
      new_diag_elem = new_diag_elem - integral_value(r,s,s,r) + integral_value(p,q,q,p)
    endif

    call get_occ_orbs(new_det_up,new_det_dn,occ_up,occ_dn)

    ! O(N) Two-body direct part: E += sum_{i in occ. but not in (r,s)} direct(i,r) + direct(i,s) - direct(i,p) - direct(i,q)
    do j=1,nup
      i = occ_up(j)
      if (i==r_in.or.i==s_in)  cycle ! only need to check r and s because 'occupied' orbitals are those occupied in new det, which includes r and s but not p and q
      new_diag_elem = new_diag_elem + integral_value(i,i,r,r) + integral_value(i,i,s,s) - integral_value(i,i,p,p) - integral_value(i,i,q,q)
    enddo
    do j=1,ndn
      i = occ_dn(j)
      if (i==r_in-norb.or.i==s_in-norb)  cycle
      new_diag_elem = new_diag_elem + integral_value(i,i,r,r) + integral_value(i,i,s,s) - integral_value(i,i,p,p) - integral_value(i,i,q,q)
    enddo

    ! O(N) Two-body exchange part: E += sum_{i in occ. but not in (r,s)} -exchange(i,r) - exchange(i,s) + exchange(i,p) + exchange(i,q)
    if (p_up.or.q_up) then
      do j=1,nup
        i = occ_up(j)
        if (i==r_in.or.i==s_in)  cycle
        if (p_up)  new_diag_elem = new_diag_elem - integral_value(i,r,r,i) + integral_value(i,p,p,i)
        if (q_up)  new_diag_elem = new_diag_elem - integral_value(i,s,s,i) + integral_value(i,q,q,i)
      enddo
    endif
    if ((.not.p_up).or.(.not.q_up)) then
      do j=1,ndn
        i = occ_dn(j)
        if (i==r_in-norb.or.i==s_in-norb)  cycle
        if (.not.p_up)  new_diag_elem = new_diag_elem - integral_value(i,r,r,i) + integral_value(i,p,p,i)
        if (.not.q_up)  new_diag_elem = new_diag_elem - integral_value(i,s,s,i) + integral_value(i,q,q,i)
      enddo
    endif

  end subroutine get_new_diag_elem

  !===========================================================================
  subroutine sort_rs_absH(arr)
    ! Sort arr of rs_absH by absH in decreasing order

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! From http://jblevins.org/mirror/amiller/qsort.f90
    ! Modified to be of type rs_absH by A Holmes, 24 Jan 2016

    implicit none
    type(rs_absH),intent(inout) :: arr(:)

    call quick_sort(1, size(arr))

    contains

      recursive subroutine quick_sort(left_end, right_end)

        integer, intent(in) :: left_end, right_end

        !     local variables
        integer             :: i, j
        type(rs_absH)          :: reference, temp
        integer, parameter  :: max_simple_sort_size = 6

        if (right_end<left_end+max_simple_sort_size) then
          ! use interchange sort for small arrs
          call interchange_sort(left_end, right_end)
        else
          ! use partition ("quick") sort
          reference = arr((left_end+ right_end)/2)
          i = left_end- 1; j = right_end + 1
          do
            ! scan arr from left end until element >= reference is found
            do
              i = i + 1
              if (.not.(arr(i)%absH > reference%absH)) exit
             !if (arr(i) >= reference) exit
            enddo
            ! scan arr from right end until element <= reference is found
            do
              j = j - 1
              if (.not.(arr(j)%absH < reference%absH)) exit
             !if (arr(j) <= reference) exit
            enddo
            if (i < j) then
              ! swap two out-of-order elements
              temp = arr(i); arr(i) = arr(j); arr(j) = temp
            elseif (i == j) then
              i = i + 1
              exit
            else
              exit
            endif
          enddo
          if (left_end<j)  call quick_sort(left_end, j)
          if (i<right_end)  call quick_sort(i, right_end)
        endif

      end subroutine quick_sort

      subroutine interchange_sort(left_end, right_end)
        integer, intent(in) :: left_end, right_end
        integer             :: i, j
        type(rs_absH)          :: temp

        do i = left_end, right_end - 1
          do j = i+1, right_end
            if (arr(i)%absH < arr(j)%absH) then
              temp = arr(i); arr(i) = arr(j); arr(j) = temp
            endif
          enddo
        enddo
      end subroutine interchange_sort

  end subroutine sort_rs_absH

  subroutine get_n_minus_1_configs(n_det,nelec,dets,nm1_dets)
  ! Get all the ways of removing 1 electron from the input configurations
  ! Useful for Sandeep's partial connections
  ! A Holmes, 10 Nov 2016

    use more_tools, only : get_occ_orbs

    integer(i8b),intent(in) :: n_det
    integer,intent(in) :: nelec
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets(:)
    type(ik_vec),intent(out) :: nm1_dets(:)
#else
    integer(ik),intent(in) :: dets(:)
    integer(ik),intent(out) :: nm1_dets(:)
#endif

    integer :: i_det,i_el
    integer,allocatable :: occ(:)

    allocate(occ(nelec))

    do i_det = 1,n_det
      call get_occ_orbs(dets(i_det),occ)
      do i_el = 1,nelec
        nm1_dets(nelec*(i_det-1)+i_el) = ibclr(dets(i_det),occ(i_el)-1)
      enddo ! i_el
    enddo ! i_det

  end subroutine get_n_minus_1_configs


  subroutine get_connected_dets_in_list(index,det_up,det_dn,dets_up,dets_dn,beta,beta_ind,alpha_m1,alpha_m1_ind,n_connected_dets,connected_indices,connected_matrix_elements,is_included)
  ! Use Sandeep's partial connections to get the matrix elements in the list represented by beta and alpha_m1
  ! (in HCI that list is the list of variational dets)
  ! A Holmes, 10 Nov 2016

    use tools, only : sort_and_merge
    use more_tools, only : get_occ_orbs, binary_search_single, &
        & binary_search_lbound, binary_search_rbound

    integer,intent(in) :: index
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
    type(ik_vec),intent(in) :: dets_up(:), dets_dn(:)
    type(ik_vec),intent(in) :: beta(:), alpha_m1(:)
    type(ik_vec) :: tmp
#else
    integer(ik),intent(in) :: det_up,det_dn
    integer(ik),intent(in) :: dets_up(:), dets_dn(:)
    integer(ik),intent(in) :: beta(:), alpha_m1(:)
    integer(ik) :: tmp
#endif
    integer,intent(in) :: beta_ind(:), alpha_m1_ind(:)
    integer,intent(out) :: n_connected_dets
    integer,intent(out) :: connected_indices(:)
    logical,intent(inout) :: is_included(:) ! Always changed back upon return.
    real(rk),intent(out) :: connected_matrix_elements(:)
    real(rk) :: elem

    integer :: i_el,i,j,k
    integer :: n_beta, n_alpha_m1
    integer :: left, right

    n_beta = size(beta)
    n_alpha_m1 = size(alpha_m1)

    ! Diagonal element first
    n_connected_dets = 1
    call hamiltonian(det_up,det_dn,det_up,det_dn,elem)
    connected_indices(1) = index
    connected_matrix_elements(1) = elem

    ! Compute all alpha excitations: iterate over all pairs of dets in each dets_dn(i)
    call binary_search_lbound(det_dn, beta, left)
    call binary_search_rbound(det_dn, beta, right)
    k = left
    if (k.ne.0) then ! It could only be 0 if part of the sparse_ham had been computed on a previous iteration
      do i = left, right
        j = beta_ind(i)
        if (j > index) then
          call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
          if (abs(elem).gt.1.d-12) then
            n_connected_dets = n_connected_dets + 1
            connected_indices(n_connected_dets) = j
            connected_matrix_elements(n_connected_dets) = elem
            is_included(j) = .true.
          endif
        endif
      enddo
    endif ! k.ne.0

    ! Compute all other excitations: iterate over all pairs of dets in each alpha_m1
    call get_occ_orbs(det_up,occ_up)
    do i_el=1,nup
      tmp = ibclr(det_up,occ_up(i_el)-1)
      call binary_search_lbound(tmp, alpha_m1, left)
      call binary_search_rbound(tmp, alpha_m1, right)
      k = left
      if (k.ne.0) then ! It could only be 0 if part of the sparse_ham had been computed on a previous iteration
        do i = left, right
          j = alpha_m1_ind(i)
          if (j > index .and. (.not. is_included(j))) then
            call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
            if (abs(elem).gt.1.d-12) then
              n_connected_dets = n_connected_dets + 1
              connected_indices(n_connected_dets) = j
              connected_matrix_elements(n_connected_dets) = elem
              is_included(j) = .true.
            endif
          endif
        enddo
      endif ! k.ne.0
    enddo

    if (time_sym) then ! Do the same thing, but with spins flipped

      ! Compute all alpha excitations: iterate over all pairs of dets in each dets_dn(i)
      call binary_search_lbound(det_up, beta, left)
      call binary_search_rbound(det_up, beta, right)
      k = left
      if (k.ne.0) then ! It could only be 0 if part of the sparse_ham had been computed on a previous iteration
        do i = left, right
          j = beta_ind(i)
          if (j > index .and. (.not. is_included(j))) then
         !if (j > index) then
            call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
            if (abs(elem).gt.1.d-12) then
              n_connected_dets = n_connected_dets + 1
              connected_indices(n_connected_dets) = j
              connected_matrix_elements(n_connected_dets) = elem
              is_included(j) = .true.
            endif
          endif
        enddo
      endif ! k.ne.0
  
      ! Compute all other excitations: iterate over all pairs of dets in each alpha_m1
      call get_occ_orbs(det_dn,occ_dn)
      do i_el=1,ndn
        tmp = ibclr(det_dn,occ_dn(i_el)-1)
        call binary_search_lbound(tmp, alpha_m1, left)
        call binary_search_rbound(tmp, alpha_m1, right)
        k = left
        if (k.ne.0) then ! It could only be 0 if part of the sparse_ham had been computed on a previous iteration
          do i = left, right
            j = alpha_m1_ind(i)
            if (j > index .and. (.not. is_included(j))) then
              call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
              if (abs(elem).gt.1.d-12) then
                n_connected_dets = n_connected_dets + 1
                connected_indices(n_connected_dets) = j
                connected_matrix_elements(n_connected_dets) = elem
                is_included(j) = .true.
              endif
            endif
          enddo
        endif ! k.ne.0
      enddo

    endif ! time_sym

    ! Now, remove duplicates by sorting and merging
    call sort_and_merge(n_connected_dets,connected_indices,connected_matrix_elements)

    ! TODO: Only compute matrix element after merge, to avoid computing expensive single excitations and diagonal elements multiple times
   !do i=1,n_connected_dets
   !  call hamiltonian(det_up,det_dn,dets_up(connected_indices(i)),dets_dn(connected_indices(i)),elem)
   !  connected_matrix_elements(i) = elem
   !enddo
    !       Also, compute matrix elements in O(N_el) time rather than O(N_el^2) time

  end subroutine get_connected_dets_in_list

    subroutine get_connected_dets_in_list_hash(index,det_up,det_dn,dets_up,dets_dn,beta_hash,alpha_m1_hash,n_connected_dets,connected_indices,connected_matrix_elements,is_included)
  ! Use Sandeep's partial connections to get the matrix elements in the list represented by beta and alpha_m1
  ! (in HCI that list is the list of variational dets)
  ! A Holmes, 10 Nov 2016
    use fhash_module__ik_int_list
    use tools, only : sort_and_merge
    use more_tools, only : get_occ_orbs, binary_search_single, &
        & binary_search_lbound, binary_search_rbound

    integer,intent(in) :: index
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
    type(ik_vec),intent(in) :: dets_up(:), dets_dn(:)
    type(ik_vec) :: tmp
#else
    integer(ik),intent(in) :: det_up,det_dn
    integer(ik),intent(in) :: dets_up(:), dets_dn(:)
    integer(ik) :: tmp
#endif
    type(fhash_type__ik_int_list),intent(inout) :: beta_hash, alpha_m1_hash
    type(int_vec) :: det_indices
    integer,intent(out) :: n_connected_dets
    integer,intent(out) :: connected_indices(:)
    logical,intent(inout) :: is_included(:) ! Always changed back upon return.
    real(rk),intent(out) :: connected_matrix_elements(:)
    real(rk) :: elem

    integer :: i_el,i,j,k

    ! Diagonal element first
    n_connected_dets = 1
    call hamiltonian(det_up,det_dn,det_up,det_dn,elem)
    connected_indices(1) = index
    connected_matrix_elements(1) = elem

    ! Compute all alpha excitations: iterate over all pairs of dets in each dets_dn(i)
    call beta_hash%get(det_dn,det_indices)

    do i=1,det_indices%n
      j = det_indices%list(i)
      if (j > index) then
        call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
        if (abs(elem).gt.1.d-12) then
          n_connected_dets = n_connected_dets + 1
          connected_indices(n_connected_dets) = j
          connected_matrix_elements(n_connected_dets) = elem
          is_included(j) = .true.
        endif
      endif
    enddo

    ! Compute all other excitations: iterate over all pairs of dets in each alpha_m1
    call get_occ_orbs(det_up,occ_up)
    do i_el=1,nup
      tmp = ibclr(det_up,occ_up(i_el)-1)
      call alpha_m1_hash%get(tmp,det_indices)
      do i = 1,det_indices%n
        j = det_indices%list(i)
        if (j > index .and. (.not. is_included(j))) then
          call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
          if (abs(elem).gt.1.d-12) then
            n_connected_dets = n_connected_dets + 1
            connected_indices(n_connected_dets) = j
            connected_matrix_elements(n_connected_dets) = elem
            is_included(j) = .true.
          endif
        endif
      enddo
    enddo

    if (time_sym) then ! Do the same thing, but with spins flipped

      ! Compute all alpha excitations: iterate over all pairs of dets in each dets_dn(i)
      call beta_hash%get(det_up,det_indices)

      do i = 1,det_indices%n
        j = det_indices%list(i)
        if (j > index .and. (.not. is_included(j))) then
          !if (j > index) then
          call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
          if (abs(elem).gt.1.d-12) then
            n_connected_dets = n_connected_dets + 1
            connected_indices(n_connected_dets) = j
            connected_matrix_elements(n_connected_dets) = elem
            is_included(j) = .true.
          endif
        endif
      enddo

      ! Compute all other excitations: iterate over all pairs of dets in each alpha_m1
      call get_occ_orbs(det_dn,occ_dn)
      do i_el=1,ndn
        tmp = ibclr(det_dn,occ_dn(i_el)-1)
        call alpha_m1_hash%get(tmp,det_indices)
        do i = 1,det_indices%n
          j = det_indices%list(i)
          if (j > index .and. (.not. is_included(j))) then
            call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
            if (abs(elem).gt.1.d-12) then
              n_connected_dets = n_connected_dets + 1
              connected_indices(n_connected_dets) = j
              connected_matrix_elements(n_connected_dets) = elem
              is_included(j) = .true.
            endif
          endif
        enddo 
      enddo

    endif ! time_sym

    ! Now, remove duplicates by sorting and merging
    call sort_and_merge(n_connected_dets,connected_indices,connected_matrix_elements)

    ! TODO: Only compute matrix element after merge, to avoid computing expensive single excitations and diagonal elements multiple times
   !do i=1,n_connected_dets
   !  call hamiltonian(det_up,det_dn,dets_up(connected_indices(i)),dets_dn(connected_indices(i)),elem)
   !  connected_matrix_elements(i) = elem
   !enddo
    !       Also, compute matrix elements in O(N_el) time rather than O(N_el^2) time

  end subroutine get_connected_dets_in_list_hash

  subroutine get_connected_dets_in_list_hash_new(index,det_up,det_dn,dets_up,dets_dn,beta_hash,alpha_hash,alpha_prime_hash,n_connected_dets,connected_indices,connected_matrix_elements,is_included)
  ! Use Sandeep's partial connections to get the matrix elements in the list represented by beta and alpha_m1
  ! (in HCI that list is the list of variational dets)
  ! A Holmes, 10 Nov 2016
    use fhash_module__ik_int_list
    use tools, only : sort_and_merge
    use more_tools, only : get_occ_orbs, binary_search_single, &
        & binary_search_lbound, binary_search_rbound

    integer,intent(in) :: index
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det_up,det_dn
    type(ik_vec),intent(in) :: dets_up(:), dets_dn(:)
    type(ik_vec) :: tmp
#else
    integer(ik),intent(in) :: det_up,det_dn
    integer(ik),intent(in) :: dets_up(:), dets_dn(:)
    integer(ik) :: tmp
#endif
    type(fhash_type__ik_int_list),intent(inout) :: beta_hash, alpha_hash, alpha_prime_hash
    type(int_vec) :: det_indices
    integer,intent(out) :: n_connected_dets
    integer,intent(out) :: connected_indices(:)
    logical,intent(inout) :: is_included(:) ! Always changed back upon return.
    real(rk),intent(out) :: connected_matrix_elements(:)
    real(rk) :: elem

    integer :: i_el,i,j,k,excite_level

    ! Diagonal element first
    n_connected_dets = 1
    call hamiltonian(det_up,det_dn,det_up,det_dn,elem)
    connected_indices(1) = index
    connected_matrix_elements(1) = elem

    ! Compute all beta excitations
    call alpha_hash%get(det_up,det_indices)
    if (index.eq.2947) write(6,*) 'det_indices',det_indices%n
    do i=1,det_indices%n
      j = det_indices%list(i)
      if(index.eq.2947) then
        write(6,*) 'j,index',j,index
      endif
      if (j > index) then
        call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
        if (abs(elem).gt.1.d-12) then
          n_connected_dets = n_connected_dets + 1
          connected_indices(n_connected_dets) = j
          connected_matrix_elements(n_connected_dets) = elem
          is_included(j) = .true.
        endif
      endif
    enddo

    ! Compute all alpha excitations
    call beta_hash%get(det_dn,det_indices)

    do i=1,det_indices%n
      j = det_indices%list(i)
      if (j > index) then
        call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
        if (abs(elem).gt.1.d-12) then
          n_connected_dets = n_connected_dets + 1
          connected_indices(n_connected_dets) = j
          connected_matrix_elements(n_connected_dets) = elem
          is_included(j) = .true.
        endif
      endif
    enddo

    ! Compute all mixed doubles
    call alpha_prime_hash%get(det_up,det_indices)
    if (index.eq.2947) write(6,*) 'det_indices_prime',det_indices%n
    do i=1,det_indices%n
      j = det_indices%list(i)
      if (j > index) then
        ! MJO Check if beta strings are singly connected
        excite_level = popcnt(iand(det_dn,not(dets_dn(j))))
        if (excite_level.eq.1) then
          call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
          if (abs(elem).gt.1.d-12) then
            n_connected_dets = n_connected_dets + 1
            connected_indices(n_connected_dets) = j
            connected_matrix_elements(n_connected_dets) = elem
            is_included(j) = .true.
          endif
        endif
      endif
    enddo

    if (time_sym) then ! Do the same thing, but with spins flipped

      ! Compute all beta excitations
      call alpha_hash%get(det_dn,det_indices)

      do i=1,det_indices%n
        j = det_indices%list(i)
        if (j > index) then
          call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
          if (abs(elem).gt.1.d-12) then
            n_connected_dets = n_connected_dets + 1
            connected_indices(n_connected_dets) = j
            connected_matrix_elements(n_connected_dets) = elem
            is_included(j) = .true.
          endif
        endif
      enddo

      ! Compute all alpha excitations
      call beta_hash%get(det_up,det_indices)

      do i=1,det_indices%n
        j = det_indices%list(i)
        if (j > index) then
          call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
          if (abs(elem).gt.1.d-12) then
            n_connected_dets = n_connected_dets + 1
            connected_indices(n_connected_dets) = j
            connected_matrix_elements(n_connected_dets) = elem
            is_included(j) = .true.
          endif
        endif
      enddo
      ! Compute all mixed doubles
      call alpha_prime_hash%get(det_dn,det_indices)
      do i=1,det_indices%n
        j = det_indices%list(i)
        if (j > index) then
          ! MJO check if alpha strings are singly connected
          excite_level = popcnt(iand(det_up,not(dets_up(j))))
          if (excite_level.eq.1) then
            call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
            if (abs(elem).gt.1.d-12) then
              n_connected_dets = n_connected_dets + 1
              connected_indices(n_connected_dets) = j
              connected_matrix_elements(n_connected_dets) = elem
              is_included(j) = .true.
            endif
          endif
        endif
        enddo

    endif ! time_sym

    ! Now, remove duplicates by sorting and merging
    ! MJO Should not be necessary
    call sort_and_merge(n_connected_dets,connected_indices,connected_matrix_elements)

    ! TODO: Only compute matrix element after merge, to avoid computing expensive single excitations and diagonal elements multiple times
   !do i=1,n_connected_dets
   !  call hamiltonian(det_up,det_dn,dets_up(connected_indices(i)),dets_dn(connected_indices(i)),elem)
   !  connected_matrix_elements(i) = elem
   !enddo
    !       Also, compute matrix elements in O(N_el) time rather than O(N_el^2) time

  end subroutine get_connected_dets_in_list_hash_new


  subroutine hamiltonian(up1,dn1,up2,dn2,elem)
  ! Get Hamiltonian matrix element
  ! Checks for time-reversal symmetry automatically
  ! A Holmes, 10 Nov 2016
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: up1,dn1,up2,dn2
#else
  integer(ik),intent(in) :: up1,dn1,up2,dn2
#endif
  real(rk),intent(out) :: elem
  integer :: excite_level

    elem = 0._rk
    if (time_sym) then
      call hamiltonian_chem_time_sym(up1, dn1, up2, dn2, elem)
    else
      call excitation_level(up1, dn1, up2, dn2, excite_level)
      if (excite_level >= 0) call hamiltonian_chem(up1, dn1, up2, dn2, excite_level, elem)
    endif

  end subroutine hamiltonian


!  subroutine update_index_mapping(n_old,old_up,old_dn,n_new,new_up,new_dn,mapping)
!  ! Update the mapping of sorted index to stored index, where:
!  ! - Sorted index means the index that the determinant would have in the new det list (includes the old dets)
!  !   if it were sorted by label
!  ! - Stored index means the index that represents how the rows of the sparse Hamiltonian are stored,
!  !   i.e., all the old index labels come first, since they were already stored
!  ! So, mapping(sorted_index) = stored_index
!  ! Assumes that before calling subroutine, mapping is correct for the set of old dets
!  ! After calling subroutine, mapping will be correct for the new dets
!  ! A Holmes, 11 Nov 2016
!  integer,intent(in) :: n_old,n_new
!#ifdef NUM_ORBITALS_GT_127
!  type(ik_vec),intent(in) :: old_up(:),old_dn(:),new_up(:),new_dn(:)
!  type(ik_vec),allocatable :: sorted_old_up(:),sorted_old_dn(:),sorted_new_up(:),sorted_new_dn(:)
!#else
!  integer(ik),intent(in) :: old_up(:),old_dn(:),new_up(:),new_dn(:)
!  integer(ik),allocatable :: sorted_old_up(:),sorted_old_dn(:),sorted_new_up(:),sorted_new_dn(:)
!#endif
!  integer,allocatable,intent(inout) :: mapping(:)
!  integer :: i,j,n
!  integer,allocatable :: old_mapping(:)
!
!    allocate(old_mapping(n_old))
!    old_mapping(:) = mapping(:)
!    deallocate(mapping)
!    allocate(mapping(n_new))
!
!    ! Sort old and new dets by label
!    allocate(sorted_old_up(n_old))
!    allocate(sorted_old_dn(n_old))
!    allocate(sorted_new_up(n_new))
!    allocate(sorted_new_dn(n_new))
!
!    sorted_old_up(:) = old_up(:)
!    sorted_old_dn(:) = old_dn(:)
!    sorted_new_up(:) = new_up(:)
!    sorted_new_dn(:) = new_dn(:)
!
!    call sort_dets(n_old,sorted_old_up,sorted_old_dn)
!    call sort_dets(n_new,sorted_new_up,sorted_new_dn)
!
!    ! Iterate simultaneously through sorted old and new,
!    ! if same det, new mapping = old mapping value
!    ! else, it's a det that wasn't there before, so increment n and use mapping = n
!    n = n_old
!    j = 1
!    do i=1,n_old
!      do while (sorted_old_up(i)>sorted_new_up(j).or.(sorted_old_up(i)==sorted_new_up(j).and.sorted_old_dn(i)>sorted_new_dn(j)))
!        j = j+1
!        if (j>n_new)  exit
!        n = n+1
!        mapping(j) = n
!      enddo
!      if (j>n_new)  exit
!      if (sorted_old_up(i)==sorted_new_up(j).and.sorted_old_dn(i)==sorted_new_dn(j)) then
!        mapping(j) = old_mapping(i)
!        n = n-1
!      endif
!    enddo
!
!  end subroutine update_index_mapping


   subroutine auto_assign_hci0_occs(nup, ndn, norb, hci0_sym, hci0_up, hci0_dn)
   ! Automatically assigns the occupations to the HCI starting det (or
   ! linear combination of dets for time_sym) in order to get the input
   ! symmetry hci0_sym
   ! Does this by first assigning 2 electrons to each of the lowest-energy
   ! orbitals, then finding single or double excitation with right symmetry
   ! and lowest energy
   ! A Holmes, 28 Apr 2017

   use mpi_routines, only : mpi_stop,mpi_barr
   use common_selected_ci, only : up,dn

   integer,intent(in) :: nup,ndn,norb,hci0_sym
  !real(rk),intent(in) :: orbital_energies(:)
#ifdef NUM_ORBITALS_GT_127
   type(ik_vec),intent(inout) :: hci0_up,hci0_dn
   type(ik_vec) :: det_up,det_dn
#else
   integer(ik),intent(inout) :: hci0_up,hci0_dn
   integer(ik) :: det_up,det_dn
#endif

   integer :: nup_remaining, ndn_remaining, i, j
   integer :: target_sym,sym,new_orb
   real(rk),allocatable :: tmp_orb_energies(:)
   real(rk) :: min_energy,det_energy

   if (hci0_up==0_ik) then ! No input HF det, so guess using 1-body diagonal terms
     hci0_up = 0_ik
     hci0_dn = 0_ik
     allocate(tmp_orb_energies(norb))
     do i=1,norb
       tmp_orb_energies(i) = integral_value(i,i,norb+1,norb+1)
     enddo
     
     ! Start by assigning electrons to lowest-energy orbitals
    
     nup_remaining = nup
     ndn_remaining = ndn
    
     sym = 1
    
     new_orb = 0
     do i=1,max(nup,ndn)
       min_energy = maxval(orbital_energies(:))
       do j=1,norb
         if (tmp_orb_energies(j)<min_energy) then
           new_orb = j
           min_energy = tmp_orb_energies(j)
         endif
       enddo
       write (6,*) "Adding electron(s) to orbital", new_orb; call flush(6)
       tmp_orb_energies(new_orb) = 1e50_rk
       if (nup_remaining>0) then
         hci0_up = ibset(hci0_up,new_orb-1)
         sym = product_table(sym, orbital_symmetries(new_orb))
         nup_remaining = nup_remaining - 1
       endif
       if (ndn_remaining>0) then
         hci0_dn = ibset(hci0_dn,new_orb-1)
         sym = product_table(sym, orbital_symmetries(new_orb))
         ndn_remaining = ndn_remaining - 1
       endif
       if (nup_remaining==0.and.ndn_remaining==0)  exit
     enddo
  
     write (6,*) "After assigning electrons to lowest-energy orbitals, HF det=",hci0_up,hci0_dn; call flush(6)
   else
     sym = det_sym(hci0_up,hci0_dn)
   endif
   if (time_sym) then
     call hamiltonian_chem_time_sym(hci0_up,hci0_dn,hci0_up,hci0_dn,det_energy)
   else
     call hamiltonian_chem(hci0_up,hci0_dn,hci0_up,hci0_dn,0,det_energy)
   endif
   write (6,*) "Initial starting det energy=",det_energy; call flush(6)
   write (6,*) "HF symmetry is",sym; call flush(6) 

   ! Now, using this HF det as a reference, take the determinant from its CISD
   ! wavefunction with lowest energy out of those with correct symmetry
   det_up = 0_ik
   det_dn = 0_ik

   do while (.not.(hci0_up==det_up.and.hci0_dn==det_dn))

     if (det_up.ne.0_ik) then
       hci0_up = det_up
       hci0_dn = det_dn
     endif

     call find_lowest_energy_det_in_cisd(hci0_up,hci0_dn,hci0_sym,det_up,det_dn,det_energy)
     write (6,*) "New starting det is:",det_up,det_dn,"with energy", det_energy
 
   enddo

  end subroutine auto_assign_hci0_occs


  subroutine find_lowest_energy_det_in_cisd(hf_up, hf_dn, sym, det_up, det_dn, det_energy)
  ! Returns the lowest-energy det of the specified symmetry from the CISD
  ! wavefunction.
  ! If there are several degenerate dets, this picks from among them arbitrarily
  ! (FIXME)
  ! A Holmes, 4 May 2017

  use common_walk, only : n_connected_dets
  use common_run, only : connected_dets_up,connected_dets_dn

#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: hf_up,hf_dn
  type(ik_vec),intent(out) :: det_up,det_dn
  type(ik_vec) :: tmp
#else
  integer(ik),intent(in) :: hf_up,hf_dn
  integer(ik),intent(out) :: det_up,det_dn
  integer(ik) :: tmp
#endif
  integer,intent(in) :: sym
  real(rk),intent(out) :: det_energy
  real(rk) :: tmp_energy,off_diag_elem
  integer :: i,j

    write (6,*) "Finding lowest energy det in CISD with symmetry",sym; call flush(6)

    det_energy = 1e50_rk

    call find_connected_dets_chem(hf_up, hf_dn, n_connected_dets, connected_dets_up, connected_dets_dn, norb, connected_matrix_elements, sym=sym)
   !call find_connected_dets_chem(hf_up, hf_dn, n_connected_dets, connected_dets_up, connected_dets_dn, norb, connected_matrix_elements) ! This doesn't work because it preserves symmetry

    do i=1,n_connected_dets ! Including HF det

      ! If wrong symmetry, go to next
      if (det_sym(connected_dets_up(i),connected_dets_dn(i)).ne.sym)  cycle

      ! Compute energy of det
      if (time_sym) then ! just calculate E the slow way for now

        ! First, check to make sure det is not a singlet if triplet state
        ! targeted!
        if (z<0) then
          if (connected_dets_up(i)==connected_dets_dn(i))  cycle
        endif
 
        call hamiltonian_chem_time_sym(connected_dets_up(i),connected_dets_dn(i),connected_dets_up(i),connected_dets_dn(i),tmp_energy)

      else

        call hamiltonian_chem(connected_dets_up(i),connected_dets_dn(i),connected_dets_up(i),connected_dets_dn(i),0,tmp_energy)
 
      endif

      ! If det is new lowest energy, keep it
      if (tmp_energy<det_energy) then
        det_energy = tmp_energy
        det_up = connected_dets_up(i)
        det_dn = connected_dets_dn(i)
        off_diag_elem = connected_matrix_elements(i)
      endif

    enddo ! i

    write (6,*) "Starting det connected to input guess by matrix element",off_diag_elem; call flush(6)

  end subroutine find_lowest_energy_det_in_cisd


  function det_sym(det_up,det_dn)
  ! Returns the symmetry of the det
  ! A Holmes, 4 May 2017

  use more_tools, only : get_occ_orbs

#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: det_up,det_dn
#else
  integer(ik),intent(in) :: det_up,det_dn
#endif

  integer :: det_sym,i

    call get_occ_orbs(det_up,det_dn,occ_up,occ_dn)

    det_sym = 1
    do i=1,nup
      det_sym = product_table(det_sym, orbital_symmetries(occ_up(i)))
    enddo
    do i=1,ndn
      det_sym = product_table(det_sym, orbital_symmetries(occ_dn(i)))
    enddo

  end function det_sym


  subroutine assign_hf_occs_by_irrep(n_irrep,irreps,irrep_occs_up,irrep_occs_dn,orbital_symmetries,hf_up,hf_dn)
  ! Assign the occupied HF orbitals using counts of electrons per irrep
  ! Goes through the orbitals in order, since they will usually be sorted by
  ! occupation number
  ! A Holmes, 16 May 2017

  integer,intent(in) :: n_irrep,irreps(:),irrep_occs_up(:),irrep_occs_dn(:),orbital_symmetries(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(out) :: hf_up,hf_dn
#else
  integer(ik),intent(out) :: hf_up,hf_dn
#endif
  integer :: i,j,k
  integer :: new_orb

    ! Iterate over up electrons
    hf_up = 0_ik
    do i=1,n_irrep
      if (irreps(i).ne.0) then
        new_orb = 0
        do j=1,irrep_occs_up(i)
          do k=new_orb+1,norb
            if (orbital_symmetries(k)==irreps(i)) then
              new_orb = k
              exit
            endif
          enddo ! k
          hf_up = ibset(hf_up,new_orb-1)
        enddo ! j
      endif
    enddo

    ! Iterate over dn electrons
    hf_dn = 0_ik
    do i=1,n_irrep
      if (irreps(i).ne.0) then
        new_orb = 0
        do j=1,irrep_occs_dn(i)
          do k=new_orb+1,norb
            if (orbital_symmetries(k)==irreps(i)) then
              new_orb = k
              exit
            endif
          enddo ! k
          hf_dn = ibset(hf_dn,new_orb-1)
        enddo ! j
      endif
    enddo

  end subroutine assign_hf_occs_by_irrep


end module chemistry
