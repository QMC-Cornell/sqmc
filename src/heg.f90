module heg
  use types, only: i4b, i8b, rk, ik, ik_vec
#ifdef NUM_ORBITALS_GT_127
  use overload
#endif
  use constants, only: pi
  use generic_sort, only : sort
  use tools, only : n_choose_k, random_int, merge_sort2_up_dn
  use more_tools, only     :    signmod, constrained_dets,fermionic_phase, &
                              & print_walker,print_walker2,print_real_matrix,print_int_matrix,print_map_on_square,          &
                              & real_symmetric_diagonalize,real_symmetric_diagonalize_ow_ham,matrix_lanczos,                &
                              & real_sym_gen_eig

  use common_run, only     :  proposal_method, tau,run_type,max_connected_dets,connected_dets_up,connected_dets_dn,connected_matrix_elements, connected_diag_elems_info
  use common_ham, only     : nelec, nup, ndn, diagonalize_ham, hamiltonian_type, norb
  use common_psi_t, only   : ndet_psi_t, dets_up_psi_t, dets_dn_psi_t, cdet_psi_t, e_trial, trial_wf_iters, norb_trial_wf, n_initiators_trial_wf, n_truncate_trial_wf
  implicit none
  save
!  private
  public
  public                    :: read_heg, system_setup_heg, off_diagonal_move_heg, find_connected_dets_heg, energy_pieces_heg, theoretical_hf_heg  !subroutines

  real(rk), parameter :: EPSILON = 1.0e-15_rk

  integer                   :: n_dim, n_sym_uniq_det_trial_wf
  real(rk)                  :: r_s, length_cell
  real(rk)                  :: cutoff_radius, cutoff_radius_trial_wf !in units of 2 pi /L
  real(rk), allocatable     :: k_vectors(:, :)
  real(rk), allocatable :: orbital_symmetries(:)
  integer, allocatable      :: filled_up(:), empty_up(:)
  integer, allocatable      :: filled_dn(:), empty_dn(:)

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


  type rs_absH
    integer :: r,s
    real(rk) :: absH
  end type rs_absH
  type int3_absH
    integer :: int3(3)
    real(rk) :: absH
  end type int3_absH

  ! For same spin:
  ! 1st dimension is the momentum difference between occupied pair (p, q).
  ! 2nd dimension is the momentum difference between p and an unoccupied orbital r.
  type(int3_absH), allocatable, target :: dtm_hb_same_spin(:, :, :, :)

  ! For opposite spin:
  ! 1st dimension is the momentum difference between occupied pair (p, r).
  type(int3_absH), allocatable, target :: dtm_hb_opposite_spin(:)

  integer, allocatable :: dtm_hb_same_spin_pq_ind(:, :, :) ! For each pq pair, the current index of pr.
  logical, allocatable :: dtm_hb_same_spin_visited(:, :, :, :, :, :)
  logical, allocatable :: dtm_hb_opposite_spin_visited(:, :, :)
  integer, allocatable :: k_vectors_rel(:, :) ! In units of 2 PI / L.
  integer, allocatable :: orb_lut(:, :, :) ! Look up table from relative k vectors to orbital id.
  integer :: n_max ! Largest integer that is smaller than the cutoff radius.
  integer :: n_diff ! Number of different k in each dimension.
  integer :: n_diff_offset ! Offset for array index when storing n_diff related items.

  real(rk) :: max_double, max_single
  integer :: max_h_pair
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
  integer :: n_core_orb
  real(rk),allocatable :: Htot_same(:,:)
  real(rk),allocatable :: Htot_opposite(:,:,:)
  real(rk) :: energy_hf, energy_madelung

contains

  !===========================================================================
  subroutine read_heg
    !---------------------------------------------------------------------------
    ! Description : Read input relevant to the HEG.
    !---------------------------------------------------------------------------
    use types, only: rk, ik
    use common_imp, only: semistochastic
    implicit none
    real(rk) kinetic_energy, exchange_energy
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec) :: hf_up,hf_dn
#else
    integer(ik) :: hf_up,hf_dn
#endif


    !read input
  ! open(5, file='input', status='old')
    read(5,*) n_dim
    write(6,'(/,''Homogeneous electron gas in'',i2,'' dimensions'')') n_dim
    read(5,*) r_s
    write(6,'(''r_s='',f8.4)') r_s
    read(5,*) nelec, nup
    ndn=nelec-nup
    write(6,'(''nelec, nup, ndn='',9i5)') nelec, nup, ndn
    read(5,*) cutoff_radius, cutoff_radius_trial_wf, n_sym_uniq_det_trial_wf
    write(6,'(''cutoff_radius, cutoff_radius_trial_wf, n_sym_uniq_det_trial_wf'',2f8.4,i4)') &
 &  cutoff_radius, cutoff_radius_trial_wf, n_sym_uniq_det_trial_wf
  ! read(5,*) tau

  ! close(5)

    call system_setup_heg()
    if(nup.eq.5 .and. ndn.eq.5) call test_connected_dets()
!   if(nup.eq.5 .and. ndn.eq.5) call test_off_diagonal_move_heg()
#ifdef NUM_ORBITALS_GT_127
        hf_up = maskr_vec(nup)
        hf_dn = maskr_vec(ndn)
#else
        hf_up = maskr(nup,ik)
        hf_dn = maskr(ndn,ik)
#endif
    write (6,*) "HF det=", hf_up,hf_dn
    call kinetic(hf_up, hf_dn, kinetic_energy)
    call exchange(hf_up, hf_dn, exchange_energy)
    energy_hf=kinetic_energy+exchange_energy
    write(6,'(/,''For this cell, HF kinetic, exchange, total energies ='',9f14.8)') kinetic_energy, exchange_energy, energy_hf
    if(n_dim.eq.3) then
      call madelung_energy()
      write(6,'(''For this cell, HF energy including Madelung ='',9f14.8)') energy_hf + energy_madelung
    endif
    call theoretical_hf_heg()

    !if (semistochastic) then
        read(5,*) trial_wf_iters
        allocate(norb_trial_wf(max(0,trial_wf_iters)))
        allocate(n_initiators_trial_wf(max(0,trial_wf_iters)))
        allocate(n_truncate_trial_wf(max(0,trial_wf_iters)))
        read(5,*) norb_trial_wf
        read(5,*) n_initiators_trial_wf
        read(5,*) n_truncate_trial_wf
        write (6,*) "trial_wf_iters",trial_wf_iters
        write (6,*) "n_initiators_trial_wf",n_initiators_trial_wf
        write (6,*) "n_truncate_trial_wf",n_truncate_trial_wf
    !endif
    read(5,*) diagonalize_ham
    write(6,'(''diagonalize_ham (heg) ='',i4)') diagonalize_ham
    call flush(6)

  end subroutine read_heg

  !===========================================================================
  subroutine system_setup_heg()
    !---------------------------------------------------------------------------
    ! Description : Define system dependent parameters.
    !               ndn :: Number of down electrons
    !               length_cell :: length of the side of the real space cube
    !               containing the electron gas
    !
    !               Perform various checks
    !               Generate k-points
    !               Generate trial wavefunction
    !
    ! Created     : F. Petruzielo, 1 Nov 2010
    !---------------------------------------------------------------------------
    use common_run, only : n_connected_dets_hf

    implicit none

    !local variables
    real(rk) :: density
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec) :: hf_up,hf_dn
#else
    integer(ik) :: hf_up,hf_dn
#endif

    if (n_dim /= 2 .and. n_dim /= 3) stop 'Only 2D and 3D are implemented'

    ndn = nelec - nup

    !determine length_cell from r_s and n_dim
    if (n_dim == 2) then
       density = 1._rk / (pi * r_s**2)
    elseif (n_dim == 3) then
       density = 3._rk / (4._rk * pi * r_s**3)
    endif
    length_cell = (nelec / density) ** (1._rk / n_dim)

    call generate_k_vectors()
    if (run_type .eq. 'none') call generate_trial_wf()

    ! allocate the connected dets vectors
    max_connected_dets = 1 + nup * (nup - 1) / 2 * (norb - nup) * (norb - nup - 1) / 2 &
                 &         + ndn * (ndn - 1) / 2 * (norb - ndn) * (norb - ndn - 1) / 2 &
                 &         + nup * ndn * (norb - ndn) * (norb - nup)
    allocate(connected_dets_up(max_connected_dets))
    allocate(connected_dets_dn(max_connected_dets))
    allocate(connected_matrix_elements(max_connected_dets))
    allocate(connected_diag_elems_info(max_connected_dets))

    ! also allocate these vectors, which are used to find connected dets
    allocate(occ_up(nup))
    allocate(occ_dn(ndn))
    allocate(filled_up(nup))
    allocate(filled_dn(ndn))
    allocate(empty_up(norb-nup))
    allocate(empty_dn(norb-ndn))

#ifdef NUM_ORBITALS_GT_127
    hf_up = maskr_vec(nup)
    hf_dn = maskr_vec(ndn)
#else
    hf_up = maskr(nup,ik)
    hf_dn = maskr(ndn,ik)
#endif

    call find_connected_dets_heg(hf_up,hf_dn,n_connected_dets_hf,connected_dets_up,connected_dets_dn)

  end subroutine system_setup_heg
  !=============================================================

  subroutine setup_efficient_heatbath_heg()
    ! Setup efficient heatbath probabilities
    ! A Holmes, Mar 2015 (moved to its own subroutine on 16 Oct 2015)

    use types, only : i8b
    use more_tools, only : setup_alias
    use mpi_routines, only:ncores,whoami,master_core,master_core_node,mpi_bsend,shmem_allocate,mpi_barr,mpi_stop

    implicit none

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
    integer :: p,q,r,s,occ_pair_ind,unocc_pair_ind

    integer, allocatable :: diff_pq(:), diff_pr(:) ! Momentum difference in units of 2 PI / L.
    integer, allocatable :: diff_pq_ind(:), diff_pr_ind(:) ! diff_pq/pr + offset.
    integer :: imp_diff_pr_ind ! Caching the index of pr pair.
    type(int3_absH), pointer :: dtm_hb_same_spin_ptr
    type(int3_absH), pointer :: dtm_hb_opposite_spin_ptr

    if (proposal_method.ne.'fast_heatbath')  return

    if (run_type.eq.'hci') then
      call my_second(1, 'allocate dtm_hb')

      n_diff_offset = 2 * n_max + 1 ! Offset from relative k to index, which is larger than 0.
      n_diff = 4 * n_max + 1 ! Number of possible k difference in each dimension.
      allocate(dtm_hb_same_spin(n_diff, n_diff, n_diff, n_diff**3))
      allocate(dtm_hb_same_spin_pq_ind(n_diff, n_diff, n_diff))
      allocate(dtm_hb_same_spin_visited(n_diff, n_diff, n_diff, n_diff, n_diff, n_diff))
      allocate(dtm_hb_opposite_spin(n_diff**3))
      allocate(dtm_hb_opposite_spin_visited(n_diff, n_diff, n_diff))
      dtm_hb_same_spin_pq_ind = 0
      dtm_hb_same_spin_visited = .false.
      dtm_hb_opposite_spin_visited = .false.

      call my_second(2, 'allocate dtm_hb')
      max_double = 0._rk

      allocate(diff_pq(3))
      allocate(diff_pr(3))
      allocate(diff_pq_ind(3))
      allocate(diff_pr_ind(3))

      ! Same spin.
      do p = 1, norb
        write (6, *) 'Processing same spin. p = ', p
        call flush(6)
        do q = p + 1, norb
          ! Get the index to dtm_hb_same_spin and dtm_hb_same_spin_pq_ind.
          diff_pq = k_vectors_rel(q, :) - k_vectors_rel(p, :)
          diff_pq_ind = diff_pq + n_diff_offset

          do r = 1, norb
            s = find_orb_id(k_vectors_rel(p, :) + k_vectors_rel(q, :) - k_vectors_rel(r, :))
            if (s < 0 .or. s < r) then
              cycle
            endif

            ! Get the index to dtm_hb_same_spin.
            diff_pr = k_vectors_rel(r, :) - k_vectors_rel(p, :)
            diff_pr_ind = diff_pr + n_diff_offset

            if (dtm_hb_same_spin_visited( &
                & diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3), &
                & diff_pr_ind(1), diff_pr_ind(2), diff_pr_ind(3))) then
              cycle
            endif

            H_ijkl = abs(double_excitation_matrix_element_no_ref_abs(p, q, r, s))

            dtm_hb_same_spin_visited( &
                & diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3), &
                & diff_pr_ind(1), diff_pr_ind(2), diff_pr_ind(3)) = &
                & .true.

            if (H_ijkl > EPSILON) then
              imp_diff_pr_ind = &
                  & dtm_hb_same_spin_pq_ind(diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3)) + 1
              dtm_hb_same_spin_pq_ind(diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3)) = &
                  & imp_diff_pr_ind
              dtm_hb_same_spin_ptr => &
                  & dtm_hb_same_spin( &
                      & diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3), imp_diff_pr_ind)
              dtm_hb_same_spin_ptr%int3 = diff_pr
              dtm_hb_same_spin_ptr%absH = H_ijkl
            endif
          enddo

          if (imp_diff_pr_ind > 0) then
            ! Sort in decreasing order of absH.
            call sort_int3_absH(dtm_hb_same_spin( &
                & diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3), 1:imp_diff_pr_ind))

            max_double = &
                & max(max_double, dtm_hb_same_spin( &
                    & diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3), 1)%absH)
          else
            dtm_hb_same_spin(diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3), 1)%absH = 0.0_rk
          endif
        enddo
      enddo

      ! Opposite spin.
      imp_diff_pr_ind = 0
      do p = 1, norb
        write (6, *) 'Processing opposite spin. p = ', p
        call flush(6)
        do q = p, norb
          do r = 1, norb
            s = find_orb_id(k_vectors_rel(p, :) + k_vectors_rel(q, :) - k_vectors_rel(r, :))
            if (s < 0) then
              cycle
            endif

            ! Get the index to dtm_hb_opposite_spin.
            diff_pr = k_vectors_rel(r, :) - k_vectors_rel(p, :)
            diff_pr_ind = diff_pr + n_diff_offset

            if (dtm_hb_opposite_spin_visited( &
                & diff_pr_ind(1), diff_pr_ind(2), diff_pr_ind(3))) then
              cycle
            endif

            H_ijkl = abs(double_excitation_matrix_element_no_ref_abs(p, q + norb, r, s + norb))

            dtm_hb_opposite_spin_visited( &
                & diff_pr_ind(1), diff_pr_ind(2), diff_pr_ind(3)) = &
                & .true.

            if (H_ijkl > EPSILON) then
              imp_diff_pr_ind = imp_diff_pr_ind + 1
              dtm_hb_opposite_spin_ptr => dtm_hb_opposite_spin(imp_diff_pr_ind)
              dtm_hb_opposite_spin_ptr%int3 = diff_pr
              dtm_hb_opposite_spin_ptr%absH = H_ijkl
            endif
          enddo

        enddo
      enddo

      if (imp_diff_pr_ind > 0) then
        call sort_int3_absH(dtm_hb_opposite_spin(1:imp_diff_pr_ind))
        max_double = max(max_double, dtm_hb_opposite_spin(1)%absH)
      else
        dtm_hb_opposite_spin(1)%absH = 0.0_rk
      endif

      do i = 1, imp_diff_pr_ind
        print *, dtm_hb_opposite_spin(i)
      enddo

      write (6,*) "Max double excitation magnitude=",max_double

      max_h_pair = norb
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
      n_core_orb = 0
      do i=n_core_orb+1,norb
        if (mod(i,2)==0)  write (6,*) nint(100*float(i-1)/float(norb)),"% done"; call flush(6)
        do j=n_core_orb+1,norb
          do k=n_core_orb+1,norb
            do l=n_core_orb+1,norb

              if (i==j .or. k==l .or. i==k .or. i==l .or. j==k .or. j==l)  cycle
              ! No time-reversal symmetry here because there are no determinants! (only orbitals ijkl)
              call hamiltonian_heg(ibset(ibset(zero_det,(i-1)),(j-1)),zero_det,ibset(ibset(zero_det,(k-1)),(l-1)),zero_det,H_ijkl)

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
                call hamiltonian_heg(ibset(zero_det,i-1),ibset(zero_det,j-1),ibset(zero_det,k-1),ibset(zero_det,l-1),H_ijkl)
                one_orbital_probabilities(i) = one_orbital_probabilities(i) + abs(H_ijkl)

                two_orbital_probabilities(i,norb+j) = two_orbital_probabilities(i,norb+j) + abs(H_ijkl)
                two_orbital_probabilities(norb+i,j) = two_orbital_probabilities(norb+i,j) + abs(H_ijkl)

                three_orbital_probabilities_opposite_spin(i,j,k) = three_orbital_probabilities_opposite_spin(i,j,k) + abs(H_ijkl)
                if (master_core_node)  four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,l)) = real(four_orbital_probabilities_opposite_spin(opposite_index(i,j,k,l)) + abs(H_ijkl))
              endif
              ! i,l have same spin, j,k have same spin
              if (.not.(i==l.or.j==k)) then
                ! No time-reversal symmetry here because there are no determinants! (only orbitals ijkl)
                call hamiltonian_heg(ibset(zero_det,i-1),ibset(zero_det,j-1),ibset(zero_det,l-1),ibset(zero_det,k-1),H_ijkl)
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


  end subroutine setup_efficient_heatbath_heg
  !=============================================================

  subroutine generate_k_vectors()
    !---------------------------------------------------------------------------
    ! Description : Create all k-vectors within a cutoff radius of
    !               cutoff_radius for a n_dim dimensional
    !               real space cube of side length_cell.
    !               Set number of spatial orbitals, norb.
    !
    ! Created     : F. Petruzielo, 26 Oct 2010
    !---------------------------------------------------------------------------
    use types, only : num_words
    implicit none

    !local variables
    integer :: i, j, k, index
    real(rk), allocatable :: values(:), temp(:, :)
    character(len=3) :: fmt
    integer :: k_ind(3) ! Index for orbital look up, equals relative k plus offset.

    n_max = int(cutoff_radius + EPSILON)

    !initially we define all k-vectors within a cube of side 2*cutoff_radius
    !later we will restrict to the sphere of radius cutoff_radius
    allocate(k_vectors(n_dim, (2*n_max + 1)**n_dim))
    allocate(values(2*n_max + 1))

    !define possible k-point values for a single coordinate
    index = 1
    do i = -n_max, n_max
       values(index) = 2 * pi / length_cell * i
       index = index + 1
    enddo

    !define k-vector n-tuples where n=n_dim
    if (n_dim == 3) then
       index = 1
       do i = -n_max, n_max
          do j = -n_max, n_max
             do k = -n_max, n_max
                k_vectors(1, index) = values(i + n_max+1)
                k_vectors(2, index) = values(j + n_max+1)
                k_vectors(3, index) = values(k + n_max+1)
                index = index + 1
             enddo
          enddo
       enddo
    elseif (n_dim == 2) then
       index = 1
       do i = -n_max, n_max
          do j = -n_max, n_max
             k_vectors(1, index) = values(i + n_max+1)
             k_vectors(2, index) = values(j + n_max+1)
             index = index + 1
          enddo
       enddo
    endif

    !sort k-vector by magnitude
    call sort(k_vectors)

    !remove k-points with magnitude larger than cutoff_radius and count number of remaining k-points
    !these remaining k-points are just the spatial orbitals
    norb = 0
    do i = 1, size(k_vectors, 2)
       if (sqrt(sum(k_vectors(:, i)**2)) > 2 * pi / length_cell * cutoff_radius + epsilon)  exit
       norb = norb + 1
    enddo

    if (norb > num_words*(bit_size(1_ik) - 1)) then
       write(6,*) "Integers used in the program cannot represent the requested number of orbitals"
       write(6,'(i5,a,i3)') norb, " > ", num_words*(bit_size(1_ik) - 1)
       call flush(6)
       stop 'Integers used in the program cannot represent the requested number of orbitals. Recompile'
    endif

    !change size of k_vectors to get rid of k-points with too large of a magnitude
    allocate(temp(n_dim, norb))
    temp = k_vectors(:, 1:norb)
    deallocate(k_vectors)
    allocate(k_vectors(n_dim, norb))
    k_vectors = temp
    deallocate(temp)

    ! Construct relative k vectors.
    allocate(k_vectors_rel(norb, 3))
    do i = 1, norb
      do j = 1, 3
        k_vectors_rel(i, j) = nint(k_vectors(j, i) * length_cell / (2 * pi))
      enddo
    enddo

    ! Build reverse lookup table.
    allocate(orb_lut(2 * n_max + 1, 2 * n_max + 1, 2 * n_max + 1))
    orb_lut = -1
    do i = 1, norb
      k_ind = k_vectors_rel(i, :) + n_max + 1
      orb_lut(k_ind(1), k_ind(2), k_ind(3)) = i
    enddo

    write(6,'(''Within cutoff_radius ='', f10.5,'' number of spatial orbitals ='',i4)') cutoff_radius, norb
    write(6,'(''K Points'')')
    write(fmt, '(i1)') n_dim
    do i = 1, norb
      write(6,'(i5,' // trim(fmt) // 'f10.4,f11.5)') i, k_vectors(:, i), sqrt(sum(k_vectors(:, i)**2))
    enddo
   call flush(6)

  end subroutine generate_k_vectors
  !=============================================================

  type(integer) function find_orb_id(k_rel) result(orb_id)
    ! Find orbital id according to relative k.
    ! Return -1 if not found.

    implicit none
    integer, intent(in) :: k_rel(3)
    integer :: k_ind(3)
    integer :: i

    do i = 1, 3
      if (k_rel(i) < -n_max .or. k_rel(i) > n_max) then
        orb_id = -1
        return
      endif
    enddo

    k_ind = k_rel + n_max + 1
    orb_id = orb_lut(k_ind(1), k_ind(2), k_ind(3))

  end function find_orb_id

  !=============================================================

  subroutine find_set_bits(det, set_bits_pos, n_set_bits_in)
    ! Find the positions (count from right) of bits that have been set to 1.

    implicit none

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det
    type(ik_vec) :: det_tmp
#else
    integer(ik), intent(in) :: det
    integer(ik) :: det_tmp
#endif

    integer, allocatable, intent(out) :: set_bits_pos(:)
    integer, optional, intent(in) :: n_set_bits_in

    integer :: n_set_bits
    integer :: i, pos

    if (present(n_set_bits_in)) then
      n_set_bits = n_set_bits_in
    else
      n_set_bits = popcnt(det)
    endif

    allocate(set_bits_pos(n_set_bits))

    det_tmp = det
    do i = 1, n_set_bits
      pos = trailz(det_tmp) + 1
      set_bits_pos(i) = pos
      det_tmp = ibclr(det_tmp, pos - 1)
    enddo
  end subroutine find_set_bits
  !=============================================================

  function get_gamma_exp(det, occ, eor_set_bits, n_eor_bits)
    ! Count number of bits set in det before the positions in eor_set_bits.

    implicit none

    integer :: get_gamma_exp

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det
#else
    integer(ik), intent(in) :: det
#endif

    integer, intent(in) :: occ(:)
    integer, intent(in) :: eor_set_bits(:), n_eor_bits
    integer :: i, ptr
    integer :: orb_id

    get_gamma_exp = 0
    ptr = 0
    do i = 1, n_eor_bits
      orb_id = eor_set_bits(i)
      if (.not. btest(det, orb_id - 1)) then
        cycle
      endif
      do while (occ(ptr + 1) < orb_id)
        ptr = ptr + 1
      enddo
      get_gamma_exp = get_gamma_exp + ptr
    enddo

  end function get_gamma_exp
  !=============================================================

  subroutine hamiltonian_heg(det_i_up, det_i_dn, det_j_up, det_j_dn, matrix_element, abs_only)

    implicit none

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up, det_i_dn, det_j_up, det_j_dn
    type(ik_vec) :: eor_up, eor_dn
#else
    integer(ik), intent(in) :: det_i_up, det_i_dn, det_j_up, det_j_dn
    integer(ik) :: eor_up, eor_dn
#endif
    real(rk), intent(out) :: matrix_element
    logical, intent(in), optional :: abs_only

    integer :: i, j, k, l, m, n
    integer :: p, q, r, s, orb_p, orb_q, orb_r, orb_s, orb_id
    integer, allocatable :: occ_i_up(:), occ_i_dn(:), occ_j_up(:), occ_j_dn(:)
    integer :: n_occ_up, n_occ_dn
    real(rk), allocatable :: momentum_change(:)
    real(rk), allocatable :: momentum_p(:), momentum_q(:), momentum_s(:)
    real(rk) :: potential_energy
    logical :: momentum_p_set, momentum_q_set, momentum_s_set
    integer :: n_eor_up, n_eor_dn
    integer, allocatable :: eor_up_set_bits(:), eor_dn_set_bits(:)
    integer :: gamma_exp

    real(rk) :: FOUR_PI

    FOUR_PI = 4._rk * PI

    matrix_element = 0._rk

    if (det_i_up == det_j_up .and. det_i_dn == det_j_dn) then
      n_occ_up = popcnt(det_i_up)
      n_occ_dn = popcnt(det_i_dn)
      call find_set_bits(det_i_up, occ_i_up, n_occ_up)
      call find_set_bits(det_i_dn, occ_i_dn, n_occ_dn)

      ! One-electron operator.
      do i = 1, n_occ_up
        p = occ_i_up(i)
        matrix_element = matrix_element + sum(k_vectors(:, p)**2) * 0.5_rk
      enddo
      do i = 1, n_occ_dn
        p = occ_i_dn(i)
        matrix_element = matrix_element + sum(k_vectors(:, p)**2) * 0.5_rk
      enddo

      ! Two-electron operator.
      potential_energy = 0._rk
      do i = 1, n_occ_up
        p = occ_i_up(i)
        do j = i + 1, n_occ_up
          q = occ_i_up(j)
          potential_energy = potential_energy &
              & + FOUR_PI / sum((k_vectors(:, p) - k_vectors(:, q))**2)
        enddo
      enddo
      do i = 1, n_occ_dn
        p = occ_i_dn(i)
        do j = i + 1, n_occ_dn
          q = occ_i_dn(j)
          potential_energy = potential_energy &
              & + FOUR_PI / sum((k_vectors(:, p) - k_vectors(:, q))**2)
        enddo
      enddo
      matrix_element = matrix_element - potential_energy / length_cell**3
    else
      ! Only two-electron operator.
      potential_energy = 0._rk

      eor_up = ieor(det_i_up, det_j_up)
      eor_dn = ieor(det_i_dn, det_j_dn)
      n_eor_up = popcnt(eor_up)
      n_eor_dn = popcnt(eor_dn)

      ! Return zero if not double excitation.
      if (n_eor_up + n_eor_dn /= 4) then
        return
      endif

      ! Obtain momentum.
      call find_set_bits(eor_up, eor_up_set_bits, n_eor_up)
      call find_set_bits(eor_dn, eor_dn_set_bits, n_eor_dn)

      allocate(momentum_change(n_dim))
      allocate(momentum_p(n_dim))
      allocate(momentum_q(n_dim))
      allocate(momentum_s(n_dim))
      momentum_change = 0._rk
      momentum_p_set = .false.
      momentum_q_set = .false.
      momentum_s_set = .false.
      do i = 1, n_eor_up
        orb_id = eor_up_set_bits(i)
        if (btest(det_i_up, orb_id - 1)) then

          momentum_change = momentum_change - k_vectors(:, orb_id)
          if (.not. momentum_p_set) then
            orb_p = orb_id
            momentum_p_set = .true.
          endif
        else
          momentum_change = momentum_change + k_vectors(:, orb_id)
          if (.not. momentum_q_set) then
            orb_q = orb_id
            momentum_q_set = .true.
          elseif (.not. momentum_s_set) then
            orb_s = orb_id
            momentum_s_set = .true.
          endif
        endif
      enddo
      do i = 1, n_eor_dn
        orb_id = eor_dn_set_bits(i)

        if (btest(det_i_dn, orb_id - 1)) then
          momentum_change = momentum_change - k_vectors(:, orb_id)
          if (.not. momentum_p_set) then
            orb_p = orb_id
            momentum_p_set = .true.
          endif
        else
          momentum_change = momentum_change + k_vectors(:, orb_id)
          if (.not. momentum_q_set) then
            orb_q = orb_id
            momentum_q_set = .true.
          elseif (.not. momentum_s_set) then
            orb_s = orb_id
            momentum_s_set = .true.
          endif
        endif
      enddo

      ! Return zero if momentum is not conserved.
      if (sum(momentum_change**2) * length_cell**2 > EPSILON) then
        return
      endif

      potential_energy = FOUR_PI / sum((k_vectors(:, orb_p) - k_vectors(:, orb_q))**2)
      if (n_eor_up /= 2) then
        ! Either two up excitation or two down excitation.
        potential_energy = potential_energy &
            & - FOUR_PI / sum((k_vectors(:, orb_p) - k_vectors(:, orb_s))**2)
      endif

      ! Calculate gamma exponential factor.
      if ((.not. present(abs_only)) .or. (abs_only .eqv. .false.)) then
        call find_set_bits(det_i_up, occ_i_up)
        call find_set_bits(det_j_up, occ_j_up)
        call find_set_bits(det_i_dn, occ_i_dn)
        call find_set_bits(det_j_dn, occ_j_dn)

        gamma_exp = get_gamma_exp(det_i_up, occ_i_up, eor_up_set_bits, n_eor_up) &
            & + get_gamma_exp(det_j_up, occ_j_up, eor_up_set_bits, n_eor_up) &
            & + get_gamma_exp(det_i_dn, occ_i_dn, eor_dn_set_bits, n_eor_dn) &
            & + get_gamma_exp(det_j_dn, occ_j_dn, eor_dn_set_bits, n_eor_dn)

        if (iand(gamma_exp, 1) == 1) then
          potential_energy = -potential_energy
        endif
      endif

      matrix_element = potential_energy / length_cell**3
    endif
  end subroutine hamiltonian_heg
  !=============================================================

  subroutine exchange(det_i_up, det_i_dn, exchange_energy)
    !---------------------------------------------------------------------------
    ! Description : Calculate exchange energy for det_i (split into up and down spin)
    !               Note det_i_up and det_i_dn are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !
    ! Created     : F. Petruzielo, 26 Oct 2010
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
    real(rk), intent(out) :: exchange_energy

    !local variables
    integer :: i, j

    exchange_energy = 0._rk

    !handle up spins
    if (n_dim == 2) then
       do i = 1, norb
          if (btest(det_i_up, i - 1)) then
             do j = i + 1, norb
                if (btest(det_i_up, j - 1)) then
                   exchange_energy = exchange_energy - 1._rk / sqrt(sum((k_vectors(:, j) - k_vectors(:,i))**2))
                endif
             enddo
          endif
       enddo
    elseif (n_dim == 3) then
       do i = 1, norb
          if (btest(det_i_up, i - 1)) then
             do j = i + 1, norb
                if (btest(det_i_up, j - 1)) then
                   exchange_energy = exchange_energy - 1._rk / sum((k_vectors(:, j) - k_vectors(:,i))**2)
                endif
             enddo
          endif
       enddo
    endif

    !handle down spins
    if (det_i_dn == det_i_up) then
       !double exchange energy
       exchange_energy = exchange_energy * 2._rk
    elseif (det_i_dn /= 0) then
       if (n_dim == 2) then
          do i = 1, norb
             if (btest(det_i_up, i - 1)) then
                do j = i + 1, norb
                   if (btest(det_i_up, j - 1)) then
                      exchange_energy = exchange_energy - 1._rk / sqrt(sum((k_vectors(:, j) - k_vectors(:,i))**2))
                   endif
                enddo
             endif
          enddo
       elseif (n_dim == 3) then
          do i = 1, norb
             if (btest(det_i_up, i - 1)) then
                do j = i + 1, norb
                   if (btest(det_i_up, j - 1)) then
                      exchange_energy = exchange_energy - 1._rk / sum((k_vectors(:, j) - k_vectors(:,i))**2)
                   endif
                enddo
             endif
          enddo
       endif
    !else system is completely spin polarized
    endif

    !scale by appropriate constant
    if (n_dim == 2) then
       exchange_energy = exchange_energy * 2._rk * pi / length_cell ** 2
    elseif (n_dim == 3) then
       exchange_energy = exchange_energy * 4._rk * pi / length_cell ** 3
    endif

    ! Hard-core Bosons
    !exchange_energy=-exchange_energy

  end subroutine exchange
  !===========================================================================

  !===========================================================================
  subroutine kinetic(det_i_up, det_i_dn, kinetic_energy)
    !---------------------------------------------------------------------------
    ! Description : Calculate kinetic energy for det_i (split into up and down spin)
    !               Note det_i_up and det_i_dn are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !
    ! Created     : F. Petruzielo, 26 Oct 2010
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
    real(rk), intent(out) :: kinetic_energy

    !local variables
    integer :: i

    kinetic_energy = 0._rk

    !handle up spins
    do i = 1, norb
       if (btest(det_i_up, i - 1)) then
          kinetic_energy = kinetic_energy + sum(k_vectors(:, i)**2)
       endif
    enddo

    !handle down spins
    if (det_i_dn == det_i_up) then
       !double kinetic energy
       kinetic_energy = kinetic_energy * 2._rk
    elseif (det_i_dn /= 0) then
       do i = 1, norb
          if (btest(det_i_dn, i - 1)) then
             kinetic_energy = kinetic_energy + sum(k_vectors(:, i)**2)
          endif
       enddo
    !else System is completely spin polarized
    endif

    !KE = 1/2 Sum_i k_i^2
    kinetic_energy = kinetic_energy / 2._rk

  end subroutine kinetic
  !===========================================================================

  !===========================================================================
  subroutine off_diagonal_coulomb(det_i_up, det_i_dn, det_j_up, det_j_dn, off_diagonal_coulomb_energy)
    !---------------------------------------------------------------------------
    ! Description : Calculate off_diagonal_coulomb energy for det_i and det_j(split into up and down spin)
    !               Note det_i_up, det_i_dn, det_j_up, det_j_dn are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !               ***Assume**** that det_i and det_j are related by a double excitation
    !
    ! Created     : F. Petruzielo, 28 Oct 2010
    !---------------------------------------------------------------------------
    implicit none

    !dummy variables
    real(rk), intent(out) :: off_diagonal_coulomb_energy
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(in) :: det_j_up
    type(ik_vec), intent(in) :: det_j_dn
    type(ik_vec) :: unique_i_up, unique_j_up, unique_i_dn, unique_j_dn
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(in) :: det_j_up
    integer(ik), intent(in) :: det_j_dn
    integer(ik) :: unique_i_up, unique_j_up, unique_i_dn, unique_j_dn
#endif

    !local variables
    integer :: i
    integer :: spin_first_i_bit, spin_second_i_bit, spin_first_j_bit, spin_second_j_bit
    integer :: first_i_bit, second_i_bit, first_j_bit, second_j_bit
    integer :: gamma_i_first, gamma_i_second, gamma_j_first, gamma_j_second
    integer :: n_occ
    logical :: found_first, found_second
    integer, parameter :: up = 1, dn = -1

    off_diagonal_coulomb_energy = 0._rk
    unique_i_up = iand(det_i_up, not(det_j_up))
    unique_j_up = iand(det_j_up, not(det_i_up))
    unique_i_dn = iand(det_i_dn, not(det_j_dn))
    unique_j_dn = iand(det_j_dn, not(det_i_dn))

    !find position of unique bits for det_i and calculate gammas
    !assume up det before dn det in ON vector
    n_occ = 0
    found_first = .false.
    found_second = .false.

    !handle up spin det first
    do i = 1, norb

       if (btest(unique_i_up, i - 1)) then
          if (.not. found_first) then
             first_i_bit = i - 1
             spin_first_i_bit = up
             gamma_i_first = (-1) ** n_occ
             found_first = .true.
          else
             second_i_bit = i - 1
             spin_second_i_bit = up
             gamma_i_second = (-1) ** n_occ
             found_second = .true.
             exit
          endif
       endif

       !increment number of occupied orbs for calculating gamma.
       if (btest(det_i_up, i - 1)) then
          n_occ = n_occ + 1
       endif

    enddo

    !handle dn spin det second
    if (.not. found_second) then
       do i = 1, norb

          if (btest(unique_i_dn, i - 1)) then
             if (.not. found_first) then
                first_i_bit = i - 1
                spin_first_i_bit = dn
                gamma_i_first = (-1) ** n_occ
                found_first = .true.
             else
                second_i_bit = i - 1
                spin_second_i_bit = dn
                gamma_i_second = (-1) ** n_occ
                found_second = .true.
                exit
             endif
          endif

          !increment number of occupied for calculating gamma.
          if (btest(det_i_dn, i - 1)) then
             n_occ = n_occ + 1
          endif

       enddo

    endif


    !find position of unique bits for det_j and calculate gammas
    !assume up det before dn det in ON vector
    n_occ = 0
    found_first = .false.
    found_second = .false.

    !handle up spin det first
    do i = 1, norb

       if (btest(unique_j_up, i - 1)) then
          if (.not. found_first) then
             first_j_bit = i - 1
             spin_first_j_bit = up
             gamma_j_first = (-1) ** n_occ
             found_first = .true.
          else
             second_j_bit = i - 1
             spin_second_j_bit = up
             gamma_j_second = (-1) ** n_occ
             found_second = .true.
             exit
          endif
       endif

       !increment number of occupied for calculating gamma.
       if (btest(det_j_up, i - 1)) then
          n_occ = n_occ + 1
       endif

    enddo

    !handle dn spin det second
    if (.not. found_second) then
       do i = 1, norb

          if (btest(unique_j_dn, i - 1)) then
             if (.not. found_first) then
                first_j_bit = i - 1
                spin_first_j_bit = dn
                gamma_j_first = (-1) ** n_occ
                found_first = .true.
             else
                second_j_bit = i - 1
                spin_second_j_bit = dn
                gamma_j_second = (-1) ** n_occ
                found_second = .true.
                exit
             endif
          endif

          !increment number of occupied for calculating gamma.
          if (btest(det_j_dn, i - 1)) then
             n_occ = n_occ + 1
          endif

       enddo

    endif

    !remember bit index starts at 0 while array index starts at 1
    if (n_dim == 2) then
       if (spin_first_i_bit == spin_first_j_bit .and. spin_second_i_bit == spin_second_j_bit) then
          off_diagonal_coulomb_energy = off_diagonal_coulomb_energy + 1._rk / sqrt(sum((k_vectors(:, first_j_bit + 1) - k_vectors(:,first_i_bit + 1))**2))
       endif
       if (spin_first_i_bit == spin_second_j_bit .and. spin_second_i_bit == spin_first_j_bit) then
          off_diagonal_coulomb_energy = off_diagonal_coulomb_energy - 1._rk / sqrt(sum((k_vectors(:, second_j_bit + 1) - k_vectors(:,first_i_bit + 1))**2))
       endif
       off_diagonal_coulomb_energy = off_diagonal_coulomb_energy * 2._rk * pi / length_cell ** 2
    elseif (n_dim == 3) then
       if (spin_first_i_bit == spin_first_j_bit .and. spin_second_i_bit == spin_second_j_bit) then
          off_diagonal_coulomb_energy = off_diagonal_coulomb_energy + 1._rk / sum((k_vectors(:, first_j_bit + 1) - k_vectors(:,first_i_bit + 1))**2)
       endif
       if (spin_first_i_bit == spin_second_j_bit .and. spin_second_i_bit == spin_first_j_bit) then
          off_diagonal_coulomb_energy = off_diagonal_coulomb_energy - 1._rk / sum((k_vectors(:, second_j_bit + 1) - k_vectors(:,first_i_bit + 1))**2) ! Fermions
         !off_diagonal_coulomb_energy = off_diagonal_coulomb_energy + 1._rk / sum((k_vectors(:, second_j_bit + 1) - k_vectors(:,first_i_bit + 1))**2) ! Hard-core Bosons
       endif
       off_diagonal_coulomb_energy = off_diagonal_coulomb_energy * 4._rk * pi / length_cell ** 3
    endif
    off_diagonal_coulomb_energy = off_diagonal_coulomb_energy * gamma_i_first * gamma_i_second * gamma_j_first * gamma_j_second ! Fermions
    !off_diagonal_coulomb_energy = -abs(off_diagonal_coulomb_energy * gamma_i_first * gamma_i_second * gamma_j_first * gamma_j_second) ! alternative for Fermions?
    !off_diagonal_coulomb_energy = -abs(off_diagonal_coulomb_energy)  ! Hard core boson - but all negative elements             ! Hard-core Bosons

  end subroutine off_diagonal_coulomb
  !===========================================================================

  !===========================================================================
  subroutine off_diagonal_move_heg(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j)
    !---------------------------------------------------------------------------
    ! Description : Perform an offdiagonal move in SQMC for heg starting at det_i
    !               which is split into det_i_up and det_i_dn. det_i has weight = +-1.
    !               since the looping over walkers on same determinant is done
    !               in general part of the code. If move fails, for any reason,
    !               weight_j is returned as 0. Otherwise, return new determinant,
    !               det_j, and its weight, H_jj (diagonal element of Hamiltonian).
    !
    ! Created     : F. Petruzielo, 2 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !dummy variables
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up, det_i_dn
    type(ik_vec), intent(out) :: det_j_up, det_j_dn
#else
    integer(ik), intent(in) :: det_i_up, det_i_dn
    integer(ik), intent(out) :: det_j_up, det_j_dn
#endif
    real(rk), intent(out) :: weight_j

    !local variables
    integer :: excite_from_1, excite_from_2, excite_to_1
    integer :: tot_spin_from !in units of 1/2
    integer :: i, i_elec, i_open, j
    real(rk), allocatable :: tot_mom_from(:), mom_to_1(:)
    integer :: success
!   real(rk) :: ran, rannyu, proposal_prob, acceptance_prob, matrix_element
    real(rk) :: proposal_prob, acceptance_prob, matrix_element

    !initialize det_j to det_i
    det_j_up = det_i_up
    det_j_dn = det_i_dn

    !pick two unique electrons to excite at random
    excite_from_1 = random_int(nelec)
    do
       excite_from_2 = random_int(nelec)
       if (excite_from_1 /= excite_from_2) exit
    enddo

    !calculate total spin of the excitation
    tot_spin_from = 0

    if (excite_from_1 > nup) then
       !spin down
       tot_spin_from = tot_spin_from - 1
    else
       tot_spin_from = tot_spin_from + 1
    endif

    if (excite_from_2 > nup) then
       !spin down
       tot_spin_from = tot_spin_from - 1
    else
       tot_spin_from = tot_spin_from + 1
    endif

    !determine total momentum of the electrons that are being excited
    !and remove excited electrons in det_j
    allocate(tot_mom_from(n_dim))
    tot_mom_from = 0._rk
    i_elec = 0

    do i = 1, norb
       if (btest(det_i_up, i - 1)) then
          i_elec = i_elec + 1
          if (i_elec == excite_from_1 .or. i_elec == excite_from_2) then
             tot_mom_from = tot_mom_from + k_vectors(:,i)
             det_j_up = ibclr(det_j_up, i-1)
          endif
       endif
    enddo

    do i = 1, norb
       if (btest(det_i_dn, i - 1)) then
          i_elec = i_elec + 1
          if (i_elec == excite_from_1 .or. i_elec == excite_from_2) then
             tot_mom_from = tot_mom_from + k_vectors(:,i)
             det_j_dn = ibclr(det_j_dn, i-1)
          endif
       endif
    enddo

!   write(6,'(''excite_from_1, excite_from_2, tot_mom_from ='',2i6,9f9.4)') excite_from_1, excite_from_2, tot_mom_from

    allocate(mom_to_1(n_dim))
    mom_to_1 = 0._rk
    i_open = 0

    !find a single spin orbital to excite to

    if (tot_spin_from == 2) then
       !have to excite to spinup orbital
       excite_to_1 = random_int(norb - nup)
       do i = 1, norb
          if (.not. btest(det_i_up, i - 1)) then
             i_open = i_open + 1
             if (i_open == excite_to_1) then
                mom_to_1 = mom_to_1 + k_vectors(:,i)
                det_j_up = ibset(det_j_up, i-1)
             endif
          endif
       enddo

       !try to find a second spinup orbital that conserves momentum
       do i = 1, norb
          if (.not. btest(ior(det_i_up, det_j_up), i - 1)) then
             success = 0
             do j = 1, n_dim
                if (abs(tot_mom_from(j) - (mom_to_1(j) + k_vectors(j,i))) < epsilon) then
                   success = success + 1
                endif
             enddo
             if (success == n_dim) then
                det_j_up = ibset(det_j_up, i-1)
                proposal_prob = 4. / (nelec * (nelec - 1) * (norb - nup))
                exit
             endif
          endif
       enddo

    elseif (tot_spin_from == 0) then
       !can excite to any open spin orbital
       excite_to_1 = random_int(2*norb - nup - ndn)
       if (excite_to_1 <= norb - nup) then
          !excite to spin-up orbital
          do i = 1, norb
             if (.not. btest(det_i_up, i - 1)) then
                i_open = i_open + 1
                if (i_open == excite_to_1) then
                   mom_to_1 = mom_to_1 + k_vectors(:,i)
                   det_j_up = ibset(det_j_up, i-1)
                endif
             endif
          enddo

          !try to find a spin_down orbital that conserves momentum
          do i = 1, norb
             if (.not. btest(ior(det_i_dn, det_j_dn), i - 1)) then
                success = 0
                do j = 1, n_dim
                   if (abs(tot_mom_from(j) - (mom_to_1(j) + k_vectors(j,i))) < epsilon) then
                      success = success + 1
                   endif
                enddo
                if (success == n_dim) then
                   det_j_dn = ibset(det_j_dn, i-1)
                   proposal_prob = 4. / (nelec * (nelec - 1) * (2*norb - nelec))
                   exit
                endif
             endif
          enddo

       else
          !excite to spin-down orbital
          excite_to_1 = excite_to_1 - (norb - nup)
          do i = 1, norb
             if (.not. btest(det_i_dn, i - 1)) then
                i_open = i_open + 1
                if (i_open == excite_to_1) then
                   mom_to_1 = mom_to_1 + k_vectors(:,i)
                   det_j_dn = ibset(det_j_dn, i-1)
                endif
             endif
          enddo

          !try to find a spinup orbital that conserves momentum
          do i = 1, norb
             if (.not. btest(ior(det_i_up, det_j_up), i - 1)) then
                success = 0
                do j = 1, n_dim
                   if (abs(tot_mom_from(j) - (mom_to_1(j) + k_vectors(j,i))) < epsilon) then
                      success = success + 1
                   endif
                enddo
                if (success == n_dim) then
                   det_j_up = ibset(det_j_up, i-1)
                   proposal_prob = 4. / (nelec * (nelec - 1) * (2*norb - nelec))
                   exit
                endif
             endif
          enddo

       endif

    elseif (tot_spin_from == -2) then
       !have to excite to spindn orbital
       excite_to_1 = random_int(norb - ndn)
       do i = 1, norb
          if (.not. btest(det_i_dn, i - 1)) then
             i_open = i_open + 1
             if (i_open == excite_to_1) then
                mom_to_1 = mom_to_1 + k_vectors(:,i)
                det_j_dn = ibset(det_j_dn, i-1)
             endif
          endif
       enddo

       !try to find a second spindn orbital that conserves momentum
       do i = 1, norb
          if (.not. btest(ior(det_i_dn, det_j_dn), i - 1)) then
             success = 0
             do j = 1, n_dim
                if (abs(tot_mom_from(j) - (mom_to_1(j) + k_vectors(j,i))) < epsilon) then
                   success = success + 1
                endif
             enddo
             if (success == n_dim) then
                det_j_dn = ibset(det_j_dn, i-1)
                proposal_prob = 4. / (nelec * (nelec - 1) * (norb - ndn))
                exit
             endif
          endif
       enddo

    endif

    weight_j = 0
   !diagonal_matrix_element = 0._rk
   !      call hamiltonian_heg(det_i_up, det_i_dn, det_i_up, det_i_dn, diagonal_matrix_element) ! for old position
    if (success == n_dim) then
       !proposed a possible move so calculate acceptance
       call hamiltonian_heg(det_i_up, det_i_dn, det_j_up, det_j_dn, matrix_element)
       acceptance_prob = tau * abs(matrix_element) / proposal_prob
!write(6,'(''det_i_up, det_i_dn, det_j_up, det_j_dn, matrix_element, acceptance_prob'',4i5,9f9.4)') det_i_up, det_i_dn, det_j_up, det_j_dn, matrix_element, acceptance_prob
! Warning: for the moment we are not saving the diagonal matrix element, so compute it regardless of whether weight_j.ne.0 or not
!      weight_j = int((acceptance_prob+rannyu())*sign(1._rk,-matrix_element))
       weight_j = acceptance_prob*sign(1._rk,-matrix_element)
!      if (weight_j .ne. 0) then
          !move was accepted so calculate diagonal element of hamiltonian for h_jj
!         call hamiltonian_heg(det_j_up, det_j_dn, det_j_up, det_j_dn, diagonal_matrix_element) ! for new position
!         write(6,'(''det_i_up, det_i_dn, diagonal_matrix_element='',2i4,d12.4)') det_i_up, det_i_dn, diagonal_matrix_element
!      endif

!      weight_j = int(acceptance_prob)
!      call random_number(ran) !to be replaced
!      if (acceptance_prob - int(acceptance_prob) > ran) then
!         weight_j = weight_j + 1
!      endif
!      if (weight_j > 0) then
!         !move was accepted so calculate diagonal element of hamiltonian for h_jj
!         call hamiltonian_heg(det_j_up, det_j_dn, det_j_up, det_j_dn, diagonal_matrix_element)
!         !set sign of walker
!         if (-matrix_element * weight_i < 0) then
!            weight_j = -weight_j
!         endif
!      endif
    endif

  end subroutine off_diagonal_move_heg
  !=============================================================

  !===========================================================================
  subroutine energy_pieces_heg(det_i_up, det_i_dn, e_mix_numerator, e_mix_denominator)
    !---------------------------------------------------------------------------
    ! Description : Calculate pieces of the local energy for det_i
    !               Numerator   is sum_j=1^{ndet_psi_t} H_ij * trial_wf_weight on det_j
    !               Denominator is trial_wf_weight on det_i
    !
    ! Created     : F. Petruzielo, 9 Nov 2010
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

    integer :: j
    real(rk) :: matrix_element

    e_mix_numerator = 0._rk
    e_mix_denominator = 0._rk
    !calculate contributions to the energy
    do j = 1, ndet_psi_t
       if (is_connected_heg(det_i_up, det_i_dn, dets_up_psi_t(j), dets_dn_psi_t(j))) then
          !det_i is connected to a det in the trial wavefunction
          call hamiltonian_heg(det_i_up, det_i_dn, dets_up_psi_t(j), dets_dn_psi_t(j), matrix_element)
          e_mix_numerator = e_mix_numerator + matrix_element * cdet_psi_t(j)
          if (dets_up_psi_t(j) == det_i_up .and. dets_dn_psi_t(j) == det_i_dn) then
             !det_i is in the trial wavefunction
             e_mix_denominator = cdet_psi_t(j)
          endif
       endif
    enddo

  end subroutine energy_pieces_heg
  !=============================================================

  !=============================================================
  subroutine generate_trial_wf()
    !---------------------------------------------------------------------------
    ! Description : Create trial wavefunction.
    !               Diagonalize H in a subset of the entire set of determinants.
    !               The space is specified by allowing all nelec electrons to excite
    !               to all spatial orbitals within cutoff_radius_trial_wf.
    !               (To be updated to allow from CAS and more generally a RAS)
    !               Then pick n_sym_uniq_det_trial_wf symmetry unique determinants.
    !               Note that the total number of determinants will be larger than
    !               this. Finally, diagonalize H in this yet smaller space to generate
    !               the trial wavefunction and weights.
    !
    ! Created     : F. Petruzielo, 5 Nov 2010
    !---------------------------------------------------------------------------
!   use common_psi_t, only: ndet_psi_t, dets_up_psi_t, dets_dn_psi_t, cdet_psi_t, e_trial
    use common_run, only: tau_multiplier
    use common_ham, only: diagonal_ham_lowest, diagonal_ham_highest
    implicit none

    !local variables
    integer :: norb_psi_t, n_init_det_up, n_init_det_dn, n_det
    integer :: count
    integer :: i, j
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec) det_up, det_dn
    type(ik_vec), allocatable :: init_dets_up(:), init_dets_dn(:)
    type(ik_vec), allocatable :: dets_up(:), dets_dn(:), temp_dets(:)
    type(ik_vec) :: temp
#else
    integer(ik) det_up, det_dn
    integer(ik), allocatable :: init_dets_up(:), init_dets_dn(:)
    integer(ik), allocatable :: dets_up(:), dets_dn(:), temp_dets(:)
    integer(ik) :: temp
#endif
!   real(rk) diagonal_ham_lowest, diagonal_ham_highest
    real(rk) :: previous_value
    real(rk) :: lowest_eigenvalue
    real(rk), allocatable :: lowest_eigenvector(:)
    logical :: valid_combination

    ! First we generate all determinants formed with nup spin-up electrons
    ! and ndn spin-down electrons in all spatial orbitals within
    ! a cutoff of cutoff_radius_trial_wf
    norb_psi_t = 0
    do i = 1, norb
       if (sqrt(sum(k_vectors(:, i)**2)) > 2._rk * pi * cutoff_radius_trial_wf / length_cell + epsilon) exit
       norb_psi_t = norb_psi_t + 1
    enddo

    n_init_det_up = int(n_choose_k(norb_psi_t, nup),i4b)
    n_init_det_dn = int(n_choose_k(norb_psi_t, ndn),i4b)

    allocate(init_dets_up(n_init_det_up))
    allocate(init_dets_dn(n_init_det_dn))

    !generate all possible up-spin determinants out of norb_psi_t orbitals and nup up-spin electrons (from Bit Twiddling Hacks)
    init_dets_up(1) = 2_ik ** nup - 1
    do i = 2, n_init_det_up
       temp = ior(init_dets_up(i-1), init_dets_up(i-1) - 1) + 1
       init_dets_up(i) = ior(temp,ishft(iand(temp, -temp) / iand(init_dets_up(i-1),-init_dets_up(i-1)),-1) - 1)
    enddo

    !generate all possible dn-spin determinants out of norb_psi_t orbitals and ndn dn-spin electrons (from Bit Twiddling Hacks)
    init_dets_dn(1) = 2_ik ** ndn - 1
    do i = 2, n_init_det_dn
       temp = ior(init_dets_dn(i-1), init_dets_dn(i-1) - 1) + 1
       init_dets_dn(i) = ior(temp,ishft(iand(temp, -temp) / iand(init_dets_dn(i-1),-init_dets_dn(i-1)),-1) - 1)
    enddo

    !only include determinants which have zero momentum
    n_det = 0
    allocate(dets_up(n_init_det_up * n_init_det_dn))
    allocate(dets_dn(n_init_det_up * n_init_det_dn))
    do i = 1, n_init_det_up
       do j = 1, n_init_det_dn
          call filter_dets(init_dets_up(i), init_dets_dn(j), valid_combination)
          if (valid_combination) then
             n_det = n_det + 1
             dets_up(n_det) = init_dets_up(i)
             dets_dn(n_det) = init_dets_dn(j)
          endif
       enddo
    enddo
    write(6,'(/,i8,'' determinants formed from orbs within cutoff_radius_trial_wf='',f9.4,'' and'',i7,'' have zero momentum'')') n_init_det_up*n_init_det_dn, cutoff_radius_trial_wf, n_det
    call flush(6)

    !cleanup
    deallocate(init_dets_up)
    deallocate(init_dets_dn)
    allocate(temp_dets(n_det))
    temp_dets = dets_up(1:n_det)
    deallocate(dets_up)
    allocate(dets_up(n_det))
    dets_up = temp_dets
    temp_dets = dets_dn(1:n_det)
    deallocate(dets_dn)
    allocate(dets_dn(n_det))
    dets_dn = temp_dets
    deallocate(temp_dets)

    !perform first diagonalization
    n_det = size(dets_up)
    allocate(lowest_eigenvector(n_det))

    call diagonalize_hamiltonian_heg(dets_up, dets_dn, lowest_eigenvector, lowest_eigenvalue)
    e_trial=lowest_eigenvalue ! We set e_trial from the first diagonalization since it is a bit more accurate
    write(6,'(/,''from initial trial wavefn. e_trial='',f12.6)') e_trial
    call flush(6)

    !pick the number of desired unique (by symmetry) determinants
    count = 0
    previous_value = 0._rk
    do i = 1, n_det
       if ( abs(abs(previous_value) - abs(lowest_eigenvector(i))) > epsilon ) then
          previous_value = lowest_eigenvector(i)
          count = count + 1
       endif
       ndet_psi_t = i
       if ( count > n_sym_uniq_det_trial_wf .and. n_sym_uniq_det_trial_wf >= 1 ) then
         ndet_psi_t = i - 1
         count = count - 1
         exit
       endif
    enddo
    if(n_sym_uniq_det_trial_wf >= 1) then
       n_sym_uniq_det_trial_wf=min(n_sym_uniq_det_trial_wf,count)
    else
       n_sym_uniq_det_trial_wf=count
    endif

    allocate(dets_up_psi_t(ndet_psi_t))
    allocate(dets_dn_psi_t(ndet_psi_t))
    dets_up_psi_t = dets_up(1:ndet_psi_t)
    dets_dn_psi_t = dets_dn(1:ndet_psi_t)

    !diagonalize in the yet smaller space to generate trial wavefunction
    allocate(cdet_psi_t(ndet_psi_t))
    call diagonalize_hamiltonian_heg(dets_up_psi_t, dets_dn_psi_t, cdet_psi_t, lowest_eigenvalue)
    write(6,'(/,''From initial trial wavefn. in smaller space of'',i6,'' CSFs'',i6,'' dets., e_trial='',f12.6)') n_sym_uniq_det_trial_wf, ndet_psi_t, lowest_eigenvalue

    if(cdet_psi_t(maxloc(abs(cdet_psi_t),1)).lt.0) cdet_psi_t=-cdet_psi_t
    write(6,'(/,a)') "Trial Wavefunction"
    write(6,'(''Number of determinants in final Psi_T with n_sym_uniq_det_trial_wf='',i3,'' determinants = '', i6)') n_sym_uniq_det_trial_wf, ndet_psi_t
    write(*,*) "    up-det     dn-det   coefficient"
    do i = 1, ndet_psi_t
       write(6,'(2i10,f10.4)') dets_up_psi_t(i), dets_dn_psi_t(i), cdet_psi_t(i)
    enddo
    call flush(6)

! Use diagonal elements of Hamiltonian to set tau.
    det_up=2_ik**nup-1
    det_dn=2_ik**ndn-1
    call hamiltonian_heg(det_up, det_dn, det_up, det_dn, diagonal_ham_lowest)
    det_up=2_ik**norb - 2_ik**(norb-nup)
    det_dn=2_ik**norb - 2_ik**(norb-ndn)
    call hamiltonian_heg(det_up, det_dn, det_up, det_dn, diagonal_ham_highest)

    if(tau.eq.0._rk) then
      tau=tau_multiplier/(diagonal_ham_highest-diagonal_ham_lowest)
      write(6,'(/,''diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau='',9f12.6)') diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau
     else
      write(6,'(/,''diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, deduced tau='',9f12.6)') diagonal_ham_lowest, diagonal_ham_highest, tau_multiplier, tau_multiplier/(diagonal_ham_highest-diagonal_ham_lowest)
      write(6,'(''Input value (actually used) of tau='',9f12.6)') tau
    endif
    call flush(6)

  end subroutine generate_trial_wf
  !=============================================================

  !=============================================================
  subroutine diagonalize_hamiltonian_heg(dets_up, dets_dn, lowest_eigenvector, lowest_eigenvalue)
    !---------------------------------------------------------------------------
    ! Description : Construct and diagonalize H in a space of dets (broken up into dets_up and dets_dn)
    !               Sort ground state eigenvector by decreasing magnitude (HF first for example).
    !               Put dets_up and dets_dn in the same order.
    !
    ! Created     : F. Petruzielo, 6 Nov 2010
    !---------------------------------------------------------------------------
    use common_run, only: tau_multiplier
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
    integer :: n_det
    integer :: len_work, info
    integer :: i, j
    real(rk) tau_optimal_deterministic, tau_optimal_stochastic
    real(rk), allocatable :: ham(:,:)
    real(rk), allocatable :: eigenvalues(:), work(:)
    character(len=2) :: fmt

    n_det = size(dets_up)
    allocate(ham(n_det, n_det))
    allocate(eigenvalues(n_det))

    !construct hamiltonian
    ham = 0._rk
    do i = 1, n_det
       call hamiltonian_heg(dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), ham(i,i))  ! diagonal element
       do j = i+1, n_det
          if (is_connected_heg(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j))) then
             !determinants are connected by a double excitation so evaluate matrix elementx
             call hamiltonian_heg(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), ham(i,j))
          else
             ham(i,j) = 0._rk
          endif
          !symmetric matrix
          ham(j,i) = ham(i,j)
       enddo
    enddo

    !diagonalize with lapack routine
    len_work = 3 * n_det -1
    allocate(work(len_work))

    write(6,'(/,''n_det, len_work='',9i8)') n_det,len_work
    call flush(6)

    ! ASK HITESH:   Here we are only using the lapack routine dsyev, but the chem module uses a lanczos
    !               diagonalization function at this point also.  Does this need to be changed?

! Write Hamiltonian matrix for testing
    if(n_det.le.1000) then
!     open(2, file='hamiltonian_new')
      write(6,'(/,''Hamiltonian matrix'')')
      write(6,'(i5,'' n_det'')') n_det
      do i=1,n_det
        write(6,'(es21.14,999es22.14)') ham(i,:)
      enddo
!     close(2)
    endif

    call dsyev('V', 'U', n_det, ham, n_det, eigenvalues, work, len_work, info)
    if (info /= 0) then
       write(6,*) "info = ", info
       stop 'Diagonalization Failed!'
    endif

    lowest_eigenvector = ham(:,1)
    lowest_eigenvalue = eigenvalues(1)
    call sort(lowest_eigenvector, dets_up, dets_dn)

    write(6,'(/,''Wavefunction has'',i9,'' determinants. Dets. with largest'',i3,'' abs coefs are:'',/,'' det_up   det_dn,  cdet'')')  n_det,min(n_det,20)
    write (fmt, '(i2)') 2*num_words
    write(6,'(' // trim(fmt) // 'b20,' // trim(fmt) // 'i15,f13.9)') (dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), lowest_eigenvector(i),i=1,min(n_det,20))
    write(6,'(/,''Lowest'',i3,'' eigenvalues='',20f12.6)') min(n_det,20),(eigenvalues(i),i=1,min(n_det,20))
    write(6,'(''Lowest, highest eigenvalues='', 9f18.9)') eigenvalues(1),eigenvalues(n_det)

! Calculate what tau_optimal_deterministic and tau_optimal_stochastic would be if the wavefn were this trial wavefn.
    tau_optimal_stochastic=1/(eigenvalues(n_det)-eigenvalues(1))
    if(n_det.ge.2) then
      tau_optimal_deterministic=2/(eigenvalues(n_det)+eigenvalues(2)-2*eigenvalues(1))
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_deterministic, tau_optimal_stochastic='',9f10.6)') &
   &  tau_optimal_deterministic, tau_optimal_stochastic
     else
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_stochastic='',9f10.6)') tau_optimal_stochastic
    endif

  end subroutine diagonalize_hamiltonian_heg
  !=============================================================

  !=============================================================
  subroutine diagonalize_sparse_hamiltonian_heg(dets_up, dets_dn, n_det, lowest_eigenvector, lowest_eigenvalue)
    !---------------------------------------------------------------------------
    ! Description : Construct and diagonalize H in a space of dets (broken up into dets_up and dets_dn)
    !               Sort ground state eigenvector by decreasing magnitude (HF first for example).
    !               Put dets_up and dets_dn in the same order.
    !
    ! Created     : RZ Lamberty, 16 Sep 2012 (same as diagonalize_hamiltonian_heg, but stores H sparsely)
    !---------------------------------------------------------------------------
    use types, only : i8b, num_words
    use more_tools, only : real_symmetric_diagonalize_ow_ham
    use common_run, only: tau_multiplier
    use tools, only : print_excitation_levels_and_wts
    implicit none

    !dummy arguments
    integer,intent(in) :: n_det
    real(rk), intent(out) :: lowest_eigenvector(:)
    real(rk), intent(out) :: lowest_eigenvalue
    !real(rk),intent(out) :: tau_out
    !logical,intent(in)   :: time_sym_on
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(inout) :: dets_up(:)
    type(ik_vec), intent(inout) :: dets_dn(:)
    type(ik_vec),allocatable    :: temp_i16_up(:),temp_i16_dn(:)!,temp_i_2(:)
#else
    integer(ik), intent(inout) :: dets_up(:)
    integer(ik), intent(inout) :: dets_dn(:)
    integer(ik),allocatable    :: temp_i16_up(:),temp_i16_dn(:)!,temp_i_2(:)
#endif

    !local variables
    real(rk), allocatable :: H_values(:)
    integer(i8b), allocatable :: H_indices(:),H_nonzero_elements(:)
    real(rk) :: highest_eigenvalue,second_lowest_eigenvalue
    real(rk) :: tau_optimal_deterministic, tau_optimal_stochastic
    real(rk),allocatable :: ham(:,:),eigenvalues(:)
    real(rk) :: matrix_element
    integer :: i,j
    integer(i8b),parameter :: max_nonzero_elements = int(1e9,i8b)
    integer(i8b) :: num_nonzero_elements
    integer,allocatable        :: iorder(:),temp_i_2(:)
    integer :: n
    character(len=2) :: fmt

    call my_second(1, 'diagonalize_sparse_hamiltonian_heg')

    if (n_det.ge.2000) then
      !construct hamiltonian

      ! ASK HITESH: Is my g_s_h_h function up to snuff?
      call generate_sparse_ham_heg(dets_up,dets_dn,H_indices,H_nonzero_elements,H_values,max_nonzero_elements,num_nonzero_elements)

      write (6,'(''ndet ='',i10,'''')') n_det
      call flush(6)

      call my_second(2, 'generate_sparse_ham_heg called from diagonalize_sparse_ham_heg')

      write (6,'(''Number of nonzero elements in sparse H='',i10,'' ='',es9.2,'' Max number of nonzero elements='',i12,'' ='',es9.2)') num_nonzero_elements, real(num_nonzero_elements), max_nonzero_elements, real(max_nonzero_elements)
      call flush(6)

      if (num_nonzero_elements>max_nonzero_elements) then
        write (6,'(''Cannot store Hamiltonian sparsely for Lanczos, so use matrix_lanczos_on_the_fly'')') ; call flush(6)
        ! sort by label so that a binary search can be performed
        n=size(dets_up)
        allocate(temp_i16_up((n+1)/2))
        allocate(temp_i16_dn((n+1)/2))
        allocate(temp_i_2((n+1)/2))
        allocate(iorder(n) )
        ! ASK HITESH:   Is this function ham type independent?
        call merge_sort2_up_dn(dets_up,dets_dn, iorder, n, temp_i16_up, temp_i16_dn, temp_i_2)
        deallocate(iorder)
        deallocate(temp_i16_up)
        deallocate(temp_i16_dn)
        deallocate(temp_i_2)
        ! ASK HITESH:   Is the version I wrote acceptable?
        call matrix_lanczos_on_the_fly(dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,highest_eigenvalue,second_lowest_eigenvalue)
      else
        write (6,*) "Using sparsely stored Hamiltonian for Lanczos"
        call flush(6)
        !write(6,'(''Number of nonzero elements in sparse Hamiltonian='',i10)') nonzero_elements
        call my_second(2, 'generate_sparse_ham_heg called from diagonalize_sparse_ham_heg')
        call flush(6)
        write (6,'(''Performing Lanczos using matrix_lanczos in more_tools, n_det='',i10)') n_det
        call flush(6)
        ! ASK HITESH:   Is this function ham type independent?  Looks like it is
        call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue)
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
            call hamiltonian_heg(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
            ham(i,j)=matrix_element
          else
            if (is_connected_heg(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j))) then
              !determinants are connected by a double excitation so evaluate matrix elementx
              call hamiltonian_heg(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
              ham(i,j) = matrix_element
              ham(j,i) = matrix_element
            endif
          endif
        enddo
      enddo
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

    call sort(lowest_eigenvector, dets_up, dets_dn) ! Sort in order of decreasing absolute coef
    write(6,'(/,''Wavefunction has'',i9,'' determinants. Dets. with largest'',i3,'' abs coefs are:'',/,''     dets_up             dets_dn        dets_up      dets_dn  excit_level cdet'')')  n_det,min(n_det,20)
    write (fmt, '(i2)') 2*num_words
    write(6,'(' // trim(fmt) // 'b20,' // trim(fmt) // 'i15,i3,f13.9)') (dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), popcnt(iand(dets_up(i),not(dets_up(1))))+popcnt(iand(dets_dn(i),not(dets_dn(1)))), lowest_eigenvector(i),i=1,min(n_det,20))
    write(6,'(/,''Lowest, highest eigenvalue='',9f18.9)') lowest_eigenvalue, highest_eigenvalue

! Calculate number of dets with various excitation levels and the sum of their squared wts.
    !call print_excitation_levels_and_wts(n_det,dets_up,dets_dn,lowest_eigenvector,norb,orbital_symmetries)

! Calculate what tau_optimal_deterministic and tau_optimal_stochastic would be if the wavefn were this trial wavefn.
    tau_optimal_stochastic=1/(highest_eigenvalue-lowest_eigenvalue)
    !tau_out=tau_multiplier*tau_optimal_stochastic
    if(n_det.ge.2) then
      tau_optimal_deterministic=2/(highest_eigenvalue+second_lowest_eigenvalue-2*lowest_eigenvalue)
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_deterministic, tau_optimal_stochastic='',9f10.6)') &
   &  tau_optimal_deterministic, tau_optimal_stochastic
     else
      write(6,'(/,''If the wavefn were this trial wf, tau_optimal_stochastic='',9f10.6)') tau_optimal_stochastic
    endif

    call my_second(2, 'diagonalize_sparse_hamiltonian_heg')

  end subroutine diagonalize_sparse_hamiltonian_heg
  !=============================================================

  ! ASK HITESH: I modelled this after generate_sparse_ham_hubbardk, not chem or chem2.  Is this okay?
  !=============================================================
  subroutine generate_sparse_ham_heg(dets_up,dets_dn,H_indices,H_nonzero_elements,H_values,max_nonzero_elems,num_nonzero_elems)
    !---------------------------------------------------------------------------
    ! Description : Compute the Hamiltonian in a given space and store it sparsely.
    !
    ! Created     : RZ Lamberty, 16 Sep 2012
    ! Modified    : A Holmes, 18 Oct 2012. Use brute force method when it is cheaper
    !               than binary search (using estimate from chemistry; TODO: find
    !               estimates specific to heg).
    !---------------------------------------------------------------------------

    use types, only : i8b
    use more_tools, only : binary_search
    use mpi_routines, only : ncores

    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
#endif
    integer(i8b),allocatable,intent(out) :: H_indices(:),H_nonzero_elements(:)
    real(rk),allocatable,intent(out) :: H_values(:)
    integer(i8b),optional,intent(in) :: max_nonzero_elems
    integer(i8b),optional,intent(out) :: num_nonzero_elems

    integer :: i,j,k,n_det
    integer(i8b) :: nonzero_elements,isparse,itotal
    real(rk) :: matrix_element
    integer,allocatable :: H_rows(:),row_count(:)
    integer :: n_connected_dets
    integer :: connected_indices(max_connected_dets)
    integer :: numops
    logical :: brute_force ! whether or not to use the brute force method (instead of binary search)

    call my_second(1, 'generate_sparse_ham_heg')

    n_det=size(dets_up)

    allocate(H_nonzero_elements(n_det))
    H_nonzero_elements(:)=0

    allocate(row_count(n_det))

    ! Count number of nonzero elements of H
    nonzero_elements=0

    call find_connected_dets_heg(dets_up(1),dets_dn(1),n_connected_dets,connected_dets_up,connected_dets_dn)

    write (6,*) "N_det=",n_det,"N_con=",n_connected_dets
    call flush(6)

    ! Determine which method to use based on n_det and n_connected_dets (in general, larger matrices are faster with binary search)
    ! Note: the factor of 2.1 comes from chemistry tests. This is probably the same factor here, but we should check!
    if (real(n_det)>2.1*real(n_connected_dets)*log(real(n_det))) then
      brute_force = .false.
    else
      brute_force = .true.
    endif

    if (brute_force) then

      do i = 1, n_det
        do j = 1, i
          if (is_connected_heg(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j))) then
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

      do i = 1, n_det
        call find_connected_dets_heg(dets_up(i),dets_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn)
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

    write(6,'(''allocating H_indices etc. arrays of size nonzero_elements='',i10,'' ='',es11.4,'' in generate_sparse_ham_hubbard'')') nonzero_elements, real(nonzero_elements)
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
            if (is_connected_heg(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j))) then
                call hamiltonian_heg(dets_up(i),dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
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

    else

      do i = 1, n_det
        if (i>1) then
          itotal=itotal+H_nonzero_elements(i-1)
          H_rows(isparse+1:itotal)=H_rows(isparse)
          isparse = itotal
        endif

        call find_connected_dets_heg(dets_up(i),dets_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn)

        do k=1,n_connected_dets
          call binary_search(connected_dets_up(k),connected_dets_dn(k),dets_up(1:n_det),dets_dn(1:n_det),j)
          connected_indices(k) = j
        enddo
        ! sort connected_indices
        call sort(n_connected_dets,connected_indices(1:n_connected_dets),numops)
        do k=1,n_connected_dets
          j = connected_indices(n_connected_dets-k+1)
          if (j>0.and.i>=j) then
            call hamiltonian_heg(dets_up(i), dets_dn(i), dets_up(j), dets_dn(j), matrix_element)
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

    call my_second(2, 'generate_sparse_ham_heg')

  end subroutine generate_sparse_ham_heg
  !=============================================================

  ! ASK HITESH: is this okay?
  !===========================================================================
  subroutine matrix_lanczos_on_the_fly(dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,highest_eigenvalue,second_lowest_eigenvalue)
    !---------------------------------------------------------------------------
    ! Created by  : RZ Lamberty, 16 Sep 2012
    !---------------------------------------------------------------------------
    implicit none

    ! Dummy
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
#endif
    real(rk),intent(out) :: lowest_eigenvector(:)
    real(rk),intent(out) :: lowest_eigenvalue
    real(rk),intent(out),optional :: highest_eigenvalue,second_lowest_eigenvalue

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
    real(rk)                   :: epsilon_convergence=1.0e-6  ! for Lanczos
    integer                    :: n

    n=size(dets_up)

    iterations=30          ! User option
    iterations=min(n,iterations)
    allocate (v(n,iterations+1))
    allocate (w(n))
    allocate(alphas(iterations+1))
    allocate(betas(iterations+1))
    w(:)=0._rk
    v(:,1)=0._rk
    v(1,1)=1._rk
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
        if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<epsilon_convergence) then
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
  !===========================================================================

  ! ASK HITESH: is this okay?
  !===========================================================================
  subroutine apply_H_on_the_fly(dets_up,dets_dn,v1,v2)
    !---------------------------------------------------------------------------
    ! Created by  : RZ Lamberty, 16 Sep 2012
    !---------------------------------------------------------------------------

    use more_tools, only : binary_search

    implicit none

    real(rk),intent(in) :: v1(:)
    real(rk),intent(out) :: v2(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#endif

    integer :: i,j,k,n_connected_dets
    real(rk) :: matrix_element

    v2(:)=0._rk
    do i=1,size(dets_up)
      ! generate connections = connected_dets_up,dn
      call find_connected_dets_heg(dets_up(i),dets_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn)
      ! loop over connections: binary search the dets list for each connection.
      do j=1,n_connected_dets
        call binary_search(connected_dets_up(j),connected_dets_dn(j),dets_up,dets_dn,k)
        if (k>0) then
          call hamiltonian_heg(dets_up(i),dets_dn(i),dets_up(k),dets_dn(k),matrix_element)
          v2(k) = v2(k) + matrix_element * v1(i)
        endif
      enddo
    enddo

  end subroutine apply_H_on_the_fly
  !===========================================================================

  !===========================================================================
  subroutine filter_dets(det_i_up, det_i_dn, valid_combination)
    !---------------------------------------------------------------------------
    ! Description : Determine if det_i is valid determinant (0 momentum)
    !               Note det_i_up, det_i_dn, are integers with their binary
    !               representation giving the occupancies of the spatial orbitals
    !
    ! Created     : F. Petruzielo, 4 Nov 2010
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

    !local
    real(rk), allocatable :: momentum(:)
    integer :: i,j

    allocate(momentum(n_dim))
    momentum = 0._rk
    valid_combination = .true.

    do j = 1, norb
       if (btest(det_i_up, j - 1)) then
          do i = 1, n_dim
             momentum(i) = momentum(i) + k_vectors(i, j)
          enddo
       endif
       if (btest(det_i_dn, j - 1)) then
          do i = 1, n_dim
             momentum(i) = momentum(i) + k_vectors(i, j)
          enddo
       endif
    enddo
    if (sum(momentum(:)**2) > epsilon ) then
       valid_combination = .false.
    endif

  end subroutine filter_dets
  !===========================================================================

  subroutine find_important_connected_dets_heg( &
      & det_up, det_dn, eps, n_connected_dets, connected_dets_up, connected_dets_dn, &
      & matrix_elements, ref_diag_elem, diag_elems_info, eps_big, matrix_elements_big)
  ! Same as find_connected_dets_chem, but instead of finding *all*
  ! connected dets, it finds all single excitations but only the
  ! important double excitations, i.e., only those with value
  ! greater than eps
  ! Requires efficient heatbath probabilities
  ! A Holmes, 22 Feb 2016

    use common_run, only : diag_elem_info, connected_diag_elems_info

    implicit none

    !dummy argument
    real(rk), intent(in) :: eps
    integer, intent(out) :: n_connected_dets
    real(rk), intent(out), optional :: matrix_elements(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_up
    type(ik_vec), intent(in) :: det_dn
    type(ik_vec), intent(out) :: connected_dets_up(:), connected_dets_dn(:)
    type(ik_vec) :: new_up, new_dn, tmp_det
#else
    integer(ik), intent(in) :: det_up
    integer(ik), intent(in) :: det_dn
    integer(ik), intent(out) :: connected_dets_up(:), connected_dets_dn(:)
    integer(ik) :: new_up, new_dn, tmp_det
#endif
    real(rk),optional,intent(in) :: ref_diag_elem
    type(diag_elem_info),optional,intent(out), target :: diag_elems_info(:)
    type(diag_elem_info), pointer :: diag_elem_info_ptr
    real(rk), optional, intent(in) :: eps_big
    real(rk), optional, intent(out) :: matrix_elements_big(:)

    integer :: p,q,r,s
    real(rk) :: matrix_element
    integer :: ipair,n_pairs
    integer :: i
    integer :: epair,hpair
    integer :: p2, q2, r_tmp
    logical :: valid_combination

    ! For memory optimized implementation in HEG.
    integer :: n_hpair ! Number of important connected dets.
    integer, allocatable :: hpairs_r(:), hpairs_s(:)
    integer :: diff_pq(3), diff_pr(3)
    integer :: diff_pq_ind(3) ! diff_pq + offset to get a value larger than 0 for index.

    ! First, add diagonal element, to be consistent with find_connected_dets_chem
    ! WARNING: Setting matrix_elements(1) to 0 since it is usually not needed for anything
    n_connected_dets = 1
    connected_dets_up(1) = det_up
    connected_dets_dn(1) = det_dn
    if (present(matrix_elements))  matrix_elements(1) = 0._rk
    if (present(matrix_elements_big))  matrix_elements_big(1) = 0._rk
    if (present(diag_elems_info))  diag_elems_info(1)%old_diag_elem = 1.e51_rk

    ! Get occ_up, occ_dn
    new_up = det_up
    new_dn = det_dn

    do p=1,nup
      occ_up(p) = trailz(new_up) + 1
      new_up = ibclr(new_up,occ_up(p)-1)
    enddo
    do p=1,ndn
      occ_dn(p) = trailz(new_dn) + 1
      new_dn = ibclr(new_dn,occ_dn(p)-1)
    enddo

    if (present(ref_diag_elem)) then
      do i=1,n_connected_dets
        diag_elems_info(i)%old_diag_elem = 1.e51_rk
      enddo
    endif

    ! Next, add double excitations
    if (eps > max_double)  return

    n_pairs = 0
    ! Same spin, up
    do p=1,nup
      do q=p+1,nup
        n_pairs = n_pairs + 1
        pairs_e1(n_pairs) = occ_up(p)
        pairs_e2(n_pairs) = occ_up(q)
      enddo
    enddo
    ! Same spin, dn
    do p=1,ndn
      do q=p+1,ndn
        n_pairs = n_pairs + 1
        pairs_e1(n_pairs) = occ_dn(p)+norb
        pairs_e2(n_pairs) = occ_dn(q)+norb
      enddo
    enddo

    ! Opposite spin
    do p=1,nup
      do q=1,ndn
        n_pairs = n_pairs + 1
        pairs_e1(n_pairs) = occ_up(p)
        pairs_e2(n_pairs) = occ_dn(q)+norb
      enddo
    enddo

    allocate(hpairs_r(norb))
    allocate(hpairs_s(norb))

    do ipair=1,n_pairs
      p = pairs_e1(ipair)
      q = pairs_e2(ipair)

      p2 = p
      q2 = q
      if (p > norb .and. q > norb) then
        p2 = p - norb
        q2 = q - norb
      elseif (p <= norb .and. q > norb .and. p > q - norb) then
        p2 = q - norb
        q2 = p + norb
      endif

      epair = int(combine_2_indices(p2, q2))

      ! Obtain orbital pairs to excite to.
      n_hpair = 0 ! Number of possible excitations.
      do hpair = 1, n_diff**3
        if (p2 <= norb .and. q2 <= norb) then
          ! Same spin.
          diff_pq = k_vectors_rel(q2, :) - k_vectors_rel(p2, :)
          diff_pq_ind = diff_pq + n_diff_offset
          if (dtm_hb_same_spin( &
              & diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3), hpair)%absH <= eps) then
            exit
          endif
          diff_pr = &
              & dtm_hb_same_spin(diff_pq_ind(1), diff_pq_ind(2), diff_pq_ind(3), hpair)%int3
          r = find_orb_id(diff_pr + k_vectors_rel(p2, :))
          if (r < 0) then
            cycle
          endif
          s = &
              & find_orb_id(k_vectors_rel(p2, :) + k_vectors_rel(q2, :) &
              & - k_vectors_rel(r, :))
          if (s < 0 .or. s < r) then
            cycle
          endif
        else
          if (dtm_hb_opposite_spin(hpair)%absH <= eps) then
            exit
          endif

          diff_pr = dtm_hb_opposite_spin(hpair)%int3
          r = find_orb_id(diff_pr + k_vectors_rel(p2, :))

          if (r < 0) then
            cycle
          endif

          s = &
              & find_orb_id(k_vectors_rel(p2, :) + k_vectors_rel(q2 - norb, :) &
              & - k_vectors_rel(r, :))
          if (s < 0) then
            cycle
          endif
          s = s + norb
        endif
        n_hpair = n_hpair + 1
        hpairs_r(n_hpair) = r
        hpairs_s(n_hpair) = s
      enddo

      ! Process each possible excitation.
      do hpair = 1, n_hpair

        r = hpairs_r(hpair)
        s = hpairs_s(hpair)

        if (p > norb .and. q > norb) then
          r = r + norb
          s = s + norb
        elseif (p <= norb .and. q > norb .and. p > q - norb) then
          r_tmp = s - norb
          s = r + norb
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

        n_connected_dets = n_connected_dets + 1
        connected_dets_up(n_connected_dets)=new_up
        connected_dets_dn(n_connected_dets)=new_dn

        if (present(matrix_elements) .or. present(matrix_elements_big)) then
          call hamiltonian_heg(det_up,det_dn,new_up,new_dn,matrix_element)
        endif
        if (present(matrix_elements)) then
          matrix_elements(n_connected_dets) = matrix_element
        endif
        if (present(matrix_elements_big)) then
          if (.not. present(eps_big)) stop 'eps_big is unavailable.'
          if (abs(matrix_element) > eps_big) then
            matrix_elements_big(n_connected_dets) = matrix_element
          else
            matrix_elements_big(n_connected_dets) = 0._rk
          endif
        endif

        ! Save info for computing new diagonal element in O(N) time
        if (present(ref_diag_elem)) then
          diag_elem_info_ptr => diag_elems_info(n_connected_dets)
          diag_elem_info_ptr%old_diag_elem = ref_diag_elem
          diag_elem_info_ptr%p = p
          diag_elem_info_ptr%q = q
          diag_elem_info_ptr%r = r
          diag_elem_info_ptr%s = s
        endif

      enddo ! hpair (r,s)

      ! if (p2 == 2 .and. q2 == 4) stop 'debug'
    enddo ! ipair (p,q)

    ! stop 'debug'

  end subroutine find_important_connected_dets_heg
  !===========================================================================

  subroutine find_connected_dets_heg(det_up, det_dn, n_connected_dets, connected_dets_up, connected_dets_dn, matrix_elements)
    implicit none

    integer, intent(inout) :: n_connected_dets
    !integer, intent(in),optional :: n_orb_include
    real(rk), intent(out), optional :: matrix_elements(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_up
    type(ik_vec), intent(in) :: det_dn
    type(ik_vec), intent(out) :: connected_dets_up(:), connected_dets_dn(:)
    type(ik_vec) :: temp_det_up, temp_det_dn
#else
    integer(ik), intent(in) :: det_up
    integer(ik), intent(in) :: det_dn
    integer(ik), intent(out) :: connected_dets_up(:), connected_dets_dn(:)
    integer(ik) :: temp_det_up, temp_det_dn
#endif

    call find_important_connected_dets_heg( &
        & det_up, det_dn, 1.e-15_rk, n_connected_dets, connected_dets_up, connected_dets_dn, matrix_elements)
  end subroutine find_connected_dets_heg
  !===========================================================================

  function is_connected_heg(det_i_up, det_i_dn, det_j_up, det_j_dn)
    !---------------------------------------------------------------------------
    ! Description : Determine if det_i is connected to det_j by the hamiltonian.
    !               In other words, is det_i=det_j or are they related by double excitation
    !
    ! Created     : F. Petruzielo, 8 Nov 2010
    !---------------------------------------------------------------------------

    !dummy arguments
    logical :: is_connected_heg
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_i_up
    type(ik_vec), intent(in) :: det_i_dn
    type(ik_vec), intent(in) :: det_j_up
    type(ik_vec), intent(in) :: det_j_dn
    type(ik_vec) :: unique_i_up, unique_j_up, unique_i_dn, unique_j_dn
#else
    integer(ik), intent(in) :: det_i_up
    integer(ik), intent(in) :: det_i_dn
    integer(ik), intent(in) :: det_j_up
    integer(ik), intent(in) :: det_j_dn
    integer(ik) :: unique_i_up, unique_j_up, unique_i_dn, unique_j_dn
#endif

    !local variables
    integer :: k
    integer :: temp_tot_uniq

    !are determinants equal
    if (det_i_up == det_j_up .and. det_i_dn == det_j_dn) then
       is_connected_heg = .true.
       return
    endif

    !are determinants related by double excitation
    unique_i_up = iand(det_i_up, not(det_j_up))
    unique_j_up = iand(det_j_up, not(det_i_up))
    unique_i_dn = iand(det_i_dn, not(det_j_dn))
    unique_j_dn = iand(det_j_dn, not(det_i_dn))

    temp_tot_uniq = 0
    do k = 1, norb
       if (btest(unique_i_up, k - 1)) then
          temp_tot_uniq = temp_tot_uniq + 1
       endif
       if (btest(unique_i_dn, k - 1)) then
          temp_tot_uniq = temp_tot_uniq + 1
       endif
    enddo
    if (temp_tot_uniq > 2) then
       is_connected_heg = .false.
       return
    endif

    temp_tot_uniq = 0
    do k = 1, norb
       if (btest(unique_j_up, k - 1)) then
          temp_tot_uniq = temp_tot_uniq + 1
       endif
       if (btest(unique_j_dn, k - 1)) then
          temp_tot_uniq = temp_tot_uniq + 1
       endif
    enddo
    if (temp_tot_uniq > 2) then
       is_connected_heg = .false.
       return
    endif

    !yes
    is_connected_heg = .true.

  end function is_connected_heg
  !===========================================================================

  !===========================================================================
  subroutine madelung_energy()
    !---------------------------------------------------------------------------
    ! Description : Calculate the Madelung energy. For the formula in 3d see
    !               Dr. Umrigar's ewald notes.
    !
    ! Created     : F. Petruzielo, 15 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !local variables
    integer :: i, j, k
    integer :: n_max, index
    real(rk), parameter :: epsilon = 1.e-10_rk
    real(rk), allocatable :: g_vectors(:, :), values(:)
    real(rk) :: kappa, g_max

    if (n_dim /= 3) then
       stop 'Madelung energy_madelung is only implemented for 3d'
    endif

    !kappa is chosen large enough that we can leave out erfc term
    kappa = 10._rk / length_cell

    !determine max n for the g-vectors to include in the calculation
    n_max = 1
    do
       g_max = (2 * pi * n_max / length_cell)
       if (4 * pi / length_cell**3 * exp(-(g_max/(2*kappa))**2) / g_max**2 < epsilon) exit !neglible contribution so exit
       n_max = n_max + 1
    enddo

    !we define all g-vectors within a cube of side 2 * 2 pi n_max / L
    !it is unneccessary to restrict to sphere because calculating a convergent sum
    allocate(g_vectors(n_dim, (2*n_max + 1)**n_dim))
    allocate(values(2*n_max + 1))

    !define possible g-point values for a single coordinate
    index = 1
    do i = -n_max, n_max
       values(index) = 2 * pi / length_cell * i
       index = index + 1
    enddo

    !define g-vector n-tuples where n=n_dim
    if (n_dim == 3) then
       index = 1
       do i = -n_max, n_max
          do j = -n_max, n_max
             do k = -n_max, n_max
                g_vectors(1, index) = values(i + n_max+1)
                g_vectors(2, index) = values(j + n_max+1)
                g_vectors(3, index) = values(k + n_max+1)
                index = index + 1
             enddo
          enddo
       enddo
    ! elseif (n_dim == 2) then
    !    index = 1
    !    do i = -n_max, n_max
    !       do j = -n_max, n_max
    !          g_vectors(1, index) = values(i + n_max+1)
    !          g_vectors(2, index) = values(j + n_max+1)
    !          index = index + 1
    !       enddo
    !    enddo
    endif

    energy_madelung = 0._rk
    do i = 1, size(g_vectors, 2)
       if (sum(g_vectors(:, i)**2) < epsilon) cycle !dont include g=0
       energy_madelung = energy_madelung + exp(- sum(g_vectors(:, i)**2)/(2*kappa)**2) / sum(g_vectors(:, i)**2)
    enddo
    energy_madelung = energy_madelung * 4 * pi / length_cell**3
    energy_madelung = energy_madelung - pi / length_cell**3 / kappa**2 - 2 * kappa / pi ** 0.5_rk
    energy_madelung = energy_madelung * nelec / 2._rk

    write(6,'(''Madelung energy ='', f10.6)') energy_madelung

  end subroutine madelung_energy
  !===========================================================================

  !===========================================================================
  subroutine theoretical_hf_heg()
    !---------------------------------------------------------------------------
    ! Description : Calculate theoretical HF results
    !
    ! Created     : F. Petruzielo, 1 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !local variables
    real(rk) :: kinetic_energy, exchange_energy, alpha, polarization

    write (6,'(/,''Theoretical infinite cell Hartree-Fock Results'')')
    write (6,'(a,i1,a,f6.4)') "n_dim = ", n_dim, ", r_s = ", r_s

    if (n_dim == 2) then
       alpha = 1._rk / sqrt(2._rk)
       kinetic_energy = 0.5 * 0.5_rk / (alpha * r_s)**2
       exchange_energy = 0.5 * (-8._rk*sqrt(2._rk) / (3*pi*r_s))
    elseif (n_dim == 3) then
       write (6,'(''r_s = (3/(4*pi*n))^(1/3);  E_F = hbar^2*k_F^2/(2m)'')')
       write (6,'(''E_tot = E_kin + E_x = (3/5)E_F -(3/4)(3n/pi)^(1/3) = 1.105/r_s^2-0.458/r_s'')')
       alpha = (4._rk / (9._rk * pi))**(1._rk / 3._rk)
       kinetic_energy = 0.5 * 0.6_rk / (alpha * r_s)**2
       exchange_energy = 0.5 * (-3._rk / (2*pi*r_s*alpha))
    endif

    polarization = (nup - ndn) / (1._rk * nelec)
    kinetic_energy = kinetic_energy * ( (1 + polarization)**((n_dim + 2._rk) / n_dim) + (1 - polarization)**((n_dim + 2._rk) / n_dim) ) / 2._rk
    exchange_energy = exchange_energy * ( (1 + polarization)**((n_dim + 1._rk) / n_dim) + (1 - polarization)**((n_dim + 1._rk) / n_dim) ) / 2._rk

    write (6,'(a,t32,f10.6)') "Kinetic energy per particle = ", kinetic_energy
    write (6,'(a,t32,f10.6)') "Exchange energy per particle = ", exchange_energy
    write (6,'(a,t32,f10.6)') "Total HF energy per particle = ", kinetic_energy + exchange_energy
    write (6,'(a,i3)') "Scaled by number of electrons = ", nelec
    write (6,'(a,t32,f10.6)') "Kinetic energy = ", kinetic_energy * nelec
    write (6,'(a,t32,f10.6)') "Exchange energy = ", exchange_energy * nelec
    write (6,'(a,t32,f10.6,/)') "Total HF energy = ", (kinetic_energy + exchange_energy) * nelec

  end subroutine theoretical_hf_heg
  !===========================================================================

  !===========================================================================
  subroutine test_connected_dets()
    !---------------------------------------------------------------------------
    ! Description : Test find_connected_dets_heg. det_up = det_dn = 31
    !               appropriate for input with 5 up and 5 dn electrons (this can be changed below)
    !
    ! Created     : F. Petruzielo, 6 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !local variables
    !integer :: i
    integer :: i, n_connected_dets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec) :: det_up, det_dn
    type(ik_vec), allocatable :: connected_dets_up(:), connected_dets_dn(:)
#else
    integer(ik) :: det_up, det_dn
    integer(ik), allocatable :: connected_dets_up(:), connected_dets_dn(:)
#endif
    real(rk), allocatable :: matrix_elements(:)

    det_up = 31
    det_dn = 31

    !call find_connected_dets_heg(det_up, det_dn, connected_dets_up, connected_dets_dn, matrix_elements)
    call find_connected_dets_heg(det_up, det_dn, n_connected_dets, connected_dets_up, connected_dets_dn, matrix_elements)

    write(6,'(/,''Test find_connected_dets_heg 31 31: found'',i6)') n_connected_dets
    !do i = 1, size(connected_dets_up)
    do i = 1, n_connected_dets
       write(6,'(4i6,es12.4)') det_up, det_dn, connected_dets_up(i), connected_dets_dn(i), matrix_elements(i)
    enddo

  end subroutine test_connected_dets
  !===========================================================================

  !===========================================================================
  subroutine test_off_diagonal_move_heg()
    !---------------------------------------------------------------------------
    ! Description : Test off-diagonal move. det_up = det_dn = 31
    !               appropriate for input with 5 up and 5 dn electrons (this can be changed below)
    !
    ! Created     : F. Petruzielo, 8 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !local variables
!   integer :: i, weight_i, weight_j
    integer :: i, weight_i
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec) :: det_i_up, det_i_dn, det_j_up, det_j_dn
#else
    integer(ik) :: det_i_up, det_i_dn, det_j_up, det_j_dn
#endif
    real(rk) :: diagonal_matrix_element, proj_energy_j, weight_j_times_trial_wf_weight_j, weight_j

    det_i_up = 31
    det_i_dn = 31
    weight_i = 1

    write(6,'(/,''Test off_diagonal_move_heg for 31 31'')')
    do i = 1, 100
!      call off_diagonal_move_heg(det_i_up, det_i_dn, weight_i, det_j_up, det_j_dn, weight_j, diagonal_matrix_element)
       call off_diagonal_move_heg(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j)
       if (abs(weight_j) > 0) then
          call energy_pieces_heg(det_j_up, det_j_dn, proj_energy_j, weight_j_times_trial_wf_weight_j)
          write(*,*) det_j_up, det_j_dn, weight_j, diagonal_matrix_element, proj_energy_j, weight_j_times_trial_wf_weight_j
       endif
    enddo

    det_i_up = 279
    det_i_dn = 91
    weight_i = 1

    write(6,'(/,''Test off_diagonal_move_heg for 279 91'')')
    do i = 1, 100
!      call off_diagonal_move_heg(det_i_up, det_i_dn, weight_i, det_j_up, det_j_dn, weight_j, diagonal_matrix_element)
       call off_diagonal_move_heg(det_i_up, det_i_dn, det_j_up, det_j_dn, weight_j)
       if (abs(weight_j) > 0) then
          call energy_pieces_heg(det_j_up, det_j_dn, proj_energy_j, weight_j_times_trial_wf_weight_j)
          write(*,*) det_j_up, det_j_dn, weight_j, diagonal_matrix_element, proj_energy_j, weight_j_times_trial_wf_weight_j
       endif
    enddo

  end subroutine test_off_diagonal_move_heg
  !===========================================================================

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
  !===========================================================================

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
  !===========================================================================

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
  !===========================================================================

  real(rk) function double_excitation_matrix_element_no_ref_abs(p,q,r,s)
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

    zero_det = 0_ik

    if (p==q.or.r==s.or.p==r.or.q==s.or.p==s.or.q==r) then
      double_excitation_matrix_element_no_ref_abs = 0._rk
      return
    endif
    if (p<=norb.and.q<=norb) then ! same spin, up
      call hamiltonian_heg(ibset(ibset(zero_det,p-1),q-1),zero_det,ibset(ibset(zero_det,r-1),s-1),zero_det,elem, .true.)
    elseif (p>norb.and.q>norb) then ! same spin, dn
      call hamiltonian_heg(zero_det,ibset(ibset(zero_det,p-norb-1),q-norb-1),zero_det,ibset(ibset(zero_det,r-norb-1),s-norb-1),elem, .true.)
    else ! opposite spin
      call hamiltonian_heg(ibset(zero_det,p-1),ibset(zero_det,q-norb-1),ibset(zero_det,r-1),ibset(zero_det,s-norb-1),elem, .true.)
    endif

    double_excitation_matrix_element_no_ref_abs = elem

  end function double_excitation_matrix_element_no_ref_abs
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
  !===========================================================================

  subroutine sort_int3_absH(arr)
    ! Sort arr of int3_absH by absH in decreasing order

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! From http://jblevins.org/mirror/amiller/qsort.f90
    ! Modified to be of type rs_absH by A Holmes, 24 Jan 2016

    implicit none
    type(int3_absH),intent(inout) :: arr(:)

    call quick_sort(1, size(arr))

    contains

      recursive subroutine quick_sort(left_end, right_end)

        integer, intent(in) :: left_end, right_end

        !     local variables
        integer             :: i, j
        type(int3_absH)          :: reference, temp
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
        type(int3_absH)          :: temp

        do i = left_end, right_end - 1
          do j = i+1, right_end
            if (arr(i)%absH < arr(j)%absH) then
              temp = arr(i); arr(i) = arr(j); arr(j) = temp
            endif
          enddo
        enddo
      end subroutine interchange_sort

  end subroutine sort_int3_absH
  !===========================================================================
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
          call hamiltonian_heg(ibset(zero_det,p-1),zero_det,ibset(zero_det,q-1),zero_det,a)
#else
          call hamiltonian_heg(ibset(0_ik,p-1),0_ik,ibset(0_ik,q-1),0_ik,a)
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
  !===========================================================================

  subroutine get_new_diag_elem_heg(info,new_det_up,new_det_dn,new_diag_elem)
  ! This is not used currently.
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
    real(rk) :: tmp
    integer,allocatable :: occ_up_tmp(:),occ_dn_tmp(:)

    allocate(occ_up_tmp(nup))
    allocate(occ_dn_tmp(ndn))

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

    call get_occ_orbs(new_det_up,new_det_dn,occ_up_tmp,occ_dn_tmp)

    ! O(N) Two-body direct part: E += sum_{i in occ. but not in (r,s)} direct(i,r) + direct(i,s) - direct(i,p) - direct(i,q)
    do j=1,nup
      i = occ_up_tmp(j)
      if (i==r_in.or.i==s_in)  cycle ! only need to check r and s because 'occupied' orbitals are those occupied in new det, which includes r and s but not p and q
      new_diag_elem = new_diag_elem + integral_value(i,i,r,r) + integral_value(i,i,s,s) - integral_value(i,i,p,p) - integral_value(i,i,q,q)
    enddo
    do j=1,ndn
      i = occ_dn_tmp(j)
      if (i==r_in-norb.or.i==s_in-norb)  cycle
      new_diag_elem = new_diag_elem + integral_value(i,i,r,r) + integral_value(i,i,s,s) - integral_value(i,i,p,p) - integral_value(i,i,q,q)
    enddo

    ! O(N) Two-body exchange part: E += sum_{i in occ. but not in (r,s)} -exchange(i,r) - exchange(i,s) + exchange(i,p) + exchange(i,q)
    if (p_up.or.q_up) then
      do j=1,nup
        i = occ_up_tmp(j)
        if (i==r_in.or.i==s_in)  cycle
        if (p_up)  new_diag_elem = new_diag_elem - integral_value(i,r,r,i) + integral_value(i,p,p,i)
        if (q_up)  new_diag_elem = new_diag_elem - integral_value(i,s,s,i) + integral_value(i,q,q,i)
      enddo
    endif
    if ((.not.p_up).or.(.not.q_up)) then
      do j=1,ndn
        i = occ_dn_tmp(j)
        if (i==r_in-norb.or.i==s_in-norb)  cycle
        if (.not.p_up)  new_diag_elem = new_diag_elem - integral_value(i,r,r,i) + integral_value(i,p,p,i)
        if (.not.q_up)  new_diag_elem = new_diag_elem - integral_value(i,s,s,i) + integral_value(i,q,q,i)
      enddo
    endif

  end subroutine get_new_diag_elem_heg
  !===========================================================================

  function integral_value(p, q, r, s, tmp)
    implicit none

    integer, intent(in) :: p, q, r, s
    integer, optional, intent(in) :: tmp
    real(rk) :: integral_value

    real(rk) :: k_pr(3), k_qs(3), k_pq_diff(3)

    if (p == q) then
      integral_value = 0._rk
      return
    endif

    if (r == norb + 1 .and. s == norb + 1) then
      integral_value = sum(k_vectors(:, p)**2) * 0.5_rk
      return
    endif

    k_pr = k_vectors(:, p) + k_vectors(:, r)
    k_qs = k_vectors(:, q) + k_vectors(:, s)
    if (sum((k_pr - k_qs)**2) * length_cell**2 > 1.0e-3_rk) then
      integral_value = 0._rk
      return
    endif

    k_pq_diff = k_vectors(:, q) - k_vectors(:, p)
    integral_value = 4._rk * pi / (sum(k_pq_diff**2) * length_cell**3)
  end function integral_value
  !===========================================================================

  subroutine diagonalize_sparse_hamiltonian_heg_excited(dets_up, dets_dn, n_det, n_states, lowest_eigenvectors, lowest_eigenvalues, time_sym_on, initial_vectors, sparse_ham)

    use types, only : i8b, num_words
    use more_tools, only : real_symmetric_diagonalize_ow_ham,davidson_sparse
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
    integer,intent(in) :: n_det,n_states
    real(rk), intent(out) :: lowest_eigenvectors(:,:)
    real(rk), intent(out) :: lowest_eigenvalues(:)
    logical,intent(in)   :: time_sym_on
    real(rk),intent(in) :: initial_vectors(:,:)
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
    real(rk) :: tau_optimal_deterministic, tau_optimal_stochastic
    logical  :: sym_on
    real(rk),allocatable :: ham(:,:),eigenvalues(:)
    real(rk) :: matrix_element
    integer :: i,j,excite_level
    logical               :: srt
    integer :: n_connected_dets
    character(len=2) :: fmt

    call my_second(1,'diagonalize_sparse_hamiltonian_heg')

    sym_on=.false.

    call generate_sparse_ham_heg_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,H_values,sym_on,hf_to_psit=.false.,sparse_ham=sparse_ham) ! Pass in old H, return new H
    write (6,'(/,''Performing Davidson using sparsely stored H and matrix_lanczos in more_tools, n_det='',i10,es10.2)') n_det, real(n_det)
    call flush(6)
    call davidson_sparse(n_det,n_states,lowest_eigenvectors,lowest_eigenvalues,H_indices,H_nonzero_elements,H_values,initial_vectors)
   !call matrix_lanczos(n_det,lowest_eigenvector,lowest_eigenvalue,H_indices,H_nonzero_elements,H_values,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
    deallocate(H_values)
    deallocate(H_indices)
    deallocate(H_nonzero_elements)

    write(6,'(/,''Wavefunction has'',i8,'' determinants'')') n_det
    write(6,'(/,''Lowest eigenvalues='',10f13.6)') lowest_eigenvalues(1:n_states)

! Calculate number of dets with various excitation levels and the sum of their squared wts.
    call print_excitation_levels_and_wts(n_det,dets_up,dets_dn,lowest_eigenvectors(:,1),norb)

    call my_second(2,'diagonalize_sparse_hamiltonian_heg')

  end subroutine diagonalize_sparse_hamiltonian_heg_excited
  !=====================================================================================================

  subroutine generate_sparse_ham_heg_upper_triangular(dets_up,dets_dn,H_indices,H_nonzero_elements,H_values,time_sym_on,hf_to_psit,max_nonzero_elems,num_nonzero_elems,sparse_ham)
    !---------------------------------------------------------------------------
    ! Description : Compute the Hamiltonian in a given space and store it sparsely, as an upper triangular matrix.
    !
    ! Created     : A Holmes, 17 Dec 2012
    ! Modified    : A Holmes, 11 Jan 2013. Input hf_to_psit replaces the first state (HF) with psi_trial.
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
#else
      integer(ik),intent(in) :: dets_up(:),dets_dn(:)
      integer(ik),allocatable :: beta(:),alpha_m1(:)
#endif
      integer(i8b),allocatable,intent(out) :: H_indices(:),H_nonzero_elements(:)
      real(rk),allocatable,intent(out) :: H_values(:)
      logical,intent(in)               :: time_sym_on
      integer(i8b),optional,intent(in) :: max_nonzero_elems
      integer(i8b),optional,intent(out) :: num_nonzero_elems
      logical,intent(in) :: hf_to_psit
      type(sparse_mat),optional,intent(inout) :: sparse_ham

      ! Local
      integer(i8b),allocatable :: tmp_H_indices(:)
      real(rk),allocatable :: tmp_H_values(:)
      integer(i8b) :: i
      integer :: k
      integer(i8b) :: j
      integer :: excite_level
      integer(i8b) :: n_det
      integer(i8b) :: nonzero_elements,old_nonzero_elements,isparse
      real(rk) :: matrix_element
      logical             :: sym_on
      integer :: n_connected_dets
      integer :: connected_indices(max_connected_dets)
      integer :: numops
      logical :: brute_force ! whether or not to use the brute force method (instead of binary search)
      logical :: use_partial_connections ! whether or not to use Sandeep's partial connections idea
      integer,allocatable :: beta_ind(:),alpha_m1_ind(:)
      integer :: isparse_old
      integer(i8b) :: n_copy
      integer(i8b) :: ndet_old,ndet_new

      call my_second(1, 'generate_sparse_ham_heg_upper_triangular')

      sym_on=.false.

      n_det=size(dets_up)

      allocate(H_nonzero_elements(n_det))
      H_nonzero_elements(:)=0

      ! Count number of nonzero elements of H
      nonzero_elements=0

      call find_connected_dets_heg(dets_up(1), dets_dn(1), n_connected_dets, connected_dets_up, connected_dets_dn)

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

      nonzero_elements = int(n_det*real(min(n_det,n_connected_dets))**0.5,i8b) ! Just a starting guess

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

          ! If old sparse ham present, then only need to create helper lists for the new dets

          ndet_old = sparse_ham%ndet
          ndet_new = n_det - ndet_old

          ! First, sort a list of (dets_dn(i),i) by dets_dn
          allocate(beta(ndet_new))
          allocate(beta_ind(ndet_new))
          beta(:) = dets_dn(ndet_old+1:n_det)
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

        else

          ! First, sort a list of (dets_dn(i),i) by dets_dn
          allocate(beta(n_det))
          allocate(beta_ind(n_det))
          beta(:) = dets_dn(:)
          do i=1,n_det
            beta_ind(i) = i
          enddo
          call sort_by_first_argument(int(n_det),beta,beta_ind)

          ! Next, generate alpha_m1, a list of all up configurations reachable by removing one up electron
          allocate(alpha_m1(n_det*nup))
          allocate(alpha_m1_ind(n_det*nup))
          call get_n_minus_1_configs(n_det,nup,dets_up,alpha_m1)
          do i=1,n_det
            alpha_m1_ind(nup*(i-1)+1:nup*i) = i
          enddo
          call sort_by_first_argument(int(n_det*nup),alpha_m1,alpha_m1_ind)

        endif

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
                write (6,'(''Reallocating from size'',i14,'' ='',es11.4,'' to size'',i14,'' ='',es11.4)') nonzero_elements, real(nonzero_elements), nonzero_elements+n_copy, real(nonzero_elements+n_copy)
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
          call get_connected_dets_in_list(int(i),dets_up(i),dets_dn(i),dets_up(1:n_det),dets_dn(1:n_det),beta,beta_ind,alpha_m1,alpha_m1_ind,n_connected_dets,connected_indices,connected_matrix_elements)

          ! Add to Hamiltonian if connected_index >= i
          do j=1,n_connected_dets
            if (connected_indices(j)>i.or.(connected_indices(j)==i.and.i>ndet_old)) then
              ! Add to H
              isparse=isparse+1_i8b
              if (isparse>nonzero_elements) then ! increment nonzero_elements by 41%
                write (6,'(''Reallocating from size'','' ='',es11.4,i14,'' to size'',i14,'' ='',es11.4)') real(nonzero_elements), (nonzero_elements), int(1.414*nonzero_elements,i8b), 1.414*nonzero_elements
                allocate(tmp_H_indices(nonzero_elements))
                allocate(tmp_H_values(nonzero_elements))
                tmp_H_indices(1_i8b:nonzero_elements) = H_indices(1_i8b:nonzero_elements)
                tmp_H_values(1_i8b:nonzero_elements) = H_values(1_i8b:nonzero_elements)
                deallocate(H_indices,H_values)
                old_nonzero_elements = nonzero_elements
                nonzero_elements = max(3,int(1.414*nonzero_elements,i8b))
                allocate(H_indices(nonzero_elements))
                allocate(H_values(nonzero_elements))
                H_indices(1_i8b:old_nonzero_elements) = tmp_H_indices(1_i8b:old_nonzero_elements)
                H_values(1_i8b:old_nonzero_elements) = tmp_H_values(1_i8b:old_nonzero_elements)
                deallocate(tmp_H_indices,tmp_H_values)
              endif
              H_indices(isparse)=connected_indices(j)
              H_values(isparse)=connected_matrix_elements(j)
              H_nonzero_elements(i)=H_nonzero_elements(i)+1
            endif
          enddo ! connected dets

        enddo ! reference dets

      else
        stop 'Only partial connection implemented.'
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

      call my_second(2, 'generate_sparse_heg_heg_upper_triangular')

  end subroutine generate_sparse_ham_heg_upper_triangular
  !=====================================================================================================

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

    integer(i8b) :: i_det
    integer :: i_el
    integer,allocatable :: occ(:)

    allocate(occ(nelec))

    do i_det = 1,n_det
      call get_occ_orbs(dets(i_det),occ)
      do i_el = 1,nelec
        nm1_dets(nelec*(i_det-1)+i_el) = ibclr(dets(i_det),occ(i_el)-1)
      enddo ! i_el
    enddo ! i_det

  end subroutine get_n_minus_1_configs

  !=====================================================================================================

  subroutine get_connected_dets_in_list(index,det_up,det_dn,dets_up,dets_dn,beta,beta_ind,alpha_m1,alpha_m1_ind,n_connected_dets,connected_indices,connected_matrix_elements)
  ! Use Sandeep's partial connections to get the matrix elements in the list represented by beta and alpha_m1
  ! (in HCI that list is the list of variational dets)
  ! A Holmes, 10 Nov 2016

    use tools, only : sort_and_merge
    use more_tools, only : get_occ_orbs,binary_search_single

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
    real(rk),intent(out) :: connected_matrix_elements(:)
    real(rk) :: elem

    integer :: i_el,i,j,k
    integer :: n_beta, n_alpha_m1

    n_beta = size(beta)
    n_alpha_m1 = size(alpha_m1)

    ! Diagonal element first
    n_connected_dets = 1
    call hamiltonian(det_up,det_dn,det_up,det_dn,elem)
    connected_indices(1) = index
    connected_matrix_elements(1) = elem

    ! Compute all alpha excitations: iterate over all pairs of dets in each dets_dn(i)
    call binary_search_single(det_dn,beta,k)
    if (k.ne.0) then ! It could only be 0 if part of the sparse_ham had been computed on a previous iteration
      i=k
      do while (det_dn==beta(i))
        j = beta_ind(i)
        if (j.ne.index) then
          call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
          if (abs(elem).gt.1.d-12) then
            n_connected_dets = n_connected_dets + 1
            connected_indices(n_connected_dets) = j
            connected_matrix_elements(n_connected_dets) = elem
          endif
        endif
        i = i-1
        if (i==0)  exit
      enddo
      i=k
      do while (det_dn==beta(i))
        j = beta_ind(i)
        if (j.ne.index) then
          call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
          if (abs(elem).gt.1.d-12) then
            n_connected_dets = n_connected_dets + 1
            connected_indices(n_connected_dets) = j
            connected_matrix_elements(n_connected_dets) = elem
          endif
        endif
        i = i+1
        if (i==n_beta+1)  exit
      enddo
    endif ! k.ne.0

    ! Compute all other excitations: iterate over all pairs of dets in each alpha_m1
    call get_occ_orbs(det_up,occ_up)
    do i_el=1,nup
      tmp = ibclr(det_up,occ_up(i_el)-1)
      call binary_search_single(tmp,alpha_m1,k)
      if (k.ne.0) then ! It could only be 0 if part of the sparse_ham had been computed on a previous iteration
        i=k
        do while (tmp==alpha_m1(i))
          j = alpha_m1_ind(i)
          if (j.ne.index) then
            call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
            if (abs(elem).gt.1.d-12) then
              n_connected_dets = n_connected_dets + 1
              connected_indices(n_connected_dets) = j
              connected_matrix_elements(n_connected_dets) = elem
            endif
          endif
          i = i-1
          if (i==0)  exit
        enddo
        i=k
        do while (tmp==alpha_m1(i))
          j = alpha_m1_ind(i)
          if (j.ne.index) then
            call hamiltonian(det_up,det_dn,dets_up(j),dets_dn(j),elem)
            if (abs(elem).gt.1.d-12) then
              n_connected_dets = n_connected_dets + 1
              connected_indices(n_connected_dets) = j
              connected_matrix_elements(n_connected_dets) = elem
            endif
          endif
          i = i+1
          if (i==n_alpha_m1+1)  exit
        enddo
      endif ! k.ne.0
    enddo

    ! Now, remove duplicates by sorting and merging
    call sort_and_merge(n_connected_dets,connected_indices,connected_matrix_elements)

    ! TODO: Only compute matrix element after merge, to avoid computing expensive single excitations and diagonal elements multiple times
   !do i=1,n_connected_dets
   !  call hamiltonian(det_up,det_dn,dets_up(connected_indices(i)),dets_dn(connected_indices(i)),elem)
   !  connected_matrix_elements(i) = elem
   !enddo
    !       Also, compute matrix elements in O(N_el) time rather than O(N_el^2) time

  end subroutine get_connected_dets_in_list

  !=====================================================================================================

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
    call excitation_level(up1, dn1, up2, dn2, excite_level)
    if (excite_level >= 0) call hamiltonian_heg(up1, dn1, up2, dn2, elem)

  end subroutine hamiltonian

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

end module heg
