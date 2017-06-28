module common_selected_ci
  use types, only : rk,ik,ik_vec,i8b, optional_rk, optional_integer

  integer :: det_sel_iters ! number of times find connected dets and truncation are performed
  integer,allocatable :: norb_det_sel(:),n_sym_uniq_det_det_sel(:)
  integer :: ndet_det_sel
  real(rk), allocatable :: cdet_det_sel(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable :: dets_up_det_sel(:), dets_dn_det_sel(:)
#else
  integer(ik), allocatable :: dets_up_det_sel(:), dets_dn_det_sel(:)
#endif
  integer :: lanczos_iters,lanczos_initiators,lanczos_truncate
  integer(i8b),parameter :: max_nonzero_elements=int(5e9,i8b) ! Maximum number of nonzero matrix elements that can be stored in a sparse matrix for Lanczos
  logical :: too_big_to_store
  real(rk) :: log_num_nonzero_elements
  ! These are for Norm's IACI algorithm which is no longer used:
  integer :: cdets,tdets
  ! These are for our HCI algorithm:
  real(rk) :: eps_var_sched(30)=0_rk ! Schedule for how eps_var changes at each HCI iteration. If omitted, it is eps_var. (typically start at 2*eps_var and go down to eps_var)
  real(rk) :: eps_var           ! Cutoff for terms in the sum for variational wf. (typically 1e-3 - 2e-5)
  real(rk) :: eps_pt            ! Cutoff for terms in the sum for 2PT energy correction (typically 1e-6 - 1e-8)
  real(rk) :: eps_pt_big = -1   ! Cutoff for terms in the deterministic part of the semistochastic 2PT energy correction (typically 1e-5 - 1e-7)
  real(rk) :: target_error      ! Maximum standard deviation in 2PT energy correction when computed stochastically (typically 5e-4)
  logical :: dump_wf_var
  type sparse_mat
    integer(i8b) :: ndet = 0
    integer(i8b),allocatable :: indices(:),nonzero_elements(:)
    real(rk),allocatable :: values(:)
  end type sparse_mat
  type(sparse_mat), save :: sparse_ham
  integer :: n_states ! number of ground and excited states to compute with HCI
  real(rk) :: eps_pt_big_energy = 1 ! Energy from a previous eps_pt_big run
  real(rk) :: n_max_connections = -1 ! It is real, to allow input in scientific notation. If <=0, it is set depending on available memory
  integer :: n_mc = -1 ! Number of monte carlo samples for stochastic PT. If <=0, it is set depending on available memory
  integer :: n_energy_batch = -1 ! Batch size for 'optimal' wf printing
  logical :: use_hash_generation = .false. ! Control flag for using hash function in Hamiltonian generation

  namelist /selected_ci/ eps_var_sched, eps_var, eps_pt, eps_pt_big, eps_pt_big_energy, target_error, dump_wf_var, n_max_connections, n_states, n_mc, n_energy_batch, use_hash_generation

  ! Automatically choose HF det of chosen symmetry
  logical :: get_auto_hf = .false.
  integer :: hf_symmetry = 0
  integer,allocatable :: up(:),dn(:)
  integer :: lz
  logical :: g = .false.
  logical :: u = .false.
  integer,allocatable :: irreps(:),irrep_occs_up(:),irrep_occs_dn(:) ! For assigning starting det's occ orbs to be a set of numbers of irreps, e.g., 3 of irrep 1, 1 of irrep 2, etc.
  integer :: n_irrep
  namelist /hf_det/ up, dn, hf_symmetry, lz, g, u, n_irrep, irreps, irrep_occs_up, irrep_occs_dn

  ! HCI Natural Orbitals
  logical :: get_natorbs = .false.
  logical :: use_pt = .false.
  namelist /natorb/  get_natorbs, use_pt

  ! Green's Functions
  logical :: get_greens_function = .false.
  integer :: n_w = 20 ! number of frequencies to compute G0_pq(w) for
  real(rk) :: w_min, w_max ! min and max frequency range
  namelist /greens_function/  get_greens_function, n_w, w_min, w_max

  ! HCI Active Space (for using a smaller active space for the variational part
  ! of HCI)
  integer :: n_var_e_up, n_var_e_dn, n_var_orbs = 0
  integer,allocatable :: var_orbs(:)
  namelist /active_space/ n_var_e_up, n_var_e_dn, n_var_orbs, var_orbs

end module common_selected_ci
